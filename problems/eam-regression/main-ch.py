#!/usr/bin/env python3
"""C/H molecular dynamics: VD3 MPI decomposition + threads per MPI rank.

Run as `srun python -u main.py` (one MPI rank / 128-core node), or with
multiple MPI ranks/node. All physical parameters retain main.py defaults.

This script assumes the corrected Multi and EAM C/H implementation and that
VD3 halo exchanges preserve all required density fields (nc, nc_mod).
"""

import argparse
from pathlib import Path
import socket
import time as walltime

import numpy as np
from mpi4py import MPI

from mdcraft.data import Atoms
from mdcraft.decomp import VD3
from mdcraft.lattice import Domain
from mdcraft.neibs import VerletList
from mdcraft.solver import Multi as Solver
from mdcraft.solver.boundary import Periodic
from mdcraft.solver.stepper import Verlet as Stepper
from mdcraft.solver.thermostat import Langevin
from mdcraft.tools import Threads

from mpi_output import dump_frame, global_thermodynamics


LOG_HEADER = (
    "%8s %15s %15s %15s %15s %15s %15s %15s %15s\n"
    % ("step", "time", "potentialEnergy", "totalEnergy", "kineticEnergy",
       "temperature", "virial", "pressure", "density")
)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="C-H EAM on multi-node MPI + threads (VD3)")
    p.add_argument("--data", type=Path, default=Path("LAYER_15_25_15.lmp"),
                   help="LAMMPS data file, read only on rank 0")
    p.add_argument("--inputs", type=Path, default=Path("inputs"),
                   help="directory with EAMW1046.dat, EAMW0644.dat, EAMW1044.dat")
    p.add_argument("--threads-per-node", type=int, default=128,
                   help="total CPU worker threads across MPI ranks on one compute node")
    p.add_argument("--threads", type=int, default=None,
                   help="threads per MPI rank (default: threads-per-node // ranks-on-node)")
    p.add_argument("--like-pairs", choices=("original", "modified"), default="original",
                   help="original uses EAMW1046/EAMW0644 as in main.py; modified uses EAMWMOD1046/EAMWMOD0644")
    p.add_argument("--temperature", type=float, default=2000.0, help="Langevin bath temperature, K")
    p.add_argument("--beta", type=float, default=0.5, help="Langevin friction, ps^-1")
    p.add_argument("--dt", type=float, default=1e-3, help="timestep, ps")
    p.add_argument("--tend", type=float, default=1200.0, help="end time, ps")
    p.add_argument("--steps", type=int, default=None, help="override tend/dt step count")
    p.add_argument("--kbuf", type=float, default=1.7, help="Verlet radius / cutoff")
    p.add_argument("--nlist-every", type=int, default=1,
                   help="maximum number of steps between rebuilds; 0=only displacement decisions")
    p.add_argument("--log-every", type=int, default=2, help="thermodynamic log frequency, in MD steps (0=off)")
    p.add_argument("--dump-every", type=int, default=1000,
                   help="MPI-gathered LAMMPS dump frequency; 0=off (original was every step)")
    p.add_argument("--log-dir", type=Path, default=Path("log"))
    p.add_argument("--dump-dir", type=Path, default=Path("lammps_dump"))
    p.add_argument("--dump-coordinate-scale", type=float, default=10.0,
                   help="multiply native MDCraft coordinates/box by this value (default nm -> Angstrom)")
    p.add_argument("--no-initial-dump", action="store_true")
    args = p.parse_args(argv)
    if args.threads_per_node < 1 or (args.threads is not None and args.threads < 1):
        p.error("thread counts must be positive")
    if args.dt <= 0 or args.tend < 0 or (args.steps is not None and args.steps < 0):
        p.error("dt must be positive, tend/steps non-negative")
    if args.kbuf <= 1 or args.nlist_every < 0 or args.log_every < 0 or args.dump_every < 0:
        p.error("kbuf > 1 and non-negative output/rebuild intervals required")
    if args.dump_coordinate_scale <= 0:
        p.error("dump-coordinate-scale must be positive")
    return args


def factor_grid(size, lengths):
    """Find factorization nx*ny*nz=size with near-isotropic subdomain sizes."""
    if size < 1 or np.any(np.asarray(lengths) <= 0):
        raise ValueError("size and box dimensions must be positive")
    best, best_error = None, float("inf")
    for nx in range(1, size + 1):
        if size % nx:
            continue
        rest = size // nx
        for ny in range(1, rest + 1):
            if rest % ny:
                continue
            nz = rest // ny
            cells = np.asarray(lengths, dtype=float) / np.array([nx, ny, nz])
            logs = np.log(cells)
            error = float(np.sum((logs - logs.mean()) ** 2))
            if error < best_error:
                best, best_error = (nx, ny, nz), error
    return best


def make_centers(size, bounds):
    """Identical deterministic initial Voronoi centers on every MPI rank."""
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    lengths = (xmax - xmin, ymax - ymin, zmax - zmin)
    nx, ny, nz = factor_grid(size, lengths)
    centers = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                centers.append([
                    xmin + (ix + 0.5) * lengths[0] / nx,
                    ymin + (iy + 0.5) * lengths[1] / ny,
                    zmin + (iz + 0.5) * lengths[2] / nz,
                ])
    assert len(centers) == size
    return centers, (nx, ny, nz)


def determine_threads(comm, requested, budget):
    """Use 128 threads/node, not 128 threads *per rank* when ranks share a node."""
    node_comm = comm.Split_type(MPI.COMM_TYPE_SHARED, key=comm.Get_rank())
    local_ranks = node_comm.Get_size()
    node_comm.Free()
    if requested is None:
        if budget % local_ranks:
            raise ValueError(
                f"{budget} threads/node cannot be evenly divided among "
                f"{local_ranks} MPI ranks on this node; use --threads explicitly"
            )
        per_rank = budget // local_ranks
    else:
        per_rank = requested
    if per_rank * local_ranks > budget:
        raise ValueError(
            f"oversubscription: {local_ranks} ranks/node x {per_rank} threads/rank "
            f"> {budget} threads/node"
        )
    return per_rank, local_ranks


def load_potentials(input_dir, like_pairs):
    from eamReader import eamReader
    if like_pairs == "modified":
        from mdcraft.solver.potential import EAMWMOD1046 as PotentialCC
        from mdcraft.solver.potential import EAMWMOD0644 as PotentialHH
    else:
        from mdcraft.solver.potential import EAMW1046 as PotentialCC
        from mdcraft.solver.potential import EAMW0644 as PotentialHH
    from mdcraft.solver.potential import EAMMOD1044 as PotentialCH

    def construct(cls, filename):
        values = eamReader(str(input_dir / filename))
        rcut, aVr, bVr, rVr, V, W, C, D, F, G = values
        return cls(Rcutoff=rcut, aVr=aVr, bVr=bVr, rVr=rVr,
                   V=V, W=W, C=C, D=D, F=F, G=G), float(rcut)

    cc, cc_cut = construct(PotentialCC, "EAMW1046.dat")
    hh, hh_cut = construct(PotentialHH, "EAMW0644.dat")
    ch, ch_cut = construct(PotentialCH, "EAMW1044.dat")
    return (cc, hh, ch), max(cc_cut, hh_cut, ch_cut), (cc_cut, hh_cut, ch_cut)


def read_initial_atoms(data_path, rcut, kbuf):
    """Root-only ASE conversion; swap ASE kinds 1(C)/0(H) -> MDCraft 0(C)/1(H)."""
    from ase.io import read
    from ase2mdcraft import convert_domain, convert_atoms

    ase_atoms = read(str(data_path), format="lammps-data")
    initial_domain = convert_domain(ase_atoms)
    bounds = (initial_domain.xmin, initial_domain.xmax,
              initial_domain.ymin, initial_domain.ymax,
              initial_domain.zmin, initial_domain.zmax)
    atoms = convert_atoms(ase_atoms)
    kinds = np.asarray(atoms.data["kind"])
    if not np.all((kinds == 0) | (kinds == 1)):
        raise ValueError("expected ASE kinds 1=carbon, 0=hydrogen; unexpected kind values")
    n_c = int(np.count_nonzero(kinds == 1))
    n_h = int(np.count_nonzero(kinds == 0))
    atoms.data["kind"] = 1 - kinds
    atoms.data["uid"] = np.arange(1, atoms.data.size + 1, dtype=np.uint32)
    atoms.data["rcut"] = rcut
    atoms.data["rns"] = kbuf * rcut
    return atoms, bounds, (n_c, n_h)


def make_solver(periodic, threads, potentials, temp, beta, dt, thermostat_on):
    options = {"boundary": periodic, "threads": threads, "nkinds": 2}
    if thermostat_on:
        options["thermostat"] = Langevin(
            beta=beta, temperature=temp, time_step=dt,
            heat_x=1, heat_y=1, heat_z=1, Ux=0, Uy=0, Uz=0,
        )
    solver = Solver(**options)
    cc, hh, ch = potentials
    solver.add_potential(0, 0, cc)
    solver.add_potential(1, 1, hh)
    solver.add_potential(0, 1, ch)
    # Pair stores a C++ reference: the caller must retain the Python owner.
    return solver, options.get("thermostat")


def format_log(step, time_ps, numbers):
    return "%8d %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n" % (
        step, time_ps, *numbers
    )


def main(argv=None):
    args = parse_args(argv)
    comm = MPI.COMM_WORLD
    rank, size = comm.Get_rank(), comm.Get_size()
    nthreads, local_ranks = determine_threads(comm, args.threads, args.threads_per_node)
    rank_layouts = comm.allgather((socket.gethostname(), local_ranks, nthreads))
    if rank == 0:
        layout = sorted(set(rank_layouts))
        print(f"MPI size={size}, threads/rank={nthreads}, threads/node budget={args.threads_per_node}; "
              f"per-node layouts: {layout}", flush=True)

    threads = Threads(nthreads, comm)  # actual MDCraft MPI overload
    actual_threads = comm.allgather(threads.count)
    if rank == 0:
        print(f"Actual MDCraft threads per rank: {actual_threads}", flush=True)
    start = walltime.perf_counter()
    potentials, rcut, cutoffs = load_potentials(args.inputs, args.like_pairs)

    atoms = Atoms(0)
    bounds = None
    root_meta = None
    init_error = None
    if rank == 0:
        try:
            atoms, bounds, root_meta = read_initial_atoms(args.data, rcut, args.kbuf)
        except Exception as exc:
            init_error = f"rank 0 failed reading {args.data}: {exc}"
    init_error = comm.bcast(init_error, root=0)
    if init_error:
        raise RuntimeError(init_error)
    bounds = comm.bcast(bounds, root=0)
    counts = comm.bcast(root_meta, root=0)

    domain = Domain(*bounds)
    for axis in range(3):
        domain.set_periodic(axis, True)
    boundary = Periodic(domain)
    volume = (bounds[1] - bounds[0]) * (bounds[3] - bounds[2]) * (bounds[5] - bounds[4])
    centers, grid = make_centers(size, bounds)
    vd = VD3(comm=comm, atoms=atoms, domain=domain, dimension=3,
             centers=centers, centroidal=0.25, measurer="size", threads=threads)
    vd.prebalancing(verbose=False)

    # EAM densities require directed full lists, including all local <-> local pairs.
    nlist_local = VerletList(atoms=vd.locals, neighbors=vd.locals,
                             domain=domain, threads=threads, half=False)
    nlist_halo = VerletList(atoms=vd.locals, neighbors=vd.aliens,
                            domain=domain, threads=threads, half=False)
    nlist_local.update()
    nlist_halo.update()

    solver, thermostat = make_solver(boundary, threads, potentials, args.temperature, args.beta, args.dt, True)
    boundary_cold = Periodic(domain)
    solver_cold, thermostat_cold = make_solver(boundary_cold, threads, potentials, args.temperature, args.beta, args.dt, False)
    stepper = Stepper(solver=solver, threads=threads)
    solver_cold.prepare(vd, nlist_local.get(), nlist_halo.get())
    solver_cold.charge(vd, nlist_local.get(), nlist_halo.get())
    solver_cold.forces(vd, nlist_local.get(), nlist_halo.get())
    solver_cold.virials(vd, nlist_local.get(), nlist_halo.get())

    logfile = None
    if rank == 0:
        if args.log_every:
            args.log_dir.mkdir(parents=True, exist_ok=True)
            logfile = (args.log_dir / "log-prepare.dat").open("w", encoding="utf-8")
            logfile.write(LOG_HEADER)
        if args.dump_every:
            args.dump_dir.mkdir(parents=True, exist_ok=True)
        print(f"Atoms total={sum(counts)} (C={counts[0]}, H={counts[1]}), "
              f"VD3 grid={grid}, cutoffs={cutoffs}, kbuf={args.kbuf}, "
              f"pair EAM={args.like_pairs}", flush=True)

    try:
        if args.log_every:
            values = global_thermodynamics(vd.locals, volume, comm)
            if rank == 0:
                logfile.write(format_log(0, 0.0, values))
                logfile.flush()
        if args.dump_every and not args.no_initial_dump:
            dump_frame(vd.locals, comm, bounds, args.dump_dir, step=0,
                       expected_atoms=sum(counts), coordinate_scale=args.dump_coordinate_scale)

        steps = args.steps if args.steps is not None else int(round(args.tend / args.dt))
        rebuilds = 0
        for step in range(1, steps + 1):
            stepper.make_step(vd, nlist_local.get(), nlist_halo.get(), args.dt)

            # Diagnose with the SAME neighbor lists used for this step's forces.
            # Rebuild for the *next* step only after virials/output, as in VD3 example.
            if args.log_every and step % args.log_every == 0:
                solver.virials(vd, nlist_local.get(), nlist_halo.get())
                values = global_thermodynamics(vd.locals, volume, comm)
                if rank == 0:
                    logfile.write(format_log(step, step * args.dt, values))
                    logfile.flush()
                    if step % max(args.log_every, 200) == 0:
                        print(f"step={step} time={step * args.dt:.5f} ps "
                              f"elapsed={walltime.perf_counter() - start:.1f} s", flush=True)
            if args.dump_every and step % args.dump_every == 0:
                dump_frame(vd.locals, comm, bounds, args.dump_dir, step=step,
                           expected_atoms=sum(counts), coordinate_scale=args.dump_coordinate_scale)

            # Both displacement checks are evaluated; MPI uses a collective decision.
            local_decision = bool(nlist_local.decide())
            halo_decision = bool(nlist_halo.decide())
            forced = args.nlist_every != 0 and step % args.nlist_every == 0
            need_rebuild = comm.allreduce(local_decision or halo_decision or forced, op=MPI.LOR)
            if need_rebuild:
                rebuilds += 1
                # Refresh ownership and halo routes before rebuilding neighbour lists.
                vd.update(verbose=False)
                nlist_local.update()
                nlist_halo.update()

        elapsed = comm.reduce(walltime.perf_counter() - start, op=MPI.MAX, root=0)
        n_rebuilt = comm.reduce(rebuilds, op=MPI.SUM, root=0)
        if rank == 0:
            print(f"Finished: ranks={size}, threads/rank={nthreads}, steps={steps}, "
                  f"max wall={elapsed:.3f} s, rebuilds/rank={n_rebuilt / size:.1f}", flush=True)
    finally:
        if logfile is not None:
            logfile.close()


if __name__ == "__main__":
    main()
