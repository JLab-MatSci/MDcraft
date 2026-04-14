import argparse
import os
import time as Time
from pathlib import Path

import numpy as np
from mpi4py import MPI

import analysis
from mdcraft.data import Atoms
from mdcraft.decomp import VD3
from mdcraft.lattice import Domain
from mdcraft.neibs import VerletList
from mdcraft.solver import Single as Solver
from mdcraft.solver.boundary import Periodic
from mdcraft.solver.potential import LJsCut
from mdcraft.solver.stepper import Verlet as Stepper
from mdcraft.tools import Threads
from read_lammps_coords import read_lammps_data


def parse_bool(value):
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "on"}:
        return True
    if text in {"0", "false", "no", "off"}:
        return False
    raise argparse.ArgumentTypeError(f"Invalid boolean value: {value}")


def parse_args():
    parser = argparse.ArgumentParser(description="Argon ball MPI run for neighbor-list validation.")
    parser.add_argument("--half", type=parse_bool, default=True, help="Use half neighbor list for local-local pairs.")
    parser.add_argument("--threads", type=int, default=1, help="Number of CPU threads per MPI rank.")
    parser.add_argument("--steps", type=int, default=3000, help="Number of MD steps after the initial prepare stage.")
    parser.add_argument("--write-log", type=parse_bool, default=False, help="Write thermodynamic output to a log file.")
    parser.add_argument("--output-every", type=int, default=20, help="Write thermodynamic output every N steps.")
    parser.add_argument("--kbuf", type=float, default=10.125 / 8.125, help="Neighbor-list buffer factor.")
    parser.add_argument("--dt", type=float, default=4e-3, help="MD timestep in ps.")
    parser.add_argument("--seed", type=int, default=12345, help="Seed for deterministic decomposition centers.")
    parser.add_argument("--log-dir", default="log", help="Directory for output logs.")
    return parser.parse_args()

def ensure_log_dir(path, rank):
    if rank == 0:
        os.makedirs(path, exist_ok=True)


def main():
    args = parse_args()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    ar_mass = 40.0
    eV_NA = 1.602176634 * 6.02214076
    a_lmp = 0.01068764262454299
    s_lmp = 3.3841043857528683
    r_lmp = 8.125

    start = Time.time()

    potential = LJsCut(
        a=a_lmp * eV_NA * 10.0,
        s=s_lmp * 0.1,
        rcut=r_lmp * 0.1,
    )

    atoms = Atoms(0)
    if rank == 0:
        data_path = Path(__file__).with_name("LMP_MEGA_BALL_Ar")
        atoms, xmin, xmax, ymin, ymax, zmin, zmax = read_lammps_data(data_path, atoms)
        atoms.data["rcut"] = potential.rcut
        atoms.data["rns"] = args.kbuf * atoms.data["rcut"]
        atoms.data["m"] = ar_mass
        atoms.data["v"]["x"] = 0.0
        atoms.data["v"]["y"] = 0.0
        atoms.data["v"]["z"] = 0.0
        atoms.data["uid"] = np.arange(atoms.data.size, dtype=np.uint32)

    xmin, xmax = -16.08947906494, 48.26843719482
    ymin, ymax = -16.08947906494, 48.26843719482
    zmin, zmax = -16.08947906494, 48.26843719482

    volume = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)

    domain = Domain(xmin, xmax, ymin, ymax, zmin, zmax)
    domain.set_periodic(0, True)
    domain.set_periodic(1, True)
    domain.set_periodic(2, True)
    periodic = Periodic(domain)

    R = 20.7
    xc = yc = zc = (xmin + xmax) / 2.0

    centers = []
    while len(centers) < size:
        x = np.random.uniform(low=-R, high=R, size=1)[0]
        y = np.random.uniform(low=-R, high=R, size=1)[0]
        z = np.random.uniform(low=-R, high=R, size=1)[0]
        if x*x + y*y + z*z < R*R:
            centers.append([x + xc, y + yc, z + zc])

    vd = VD3(
        comm=comm,
        atoms=atoms,
        domain=domain,
        dimension=3,
        centers=centers,
        centroidal=0.25,
        measurer="size",
    )
    vd.prebalancing(verbose=False)

    solver = Solver(
        potential=potential,
        boundary=periodic,
    )
    solver_cold = Solver(
        potential=potential,
        boundary=periodic,
    )
    stepper = Stepper(
        solver=solver,
    )

    if args.write_log:
        ensure_log_dir(args.log_dir, rank)
        comm.Barrier()

        if rank == 0:
            logfile = open(Path(args.log_dir) / "log-prepare.dat", "w", encoding="utf-8")
            logfile.write(
                "%06s %15s %15s %15s %15s %15s %15s %15s %15s\n"
                % (
                    "step",
                    "time",
                    "potentialEnergy",
                    "totalEnergy",
                    "kineticEnergy",
                    "temperature",
                    "virial",
                    "pressure",
                    "density",
                )
            )
        else:
            logfile = None

    nlist1 = VerletList(
        atoms=vd.locals,
        neighbors=vd.locals,
        domain=domain,
        half=args.half,
    )
    nlist2 = VerletList(
        atoms=vd.locals,
        neighbors=vd.aliens,
        domain=domain,
        half=False,
    )

    nlist1.update()
    nlist2.update()

    solver_cold.prepare(vd, nlist1.get(), nlist2.get())
    solver_cold.forces(vd, nlist1.get(), nlist2.get())
    solver_cold.virials(vd, nlist1.get(), nlist2.get())

    energy_potential, energy_total, kinetic_energy, temperature, virial, pressure, density = (
        analysis.averageValues_mpi(vd.locals, volume)
    )

    if args.write_log:
        if rank == 0:
            logfile.write(
                "%06d %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n"
                % (
                    0,
                    0.0,
                    energy_potential,
                    energy_total,
                    kinetic_energy,
                    temperature,
                    virial,
                    pressure,
                    density,
                )
            )
            logfile.flush()
            print(
                "initial: half=%s threads=%d ranks=%d Ep=% .8e virial=% .8e"
                % (args.half, args.threads, size, energy_potential, virial),
                flush=True,
            )

    step = 0
    time_ps = 0.0
    rebuilds = 0

    while step < args.steps:
        stepper.make_step(vd, nlist1.get(), nlist2.get(), args.dt)

        time_ps += args.dt
        step += 1

        if args.write_log:
            if step % args.output_every == 0:
                solver.virials(vd, nlist1.get(), nlist2.get())

                energy_potential, energy_total, kinetic_energy, temperature, virial, pressure, density = (
                    analysis.averageValues_mpi(vd.locals, volume)
                )

                if rank == 0:
                    logfile.write(
                        "%06d %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n"
                        % (
                            step,
                            time_ps,
                            energy_potential,
                            energy_total,
                            kinetic_energy,
                            temperature,
                            virial,
                            pressure,
                            density,
                        )
                    )
                    logfile.flush()

        local_rebuild = nlist1.decide() or nlist2.decide()
        need_rebuild = comm.allreduce(local_rebuild, op=MPI.LOR)
        # need_rebuild = (step % 20 == 0)
        if need_rebuild:
            rebuilds += 1
            nlist1.update()
            nlist2.update()

    elapsed = Time.time() - start
    max_elapsed = comm.reduce(elapsed, op=MPI.MAX, root=0)
    total_rebuilds = comm.reduce(rebuilds, op=MPI.SUM, root=0)

    if rank == 0:
        print(
            "summary: half=%s threads=%d ranks=%d steps=%d wall=%.6f rebuilds=%d"
            % (args.half, args.threads, size, args.steps, max_elapsed, total_rebuilds / size),
            flush=True,
        )
        if args.write_log:
            logfile.close()


if __name__ == "__main__":
    main()
