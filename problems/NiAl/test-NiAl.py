#!/usr/bin/env python3
"""Public NiAl reproduction of the C/H Multi/Langevin workflow.

python test-NiAl.py --serial --threads 2
mpiexec -n 2 python test-NiAl.py --threads 2
Use a matching serial/MPI build of MDcraft. Outputs go to --output.
"""
import argparse
import gc
from pathlib import Path
import weakref
import numpy as np
from mdcraft.data import Atoms
from mdcraft.lattice import Domain
from mdcraft.neibs import VerletList
from mdcraft.solver import Multi
from mdcraft.solver.boundary import Periodic
from mdcraft.solver.potential import EAM, Pair
from mdcraft.solver.stepper import Verlet
from mdcraft.solver.thermostat import Langevin
from mdcraft.tools import Threads


def make_solver(boundary, threads, potentials, temperature, beta, dt, thermostat_on):
    # C++ stores references: retain the thermostat beyond this factory call.
    thermostat = Langevin(beta, temperature, dt) if thermostat_on else None
    options = dict(boundary=boundary, threads=threads, nkinds=2)
    if thermostat is not None:
        options['thermostat'] = thermostat
    solver = Multi(**options)
    for i, j, potential in ((0, 0, potentials[0]), (1, 1, potentials[1]), (0, 1, potentials[2])):
        solver.add_potential(i, j, potential)
    return solver, thermostat


def initial_atoms(cells, cutoff, kbuf, temperature):
    basis = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
    positions = np.concatenate([basis + [x, y, z] for x in range(cells)
                                for y in range(cells) for z in range(cells)]) * .405
    length = cells * .405
    atoms = Atoms(len(positions))
    data = atoms.data
    data['uid'] = np.arange(1, len(positions) + 1)
    data['kind'] = (data['uid'] % 5 == 0).astype(np.uint16)  # 0=Al, 1=Ni
    data['m'] = np.where(data['kind'] == 0, 26.9815384, 58.6934)
    data['rcut'], data['rns'] = cutoff, kbuf * cutoff
    rng = np.random.default_rng(2026)
    velocities = rng.normal(size=(len(positions), 3)) * np.sqrt(.008314462175 * temperature / data['m'])[:, None]
    velocities -= np.sum(velocities * data['m'][:, None], axis=0) / data['m'].sum()
    for i, axis in enumerate('xyz'):
        data['r'][axis] = positions[:, i] - length / 2
        data['v'][axis] = velocities[:, i]
    return atoms, (-length / 2, length / 2) * 3


def lifetime_regression(boundary, threads, potentials, dt):
    # Mirrors the user's original factory: options and Langevin are local.
    options = dict(boundary=boundary, threads=threads, nkinds=2,
                   thermostat=Langevin(.5, 300, dt))
    reference = weakref.ref(options['thermostat'])
    solver = Multi(**options)
    del options
    gc.collect()
    assert reference() is not None, 'Multi holds a dangling thermostat: missing keep_alive'
    for i, j, p in ((0, 0, potentials[0]), (1, 1, potentials[1]), (0, 1, potentials[2])):
        solver.add_potential(i, j, p)
    return solver, reference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--serial', action='store_true', help='use a build with MPI disabled')
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--steps', type=int, default=20)
    parser.add_argument('--cells', type=int, default=4)
    parser.add_argument('--temperature', type=float, default=300)
    parser.add_argument('--dt', type=float, default=.0001)
    parser.add_argument('--beta', type=float, default=.5)
    parser.add_argument('--kbuf', type=float, default=1.25)
    parser.add_argument('--output', type=Path, default=Path('nial-output'))
    parser.add_argument('--lifetime-regression', action='store_true')
    args = parser.parse_args()
    if args.threads < 1 or args.steps < 1 or args.cells < 4 or args.dt <= 0 or args.kbuf <= 1 or args.temperature < 0:
        parser.error('positive threads/steps/dt, cells>=4, kbuf>1, temperature>=0 required')
    rank, size, comm = 0, 1, None
    if not args.serial:
        from mpi4py import MPI
        from mdcraft.decomp import VD3
        from mpi_output import global_thermodynamics, dump_frame
        comm = MPI.COMM_WORLD
        rank, size = comm.rank, comm.size
    threads = Threads(args.threads) if args.serial else Threads(args.threads, comm)
    directory = Path(__file__).resolve().parent
    potentials = (EAM(str(directory / 'Al.eam.alloy')), EAM(str(directory / 'Ni.eam.alloy')),
                  Pair(str(directory / 'NiAlpair.eam.alloy')))
    cutoff = max(p.rcut for p in potentials)
    atoms, bounds = initial_atoms(args.cells, cutoff, args.kbuf, args.temperature)
    if rank:
        atoms = Atoms(0)
    domain = Domain(*bounds)
    for axis in range(3):
        domain.set_periodic(axis, True)
    # A boundary binds itself to one solver; avoid sharing it with solver_cold.
    boundary, cold_boundary = Periodic(domain), Periodic(domain)
    solver, thermostat = make_solver(boundary, threads, potentials, args.temperature, args.beta, args.dt, True)
    cold, _ = make_solver(cold_boundary, threads, potentials, args.temperature, args.beta, args.dt, False)
    if args.lifetime_regression:
        solver, thermostat_reference = lifetime_regression(boundary, threads, potentials, args.dt)
        thermostat = None  # binding must retain the factory-created object
    stepper = Verlet(solver=solver, threads=threads)
    vd = None
    if comm is not None:
        vd = VD3(comm=comm, atoms=atoms, domain=domain, threads=threads, dimension=3, measurer='size')
        vd.prebalancing(verbose=False)
        vd.exchange()
        atoms = vd.locals
    local = VerletList(atoms=atoms, neighbors=atoms, domain=domain, threads=threads, half=False)
    halo = None if vd is None else VerletList(atoms=atoms, neighbors=vd.aliens, domain=domain, threads=threads, half=False)
    local.update()
    if halo is not None:
        halo.update()
    def arguments():
        return (atoms, atoms, local.get()) if vd is None else (vd, local.get(), halo.get())
    cold.prepare(*arguments())
    cold.charge(*arguments())  # EAM density must precede forces and virials
    cold.forces(*arguments())
    cold.virials(*arguments())
    if rank == 0:
        args.output.mkdir(parents=True, exist_ok=True)
        print(f'ranks={size} requested_threads={args.threads} actual_threads={threads.count}', flush=True)
    if comm is not None:
        comm.Barrier()
    for step in range(1, args.steps + 1):
        step_args = (atoms, local.get(), args.dt) if vd is None else (vd, local.get(), halo.get(), args.dt)
        stepper.make_step(*step_args)
        data = atoms.data
        finite = all(np.isfinite(data[field][axis]).all() for field in ('r', 'v', 'f') for axis in 'xyz')
        if comm is not None:
            finite = comm.allreduce(finite, op=MPI.LAND)
        if not finite:
            raise RuntimeError(f'nonfinite state at step {step}')
        # Match main.py: output first, rebuild lists for the next step afterwards.
        solver.virials(*arguments())
        local.update()
        if halo is not None:
            halo.update()
    count = len(atoms.data)
    if comm is not None:
        count = comm.allreduce(count, op=MPI.SUM)
        values = global_thermodynamics(atoms, (bounds[1]-bounds[0])**3, comm)
        if not np.isfinite(values).all():
            raise RuntimeError('nonfinite thermodynamics')
        dump_frame(atoms, comm, bounds, args.output, args.steps, 4 * args.cells**3)
    if count != 4 * args.cells**3:
        raise RuntimeError(f'atom count changed: {count}')
    if rank == 0:
        print(f'PASS atoms={count} steps={args.steps} lifetime_regression={args.lifetime_regression}', flush=True)


if __name__ == '__main__':
    main()
