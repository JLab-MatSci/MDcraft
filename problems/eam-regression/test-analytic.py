#!/usr/bin/env python3
"""Exercise every supplied analytic EAM variant with public synthetic coefficients.

The coefficients test implementation consistency; they are NOT a C/H force field.
Run serial for finite-difference forces, or with mpiexec for VD3/halo execution.
"""
import argparse
import json
import numpy as np
from mdcraft.data import Atoms
from mdcraft.lattice import Domain
from mdcraft.neibs import VerletList
from mdcraft.solver import Multi
from mdcraft.solver.boundary import Periodic
from mdcraft.solver.stepper import Verlet
from mdcraft.solver import potential as pots
from mdcraft.tools import Threads


def coefficients():
    rng = np.random.default_rng(2026)
    arrays = {name: rng.uniform(-.001, .001, 12) for name in 'VWCDFG'}
    for values in arrays.values(): values[:2] = 0.0
    arrays['V'][2] = .05
    arrays['C'][2] = .1
    arrays['F'][2] = .02
    return arrays


def construct(name, arrays=None):
    cls = getattr(pots, name)
    arrays = coefficients() if arrays is None else arrays
    keys = 'VW' if name == 'EAMV1622' else 'VWCDFG'
    return cls(Rcutoff=1., aVr=.001, bVr=.01, rVr=.01,
               **{key: arrays[key] for key in keys})


def check_bindings():
    tabulated = pots.EAM(Rcutoff=1., pairPotential=np.zeros(8), chargeDensity=np.ones(8),
                         embeddingEnergy=np.linspace(0, 1, 8), NrPot=8, drPot=.2, NrEmb=8, drEmb=.2)
    rho = np.array([-1., 0., .3, 1.4, 100.])
    np.testing.assert_allclose(tabulated.embeddingEnergy(rho), rho / 1.4, rtol=1e-12, atol=1e-12)
    try:
        pots.EAM(Rcutoff=1., pairPotential=np.zeros(3), chargeDensity=np.ones(8),
                 embeddingEnergy=np.zeros(8), NrPot=8, drPot=.2, NrEmb=8, drEmb=.2)
    except ValueError:
        pass
    else:
        raise AssertionError('EAM accepted an undersized radial table')
    names = ['EAMW1046' , 'EAMW0644', 'EAMWMOD1046', 'EAMWMOD0644', 'EAMMOD1044', 'EAMV1622']
    for name in names:
        p = construct(name)
        assert np.isfinite(p.value(.3)), name
        bad = coefficients(); bad['W'] = np.zeros(11)
        try:
            construct(name, bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f'{name} accepted mismatched arrays')
        bad = coefficients(); bad['V'] = np.zeros(24)[::2]
        try:
            construct(name, bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f'{name} accepted a strided coefficient array')


def run(model, args, comm):
    rank = 0 if comm is None else comm.rank
    threads = Threads(args.threads) if comm is None else Threads(args.threads, comm)
    names = {'original': ('EAMW1046','EAMW0644'),
             'modified': ('EAMWMOD1046','EAMWMOD0644'),
             'pair': ('EAMV1622','EAMV1622')}[model]
    potentials = [construct(name) for name in (*names, 'EAMMOD1044')]
    atoms = Atoms(4 if rank == 0 else 0)
    if rank == 0:
        data = atoms.data
        data['uid'] = np.arange(1,5)
        data['kind'] = [0,1,0,1]
        data['m'] = [12.,1.,12.,1.]
        data['rcut'], data['rns'] = 1., 1.25
        for axis in 'xyz':
            data['r'][axis] = 0; data['v'][axis] = 0
        data['r']['x'] = [-.4,-.1,.2,.5]
    domain = Domain(-2,2,-2,2,-2,2)
    for axis in range(3): domain.set_periodic(axis,True)
    boundary = Periodic(domain)
    solver = Multi(boundary=boundary, nkinds=2, threads=threads)
    for i,j,p in [(0,0,potentials[0]),(1,1,potentials[1]),(0,1,potentials[2])]:
        solver.add_potential(i,j,p)
    vd = None
    if comm is not None:
        from mdcraft.decomp import VD3
        centers = [[-2 + (i + .5) * 4 / comm.size, 0, 0] for i in range(comm.size)]
        vd = VD3(comm=comm, atoms=atoms, domain=domain, threads=threads, dimension=1, centers=centers)
        vd.exchange(); atoms = vd.locals
    local = VerletList(atoms=atoms,neighbors=atoms,domain=domain,threads=threads,half=False)
    halo = None if vd is None else VerletList(atoms=atoms,neighbors=vd.aliens,domain=domain,threads=threads,half=False)
    local.update()
    if halo is not None: halo.update()
    def inputs(): return (atoms,atoms,local.get()) if vd is None else (vd,local.get(),halo.get())
    def energy():
        solver.prepare(*inputs()); solver.charge(*inputs()); solver.forces(*inputs()); solver.virials(*inputs())
        value = float(np.sum(atoms.data['Ep']+atoms.data['Em']))
        return value if comm is None else comm.allreduce(value)
    initial = energy()
    forces = atoms.data['f']['x'].copy()
    solver.virials(*inputs())
    repeated = float(np.sum(atoms.data['Ep']+atoms.data['Em']))
    if comm is not None: repeated = comm.allreduce(repeated)
    np.testing.assert_allclose(repeated, initial, rtol=1e-12, atol=1e-12)
    if comm is None:
        virial = float(np.sum(atoms.data['V']['xx'] + atoms.data['V']['yy'] + atoms.data['V']['zz']))
        coordinates = {axis: atoms.data['r'][axis].copy() for axis in 'xyz'}
        h = 1e-6
        for axis in 'xyz': atoms.data['r'][axis] = coordinates[axis] * (1 + h)
        ep = energy()
        for axis in 'xyz': atoms.data['r'][axis] = coordinates[axis] * (1 - h)
        em = energy()
        for axis in 'xyz': atoms.data['r'][axis] = coordinates[axis]
        np.testing.assert_allclose(virial, (ep-em)/(2*h), rtol=2e-5, atol=1e-7,
                                   err_msg=f'{model} virial/strain gradient')
        energy()
        for i in range(len(atoms.data)):
            x = float(atoms.data['r']['x'][i]); h=1e-6
            atoms.data['r']['x'][i]=x+h; ep=energy()
            atoms.data['r']['x'][i]=x-h; em=energy()
            atoms.data['r']['x'][i]=x
            np.testing.assert_allclose(forces[i], -(ep-em)/(2*h),rtol=2e-5,atol=1e-7,
                                       err_msg=f'{model} force/energy gradient atom {i}')
    energy()
    stepper = Verlet(solver=solver,threads=threads)
    for step in range(args.steps):
        stepper.make_step(*( (atoms,local.get(),1e-5) if vd is None else (vd,local.get(),halo.get(),1e-5)))
        for field in ('r','v','f'):
            for axis in 'xyz':
                assert np.isfinite(atoms.data[field][axis]).all(), (model,step,field)
    final = energy()
    owned = [(int(a['uid']),float(a['r']['x']),float(a['f']['x'])) for a in atoms.data]
    if comm is not None: owned = [a for chunk in comm.allgather(owned) for a in chunk]
    if rank == 0:
        assert len(owned)==4 and len({a[0] for a in owned})==4
        print(json.dumps(dict(model=model,initial_energy=initial,final_energy=final,
                              atoms=sorted(owned),threads=threads.count)),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--serial',action='store_true');p.add_argument('--threads',type=int,default=1)
    p.add_argument('--steps',type=int,default=3);args=p.parse_args()
    comm=None
    if not args.serial:
        from mpi4py import MPI
        comm=MPI.COMM_WORLD
    check_bindings()
    for model in ('original','modified','pair'):run(model,args,comm)


if __name__=='__main__':main()
