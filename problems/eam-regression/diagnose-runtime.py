#!/usr/bin/env python3
"""Read-only installation report and Multi construction stress test; no private input."""
import argparse
import gc
import importlib
import json
from pathlib import Path
import sys
import numpy as np

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--serial',action='store_true')
p.add_argument('--iterations',type=int,default=1000)
p.add_argument('--threads',type=int,default=1)
a=p.parse_args()
if a.iterations<1 or a.threads<1:p.error('positive iterations and threads required')
comm=None;rank=0
if not a.serial:
    from mpi4py import MPI
    comm=MPI.COMM_WORLD;rank=comm.rank
from mdcraft.data import Atoms
from mdcraft.tools import Threads
from mdcraft.solver import Multi
from mdcraft.solver.thermostat import Langevin
from mdcraft.solver.boundary import Periodic
from mdcraft.lattice import Domain
from mdcraft.solver import potential
modules={}
for name in ('mdcraft.data._mdcraft_data','mdcraft.tools._mdcraft_tools',
             'mdcraft.solver._mdcraft_solver','mdcraft.solver.potential._mdcraft_solver_potential',
             'mdcraft.solver.thermostat._mdcraft_solver_thermostat'):
    module=importlib.import_module(name);modules[name]=str(Path(module.__file__).resolve())
if comm is not None:
    module=importlib.import_module('mdcraft.decomp._mdcraft_decomp');modules[module.__name__]=str(Path(module.__file__).resolve())
atom=Atoms(1);dtype=atom.data.dtype
maps=Path('/proc/self/maps')
libraries=sorted({line.split()[-1] for line in maps.read_text().splitlines()
                  if '/' in line and any(word in line for word in ('libmdcraft','libmpi','libmpich'))}) if maps.exists() else []
report=dict(rank=rank,python=sys.executable,modules=modules,libraries=libraries,
            atom_itemsize=dtype.itemsize,atom_fields={name:int(dtype.fields[name][1]) for name in dtype.names})
if comm is not None:report['mpi']=MPI.Get_library_version()
print(json.dumps(report,indent=2),flush=True)
threads=Threads(a.threads) if comm is None else Threads(a.threads,comm)
domain=Domain(-2,2,-2,2,-2,2);boundary=Periodic(domain)
coeff={name:np.zeros(12) for name in 'VWCDFG'}
coeff['V'][2]=.05;coeff['C'][2]=.1;coeff['F'][2]=.02
potentials=[getattr(potential,name)(Rcutoff=1.,aVr=.001,bVr=.01,rVr=.01,**coeff)
            for name in ('EAMW1046','EAMW0644','EAMMOD1044')]
for i in range(a.iterations):
    thermostat=Langevin(.5,300,.001)
    solver=Multi(boundary=boundary,thermostat=thermostat,nkinds=2,threads=threads)
    for x,y,pot in ((0,0,potentials[0]),(1,1,potentials[1]),(0,1,potentials[2])):
        solver.add_potential(x,y,pot)
    del solver,thermostat
    if i%100==0:gc.collect()
print(f'PASS rank={rank}: {a.iterations} Multi factories, no atoms/input file/decomposition',flush=True)
