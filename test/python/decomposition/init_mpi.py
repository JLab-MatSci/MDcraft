from mpi4py import MPI
import numpy as np

print(MPI.Get_library_version())

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
print("load mpi4py", rank, size)

from mdcraft.data import Atoms
from mdcraft.tools import Threads
from mdcraft.lattice import Domain
from mdcraft.decomp import VD3

domain = Domain(
	xmin = -5, xmax = 5,
	ymin = -5, ymax = 5,
	zmin = -5, zmax = 5
)

X, Y, Z = np.meshgrid(
	np.linspace(domain.xmin, domain.xmax, 11),
	np.linspace(domain.ymin, domain.ymax, 11),
	np.linspace(domain.zmin, domain.zmax, 11),
)
atoms = Atoms(X.size)

atoms.data["r"]["x"] = X.ravel()
atoms.data["r"]["y"] = Y.ravel()
atoms.data["r"]["z"] = Z.ravel()

decomp = VD3(
	comm = comm, 
	atoms = atoms,
	domain = domain
)

print("decomposition size", decomp.size)
print("decomposition rank", decomp.rank)
