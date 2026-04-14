import numpy as np
from mdcraft.data import Atoms
from mdcraft.lattice import Domain

import matplotlib.pyplot as plt

domain = Domain(
	xmin = -10, xmax = 10,
	ymin = -5,  ymax = 5,
	zmin = -5,  zmax = 5
)

domain.set_periodic(0, True)
domain.set_periodic(1, True)

print(domain.xmin)
print(domain.ymin)
print(domain.zmin)
print(domain.xmax)
print(domain.ymax)
print(domain.zmax)
print(domain.xsize)
print(domain.ysize)
print(domain.zsize)
print(domain.xperiodic)
print(domain.yperiodic)
print(domain.zperiodic)

X, Y, Z = np.meshgrid(
	np.linspace(-20, 20, 41),
	np.linspace( -6,  6, 13),
	np.linspace( -6,  6, 13),
)
atoms = Atoms(X.size)

atoms.data["r"]["x"] = X.ravel()
atoms.data["r"]["y"] = Y.ravel()
atoms.data["r"]["z"] = Z.ravel()

plt.scatter(atoms.data["r"]["x"], atoms.data["r"]["y"], color="black")

ibelong = np.where(domain.belong(atoms))
plt.scatter(atoms.data["r"]["x"][ibelong], atoms.data["r"]["y"][ibelong], color="red")

plt.show()

distance = domain.nearest_distance(atoms, 1)
plt.scatter(atoms.data["r"]["y"], distance)
plt.show()

