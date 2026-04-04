import numpy as np

from mdcraft.data import Atoms

class Mesh:

	def __init__(self, Lx, Ly, Lz, D, kind='pcc'):
		self.Lx = Lx
		self.Ly = Ly
		self.Lz = Lz
		self.D  = D
		self.kind = kind

	def setup(self, \
		Bxmin = None, Bymin = None, Bzmin = None, \
		Bxmax = None, Bymax = None, Bzmax = None):
		# adjust sizes for PBC
		D  = self.D
		Lx = self.Lx
		Ly = self.Ly
		Lz = self.Lz 
		Nx = int(Lx/D)
		Ny = int(Ly/D)
		Nz = int(Lz/D)
		Lx = Nx*D
		Ly = Ny*D
		Lz = Nz*D
		self.Lx = Lx
		self.Ly = Ly 
		self.Lz = Lz

		if  Bxmin == None:
			Bxmin = -0.5*Lx
		if  Bxmax == None:
			Bxmax =  0.5*Lx

		if  Bymin == None:
			Bymin = -0.5*Ly
		if  Bymax == None:
			Bymax =  0.5*Ly

		if  Bzmin == None:
			Bzmin = -0.5*Lz
		if  Bzmax == None:
			Bzmax =  0.5*Lz

		Bx = Bxmax - Bxmin
		By = Bymax - Bymin
		Bz = Bzmax - Bzmin

		if Bx >= Lx:
			Bxmin = -0.5*Lx
			Bxmax =  0.5*Lx
		else:
			Nx = int(Bx/D)
			Bxmin = Bxmin + 0.5*(Bx - Nx*D)
			Bxmax = Bxmin + Nx*D
		if By >= Ly:
			Bymin = -0.5*Ly
			Bymax =  0.5*Ly
		else:
			Ny = int(By/D)
			Bymin = Bymin + 0.5*(By - Ny*D)
			Bymax = Bymin + Ny*D
		if Bz >= Lz:
			Bzmin = -0.5*Lz
			Bzmax =  0.5*Lz
		else:
			Nz = int(Bz/D)
			Bzmin = Bzmin + 0.5*(Bz - Nz*D)
			Bzmax = Bzmin + Nz*D
		print("MDbox sizes: Lx = ", Lx, ", Ly = ", Ly, ", Lz = ", Lz)
		print("Bxmin = ", Bxmin, ", Bymin = ", Bymin, ", Bzmin = ", Bzmin)
		print("Bxmax = ", Bxmax, ", Bymax = ", Bymax, ", Bzmax = ", Bzmax)
		atoms = Atoms(0)
		# evaluate mesh
		if self.kind == 'pcc':
			X, Y, Z = np.meshgrid(
			    np.linspace(Bxmin + 0.5*D, Bxmax - 0.5*D, Nx),
			    np.linspace(Bymin + 0.5*D, Bymax - 0.5*D, Ny),
			    np.linspace(Bzmin + 0.5*D, Bzmax - 0.5*D, Nz),
			    indexing = 'ij'
			)
			atoms.resize(X.size)
			atoms.data["r"]["x"] = X.ravel()
			atoms.data["r"]["y"] = Y.ravel()
			atoms.data["r"]["z"] = Z.ravel()

		if self.kind == 'bcc':
			X, Y, Z = np.meshgrid(
			    np.linspace(Bxmin + 0.25*D, Bxmax - 0.75*D, Nx),
			    np.linspace(Bymin + 0.25*D, Bymax - 0.75*D, Ny),
			    np.linspace(Bzmin + 0.25*D, Bzmax - 0.75*D, Nz),
			    indexing = 'ij'
			)
			atoms.resize(2*X.size)
			atoms.data["r"]["x"][0::2] = X.ravel()
			atoms.data["r"]["y"][0::2] = Y.ravel()
			atoms.data["r"]["z"][0::2] = Z.ravel()
			atoms.data["r"]["x"][1::2] = X.ravel() + 0.5*D
			atoms.data["r"]["y"][1::2] = Y.ravel() + 0.5*D
			atoms.data["r"]["z"][1::2] = Z.ravel() + 0.5*D

		if self.kind == 'fcc':
			X, Y, Z = np.meshgrid(
			    np.linspace(Bxmin + 0.25*D, Bxmax - 0.75*D, Nx),
			    np.linspace(Bymin + 0.25*D, Bymax - 0.75*D, Ny),
			    np.linspace(Bzmin + 0.25*D, Bzmax - 0.75*D, Nz),
			    indexing = 'ij'
			)
			atoms.resize(4*X.size)
			atoms.data["r"]["x"][0::4] = X.ravel()
			atoms.data["r"]["y"][0::4] = Y.ravel()
			atoms.data["r"]["z"][0::4] = Z.ravel()
			atoms.data["r"]["x"][1::4] = X.ravel() + 0.5*D
			atoms.data["r"]["y"][1::4] = Y.ravel() + 0.5*D
			atoms.data["r"]["z"][1::4] = Z.ravel()
			atoms.data["r"]["x"][2::4] = X.ravel() + 0.5*D
			atoms.data["r"]["y"][2::4] = Y.ravel()
			atoms.data["r"]["z"][2::4] = Z.ravel() + 0.5*D
			atoms.data["r"]["x"][3::4] = X.ravel() 
			atoms.data["r"]["y"][3::4] = Y.ravel() + 0.5*D
			atoms.data["r"]["z"][3::4] = Z.ravel() + 0.5*D

		atoms.data["kind"] = 1

		self.Bxmin = Bxmin
		self.Bxmax = Bxmax
		self.Bymin = Bymin
		self.Bymax = Bymax
		self.Bzmin = Bzmin
		self.Bzmax = Bzmax

		print("Total number of atoms: ", atoms.data.size)

		return atoms

	def box_sizes(self):
		return self.Lx, self.Ly, self.Lz

	def sample_sizes(self):
		return self.Bxmin, self.Bymin, self.Bzmin, \
		       self.Bxmax, self.Bymax, self.Bzmax

	# def load_from_hdf(self, hdf_path):

	# 	Lx = self.Lx
	# 	Ly = self.Ly
	# 	Lz = self.Lz

	# 	load_file = hdf.file(name=hdf_path, access='r')
	# 	numpy_array = load_file["storage"][:]
	# 	atoms = Atoms(numpy_array)
	# 	box = Box(-Lx/2, -Ly/2, -Lz/2,
	# 		       Lx/2,  Ly/2,  Lz/2)
	# 	return box, atoms
