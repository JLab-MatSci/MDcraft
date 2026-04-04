import numpy as np
from mpi4py import MPI
from scipy.ndimage import gaussian_filter1d

Kb       = 0.008314462175# kJ/mol/K
m_a      = 1.6605390666e-24  # g
N_a      = 6.02214076 * 1e23
P_factor = 100*10**4/N_a # from Molecular dynamic units to kPa

def averageValues_mpi(atoms, volume):
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	size = comm.Get_size()

	local_px = np.sum(atoms.data["m"] * atoms.data["v"]["x"])
	local_py = np.sum(atoms.data["m"] * atoms.data["v"]["y"])
	local_pz = np.sum(atoms.data["m"] * atoms.data["v"]["z"])
	local_mass = np.sum(atoms.data["m"])

	local_mom_buf = np.array([local_px, local_py, local_pz, local_mass])
	global_mom_buf = np.zeros(4)

	comm.Allreduce(local_mom_buf, global_mom_buf, op=MPI.SUM)

	# if rank == 0:
	total_px, total_py, total_pz, total_mass = global_mom_buf

	v_cm_x = total_px / total_mass
	v_cm_y = total_py / total_mass
	v_cm_z = total_pz / total_mass

	dv_x = atoms.data["v"]["x"] - v_cm_x
	dv_y = atoms.data["v"]["y"] - v_cm_y
	dv_z = atoms.data["v"]["z"] - v_cm_z

	V2_thermal = dv_x**2 + dv_y**2 + dv_z**2
	
	local_2Ek = np.sum(V2_thermal * atoms.data["m"])

	local_Ep = np.sum(atoms.data["Ep"])
	try:
		local_Ep += np.sum(atoms.data["Em"])
	except (ValueError, KeyError):
		pass
	
	local_virial = np.sum(atoms.data["V"]["xx"] + atoms.data["V"]["yy"] + atoms.data["V"]["zz"])
	local_N = float(atoms.data["m"].size)

	# print(f"local_N = {local_N} at rank {rank}")

	local_stats = np.array([local_2Ek, local_Ep, local_virial, local_mass, local_N])
	global_stats = np.zeros(5)

	comm.Allreduce(local_stats, global_stats, op=MPI.SUM)

	two_kinetic_energies = global_stats[0]
	energy_potential = global_stats[1]
	virial = global_stats[2]
	total_mass_sys = global_stats[3]
	n_atoms = global_stats[4]

	vol_global = volume * 1e-21 
	
	kinetic_energy = two_kinetic_energies / 2.0
	energy_total = kinetic_energy + energy_potential
	
	temperature = two_kinetic_energies / (3 * Kb * n_atoms)
	
	pressure = (two_kinetic_energies - virial) / (3.0 * vol_global * N_a)
	
	density = m_a * total_mass_sys / vol_global

	return energy_potential, \
		energy_total, \
		two_kinetic_energies/2, \
		temperature, \
		virial, \
		pressure, \
		density

def averageValues(atoms, volume):
	volume *= 1e-21
	vMCx = np.mean(atoms.data["v"]["x"])
	vMCy = np.mean(atoms.data["v"]["y"])
	vMCz = np.mean(atoms.data["v"]["z"])
	V2 = (atoms.data["v"]["x"]-vMCx)**2 + \
		 (atoms.data["v"]["y"]-vMCy)**2 + \
		 (atoms.data["v"]["z"]-vMCz)**2
	energy_potential = np.sum(atoms.data["Ep"])
	try:
		energy_potential += np.sum(atoms.data["Em"])
	except ValueError:
		pass
	two_kinetic_energies = np.sum(V2*atoms.data["m"])
	energy_total = two_kinetic_energies/2 + energy_potential
	temperature = two_kinetic_energies/Kb/3/atoms.data.size
	virial = np.sum(atoms.data["V"]["xx"] + \
		          atoms.data["V"]["yy"] + \
		          atoms.data["V"]["zz"])		
	pressure = (two_kinetic_energies - virial)/3.0/volume/N_a
	density = m_a*np.sum(atoms.data["m"])/volume
	return energy_potential, \
	       energy_total, \
	       two_kinetic_energies/2, \
		   temperature, \
		   virial, \
		   pressure, \
		   density


class Profile:

	def __init__(self, xmin, xmax, nbins, Ly, Lz):
		self.nbins = nbins
		self.bin_edges = np.linspace(xmin, xmax, nbins + 1)
		self.dx   = self.bin_edges[1] - self.bin_edges[0]
		self.bins = self.bin_edges[:-1] + 0.5*self.dx
		self.Vbin = self.dx*Ly*Lz*1e-21 # cm3

	def analyze(self, atoms, smooth = -1):
		data = np.zeros(self.nbins, dtype=[
			('concentration', 'f8'), # 1/cm3
			('density', 'f8'),       # g/cm3
			('temperature', 'f8'),   # K
			('pressure', 'f8'),      # kBar
			('velocity', 'f8')       # km/s
		])
		# evaluate concentration
		counts, tmp = np.histogram(
			atoms.data["r"]["x"], 
			bins    = self.bin_edges
		)
		data['concentration'] = counts/self.Vbin
		# evaluate density
		m, np.histogram(
			atoms.data["r"]["x"], 
			bins    = self.bin_edges,
			weights = atoms.data["m"]
		)
		data['density'] = m*mVbin
		# evaluate temperature and velocity
		binid = np.digitize(atoms.data["r"]["x"], bins = self.bin_edges, right=False)
		V2 = np.zeros(atoms.data.size, dtype=np.float64)
		for i in range(1, self.nbins + 1):
			idx = np.where(binid == i)
			# m calong x axis at bin i
			vMCx = 0.0
			if idx[0].size > 0:
				vMCx = np.sum(atoms.data["m"][idx]*atoms.data["v"]["x"][idx])
				vMCx /= np.sum(atoms.data["m"][idx])
			V2[idx] = (atoms.data["v"]["x"][idx] - vMCx)**2 + \
			           atoms.data["v"]["y"][idx]**2 + \
			           atoms.data["v"]["z"][idx]**2
			data['velocity'][i - 1] = vMCx

		T, tmp = np.histogram(
			atoms.data["r"]["x"], 
			bins    = self.bin_edges,
			weights = V2*atoms.data["m"]
		)
		ige0 = np.where(counts > 0)
		data['temperature'][ige0] = T[ige0]/(counts[ige0] * Kb * 3)

		virials = atoms.data["V"]["xx"] + \
		          atoms.data["V"]["yy"] + \
		          atoms.data["V"]["zz"]
		V, tmp = np.histogram(
			atoms.data["r"]["x"], 
			bins    = self.bin_edges,
			weights = virials
		)
		data['pressure'] = P_factor*(counts*Kb*T - V/3.0)/self.Vbin
		if smooth > 0:
			data['concentration'] = gaussian_filter1d(data['concentration'], sigma = smooth)
			data['density']       = gaussian_filter1d(data['density'],       sigma = smooth)
			data['temperature']   = gaussian_filter1d(data['temperature'],   sigma = smooth)
			data['pressure']      = gaussian_filter1d(data['pressure'],      sigma = smooth)
			data['velocity']      = gaussian_filter1d(data['velocity'],      sigma = smooth)
		return data
