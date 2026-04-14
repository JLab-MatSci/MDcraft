import numpy as np
import math, sys, os, time as Time
import h5py as hdf
from mdcraft.tools import Threads
threads = Threads(4)
# data 
from mdcraft.data import Atoms
# output
from mdcraft.neibs import VerletList
# geometry
from mdcraft.lattice import Domain
# solver
from mdcraft.solver.potential  import EAM
from mdcraft.solver.potential  import Pair
from mdcraft.solver.thermostat import Langevin
from mdcraft.solver.boundary   import Periodic
from mdcraft.solver            import Multi as Solver
from mdcraft.solver.stepper    import Verlet as Stepper

kBuf	 = 1.25           # buffer size
NCbuff   = 10             # nlist update steps
Al_mass  = 26.9815384     # g/mol
Ni_mass  = 58.6934        # g/mol
Al_part  = 0.8            # % of aluminum atoms
Ni_part  = 1.0 - Al_part  # % of nickel atoms
T0	     = 300.           # K
Kb	     = 0.008314462175 # kJ/mol/K
N_a      = 6.02214076e23  # avogadro

NCbuff   = 10           # nlist update steps
NCbuffMW = 5*NCbuff     # nlist update
NCout    = 100          # output profiles
tend	 = 12.0         # final time (ps)
dt	     = 4e-3         # timestep (ps)
NCsave   = int(tend/dt) # save atoms
# load potential
Al_potential = EAM('Al.eam.alloy')
Ni_potential = EAM('Al.eam.alloy')
NiAl_potential = Pair('NiAlpair.eam.alloy')
# setup atoms mesh
from mesh import Mesh
mesh = Mesh(
	Lx   = 8.,     # nm
	Ly   = 8.,     # nm
	Lz   = 8.,     # nm
	D    = 0.405,  # nm, from potential file
	kind = 'fcc'
)
atoms = mesh.setup()

rndindex = np.random.choice([0, 1], atoms.data.size, p=[Al_part, Ni_part])
iAl = np.where(rndindex == 0)
iNi = np.where(rndindex == 1)

atoms.data["m"][iAl] = Al_mass
atoms.data["m"][iNi] = Ni_mass
atoms.data["rcut"][iAl] = Al_potential.rcut 
atoms.data["rcut"][iNi] = Ni_potential.rcut 

# initial velocity distribution
mu = 0.0

sigma = (1.2*Kb*T0/Al_mass)**0.5
Vx = np.random.normal(mu, sigma, iAl[0].size)
Vy = np.random.normal(mu, sigma, iAl[0].size)
Vz = np.random.normal(mu, sigma, iAl[0].size)
atoms.data["v"]["x"][iAl] = Vx
atoms.data["v"]["y"][iAl] = Vy
atoms.data["v"]["z"][iAl] = Vz

sigma = (1.2*Kb*T0/Ni_mass)**0.5
Vx = np.random.normal(mu, sigma, iNi[0].size)
Vy = np.random.normal(mu, sigma, iNi[0].size)
Vz = np.random.normal(mu, sigma, iNi[0].size)
atoms.data["v"]["x"][iNi] = Vx
atoms.data["v"]["y"][iNi] = Vy
atoms.data["v"]["z"][iNi] = Vz

V_comx = np.sum(atoms.data["m"]*atoms.data["v"]["x"])/np.sum(atoms.data["m"])
V_comy = np.sum(atoms.data["m"]*atoms.data["v"]["y"])/np.sum(atoms.data["m"])
V_comz = np.sum(atoms.data["m"]*atoms.data["v"]["z"])/np.sum(atoms.data["m"])
atoms.data["v"]["x"] -= V_comx
atoms.data["v"]["y"] -= V_comy
atoms.data["v"]["z"] -= V_comz

# setup boundary conditions
Lx, Ly, Lz = mesh.box_sizes()
volume = Lx * Ly * Lz
Bxmin, Bymin, Bzmin, \
Bxmax, Bymax, Bzmax = mesh.sample_sizes()
domain = Domain(
	Bxmin, Bxmax,
	Bymin, Bymax,
	Bzmin, Bzmax
)
domain.set_periodic(0, True)
domain.set_periodic(1, True)
domain.set_periodic(2, True)
periodic = Periodic(domain)

# setup thermostat
thermostat = Langevin(
	beta		= 0.5, # ps^-1 
	temperature = T0,  # K 
	time_step   = dt,  # ps
	heat_x	    = 1,   # x axis heating
	heat_y	    = 1,   # y axis heating
	heat_z	    = 1,   # z axis heating
	Ux	        = 0,   # flow velocity x
	Uy	        = 0,   # flow velocity y
	Uz	        = 0    # flow velocity z
)
# solver with thermostat
solver = Solver(
	boundary   = periodic, 
	# thermostat = thermostat, # uncomment to use
	nkinds     = 2,
	threads    = threads
)
solver.add_potential(0, 0, Al_potential)
solver.add_potential(1, 1, Ni_potential)
solver.add_potential(0, 1, NiAl_potential)
# solver without thermostat
solver_cold = Solver(
	boundary   = periodic, 
	nkinds     = 2,
	threads    = threads
)
solver_cold.add_potential(0, 0, Al_potential)
solver_cold.add_potential(1, 1, Ni_potential)
solver_cold.add_potential(0, 1, NiAl_potential)
# setup MD stepper
stepper = Stepper(
	solver  = solver,
	threads = threads
)
# setup neighbors list
atoms.data["rns"] = kBuf*atoms.data["rcut"]
nlist = VerletList(
	atoms     = atoms,
	neighbors = atoms,
	threads   = threads,
	domain    = domain
)
nlist.update()

# init profile for analysis
import analysis

nbins = 1
profile = analysis.Profile(
	xmin = Bxmin, xmax = Bxmax, 
	Ly   = Bymax - Bymin, 
	Lz   = Bzmax - Bzmin, 
	nbins = nbins
)

# prepare output
if not os.path.exists("./log/"):
	os.makedirs("./log/")
if not os.path.exists("./atoms/"):
	os.makedirs("./atoms/")

logfile = open("./log/log-prepare.dat","w")
logfile.write(
	"%06s %15s %15s %15s %15s %15s %15s %15s %15s\n" % (
		"step", "time", "potentialEnergy", "totalEnergy",
		"kineticEnergy", "temperature", "virial", "pressure", "density"
	)
)

# prepare main cycle
ICbuff   = 0
ICout    = 0
ICsave   = 0
ICwait   = 0
ICskip   = 0

time	 = 0.0
step	 = 0

start = Time.time()

# Calculate initial forces and energy 
solver_cold.prepare(atoms, atoms, nlist.get())
solver_cold.forces(atoms,  atoms, nlist.get())
solver_cold.virials(atoms, atoms, nlist.get())

# start main cycle
while step < NCsave:

	stepper.make_step(atoms, nlist.get(), dt)

	time += dt
	step += 1

	ICbuff = step % NCbuff
	ICout  = step % NCout
	ICsave = step % NCsave

	if ICout  == 0:
		solver.virials(atoms, atoms, nlist.get())
		energy_potential, \
		energy_total, \
		kinetic_energy, \
		temperature, \
		virial, \
		pressure, \
		density = analysis.averageValues(atoms, volume)
						
		elapsed = Time.time() - start
		logfile.write("%06d %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e\n" %(
			step, time, 
			energy_potential, 
			energy_total, 
			kinetic_energy, 
			temperature, 
			virial, 
			pressure, 
			density
		))
		print("step: %06d, phys_time: %6.2f ps, sys_time: %6.2f s" % (step, time, elapsed))
		logfile.flush()

	if ICsave == 0:
		save_file = hdf.File("atoms/prep%06d.h5" % step, 'w')
		save_file["atoms"] = atoms.data
		save_file["domain"] = np.array([Bxmin, Bxmax, Bymin, Bymax, Bzmin, Bzmax])
		save_file.close()
		
	if ICbuff == 0:
		nlist.update()

logfile.close()
