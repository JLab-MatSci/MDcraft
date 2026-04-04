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
from mdcraft.solver.thermostat import Langevin
from mdcraft.solver.boundary   import Periodic
from mdcraft.solver            import Single as Solver
from mdcraft.solver.stepper    import Verlet as Stepper

kBuf	 = 1.25           # buffer size
NCbuff   = 10             # nlist update steps
Al_mass  = 26.9815384     # g/mol
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
potential = EAM('Al2009.eam.alloy')
# setup atoms mesh
from mesh import Mesh
mesh = Mesh(
	Lx   = 8.,     # nm
	Ly   = 8.,     # nm
	Lz   = 8.,     # nm
	D    = 0.4032, # nm, from potential file
	kind = 'fcc'
)
atoms = mesh.setup()
atoms.data["rcut"] = potential.rcut 
atoms.data["m"] = Al_mass

# initial velocity distribution
mu = 0.0
sigma = (1.2*Kb*T0/Al_mass)**0.5
Vx = np.random.normal(mu, sigma, atoms.data.size)
Vy = np.random.normal(mu, sigma, atoms.data.size)
Vz = np.random.normal(mu, sigma, atoms.data.size)
V_comx = np.sum(Vx)/atoms.data.size
V_comy = np.sum(Vy)/atoms.data.size
V_comz = np.sum(Vz)/atoms.data.size
atoms.data["v"]["x"] = Vx - V_comx
atoms.data["v"]["y"] = Vy - V_comy
atoms.data["v"]["z"] = Vz - V_comz

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
	potential  = potential, 
	boundary   = periodic, 
	# thermostat = thermostat, # uncomment to use
	threads    = threads
)
# solver without thermostat
solver_cold = Solver(
	potential  = potential, 
	boundary   = periodic, 
	threads    = threads
)
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
import analysis_serial

nbins = 1
profile = analysis_serial.Profile(
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
		density = analysis_serial.averageValues(atoms, volume)
						
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
