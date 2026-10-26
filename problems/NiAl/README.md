# NiAl Multi / Langevin regression

`test-NiAl.py` is a small deterministic FCC Al/Ni setup using the public potential
tables in this directory. It follows the C/H factory, cold initialization, full
neighbor lists, VD3 decomposition and Verlet workflow without private inputs.
The default is 256 atoms and 20 steps; this is a diagnostic, not an equilibrated
production simulation or a force-accuracy reference.

After rebuilding and installing the Python bindings into the selected environment:

```bash
cd problems/NiAl
# Requires an MPI-disabled MDcraft build:
python test-NiAl.py --serial --threads 2 --output /tmp/nial-serial
# Requires an MPI-enabled build and matching mpi4py runtime:
mpiexec -n 1 python test-NiAl.py --threads 2 --output /tmp/nial-mpi1
mpiexec -n 2 python test-NiAl.py --threads 2 --output /tmp/nial-mpi2
```

Add `--lifetime-regression` to exercise a factory that drops its local Python
thermostat reference. The test asserts that `Multi` keeps the thermostat alive
before executing forces. This requires the `keep_alive` policies in `_solver.cxx`.
Ordinary runs additionally retain the thermostat explicitly in Python.

The test validates finite atom state and conserved global atom count. MPI runs
also calculate thermodynamics and write a gathered LAMMPS snapshot using
`mpi_output.py`. Actual worker counts are printed because `Threads` may clamp
the requested count. Use distinct output directories for simultaneous runs.

When `Atom` changes, rebuild all C++ libraries and Python modules together.
Do not combine modules compiled against different layouts of `Atom`.
