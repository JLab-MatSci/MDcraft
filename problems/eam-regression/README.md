# NiAl Multi / Langevin regression

`test-NiAl.py` is a small deterministic FCC Al/Ni setup using the public potential
tables in ../NiAl. It follows the C/H factory, cold initialization, full
neighbor lists, VD3 decomposition and Verlet workflow without private inputs.
The default is 256 atoms and 20 steps; this is a diagnostic, not an equilibrated
production simulation or a force-accuracy reference.

After rebuilding and installing the Python bindings into the selected environment:

```bash
cd problems/eam-regression
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

## All analytic EAM variants

`test-analytic.py` constructs EAMW1046, EAMW0644, EAMWMOD1046,
EAMWMOD0644, EAMMOD1044 and EAMV1622 from small **synthetic public** coefficient
arrays. They are implementation fixtures, not calibrated C/H parameters.
The serial run checks forces against central differences of total energy,
repeated virial/energy evaluation, array validation and short MD trajectories.
MPI runs use the same coordinates and coefficients with full local/halo lists.

```bash
python test-analytic.py --serial --threads 2
mpiexec -n 1 python test-analytic.py --threads 2
mpiexec -n 2 python test-analytic.py --threads 2
```

The programs print JSON energies and per-UID coordinates/forces so serial and
MPI results can be compared. Thermostat-free synthetic fixtures should agree
across rank and thread counts within floating-point roundoff. Private C/H
coefficients are not included, and these checks do not validate their calibration.

C++ CTest regressions cover periodic density counting, table endpoints, invalid
input handling and once-per-step MPI thermostat application. `ctest` also runs
the analytic Python test and the NiAl lifetime regression after installation.


## User-supplied C/H inputs

`main-ch.py` is the corrected MPI C/H driver. It retains the thermostat,
uses separate hot/cold boundary objects, computes the initial densities, and
prints actual worker counts. It shares the generic `mpi_output.py` helper.
It requires external `eamReader.py`, `ase2mdcraft.py`, LAMMPS initial data and
coefficient files (`inputs/EAMW1046.dat`, `EAMW0644.dat`, `EAMW1044.dat`).
These private inputs are not included and this C/H production driver is not
part of the self-contained regression suite.

```bash
mpiexec -n 2 python main-ch.py --threads 2 --steps 20 --data /path/to/system.lmp --inputs /path/to/inputs
```

## Validation of this integration

Both MPI-disabled and MPI-enabled builds were tested with AddressSanitizer and
UndefinedBehaviorSanitizer: 30 and 31 CTest cases passed, respectively. Leak
checking was disabled for the sanitizer runs. Eighteen separate Python runs
covered both scripts at requested thread counts 1, 2 and 4, in serial and with
1/2 MPI ranks; analytic energies and per-UID states agreed with the serial
reference (rtol 1e-10, atol 1e-12). The two-rank run used 3 actual threads when
4 were requested, as reported by the thread pool. A standalone ThreadSanitizer
stress test of the shared Langevin RNG also passed.

These results validate implementation consistency with the supplied public
and synthetic fixtures; they do not validate a private C/H parameter fit.

## Diagnosing cluster heap corruption

`malloc(): invalid size (unsorted)` at `Multi::add_potential` can be a delayed
report of an earlier invalid write. Do not treat the vector allocation itself
as the established source of corruption. First run without private data:

```bash
mpiexec -n 1 python diagnose-runtime.py --threads 1
mpiexec -n 2 python diagnose-runtime.py --threads 2
```

The script prints the interpreter, actual extension/shared-library paths, Atom
dtype offsets and MPI implementation, then creates 1000 Multi factories using
synthetic coefficients. It does not read a data file or create a decomposition.
This is a localization aid, not a complete ABI-compatibility proof. In this
branch's tested build Atom has itemsize 224 and nc_mod has two double entries;
platform/compiler differences must be evaluated against matching native builds.

Check that every rank imports the intended package and no older libmdcraft is
loaded from another prefix through LD_LIBRARY_PATH. After changes to Atom,
Pair/ISolver, potential classes, compiler flags or MPI, use a new build directory
and reinstall all native libraries AND Python extensions together. The MPI used
by mpi4py must match the one used to build MDcraft.

If factory stress passes but main-ch.py fails, capture the full native backtrace
and determine whether the first failure occurs before or after data conversion,
VD3 initialization, initial forces or the first integration step. For a debug
build, use AddressSanitizer/UndefinedBehaviorSanitizer (including Python extension
modules), or Valgrind to locate the first invalid write. Fixing the final malloc
failure site without that evidence can only hide the underlying error.
