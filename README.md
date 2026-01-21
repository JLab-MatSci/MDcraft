# Building and Installing MDcraft

## Cmake install

Load the install files for **CMake** of the needed version from the official web-site [https://cmake.org/download](https://cmake.org/download). **CMake**  can be built from the source code, but better use the ready-to-use binaries for your specific platform.

### macOS

Load the file [cmake-x.y.z-Darwin-x86_64.dmg](https://cmake.org/download) where x.y.z is the actual version number. After opening, the standard installation window (in macOS) must emerge.

![CMake Install macOS](./images/CMake_install_macOS.png)

Move **CMake** into Applications. To make **CMake** available from the command line interface (**cmake**, **сcmake**, **ctest** and **cpack** commands) write `/usr/local/bin`:
```bash
sudo "/Applications/CMake.app/Contents/bin/cmake-gui" --install
```
where the `/usr/local/bin` path is the default one for the binaries, but if you need another option, do:
```bash
sudo "/Applications/CMake.app/Contents/bin/cmake-gui" --install=/your/cmake/path/bin
```
This case, do not forget to add `/your/cmake/path/bin` into the environment **PATH** variable:
```bash
export PATH=/your/cmake/path/bin:${PATH}
```

### Linux

Load the file [cmake-x.y.z-Linux-x86_64.sh](https://cmake.org/download) where x.y.z is the actual version number. Go into the folder, where the file was loaded and change the execution rights accordingly:
```bash
chmod +x cmake-x.y.z-Linux-x86_64.sh
```
Next, run the script, set the installation path **prefix** (the default is `/usr/local`), the super-user rights may be needed (sudo):
```bash
sudo ./cmake-x.y.z-Linux-x86_64.sh --prefix=/your/cmake/path
```
Accept the license *y*:
```
Do you accept the license? [yN]: y
```
and install into the given path **prefix** without any additional subdirectory:
```
Do you want to include the subdirectory cmake-x.y.z-Linux-x86_64?
Saying no will install in: "/your/cmake/path" [Yn]: n
```
The build system will store all the needed files into `/your/cmake/path/bin`, `/your/cmake/path/doc`, `/your/cmake/path/man`, `/your/cmake/path/share`. Do not forget to add `/your/cmake/path/bin` into the environment **PATH** variable:
```bash
export PATH=/your/cmake/path/bin:${PATH}
```

## Eigen library

Load the source code of **Eigen** from the official web-site [https://eigen.tuxfamily.org](https://eigen.tuxfamily.org) or from GitLab [https://gitlab.com/libeigen/eigen](https://gitlab.com/libeigen/eigen). **Eigen** is a С++ header-only library. But we higly recommend to install it via **CMake**.

For now suppose the `eigen-3.4.0.tar.gz` version is loaded. Extract the files:
```bash
tar -zxf eigen-3.4.0.tar.gz
```
Make an additional folder to build the library and go into it:
```bash
mkdir eigen-build
cd eigen-build
```

Run the configuration **CMake** step from the command line:
```bash
ccmake ../eigen-3.4.0
```
Run the initial configuration step with pressing the [c] button. You must see the following message:
```
Configured Eigen 3.4.0

 Available targets (use: make TARGET):
 ---------+--------------------------------------------------------------
 Target   |   Description
 ---------+--------------------------------------------------------------
 install  | Install Eigen. Headers will be installed to:
          |     <CMAKE_INSTALL_PREFIX>/<INCLUDE_INSTALL_DIR>
          |   Using the following values:
          |     CMAKE_INSTALL_PREFIX: /usr/local
          |     INCLUDE_INSTALL_DIR:  include/eigen3
          |   Change the install location of Eigen headers using:
          |     cmake . -DCMAKE_INSTALL_PREFIX=yourprefix
          |   Or:
          |     cmake . -DINCLUDE_INSTALL_DIR=yourdir
 doc      | Generate the API documentation, requires Doxygen & LaTeX
 check    | Build and run the unit-tests. Read this page:
          |   http://eigen.tuxfamily.org/index.php?title=Tests
 blas     | Build BLAS library (not the same thing as Eigen)
 uninstall| Remove files installed by the install target
 ---------+--------------------------------------------------------------

 Configuring done
```
Next, press [e]. Now choose the following options:
```
BUILD_TESTING                   OFF
CMAKEPACKAGE_INSTALL_DIR        share/cmake/eigen/
CMAKE_BUILD_TYPE                Release
CMAKE_INSTALL_PREFIX            /opt/eigen # an example
EIGEN_BUILD_DOC                 OFF
```
Again, do the configuration steps with pressing the [c] button. Finally, run the generation step with pressig the [g] button. Now run
```bash
sudo make install
```

## Python installation

### Standard python3 and NumPy installation

We work with python 3.x, so install the corresponding packages. To start, make sure that **python3** with all the necessary header files is installed
```bash
sudo apt install python3-dev
```
then install **NumPy**, **SciPy** and **Matplotlib**:
```bash
sudo apt install python3-numpy python3-scipy python3-matplotlib
```
Now see the installed python3 version:
```bash
python3 --version
# Python 3.x.y
```

### Anaconda installation

If it is impossible to install the interpreter using a package manager, it is recommended to use the Anaconda distribution https://www.anaconda.com/products/individual Linux 64-Bit (x86) Installer. This package contains the Python interpreter and all necessary libraries.

## pybind11 installation

To build python modules from the C++ source code we used **pybind11** library: 
```bash
git clone https://github.com/pybind/pybind11
cd pybind11
git checkout v2.10.4 # переключаемся в стабильную версию, например
```
Create some working folder to build this library, enter it and run **CMake**:
```bash
mkdir pybind11-build
cd pybind11-build
ccmake ../pybind11
```
Follow the instructions described above for eigen library and set the following cache variables as:
```
 CMAKE_INSTALL_PREFIX             /path/to/pybind11/install # например: /opt/pybind11
 CMAKE_BUILD_TYPE                 Release
 PYBIND11_TEST                    OFF
 PYBIND11_CPP_STANDARD            -std=c++11
```
You must see one more variable to set the same python version as you installed above:
```
 PYBIND11_PYTHON_VERSION          3.x
```
After the configuration step, again run the generation with [g] button and then run
```bash
$ sudo make install
``` 
Check whether the folders *include* and *share* are correctly created in the *CMAKE_INSTALL_PREFIX*.

## Build and install MDcraft

Create the folder to build **MDcraft**, and launch **CMake**:
```bash
mkdir mdcraft-build
cd mdcraft-build
ccmake /path/to/src/mdcraft
```
Again do the similar configuration steps with pressing [c] button. The following options are presented:
```
 CMAKE_BUILD_TYPE                Release # choose the build type
 CMAKE_CXX_STANDARD              17      # choose the C++ version (11, 14 or 17, better use 17)
 CMAKE_INSTALL_PREFIX            /home/username/soft/mdcraft # or /Users/username/soft/mdcraft in mac
 Eigen3_DIR                      /opt/eigen/share/cmake/eigen3 # from Eigen install dir
 pybind11_DIR                    /opt/pybind11/share/cmake/pybind11 # from Pybind11 install dir
 PYBIND11_PYTHON_VERSION         3.x  # the actual Python version
 mdcraft_ENABLE_BPNN             ON   # enable BPNN potentials
 mdcraft_ENABLE_DeePMD           ON   # enable DeePMD potentials
 mdcraft_ENABLE_MLIP4            ON   # enable MTP potentials
 mdcraft_ENABLE_MPI              ON   # MPI support
 mdcraft_ENABLE_PROBLEMS         ON   # install some problems
 mdcraft_ENABLE_PYTHON           ON   # Python-interface for mdcraft
 mdcraft_ENABLE_TESTS            ON   # install tests
 mdcraft_THREADS_TYPE            STD  # choose the type of threads (either from std C++ lib or TBB lib)
```
After configuration and generation steps build and install the project:

```bash
make -j4 install
```

Check the tests which must appear in the folder for test, e.g.:
```bash
cd /home/username/soft/mdcraft/test/solver
```

Set the dynamic libraries paths in linux:
```bash
export LD_LIBRARY_PATH=/home/username/soft/mdcraft/lib:${LD_LIBRARY_PATH}
```
or in mac:
```bash
export DYLD_LIBRARY_PATH=/Users/username/soft/mdcraft/lib:${DYLD_LIBRARY_PATH}
```

Set the path to the python modules:
```bash
export PYTHONPATH=/home/username/soft/mdcraft/python:${PYTHONPATH}
```
or in mac:
```bash
export PYTHONPATH=/Users/username/soft/mdcraft/lib:${PYTHONPATH}
```

Now run some test script
```bash
python test-LJs.py
```
