if(NOT Python_EXECUTABLE)
    find_package(Python 3 REQUIRED COMPONENTS Interpreter)
endif()

if(Python_EXECUTABLE)
    execute_process(COMMAND "${Python_EXECUTABLE}" -c "import mpi4py; print(mpi4py.__path__[0])"
                    OUTPUT_VARIABLE MPI4PY_PATH
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    RESULT_VARIABLE _mpi4py_result)

    if(_mpi4py_result EQUAL 0)
        set(MPI4PY_FOUND TRUE CACHE BOOL "Whether mpi4py was found")
        set(MPI4PY_INCLUDE_DIR "${MPI4PY_PATH}/include" CACHE PATH "Path to mpi4py installation")
        # You might need to add logic to find MPI libraries if your project links against them
        # find_package(MPI REQUIRED)
        # set(MPI4PY_LIBRARIES ${MPI_LIBRARIES})
    else()
        set(MPI4PY_FOUND FALSE CACHE BOOL "Whether mpi4py was found")
        if(MPI4PY_FIND_REQUIRED)
            message(FATAL_ERROR
                "mpi4py was not found for Python interpreter '${Python_EXECUTABLE}'. "
                "For the automated Linux flow use ENABLE_MPI=ON ./scripts/bootstrap_linux.sh, "
                "or install it manually with ${Python_EXECUTABLE} -m pip install mpi4py.")
        endif()
    endif()
else()
    set(MPI4PY_FOUND FALSE CACHE BOOL "Whether mpi4py was found")
    if(MPI4PY_FIND_REQUIRED)
        message(FATAL_ERROR "Python interpreter was not resolved before searching for mpi4py.")
    endif()
endif()

mark_as_advanced(MPI4PY_INCLUDE_DIR)
