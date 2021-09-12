#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2016 by the GROMACS development team.
# Copyright (c) 2019,2020,2021, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.

if (GMX_MPI)
    if (GMX_THREAD_MPI)
        message(STATUS "MPI is not compatible with thread-MPI. Disabling thread-MPI.")
        set(GMX_THREAD_MPI OFF CACHE BOOL
            "Build a thread-MPI-based multithreaded version of GROMACS (not compatible with MPI)" FORCE)
    endif ()
    set(GMX_LIB_MPI 1)
else ()
    set(GMX_LIB_MPI 0)
endif ()

# Manage the MPI setup.
# Note that we may want to execute tests or Python with MPI,
# even if we are not using an MPI-enabled GROMACS build.
#TODO(#3672): find_package(MPI COMPONENTS ...)
find_package(MPI)
if (GMX_LIB_MPI)
    if (NOT MPI_CXX_FOUND)
        message(FATAL_ERROR
                "MPI support requested, but no suitable MPI compiler found. Either set the "
                "MPI_CXX_COMPILER to the MPI compiler wrapper (often called mpicxx or mpic++), "
                "set CMAKE_CXX_COMPILER to a default-MPI-enabled compiler, "
                "or set the variables reported missing for MPI_CXX above.")
    elseif (MPI_CXX_VERSION VERSION_LESS 2.0)
        message(FATAL_ERROR "MPI version 2.0 or higher is required. Please update your MPI library.")
    endif ()
    # TODO(#3672): Replace usage of `MPI_LINKER_FLAGS` with MPI_LINK_FLAGS,
    #   MPI_<lang>_LINK_FLAGS, or target_link_libraries(... MPI::MPI_CXX)
    set(MPI_LINKER_FLAGS ${MPI_CXX_LINK_FLAGS})
    separate_arguments(MPI_CXX_LINK_FLAGS)
    #TODO(#3672, #3776): These should be acquired through the MPI::MPI_CXX target.
    include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})
    list(APPEND GMX_COMMON_LIBRARIES ${MPI_CXX_LIBRARIES})
endif ()

if (GMX_LIB_MPI)
    # Test for and warn about unsuitable MPI versions
    # Find path of the mpi compilers
    if (${MPI_CXX_FOUND})
        get_filename_component(_mpi_c_compiler_path "${MPI_CXX_COMPILER}" PATH)
        get_filename_component(_mpiexec_path "${MPIEXEC_EXECUTABLE}" PATH)
    else ()
        get_filename_component(_cmake_c_compiler_path "${CMAKE_C_COMPILER}" PATH)
        get_filename_component(_cmake_cxx_compiler_path "${CMAKE_CXX_COMPILER}" PATH)
    endif ()

    # Execute the ompi_info binary with the full path of the compiler wrapper
    # found, otherwise we run the risk of false positives.
    find_file(MPI_INFO_BIN ompi_info
              HINTS ${_mpi_c_compiler_path} ${_mpiexec_path}
              ${_cmake_c_compiler_path} ${_cmake_cxx_compiler_path}
              NO_DEFAULT_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)
    if (MPI_INFO_BIN)
        exec_program(${MPI_INFO_BIN}
                     ARGS -v ompi full
                     OUTPUT_VARIABLE OPENMPI_TYPE
                     RETURN_VALUE OPENMPI_EXEC_RETURN)
        if (OPENMPI_EXEC_RETURN EQUAL 0)
            string(REGEX REPLACE ".*Open MPI: \([0-9]+\\.[0-9]*\\.?[0-9]*\).*" "\\1" OPENMPI_VERSION ${OPENMPI_TYPE})
            if (OPENMPI_VERSION VERSION_LESS "1.4.1")
                MESSAGE(WARNING
                        "CMake found OpenMPI version ${OPENMPI_VERSION} on your system. "
                        "There are known problems with GROMACS and OpenMPI version < 1.4.1. "
                        "Please consider updating your OpenMPI if your MPI wrapper compilers "
                        "are using the above OpenMPI version.")
            endif ()
            if (OPENMPI_VERSION VERSION_EQUAL "1.8.6")
                MESSAGE(WARNING
                        "CMake found OpenMPI version ${OPENMPI_VERSION} on your system. "
                        "This OpenMPI version is known to leak memory with GROMACS,"
                        "please update to a more recent version. ")
            endif ()
            unset(OPENMPI_VERSION)
            unset(OPENMPI_TYPE)
            unset(OPENMPI_EXEC_RETURN)
        endif ()
    endif ()
    unset(MPI_INFO_BIN CACHE)

    # Execute the mpiname binary with the full path of the compiler wrapper
    # found, otherwise we run the risk of false positives.
    find_file(MPINAME_BIN mpiname
              HINTS ${_mpi_c_compiler_path}
              ${_cmake_c_compiler_path} ${_cmake_cxx_compiler_path}
              NO_DEFAULT_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)
    if (MPINAME_BIN)
        exec_program(${MPINAME_BIN}
                     ARGS -n -v
                     OUTPUT_VARIABLE MVAPICH2_TYPE
                     RETURN_VALUE MVAPICH2_EXEC_RETURN)
        if (MVAPICH2_EXEC_RETURN EQUAL 0)
            string(REGEX MATCH "MVAPICH2" MVAPICH2_NAME ${MVAPICH2_TYPE})
            # Want to check for MVAPICH2 in case some other library supplies mpiname
            string(REGEX REPLACE "MVAPICH2 \([0-9]+\\.[0-9]*[a-z]?\\.?[0-9]*\)" "\\1" MVAPICH2_VERSION ${MVAPICH2_TYPE})
            if (${MVAPICH2_NAME} STREQUAL "MVAPICH2" AND MVAPICH2_VERSION VERSION_LESS "1.5")
                # This test works correctly even with 1.5a1
                MESSAGE(WARNING
                        "CMake found MVAPICH2 version ${MVAPICH2_VERSION} on your system. "
                        "There are known problems with GROMACS and MVAPICH2 version < 1.5. "
                        "Please consider updating your MVAPICH2 if your MPI wrapper compilers "
                        "are using the above MVAPICH2 version.")
            endif ()
            unset(MVAPICH2_VERSION)
            unset(MVAPICH2_NAME)
            unset(MVAPICH2_TYPE)
            unset(MVAPICH2_EXEC_RETURN)
        endif ()
    endif ()
    unset(MPINAME_BIN CACHE)
endif () # end if(GMX_LIB_MPI)

# Look for MPI process launchers that may be missed, especially if we didn't
# need to find the full MPI library build system support.
if (NOT MPIEXEC_EXECUTABLE)
    find_program(MPIEXEC
                 NAMES mpiexec mpirun lamexec srun aprun poe
                 HINTS ${MPI_HOME} $ENV{MPI_HOME}
                 PATH_SUFFIXES bin
                 DOC "Executable for running MPI programs.")

    set(MPIEXEC_EXECUTABLE "$MPIEXEC")
    set(MPIEXEC_NUMPROC_FLAG "-np" CACHE STRING "Flag used by MPI to specify the number of processes for MPIEXEC; the next option will be the number of processes.")
    set(MPIEXEC_PREFLAGS "" CACHE STRING "These flags will be directly before the executable that is being run by MPIEXEC.")
    set(MPIEXEC_POSTFLAGS "" CACHE STRING "These flags will come after all flags given to MPIEXEC.")
    set(MPIEXEC_MAX_NUMPROCS "2" CACHE STRING "Maximum number of processors available to run MPI applications.")
    mark_as_advanced(MPIEXEC MPIEXEC_NUMPROC_FLAG MPIEXEC_PREFLAGS MPIEXEC_POSTFLAGS MPIEXEC_MAX_NUMPROCS)
endif ()
