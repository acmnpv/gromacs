#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2018 by the GROMACS development team.
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

# OpenCL required version: 1.2 or newer
set(REQUIRED_OPENCL_MIN_VERSION_MAJOR 1)
set(REQUIRED_OPENCL_MIN_VERSION_MINOR 2)
set(REQUIRED_OPENCL_MIN_VERSION ${REQUIRED_OPENCL_MIN_VERSION_MAJOR}.${REQUIRED_OPENCL_MIN_VERSION_MINOR})

set(GMX_GPU_OPENCL ON)

if(GMX_DOUBLE)
    message(FATAL_ERROR "OpenCL acceleration is not available in double precision")
endif()

# for some reason FindOpenCL checks CUDA_PATH but not CUDA_HOME
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH};$ENV{CUDA_HOME})
find_package(OpenCL)

if (OpenCL_FOUND)
    if (OpenCL_VERSION_STRING VERSION_LESS REQUIRED_OPENCL_MIN_VERSION)
        message(FATAL_ERROR "OpenCL " "${OpenCL_VERSION_STRING}" " is not supported. OpenCL version " "${REQUIRED_OPENCL_MIN_VERSION}" " or newer is required.")
        return ()
    endif()
else ()
    message(FATAL_ERROR "OpenCL not found.")
    return()
endif()

# Tell compiler to hide warnings for comments caused by cl_gl_ext.h on Linux
if (UNIX)
    set(OpenCL_DEFINITIONS ${OpenCL_DEFINITIONS} " -Wno-comment")
endif()

# Yes Virginia, Darwin kernel version 14.4 corresponds to OS X 10.4.
if(APPLE AND ${CMAKE_SYSTEM_VERSION} VERSION_LESS "14.4")
        message(WARNING "OS X prior to version 10.10.4 produces incorrect AMD OpenCL code at runtime. You will not be able to use AMD GPUs on this host unless you upgrade your operating system.")
endif()

add_definitions(${OpenCL_DEFINITIONS})

include_directories(SYSTEM ${OpenCL_INCLUDE_DIRS})

# Ensure the OpenCL implementation is 64-bit, because we only support that;
# Note that this can only be revised if the cl_mem size assumptions made
# (originally in pme-gpu-types.h) are relieved.
if (NOT CMAKE_SIZEOF_VOID_P EQUAL 8)
    message(FATAL_ERROR "The OpenCL implementation is only supported on 64-bit platforms.")
endif()

set(GMX_INSTALL_OCLDIR       ${GMX_INSTALL_GMXDATADIR}/opencl)
