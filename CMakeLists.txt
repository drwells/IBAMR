## ---------------------------------------------------------------------
##
## Copyright (c) 2020 - 2020 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

#
# config info
#
CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)
PROJECT(IBAMR
  VERSION 1.0.0
  DESCRIPTION "Software infrastructure for the IB method with adaptively-refined grids"
  HOMEPAGE_URL "https://ibamr.github.io"
  # include C so that we can link against C libraries (e.g., MPI::MPI_C) easily
  LANGUAGES C CXX Fortran)
SET(IBAMR_DIMENSIONS "2" "3")

#
# -1. Hardcode some things that will ultimately be configuration options:
#
SET(SAMRAI_PREFIX "/home/drwells/Applications/samrai-2.4.4")
SET(PETSC_PREFIX "/home/drwells/Applications/petsc-3.13.0/x86_64-opt")
SET(LIBMESH_PREFIX "/home/drwells/Applications/libmesh-dev")

#
# 0. set the stage:
#
MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")
include(GNUInstallDirs)
INCLUDE(CMakePackageConfigHelpers)
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)

# we need this since SAMRAI is usually statically linked
SET(CMAKE_POSITION_INDEPENDENT_CODE ON)

#
# 1: manage mandatory dependencies:
#
FIND_PACKAGE(MPI REQUIRED)

# a couple of dependencies are also bundled
FIND_PACKAGE(Boost 1.57 QUIET)
IF(Boost_FOUND)
  SET(IBAMR_USE_BUNDLED_Boost FALSE)
  MESSAGE(STATUS "Found external boost ${Boost_VERSION} at ${Boost_DIR}")
ELSE()
  SET(IBAMR_USE_BUNDLED_Boost TRUE)
  MESSAGE(STATUS "Setting up boost as a bundled dependency")
  ADD_LIBRARY(BUNDLED_Boost INTERFACE)
  TARGET_INCLUDE_DIRECTORIES(
    BUNDLED_Boost
    INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/ibtk/contrib/boost>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/contrib/boost>)
  INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/ibtk/contrib/boost DESTINATION include/contrib)
ENDIF()

FIND_PACKAGE(Eigen3 3.2.5 QUIET)
IF(Eigen3_FOUND)
  SET(IBAMR_USE_BUNDLED_Eigen FALSE)
  MESSAGE(STATUS "Found external Eigen ${Eigen_VERSION} at ${Eigen_DIR}")
ELSE()
  SET(IBAMR_USE_BUNDLED_Eigen TRUE)
  MESSAGE(STATUS "Setting up Eigen as a bundled dependency")
  ADD_LIBRARY(BUNDLED_Eigen INTERFACE)
  TARGET_INCLUDE_DIRECTORIES(
    BUNDLED_Eigen
    INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/ibtk/contrib/eigen>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/contrib/eigen>)
  INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/ibtk/contrib/eigen/Eigen DESTINATION include/contrib)
  INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/ibtk/contrib/eigen/unsupported DESTINATION include/contrib)
ENDIF()

SET(HDF5_FIND_DEBUG TRUE)
FIND_PACKAGE(HDF5 REQUIRED COMPONENTS C HL)

#
# hypre:
#
MESSAGE(STATUS "")
MESSAGE(STATUS "Setting up HYPRE")
FIND_PATH(HYPRE_INCLUDE_DIRS NAMES HYPRE.h HINTS ${HYPRE_DIR}/include /usr/include/hypre)
MESSAGE(STATUS "HYPRE_INCLUDE_DIRS: ${HYPRE_INCLUDE_DIRS}")
FIND_LIBRARY(HYPRE_LIBRARIES NAMES HYPRE HINTS ${HYPRE_DIR}/lib)
MESSAGE(STATUS "HYPRE_LIBRARIES: ${HYPRE_LIBRARIES}")
# samrai 2.4.4 does not yet have a CMake module
SET(_samrai_library_suffixes "algs" "appu" "geom" "hier"
  "math_std" "mesh" "pdat_std" "solv" "xfer")

# SAMRAI consists of about a dozen libraries (see the above list) as well as
# header files for both 2D and 3D. Rather than setting up each as a target, set
# up our own 'header-only' target libraries that depend on SAMRAI libraries to
# encapsulate them neatly:
ADD_LIBRARY(SAMRAI2d INTERFACE IMPORTED)
ADD_LIBRARY(SAMRAI3d INTERFACE IMPORTED)

# there is also just libSAMRAI.a
ADD_LIBRARY(SAMRAI UNKNOWN IMPORTED)
FIND_LIBRARY(SAMRAI_path SAMRAI HINTS "${SAMRAI_PREFIX}/lib")
MESSAGE(STATUS "found samrai library ${SAMRAI_path}")
SET_TARGET_PROPERTIES(SAMRAI PROPERTIES IMPORTED_LOCATION ${SAMRAI_path})

TARGET_LINK_LIBRARIES(SAMRAI2d INTERFACE SAMRAI)
TARGET_LINK_LIBRARIES(SAMRAI3d INTERFACE SAMRAI)

FOREACH(_d ${IBAMR_DIMENSIONS})
  FOREACH(_suffix ${_samrai_library_suffixes})
    SET(_lib_name "SAMRAI${_d}d_${_suffix}")
    ADD_LIBRARY(
      ${_lib_name}
      UNKNOWN # library type (though its currently always static)
      IMPORTED)

    FIND_LIBRARY("${_lib_name}_path"
      ${_lib_name}
      HINTS
      "${SAMRAI_PREFIX}/lib"
      )
    MESSAGE(STATUS "found samrai library ${${_lib_name}_path}")

    SET_TARGET_PROPERTIES(
      ${_lib_name}
      PROPERTIES
      IMPORTED_LOCATION
      "${${_lib_name}_path}"
      )

    # add the library we just set up to the interface target:
    TARGET_LINK_LIBRARIES("SAMRAI${_d}d" INTERFACE ${_lib_name})
  ENDFOREACH()
ENDFOREACH()

#
# 1.25: manage mandatory dependencies that have pkgconfig information:
#
FIND_PACKAGE(PkgConfig REQUIRED)
PKG_CHECK_MODULES(MUPARSER REQUIRED IMPORTED_TARGET GLOBAL muparser)
LIST(PREPEND CMAKE_PREFIX_PATH "${PETSC_PREFIX}")
PKG_CHECK_MODULES(PETSC IMPORTED_TARGET GLOBAL PETSc)
LIST(POP_FRONT CMAKE_PREFIX_PATH)

#
# 1.5: manage optional dependencies
#
LIST(PREPEND CMAKE_PREFIX_PATH "${LIBMESH_PREFIX}")
PKG_CHECK_MODULES(LIBMESH IMPORTED_TARGET GLOBAL libmesh)
# todo - set up variables correctly
SET(IBAMR_WITH_LIBMESH TRUE)
LIST(POP_FRONT CMAKE_PREFIX_PATH)
#
# there are several issues with libmesh's pkgconfig file:
# 1. the link order is wrong, e.g., -L/usr/lib/openmpi comes after -lmpi (this will be fixed in 1.6)
# 2. cflags contains unnecessary flags (e.g., -felide-constructors) that should
#    be in Cflags_other (this will be fixed in 1.6)
# 3. it assumes -I${includedir}/contrib always exists (this will be fixed in 1.6)
#
# It may be better in the long run to write our own wrapper that avoids these
# problems.
#


#
# 2: IBAMR-specific configuration:
#

# Macro to process m4 files to generate Fortran
INCLUDE(CMakeParseArguments)
MACRO(IBAMR_PROCESS_M4)
  SET(options)
  SET(singleValueArgs NAME)
  SET(multiValueArgs)

  CMAKE_PARSE_ARGUMENTS(arg "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})
  GET_FILENAME_COMPONENT(_in_directory "${CMAKE_CURRENT_SOURCE_DIR}/${arg_NAME}.m4" DIRECTORY)
  GET_FILENAME_COMPONENT(_out_directory "${CMAKE_CURRENT_BINARY_DIR}/${arg_NAME}.f" DIRECTORY)
  FILE(MAKE_DIRECTORY ${_out_directory})

  SET(_args "-DSAMRAI_FORTDIR=${SAMRAI_PREFIX}/include;-DCURRENT_SRCDIR=${_in_directory};-DTOP_SRCDIR=${CMAKE_SOURCE_DIR}")
  SET(_input "${CMAKE_CURRENT_SOURCE_DIR}/${arg_NAME}.m4")
  SET(_output "${CMAKE_CURRENT_BINARY_DIR}/${arg_NAME}")

  # MESSAGE(STATUS "input:  ${_input}")
  # MESSAGE(STATUS "output: ${_output}")

  ADD_CUSTOM_COMMAND(
    OUTPUT ${_output}
    DEPENDS ${_input}
    COMMAND m4 ${_args} ${_input} > ${_output}
    VERBATIM)

  SET_SOURCE_FILES_PROPERTIES(${_output} PROPERTIES GENERATED true)
  SET_SOURCE_FILES_PROPERTIES(${_output} PROPERTIES LANGUAGE Fortran)
ENDMACRO()

# Macro to setup an IBAMR target library with all common features (i.e.,
# everything but the source files)
FUNCTION(IBAMR_SETUP_TARGET_LIBRARY target_library)
  MESSAGE(STATUS "setting up target ${target_library}")
  # we need C++11 (todo: libMesh assumes C++14 so this is temporarily wrong)
  SET_PROPERTY(TARGET ${target_library} PROPERTY CXX_STANDARD 14)
  # we and our users will use these MPI functions so make the interface public:
  TARGET_LINK_LIBRARIES(${target_library} PUBLIC MPI::MPI_C)
  IF(${IBAMR_WITH_LIBMESH})
    TARGET_LINK_LIBRARIES(${target_library} PUBLIC MPI::MPI_CXX)
  ENDIF()
  # same for HDF5:
  # todo: no nice target support yet for HDF5
  TARGET_LINK_LIBRARIES(${target_library} PUBLIC "${HDF5_LIBRARIES}")
  TARGET_LINK_LIBRARIES(${target_library} PUBLIC "${HDF5_HL_LIBRARIES}")
  TARGET_INCLUDE_DIRECTORIES(${target_library} PUBLIC "${HDF5_INCLUDE_DIRS}")
  # similar for hypre
  # todo: no nice target support yet for hypre
  TARGET_LINK_LIBRARIES(${target_library} PUBLIC "${HYPRE_LIBRARIES}")
  TARGET_INCLUDE_DIRECTORIES(${target_library} PUBLIC "${HYPRE_INCLUDE_DIRS}")
  # we only need boost headers
  IF(IBAMR_USE_BUNDLED_Boost)
    TARGET_LINK_LIBRARIES(${target_library} PUBLIC BUNDLED_Boost)
  ELSE()
    TARGET_LINK_LIBRARIES(${target_library} PUBLIC Boost::headers)
  ENDIF()
  # Eigen only has one target:
  IF(IBAMR_USE_BUNDLED_Eigen)
    TARGET_LINK_LIBRARIES(${target_library} PUBLIC BUNDLED_Eigen)
  ELSE()
    TARGET_LINK_LIBRARIES(${target_library} PUBLIC Eigen3::Eigen)
  ENDIF()
  # muparser:
  TARGET_LINK_LIBRARIES(${target_library} PUBLIC PkgConfig::MUPARSER)
  # PETSc:
  TARGET_LINK_LIBRARIES(${target_library} PUBLIC PkgConfig::PETSC)
  # libMesh:
  IF(IBAMR_WITH_LIBMESH)
    TARGET_LINK_LIBRARIES(${target_library} PUBLIC PkgConfig::LIBMESH)
  ENDIF()

  # Figure out if we are a 2D or a 3D library:
  SET(_2_location = -1)
  SET(_3_location = -1)
  STRING(FIND "${target_library}" "2" _2_location)
  STRING(FIND "${target_library}" "3" _3_location)

  # We cannot add SAMRAI to the link interface for our own libraries since
  # SAMRAI doesn't necessarily compile with -fPIC: we will instead put it in a
  # link interface for executables.
  TARGET_INCLUDE_DIRECTORIES(${target_library} PUBLIC "${SAMRAI_PREFIX}/include")

  # add the right preprocessor flag:
  IF(NOT ${_2_location} STREQUAL "-1")
    MESSAGE(STATUS "Adding flag -DNDIM=2 to target ${target_library}")
    TARGET_COMPILE_OPTIONS(${target_library} PUBLIC -DNDIM=2)
  ELSEIF(NOT ${_3_location} STREQUAL "-1")
    MESSAGE(STATUS "Adding flag -DNDIM=3 to target ${target_library}")
    TARGET_COMPILE_OPTIONS(${target_library} PUBLIC -DNDIM=3)
  ELSE()
    MESSAGE(STATUS "no NDIM preprocessor constant added to ${target_library}")
  ENDIF()
ENDFUNCTION()

#
# IBTK and IBAMR are compiled in 2D and 3D, but the headers are only installed
# once - hence set up the headers as separate targets
#
ADD_LIBRARY(IBTKHeaders INTERFACE)
TARGET_INCLUDE_DIRECTORIES(
  IBTKHeaders
  INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/ibtk/include>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/ibtk>)

ADD_LIBRARY(IBAMRHeaders INTERFACE)
TARGET_INCLUDE_DIRECTORIES(
  IBAMRHeaders
  INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/ibamr>)
# TODO: there is much more to do here for packaging - see
# https://dominikberner.ch/cmake-interface-lib/
# http://www.mariobadr.com/creating-a-header-only-library-with-cmake.html
#
# in particular: we presently export no cmake information (no targets are
# exported)
INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ibamr DESTINATION include)
INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/ibtk/include/ibtk DESTINATION include)

#
# 3. Set up specific targets for executables and libraries that want to link to
# IBAMR. As noted above, we have to do this since SAMRAI might not be compiled
# with -fPIC: i.e., if we added SAMRAI as a link dependency to libIBTK and
# libIBTK as a link dependency to libFoo, then the linker would encounter errors
# as a result. With executables we can use the full link interface.
#
ADD_LIBRARY(EXECUTABLE2d INTERFACE)
ADD_LIBRARY(EXECUTABLE3d INTERFACE)
TARGET_LINK_LIBRARIES(EXECUTABLE2d INTERFACE IBAMR2d IBTK2d SAMRAI2d)
TARGET_LINK_LIBRARIES(EXECUTABLE3d INTERFACE IBAMR3d IBTK3d SAMRAI3d)

ADD_LIBRARY(LIBRARY2d INTERFACE)
ADD_LIBRARY(LIBRARY3d INTERFACE)
# no SAMRAI link instructions:
TARGET_LINK_LIBRARIES(LIBRARY2d INTERFACE IBAMR2d IBTK2d)
TARGET_LINK_LIBRARIES(LIBRARY3d INTERFACE IBAMR3d IBTK3d)

# proceed to compilation units:
ADD_SUBDIRECTORY(ibtk)
ADD_SUBDIRECTORY(src)

ADD_SUBDIRECTORY(tests)
