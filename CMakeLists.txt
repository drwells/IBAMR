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
# 0. set the stage:
#
MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")
include(GNUInstallDirs)
INCLUDE(CMakePackageConfigHelpers)
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)

# we need this since SAMRAI is usually statically linked
SET(CMAKE_POSITION_INDEPENDENT_CODE ON)

# print out the relevant command line arguments to make debugging easier
#
# With modern CMake, each find module uses the environment variable ending in
# _ROOT to look for an installation.
SET(_roots "Boost_ROOT" "Eigen3_ROOT" "HDF5_ROOT" "HYPRE_ROOT" "LIBMESH_ROOT"
  "MUPARSER_ROOT" "MPI_ROOT" "PETSC_ROOT" "SAMRAI_ROOT" "SILO_ROOT" "ZLIB_ROOT")
FOREACH(_root ${_roots})
  IF("${${_root}}" STREQUAL "")
    MESSAGE(STATUS "${_root} was not provided to CMake: default search paths will be used.")
  ELSE()
    MESSAGE(STATUS "${_root}=${${_root}}")
  ENDIF()
ENDFOREACH()

#
# 1: manage mandatory dependencies:
#
MESSAGE(STATUS "")
MESSAGE(STATUS "Setting up MPI")
FIND_PACKAGE(MPI REQUIRED)
MESSAGE(STATUS "MPI_C_INCLUDE_DIRS: ${MPI_C_INCLUDE_DIRS}")
MESSAGE(STATUS "MPI_C_LIBRARIES: ${MPI_C_LIBRARIES}")

#
# Boost, which may be bundled:
#
MESSAGE(STATUS "")
MESSAGE(STATUS "Setting up Boost")
FIND_PACKAGE(Boost 1.57 QUIET)
IF(Boost_FOUND)
  SET(IBAMR_USE_BUNDLED_Boost FALSE)
  MESSAGE(STATUS "Found external boost ${Boost_VERSION} at ${Boost_INCLUDE_DIRS}")
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

#
# Eigen, which may be bundled:
#
MESSAGE(STATUS "")
MESSAGE(STATUS "Setting up Eigen")
FIND_PACKAGE(Eigen3 3.2.5 QUIET)
IF(Eigen3_FOUND)
  SET(IBAMR_USE_BUNDLED_Eigen FALSE)
  MESSAGE(STATUS "Found external Eigen ${Eigen3_VERSION} at ${EIGEN3_INCLUDE_DIRS}")
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

#
# muParser, which may be bundled:
#
MESSAGE(STATUS "")
MESSAGE(STATUS "Setting up muParser")
FIND_PATH(MUPARSER_INCLUDE_DIRS NAMES muParser.h HINTS ${MUPARSER_ROOT}/include /usr/include/)
FIND_LIBRARY(MUPARSER_LIBRARIES NAMES muparser HINTS ${MUPARSER_ROOT}/lib /usr/lib)
IF("${MUPARSER_INCLUDE_DIRS}" STREQUAL "MUPARSER_INCLUDE_DIRS-NOTFOUND" OR
   "${MUPARSER_LIBRARIES}" STREQUAL "MUPARSER_LIBRARIES-NOTFOUND")
  SET(IBAMR_USE_BUNDLED_muParser TRUE)
  MESSAGE(STATUS "Setting up muParser as a bundled dependency")
  ADD_LIBRARY(BUNDLED_muParser SHARED)
  TARGET_INCLUDE_DIRECTORIES(
    BUNDLED_muParser
    PUBLIC $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/ibtk/contrib/muparser/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/contrib/muparser>)

  TARGET_SOURCES(BUNDLED_muParser
    PRIVATE
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/src/muParser.cpp
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/src/muParserBase.cpp
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/src/muParserBytecode.cpp
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/src/muParserCallback.cpp
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/src/muParserDLL.cpp
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/src/muParserError.cpp
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/src/muParserInt.cpp
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/src/muParserTest.cpp
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/src/muParserTokenReader.cpp)

  INSTALL(FILES
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParser.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserBase.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserBytecode.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserCallback.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserDLL.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserDef.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserError.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserFixes.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserInt.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserStack.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserTemplateMagic.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserTest.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserToken.h
    ${CMAKE_SOURCE_DIR}/ibtk/contrib/muparser/include/muParserTokenReader.h
    DESTINATION include/contrib/muparser)
  INSTALL(TARGETS BUNDLED_muParser)
ELSE()
  SET(IBAMR_USE_BUNDLED_muParser FALSE)
  MESSAGE(STATUS "MUPARSER_INCLUDE_DIRS: ${MUPARSER_INCLUDE_DIRS}")
  MESSAGE(STATUS "MUPARSER_LIBRARIES: ${MUPARSER_LIBRARIES}")
ENDIF()

#
# HDF5:
#
MESSAGE(STATUS "")
MESSAGE(STATUS "Setting up HDF5")
SET(HDF5_FIND_DEBUG TRUE)
FIND_PACKAGE(HDF5 REQUIRED COMPONENTS C HL)

#
# hypre:
#
MESSAGE(STATUS "")
MESSAGE(STATUS "Setting up HYPRE")
FIND_PATH(HYPRE_INCLUDE_DIRS REQUIRED NAMES HYPRE.h HINTS ${HYPRE_ROOT}/include /usr/include/hypre)
MESSAGE(STATUS "HYPRE_INCLUDE_DIRS: ${HYPRE_INCLUDE_DIRS}")
FIND_LIBRARY(HYPRE_LIBRARIES REQUIRED NAMES HYPRE HINTS ${HYPRE_ROOT}/lib /usr/lib)
MESSAGE(STATUS "HYPRE_LIBRARIES: ${HYPRE_LIBRARIES}")

#
# SAMRAI:
#
# samrai 2.4.4 does not yet have a CMake module
MESSAGE(STATUS "")
MESSAGE(STATUS "Setting up SAMRAI")
SET(_samrai_library_suffixes "algs" "appu" "geom" "hier"
  "math_std" "mesh" "pdat_std" "solv" "xfer")

# SAMRAI consists of about a dozen libraries (see the above list) as well as
# header files for both 2D and 3D. Rather than setting up each as a target, set
# up our own 'header-only' target libraries that depend on SAMRAI libraries to
# encapsulate them neatly:
ADD_LIBRARY(SAMRAI2d INTERFACE IMPORTED)
ADD_LIBRARY(SAMRAI3d INTERFACE IMPORTED)

ADD_LIBRARY(SAMRAI UNKNOWN IMPORTED)
FIND_LIBRARY(SAMRAI_path SAMRAI REQUIRED HINTS "${SAMRAI_ROOT}/lib")

# samrai version info:
FILE(STRINGS "${SAMRAI_ROOT}/include/SAMRAI_config.h" SAMRAI_VERSION_MAJOR_LINE
  REGEX "#define.*SAMRAI_VERSION_MAJOR")
STRING(REGEX REPLACE "^.*SAMRAI_VERSION_MAJOR.* ([0-9]+).*" "\\1"
  SAMRAI_VERSION_MAJOR "${SAMRAI_VERSION_MAJOR_LINE}"
  )
FILE(STRINGS "${SAMRAI_ROOT}/include/SAMRAI_config.h" SAMRAI_VERSION_MINOR_LINE
  REGEX "#define.*SAMRAI_VERSION_MINOR")
STRING(REGEX REPLACE "^.*SAMRAI_VERSION_MINOR.* ([0-9]+).*" "\\1"
  SAMRAI_VERSION_MINOR "${SAMRAI_VERSION_MINOR_LINE}"
  )
FILE(STRINGS "${SAMRAI_ROOT}/include/SAMRAI_config.h" SAMRAI_VERSION_PATCHLEVEL_LINE
  REGEX "#define.*SAMRAI_VERSION_PATCHLEVEL")
STRING(REGEX REPLACE "^.*SAMRAI_VERSION_PATCHLEVEL.* ([0-9]+).*" "\\1"
  SAMRAI_VERSION_PATCHLEVEL "${SAMRAI_VERSION_PATCHLEVEL_LINE}"
  )
SET(SAMRAI_VERSION_STRING "${SAMRAI_VERSION_MAJOR}.${SAMRAI_VERSION_MINOR}\
.${SAMRAI_VERSION_PATCHLEVEL}")
MESSAGE(STATUS "Found SAMRAI ${SAMRAI_VERSION_STRING} at ${SAMRAI_ROOT}")
SET_TARGET_PROPERTIES(SAMRAI PROPERTIES IMPORTED_LOCATION ${SAMRAI_path})

# samrai libs:
TARGET_LINK_LIBRARIES(SAMRAI2d INTERFACE SAMRAI)
TARGET_LINK_LIBRARIES(SAMRAI3d INTERFACE SAMRAI)

FOREACH(_d ${IBAMR_DIMENSIONS})
  FOREACH(_suffix ${_samrai_library_suffixes})
    SET(_lib_name "SAMRAI${_d}d_${_suffix}")
    ADD_LIBRARY( ${_lib_name} UNKNOWN IMPORTED)

    FIND_LIBRARY("${_lib_name}_path" ${_lib_name} REQUIRED
      HINTS "${SAMRAI_ROOT}/lib")
    IF("${${_lib_name}_path}" STREQUAL "${_lib_name}-NOTFOUND")
      MESSAGE(FATAL_ERROR("\
Unable to find required library ${_lib_name} in directory ${SAMRAI_ROOT}/lib"))
    ENDIF()

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
MESSAGE(STATUS "")
FIND_PACKAGE(PkgConfig REQUIRED)

MESSAGE(STATUS "")
MESSAGE(STATUS "Setting up PETSc")
LIST(PREPEND CMAKE_PREFIX_PATH "${PETSC_ROOT}")
PKG_CHECK_MODULES(PETSC REQUIRED IMPORTED_TARGET GLOBAL PETSc)
LIST(POP_FRONT CMAKE_PREFIX_PATH)

#
# 1.5: manage optional dependencies
#
MESSAGE(STATUS "")
MESSAGE(STATUS "Setting up libMesh")
LIST(PREPEND CMAKE_PREFIX_PATH "${LIBMESH_ROOT}")
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

MESSAGE(STATUS "")
MESSAGE(STATUS "Setting up Silo")
SET(IBAMR_WITH_SILO FALSE)
FIND_PATH(SILO_INCLUDE_DIRS NAMES silo.h HINTS ${SILO_ROOT}/include /usr/include/)
FIND_LIBRARY(SILO_LIBRARIES NAMES silo siloh5 HINTS ${SILO_ROOT}/lib /usr/lib)
IF(NOT "${SILO_LIBRARIES}" STREQUAL "SILO_LIBRARIES-NOTFOUND" AND
    NOT "${SILO_INCLUDE_DIRS}" STREQUAL "SILO_INCLUDE_DIRS-NOTFOUND")
  MESSAGE(STATUS "SILO_INCLUDE_DIRS: ${SILO_INCLUDE_DIRS}")
  MESSAGE(STATUS "SILO_LIBRARIES: ${SILO_LIBRARIES}")
  SET(IBAMR_WITH_SILO TRUE)

  # Silo depends on zlib
  FIND_PACKAGE(ZLIB REQUIRED)
ENDIF()


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

  SET(_args "-DSAMRAI_FORTDIR=${SAMRAI_ROOT}/include;-DCURRENT_SRCDIR=${_in_directory};-DTOP_SRCDIR=${CMAKE_SOURCE_DIR}")
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
  TARGET_COMPILE_FEATURES(${target_library} PUBLIC cxx_std_14)
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
  IF(IBAMR_USE_BUNDLED_muParser)
    TARGET_LINK_LIBRARIES(${target_library} PUBLIC BUNDLED_muParser)
  ELSE()
    # todo: no nice target support yet for muParser
    TARGET_LINK_LIBRARIES(${target_library} PUBLIC "${MUPARSER_LIBRARIES}")
    TARGET_INCLUDE_DIRECTORIES(${target_library} PUBLIC "${MUPARSER_INCLUDE_DIRS}")
  ENDIF()
  # PETSc:
  TARGET_LINK_LIBRARIES(${target_library} PUBLIC PkgConfig::PETSC)
  # libMesh:
  IF(IBAMR_WITH_LIBMESH)
    TARGET_LINK_LIBRARIES(${target_library} PUBLIC PkgConfig::LIBMESH)
  ENDIF()
  # Silo:
  IF(IBAMR_WITH_SILO)
    TARGET_LINK_LIBRARIES(${target_library} PRIVATE "${SILO_LIBRARIES}")
    TARGET_INCLUDE_DIRECTORIES(${target_library} PRIVATE "${SILO_INCLUDE_DIRS}")
    TARGET_LINK_LIBRARIES(${target_library} PRIVATE ZLIB::ZLIB)
  ENDIF()

  # Figure out if we are a 2D or a 3D library:
  SET(_2_location = -1)
  SET(_3_location = -1)
  STRING(FIND "${target_library}" "2" _2_location)
  STRING(FIND "${target_library}" "3" _3_location)

  # We cannot add SAMRAI to the link interface for our own libraries since
  # SAMRAI doesn't necessarily compile with -fPIC: we will instead put it in a
  # link interface for executables.
  TARGET_INCLUDE_DIRECTORIES(${target_library} PUBLIC "${SAMRAI_ROOT}/include")

  # add the right preprocessor flag:
  IF(NOT ${_2_location} STREQUAL "-1")
    MESSAGE(STATUS "Adding flag -DNDIM=2 to target ${target_library}")
    TARGET_COMPILE_OPTIONS(${target_library} PUBLIC -DNDIM=2)
  TARGET_LINK_LIBRARIES(${target_library} PUBLIC SAMRAI2d)
  ELSEIF(NOT ${_3_location} STREQUAL "-1")
    MESSAGE(STATUS "Adding flag -DNDIM=3 to target ${target_library}")
    TARGET_COMPILE_OPTIONS(${target_library} PUBLIC -DNDIM=3)
  TARGET_LINK_LIBRARIES(${target_library} PUBLIC SAMRAI3d)
  ELSE()
    MESSAGE(STATUS "no NDIM preprocessor constant added to ${target_library}")
  ENDIF()
ENDFUNCTION()

MESSAGE(STATUS "")
MESSAGE(STATUS "IBAMR dependencies have been successfully set up.")

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
