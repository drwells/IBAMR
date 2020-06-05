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

#
# 1. manage mandatory dependencies:
#
FIND_PACKAGE(Boost 1.57 REQUIRED)
FIND_PACKAGE(MPI REQUIRED)

SET(HDF5_FIND_DEBUG TRUE)
FIND_PACKAGE(HDF5 REQUIRED COMPONENTS C HL)

#
# 1.5. manage optional dependencies
#

#
# 2. IBAMR-specific configuration:
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

  # todo - un-hardcode this when we properly support SAMRAI
  SET(_args "-DSAMRAI_FORTDIR=/home/drwells/Applications/samrai-2.4.4/include;-DCURRENT_SRCDIR=${_in_directory};-DTOP_SRCDIR=${CMAKE_SOURCE_DIR}")
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
  # we need C++11
  SET_PROPERTY(TARGET ${target_library} PROPERTY CXX_STANDARD 11)
  # we and our users will use these MPI functions so make the interface public:
  TARGET_LINK_LIBRARIES(${target_library} PUBLIC MPI::MPI_C)
  # same for HDF5:
  # todo: no nice target support yet for HDF5
  TARGET_LINK_LIBRARIES(${target_library} PUBLIC "${HDF5_LIBRARIES}")
  TARGET_INCLUDE_DIRECTORIES(${target_library} PUBLIC "${HDF5_INCLUDE_DIRS}")
  # we only need boost headers
  TARGET_LINK_LIBRARIES(${target_library} PUBLIC Boost::headers)
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
#
# in particular: we presently export no cmake information (no targets are
# exported)
INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ibamr DESTINATION include)
INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/ibtk/include/ibtk DESTINATION include)

# proceed to compilation units:
ADD_SUBDIRECTORY(ibtk)
ADD_SUBDIRECTORY(src)
