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

ADD_LIBRARY(IBTK2d SHARED)
ADD_LIBRARY(IBTK3d SHARED)

# 'base' names (without, 2d.f.m4 or 3d.f.m4) for all the fortran that is
# explicitly compiled
SET(FORTRAN_SRC_BASE
  boundary/cf_interface/fortran/linearcfinterpolation
  boundary/cf_interface/fortran/quadcfinterpolation
  boundary/physical_boundary/fortran/cartphysbdryop

  coarsen_ops/fortran/cubiccoarsen
  coarsen_ops/fortran/rt0coarsen
  # lagrangian/fortran/lagrangian_delta is included in lagrangian_interaction
  lagrangian/fortran/lagrangian_interaction

  math/fortran/curl
  math/fortran/div
  math/fortran/flux
  math/fortran/grad
  math/fortran/graddetect
  math/fortran/interp
  math/fortran/laplace
  math/fortran/miscmath
  math/fortran/rot
  math/fortran/strain
  math/fortran/vclaplace

  refine_ops/fortran/cart_side_refine
  refine_ops/fortran/divpreservingrefine
  # refine_ops/fortran/cart_side_helpers is included in cart_side_refine

  solvers/impls/fortran/patchsmoothers
  )

FOREACH(_d ${IBAMR_DIMENSIONS})
  FOREACH(_fortran_root_name ${FORTRAN_SRC_BASE})
    IBAMR_PROCESS_M4(NAME ${_fortran_root_name}${_d}d.f)
  ENDFOREACH()

  SET("FORTRAN_GENERATED_SRC${_d}D" "")
  FOREACH(_fortran_root_name ${FORTRAN_SRC_BASE})
    LIST(APPEND "FORTRAN_GENERATED_SRC${_d}D" "${CMAKE_CURRENT_BINARY_DIR}/${_fortran_root_name}${_d}d.f")
  ENDFOREACH()
ENDFOREACH()

TARGET_SOURCES(IBTK2d PUBLIC ${FORTRAN_GENERATED_SRC2D})
TARGET_SOURCES(IBTK3d PUBLIC ${FORTRAN_GENERATED_SRC3D})

TARGET_LINK_LIBRARIES(IBTK2d PUBLIC IBTKHeaders)
TARGET_LINK_LIBRARIES(IBTK3d PUBLIC IBTKHeaders)

IBAMR_SETUP_TARGET_LIBRARY(IBTK2d)
IBAMR_SETUP_TARGET_LIBRARY(IBTK3d)

INSTALL(TARGETS IBTK2d IBTK3d COMPONENT library)
