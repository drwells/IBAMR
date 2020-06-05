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

ADD_LIBRARY(IBAMR2d SHARED)
ADD_LIBRARY(IBAMR3d SHARED)

# 'base' names (without, 2d.f.m4 or 3d.f.m4) for all the fortran that is
# explicitly compiled
SET(FORTRAN_SRC_BASE
  adv_diff/fortran/adv_diff_consdiff
  adv_diff/fortran/adv_diff_wp_convective_op
  advect/fortran/advect_centered_derivatives
  advect/fortran/advect_detect
  advect/fortran/advect_diff
  advect/fortran/advect_helpers
  advect/fortran/advect_predictors
  advect/fortran/advect_stable
  complex_fluids/fortran/div_tensor
  complex_fluids/fortran/log_upper_convective_op
  complex_fluids/fortran/sqrt_upper_convective_op
  complex_fluids/fortran/upper_convective_op
  level_set/fortran/levelsetops
  navier_stokes/fortran/navier_stokes_bdryop
  navier_stokes/fortran/navier_stokes_divsource
  navier_stokes/fortran/navier_stokes_stabledt
  navier_stokes/fortran/navier_stokes_staggered_derivatives
  navier_stokes/fortran/navier_stokes_staggered_helpers
  navier_stokes/fortran/navier_stokes_stochastic_forcing
  navier_stokes/fortran/navier_stokes_surface_tension_forcing
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

TARGET_SOURCES(IBAMR2d PUBLIC ${FORTRAN_GENERATED_SRC2D})
TARGET_SOURCES(IBAMR3d PUBLIC ${FORTRAN_GENERATED_SRC3D})

TARGET_LINK_LIBRARIES(IBAMR2d PUBLIC IBAMRHeaders IBTKHeaders)
TARGET_LINK_LIBRARIES(IBAMR3d PUBLIC IBAMRHeaders IBTKHeaders)

IBAMR_SETUP_TARGET_LIBRARY(IBAMR2d)
IBAMR_SETUP_TARGET_LIBRARY(IBAMR3d)

INSTALL(TARGETS IBAMR2d IBAMR3d COMPONENT library)
