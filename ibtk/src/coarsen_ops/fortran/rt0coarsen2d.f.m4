c ---------------------------------------------------------------------
c
c Copyright (c) 2015 - 2017 by the IBAMR developers
c All rights reserved.
c
c This file is part of IBAMR.
c
c IBAMR is free software and is distributed under the 3-clause BSD
c license. The full text of the license can be found in the file
c COPYRIGHT at the top level directory of IBAMR.
c
c ---------------------------------------------------------------------

define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
define(coarsen_index,`dnl
if ($1.lt.0) then
            $2=($1+1)/$3-1
         else
            $2=$1/$3
         endif
')dnl'c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Coarsen side-centered data via the adjoint of RT0 interpolation.
c
c     NOTE: Values at physical boundaries and coarse-fine interfaces
c     will need to be corrected.  This routine uses the correct stencil
c     only *away* from such boundaries.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scrt0coarsen2d(
     &     U_crse0,U_crse1,U_crse_gcw,
     &     U_fine0,U_fine1,U_fine_gcw,
     &     ilowerc0,iupperc0,
     &     ilowerc1,iupperc1,
     &     ilowerf0,iupperf0,
     &     ilowerf1,iupperf1,
     &     ratio_to_coarser,
     &     coarse_box_lower,coarse_box_upper)
c
      implicit none
c
c     Input.
c
      INTEGER ilowerc0,iupperc0
      INTEGER ilowerc1,iupperc1
      INTEGER ilowerf0,iupperf0
      INTEGER ilowerf1,iupperf1
      INTEGER U_crse_gcw,U_fine_gcw

      INTEGER ratio_to_coarser(0:NDIM-1)

      INTEGER coarse_box_lower(0:NDIM-1), coarse_box_upper(0:NDIM-1)

      REAL U_fine0(
     &     SIDE2d0(ilowerf,iupperf,U_fine_gcw)
     &     )
      REAL U_fine1(
     &     SIDE2d1(ilowerf,iupperf,U_fine_gcw)
     &     )
c
c     Input/Output.
c
      REAL U_crse0(
     &     SIDE2d0(ilowerc,iupperc,U_crse_gcw)
     &     )
      REAL U_crse1(
     &     SIDE2d1(ilowerc,iupperc,U_crse_gcw)
     &     )
c
c     Local variables.
c
      INTEGER i,i_c,i_f,j,j_c,j_f
      INTEGER axis,d,istencil_lower(0:NDIM-1),istencil_upper(0:NDIM-1)
      REAL w,w_fac
c
c     Treat x components:
c
      axis = 0
      do d = 0,NDIM-1
         if (d .eq. axis) then
            istencil_lower(d) = -ratio_to_coarser(d)+1
            istencil_upper(d) = +ratio_to_coarser(d)-1
         else
            istencil_lower(d) = 0
            istencil_upper(d) = ratio_to_coarser(d)-1
         endif
      enddo

      w_fac = 0.d0
      do j = istencil_lower(1),istencil_upper(1)
         do i = istencil_lower(0),istencil_upper(0)
            w = 1.d0 - dabs(dble(i))/dble(ratio_to_coarser(axis))
            w_fac = w_fac + w
         enddo
      enddo
      w_fac = 1.d0/w_fac

      do j_c = coarse_box_lower(1),coarse_box_upper(1)
         do i_c = coarse_box_lower(0),coarse_box_upper(0)+1
            i_f = i_c*ratio_to_coarser(0)
            j_f = j_c*ratio_to_coarser(1)
            U_crse0(i_c,j_c) = 0.d0
            do j = istencil_lower(1),istencil_upper(1)
               do i = istencil_lower(0),istencil_upper(0)
                  w = 1.d0 - dabs(dble(i))/dble(ratio_to_coarser(axis))
                  U_crse0(i_c,j_c) = U_crse0(i_c,j_c) +
     &                 w*w_fac*U_fine0(i_f+i,j_f+j)
               enddo
            enddo
         enddo
      enddo
c
c     Treat y components:
c
      axis = 1
      do d = 0,NDIM-1
         if (d .eq. axis) then
            istencil_lower(d) = -ratio_to_coarser(d)+1
            istencil_upper(d) = +ratio_to_coarser(d)-1
         else
            istencil_lower(d) = 0
            istencil_upper(d) = ratio_to_coarser(d)-1
         endif
      enddo

      w_fac = 0.d0
      do j = istencil_lower(1),istencil_upper(1)
         do i = istencil_lower(0),istencil_upper(0)
            w = 1.d0 - dabs(dble(j))/dble(ratio_to_coarser(axis))
            w_fac = w_fac + w
         enddo
      enddo
      w_fac = 1.d0/w_fac

      do j_c = coarse_box_lower(1),coarse_box_upper(1)+1
         do i_c = coarse_box_lower(0),coarse_box_upper(0)
            i_f = i_c*ratio_to_coarser(0)
            j_f = j_c*ratio_to_coarser(1)
            U_crse1(i_c,j_c) = 0.d0
            do j = istencil_lower(1),istencil_upper(1)
               do i = istencil_lower(0),istencil_upper(0)
                  w = 1.d0 - dabs(dble(j))/dble(ratio_to_coarser(axis))
                  U_crse1(i_c,j_c) = U_crse1(i_c,j_c) +
     &                 w*w_fac*U_fine1(i_f+i,j_f+j)
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Fix up RT0 coarsening along a boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scrt0coarsenbdry2d(
     &     U_crse0,U_crse1,U_crse_gcw,
     &     U_fine0,U_fine1,U_fine_gcw,
     &     ilowerc0,iupperc0,
     &     ilowerc1,iupperc1,
     &     ilowerf0,iupperf0,
     &     ilowerf1,iupperf1,
     &     ratio_to_coarser,
     &     coarse_box_lower,coarse_box_upper,
     &     bbox_ilowerf0,bbox_iupperf0,
     &     bbox_ilowerf1,bbox_iupperf1,
     &     axis,upperlower)
c
      implicit none
c
c     Input.
c
      INTEGER ilowerc0,iupperc0
      INTEGER ilowerc1,iupperc1
      INTEGER ilowerf0,iupperf0
      INTEGER ilowerf1,iupperf1
      INTEGER U_crse_gcw,U_fine_gcw

      INTEGER ratio_to_coarser(0:NDIM-1)

      INTEGER coarse_box_lower(0:NDIM-1), coarse_box_upper(0:NDIM-1)

      INTEGER bbox_ilowerf0,bbox_iupperf0
      INTEGER bbox_ilowerf1,bbox_iupperf1
      INTEGER axis,upperlower

      REAL U_fine0(
     &     SIDE2d0(ilowerf,iupperf,U_fine_gcw)
     &     )
      REAL U_fine1(
     &     SIDE2d1(ilowerf,iupperf,U_fine_gcw)
     &     )
c
c     Input/Output.
c
      REAL U_crse0(
     &     SIDE2d0(ilowerc,iupperc,U_crse_gcw)
     &     )
      REAL U_crse1(
     &     SIDE2d1(ilowerc,iupperc,U_crse_gcw)
     &     )
c
c     Local variables.
c
      INTEGER i,i_c,i_f,j,j_c,j_f
      INTEGER bbox_ilowerc(0:NDIM-1),bbox_iupperc(0:NDIM-1)
      INTEGER d,istencil_lower(0:NDIM-1),istencil_upper(0:NDIM-1)
      INTEGER ibdryc(0:NDIM-1)
      REAL w,w_fac
c
c     Setup boundary box extents.
c
      coarsen_index(bbox_ilowerf0, bbox_ilowerc(0), ratio_to_coarser(0))
      coarsen_index(bbox_ilowerf1, bbox_ilowerc(1), ratio_to_coarser(1))
      coarsen_index(bbox_iupperf0, bbox_iupperc(0), ratio_to_coarser(0))
      coarsen_index(bbox_iupperf1, bbox_iupperc(1), ratio_to_coarser(1))

      do d = 0,NDIM-1
         bbox_ilowerc(d) = max(coarse_box_lower(d), bbox_ilowerc(d))
         bbox_iupperc(d) = min(coarse_box_upper(d), bbox_iupperc(d))
      enddo
c
c     Setup stencil indices.
c
      do d = 0,NDIM-1
         if (d .eq. axis) then
            if (upperlower .eq. 0)  then
               istencil_lower(d) = 0
               istencil_upper(d) = ratio_to_coarser(d)-1
               ibdryc(d) = bbox_iupperc(d)+1
            else if (upperlower .eq. 1)  then
               istencil_lower(d) = -ratio_to_coarser(d)+1
               istencil_upper(d) = 0
               ibdryc(d) = bbox_ilowerc(d)
            endif
         else
            istencil_lower(d) = 0
            istencil_upper(d) = ratio_to_coarser(d)-1
         endif
      enddo
c
c     Set values along the boundary.
c
      if (axis .eq. 0) then
c
c     Treat x boundaries:
c
         w_fac = 0.d0
         do j = istencil_lower(1),istencil_upper(1)
            do i = istencil_lower(0),istencil_upper(0)
               w = 1.d0 - dabs(dble(i))/dble(ratio_to_coarser(axis))
               w_fac = w_fac + w
            enddo
         enddo
         w_fac = 1.d0/w_fac

         do j_c = bbox_ilowerc(1),bbox_iupperc(1)
            i_c = ibdryc(0)
            i_f = i_c*ratio_to_coarser(0)
            j_f = j_c*ratio_to_coarser(1)
            U_crse0(i_c,j_c) = 0.d0
            do j = istencil_lower(1),istencil_upper(1)
               do i = istencil_lower(0),istencil_upper(0)
                  w = 1.d0 - dabs(dble(i))/dble(ratio_to_coarser(axis))
                  U_crse0(i_c,j_c) = U_crse0(i_c,j_c) +
     &               w*w_fac*U_fine0(i_f+i,j_f+j)
               enddo
            enddo
         enddo

      else if (axis .eq. 1) then
c
c     Treat y boundaries:
c
         w_fac = 0.d0
         do j = istencil_lower(1),istencil_upper(1)
            do i = istencil_lower(0),istencil_upper(0)
               w = 1.d0 - dabs(dble(j))/dble(ratio_to_coarser(axis))
               w_fac = w_fac + w
            enddo
         enddo
         w_fac = 1.d0/w_fac

         j_c = ibdryc(1)
         do i_c = bbox_ilowerc(0),bbox_iupperc(0)
            i_f = i_c*ratio_to_coarser(0)
            j_f = j_c*ratio_to_coarser(1)
            U_crse1(i_c,j_c) = 0.d0
            do j = istencil_lower(1),istencil_upper(1)
               do i = istencil_lower(0),istencil_upper(0)
                  w = 1.d0 - dabs(dble(j))/dble(ratio_to_coarser(axis))
                  U_crse1(i_c,j_c) = U_crse1(i_c,j_c) +
     &               w*w_fac*U_fine1(i_f+i,j_f+j)
               enddo
            enddo
         enddo

      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
