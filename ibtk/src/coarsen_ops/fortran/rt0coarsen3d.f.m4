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

define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Coarsen side-centered data via the adjoint of RT0 interpolation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scrt0coarsen3d(
     &     U_crse0,U_crse1,U_crse2,U_crse_gcw,
     &     U_fine0,U_fine1,U_fine2,U_fine_gcw,
     &     ilowerc0,iupperc0,
     &     ilowerc1,iupperc1,
     &     ilowerc2,iupperc2,
     &     ilowerf0,iupperf0,
     &     ilowerf1,iupperf1,
     &     ilowerf2,iupperf2,
     &     ratio_to_coarser,
     &     coarse_box_lower,coarse_box_upper,
     &     physical_bdry_flags)
c
      implicit none
c
c     Input.
c
      INTEGER ilowerc0,iupperc0
      INTEGER ilowerc1,iupperc1
      INTEGER ilowerc2,iupperc2
      INTEGER ilowerf0,iupperf0
      INTEGER ilowerf1,iupperf1
      INTEGER ilowerf2,iupperf2
      INTEGER U_crse_gcw,U_fine_gcw

      INTEGER ratio_to_coarser(0:NDIM-1)

      INTEGER coarse_box_lower(0:NDIM-1), coarse_box_upper(0:NDIM-1)

      INTEGER physical_bdry_flags(0:NDIM-1,0:1)

      REAL U_fine0(
     &     SIDE3d0(ilowerf,iupperf,U_fine_gcw)
     &     )
      REAL U_fine1(
     &     SIDE3d1(ilowerf,iupperf,U_fine_gcw)
     &     )
      REAL U_fine2(
     &     SIDE3d2(ilowerf,iupperf,U_fine_gcw)
     &     )
c
c     Input/Output.
c
      REAL U_crse0(
     &     SIDE3d0(ilowerc,iupperc,U_crse_gcw)
     &     )
      REAL U_crse1(
     &     SIDE3d1(ilowerc,iupperc,U_crse_gcw)
     &     )
      REAL U_crse2(
     &     SIDE3d2(ilowerc,iupperc,U_crse_gcw)
     &     )
c
c     Local variables.
c
      INTEGER i,i_c,i_f,j,j_c,j_f,k,k_c,k_f
      REAL w0,w1,w2,w_fac
c
c     Coarsen data.
c

c
c     Treat x components away from physical boundaries:
c
      w_fac = 0.d0
      do k = 0,ratio_to_coarser(2)-1
         do j = 0,ratio_to_coarser(1)-1
            do i = -ratio_to_coarser(0)+1,ratio_to_coarser(0)-1
               w0 = 1.d0 - dabs(dble(i))/dble(ratio_to_coarser(0))
               w1 = 1.d0
               w2 = 1.d0
               w_fac = w_fac + w0*w1*w2
            enddo
         enddo
      enddo
      w_fac = 1.d0/w_fac

      do k_c = coarse_box_lower(2),coarse_box_upper(2)
         do j_c = coarse_box_lower(1),coarse_box_upper(1)
            do i_c = coarse_box_lower(0),coarse_box_upper(0)+1
               i_f = i_c*ratio_to_coarser(0)
               j_f = j_c*ratio_to_coarser(1)
               k_f = k_c*ratio_to_coarser(2)
               U_crse0(i_c,j_c,k_c) = 0.d0
               do k = 0,ratio_to_coarser(2)-1
                  do j = 0,ratio_to_coarser(1)-1
                     do i = -ratio_to_coarser(0)+1,ratio_to_coarser(0)-1
                        w0 = 1.d0 - dabs(dble(i))
     &                       /dble(ratio_to_coarser(0))
                        w1 = 1.d0
                        w2 = 1.d0
                        U_crse0(i_c,j_c,k_c) = U_crse0(i_c,j_c,k_c) +
     &                       w0*w1*w2*w_fac*U_fine0(i_f+i,j_f+j,k_f+k)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
c     Treat x components along physical boundaries:
c
      if (physical_bdry_flags(0,0) .eq. 1) then ! lower x patch bdry is on physical bdry
         w_fac = 0.d0
         do k = 0,ratio_to_coarser(2)-1
            do j = 0,ratio_to_coarser(1)-1
               do i = 0,ratio_to_coarser(0)-1
                  w0 = 1.d0 - dabs(dble(i))/dble(ratio_to_coarser(0))
                  w1 = 1.d0
                  w2 = 1.d0
                  w_fac = w_fac + w0*w1*w2
               enddo
            enddo
         enddo
         w_fac = 1.d0/w_fac

         do k_c = coarse_box_lower(2),coarse_box_upper(2)
            do j_c = coarse_box_lower(1),coarse_box_upper(1)
               i_c = ilowerc0
               i_f = i_c*ratio_to_coarser(0)
               j_f = j_c*ratio_to_coarser(1)
               k_f = k_c*ratio_to_coarser(2)
               U_crse0(i_c,j_c,k_c) = 0.d0
               do k = 0,ratio_to_coarser(2)-1
                  do j = 0,ratio_to_coarser(1)-1
                     do i = 0,ratio_to_coarser(0)-1
                        w0 = 1.d0 - dabs(dble(i))
     &                       /dble(ratio_to_coarser(0))
                        w1 = 1.d0
                        w2 = 1.d0
                        U_crse0(i_c,j_c,k_c) = U_crse0(i_c,j_c,k_c) +
     &                       w0*w1*w2*w_fac*U_fine0(i_f+i,j_f+j,k_f+k)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif

      if (physical_bdry_flags(0,1) .eq. 1) then ! upper x patch bdry is on physical bdry
         w_fac = 0.d0
         do k = 0,ratio_to_coarser(2)-1
            do j = 0,ratio_to_coarser(1)-1
               do i = -ratio_to_coarser(0)+1,0
                  w0 = 1.d0 - dabs(dble(i))/dble(ratio_to_coarser(0))
                  w1 = 1.d0
                  w2 = 1.d0
                  w_fac = w_fac + w0*w1*w2
               enddo
            enddo
         enddo
         w_fac = 1.d0/w_fac

         do k_c = coarse_box_lower(2),coarse_box_upper(2)
            do j_c = coarse_box_lower(1),coarse_box_upper(1)
               i_c = iupperc0
               i_f = i_c*ratio_to_coarser(0)
               j_f = j_c*ratio_to_coarser(1)
               k_f = k_c*ratio_to_coarser(2)
               U_crse0(i_c,j_c,k_c) = 0.d0
               do k = 0,ratio_to_coarser(2)-1
                  do j = 0,ratio_to_coarser(1)-1
                     do i = -ratio_to_coarser(0)+1,0
                        w0 = 1.d0 - dabs(dble(i))
     &                       /dble(ratio_to_coarser(0))
                        w1 = 1.d0
                        w2 = 1.d0
                        U_crse0(i_c,j_c,k_c) = U_crse0(i_c,j_c,k_c) +
     &                       w0*w1*w2*w_fac*U_fine0(i_f+i,j_f+j,k_f+k)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
c
c     Treat y components away from physical boundaries:
c
      w_fac = 0.d0
      do k = 0,ratio_to_coarser(2)-1
         do j = -ratio_to_coarser(1)+1,ratio_to_coarser(1)-1
            do i = 0,ratio_to_coarser(0)-1
               w0 = 1.d0
               w1 = 1.d0 - dabs(dble(j))/dble(ratio_to_coarser(1))
               w2 = 1.d0
               w_fac = w_fac + w0*w1*w2
            enddo
         enddo
      enddo
      w_fac = 1.d0/w_fac

      do k_c = coarse_box_lower(2),coarse_box_upper(2)
         do j_c = coarse_box_lower(1),coarse_box_upper(1)+1
            do i_c = coarse_box_lower(0),coarse_box_upper(0)
               i_f = i_c*ratio_to_coarser(0)
               j_f = j_c*ratio_to_coarser(1)
               k_f = k_c*ratio_to_coarser(2)
               U_crse1(i_c,j_c,k_c) = 0.d0
               do k = 0,ratio_to_coarser(2)-1
                  do j = -ratio_to_coarser(1)+1,ratio_to_coarser(1)-1
                     do i = 0,ratio_to_coarser(0)-1
                        w0 = 1.d0
                        w1 = 1.d0 - dabs(dble(j))
     &                       /dble(ratio_to_coarser(1))
                        w2 = 1.d0
                        U_crse1(i_c,j_c,k_c) = U_crse1(i_c,j_c,k_c) +
     &                       w0*w1*w2*w_fac*U_fine1(i_f+i,j_f+j,k_f+k)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
c     Treat y components along physical boundaries:
c
      if (physical_bdry_flags(1,0) .eq. 1) then ! lower y patch bdry is on physical bdry
         w_fac = 0.d0
         do k = 0,ratio_to_coarser(2)-1
            do j = 0,ratio_to_coarser(1)-1
               do i = 0,ratio_to_coarser(0)-1
                  w0 = 1.d0
                  w1 = 1.d0 - dabs(dble(j))/dble(ratio_to_coarser(1))
                  w2 = 1.d0
                  w_fac = w_fac + w0*w1*w2
               enddo
            enddo
         enddo
         w_fac = 1.d0/w_fac

         do k_c = coarse_box_lower(2),coarse_box_upper(2)
            j_c = ilowerc1
            do i_c = coarse_box_lower(0),coarse_box_upper(0)
               i_f = i_c*ratio_to_coarser(0)
               j_f = j_c*ratio_to_coarser(1)
               k_f = k_c*ratio_to_coarser(2)
               U_crse1(i_c,j_c,k_c) = 0.d0
               do k = 0,ratio_to_coarser(2)-1
                  do j = 0,ratio_to_coarser(1)-1
                     do i = 0,ratio_to_coarser(0)-1
                        w0 = 1.d0
                        w1 = 1.d0 - dabs(dble(j))
     &                       /dble(ratio_to_coarser(1))
                        w2 = 1.d0
                        U_crse1(i_c,j_c,k_c) = U_crse1(i_c,j_c,k_c) +
     &                       w0*w1*w2*w_fac*U_fine1(i_f+i,j_f+j,k_f+k)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif

      if (physical_bdry_flags(1,1) .eq. 1) then ! upper y patch bdry is on physical bdry
         w_fac = 0.d0
         do k = 0,ratio_to_coarser(2)-1
            do j = -ratio_to_coarser(1)+1,0
               do i = 0,ratio_to_coarser(0)-1
                  w0 = 1.d0
                  w1 = 1.d0 - dabs(dble(j))/dble(ratio_to_coarser(1))
                  w2 = 1.d0
                  w_fac = w_fac + w0*w1*w2
               enddo
            enddo
         enddo
         w_fac = 1.d0/w_fac

         do k_c = coarse_box_lower(2),coarse_box_upper(2)
            j_c = iupperc1
            do i_c = coarse_box_lower(0),coarse_box_upper(0)
               i_f = i_c*ratio_to_coarser(0)
               j_f = j_c*ratio_to_coarser(1)
               k_f = k_c*ratio_to_coarser(2)
               U_crse1(i_c,j_c,k_c) = 0.d0
               do k = 0,ratio_to_coarser(2)-1
                  do j = -ratio_to_coarser(1)+1,0
                     do i = 0,ratio_to_coarser(0)-1
                        w0 = 1.d0
                        w1 = 1.d0 - dabs(dble(j))
     &                       /dble(ratio_to_coarser(2))
                        w2 = 1.d0
                        U_crse1(i_c,j_c,k_c) = U_crse1(i_c,j_c,k_c) +
     &                       w0*w1*w2*w_fac*U_fine2(i_f+i,j_f+j,k_f+k)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
c
c     Treat z components away from physical boundaries:
c
      w_fac = 0.d0
      do k = -ratio_to_coarser(2)+1,ratio_to_coarser(2)-1
         do j = 0,ratio_to_coarser(1)-1
            do i = 0,ratio_to_coarser(0)-1
               w0 = 1.d0
               w1 = 1.d0
               w2 = 1.d0 - dabs(dble(k))/dble(ratio_to_coarser(2))
               w_fac = w_fac + w0*w1*w2
            enddo
         enddo
      enddo
      w_fac = 1.d0/w_fac

      do k_c = coarse_box_lower(2),coarse_box_upper(2)+1
         do j_c = coarse_box_lower(1),coarse_box_upper(1)
            do i_c = coarse_box_lower(0),coarse_box_upper(0)
               i_f = i_c*ratio_to_coarser(0)
               j_f = j_c*ratio_to_coarser(1)
               k_f = k_c*ratio_to_coarser(2)
               U_crse2(i_c,j_c,k_c) = 0.d0
               do k = -ratio_to_coarser(2)+1,ratio_to_coarser(2)-1
                  do j = 0,ratio_to_coarser(1)-1
                     do i = 0,ratio_to_coarser(0)-1
                        w0 = 1.d0
                        w1 = 1.d0
                        w2 = 1.d0 - dabs(dble(k))
     &                       /dble(ratio_to_coarser(2))
                        U_crse2(i_c,j_c,k_c) = U_crse2(i_c,j_c,k_c) +
     &                       w0*w1*w2*w_fac*U_fine2(i_f+i,j_f+j,k_f+k)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
c     Treat z components along physical boundaries:
c
      if (physical_bdry_flags(2,0) .eq. 1) then ! lower z patch bdry is on physical bdry
         w_fac = 0.d0
         do k = 0,ratio_to_coarser(2)-1
            do j = 0,ratio_to_coarser(1)-1
               do i = 0,ratio_to_coarser(0)-1
                  w0 = 1.d0
                  w1 = 1.d0
                  w2 = 1.d0 - dabs(dble(k))/dble(ratio_to_coarser(2))
                  w_fac = w_fac + w0*w1*w2
               enddo
            enddo
         enddo
         w_fac = 1.d0/w_fac

         k_c = ilowerc2
         do j_c = coarse_box_lower(1),coarse_box_upper(1)
            do i_c = coarse_box_lower(0),coarse_box_upper(0)
               i_f = i_c*ratio_to_coarser(0)
               j_f = j_c*ratio_to_coarser(1)
               k_f = k_c*ratio_to_coarser(2)
               U_crse2(i_c,j_c,k_c) = 0.d0
               do k = 0,ratio_to_coarser(2)-1
                  do j = 0,ratio_to_coarser(1)-1
                     do i = 0,ratio_to_coarser(0)-1
                        w0 = 1.d0
                        w1 = 1.d0
                        w2 = 1.d0 - dabs(dble(k))
     &                       /dble(ratio_to_coarser(2))
                        U_crse2(i_c,j_c,k_c) = U_crse2(i_c,j_c,k_c) +
     &                       w0*w1*w2*w_fac*U_fine2(i_f+i,j_f+j,k_f+k)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif

      if (physical_bdry_flags(2,1) .eq. 1) then ! upper z patch bdry is on physical bdry
         w_fac = 0.d0
         do k = -ratio_to_coarser(2)+1,0
            do j = 0,ratio_to_coarser(1)-1
               do i = 0,ratio_to_coarser(0)-1
                  w0 = 1.d0
                  w1 = 1.d0
                  w2 = 1.d0 - dabs(dble(k))/dble(ratio_to_coarser(2))
                  w_fac = w_fac + w0*w1*w2
               enddo
            enddo
         enddo
         w_fac = 1.d0/w_fac

         k_c = iupperc2
         do j_c = coarse_box_lower(1),coarse_box_upper(1)
            do i_c = coarse_box_lower(0),coarse_box_upper(0)
               i_f = i_c*ratio_to_coarser(0)
               j_f = j_c*ratio_to_coarser(1)
               k_f = k_c*ratio_to_coarser(2)
               U_crse2(i_c,j_c,k_c) = 0.d0
               do k = -ratio_to_coarser(2)+1,0
                  do j = 0,ratio_to_coarser(1)-1
                     do i = 0,ratio_to_coarser(0)-1
                        w0 = 1.d0
                        w1 = 1.d0
                        w2 = 1.d0 - dabs(dble(k))
     &                       /dble(ratio_to_coarser(2))
                        U_crse2(i_c,j_c,k_c) = U_crse2(i_c,j_c,k_c) +
     &                       w0*w1*w2*w_fac*U_fine2(i_f+i,j_f+j,k_f+k)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
