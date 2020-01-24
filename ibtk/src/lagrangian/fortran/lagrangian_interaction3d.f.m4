c ---------------------------------------------------------------------
c
c Copyright (c) 2011 - 2019 by the IBAMR developers
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
dnl Define some macros that help us unroll interpolation loops. These assume
dnl that we are accumulating into a 2D array V from a 4D array u which are both
dnl indexed first by depth (and then by point number s, which is specified in an
dnl outer loop). Assume d, ic0, ic1, and ic2 are available as inner loop indicies.
dnl
dnl The arguments are the lower and upper bounds in z, lower and upper bounds in
dnl y, and lower and upper bounds in x.
define(INTERPOLATE_INNER_3D,
          ` do d = 0,depth-1
               V(d, s) = 0.d0
               do ic2 = $1,$2
                  do ic1 = $3,$4
                     do ic0 = $5,$6
                        V(d,s) = V(d,s)
     &                       +w(0,ic0-$5)
     &                       *w(1,ic1-$3)
     &                       *w(2,ic2-$1)
     &                       *u(ic0,ic1,ic2,d)
                     enddo
                  enddo
               enddo
            enddo')
dnl Same arguments as before, but the seventh argument is the width of the
dnl stencil (e.g., 3 for bspline 3). The first branch is a hotter code path
dnl since when we are not at a boundary the number of inner loop iterations
dnl is known. Exposing this to the compiler helps generate code which speeds
dnl up the subroutine by about 25%.
define(INTERPOLATE_3D_SPECIALIZE_FIXED_WIDTH,
`   if ($2 - $1 == ($7 - 1) .and.
     &       $4 - $3 == ($7 - 1) .and.
     &       $6 - $5 == ($7 - 1)) then
           INTERPOLATE_INNER_3D($1, ($1 + $7 - 1),
                                $3, ($3 + $7 - 1),
                                $5, ($5 + $7 - 1))
         else
           INTERPOLATE_INNER_3D($1, $2,
                                $3, $4,
                                $5, $6)
         endif')dnl
define(SPREAD_INNER_3D,
          ` do d = 0,depth-1
               do ic2 = $1,$2
                  do ic1 = $3,$4
                     do ic0 = $5,$6
                        u(ic0,ic1,ic2,d) = u(ic0,ic1,ic2,d) + (
     &                       +w(0,ic0-$5)
     &                       *w(1,ic1-$3)
     &                       *w(2,ic2-$1)
     &                       *V(d,s)/(dx(0)*dx(1)*dx(2)))
                     enddo
                  enddo
               enddo
            enddo')
dnl Same arguments as before, but the seventh argument is the width of the
dnl stencil (e.g., 3 for bspline 3). The first branch is a hotter code path
dnl since when we are not at a boundary the number of inner loop iterations
dnl is known. Exposing this to the compiler helps generate code which speeds
dnl up the subroutine by about 25%.
define(SPREAD_3D_SPECIALIZE_FIXED_WIDTH,
`   if ($2 - $1 == ($7 - 1) .and.
     &       $4 - $3 == ($7 - 1) .and.
     &       $6 - $5 == ($7 - 1)) then
           SPREAD_INNER_3D($1, ($1 + $7 - 1),
                           $3, ($3 + $7 - 1),
                           $5, ($5 + $7 - 1))
         else
           SPREAD_INNER_3D($1, $2,
                           $3, $4,
                           $5, $6)
         endif')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl

c     this is a Fortran include, not an m4 include
      include 'lagrangian_delta.f'

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     piecewise constant delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_constant_interp3d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
      INTEGER d,l,s
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Use the piecewise constant delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         ic0 = NINT((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0)-0.5d0)+ilower0
         ic1 = NINT((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1)-0.5d0)+ilower1
         ic2 = NINT((X(2,s)+Xshift(2,l)-x_lower(2))/dx(2)-0.5d0)+ilower2
c
c     Interpolate u onto V.
c
         do d = 0,depth-1
            V(d,s) = u(ic0,ic1,ic2,d)
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the
c     piecewise constant delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_constant_spread3d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ic0,ic1,ic2
      INTEGER d,l,s
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Use the piecewise constant delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         ic0 = NINT((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0)-0.5d0)+ilower0
         ic1 = NINT((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1)-0.5d0)+ilower1
         ic2 = NINT((X(2,s)+Xshift(2,l)-x_lower(2))/dx(2)-0.5d0)+ilower2
c
c     Spread V onto u.
c
         do d = 0,depth-1
            u(ic0,ic1,ic2,d) = u(ic0,ic1,ic2,d)
     &           + V(d,s)/(dx(0)*dx(1)*dx(2))
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     discontinuous linear delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_discontinuous_linear_interp3d(
     &     dx,x_lower,x_upper,depth,axis,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Input.
c
      INTEGER depth,axis
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_cell(0:NDIM-1),X_shifted(0:NDIM-1),w(0:NDIM-1,0:1)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use the discontinuous linear delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Shift the point by X_shift.
c
         do d = 0,NDIM-1
            X_shifted(d) = X(d,s)+Xshift(d,l)
         enddo
c
c     Determine the interpolation stencils and weights.
c
         do d = 0,NDIM-1
            ic_center(d) = ilower(d) +
     &           NINT((X_shifted(d)-x_lower(d))/dx(d)-0.5d0)
            X_cell(d) = x_lower(d) +
     &           (dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)

            if ( d.eq.axis ) then
               if ( X_shifted(d).lt.X_cell(d) ) then
                  ic_lower(d) = ic_center(d)-1
                  ic_upper(d) = ic_center(d)
                  w(d,0) = (X_cell(d)-X_shifted(d))/dx(d)
                  w(d,1) = 1.d0 - w(d,0)
               else
                  ic_lower(d) = ic_center(d)
                  ic_upper(d) = ic_center(d)+1
                  w(d,0) = 1.d0 + (X_cell(d)-X_shifted(d))/dx(d)
                  w(d,1) = 1.d0 - w(d,0)
               endif
            else
               w(d,0) = 1.d0
               ic_lower(d) = ic_center(d)
               ic_upper(d) = ic_center(d)
            endif

            ic_lower(d) = max(ic_lower(d), ilower(d) - nugc(d))
            ic_upper(d) = max(ic_upper(d), iupper(d) + nugc(d))
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_3D_SPECIALIZE_FIXED_WIDTH(ic_lower(2), ic_upper(2),
                                               ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               2)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the
c     discontinuous linear delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_discontinuous_linear_spread3d(
     &     dx,x_lower,x_upper,depth,axis,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u)
c
      implicit none
c
c     Input.
c
      INTEGER depth,axis
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER ic_trimmed_lower(0:NDIM-1),ic_trimmed_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_cell(0:NDIM-1),X_shifted(0:NDIM-1),w(0:NDIM-1,0:1)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use the discontinuous linear delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Shift the point by X_shift.
c
         do d = 0,NDIM-1
            X_shifted(d) = X(d,s)+Xshift(d,l)
         enddo
c
c     Determine the interpolation stencils and weights.
c
         do d = 0,NDIM-1
            if ( d.eq.axis ) then
               ic_center(d) = ilower(d) +
     &              NINT((X_shifted(d)-x_lower(d))/dx(d)-0.5d0)
               X_cell(d) = x_lower(d) +
     &              (dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)

               if ( X_shifted(d).lt.X_cell(d) ) then
                  ic_lower(d) = ic_center(d)-1
                  ic_upper(d) = ic_center(d)
                  w(d,0) = (X_cell(d)-X_shifted(d))/dx(d)
                  w(d,1) = 1.d0 - w(d,0)
               else
                  ic_lower(d) = ic_center(d)
                  ic_upper(d) = ic_center(d)+1
                  w(d,0) = 1.d0 + (X_cell(d)-X_shifted(d))/dx(d)
                  w(d,1) = 1.d0 - w(d,0)
               endif
            else
               w(d,0) = 1.d0
               ic_lower(d) = ic_center(d)
               ic_upper(d) = ic_center(d)
            endif

            ic_trimmed_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_trimmed_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
         enddo
c
c     Spread V onto u.
c
         do d = 0,depth-1
            do ic2 = ic_trimmed_lower(2),ic_trimmed_upper(2)
               do ic1 = ic_trimmed_lower(1),ic_trimmed_upper(1)
                  do ic0 = ic_trimmed_lower(0),ic_trimmed_upper(0)
                     u(ic0,ic1,ic2,d) = u(ic0,ic1,ic2,d)+(
     &                    w(0,ic0-ic_lower(0))*
     &                    w(1,ic1-ic_lower(1))*
     &                    w(2,ic2-ic_lower(2))*
     &                    V(d,s)/(dx(0)*dx(1)*dx(2)))
                  enddo
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     piecewise linear delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_linear_interp3d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_cell(0:NDIM-1),X_shifted(0:NDIM-1),w(0:NDIM-1,0:1)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use the piecewise linear delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Shift the point by X_shift.
c
         do d = 0,NDIM-1
            X_shifted(d) = X(d,s)+Xshift(d,l)
         enddo
c
c     Determine the interpolation stencils and weights.
c
         do d = 0,NDIM-1
            ic_center(d) = ilower(d) +
     &           NINT((X_shifted(d)-x_lower(d))/dx(d)-0.5d0)
            X_cell(d) = x_lower(d) +
     &           (dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)

            if ( X_shifted(d).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)
               w(d,0) = (X_cell(d)-X_shifted(d))/dx(d)
               w(d,1) = 1.d0 - w(d,0)
            else
               ic_lower(d) = ic_center(d)
               ic_upper(d) = ic_center(d)+1
               w(d,0) = 1.d0 + (X_cell(d)-X_shifted(d))/dx(d)
               w(d,1) = 1.d0 - w(d,0)
            endif
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_3D_SPECIALIZE_FIXED_WIDTH(ic_lower(2), ic_upper(2),
                                               ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               2)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the
c     piecewise linear delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_linear_spread3d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_cell(0:NDIM-1),X_shifted(0:NDIM-1),w(0:NDIM-1,0:1)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use the piecewise linear delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Shift the point by X_shift.
c
         do d = 0,NDIM-1
            X_shifted(d) = X(d,s)+Xshift(d,l)
         enddo
c
c     Determine the interpolation stencils and weights.
c
         do d = 0,NDIM-1
            ic_center(d) = ilower(d) +
     &           NINT((X_shifted(d)-x_lower(d))/dx(d)-0.5d0)
            X_cell(d) = x_lower(d) +
     &           (dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)

            if ( X_shifted(d).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)
               w(d,0) = (X_cell(d)-X_shifted(d))/dx(d)
               w(d,1) = 1.d0 - w(d,0)
            else
               ic_lower(d) = ic_center(d)
               ic_upper(d) = ic_center(d)+1
               w(d,0) = 1.d0 + (X_cell(d)-X_shifted(d))/dx(d)
               w(d,1) = 1.d0 - w(d,0)
            endif

            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
         enddo
c
c     Spread V onto u.
c
         SPREAD_3D_SPECIALIZE_FIXED_WIDTH(
             ic_lower(2), ic_upper(2),
             ic_lower(1), ic_upper(1),
             ic_lower(0), ic_upper(0),
             2)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     piecewise cubic delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_cubic_interp3d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_piecewise_cubic_delta
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_cell(0:NDIM-1),w(0:NDIM-1,0:3)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use the piecewise cubic delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell.
c
         do d = 0,NDIM-1
c
c     Determine the Cartesian cell in which X(s) is located.
c
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d) +
     &           (dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-2
               ic_upper(d) = ic_center(d)+1
            else
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)+2
            endif
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c     
c     Compute the interpolation weights.
c     
            do ic0 = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic0-ilower(d))+0.5d0)*dx(d)
               w(d,ic0-ic_lower(d)) =
     &              lagrangian_piecewise_cubic_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_3D_SPECIALIZE_FIXED_WIDTH(ic_lower(2), ic_upper(2),
                                               ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               4)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the
c     piecewise cubic delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_cubic_spread3d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_piecewise_cubic_delta
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_cell(0:NDIM-1),w(0:NDIM-1,0:3)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use the piecewise cubic delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c     
c     Determine the Cartesian cell in which X(s) is located.
c     
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d) +
     &           (dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
c     
c     Determine the spreading stencil corresponding to the position of
c     X(s) within the cell.
c     
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-2
               ic_upper(d) = ic_center(d)+1
            else
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)+2
            endif
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c     
c     Compute the spreading weights.
c     
            do ic0 = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic0-ilower(d))+0.5d0)*dx(d)
               w(d,ic0-ic_lower(d)) =
     &              lagrangian_piecewise_cubic_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Spread V onto u.
c
         SPREAD_3D_SPECIALIZE_FIXED_WIDTH(
             ic_lower(2), ic_upper(2),
             ic_lower(1), ic_upper(1),
             ic_lower(0), ic_upper(0),
             4)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the IB
c     3-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_3_interp3d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_ib_3_delta
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_cell(0:NDIM-1),w(0:NDIM-1,0:2)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use the IB 3-point delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d = 0,NDIM-1
c     
c     Determine the Cartesian cell in which X(s) is located.
c     
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d) +
     &           (dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
c     
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell.
c     
            ic_lower(d) = ic_center(d)-1
            ic_upper(d) = ic_center(d)+1
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c     
c     Compute the interpolation weights.
c     
            do ic0 = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic0-ilower(d))+0.5d0)*dx(d)
               w(d,ic0-ic_lower(d)) =
     &              lagrangian_ib_3_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_3D_SPECIALIZE_FIXED_WIDTH(ic_lower(2), ic_upper(2),
                                               ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               3)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the IB
c     3-point delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_3_spread3d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_ib_3_delta
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_cell(0:NDIM-1),w(0:NDIM-1,0:2)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use the IB 3-point delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d = 0,NDIM-1
c     
c     Determine the Cartesian cell in which X(s) is located.
c     
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
c     
c     Determine the spreading stencil corresponding to the position of
c     X(s) within the cell.
c     
            ic_lower(d) = ic_center(d)-1
            ic_upper(d) = ic_center(d)+1
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c     
c     Compute the spreading weights.
c     
            do ic0 = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)
     &              +(dble(ic0-ilower(d))+0.5d0)*dx(d)
               w(d,ic0-ic_lower(d)) =
     &              lagrangian_ib_3_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Spread V onto u.
c
         SPREAD_3D_SPECIALIZE_FIXED_WIDTH(
             ic_lower(2), ic_upper(2),
             ic_lower(1), ic_upper(1),
             ic_lower(0), ic_upper(0),
             3)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the IB
c     4-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_4_interp3d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_o_dx,q,r
      REAL w(0:NDIM-1,0:3)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use the IB 4-point delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell and compute the interpolation weights.
c
         do d = 0,NDIM-1
            X_o_dx = (X(d,s)+Xshift(d,l)-x_lower(d))/dx(d)
            ic_lower(d) = NINT(X_o_dx)+ilower(d)-2
            ic_upper(d) = ic_lower(d) + 3
            ic_lower(d) = max(ic_lower(d), ilower(d) - nugc(d));
            ic_upper(d) = min(ic_upper(d), iupper(d) + nugc(d));

            r = X_o_dx - ((ic_lower(d)+1-ilower(d))+0.5d0)
            q = sqrt(1.d0+4.d0*r*(1.d0-r))
            w(d, 0) = 0.125d0*(3.d0-2.d0*r-q)
            w(d, 1) = 0.125d0*(3.d0-2.d0*r+q)
            w(d, 2) = 0.125d0*(1.d0+2.d0*r+q)
            w(d, 3) = 0.125d0*(1.d0+2.d0*r-q)
         enddo

c     
c     Interpolate u onto V.
c     
         INTERPOLATE_3D_SPECIALIZE_FIXED_WIDTH(ic_lower(2), ic_upper(2),
                                               ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               4)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the IB
c     4-point delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_4_spread3d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_o_dx,q,r
      REAL w(0:NDIM-1,0:3)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use the IB 4-point delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell and compute the interpolation weights.
c
         do d=0,NDIM-1
            X_o_dx = (X(d,s)+Xshift(d,l)-x_lower(d))/dx(d)
            ic_lower(d) = NINT(X_o_dx)+ilower(d)-2
            ic_upper(d) = ic_lower(d) + 3
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))

            r = X_o_dx - ((ic_lower(d)+1-ilower(d))+0.5d0)
            q = sqrt(1.d0+4.d0*r*(1.d0-r))
            w(d, 0) = 0.125d0*(3.d0-2.d0*r-q)
            w(d, 1) = 0.125d0*(3.d0-2.d0*r+q)
            w(d, 2) = 0.125d0*(1.d0+2.d0*r+q)
            w(d, 3) = 0.125d0*(1.d0+2.d0*r-q)
         enddo

c
c     Spread V onto u.
c
         SPREAD_3D_SPECIALIZE_FIXED_WIDTH(ic_lower(2), ic_upper(2),
                                          ic_lower(1), ic_upper(1),
                                          ic_lower(0), ic_upper(0),
                                          4)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the IB
c     4-point delta function that has been broadened to have a support
c     of 8 meshwidths.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_4_w8_interp3d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_o_dx,q,r
      REAL w(0:NDIM-1,0:7)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use a broadened version of the IB 4-point delta function to
c     interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)

         do d=0,NDIM-1
c     
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell and compute the interpolation weights.
c     
            X_o_dx = (X(d,s)+Xshift(d,l)-x_lower(d))/dx(d)
            ic_lower(d) = NINT(X_o_dx)+ilower(d)-4
            ic_upper(d) = ic_lower(d) + 7
            ic_lower(d) = max(ic_lower(d), ilower(d) - nugc(d))
            ic_upper(d) = min(ic_upper(d), iupper(d) + nugc(d))

            r = 0.5d0*(X_o_dx - ((ic_lower(d)+3-ilower(d))+0.5d0))
            q = sqrt(1.d0+4.d0*r*(1.d0-r))
            w(d, 1) = 0.0625d0*(3.d0-2.d0*r-q)
            w(d, 3) = 0.0625d0*(3.d0-2.d0*r+q)
            w(d, 5) = 0.0625d0*(1.d0+2.d0*r+q)
            w(d, 7) = 0.0625d0*(1.d0+2.d0*r-q)
            r = r+0.5d0
            q = sqrt(1.d0+4.d0*r*(1.d0-r))
            w(d, 0) = 0.0625d0*(3.d0-2.d0*r-q)
            w(d, 2) = 0.0625d0*(3.d0-2.d0*r+q)
            w(d, 4) = 0.0625d0*(1.d0+2.d0*r+q)
            w(d, 6) = 0.0625d0*(1.d0+2.d0*r-q)
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_3D_SPECIALIZE_FIXED_WIDTH(ic_lower(2),ic_upper(2),
                                               ic_lower(1),ic_upper(1),
                                               ic_lower(0),ic_upper(0),
                                               8)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the IB
c     4-point delta function that has been broadened to have a support
c     of 8 meshwidths.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_4_w8_spread3d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_o_dx,q,r
      REAL w(0:NDIM-1,0:7)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use a broadened version of the IB 4-point delta function to spread
c     V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell and compute the interpolation weights.
c
            X_o_dx = (X(d,s)+Xshift(d,l)-x_lower(d))/dx(d)
            ic_lower(d) = NINT(X_o_dx)+ilower(d)-4
            ic_upper(d) = ic_lower(d) + 7
            ic_lower(d) = max(ic_lower(d), ilower(d) - nugc(d))
            ic_upper(d) = min(ic_upper(d), iupper(d) + nugc(d))

            r = 0.5d0*(X_o_dx - ((ic_lower(d)+3-ilower(d))+0.5d0))
            q = sqrt(1.d0+4.d0*r*(1.d0-r))
            w(d, 1) = 0.0625d0*(3.d0-2.d0*r-q)
            w(d, 3) = 0.0625d0*(3.d0-2.d0*r+q)
            w(d, 5) = 0.0625d0*(1.d0+2.d0*r+q)
            w(d, 7) = 0.0625d0*(1.d0+2.d0*r-q)
            r = r+0.5d0
            q = sqrt(1.d0+4.d0*r*(1.d0-r))
            w(d, 0) = 0.0625d0*(3.d0-2.d0*r-q)
            w(d, 2) = 0.0625d0*(3.d0-2.d0*r+q)
            w(d, 4) = 0.0625d0*(1.d0+2.d0*r+q)
            w(d, 6) = 0.0625d0*(1.d0+2.d0*r-q)
         enddo
c
c     Spread V onto u.
c
         SPREAD_3D_SPECIALIZE_FIXED_WIDTH(ic_lower(2),ic_upper(2),
                                          ic_lower(1),ic_upper(1),
                                          ic_lower(0),ic_upper(0),
                                          8)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     5-point IB delta with three continuous derivatives
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_5_interp3d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Functions.
c
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)

c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      REAL r
      REAL phi
      REAL K
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_cell(0:NDIM-1),w(0:NDIM-1,0:4)

      PARAMETER (K = (38.0d0 - sqrt(69.0d0))/60.0d0)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use a 5-point IB delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c     
c     Determine the Cartesian cell in which X(s) is located.
c     
            ic_center(d) = floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d) +
     &           (dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
c     
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell.
c     
            ic_lower(d) = ic_center(d)-2
            ic_upper(d) = ic_center(d)+2
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c     
c     Compute the interpolation weights.
c     
            ic0 = ic_center(d)
            X_cell(d) = x_lower(d)+(dble(ic0-ilower(d))+0.5d0)*dx(d)
            r = (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d)
            phi = (136.0d0 - 40.0d0*K - 40.0d0*r**2 + sqrt(2.0d0)
     &           *sqrt(3123.0d0 - 6840.0d0*K + 3600.0d0*(K**2)
     &           - 12440.0d0*(r**2) + 25680.0d0*K*(r**2)
     &           - 12600.0d0*(K**2)*(r**2) + 8080.0d0*(r**4)
     &           - 8400.0d0*K*(r**4) - 1400.0d0*(r**6)))
     &           /280.0d0

            w(d,0) = (1.0d0/12.0d0) * (-2.0d0 + 2.0d0*phi + 2.0d0*K +
     &           r - 3.0d0*K*r + 2.0d0*r**2 - r**3)
            w(d,1) = (1.0d0/ 6.0d0) * ( 4.0d0 - 4.0d0*phi -       K -
     &           4.0d0*r + 3.0d0*K*r -       r**2 + r**3)
            w(d,2) = phi
            w(d,3) = (1.0d0/ 6.0d0) * ( 4.0d0 - 4.0d0*phi -       K +
     &           4.0d0*r - 3.0d0*K*r -       r**2 - r**3)
            w(d,4) = (1.0d0/12.0d0) * (-2.0d0 + 2.0d0*phi + 2.0d0*K -
     &           r + 3.0d0*K*r + 2.0d0*r**2 + r**3)
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_3D_SPECIALIZE_FIXED_WIDTH(ic_lower(2), ic_upper(2),
                                               ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               5)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the
c     5-point IB delta with three continuous derivatives
c     using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_5_spread3d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u)
c
      implicit none
c
c     Functions.
c
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      REAL r
      REAL phi
      REAL K
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_cell(0:NDIM-1),w(0:NDIM-1,0:4)

      PARAMETER (K = (38.0d0 - sqrt(69.0d0))/60.0d0)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use a 5-point IB delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c     
c     Determine the Cartesian cell in which X(s) is located.
c     
            ic_center(d) = floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d)+(dble(ic_center(d)-ilower(d))+0.5d0)
     &           * dx(d)
c     
c     Determine the spreading stencil corresponding to the position of
c     X(s) within the cell.
c     
            ic_lower(d) = ic_center(d)-2
            ic_upper(d) = ic_center(d)+2
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c     
c     Compute the spreading weights.
c     
            ic0 = ic_center(d)
            X_cell(d) = x_lower(d)+(dble(ic0-ilower(d))+0.5d0)*dx(d)
            r = (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d)
            phi = (136.0d0 - 40.0d0*K - 40.0d0*r**2 + sqrt(2.0d0)
     &           *sqrt(3123.0d0 - 6840.0d0*K + 3600.0d0*(K**2)
     &           - 12440.0d0*(r**2) + 25680.0d0*K*(r**2)
     &           - 12600.0d0*(K**2)*(r**2) + 8080.0d0*(r**4)
     &           - 8400.0d0*K*(r**4) - 1400.0d0*(r**6)))
     &           /280.0d0

            w(d,0) = (1.0d0/12.0d0) * (-2.0d0 + 2.0d0*phi + 2.0d0*K +
     &           r - 3.0d0*K*r + 2.0d0*r**2 - r**3)
            w(d,1) = (1.0d0/ 6.0d0) * ( 4.0d0 - 4.0d0*phi -       K -
     &           4.0d0*r + 3.0d0*K*r -       r**2 + r**3)
            w(d,2) = phi
            w(d,3) = (1.0d0/ 6.0d0) * ( 4.0d0 - 4.0d0*phi -       K +
     &           4.0d0*r - 3.0d0*K*r -       r**2 - r**3)
            w(d,4) = (1.0d0/12.0d0) * (-2.0d0 + 2.0d0*phi + 2.0d0*K -
     &           r + 3.0d0*K*r + 2.0d0*r**2 + r**3)

         enddo

c
c     Spread V onto u.
c
         SPREAD_3D_SPECIALIZE_FIXED_WIDTH(
             ic_lower(2), ic_upper(2),
             ic_lower(1), ic_upper(1),
             ic_lower(0), ic_upper(0),
             5)

c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the IB
c     6-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_6_interp3d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_o_dx,r,alpha,beta,gamma,discr,K
      REAL pm3,pm2,pm1,p,pp1,pp2
      REAL w(0:NDIM-1,0:5)

      PARAMETER (K = (59.d0/60.d0)*(1.d0-sqrt(1.d0-(3220.d0/3481.d0))))
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use the IB 6-point delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell and compute the interpolation weights.
c
         do d=0,NDIM-1
            X_o_dx = (X(d,s)+Xshift(d,l)-x_lower(d))/dx(d)
            ic_lower(d) = NINT(X_o_dx)+ilower(d)-3
            ic_upper(d) = ic_lower(d) + 5
            ic_lower(d) = max(ic_lower(d), ilower(d) - nugc(d))
            ic_upper(d) = min(ic_upper(d), iupper(d) + nugc(d))
            r = 1.d0 - X_o_dx + ((ic_lower(d)+2-ilower(d))+0.5d0)

            alpha = 28.d0
            beta = (9.d0/4.d0)-(3.d0/2.d0)*(K+r**2)
     $           +((22.d0/3.d0)-7.d0*K)*r-(7.d0/3.d0)*r**3
            gamma = (1.d0/4.d0)*(((161.d0/36.d0)
     $           -(59.d0/6.d0)*K+5.d0*K**2)*(1.d0/2.d0)*r**2
     $           + (-(109.d0/24.d0)+5.d0*K)*(1.d0/3.d0)*r**4
     $           + (5.d0/18.d0)*r**6 )
            discr = beta**2-4.d0*alpha*gamma

            pm3 = (-beta+sign(1.d0,(3.d0/2.d0)-K)*sqrt(discr))
     $           /(2.d0*alpha)
            pm2 =  -3.d0*pm3 - (1.d0/16.d0) + (1.d0/8.d0)*(K+r**2)
     $            + (1.d0/12.d0)*(3.d0*K-1.d0)*r + (1.d0/12.d0)*r**3
            pm1 =   2.d0*pm3 + (1.d0/4.d0) +  (1.d0/6.d0)*(4.d0-3.d0*K)
     $           *r - (1.d0/6.d0)*r**3
            p   =   2.d0*pm3 + (5.d0/8.d0)  - (1.d0/4.d0)*(K+r**2)
            pp1 =  -3.d0*pm3 + (1.d0/4.d0) -  (1.d0/6.d0)*(4.d0-3.d0*K)
     $           *r + (1.d0/6.d0)*r**3
            pp2 =        pm3 - (1.d0/16.d0) + (1.d0/8.d0)*(K+r**2)
     $            - (1.d0/12.d0)*(3.d0*K-1.d0)*r - (1.d0/12.d0)*r**3

            w(d, 0) = pm3
            w(d, 1) = pm2
            w(d, 2) = pm1
            w(d, 3) = p
            w(d, 4) = pp1
            w(d, 5) = pp2
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_3D_SPECIALIZE_FIXED_WIDTH(ic_lower(2), ic_upper(2),
                                               ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               6)

c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the IB
c     6-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_6_spread3d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u)
c
      implicit none
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_o_dx,r,alpha,beta,gamma,discr,K
      REAL pm3,pm2,pm1,p,pp1,pp2
      REAL w(0:NDIM-1,0:5)

      PARAMETER (K = (59.d0/60.d0)*(1.d0-sqrt(1.d0-(3220.d0/3481.d0))))
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use a broadened version of the IB 4-point delta function to spread
c     V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell and compute the interpolation weights.
c
         do d=0,NDIM-1
            X_o_dx = (X(d,s)+Xshift(d,l)-x_lower(d))/dx(d)
            ic_lower(d) = NINT(X_o_dx)+ilower(d)-3
            ic_upper(d) = ic_lower(d) + 5
            ic_lower(d) = max(ic_lower(d), ilower(d) - nugc(d))
            ic_upper(d) = min(ic_upper(d), iupper(d) + nugc(d))
            r = 1.d0 - X_o_dx + ((ic_lower(d)+2-ilower(d))+0.5d0)

            alpha = 28.d0
            beta = (9.d0/4.d0)-(3.d0/2.d0)*(K+r**2)
     $           +((22.d0/3.d0)-7.d0*K)*r-(7.d0/3.d0)*r**3
            gamma = (1.d0/4.d0)*( ((161.d0/36.d0)-(59.d0/6.d0)*K
     $           +5.d0*K**2)*(1.d0/2.d0)*r**2 + (-(109.d0/24.d0)+5.d0*K)
     $           *(1.d0/3.d0)*r**4+ (5.d0/18.d0)*r**6 )
            discr = beta**2-4.d0*alpha*gamma

            pm3 = (-beta+sign(1.d0,(3.d0/2.d0)-K)*sqrt(discr))
     $           /(2.d0*alpha)
            pm2 =  -3.d0*pm3 - (1.d0/16.d0) + (1.d0/8.d0)*(K+r**2) 
     $           + (1.d0/12.d0)*(3.d0*K-1.d0)*r + (1.d0/12.d0)*r**3
            pm1 =   2.d0*pm3 + (1.d0/4.d0) +  (1.d0/6.d0)*(4.d0-3.d0*K)
     $           *r -(1.d0/6.d0)*r**3
            p   =   2.d0*pm3 + (5.d0/8.d0)  - (1.d0/4.d0)*(K+r**2)
            pp1 =  -3.d0*pm3 + (1.d0/4.d0) -  (1.d0/6.d0)*(4.d0-3.d0*K)
     $           *r +(1.d0/6.d0)*r**3
            pp2 =        pm3 - (1.d0/16.d0) + (1.d0/8.d0)*(K+r**2) 
     $           - (1.d0/12.d0)*(3.d0*K-1.d0)*r - (1.d0/12.d0)*r**3

            w(d, 0) = pm3
            w(d, 1) = pm2
            w(d, 2) = pm1
            w(d, 3) = p
            w(d, 4) = pp1
            w(d, 5) = pp2
         enddo
c
c     Spread V onto u.
c
         SPREAD_3D_SPECIALIZE_FIXED_WIDTH(ic_lower(2), ic_upper(2),
                                          ic_lower(1), ic_upper(1),
                                          ic_lower(0), ic_upper(0),
                                          6)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     3-point B-spline delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_3_interp3d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_3_delta
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_cell(0:NDIM-1),w(0:NDIM-1,0:2)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use a 3-point B-spline function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c     
c     Determine the Cartesian cell in which X(s) is located.
c     
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
c     
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell.
c     
            ic_lower(d) = ic_center(d)-1
            ic_upper(d) = ic_center(d)+1
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c     
c     Compute the interpolation weights.
c     
            do ic0 = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic0-ilower(d))+0.5d0)*dx(d)
               w(d,ic0-ic_lower(d)) =
     &              lagrangian_bspline_3_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_3D_SPECIALIZE_FIXED_WIDTH(ic_lower(2), ic_upper(2),
                                               ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               3)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the 3-point
c     B-spline delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_3_spread3d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_3_delta
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_cell(0:NDIM-1),w(0:NDIM-1,0:2)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use a 3-point B-spline function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c     
c     Determine the Cartesian cell in which X(s) is located.
c     
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
c     
c     Determine the spreading stencil corresponding to the position of
c     X(s) within the cell.
c     
            ic_lower(d) = ic_center(d)-1
            ic_upper(d) = ic_center(d)+1
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c     
c     Compute the spreading weights.
c     
            do ic0 = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic0-ilower(d))+0.5d0)*dx(d)
               w(d,ic0-ic_lower(d)) =
     &              lagrangian_bspline_3_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Spread V onto u.
c
         SPREAD_3D_SPECIALIZE_FIXED_WIDTH(ic_lower(2), ic_upper(2),
                                          ic_lower(1), ic_upper(1),
                                          ic_lower(0), ic_upper(0),
                                          3)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     4-point B-spline delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_4_interp3d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_4_delta
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_cell(0:NDIM-1),w(0:NDIM-1,0:3)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use a 4-point B-spline function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c     
c     Determine the Cartesian cell in which X(s) is located.
c     
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d) +
     &           (dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
c     
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell.
c     
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-2
               ic_upper(d) = ic_center(d)+1
            else
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)+2
            endif


            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c     
c     Compute the interpolation weights.
c     
            do ic0 = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic0-ilower(d))+0.5d0)*dx(d)
               w(d,ic0-ic_lower(d)) =
     &              lagrangian_bspline_4_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_3D_SPECIALIZE_FIXED_WIDTH(ic_lower(2), ic_upper(2),
                                               ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               4)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the 4-point
c     B-spline delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_4_spread3d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_4_delta
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_cell(0:NDIM-1),w(0:NDIM-1,0:3)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use a 4-point B-spline function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c     
c     Determine the Cartesian cell in which X(s) is located.
c     
            ic_center(d) = floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d) +
     &           (dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
c     
c     Determine the spreading stencil corresponding to the position of
c     X(s) within the cell.
c     
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-2
               ic_upper(d) = ic_center(d)+1
            else
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)+2
            endif

            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c     
c     Compute the spreading weights.
c     
            do ic0 = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic0-ilower(d))+0.5d0)*dx(d)
               w(d,ic0-ic_lower(d)) =
     &              lagrangian_bspline_4_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Spread V onto u.
c
         SPREAD_3D_SPECIALIZE_FIXED_WIDTH(
             ic_lower(2), ic_upper(2),
             ic_lower(1), ic_upper(1),
             ic_lower(0), ic_upper(0),
             4)

c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     5-point B-spline delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_5_interp3d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_5_delta
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_cell(0:NDIM-1),w(0:NDIM-1,0:4)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use a 5-point B-spline function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         do d=0,NDIM-1
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
c     
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell.
c     
            ic_lower(d) = ic_center(d)-2
            ic_upper(d) = ic_center(d)+2

            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c     
c     Compute the interpolation weights.
c     
            do ic0 = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic0-ilower(d))+0.5d0)*dx(d)
               w(d,ic0-ic_lower(d)) =
     &              lagrangian_bspline_5_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_3D_SPECIALIZE_FIXED_WIDTH(ic_lower(2), ic_upper(2),
                                               ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               5)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the 5-point
c     B-spline delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_5_spread3d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_5_delta
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_cell(0:NDIM-1),w(0:NDIM-1,0:4)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use a 5-point B-spline function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c     
c     Determine the Cartesian cell in which X(s) is located.
c     
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
c     
c     Determine the spreading stencil corresponding to the position of
c     X(s) within the cell.
c     
            ic_lower(d) = ic_center(d)-2
            ic_upper(d) = ic_center(d)+2
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c     
c     Compute the spreading weights.
c     
            do ic0 = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic0-ilower(d))+0.5d0)*dx(d)
               w(d,ic0-ic_lower(d)) =
     &              lagrangian_bspline_5_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Spread V onto u.
c
         SPREAD_3D_SPECIALIZE_FIXED_WIDTH(
             ic_lower(2), ic_upper(2),
             ic_lower(1), ic_upper(1),
             ic_lower(0), ic_upper(0),
             5)

c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     6-point B-spline delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_6_interp3d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_6_delta
c
c     Input.
c
      INTEGER depth
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2
      INTEGER nindices

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_cell(0:NDIM-1),w(0:NDIM-1,0:5)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use a 6-point B-spline function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,NDIM-1
c     
c     Determine the Cartesian cell in which X(s) is located.
c     
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d) +
     &           (dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
c     
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell.
c     
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-3
               ic_upper(d) = ic_center(d)+2
            else
               ic_lower(d) = ic_center(d)-2
               ic_upper(d) = ic_center(d)+3
            endif

            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c     
c     Compute the interpolation weights.
c     
            do ic0 = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic0-ilower(d))+0.5d0)*dx(d)
               w(d,ic0-ic_lower(d)) =
     &              lagrangian_bspline_6_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Interpolate u onto V.
c
         INTERPOLATE_3D_SPECIALIZE_FIXED_WIDTH(ic_lower(2), ic_upper(2),
                                               ic_lower(1), ic_upper(1),
                                               ic_lower(0), ic_upper(0),
                                               6)
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the 6-point
c     B-spline delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_6_spread3d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,ilower2,iupper2,
     &     nugc0,nugc1,nugc2,
     &     u)
c
      implicit none
c
c     Functions.
c
      REAL lagrangian_bspline_6_delta
c
c     Input.
c
      INTEGER depth
      INTEGER nindices
      INTEGER ilower0,iupper0,ilower1,iupper1,ilower2,iupper2
      INTEGER nugc0,nugc1,nugc2

      INTEGER indices(0:nindices-1)

      REAL Xshift(0:NDIM-1,0:nindices-1)

      REAL dx(0:NDIM-1),x_lower(0:NDIM-1),x_upper(0:NDIM-1)
      REAL u(CELL3dVECG(ilower,iupper,nugc),0:depth-1)
      REAL X(0:NDIM-1,0:*)
c
c     Input/Output.
c
      REAL V(0:depth-1,0:*)
c
c     Local variables.
c
      INTEGER ilower(0:NDIM-1),iupper(0:NDIM-1)
      INTEGER ic0,ic1,ic2
      INTEGER ic_center(0:NDIM-1),ic_lower(0:NDIM-1),ic_upper(0:NDIM-1)
      INTEGER d,l,s,nugc(0:NDIM-1)

      REAL X_cell(0:NDIM-1),w(0:NDIM-1,0:5)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1
      ilower(2) = ilower2

      iupper(0) = iupper0
      iupper(1) = iupper1
      iupper(2) = iupper2

      nugc(0) = nugc0
      nugc(1) = nugc1
      nugc(2) = nugc2
c
c     Use a 6-point B-spline function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         do d = 0,NDIM-1
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d) +
     &           (dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
c     
c     Determine the spreading stencil corresponding to the position of
c     X(s) within the cell.
c     
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-3
               ic_upper(d) = ic_center(d)+2
            else
               ic_lower(d) = ic_center(d)-2
               ic_upper(d) = ic_center(d)+3
            endif

            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c     
c     Compute the spreading weights.
c     
            do ic0 = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic0-ilower(d))+0.5d0)*dx(d)
               w(d,ic0-ic_lower(d)) =
     &              lagrangian_bspline_6_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Spread V onto u.
c
         SPREAD_3D_SPECIALIZE_FIXED_WIDTH(ic_lower(2), ic_upper(2),
                                          ic_lower(1), ic_upper(1),
                                          ic_lower(0), ic_upper(0),
                                          6)

c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
