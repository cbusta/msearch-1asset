module lib_integra

! Module with subroutines for numerical integration
! By: Christian Bustamante
! Email: cbustaam@gmail.com
! Last modified: 27 Nov, 2019, 19:11
! v 1.10
!
! To compile this module:
! $ ifort -c Lib_Integra.f90 Lib_Kind.f90
! $ gfortran -c Lib_Integra.f90 Lib_Kind.f90
! To use it put this just after the -program- statement:
! use lib_integra
!
! Subroutines:
! (Integra) 01. gauleg     - Numerical integration using Gauss-Legendre quadrature
!
! Note: depends on Lib_Kind.f90

use lib_kind
use lib_basic
use lib_smalapack
implicit none
private
public :: gauleg

contains


    subroutine gauleg(n,a,b, x, w)
        ! Numerical integration with Gauss-Legendre quarature (finite intervals)
        ! Adapted from "Numerical Recipes in Fortran", pp. 145
        ! Compute x(i) and w(i), i = 1,n  Legendre ordinates and weights on interval
        ! -1.0 to 1.0 (length is 2.0). Use ordinates and weights for Gauss Legendre integration
        ! Usage:
        ! call gauleg(n,a,b, x,w)
        ! Inputs:
        ! n      - number of points for quadrature
        ! a,b    - upper and lower bounds for integration
        ! Output:
        ! x     - abscissas for Gauss-Legendre quarature
        ! w     - weights for Gauss-Legendre quarature
        integer,intent(in)   :: n
        real(dp),intent(in)  :: a,b
        real(dp),intent(out) :: x(n), w(n)
        integer              :: i, j, m
        real(dp)             :: p1, p2, p3, pp, xl, xm, z, z1
        real(dp), parameter  :: eps = 3.d-14
        m  = (n+1)/2
        xm = 0.5d0*(b+a)
        xl = 0.5d0*(b-a)
        do i = 1,m
            z  = cos(3.141592654d0*(i-0.25d0)/(n+0.5d0))
            z1 = 0.0
            do while (abs(z-z1) > eps)
                p1 = 1.0d0
                p2 = 0.0d0
                do j = 1,n
                    p3 = p2
                    p2 = p1
                    p1 = ((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
                end do
                pp = n*(z*p1-p2)/(z*z-1.0d0)
                z1 = z
                z  = z1 - p1/pp
            end do
            x(i)     = xm - xl*z
            x(n+1-i) = xm + xl*z
            w(i)     = (2.0d0*xl)/((1.0d0-z*z)*pp*pp)
            w(n+1-i) = w(i)
        end do
    end subroutine gauleg
 
 
end module lib_integra