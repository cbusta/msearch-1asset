module lib_smalapack

! Module to ease the usage of certain LAPACK's subroutines
! By: Christian Bustamante
! Email: cbustaam@gmail.com
! Last modified: 27 Nov, 2019, 19:20
! v 1.10
!
! To compile this module:
! $ ifort -c Lib_Smalapack.f90 Lib_Kind.f90 -lblas -llapack
! $ gfortran -c Lib_Smalapack.f90 Lib_Kind.f90 -lblas -llapack
! The Lapack pachage in Ubuntu is [liblapack-dev]
! To use it put this just after the -program- statement:
! use lib_smalapck
!
! Subroutines:
! (Smalapack) 01. invlap        - Computes the inverse of a matrix using a LU factorization.
! (Smalapack) 02. solve_linsys  - Solves a system of linear equations Ax=b, for a general matrix A.
! (Smalapack) 03. lufactor      - Computes the LU-factorization of a general matrix N-by-M.
!
! Note: depends on Lib_Kind
! When using Intel Fortran in Visual Studio, this module requires to activate the MKL library in 
! Project -> Properties -> Fortran -> Libraries -> Use Intel Math Kernel Library -or-
! Manually: include (at compilation) /Qmkl:parallel -or- /Qmkl:sequential -or- /Qmkl:cluster, depending on the case.
! If explicit interfaces needed, include mkl.fi, lapack.f90

use lib_kind
implicit none
private
public :: invlap, solve_linsys, lufactor

contains


    function invlap(A) result(Ainv)
        ! Computes the inverse of a matrix using a L-U factorization.
        ! This function is more efficient than -inv- in the -lib_basic- module
        ! Usage:
        ! call Ainv = invlap(A)
        ! Inputs:
        ! A     - squared matrix to be inverted
        ! Output:
        ! Ainv  - inverse of A
        real(dp),intent(in)  :: A(:,:)
        real(dp)             :: Ainv(size(A,1),size(A,2))
        integer              :: n,info
        integer              :: ipiv(size(A,1))
        real(dp)             :: work(size(A,1))
        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)
        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges
        call dgetrf(n,n,Ainv,n,ipiv,info)
        if (info /= 0) then
            print *, 'Error (invlap): Matrix is numerically singular!'
        end if
        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call dgetri(n, Ainv, n, ipiv, work, n, info)
        if (info /= 0) then
            print *, 'Error (invlap): Matrix inversion failed!'
        end if
    end function invlap
  
  
    function solve_linsys(A,b) result(x)
        ! Solves a system of linear equations with an LU-factored
        ! square coefficient matrix.
        ! This is more efficient than x = matmul(invlap(A),b)
        ! Usage:
        ! x = solve_linsys(A,b)
        ! Inputs:
        ! A     - coefficient matrix (squared), n-by-n
        ! b     - n vector of constants, n
        ! Output:
        ! x     - vector with the solution to Ax=b
        integer,intent(in)   :: A(:,:),b(:)
        real(dp)             :: x(size(A,1))
        integer              :: na,nb,info
        integer              :: ipiv(size(A,1))
        real(dp)             :: Ainv(size(A,1),size(A,2))
        ! Check that the dimensions are right
        if (size(A,1) /= size(A,2)) then
            print *, 'Error (solve_linsys): Matrix A is not squared!'
        end if
        if (size(A,1) /= size(b,1)) then
            print *, 'Error (solve_linsys): Number of LHS and RHS equations do not match!'
        end if
        ! Store A in Ainv and b in x to prevent it from being overwritten
        ! by LAPACK
        Ainv = A
        x    = b
        na = size(A,1)
        nb = size(B,1)
        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges
        call dgetrf(na,na,Ainv,na,ipiv,info)
        if (info /= 0) then
            print *, 'Error (solve_linsys): Matrix is numerically singular!'
        end if
        ! DGESV solves a linear system of equations with a LU-factored
        ! matrix A
        call dgesv(na,1,Ainv,na,ipiv,x,nb,info)
    end function solve_linsys
  
  
    function lufactor(A) result(luA)
        ! Computes the LU factorization for a general matrix M-by-N
        ! Usage:
        ! luA = lufactor(A)
        ! Inputs:
        ! A     - matrix to be LU-factored
        ! Output:
        ! luA   - LU factorization of A
        integer,intent(in)   :: A(:,:)
        real(dp)             :: luA(size(A,1),size(A,2))
        integer              :: m,n,info
        integer              :: ipiv(min(size(A,1),size(A,2)))
        real(dp)             :: Ainv(size(A,1),size(A,2))
        ! Store A in luA to prevent it from being overwritten by LAPACK
        luA = A
        m = size(A,1)
        n = size(A,2)
        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges
        call dgetrf(m,n,Ainv,m,ipiv,info)
    end function lufactor
    
  
end module lib_smalapack