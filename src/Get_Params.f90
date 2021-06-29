module get_params

! This module retrieves model parameters and grids. This includes a header that declares
! contants and preallocate vectors, and two other subroutines that generate grids.
!
! More money for some
! Christian Bustamante
! Last modified: 15 May 2020

!-------------------------------------------------------------------------
! Retrieving Model Parameters and Grids
!-------------------------------------------------------------------------
use lib_kind
use lib_basic
implicit none

! Lengths
integer,parameter :: r   = 48
integer,parameter :: Nm  = r+2                ! Number of grid points for value fns and to solve bargaining problem
integer,parameter :: Ng  = 500                ! Number of grid points for the distributions

! Nash bargaining?
integer,parameter  :: dnash  = 0

! Parameter (matching and money)
real(dp),parameter :: theta  = 0.99           ! Buyer's bargaining power (when using Nash bargaining)
real(dp),parameter :: alpha  = 1.0            ! Probability of being matched
real(dp),parameter :: sigma  = 0.5            ! Probability of single coincidence (halved)
real(dp),parameter :: delta  = 0.0            ! Probability of double coincidence
real(dp),parameter :: pyear  = 1.0            ! Model periods per year
real(dp),parameter :: ipyear = 1.0/pyear      ! Inverse of model periods per year
real(dp),parameter :: rreal  = 1.04**ipyear   ! Real interest rate
real(dp),parameter :: beta   = 1.0/rreal      ! Discount factor
real(dp),parameter :: mum    = 1.02**ipyear-1 ! Rate for lump sum transfers of money

! Parameters for DM
real(dp),parameter :: eta    = 0.99           ! DM: consumption risk aversion parameter (kindof)
real(dp),parameter :: b_u    = 1e-4           ! DM: needed to guarantee u(c1=0) = 0
real(dp),parameter :: Bdm    = 4.0            ! DM: scale of disutility of labor (to target avg markup?)
real(dp),parameter :: nu     = 1.5            ! DM: curvature of disutility of labor

! Parameters for CM
integer,parameter  :: separa = 1              ! Is the utility function in CM addt. separable?
real(dp),parameter :: gamma  = 1.0            ! Risk aversion parameter
real(dp),parameter :: kappa  = 0.6            ! Scale parameter of utulity function (can go well above 1, ie 2.5 works, but not close to 0)
!real(dp),parameter :: chi    = 0.5*0          ! Inverse of Frisch elasticity of labor supply

! Grids for M
real(dp),parameter :: mlow   = 0.0001         ! Lowest possible level of real money holdings
real(dp),parameter :: mhigh  = 5.0            ! Highest possible level of real money holdings

! Tolerance levels
real(dp),parameter :: tol_pp = 1e-8           ! Tolerance for the bargaining problem
real(dp),parameter :: tol_cm = 1e-7           ! Tolerance for the CM solution (Golden)
real(dp),parameter :: tol_w  = tol_cm*10      ! Tolerance for the value fn solution
real(dp),parameter :: tol_mu = tol_w *10      ! Tolerance for the conv. of distributions
real(dp),parameter :: tol_ph = tol_mu*10      ! Tolerance for conv. of price of money
real(dp),parameter :: tol_nb = tol_ph*10      ! Tolerance for conv. of ToT

! Grids and matrices
real(dp)           :: KnotsM(Nm),Grid_M(Ng)


contains


    subroutine get_grids(KnotsM,Grid_M,Nm,Ng,mlow,mhigh)
        integer,intent(in)   :: Nm,Ng
        real(dp),intent(in)  :: mlow,mhigh
        real(dp),intent(out) :: KnotsM(Nm),Grid_M(Ng)
        ! Log-linearly spaced grid for M
        call logspace(mlow,mhigh,Nm,KnotsM)
        ! Linearly spaced grid for the distribution
        call logspace(mlow,mhigh,Ng,Grid_M)
    end subroutine get_grids


end module get_params