module bargaining
! This module contains the subroutines needed to solve the "take it or leave it problem",
! i.e., Nash bargaining with bargaining power of 0 for the seller.
! This subroutines includes one that calls an optimizer for each point in the grid
! (mb,ms), and the one that evaluates the objective function for provided (d,q)
!
! More money for some
! Christian Bustamante
! Last modified: 13 May 2020

use lib_kind
use lib_basic
use par_pass
use lib_solver
use lib_interp
use tools_fnutil
implicit none

contains


    subroutine solve_takeit(KnotsM,Coeffw,toln,PARAM_PP,Qsc,Dsc,Obj)
        ! --------------------------------------------------------------------------
        ! Solving the price posting problem
        ! Agent i = buyer, j = seller
        ! --------------------------------------------------------------------------
        real(dp), intent(in) :: KnotsM(:),Coeffw(:,:)
        real(dp), intent(in) :: toln
        real(dp), intent(out):: Qsc(size(KnotsM),size(KnotsM))
        real(dp), intent(out):: Dsc(size(KnotsM),size(KnotsM))
        real(dp),intent(out) :: Obj(size(KnotsM),size(KnotsM))
        integer              :: i,j,Nm,order,r
        real(dp)             :: d_star,q_star,f_star,d_lb,d_ub
        real(dp)             :: mb,ms,ms_end,Wsp_App,Ws0_App
        real(dp)             :: Bdm,nu
        real(dp)             :: muv
        type(parpass_g), intent(inout):: PARAM_PP

        order  = 2
        Nm     = size(KnotsM)
        r      = Nm-2
        Bdm    = PARAM_PP%Bdm_p
        nu     = PARAM_PP%nu_p
        muv    = PARAM_PP%mu_p
        !$omp parallel do private(i,j,mb,ms,d_lb,d_ub,d_star,f_star,ms_end,Wsp_App,Ws0_App,q_star) firstprivate(PARAM_PP)
        do i = 1,Nm
            do j = 1,Nm
                mb = KnotsM(i)
                ms = KnotsM(j)
                if (mb == 0.0) then
                    Qsc(i,j) = 0.0
                    Dsc(i,j) = 0.0
                    Obj(i,j) = 0.0
                else
                    PARAM_PP%mb_p = mb
                    PARAM_PP%ms_p = ms
                    d_lb = 0.0
                    d_ub = mb
                    call golden(d_lb,d_ub,toln,obj_takeit,PARAM_PP, d_star,f_star)
                    ms_end = ms + d_star
                    call eval_spline(Coeffw,KnotsM,r,muv + ms_end,order, Wsp_App)
                    call eval_spline(Coeffw,KnotsM,r,muv + ms,    order, Ws0_App)
                    q_star = ((1.0_dp/Bdm)*(Wsp_App-Ws0_App))**(1.0_dp/nu)
                    Qsc(i,j) = q_star
                    Dsc(i,j) = d_star
                    Obj(i,j) =-f_star
                end if
            end do
        end do
        !$omp end parallel do
    end subroutine solve_takeit




    ! ------------------------------------------------------------------
    ! Objective function
    ! ------------------------------------------------------------------
    real(dp) function obj_takeit(dx,PARAM_PP) result(f)
        real(dp),intent(in)            :: dx
        type(parpass_g), intent(inout) :: PARAM_PP
        integer                        :: r,order
        real(dp)                       :: Bdm,nu,b_u,eta
        real(dp)                       :: mb,ms,mb_end,ms_end
        real(dp)                       :: Wsp_App,Ws0_App,Wbp_App
        real(dp)                       :: qx,ux,muv
        real(dp),allocatable           :: KnotsM(:),Coeffw(:,:)
        ! Unpacking parameters
        r      = PARAM_PP%r_p
        Bdm    = PARAM_PP%Bdm_p
        nu     = PARAM_PP%nu_p
        b_u    = PARAM_PP%b_u_p
        eta    = PARAM_PP%eta_p
        mb     = PARAM_PP%mb_p
        ms     = PARAM_PP%ms_p
        muv    = PARAM_PP%mu_p
        order  = 2
        allocate (KnotsM(r+2),Coeffw(r+1,4))
        KnotsM = PARAM_PP%M_p
        Coeffw = PARAM_PP%Coeffw_p
        ms_end = ms + dx
        mb_end = mb - dx
        ! Continuation values
        call eval_spline(Coeffw,KnotsM,r,muv + ms_end,order,Wsp_App)
        call eval_spline(Coeffw,KnotsM,r,muv + ms,    order,Ws0_App)
        call eval_spline(Coeffw,KnotsM,r,muv + mb_end,order,Wbp_App)
        ! Implied q from the binding PC
        qx = ((1.0_dp/Bdm)*(Wsp_App-Ws0_App))**(1.0_dp/nu)
        ux = Util_DM(qx,b_u,eta)
        ! Objective
        f  = ux + Wbp_App
        f  = -f
    end function obj_takeit


end module bargaining