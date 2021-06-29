    module tools_value_fns

    ! This module contains the subroutines to evaluate the value functions in the decentralized and
    ! centralized market. In the centralized market case, we need to solve for the optimal decision rules.
    ! In the decentralized market case, we only evaluate V using the ToT computed in Barganing.f90
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


    subroutine compute_ev(KnotsM,Qsc_BFix,Qsc_SFix,Dsc_BFix,Dsc_SFix,Coeffw,Grid_M,Dist_F,alpha,sigma,b_u,eta,Bdm,nu,muv,V)
        ! --------------------------------------------------------------------------
        ! Compute expected V(m)
        ! This function will be used as continuation value when solving CM problem
        ! --------------------------------------------------------------------------
        real(dp), intent(in) :: KnotsM(:),Coeffw(:,:),Grid_M(:),Dist_F(:)
        real(dp), intent(in) :: Qsc_BFix(:,:),Dsc_BFix(:,:),Qsc_SFix(:,:),Dsc_SFix(:,:)
        real(dp), intent(in) :: alpha,sigma,b_u,eta,Bdm,nu,muv
        real(dp), intent(out):: V(Nm)
        real(dp)             :: m_end,W_App,Aux_Vm
        integer              :: i,j,order

        order = 2
        !$omp parallel do private(i,Aux_Vm,j,m_end,W_App)
        do i = 1,Nm
            Aux_Vm = 0.0
            ! Value of beign a buyer (single coincidence)
            do j = 1,Ng
                m_end = muv + KnotsM(i) - Dsc_BFix(i,j)
                call eval_spline(Coeffw(:,:),KnotsM,r,m_end,order, W_App)
                Aux_Vm = Aux_Vm + alpha*sigma * (Util_DM(Qsc_BFix(i,j),b_u,eta) + W_App) * Dist_F(j)
            end do
            ! Value of beign a seller (single coincidence)
            do j = 1,Ng
                m_end = muv + KnotsM(i) + Dsc_SFix(j,i)
                call eval_spline(Coeffw(:,:),KnotsM,r,m_end,order, W_App)
                Aux_Vm = Aux_Vm + alpha*sigma * (Cost_HDM(Qsc_SFix(j,i),Bdm,nu) + W_App) * Dist_F(j)
            end do
            ! Value of not being matched
            m_end = muv + KnotsM(i)
            call eval_spline(Coeffw(:,:),KnotsM,r,m_end,order, W_App)
            Aux_Vm = Aux_Vm + (1.0_dp-2.0_dp*alpha*sigma)*W_App
            ! Total value
            V(i) = Aux_Vm
        end do
        !$omp end parallel do
    end subroutine compute_ev


    subroutine solve_cm(KnotsM,V,phim,kappa,chi,tolw,PARAM_CM,TW,Gm,Hstar,Cstar)
        ! --------------------------------------------------------------------------
        ! Solving the Centralized Market
        ! --------------------------------------------------------------------------
        real(dp), intent(in) :: KnotsM(:),V(:)
        real(dp), intent(in) :: phim,kappa,chi,tolw
        real(dp), intent(out):: TW(Nm),Gm(Nm),Hstar(Nm),Cstar(Nm)
        real(dp)             :: Coeffv(Nm-1,4)
        real(dp)             :: m_lb,m_ub,m_star,w_star,m0,muv
        integer              :: i
        type(parpass_g), intent(inout) :: PARAM_CM

        ! Fitting a splines to V
        muv = PARAM_CM%mu_p
        call spline_nak(KnotsM,V,r,Coeffv(:,:))
        PARAM_CM%Coeffv_p = Coeffv(:,:)

        ! !$omp parallel do private(i,m0,m_lb,m_ub,m_star,w_star) firstprivate(PARAM_CM)
        do i = 1,Nm
            m0    = KnotsM(i)
            PARAM_CM%m0_p = m0
            m_lb  = KnotsM(1)
            m_ub  = KnotsM(r+2)
            call golden(m_lb,m_ub,tolw,OBJ_CM,PARAM_CM, m_star,w_star)
            TW(i) = -w_star
            Gm(i) =  m_star
            Hstar(i) = PARAM_CM%hstar_p
            Cstar(i) = Hstar(i) + phim*(m0 - m_star*(1.0_dp+muv))
        end do
        ! !$omp end parallel do
    end subroutine solve_cm


    ! --------------------------------------------------------------------------------------
    ! Objective functions
    ! --------------------------------------------------------------------------------------
    real(dp) function obj_cm(mp,PARAM_CM) result(obj)
        real(dp),intent(in)            :: mp
        type(parpass_g), intent(inout) :: PARAM_CM
        integer                        :: r,order
        real(dp)                       :: gamma,kappa,chi,beta,phim,muv,m0,hstar,c,V_App
        real(dp)                       :: hl,hr,tolh,fn_sol
        real(dp)                       :: fn_soll, fn_solh
        real(dp)                       :: m_end
        real(dp),allocatable           :: KnotsM(:),Coeff(:,:)
        type(parpass_nl)               :: PARAM_HS
        gamma = PARAM_CM%gamma_p
        kappa = PARAM_CM%kappa_p
        chi   = PARAM_CM%chi_p
        beta  = PARAM_CM%beta_p
        phim  = PARAM_CM%phim_p
        muv   = PARAM_CM%mu_p
        m0    = PARAM_CM%m0_p
        r     = PARAM_CM%r_p
        order = 2
        allocate (KnotsM(r+2),Coeff(r+1,4))
        KnotsM = PARAM_CM%M_p
        Coeff  = PARAM_CM%Coeffv_p
        m_end  = mp
        call eval_spline(Coeff,KnotsM,r,m_end,order, V_App)
        ! Optimal labor
        PARAM_HS%gamma_p = gamma
        PARAM_HS%kappa_p = kappa
        PARAM_HS%chi_p   = chi
        PARAM_HS%phim_p  = phim
        PARAM_HS%mum_p   = muv
        PARAM_HS%m0_p    = m0
        PARAM_HS%mf_p    = mp
        hl   = 0.0001_dp
        hr   = 15.0_dp
        tolh = 1.0e-5_dp
        fn_soll = obj_hstar(hl,PARAM_HS)
        fn_solh = obj_hstar(hr,PARAM_HS)

        if (abs(chi) < 1.0e-4_dp) then
            ! Analytical solution
            hstar  = (1.0_dp/kappa) - phim*(m0-mp*(1.0_dp+muv))
            fn_sol = 0.0_dp
        else
            ! Numerical solution
            call bisect(hl,hr,tolh,obj_hstar,PARAM_HS, hstar,fn_sol)
            hstar = max(0.0_dp, hstar)
        end if
        c     = hstar + phim*(m0 - mp*(1.0_dp+muv))
        obj   = util_cm(c,hstar,gamma,kappa,chi) + beta*V_App
        obj   =-obj
        PARAM_CM%hstar_p = hstar
    end function obj_cm



    real(dp) function obj_hstar(h,PARAM_HS) result(obj)
        real(dp),intent(in)            :: h
        type(parpass_nl), intent(inout):: PARAM_HS
        real(dp)                       :: gamma,kappa,chi,phim,muv,m0,mf
        gamma = PARAM_HS%gamma_p
        kappa = PARAM_HS%kappa_p
        chi   = PARAM_HS%chi_p
        phim  = PARAM_HS%phim_p
        muv   = PARAM_HS%mum_p
        m0    = PARAM_HS%m0_p
        mf    = PARAM_HS%mf_p
        ! Objective function
        obj   = kappa*(h**(chi)) * ((h + phim*(m0-mf*(1.0_dp+muv)))**gamma) - 1.0_dp
    end function obj_hstar


end module tools_value_fns