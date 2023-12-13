module solve_welfare
! This module contains subroutines to compute th welfare cost of inflation
!
! Part of "The Long-Run Redistributive Effects of Monetary Policy"
! Christian Bustamante
! Last modified: 13 May 2020
    
use lib_kind
use lib_basic
use lib_interp
use lib_rwhdf5
use get_params
use par_pass
use lib_solver
use tools_fnutil
use tools_value_fns
use tools_distributions
use tools_extra
use bargaining
use solve_steady
implicit none
   
   
contains


    subroutine find_comp(chi_bench,muv_bench, gcomp)
        real(dp),intent(in)  :: chi_bench,muv_bench
        real(dp),intent(out) :: gcomp
        real(dp)             :: Qsc_Int(Ng,Ng),Dsc_Int(Ng,Ng),Gm_Dist(Ng),Dist_F(Ng),Dist_G(Ng),phim
        real(dp)             :: muv_zero,gcomp0
        real(dp)             :: WTotal_bench,WTotal_0
        real(dp)             :: diff,gcomp_l,gcomp_r,gcomp_m,fobj_l,fobj_m,gcomp_sol,fobj_sol

        ! Solve steady state for mu = 10%
        muv_zero  = 1.00_dp**ipyear-1.0_dp
        gcomp0    = 1.0_dp
        call get_steady(chi_bench,muv_bench, Dsc_Int,Qsc_Int,Gm_Dist,Dist_F,Dist_G,phim)

        ! Get lifetime value for mu = 10%
        call get_totalw_comp(chi_bench,muv_bench, Dsc_Int,Qsc_Int,Gm_Dist,Dist_F,Dist_G,phim, gcomp0, WTotal_bench)

        ! Solve steady state for mu = 0%
        call get_steady(chi_bench,muv_zero, Dsc_Int,Qsc_Int,Gm_Dist,Dist_F,Dist_G,phim)

        ! Get lifetime value for mu = 0%
        call get_totalw_comp(chi_bench,muv_zero, Dsc_Int,Qsc_Int,Gm_Dist,Dist_F,Dist_G,phim, gcomp0, WTotal_0)

        print *, WTotal_bench, WTotal_0

        ! Bisecting the compensation
        diff = 2.0_dp
        gcomp_l = 0.85_dp
        gcomp_r = 1.125_dp
        call get_totalw_comp(chi_bench,muv_zero, Dsc_Int,Qsc_Int,Gm_Dist,Dist_F,Dist_G,phim, gcomp_l, WTotal_0)
        fobj_l  = WTotal_0 - WTotal_bench
        do while(diff > 1.0e-5_dp)
            gcomp_m  = (gcomp_l+gcomp_r)/2.0_dp
            call get_totalw_comp(chi_bench,muv_zero, Dsc_Int,Qsc_Int,Gm_Dist,Dist_F,Dist_G,phim, gcomp_m, WTotal_0)
            fobj_m  = WTotal_0 - WTotal_bench
            if (fobj_l*fobj_m < 0.0_dp) then
                gcomp_r = gcomp_m
            else
                gcomp_l = gcomp_m
                fobj_l  = fobj_m
            end if
            diff = abs(gcomp_l-gcomp_r)
            print *, fobj_m, gcomp_m
        end do
        gcomp_sol = gcomp_m
        fobj_sol  = fobj_m

        ! Solution
        gcomp = gcomp_sol
        print *, gcomp, fobj_sol
    end subroutine find_comp




    subroutine get_totalw_comp(chi,muv, Dsc_Int,Qsc_Int,Gm_Dist,Dist_F,Dist_G,phim, gcomp, WTotal)
        real(dp),intent(in)     :: chi,muv,gcomp
        real(dp),intent(in)     :: Qsc_Int(Ng,Ng),Dsc_Int(Ng,Ng),Gm_Dist(Ng),Dist_F(Ng),Dist_G(Ng),phim
        real(dp),intent(out)    :: WTotal
        integer           :: i,j
        real(dp)          :: hl,hr,tolh,fn_soll,fn_solh,fn_sol,hstari,cstari
        real(dp)          :: VWelf,WWelf
        type(parpass_nl)  :: PARAM_HS

        ! Retrieving the rest of parameters and grids
        call get_grids(KnotsM,Grid_M,Nm,Ng,mlow,mhigh)

        !-----------------------------------------------------------------------------
        ! Computing expected lifetime utility
        !-----------------------------------------------------------------------------
        ! Value of DM
        VWelf = 0.0_dp
        do i = 1,Ng
            do j = 1,Ng
                VWelf = VWelf + alpha*0.5_dp*(Util_DM(Qsc_Int(i,j) * gcomp,b_u,eta) &
                      - Cost_HDM(Qsc_Int(i,j),Bdm,nu)) * Dist_F(i)*Dist_F(j)
            end do
        end do

        ! Value of CM
        PARAM_HS%gamma_p = gamma
        PARAM_HS%kappa_p = kappa
        PARAM_HS%chi_p   = chi
        PARAM_HS%phim_p  = phim
        PARAM_HS%mum_p   = muv
        WWelf = 0.0_dp
        do i = 1,Ng
            ! Optimal labor
            PARAM_HS%m0_p    = Grid_M(i)
            PARAM_HS%mf_p    = Gm_Dist(i)
            hl   = 0.0001_dp
            hr   = 15.0_dp
            tolh = 1e-5_dp
            fn_soll = OBJ_HSTAR(hl,PARAM_HS)
            fn_solh = OBJ_HSTAR(hr,PARAM_HS)
            if (abs(chi) < 1.0e-4_dp) then
                ! Analytical solution
                hstari = (1.0_dp/kappa) - phim*(Grid_M(i)-Gm_Dist(i)*(1.0_dp+muv))
                fn_sol = 0.0_dp
            else
                ! Numerical solution
                call bisect(hl,hr,tolh,OBJ_HSTAR,PARAM_HS, hstari,fn_sol)
                hstari = max(0.0_dp, hstari)
            end if
            ! Optimal consumption
            cstari = hstari + phim*(Grid_M(i)-Gm_Dist(i)*(1.0_dp+muv))
            ! Adding period utilities
            WWelf = WWelf + Util_CM(cstari * gcomp,hstari,gamma,kappa,chi) * Dist_G(i)
        end do

        ! Total lifetime value
        WTotal = (1.0_dp/(1.0_dp-beta)) * (VWelf + WWelf)
    end subroutine get_totalw_comp



    subroutine get_steady(chi,muv, Dsc_Int,Qsc_Int,Gm_Dist,Dist_F,Dist_G,phim)
        real(dp),intent(in)     :: chi,muv
        real(dp),intent(out)    :: Qsc_Int(Ng,Ng),Dsc_Int(Ng,Ng),Gm_Dist(Ng),Dist_F(Ng),Dist_G(Ng),phim
        real(dp)          :: phim1,error_phim,adj_phim
        integer           :: i,j,index,fneval
        integer           :: iter_w,iter_mu,iter_ph,iter_nb
        real(dp)          :: diff_w,diff_mu,diff_ph,diff_nb
        real(dp)          :: V(Nm),TV(Nm),W(Nm),TW(Nm),Gm(Nm),Hstar(Nm),Cstar(Nm)
        real(dp)          :: Coeffw(r+1,4),Coeffv(r+1,4)
        real(dp)          :: wint,adj
        real(dp)          :: MBar
        real(dp)          :: Dist_Fp(Ng)
        real(dp)          :: M_Mean_F,M_Mean_G,M_Var_F,M_Var_G
        real(dp)          :: m_lb,m_ub,m_star,w_star
        integer           :: howard,n_how
        type(parpass_g)   :: PARAM_PP
        type(parpass_nl)  :: PARAM_PHIM
        real(dp)          :: begin_nb
        integer           :: hours,mins,secs
        character(60)     :: out_folder,out_fname,out_data_file
        integer(hid_t)    :: fileid
        real(dp)          :: Qsc0(Nm,Nm),Dsc0(Nm,Nm),Objsc0(Nm,Nm)
        real(dp)          :: Qsc1(Nm,Nm),Dsc1(Nm,Nm),Objsc1(Nm,Nm)
        real(dp)          :: Coeffq2(r+1,4,4,r+1),Coeffd2(r+1,4,4,r+1)
        real(dp)          :: Qsc_BFix(Nm,Ng),Dsc_BFix(Nm,Ng)
        real(dp)          :: Qsc_SFix(Ng,Nm),Dsc_SFix(Ng,Nm)
        real(dp)          :: phim_l,phim_h,phim_sol,f_sol

        ! Retrieving the rest of parameters and grids
        call get_grids(KnotsM,Grid_M,Nm,Ng,mlow,mhigh)

        ! Initial guess for distribution F
        MBar   = 1.0_dp  ! Guess for aggregate (fixed) supply of money
        phim   = 1.5_dp  ! Guess for the price of money
        Dist_F = 0.0_dp
        index  = gridlookup(Grid_M,MBar)
        wint   = (MBar-Grid_M(index+1))/(Grid_M(index)-Grid_M(index+1))
        Dist_F(index)   = wint
        Dist_F(index+1) = 1.0_dp-wint

        ! Format labels
700     format(1x,a,i10,a,f12.8,a)      ! Iter + Diff
705     format(1x,a)                    ! String
710     format(1x,a,f10.06,a,f10.6,a)   ! Moments

        ! Initial guess for the terms of trade
        ! Guess for W
        W = KnotsM**(8.0_dp/10.0_dp)
        ! Fitting a spline to W
        call spline_nak(KnotsM,W,r,Coeffw(:,:))


        !-----------------------------------------------------------------------------
        ! Solving the DM problem
        ! Rows: buyer, Cols: seller
        !-----------------------------------------------------------------------------
        ! Take it or leave it
        PARAM_PP%r_p   = r
        PARAM_PP%Bdm_p = Bdm
        PARAM_PP%nu_p  = nu
        PARAM_PP%b_u_p = b_u
        PARAM_PP%eta_p = eta
        PARAM_PP%M_p   = KnotsM
        PARAM_PP%Coeffw_p = Coeffw
        PARAM_PP%mu_p     = muv
        call solve_takeit(KnotsM,Coeffw,tol_pp,PARAM_PP,Qsc0,Dsc0,Objsc0)


        !-----------------------------------------------------------------------------
        ! Solving the model
        ! Loop terms of trade
        !   Loop money market claring (price)
        !     Loop decision rules / value function
        !     Loop distribution
        !----------------------------------------------------------------------------
        diff_nb = 2.0_dp*tol_nb
        iter_nb = 0
        call tick(begin_nb)
        do while (diff_nb > tol_nb)
            iter_nb = iter_nb + 1

            !-----------------------------------------------------------------------------
            ! Bivariate splines for terms of trade
            !-----------------------------------------------------------------------------
            call spline2d_nak(KnotsM,KnotsM,Qsc0,r,r,Coeffq2)
            call spline2d_nak(KnotsM,KnotsM,Dsc0,r,r,Coeffd2)
            ! Evaluating on the finner grid for Dist (for both buyer and seller)
            do i = 1,Ng
                do j = 1,Ng
                    call eval_spline2d(Coeffd2,KnotsM,KnotsM,r,r,Grid_M(i),Grid_M(j), Dsc_Int(i,j))
                    call eval_spline2d(Coeffq2,KnotsM,KnotsM,r,r,Grid_M(i),Grid_M(j), Qsc_Int(i,j))
                end do
            end do
            ! Fixing buyer's or seller's money (useful when computing V, knots for V vs grid for Dist)
            do i = 1,Nm
                do j = 1,Ng
                    call eval_spline2d(Coeffd2,KnotsM,KnotsM,r,r,KnotsM(i),Grid_M(j), Dsc_BFix(i,j))
                    call eval_spline2d(Coeffq2,KnotsM,KnotsM,r,r,KnotsM(i),Grid_M(j), Qsc_BFix(i,j))
                    call eval_spline2d(Coeffd2,KnotsM,KnotsM,r,r,Grid_M(j),KnotsM(i), Dsc_SFix(j,i))
                    call eval_spline2d(Coeffq2,KnotsM,KnotsM,r,r,Grid_M(j),KnotsM(i), Qsc_SFix(j,i))
                end do
            end do


            !-----------------------------------------------------------------------------
            ! Solving for the equilibrium phim
            !-----------------------------------------------------------------------------
            ! Passing
            PARAM_PHIM%b_chi      = chi
            PARAM_PHIM%b_mu       = muv
            PARAM_PHIM%b_W        = W
            PARAM_PHIM%b_Qsc_BFix = Qsc_BFix
            PARAM_PHIM%b_Qsc_SFix = Qsc_SFix
            PARAM_PHIM%b_Dsc_BFix = Dsc_BFix
            PARAM_PHIM%b_Dsc_SFix = Dsc_SFix
            PARAM_PHIM%b_Dsc_Int  = Dsc_Int
            PARAM_PHIM%b_Dist_F   = Dist_F

            ! Bisectinh
            phim_l = 0.5_dp
            phim_h = 2.5_dp
            call bisect(phim_l,phim_h,tol_ph,obj_phim,PARAM_PHIM, phim_sol,f_sol)
            phim = phim_sol

            ! Extracting solution
            W  = PARAM_PHIM%b_W
            V  = PARAM_PHIM%b_V
            Gm = PARAM_PHIM%b_Gm
            Hstar   = PARAM_PHIM%b_Hstar
            Cstar   = PARAM_PHIM%b_Cstar
            Dist_F  = PARAM_PHIM%b_Dist_F
            Dist_G  = PARAM_PHIM%b_Dist_G
            Gm_Dist = PARAM_PHIM%b_Gm_Dist

            ! More
            M_Mean_F = PARAM_PHIM%b_M_Mean_F
            M_Mean_G = PARAM_PHIM%b_M_Mean_G
            M_Var_F  = PARAM_PHIM%b_M_Var_F
            M_Var_G  = PARAM_PHIM%b_M_Var_G


            !-----------------------------------------------------------------------------
            ! 4. Solving the Nash bargaining (Rows: buyer, Cols: seller)
            !-----------------------------------------------------------------------------
            call spline_nak(KnotsM,W,r,Coeffw(:,:))
            PARAM_PP%Coeffw_p = Coeffw
            call solve_takeit(KnotsM,Coeffw,tol_pp,PARAM_PP,Qsc1,Dsc1,Objsc1)

            ! Checking convergence
            diff_nb = maxval(abs(Dsc1-Dsc0))*1
            adj     = 1.00_dp/2.0_dp
            Dsc0    = adj*Dsc1 + (1.0_dp-adj)*Dsc0
            Qsc0    = adj*Qsc1 + (1.0_dp-adj)*Qsc0
            call print_status(iter_w,iter_mu,iter_ph,iter_nb,diff_w,diff_mu,diff_ph/1000,diff_nb,M_Mean_F,M_Mean_G,phim,chi,muv)

        end do
        ! Report
        call print_status(iter_w,iter_mu,iter_ph,iter_nb,diff_w,diff_mu,diff_ph/1000,diff_nb,M_Mean_F,M_Mean_G,phim,chi,muv)

        ! Clock
        call tock_hms(begin_nb,hours,mins,secs)
        print "(1x,a,i2.2,a,i2.2,a,i2.2)", "Time elapsed: ",hours,":",mins,":",secs
        write(*, 705), ' '

    end subroutine get_steady



end module solve_welfare