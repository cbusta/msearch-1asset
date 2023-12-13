module solve_steady

! This module containt the subroutine that solves the stationary equilibrium for given
! values of (chi,mu), i.e. the constant Frisch elasticity of labor supply and the money
! growth rate.
! Modified Lagos-Wright model that allows for non-quasilinear preferences.
!
! Part of "The Long-Run Redistributive Effects of Monetary Policy"
! Christian Bustamante
! Last modified: 14 May 2020

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
implicit none

contains


    subroutine solve_4chimu(chi,chi_str,muv,mu_str)
        ! Solves the stationary equilibrium
        real(dp),intent(in)     :: chi,muv
        character(3),intent(in) :: chi_str,mu_str
        real(dp)          :: phim
        integer           :: i,j,index
        integer           :: iter_w,iter_mu,iter_ph,iter_nb
        real(dp)          :: diff_w,diff_mu,diff_ph,diff_nb
        real(dp)          :: V(Nm),W(Nm),Gm(Nm),Hstar(Nm),Cstar(Nm)
        real(dp)          :: Coeffw(r+1,4)
        real(dp)          :: wint,adj
        real(dp)          :: MBar
        real(dp)          :: Dist_F(Ng),Dist_G(Ng)
        real(dp)          :: M_Mean_F,M_Mean_G,M_Var_F,M_Var_G
        real(dp)          :: Gm_Dist(Ng)
        type(parpass_g)   :: PARAM_PP
        type(parpass_nl)  :: PARAM_PHIM
        integer           :: begin_nb,hours,mins,secs,timec(3)
        integer           :: vdate(8)
        character(60)     :: out_folder,out_fname,out_data_file
        integer(hid_t)    :: fileid
        integer(4)        :: hostnm, status
        character(20)     :: machname
        character(8)      :: dstr
        character(6)      :: tstr
        real(dp)          :: Qsc0(Nm,Nm),Dsc0(Nm,Nm),Objsc0(Nm,Nm)
        real(dp)          :: Qsc1(Nm,Nm),Dsc1(Nm,Nm),Objsc1(Nm,Nm)
        real(dp)          :: Coeffq2(r+1,4,4,r+1),Coeffd2(r+1,4,4,r+1)
        real(dp)          :: Qsc_Int(Ng,Ng),Dsc_Int(Ng,Ng)
        real(dp)          :: Qsc_BFix(Nm,Ng),Dsc_BFix(Nm,Ng)
        real(dp)          :: Qsc_SFix(Ng,Nm),Dsc_SFix(Ng,Nm)
        real(dp)          :: Psc(Nm,Nm),Grid_P(Ng),Dist_P(Ng)
        real(dp)          :: Mkup(Nm,Nm),Grid_Mkup(Ng),Dist_Mkup(Ng)
        real(dp)          :: dx,qx,Dv,mtild,htild,Dwm,Mc
        real(dp)          :: CAgg_DM_Nom,CAgg_DM_Real,CAgg_CM_Nom,CAgg_CM_Real,HAgg_CM
        real(dp)          :: YAgg_Nom,YAgg_Real,Velo,inom
        real(dp)          :: Mean_P,Var_P,Mean_Mkup,Var_Mkup
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
        call tick_wc(begin_nb)
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
                    call eval_spline2d(Coeffd2,KnotsM,KnotsM,r,r,KnotsM(i),Grid_M(j),Dsc_BFix(i,j))
                    call eval_spline2d(Coeffq2,KnotsM,KnotsM,r,r,KnotsM(i),Grid_M(j),Qsc_BFix(i,j))
                    call eval_spline2d(Coeffd2,KnotsM,KnotsM,r,r,Grid_M(j),KnotsM(i),Dsc_SFix(j,i))
                    call eval_spline2d(Coeffq2,KnotsM,KnotsM,r,r,Grid_M(j),KnotsM(i),Qsc_SFix(j,i))
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

            ! Bisecting
            phim_l = 0.5
            phim_h = 2.5
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
            diff_nb = maxval(abs(Dsc1-Dsc0))
            adj     = 1.00_dp/2_dp
            Dsc0    = adj*Dsc1 + (1.0_dp-adj)*Dsc0
            Qsc0    = adj*Qsc1 + (1.0_dp-adj)*Qsc0
            call print_status(iter_w,iter_mu,iter_ph,iter_nb,diff_w,diff_mu,diff_ph/1000,diff_nb,M_Mean_F,M_Mean_G,phim,chi,muv)

        end do
        ! Report
        call print_status(iter_w,iter_mu,iter_ph,iter_nb,diff_w,diff_mu,diff_ph/1000,diff_nb,M_Mean_F,M_Mean_G,phim,chi,muv)


        !-----------------------------------------------------------------------------
        ! Computing some equilibrium results
        !-----------------------------------------------------------------------------

        ! Computing prices per unit in the DM
        Psc = Dsc0/Qsc0
        call compute_distp(Grid_M,KnotsM,Dist_F,Psc,alpha,sigma,Ng,Nm, Grid_P,Dist_P)

        ! Markup (price over marginal cost)
        ! Note that price is defined as d, nor as price per unit since mc are over q, not 1
        do i = 1,Nm
            do j = 1,Nm
                dx    = Dsc0(i,j)
                qx    = Qsc0(i,j)
                Dv    = MDisu_HDM(qx,Bdm,nu)
                mtild = KnotsM(j) + muv + dx
                htild = linterpx(Hstar,KnotsM,mtild)
                Dwm   = kappa*(htild**chi)*phim
                Mc    = Dv*(1.0_dp/Dwm)
                Mkup(i,j) = dx/Mc
            end do
        end do
        ! Distribution of Mkups (can use same subroutine as for P)
        call compute_distp(Grid_M,KnotsM,Dist_F,Mkup,alpha,sigma,Ng,Nm, Grid_Mkup,Dist_Mkup)


        ! ------------------------------------------------
        ! Computing some macro aggregates
        ! ------------------------------------------------
        ! Steady State
        ! Consumption
        call agg_dm(Dsc_Int,Dist_F,alpha,sigma, CAgg_DM_Nom)
        call agg_cm(Gm_Dist,Dist_G,Grid_M,kappa,chi,phim,muv, CAgg_CM_Real,HAgg_CM)
        ! Total output (in terms of numerarie)
        CAgg_CM_Nom  = (1.0_dp/phim)*CAgg_CM_Real
        CAgg_DM_Real = phim*CAgg_DM_Nom
        YAgg_Nom  = CAgg_DM_Nom + CAgg_CM_Nom
        YAgg_Real = Yagg_Nom*phim
        ! Interest rate
        inom = (1.0_dp/beta)*(1.0_dp+muv) - 1.0_dp
        ! Velocity
        Velo = YAgg_Nom/MBar

        ! Moments of distributions
        ! Prices
        call momdist1v(Dist_P,Grid_P, Mean_P,Var_P)
        ! Markups
        call momdist1v(Dist_Mkup,Grid_Mkup, Mean_Mkup,Var_Mkup)


        ! Clock
        call tock_hms_wc(begin_nb,hours,mins,secs)
        print "(1x,a,i2.2,a,i2.2,a,i2.2)", "Time elapsed: ",hours,":",mins,":",secs
        write(*, 705), ' '


        ! ------------------------------------------------
        ! Exporting arrays
        ! ------------------------------------------------
        out_folder    = 'out'
        out_data_file = trim(out_folder)//'/OutData_Steady_Chi'//chi_str//'_Mu'//mu_str//'.h5'
        call system( 'mkdir '// trim(out_folder))

        ! Create a new file using the default properties
        call hdf5_openf(out_data_file,fileid)

        ! Writing datasets to h5 file
        call hdf5_write(KnotsM,      fileid,'KnotsM')
        call hdf5_write(Grid_M,      fileid,'Grid_M')
        call hdf5_write(W,           fileid,'W')
        call hdf5_write(V,           fileid,'V')
        call hdf5_write(Gm,          fileid,'Gm')
        call hdf5_write(Hstar,       fileid,'Hstar')
        call hdf5_write(Cstar,       fileid,'Cstar')
        call hdf5_write(Dist_F,      fileid,'Dist_F')
        call hdf5_write(Dist_G,      fileid,'Dist_G')
        call hdf5_write(Dsc0,        fileid,'Dsc0')
        call hdf5_write(Qsc0,        fileid,'Qsc0')
        call hdf5_write(Psc,         fileid,'Psc')
        call hdf5_write(Mkup,        fileid,'Mkup')
        call hdf5_write(Grid_P,      fileid,'Grid_P')
        call hdf5_write(Grid_Mkup,   fileid,'Grid_Mkup')
        call hdf5_write(Dist_P,      fileid,'Dist_P')
        call hdf5_write(Dist_Mkup,   fileid,'Dist_Mkup')
        call hdf5_write(muv,         fileid,'Inf_Rate')
        call hdf5_write(CAgg_DM_Real,fileid,'CAgg_DM_Real')
        call hdf5_write(CAgg_CM_Real,fileid,'CAgg_CM_Real')
        call hdf5_write(HAgg_CM,     fileid,'HAgg_CM')
        call hdf5_write(YAgg_Real,   fileid,'YAgg_Real')
        call hdf5_write(MBar,        fileid,'MBar')
        call hdf5_write(Velo,        fileid,'Velo')
        call hdf5_write(phim,        fileid,'phim')
        call hdf5_write(beta,        fileid,'Par_Beta')
        call hdf5_write(rreal,       fileid,'Par_rReal')
        call hdf5_write(ipyear,      fileid,'Par_pYear')

        ! Close file
        call hdf5_closef(fileid)



        ! ------------------------------------------------
        ! Reporting file
        ! ------------------------------------------------
91      format(1x,a,f12.8)      ! Text + Real
93      format(1x,a,f8.4)       ! Text + Parameter value
95      format(1x,a)            ! String
97      format(1x,a,a)          ! String (x2)
        call itime(timec)
        call date_and_time(values=vdate)

        status = hostnm(machname)
        call date_str(dstr)
        call time_str(tstr)
        out_fname  = trim(trim(out_folder)//'/_Results_Steady_'//chi_str//'_'//mu_str//'_'//dstr//'.log')
        open(7,file= trim(out_fname))
        write(7,95) , '-------------------------------------------------------'
        write(7,95) , 'Search Theoretic Model of Money'
        write(7,95) , 'Benchmark Modification of Lagos & Wright (2005)'
        write(7,93) , 'Growth rate of money (%) = ', 100*muv
        write(7,"(1x,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i4)"), &
                      'Results compiled at ',timec(1),':',timec(2),':',timec(3), &
                      ' on ',vdate(3),'/',vdate(2),'/',vdate(1)
        write(7,97) , 'Machine name: ', machname
        write(7,95) , '-------------------------------------------------------'
        write(7,95) , ''
        write(7,95) , ''
        write(7,95) , '-------------------------------------------------------'
        write(7,95) , 'Steady state values'
        write(7,95) , '-------------------------------------------------------'
        write(7,91) , 'Total output             = ', YAgg_Real
        write(7,91) , 'Output in DM             = ', CAgg_DM_Real
        write(7,91) , 'Output in CM             = ', HAgg_CM
        write(7,91) , 'Consumption in CM        = ', CAgg_CM_Real
        write(7,91) , 'Price of money (phiM)    = ', phim
        write(7,91) , 'Inflation rate           = ', muv
        write(7,91) , 'Velocity of money        = ', Velo
        write(7,91) , 'Real interest rate       = ', rreal
        write(7,95) , ''
        write(7,95) , ''
        write(7,95) , '-------------------------------------------------------'
        write(7,95) , 'Moments of the distribution of money'
        write(7,95) , '-------------------------------------------------------'
        write(7,91) , 'Mean of money at DM      = ', M_Mean_F
        write(7,91) , 'Mean of money at CM      = ', M_Mean_G
        write(7,91) , 'Std. dev. of money at DM = ', sqrt(M_Var_F)
        write(7,91) , 'Std. dev. of money at CM = ', sqrt(M_Var_G)
        write(7,95) , ''
        write(7,95) , ''
        write(7,95) , '-------------------------------------------------------'
        write(7,95) , 'Moments of the distribution of prices and markups'
        write(7,95) , '-------------------------------------------------------'
        write(7,91) , 'Average price in the DM   = ', Mean_P
        write(7,91) , 'Std. dev. of prices in DM = ', sqrt(Var_P)
        write(7,91) , 'Average markup in the DM  = ', Mean_Mkup
        write(7,91) , 'Std. dev. of mkup in DM   = ', sqrt(Var_Mkup)
        write(7,95) , ''
        write(7,95) , ''
        write(7,95) , '-------------------------------------------------------'
        write(7,95) , 'Parameters'
        write(7,95) , '-------------------------------------------------------'
        write(7,93) , 'Prob. of being matched  (alpha)       = ', alpha
        write(7,93) , 'Prob. of single coincidence (sigma)   = ', sigma
        write(7,93) , 'Model periods per year                = ', pyear
        write(7,95) , '-------------------------------------------------------'
        write(7,93) , 'Risk aversion in DM (eta)             = ', eta
        write(7,93) , 'Scale disutility labor in DM (B)      = ', Bdm
        write(7,93) , 'Curvature disutility of labor (nu)    = ', nu
        write(7,93) , 'Risk aversion in CM (gamma)           = ', gamma
        write(7,93) , 'Scale disutility labor in CM (kappa)  = ', kappa
        write(7,93) , 'Inv. Frisch elasticity in CM (chi)    = ', chi
        write(7,95) , '-------------------------------------------------------'
        write(7,93) , 'Number of grid points for money       = ', DBLE(Nm)
        write(7,93) , 'Number of grid points for dist money  = ', DBLE(Ng)
        write(7,95) , ''
        write(7,95) , ''
        write(7,"(a,i2.2,a,i2.2,a,i2.2)"), 'Time elapsed: ',hours,':',mins,':',secs
        close(7)

    end subroutine solve_4chimu



    real(dp) function obj_phim(phim,PARAM_PHIM) result(fobj)
        ! Solves the equilibrium in the money market, *given* some terms of trade
        real(dp),intent(in)            :: phim
        type(parpass_nl),intent(inout) :: PARAM_PHIM
        integer           :: i
        integer           :: iter_w,iter_mu
        real(dp)          :: diff_w,diff_mu
        real(dp)          :: V(Nm),W(Nm),TW(Nm),Gm(Nm),Hstar(Nm),Cstar(Nm)
        real(dp)          :: Coeffw(r+1,4),Coeffv(r+1,4)
        real(dp)          :: adj,MBar
        real(dp)          :: Dist_F(Ng),Dist_G(Ng),Dist_Fp(Ng)
        real(dp)          :: M_Mean_F,M_Mean_G,M_Var_F,M_Var_G
        real(dp)          :: m_lb,m_ub,m_star,w_star,Gm_Dist(Ng)
        integer           :: howard,n_how
        type(parpass_g)   :: PARAM_CM
        real(dp)          :: Dsc_Int(Ng,Ng)
        real(dp)          :: Qsc_BFix(Nm,Ng),Dsc_BFix(Nm,Ng)
        real(dp)          :: Qsc_SFix(Ng,Nm),Dsc_SFix(Ng,Nm)
        real(dp)          :: chi,muv

        ! Format labels
700     format(1x,a,i10,a,f12.8,a)      ! Iter + Diff

        ! Extracting
        chi = PARAM_PHIM%b_chi
        muv = PARAM_PHIM%b_mu
        W   = PARAM_PHIM%b_W
        Qsc_BFix = PARAM_PHIM%b_Qsc_BFix
        Qsc_SFix = PARAM_PHIM%b_Qsc_SFix
        Dsc_BFix = PARAM_PHIM%b_Dsc_BFix
        Dsc_SFix = PARAM_PHIM%b_Dsc_SFix
        Dsc_Int  = PARAM_PHIM%b_Dsc_Int
        Dist_F   = PARAM_PHIM%b_Dist_F

        ! Retrieving the rest of parameters and grids
        call get_grids(KnotsM,Grid_M,Nm,Ng,mlow,mhigh)

        ! Passing parameters (centralized market)
        PARAM_CM%gamma_p = gamma
        PARAM_CM%kappa_p = kappa
        PARAM_CM%chi_p   = chi
        PARAM_CM%beta_p  = beta
        PARAM_CM%r_p     = r
        PARAM_CM%phim_p  = phim
        PARAM_CM%mu_p    = muv
        PARAM_CM%M_p     = KnotsM

        !-----------------------------------------------------------------------------
        ! 1. For a given guess for Dist_F (the intermediate Dist_G isn't important here)
        !    and (Q*,D*) solve V and W
        !-----------------------------------------------------------------------------
        diff_w = 2.0_dp*tol_w
        iter_w = 0
        do while (diff_w > tol_w)
            iter_w = iter_w + 1

            ! Fitting a spline to W
            call spline_nak(KnotsM,W,r,Coeffw(:,:))

            ! Compute expected V(m)
            call compute_ev(KnotsM,Qsc_BFix,Qsc_SFix,Dsc_BFix,Dsc_SFix,Coeffw,Grid_M,Dist_F,alpha,sigma,b_u,eta,Bdm,nu,muv, V)

            ! Solve the CM
            ! Using Howard's improvement algorithm (update decision rules every n iter)
            PARAM_CM%phim_p = phim
            howard = 1
            n_how  = 25
            if (howard == 1) then
                if ((iter_w == 1) .or. (mod(iter_w,n_how) == 0)) then
                    call solve_cm(KnotsM,V,phim,kappa,chi,tol_cm,PARAM_CM,TW,Gm,Hstar,Cstar)
                else
                    call spline_nak(KnotsM,V,r,Coeffv(:,:))
                    PARAM_CM%Coeffv_p = Coeffv
                    do i = 1,Nm
                        PARAM_CM%m0_p = KnotsM(i)
                        TW(i) = -obj_cm(Gm(i),PARAM_CM)
                    end do
                end if
            else
                call solve_cm(KnotsM,V,phim,kappa,chi,tol_cm,PARAM_CM,TW,Gm,Hstar,Cstar)
            end if

            ! Checking convergence
            diff_w = maxval(abs(TW-W))
            adj    = 1.0_dp
            W      = adj*TW + (1_dp-adj)*W
            if (mod(iter_w,10000) == 0) then
                write(*, 700), '[iter_w ,Diff_w ] = [',iter_w,', ',diff_w,']'
            end if
        end do

        !-----------------------------------------------------------------------------
        ! 2. Updating the distributions
        !-----------------------------------------------------------------------------

        ! Getting the policy for finner grid
        call spline_nak(KnotsM,V,r,Coeffv(:,:))
        PARAM_CM%phim_p   = phim
        PARAM_CM%Coeffv_p = Coeffv(:,:)
        do i = 1,Ng
            PARAM_CM%m0_p = Grid_M(i)
            m_lb  = KnotsM(1)
            m_ub  = min(Grid_M(Ng), (1.0_dp/(1.0_dp+muv)) * ((1.0_dp/phim) + Grid_M(i)) - 1.0e-4_dp)
            call golden(m_lb,m_ub,tol_w,obj_cm,PARAM_CM, m_star,w_star)
            Gm_Dist(i)  = m_star
        end do

        ! Finding stationary distribution
        diff_mu = 2.0_dp*tol_mu
        iter_mu = 0
        do while (diff_mu > tol_mu)
            iter_mu = iter_mu + 1

            call compute_distg(Grid_M,Dist_F,KnotsM,Dsc_Int,alpha,sigma,muv, Dist_G)
            call compute_distf(Grid_M,Dist_G,Gm_Dist, Dist_Fp)
            diff_mu = maxval(abs(Dist_Fp-Dist_F))*1
            adj     = 1.0_dp
            Dist_F  = adj*Dist_Fp + (1.0_dp-adj)*Dist_F
        end do

        !-----------------------------------------------------------------------------
        ! 3. Solving for market clearing in money market
        !-----------------------------------------------------------------------------
        ! Mean and variance of money holdings before DM
        M_Mean_F = 0.0_dp
        M_Mean_G = 0.0_dp
        M_Var_F  = 0.0_dp
        M_Var_G  = 0.0_dp
        do i = 1,Ng
            M_Mean_F = M_Mean_F + Dist_F(i)*(Grid_M(i))
            M_Mean_G = M_Mean_G + Dist_G(i)*(Grid_M(i))
            M_Var_F  = M_Var_F  + Dist_F(i)*(Grid_M(i)**2.0_dp)
            M_Var_G  = M_Var_G  + Dist_G(i)*(Grid_M(i)**2.0_dp)
        end do
        M_Var_F = M_Var_F - M_Mean_F**2.0_dp
        M_Var_G = M_Var_G - M_Mean_G**2.0_dp
        MBar    = M_Mean_F

        ! Objective
        fobj = MBar - 1.0_dp

        ! Exporting
        PARAM_PHIM%b_W  = W
        PARAM_PHIM%b_V  = V
        PARAM_PHIM%b_Gm = Gm
        PARAM_PHIM%b_Hstar   = Hstar
        PARAM_PHIM%b_Cstar   = Cstar
        PARAM_PHIM%b_Dist_F  = Dist_F
        PARAM_PHIM%b_Dist_G  = Dist_G
        PARAM_PHIM%b_Gm_Dist = Gm_Dist

        ! More
        PARAM_PHIM%b_M_Mean_F = M_Mean_F
        PARAM_PHIM%b_M_Mean_G = M_Mean_G
        PARAM_PHIM%b_M_Var_F  = M_Var_F
        PARAM_PHIM%b_M_Var_G  = M_Var_G

    end function obj_phim




end module solve_steady
