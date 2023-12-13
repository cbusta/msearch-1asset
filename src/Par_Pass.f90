module par_pass
! This module contains the derived data type variables used to pass arguments into 
! the Lib_Solver subroutines
!
! Part of "The Long-Run Redistributive Effects of Monetary Policy"
! Christian Bustamante
! Last modified: 14 May 2020
    
    use lib_kind
    use get_params
    implicit none

    type parpass_g      ! Parameter passing structure (Golden)
        ! Common to all
        integer  :: r_p
        real(dp) :: M_p(r+2)
        ! When solving the CM
        real(dp) :: gamma_p,kappa_p,chi_p,phim_p,m0_p,beta_p
        real(dp) :: Coeffv_p(r+1,4)
        real(dp) :: hstar_p
        ! When solving take it or leave it
        real(dp) :: Bdm_p,nu_p,b_u_p,eta_p,mb_p,ms_p
        real(dp) :: Coeffw_p(r+1,4)
        real(dp) :: mu_p,gcomp_p
    end type parpass_g

    type parpass_nl     ! Parameter passing structure (Bisect)
        ! Price of money in CM
        integer  :: iter_w_p,iter_mu_p,iter_nb_p,fneval_p
        real(dp) :: diff_w_p,diff_mu_p,diff_nb_p
        real(dp) :: Qsc_Int_p(Ng,Ng),Dsc_Int_p(Ng,Ng)
        real(dp) :: Qsc_BFix_p(Nm,Ng),Dsc_BFix_p(Nm,Ng)
        real(dp) :: Qsc_SFix_p(Ng,Nm),Dsc_SFix_p(Ng,Nm)
        real(dp) :: W_p(Nm),V_p(Nm),Gm_p(Nm),Cstar_p(Nm),Hstar_p(Nm)
        real(dp) :: Gm_Dist_p(Ng),Dist_F_p(Ng),Dist_G_p(Ng)
        real(dp) :: M_Mean_F_p,M_Mean_G_p

        ! Maximizing utility in CM
        real(dp) :: Wsp_p,Ws0_p
        real(dp) :: Bdm_p,nu_p,eta_p,b_u_p

        ! Choosing optimal labor in CM
        real(dp) :: gamma_p,kappa_p,chi_p,phim_p,mum_p
        real(dp) :: m0_p,mf_p
        integer  :: flag

        ! When bisecting phim
        real(dp) :: b_chi,b_mu
        real(dp) :: b_W(Nm),b_V(Nm),b_Gm(Nm),b_Cstar(Nm),b_Hstar(Nm)
        real(dp) :: b_Dist_F(Ng),b_Dist_G(Ng),b_Gm_Dist(Ng)
        real(dp) :: b_Qsc_BFix(Nm,Ng),b_Dsc_BFix(Nm,Ng)
        real(dp) :: b_Qsc_SFix(Ng,Nm),b_Dsc_SFix(Ng,Nm),b_Dsc_Int(Ng,Ng)
        real(dp) :: b_M_Mean_F,b_M_Mean_G,b_M_Var_F,b_M_Var_G

    end type parpass_nl

    type parpass_a      ! Parameter passing structure (Amoeba)
        ! ---
    end type parpass_a

end module par_pass