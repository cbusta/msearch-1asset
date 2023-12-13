module tools_extra
! This module contains additional subroutines to report the current status and perform
! aggregation
!
! Part of "The Long-Run Redistributive Effects of Monetary Policy"
! Christian Bustamante
! Last modified: 15 May 2020

use lib_kind
use lib_basic
use par_pass
use lib_solver
use lib_interp
use tools_fnutil
use tools_value_fns
implicit none


contains


    subroutine print_status(iter_w,iter_mu,iter_ph,iter_nb,diff_w,diff_mu,diff_ph,diff_nb,M_Mean_F,M_Mean_G,phim_now,chi,muval)
        integer, intent(in)  :: iter_w,iter_mu,iter_ph,iter_nb
        real(dp), intent(in) :: diff_w,diff_mu,diff_ph,diff_nb,M_Mean_F,M_Mean_G,phim_now,chi,muval
        ! Format labels
100     format(1x,a,i10,a,f12.8,a)      ! Iter + Diff
105     format(1x,a)                    ! String
110     format(1x,a,f10.6,a,f10.6,a)    ! Moments
115     format(1x,a,f10.6)              ! String + Real
        write(*, 105), ' '
        write(*, 105), '-------------------------------------------------------'
        write(*, 105), 'Current status '
        write(*, 115), 'Price of money = ', phim_now
        write(*, 105), '-------------------------------------------------------'
        write(*, 110), '[mu,     chi    ] = [',muval,   ', ',chi    ,'  ]'
        write(*, 100), '[iter_nb,Diff_nb] = [',iter_nb ,', ',diff_nb,  ']'
        write(*, 110), '[Mean_M, Mean_G ] = [',M_Mean_F,', ',M_Mean_G,'  ]'
        write(*, 105), '-------------------------------------------------------'
        write(*, 105), ' '
    end subroutine print_status


    subroutine agg_dm(Dsc_Int,Dist_F,alpha,sigma, CAgg_DM)
        real(dp),intent(in) :: Dsc_Int(:,:),Dist_F(:),alpha,sigma
        real(dp),intent(out):: CAgg_DM
        integer             :: i,j,Ng
        real(dp)            :: cij_b,cij_s,prob
        CAgg_DM = 0.0_dp
        Ng = size(Dist_F)
        do i = 1,Ng
            if (Dist_F(i)  >  0.0_dp) then
                do j = 1,Ng
                    cij_b = Dsc_Int(i,j)
                    cij_s = Dsc_Int(j,i)
                    prob  = Dist_F(i)*Dist_F(j)/2.0_dp
                    CAgg_DM = CAgg_DM + (alpha*sigma)*prob*cij_b
                    CAgg_DM = CAgg_DM + (alpha*sigma)*prob*cij_s
                end do
            end if
        end do
    end subroutine agg_dm


    subroutine agg_cm(Gm_Dist,Dist_G,Grid_M,kappa,chi,phim,muv, CAgg_CM,HAgg_CM)
        real(dp),intent(in) :: Gm_Dist(:),Dist_G(:),Grid_M(:),kappa,chi,phim,muv
        real(dp),intent(out):: CAgg_CM,HAgg_CM
        real(dp)            :: MAgg_CM,AAgg_CM,Delt_CM
        integer             :: i,Ng
        real(dp)            :: m0,mstar,hstar,cstar,delta
        real(dp)            :: hl,hr,tolh,fn_sol
        type(parpass_nl)    :: PARAM_HS
        Ng      = size(Dist_G)
        CAgg_CM = 0.0_dp
        HAgg_CM = 0.0_dp
        MAgg_CM = 0.0_dp
        AAgg_CM = 0.0_dp
        Delt_CM = 0.0_dp
        do i = 1,Ng
            if (Dist_G(i)  >  0.0_dp) then
                m0    = Grid_M(i)
                mstar = Gm_Dist(i)

                ! Optimal labor
                PARAM_HS%kappa_p = kappa
                PARAM_HS%chi_p   = chi
                PARAM_HS%phim_p  = phim
                PARAM_HS%mum_p   = muv
                PARAM_HS%m0_p    = m0
                PARAM_HS%mf_p    = mstar
                hl   = 0.001_dp
                hr   = 15.0_dp
                tolh = 1.0e-5_dp
                if (abs(chi) < 1e-4) then
                    ! Analytical solution
                    hstar  = (1.0_dp/kappa) - phim*(m0-mstar*(1.0_dp+muv))
                    fn_sol = 0.0_dp
                else
                    ! Numerical solution
                    call bisect(hl,hr,tolh,obj_hstar,PARAM_HS, hstar,fn_sol)
                    hstar = max(0.0_dp, hstar)
                end if
                cstar = hstar + phim*(m0 - mstar*(1.0_dp+muv))
                delta = m0-mstar*(1.0_dp+muv)
                ! Aggregating
                CAgg_CM = CAgg_CM + Dist_G(i)*cstar
                HAgg_CM = HAgg_CM + Dist_G(i)*hstar
                MAgg_CM = MAgg_CM + Dist_G(i)*m0
                AAgg_CM = AAgg_CM + Dist_G(i)*mstar
                Delt_CM = Delt_CM + Dist_G(i)*delta
            end if
        end do
        !print *, CAgg_CM-HAgg_CM,Delt_CM,MAgg_CM,AAgg_CM
    end subroutine agg_cm


    subroutine momdist1v(Mu,Grid, MeanD,VarD)
        ! Compute first two moments of an univariate distribution definined on a grid
        real(dp),intent(in)  :: Mu(:),Grid(:)
        real(dp),intent(out) :: MeanD,VarD
        integer              :: N,i
        N = size(Grid)
        MeanD = 0.0_dp
        VarD  = 0.0_dp
        do i = 1,N
            MeanD = MeanD + Mu(i)*(Grid(i))
            VarD  = VarD  + Mu(i)*(Grid(i)**2.0_dp)
        end do
        VarD = VarD - MeanD**2.0_dp
    end subroutine momdist1v



end module tools_extra