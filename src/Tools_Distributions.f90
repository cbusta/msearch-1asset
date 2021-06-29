module tools_distributions

! This module contains subroutines to find the stationary distributions of money and bonds on
! both the centralized and decentralized markets. It also computes the distribution of prices.
!
! More money for some
! Christian Bustamante
! Last modified: 13 May 2020

use lib_kind
use lib_basic
use lib_interp
implicit none


contains


    subroutine compute_distg(Grid_M,Dist_F,KnotsM,Dsc_Int,alpha,sigma,muv, Dist_G)
        ! --------------------------------------------------------------------------
        ! Computing Dist_G(m)
        ! Using linear interpolation over F
        ! --------------------------------------------------------------------------
        real(dp), intent(in) :: Grid_M(:),Dist_F(:),KnotsM(:),Dsc_Int(:,:)
        real(dp), intent(in) :: alpha,sigma,muv
        real(dp), intent(out):: Dist_G(size(Grid_M))
        real(dp)             :: m_end,wint,excess
        integer              :: i,j,Ng,index

        Ng = size(Grid_M)
        Dist_G = 0.0_dp
        excess = 0.0_dp
        do i = 1,Ng
            if (Dist_F(i) > 0.0_dp) then

                ! Case 1: Single-coincidence as buyer
                do j = 1,Ng
                    m_end = muv + Grid_M(i) - Dsc_Int(i,j)
                    if (m_end < Grid_M(1))  m_end = Grid_M(1)  ! Forcing dist. to remain inside grid
                    if (m_end > Grid_M(Ng)) m_end = Grid_M(Ng) ! Forcing dist. to remain inside grid
                    index = gridlookup(Grid_M,m_end)
                    if (index < Ng) then
                        wint  = (Grid_M(index+1)-m_end)/(Grid_M(index+1)-Grid_M(index))
                    else
                        wint  = 0.0_dp
                        index = index-1
                    end if
                    Dist_G(index)   = Dist_G(index)   + alpha*sigma*Dist_F(i)*Dist_F(j)*wint
                    Dist_G(index+1) = Dist_G(index+1) + alpha*sigma*Dist_F(i)*Dist_F(j)*(1.0_dp-wint)
                end do

                ! Case 2: Single-coincidence as seller
                do j = 1,Ng
                    m_end = muv + Grid_M(i) + Dsc_Int(j,i)
                    if (m_end < Grid_M(1))  m_end = Grid_M(1)  ! Forcing dist. to remain inside grid
                    if (m_end > Grid_M(Ng)) m_end = Grid_M(Ng) ! Forcing dist. to remain inside grid
                    index = gridlookup(Grid_M,m_end)
                    if (index < Ng) then
                        wint = (Grid_M(index+1)-m_end)/(Grid_M(index+1)-Grid_M(index))
                    else
                        wint  = 0.0_dp
                        index = index-1
                    end if
                    Dist_G(index)   = Dist_G(index)   + alpha*sigma*Dist_F(i)*Dist_F(j)*wint
                    Dist_G(index+1) = Dist_G(index+1) + alpha*sigma*Dist_F(i)*Dist_F(j)*(1.0_dp-wint)
                end do

                ! Case 3: Unamtched or no coincidence
                m_end = muv + Grid_M(i)
                index = gridlookup(Grid_M,m_end)
                if (index < Ng) then
                    wint  = (Grid_M(index+1)-m_end)/(Grid_M(index+1)-Grid_M(index))
                else
                    wint  = 0.0_dp
                    index = index-1
                end if
                Dist_G(index)   = Dist_G(index)   + (1.0_dp-2.0_dp*alpha*sigma)*Dist_F(i)*wint
                Dist_G(index+1) = Dist_G(index+1) + (1.0_dp-2.0_dp*alpha*sigma)*Dist_F(i)*(1.0_dp-wint)

            end if
        end do
        !Dist_G = Dist_G/sum(Dist_G)
    end subroutine compute_distg



    subroutine compute_distf(Grid_M,Dist_G,Gm_Dist,Dist_Fp)
        ! --------------------------------------------------------------------------
        ! Computing  / Updating Dist_F
        ! --------------------------------------------------------------------------
        real(dp), intent(in) :: Grid_M(:),Dist_G(:),Gm_Dist(:)
        real(dp), intent(out):: Dist_Fp(size(Grid_M))
        real(dp)             :: gmpol,m_end,wint
        integer              :: Ng,i,index

        Ng = size(Grid_M)
        Dist_Fp = 0.0_dp
        do i = 1,Ng
            if (Dist_G(i) > 0.0_dp) then
                gmpol = Gm_Dist(i)
                m_end = gmpol
                m_end = min(Grid_M(Ng),max(Grid_M(1),m_end))
                index = gridlookup(Grid_M,m_end)
                if (index < Ng) then
                    wint  = (Grid_M(index+1)-m_end)/(Grid_M(index+1)-Grid_M(index))
                    Dist_Fp(index)   = Dist_Fp(index)   + Dist_G(i)*wint
                    Dist_Fp(index+1) = Dist_Fp(index+1) + Dist_G(i)*(1.0_dp-wint)
                else
                    Dist_Fp(index)   = Dist_Fp(index)   + Dist_G(i)
                end if
            end if
        end do
        Dist_Fp = Dist_Fp/sum(Dist_Fp)
    end subroutine compute_distf



    subroutine compute_distp(Grid_M,KnotsM,Dist_F,Psc,alpha,sigma,Ng,Nm, Grid_P,Dist_P)
        real(dp),intent(in) :: Grid_M(:),KnotsM(:),Dist_F(:),Psc(:,:)
        real(dp),intent(in) :: alpha,sigma
        integer,intent(in)  :: Ng,Nm
        real(dp),intent(out):: Grid_P(Ng),Dist_P(Ng)
        integer             :: i,j,index
        real(dp)            :: plow,phigh,price,wint
        real(dp)            :: Coeffp2(Nm-1,4,4,Nm-1)
        real(dp)            :: Psc_Int(Ng,Ng)
        ! --------------------------------------------------------------------------
        ! Computing Grid_P
        ! --------------------------------------------------------------------------
        Dist_P  = 0.0
        plow  = minval(Psc)
        phigh = maxval(Psc)
        call linspace(plow,phigh,Ng,Grid_P)

        call spline2d_nak(KnotsM,KnotsM,Psc,Nm-2,Nm-2,Coeffp2)
        Psc_Int = 0.0
        do i = 1,Ng
            do j = 1,Ng
                call eval_spline2d(Coeffp2,KnotsM,KnotsM,Nm-2,Nm-2,Grid_M(i),Grid_M(j), Psc_Int(i,j))
            end do
        end do

        do i = 1,Ng
            if (Dist_F(i)  >  0.0_dp) then
                ! Case 1: Unmatched
                ! Does no applies to compute prices

                ! Case 2: Single-coincidence as buyer
                do j = 1,Ng
                    price = Psc_Int(i,j)
                    index = gridlookup(Grid_P,price)
                    if (index  <  Ng) then
                        wint  = (Grid_P(index+1)-price)/(Grid_P(index+1)-Grid_P(index))
                    else
                        wint = 0.0_dp
                        index = index-1
                    end if
                    Dist_P(index)   = Dist_P(index)   + wint*alpha*sigma*Dist_F(i)*Dist_F(j)
                    Dist_P(index+1) = Dist_P(index+1) + (1.0_dp-wint)*alpha*sigma*Dist_F(i)*Dist_F(j)
                end do

                ! Case 3: Single-coincidence as seller
                do j = 1,Ng
                    price = Psc_Int(j,i)
                    index = gridlookup(Grid_P,price)
                    if (index  <  Ng) then
                        wint  = (Grid_P(index+1)-price)/(Grid_P(index+1)-Grid_P(index))
                    else
                        wint = 0.0_dp
                        index = index-1
                    end if
                    Dist_P(index)   = Dist_P(index)   + wint*alpha*sigma*Dist_F(i)*Dist_F(j)
                    Dist_P(index+1) = Dist_P(index+1) + (1.0_dp-wint)*alpha*sigma*Dist_F(i)*Dist_F(j)
                end do

                ! Case 4: No coincidence
                ! Does no applies to compute prices
            end if
        end do
    end subroutine compute_distp



end module tools_distributions