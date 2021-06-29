module tools_fnutil
! This module contains the several functions that evaluate instantenous utility and 
! costs functions (including some marginal utilities and costs)
!
! More money for some
! Christian Bustamante
! Last modified: 13 May 2020

use lib_kind
implicit none

contains

    real(dp) function Util_DM(qc,b_u,eta)
        real(dp),intent(in) :: qc,b_u,eta
        Util_DM = (1.0_dp/(1.0_dp-eta))*(((qc+b_u)**(1.0_dp-eta))-(b_u**(1.0_dp-eta)))
    end function Util_DM


    real(dp) function Dq_DM(qc,b_u,eta)
        real(dp),intent(in) :: qc,b_u,eta
        Dq_DM = (qc+b_u)**(-eta)
    end function Dq_DM


    real(dp) function Cost_HDM(qh,B,nu)
        real(dp),intent(in) :: qh,B,nu
        Cost_HDM = -B*(qh**nu)
    end function Cost_HDM


    real(dp) function Util_CM(c,h,gamma,kappa,chi)
        real(dp),intent(in) :: c,h,gamma,kappa,chi
        if (abs((gamma-1.0_dp) < 1.0e-4_dp)) then
            Util_CM = log(c) - (kappa*(h**(1.0_dp+chi))/(1.0_dp+chi))
        else
            Util_CM = (1.0_dp/(1.0_dp-gamma))*(c**(1.0_dp-gamma)) - (kappa*(h**(1.0_dp+chi))/(1.0_dp+chi))
        end if
    end function Util_CM


    real(dp) function Cost_HCM(H,D)
        real(dp),intent(in) :: H,D
        Cost_HCM = D*LOG(1.0-H)
    end function Cost_HCM


    real(dp) function MDisu_HDM(qh,B,nu)
        real(dp),intent(in) :: qh,B,nu
        MDisu_HDM = nu*B*(qh**(nu-1.0_dp))
    end function MDisu_HDM


end module tools_fnutil