module lib_stats

! Module to implement some statistical functions.
! By: Christian Bustamante
! Email: mail@cbustamante.co
! Last modified: 27 Nov, 2019, 20:46
! v 1.30
!
! To compile this module:
! $ ifort -c Lib_Stats.f90 Lib_Basic.f90 Lib_Kind.f90
! $ gfortran -c Lib_Stats.f90 Lib_Basic.f90 Lib_Kind.f90
! To use it put this just after the -program- statement:
! use lib_stats
!
! Subroutines:
! (Basic) 01. mean          - Computes the mean of a vector
! (Basic) 02. stdev         - Computes the standard deviation of a vector
! (Basic) 03. variance      - Computes the variance of a vector
! (Basic) 04. cov           - Computes the covariance between two vectors
! (Basic) 05. corr          - Computes the correlation coefficient between to vectos
! (Basic) 06. ols           - Computes the estimation of the OLS coefficients
! (Basic) 07. wgt_pctile    - Computes the specified percentiles for weighted data
!
! Note: depends on Lib_Kind and Lib_Basic

use lib_kind
use lib_basic
implicit none
private
public :: mean, stdev, variance, cov, corr, ols, wgt_pctile

contains


    function mean(X) result(mu)
        ! Computes the mean of a vector
        ! Usage:
        ! mu = mean(X)
        ! Inputs:
        ! X     - a one dimensional vector
        ! Output:
        ! mu    - mean of X
        real(dp),intent(in)  :: X(:)
        real(dp)             :: mu
        integer              :: n
        n   = size(X)
        mu  = sum(X)/n
    end function mean


    function stdev(X) result(sig)
        ! Computes the standard deviation of a vector
        ! Usage:
        ! sig = stdev(X)
        ! Inputs:
        ! X     - a one dimensional vector
        ! Output:
        ! sig   - standard deviation of X
        real(dp),intent(in)  :: X(:)
        real(dp)             :: sig
        real(dp)             :: mu
        integer              :: n
        n   = size(X)
        mu  = sum(X)/n
        sig = sqrt(sum((X-mu)**2)/(n-1.0_dp))
    end function stdev
  
  
    function variance(X) result(sig2)
        ! Computes the variance of a vector
        ! Usage:
        ! sig = variance(X)
        ! Inputs:
        ! X     - a one dimensional vector
        ! Output:
        ! sig2  - standard deviation of X
        real(dp),intent(in)  :: X(:)
        real(dp)             :: sig2
        real(dp)             :: mu
        integer              :: n
        n    = size(X)
        mu   = sum(X)/n
        sig2 = sum((X-mu)**2)/(n-1.0_dp)
    end function variance
  
  
    function cov(X,Y) result(covxy)
        ! Computes the covariance between two vectors
        ! Usage:
        ! covxy = cov(X,Y)
        ! Inputs:
        ! X     - a one dimensional vector of size n
        ! Y     - a one dimensional vector of size n
        ! Output:
        ! covxy - covariance between X and Y
        real(dp),intent(in)  :: X(:),Y(:)
        real(dp)             :: covxy
        real(dp)             :: mux,muy
        integer              :: n
        n = size(X)
        if (size(Y) /= n) print *, 'Error (cov): vectors do not have the same length'
        mux   = sum(X)/n
        muy   = sum(Y)/n
        covxy = sum((X-mux)*(Y-muy))/(n-1.0_dp)
    end function cov
  
  
    function corr(X,Y) result(corrxy)
        ! Computes the corelation coefficient between two vectors
        ! Usage:
        ! corrxy = corr(X,Y)
        ! Inputs:
        ! X      - a one dimensional vector of size n
        ! Y      - a one dimensional vector of size n
        ! Output:
        ! corrxy - correlation between X and Y
        real(dp),intent(in)  :: X(:),Y(:)
        real(dp)             :: corrxy
        real(dp)             :: covxy,stdx,stdy
        integer              :: n
        n = size(X)
        if (size(Y) /= n) print *, 'Error (corr): vectors do not have the same length'
        covxy  = cov(X,Y)
        stdx   = stdev(X)
        stdy   = stdev(Y)
        corrxy = covxy/(stdx*stdy)
    end function corr
  
  
    subroutine ols(X,Y,Cons, Beta,Stde,R2)
        ! Computes the coefficients of the OLS estimation of Y=Beta'X
        ! Usage:
        ! call ols(X,Y,Cons, Beta,Stde,R2)
        ! Inputs:
        ! X     - matrix with the independent variables observations
        ! Y     - vector with the dependent variable observations
        ! Cons  - if equal to 1, include an intercept in the regression
        ! Output:
        ! Beta  - estimated OLS coefficients
        ! Stde  - variance covariance matrix of beta
        ! R2    - r-squared
        real(dp),intent(in)              :: X(:,:),Y(:,:)
        integer,intent(in)               :: Cons
        real(dp),allocatable,intent(out) :: Beta(:,:),Stde(:,:)
        real(dp),intent(out)             :: R2
        integer                          :: rowx,rowy,colx,coly,n,k
        real(dp),allocatable             :: X_Mat(:,:),Yhat(:,:),Resid(:,:)
        real(dp)                         :: Ybar,sighat,SH(1,1),ess,tss
        rowx = size(X,1)
        colx = size(X,2)
        rowy = size(Y,1)
        coly = size(Y,2)
        n    = rowx
        if (rowx /= rowy) print *, 'Error (ols): Number of observarions doesnt match'
        if (coly /= 1)    print *, 'Error (ols): Wrong number of columns in Y'
        if (Cons==1) then
            k = colx+1
            allocate(X_Mat(n,k))
            X_Mat(:,1) = 1.0_dp
            X_Mat(:,2:k) = X
        else
            k = colx
            allocate(X_Mat(n,k))
            X_Mat(:,1:k) = X
        end if
        allocate(Beta(k,1))
        allocate(Stde(k,k))
        allocate(Yhat(n,1))
        allocate(Resid(n,1))
        Ybar   = sum(Y)/n
        Beta   = matmul(inv(matmul(transpose(X_Mat),X_Mat)), &
                 matmul(transpose(X_Mat),Y))
        Yhat   = matmul(X_Mat,Beta)
        ess    = sum((Yhat-Ybar)**2)
        tss    = sum((Y-Ybar)**2)
        Resid  = Y - Yhat
        SH     = sum(Resid**2)/(n-k)
        sighat = SH(1,1)
        Stde   = sighat*inv(matmul(transpose(X_Mat),X_Mat))
        R2     = ess/tss
    end subroutine ols
  
  
    subroutine wgt_pctile(x,wgt,P, pctl)
        ! Computes the specified percentiles for weighted data
        ! Usage:
        ! call wgt_pctile(x,wgt,P, pctl)
        ! Inputs:
        ! x      - vector with the data
        ! wgt    - vector with the weights (same length as x)
        ! P      - scalar with the percentile to compute in [0,1]
        ! Output:
        ! pctl   - value of percentile
        real(dp),intent(in) :: x(:),wgt(:),P
        real(dp),intent(out):: pctl
        integer             :: nX,k
        real(dp)            :: wscal(size(x,1))
        real(dp)            :: wsnord(size(x,1),2),A(size(x,1),2)
        real(dp)            :: xs(size(x,1)),ws(size(x,1))
        real(dp)            :: cumws(size(x,1))
        real(dp)            :: wint,xmid
        nX    = size(x,1)
        wscal = wgt/sum(wgt)
        wsnord(:,1) = x
        wsnord(:,2) = wscal
        call sort_rows(wsnord,1,A)!? 1 vs 2?
        xs   = A(:,1)
        ws   = A(:,2)
        call cumsum(ws,cumws)
        pctl = 0.0_dp
        k = gridlookup(cumws,P)
        if (cumws(k) == P) k = k-1
        if ((k < nX) .and. (k > 1)) then
            wint = (P-cumws(k))/(cumws(k+1)-cumws(k))
            xmid = (xs(k+1)-xs(k))
            pctl = xs(k) + wint * xmid
        elseif (k == nX) then
            pctl = xs(k)
        elseif (k == 0) then
            pctl = xs(1)
        end if
    end subroutine wgt_pctile
  
  
end module lib_stats