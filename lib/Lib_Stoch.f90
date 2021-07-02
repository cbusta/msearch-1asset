module lib_stoch

! Module to implement some numerical routines that require probability distributions.
! By: Christian Bustamante
! Email: mail@cbustamante.co
! Last modified: 03 Dic, 2019, 17:17
! v 1.60
!
! To compile this module:
! $ ifort -c Lib_Stoch.f90 Lib_Basic.f90 Lib_Kind.f90
! $ gfortran -c Lib_Stoch.f90 Lib_Basic.f90 Lib_Kind.f90
! To use it put this just after the -program- statement:
! use lib_stoch
!
! Subroutines:
! (Stoch) 01. cdfnor        - Cummulative distribution function of the standard normal distribution
! (Stoch) 02. paretocdf     - Cummulative distribution function of th Pareto distrubution
! (Stoch) 03. erf           - Error funcion (used by -cdfnor-)
! (Stoch) 04. tauchen       - Tauchen's method to discretize an AR(1)
! (Stoch) 05. rouwenhorst   - Rouwenhorst's method to discretize an AR(1). Works better for persitent processes (requires less grid points)
! (Stoch) 06. markovss      - Computes the ergodic distribution of a markov process (uses eigecvectors)
! (Stoch) 07. markovss_old  - Computes the ergodic distribution of a markov process (inefficient)
! (Stoch) 08. limitdist     - Computes the ergodic distribution of a markov process (faster!)
! (Stoch) 09. simul_markov  - Simulates a Markov chain
! (Stoch) 10. ranpareto     - Generates random numbers from a Pareto distribution
! (Stoch) 11. set_seed      - Sets the seed of -random_number- with a single index (integer)
!
! Note: depends on Lib_Kind and on some routines on Lib_Basic.f90

use lib_kind
use lib_basic
implicit none
private
public :: cdfnor, erff, tauchen, rouwenhorst, markovss, markovss_old, limitdist, simul_markov, &
          ranpareto, paretocdf, set_seed, paretodist

contains


    function cdfnor(X,mu,sig)
        ! Normal cumulative distribution function (mean mu and standard deviation sigma)
        ! By J. R. M. Hosking (IBM Research Division)
        ! Version 3, August 1996
        ! Usage:
        ! y = cdfnor(X,mu,sigma)
        ! Inputs:
        ! X      - point at which cumulative probability is computed
        ! mu     - mean of the distribution
        ! sigma  - standard deviation of the distribution
        ! Output:
        ! y      - cumulative probability, Pr(z<=X)
        real(dp) :: X,mu,sig,HALF,CDFNOR,RTHALF
        data HALF/0.5D0/,RTHALF/0.707106781186547524D0/
        cdfnor = HALF+HALF*erff((X-mu)/sig*RTHALF)
        return
    end function cdfnor


    function erff(X)
        ! Error function (encountered in integrating the normal distribution)
        ! By  J. R. M. Hosking (IBM Research Division)
        ! Version 3, August 1996
        ! Based on algorithm 5666, J.F.Hart et al. (1968) 'Computer  Approximations'
        ! Accurate to 15 decimal places
        real(dp) :: ZERO,ONE,TWO,THREE,FOUR,P65
        real(dp) :: P0,P1,P2,P3,P4,P5,P6,Q0,Q1,Q2,Q3,Q4,Q5,Q6,Q7,C1,C2,BIG
        real(dp) :: XX,EXPNTL,ZZ,ERFF,X,CRIT
        data ZERO/0D0/,ONE/1D0/,TWO/2D0/,THREE/3D0/,FOUR/4D0/,P65/0.65D0/
        ! Coefficients of rational-function approximation
        data P0,P1,P2,P3,P4,P5,P6/0.2202068679123761D3,0.2212135961699311D3,0.1120792914978709D3,0.3391286607838300D2, &
                                  0.6373962203531650D1,0.7003830644436881D0,0.3526249659989109D-1/
        data Q0,Q1,Q2,Q3,Q4,Q5,Q6,Q7/0.4404137358247522D3,0.7938265125199484D3,0.6373336333788311D3,0.2965642487796737D3, &
                                  0.8678073220294608D2,0.1606417757920695D2,0.1755667163182642D1,0.8838834764831844D-1/
        ! C1 is sqrt(2), C2 is sqrt(2/PI)
        ! BIG is the point at which erff=1 to machine precision
        data C1/1.414213562373095D0/
        data C2/7.978845608028654D-1/
        data BIG/6.25D0/,CRIT/5D0/
        erff = ZERO
        if (abs(X) < 1.0D-7) return
        XX = abs(X)
        if (XX > BIG) goto 20
        EXPNTL = exp(-X*X)
        ZZ = abs(X*C1)
        if (XX > CRIT) goto 10
        erff = EXPNTL*((((((P6*ZZ+P5)*ZZ+P4)*ZZ+P3)*ZZ+P2)*ZZ+P1)*ZZ+P0 ) &
             / (((((((Q7*ZZ+Q6)*ZZ+Q5)*ZZ+Q4)*ZZ+Q3)*ZZ+Q2)*ZZ+Q1)*ZZ+Q0)
        if (X > ZERO) erff = ONE-TWO*erff
        if (X < ZERO) erff = TWO*erff-ONE
        return
    10  erff = EXPNTL*C2/(ZZ+ONE/(ZZ+TWO/(ZZ+THREE/(ZZ+FOUR/(ZZ+P65)))))
        if(X > ZERO) erff=ONE-erff
        if(X < ZERO) erff=erff-ONE
        return
    20  erff = ONE
        if(X < ZERO) erff = -ONE
        return
    end function erff
  
  
    subroutine tauchen(muz,sige,rho,width,Nx, X,Pi) 
        ! Implements the Tauchen's method to discretize an AR(1) 
        ! process with a Markov chain by computing the transition
        ! probabilities.
        ! Usage: 
        ! call tauchen(muz,sige,rho,width,Nx,X,Pi)
        ! Inputs:
        ! muz    - unconditional mean of z
        ! sige   - standard deviation of innovations to z
        ! rho    - persistence of stochastic process z
        ! width  - dist. grid's bounds for log(z): [mu-width*sigz,mu+width*sigz]
        ! Nx     - the number of desired grid points
        ! Output:
        ! X      - grid for the stochastic process log(z)
        ! Pi     - transition probabilities
        integer              :: i,j
        real(dp)             :: Phi1,Phi2,w,sigz,x1,xN,czer,ci
        integer,intent(in)   :: Nx
        real(dp),intent(in)  :: muz,sige,rho,width
        real(dp),intent(out) :: X(:),Pi(:,:)
        czer = 0.0_dp
        sigz = sqrt(sige**2/(1.0_dp-rho**2))
        x1   = muz - sigz*width
        xN   = muz + sigz*width
        call linspace(x1,xN,Nx,X)
        w = (X(Nx)-X(1))/real(Nx-1,dp)
        do i = 1,Nx
            ci = rho*X(i) + (1-rho)*muz
            Pi(i,1) = cdfnor(X(1)-ci+(w/2.0_dp),czer,sige)
            do j = 2,Nx-1
                Phi1 = cdfnor(X(j)-ci+(w/2.0_dp),czer,sige)
                Phi2 = cdfnor(X(j)-ci-(w/2.0_dp),czer,sige)
                Pi(i,j) = Phi1 - Phi2
            end do
            Pi(i,Nx) = 1.0_dp - cdfnor(X(Nx)-ci-(w/2.0_dp),czer,sige)
        end do
    end subroutine tauchen
  
  
    subroutine rouwenhorst(rho,sige,Nz, Z,Pi)
        ! Implements the Rouwenhorst's method to discretize an AR(1)
        ! process with a Markov chain by computing its transition
        ! probabilities.
        ! Usage:
        ! call rouwenhorst(rho,sige,Nz,Z,Pi)
        ! Inputs:
        ! rho    - persistence of stochastic process z
        ! sige   - standard deviation of innovations to z
        ! Nz     - the number of desired grid points
        ! Output:
        ! Z      - grid for the stochastic process z
        ! Pi     - transition probabilities
        integer,intent(in)  :: Nz
        real(dp),intent(in) :: rho,sige
        real(dp),intent(out):: Z(Nz),Pi(Nz,Nz)
        integer             :: i
        real(dp)            :: sigz,p,q,psiz,hlag(Nz,Nz),h(Nz,Nz)
        sigz = sqrt(sige**2/(1.0_dp-rho**2))
        psiz = sqrt(real(Nz,dp)-1.0_dp)*sigz
        call linspace(-psiz,psiz,Nz,Z)
        p = (1.0_dp+rho)/2.0_dp
        q = p
        hlag = 0.0_dp
        hlag(1,1) = 1.0_dp
        do i = 2,Nz
            h = 0.0_dp
            h(1:i-1,1:i-1) = p*hlag(1:i-1, 1:i-1)
            h(1:i-1,2:i)   = h(1:i-1,2:i) + (1.0-p)*hlag(1:i-1, 1:i-1)
            h(2:i,1:i-1)   = h(2:i,1:i-1) + (1.0-q)*hlag(1:i-1, 1:i-1)
            h(2:i,2:i)     = h(2:i,2:i) + q*hlag(1:i-1, 1:i-1)
            h(2:i-1,1:i)   = h(2:i-1,1:i)/2.0_dp
            hlag(1:i, 1:i) = h(1:i, 1:i)
        end do
        Pi = hlag
    end subroutine rouwenhorst

  
    subroutine markovss(Pi,N, Pi_Erg)
        ! Computes the invariant ditribution of a Markov chain
        ! Uses the power method (with scaling) to approximate eigen values.
        ! Usage: 
        ! call markovss(Pi,N, Pi_Erg)
        ! Inputs:
        ! Pi     - Markov transition probability matrix
        ! N      - number of states
        ! Output:
        ! Pi_Erg - matrix of invariant distribution
        real(dp),intent(in)  :: Pi(:,:)
        integer,intent(in)   :: N
        real(dp),intent(out) :: Pi_Erg(N)
        real(dp)             :: V(N),a
        call findmaxeig(transpose(Pi),1000, V,a)
        Pi_Erg = V/sum(V)
    end subroutine markovss
    
    
    subroutine markovss_old(Pi,N, Pi_Ss)
        ! Computes the invariant ditribution of a Markov chain
        ! [This is clearly not the most efficient way of doing it]
        ! Usage:
        ! call markovss_old(Pi,N, Pi_Ss)
        ! Inputs:
        ! Pi     - Markov transition probability matrix
        ! N      - number of states
        ! Output:
        ! Pi_Ss  - matrix of invariant distribution
        real(dp),intent(in)  :: Pi(:,:)
        integer,intent(in)   :: N
        real(dp),intent(out) :: Pi_Ss(N)
        integer              :: i
        real(dp)             :: Pi_Temp(N,N)
        Pi_Temp = Pi
        do i = 1,1000
            Pi_Temp = matmul(Pi_Temp,Pi)
        end do
        do i = 1,N
            Pi_Ss(i) = Pi_Temp(i,i)
        end do
    end subroutine markovss_old
  
  
    subroutine limitdist(Pi,N, Pi_Lim)
        ! Computes the stationary ditribution of a Markov chain
        ! Uses the state-space reduction method to compute it.
        ! Note: More efficient than markovss_old
        ! Usage:
        ! call limitdist(Pi,N, Pi_Lim)
        ! Inputs:
        ! Pi     - Markov transition probability matrix
        ! N      - number of states
        ! Output:
        ! Pi_Lim  - matrix of invariant distribution
        real(dp),intent(in)  :: Pi(:,:)
        integer,intent(in)   :: N
        real(dp),intent(out) :: Pi_Lim(N)
        integer              :: n0,n1,n2,j,j1
        real(dp)             :: Pi0(N,N),s
        Pi_Lim = 0.0_dp
        Pi0 = Pi
        n0  = N
        do while (n0 > 1)
            n1 = n0-1
            s  = sum(Pi0(n0,1:n1))
            Pi0(1:n1,n0) = Pi0(1:n1,n0)/s
            n2 = n1
            do while (n2 > 0)
                Pi0(1:n1,n2) = Pi0(1:n1,n2) + Pi0(1:n1,n0)*Pi(n0,n2)
                n2 = n2-1
            end do
            n0 = n0-1
        end do
        ! Backtracking
        Pi_Lim(1) = 1
        j = 2
        do while (j <= N)
            j1 = j-1
            Pi_Lim(j) = dot_product(Pi_Lim(1:j1),(Pi0(1:j1,j)))
            j = j+1
        end do
        Pi_Lim = Pi_Lim/(sum(Pi_Lim))
    end subroutine limitdist

  
    subroutine simul_markov(Pi,z0,T,rand, zindex)
        ! Simulates a Markov chain.
        ! Usage:
        ! call simul_markov(Pi,z0,T,rand, zindex)
        ! Inputs:
        ! Pi     - transition matrix
        ! z0     - initial state
        ! T      - number of period to simulate
        ! rand   - random seed if rand=1
        ! Output:
        ! zindex - simulated indices for z
        real(dp),intent(in)         :: Pi(:,:)
        integer,intent(in)          :: z0,T
        integer,optional,intent(in) :: rand
        real(dp),allocatable        :: cumPi(:,:),rowpi(:),diff(:),shocks(:)
        integer                     :: Nz,i,time,numpos
        integer,intent(out)         :: zindex(T)
        integer                     :: n = 12, clock, ic
        integer                     :: s, counter
        integer,dimension(:),allocatable :: iseed
        Nz = size(Pi,1)
        allocate(cumPi(Nz,Nz))
        allocate(rowpi(Nz))
        allocate(diff(Nz))
        cumPi(:,1) = Pi(:,1)
        do i = 2,Nz
            cumPi(:,i) = cumPi(:,i-1) + Pi(:,i)
        end do
        zindex = z0
        allocate(shocks(T))
        if (rand == 1) then
            ! Choose random seed from system's clock
            allocate(iseed(n))
            call random_seed(size = n)
            call system_clock(count = clock)
            iseed = clock + 37 * [(ic, ic = 0,n-1)]
            call random_seed(put = iseed)
        else
            ! No random seed selected: seed will be 1
            call random_seed(size = s)
            allocate(iseed(s))
            iseed = 1
            call random_seed(put = iseed)
        end if
        call random_number(shocks)
        do time = 1,T-1
            rowpi   = cumPi(zindex(time),1:Nz)
            diff    = rowpi - shocks(time)
            counter = 0
            do i = 1,Nz
                if (diff(i) > 0.0_dp) counter = counter+1
            end do
            numpos  = counter
            zindex(time+1) = size(Pi,1)+1-numpos
        end do
    end subroutine simul_markov
  
  
    subroutine ranpareto(xm,a,nr,rand, rvec)
        ! Generates a random number/vector from the
        ! Pareto distrubution with parameters (xm,a)
        ! Usage:
        ! call ranpareto(xm,a,nr,rand, rvec)
        ! Inputs:
        ! xm     - scale parameter (lower bound)
        ! a      - shape parameter (curvature)
        ! nr     - dimension of random vector to be generated
        ! rand   - random seed if rand=1
        ! Output:
        ! rvec   - pareto random vector
        real(dp),intent(in)  :: xm,a
        integer,intent(in)   :: rand,nr
        real(dp),intent(out) :: rvec(nr)
        real(dp)             :: uvec(nr)
        integer              :: n = 12, clock, ic
        integer              :: s
        integer,dimension(:),allocatable :: iseed
        if (rand == 1) then
            allocate(iseed(n))
            call random_seed(size = n)
            call system_clock(count = clock)
            iseed = clock + 37 * [(ic, ic = 0,n-1)]
            call random_seed(put = iseed)
        else
            call random_seed(size = s)
            allocate(iseed(s))
            iseed = 1 ! no random seed selected: seed will be 1
            call random_seed(put = iseed)
        end if
        call random_number(uvec)
        rvec = xm / (uvec**(1.0_dp/a))
    end subroutine ranpareto
  
  
    function paretocdf(xm, a, xvalue) result (prob)
        ! Computes the CDF at xvalue of a Pareto distrubution with
        ! parameters (xm,a)
        ! Usage:
        ! prob = paretocdf(rnum,xm,a,rand)
        ! Inputs:
        ! xm     - scale parameter (lower bound)
        ! a      - shape parameter (curvature)
        ! Output:
        ! prob   - pareto random number
        real(dp),intent(in) :: xm, a, xvalue
        real(dp)            :: xval1, prob
        xval1 = max(xvalue, xm)
        prob  = 1.0_dp - (xm/xval1)**a
    end function paretocdf
  
  
    subroutine set_seed(index)
        ! Determines the dimension of the seed vector and sets the seed of
        ! the random number generator in a deterministic way.
        ! The function must be called before calling RANDOM_NUMBER
        ! Usage:
        ! call set_seed(index)
        ! Inputs:
        ! index  - choosen seed
        integer, intent(in) :: index
        integer             :: i, n
        integer, dimension(:), allocatable :: seed
        call random_seed(size = n)
        allocate(seed(n))
        seed = index + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(put = seed)
    end subroutine set_seed
    
    
    subroutine paretodist(em, alpha, evec, measure)
        ! Discretised pareto distribution with parameters em, the lower
        ! bound, and curvature alpha. The discrete supports evec has neps size.
        integer              :: neps
        real(dp),intent(in)  :: em, alpha, evec(:)
        integer              :: neps2, i
        real(dp)             :: evec2(size(evec)-1), prob, sumprob
        real(dp),intent(out) :: measure(size(evec))
        ! Create a grid of interlaced points
        neps  = size(evec)
        neps2 = neps-1
        evec2 = 0.0_dp
        do i = 1,neps2
            evec2(i) = (evec(i) + evec(i+1))/2.0_dp
        end do
        ! Allocate probability at endspoints of the evec using the cdf at the
        ! nearest interior point.
        prob = paretocdf(em, alpha, evec2(neps2))
        measure(neps) = 1.0_dp - prob
        prob = paretocdf(em, alpha, evec2(1))
        measure(1) = prob
        sumprob = measure(1)
        do i = 2, neps-1
            prob = paretocdf(em, alpha, evec2(i))
            measure(i) = prob - sumprob
            sumprob = sumprob + measure(i)
        end do
    end subroutine paretodist

  
end module lib_stoch