module lib_solver

! Module to implement some root finding and minimization routines.
! By: Christian Bustamante
! Email: mail@cbustamante.co
! Last modified: 13 Mar 2020, 16:15
! v 2.0
!
! To compile this module:
! $ ifort -c Lib_Solver.f90 Par_Pass.f90 Lib_Kind.f90 Lib_Rwtxt.f90
! $ gfortran -c Lib_Solver.f90 Par_Pass.f90 Lib_Kind.f90 Lib_Rwtxt.f90
! To use it put this just after the -program- statement:
! use lib_solver
!
! Subroutines:
! (Solver) 01. bisect     - Univariate bisection method to find a zero
! (Solver) 02. nl1secant  - Univariate implementation of the root-finding secant method
! (Solver) 03. hybisecant - Hybrid bisection-secant method for univariate root-finding problems
! (Solver) 04. broyden    - Broyden's method to solve a system of non-linear equations (w/ backtrack and restarting possible)
! (Solver) 05. hybroyden  - Powell's hybrid method to solve systems of non-linear equations
! (Solver) 06. fdjac      - Computes the numerical derivative of a system of equations 
! (Solver) 07. golden     - Univariate golden section search method to find a min
! (Solver) 08. amoeba     - Performs the Nelder-Mead optimization (minimization) search
!
! Note 1: depends on Lib_Kind and Lib_Rwtxt (just to export results), AND...
! Note 2: requires to declare the module -par_pass- defining the following derived data types (they can be empty):
!   parpass_g  (Golden section search)
!   parpass_nl (Non-linear equations - univariate and multivariate)
!   parpass_a  (Nelder-Mead minimization)

use lib_kind
use lib_basic
use lib_rwtxt
use par_pass
use lib_smalapack
implicit none
private
public :: bisect, nl1secant, hybisecant, broyden, hybroyden, fdjac, fdjac_par, golden, amoeba

contains
    

    subroutine bisect(xl,xr,tol,FN,PARAM, x_sol,f_sol)
        ! Implements the bisection method.
        ! Usage:
        ! call bisect(xl,xr,tol,FN,PARAM, x_sol,f_sol)
        ! Inputs:
        ! xl,xr - left and right edges of the search interval
        ! tol   - stoping criteria
        ! FN    - name of the function to be solved, must be a (scalar) double precision function
        ! PARAM - type(parpass_nl) structure with the parameters of FN
        ! Output:
        ! x_sol - solution to FN(X,PARAM)=0
        ! f_sol - value of FN at x_sol
        ! Note:
        ! FN must have a dummy argument PARAM.
        ! Requires module -par_pass- defining the derived data type -parpass_nl-
        real(dp), intent(inout)       :: xl,xr
        real(dp), intent(in)          :: tol
        type(parpass_nl),intent(inout):: PARAM
        real(dp), intent(out)         :: x_sol,f_sol
        real(dp)                      :: diff,xm,fl,fm,fr
        interface
            real(dp) function FN(y,PARAM)
                use par_pass
                real(dp), intent(in)           :: y
                type(parpass_nl), intent(inout):: PARAM
            end function FN
        end interface
        diff = 2.0_dp*tol
        xm   = (xl+xr)/2.0_dp
        fl   = FN(xl,PARAM)
        fr   = FN(xr,PARAM)
        if (fl*fr > 0.0_dp) then
            print *, "Error (bisect): your bounds do not bracket a root"
            print *, 'x', xl, xr
            print *, 'f', fl, fr
            return
        end if
        fm   = FN(xm,PARAM)
        do while(diff > tol)
            if (fl*fm < 0.0_dp) then
                xr = xm
            else
                xl = xm
                fl = fm
            end if
            diff = abs(xl-xr)
            xm   = (xl+xr)/2.0_dp
            fm   = FN(xm,PARAM)
        end do
        x_sol = xm
        f_sol = fm
    end subroutine bisect
  
  
    subroutine nl1secant(x0,x1,tol,FN,PARAM, x_sol,f_sol)
        ! Implements the secant method for solving one-dimensional root
        ! finding problems
        ! Usage:
        ! call nl1secant(x0,x1,tol,FN,PARAM, x_sol,f_sol)
        ! Inputs:
        ! x0,x1 - initial guesses for x, with x0/=x1 (inout)
        ! tol   - stoping criteria
        ! FN    - name of the function to be solved, must be a (scalar) double precision function
        ! PARAM - type(parpass_nl) structure with the parameters of FN
        ! Output:
        ! x_sol - solution to FN(X,PARAM)=0
        ! f_sol - value of FN at x_sol
        ! Note:
        ! FN must have a dummy argument PARAM.
        ! Requires module -par_pass- defining the derived data type -parpass_nl-
        real(dp), intent(inout)       :: x0,x1
        real(dp), intent(in)          :: tol
        type(parpass_nl),intent(inout):: PARAM
        real(dp), intent(out)         :: x_sol,f_sol
        real(dp)                      :: diff,d,x,f,f0,f1
        interface
            real(dp) function FN(y,PARAM)
                use par_pass
                real(dp), intent(in)           :: y
                type(parpass_nl), intent(inout):: PARAM
            end function FN
        end interface
        !Initializing
        f0 = FN(x0,PARAM)
        f1 = FN(x1,PARAM)
        ! Secant method approximates the derivate of f at xk with the slope
        ! between xk and xk-1
        diff = 2.0_dp*tol
        do while (diff > tol)
            d = (f1-f0)/(x1-x0)
            x = x1 - f1/d
            f = FN(x,PARAM)
            ! Convergence criterion
            diff = abs(x-x1)
            ! Setting up next step
            x0 = x1
            x1 = x
            f0 = f1
            f1 = f
        end do
        x_sol = x1
        f_sol = f1
    end subroutine nl1secant

  
    subroutine hybisecant(xl,xr,tol,FN,PARAM, x_sol,f_sol)
        ! Implements a hybrid secant - bisection method for solving
        ! one-dimensional root finding problems
        ! It approximates the derivative using a secant and reverts to
        ! bisection if progress is insufficient.
        ! Bounds that bracket a zero are required
        ! Usage:
        ! call hybisecant(x0,x1,tol,FN,PARAM, x_sol,f_sol)
        ! Inputs:
        ! xl,xr - left and right edges of the search interval (inout)
        ! tol   - stoping criteria
        ! FN    - name of the function to be solved, must be a (scalar) double precision function
        ! PARAM - type(parpass_nl) structure with the parameters of FN
        ! Output:
        ! x_sol - solution to FN(X,PARAM)=0
        ! f_sol - value of FN at x_sol
        ! Note:
        ! FN must have a dummy argument PARAM.
        ! Requires module -par_pass- defining the derived data type -parpass_nl-
        real(dp), intent(inout)       :: xl,xr
        real(dp), intent(in)          :: tol
        type(parpass_nl),intent(inout):: PARAM
        real(dp), intent(out)         :: x_sol,f_sol
        integer                       :: iters
        integer,parameter             :: maxiters = 10
        real(dp)                      :: diffs,d
        real(dp)                      :: f,fl,fr,f0,flp,frp
        real(dp)                      :: x,x0,xrp,xlp
        interface
            real(dp) function FN(y,PARAM)
                use par_pass
                real(dp), intent(in)           :: y
                type(parpass_nl), intent(inout):: PARAM
            end function FN
        end interface
        !Initializing
        fl = FN(xl,PARAM)
        fr = FN(xr,PARAM)
        ! Is the zero bracketed?
        if (fl*fr > 0.0_dp) then
            if (abs(fl) < abs(fr)) then
                x = xl
                f = fl
            else
                x = xr
                f = fr
            end if
        else
            ! Hybrid method
            d   = (fr-fl)/(xr-xl)
            x0  = xr
            f0  = fr
            xlp = xl
            xrp = xr
            flp = fl
            frp = fr
            ! Using secant method (up to a max number of iter)
            iters = 0
            diffs = abs(xrp-xlp)
            do while ((diffs > tol) .and. (abs(f0) > tol) .and. (iters < maxiters))
                iters = iters + 1
                x = x0 - f0/d
                if ((x < xl) .or. (x > xr)) x = (xl+xr)/2.0_dp
                f  = FN(x,PARAM)
                d  = (f-f0)/(x-x0)
                x0 = x
                f0 = f
                if (f0*flp > 0.0_dp) then
                    xlp = x
                    flp = f
                else
                    xrp = x
                    frp = f
                end if
                diffs = abs(xrp-xlp)
            end do
            ! if secant dind't work, use bisection
            if ((diffs > tol) .and. (abs(f0) > tol) .and. (iters >= maxiters)) then
                call bisect(xlp,xrp,tol,FN,PARAM, x,f)
            end if
        end if
        x_sol = x
        f_sol = f
    end subroutine hybisecant
  
  
    subroutine broyden(x0,invD0,tol,FN,PARAM,backtrack,restart,  x_sol,f_sol,invD_sol)
        ! Solves a system of non-linear equations using Broyden's method with 
        ! an initial guess for the inverse of the Jacobian matrix. It also 
        ! allows for line search backtracking, and Jacobian restarting
        ! Parts of this subroutine are adapted from -broydn- in "Numerical 
        ! Recipies in Fortran 90"
        ! Usage:
        ! call broyden(x0,invD0,tol,FN,PARAM,backtrack,restart, x_sol,f_sol,invD_sol)
        ! Inputs:
        ! x0        - initial value for x
        ! invD      - initial guess for the inverse of Jacobian
        ! tol       - stoping criteria
        ! FN        - name of the subroutine with the function to be solved
        ! PARAM     - type(PARPASS_NL) structure with the parameters of FN
        ! backtrack - set to .true. for line search backtrack (if not, Broyden takes a full Newton step)
        ! restart   - set to .true. if allowing for Jacobian restarting
        ! Output:
        ! x_sol     - solution to FN(x,PARAM)=0
        ! f_sol     - value of FN(x_opt,PARAM)
        ! invD_sol  - current inverse Jacobian upon exit
        ! Note:
        ! FN must have a dummy argument PARAM.
        ! Requires module -par_pass- defining the derived data type -parpass_nl-
        real(dp),intent(in)           :: x0(:),invD0(:,:),tol
        type(parpass_nl),intent(inout):: PARAM
        logical, intent(in)           :: backtrack,restart
        real(dp),intent(out)          :: x_sol(size(x0)),f_sol(size(x0))
        real(dp),intent(out)          :: invD_sol(size(x0),size(x0))
        real(dp)                      :: invDk(size(x0),size(x0))
        real(dp)                      :: invDk1(size(x0),size(x0))
        real(dp)                      :: newDk(size(x0),size(x0))
        real(dp)                      :: xk(size(x0),1),xk1(size(x0),1)
        real(dp)                      :: dk(size(x0),1),uk(size(x0),1)
        real(dp)                      :: fk(size(x0),1),fk1(size(x0),1)
        real(dp)                      :: num(size(x0),size(x0)),den(1,1)
        real(dp)                      :: diffx,fmxabs
        integer                       :: k,maxiter
        real(dp)                      :: fknorm,fk1norm,xconv
        logical                       :: check,restnow
        interface
            subroutine FN(outf,y,PARAM)
                use lib_kind
                use par_pass
                real(dp),intent(in)           :: y(:)
                type(parpass_nl),intent(inout):: PARAM
                real(dp),intent(out)          :: outf(:)
            end subroutine FN
        end interface
        maxiter = 50
        restnow = .true.
        ! Initializing
        k = 0
        xk(:,1) = x0
        invDk   = invD0
        call FN(fk(:,1),xk(:,1),PARAM)
        ! Broyden
        diffx = 2.0_dp*tol
        print "(a)", 'Broyden'
        do while (k <= maxiter)
            k  = k + 1            
            ! Using line search with backtrack
            fknorm = 0.5_dp*dot_product(fk(:,1),fk(:,1))
            if (backtrack) then
                ! If using backtrack
                call lnsrch(xk(:,1),fk(:,1),fknorm,invDk,FN,PARAM, xk1(:,1),fk1(:,1),fk1norm,xconv,check)
            else
                ! Taking a full Newton step
                xk1 = xk - matmul(invDk,fk)
                call FN(fk1(:,1),xk1(:,1),PARAM)
                fk1norm = 0.5_dp*dot_product(fk1(:,1),fk1(:,1))
            end if            
            ! Convergence
            diffx  = maxval(abs(xk1-xk)/max(abs(xk),1.0_dp))
            fmxabs = maxval(abs(fk1))
            if (backtrack) then
                ! If using backtrack
                if (fmxabs < tol) then 
                    ! Test for convergence on function values
                    check = .false.
                    exit
                end if
                if (check) then 
                    ! True if line search failed to find a new x
                    if (restnow .or. (xconv < tol)) exit
                    ! If restnow is true we have failure: We have already tried reinitializing the Jacobian.
                    ! The other test is for gradient of f zero, i.e., spurious convergence.
                    ! Try reinitializing the Jacobian.
                    restnow = .true. 
                    if (restnow .and. restart) then
                        print *, '[Broyden]: Restarting. Recomputing Jacobian'
                        call fdjac_par(fk(:,1),xk(:,1),FN,PARAM, newDk)
                        invDk = invlap(newDk)
                    else
                        print *, '[Broyden]: No further progress can be made without restarting Jacobian'
                    end if                    
                else 
                    ! Successful step. Will use Broyden update for next step.
                    restnow = .false. 
                    ! Test for convergence on δx.
                    if (diffx < tol) exit
                end if                
            else
                ! Test for convergence on function values
                if (fmxabs < tol) exit
                ! Test for convergence on dk
                if (diffx < tol) exit
            end if

            ! Update invDk1, dk is n x 1, uk is (nxn)*(nx1)
            if ((.not.(restnow .and. restart)) .or. (.not.(backtrack))) then
                ! Do this if Jacobian has not been restarted
                dk  = xk1 - xk
                uk  = matmul(invDk, (fk1 - fk))
                den = matmul(transpose(dk),uk)
                num = matmul((dk - uk),matmul(transpose(dk),invDk))
                invDk1 = invDk + num/den(1,1)
                ! Updating for new iteration
                xk    = xk1
                invDk = invDk1
                fk    = fk1
                ! Progress
                print "(a,i5,a,f12.8,a,f12.8,a,f12.8)", 'Iter = ', k,', Diffx = ', diffx,', Norm = ', fmxabs,', SSR = ', fk1norm
            else
                print "(a)", '[Broyden]: Jacobian restarted!'
            end if
        end do
        if ((k > maxiter) .or. (check)) print "(a)", '[Broyden]: Your search has failed!'
        if (.not.(check)) then
            print "(a)", ' '
            print "(a)", '[Broyden]: Convergence achieved!'
            print "(a,f12.8,a,f12.8)", 'Root found with: Norm = ', fmxabs,', SSR = ', fk1norm
            print "(a)", ' '
        end if
        x_sol = xk1(:,1)
        f_sol = fk1(:,1)
        invD_sol = invDk
    end subroutine broyden


    subroutine lnsrch(xold,fvec_old,fnorm_old,invD,FN,PARAM, x,fvec,fnorm,xconv,check)
        ! Line search with backtracking. Adapted from -lnsrch- in "Numerical 
        ! Recipies in Fortran 90".
        ! Given an N-dimensional point xold, the value of a N-vector fuction, 
        ! and its norm (SSR/2), and the Jacobian at that same point, finds a 
        ! new point x along the Newtown direction from xold where the function 
        ! FN has decreased "sufficiently". 
        ! The output quantity check is false on a normal exit. It is true when 
        ! x is too close to xold. 
        ! Usage:
        ! call lnsrch(xold,fvec_old,fnorm_old,invD,FN,PARAM, x,fvec,fnorm,xconv,check)
        ! Inputs:
        ! xold      - initial value for x
        ! fvec_old  - function FN evaluated at xold
        ! fnorm_old - norm of fvec_old (SSR/2)    
        ! invD      - inverse of Jacobian
        ! FN        - name of the subroutine with the function to be solved
        ! PARAM     - type(PARPASS_NL) structure with the parameters of FN
        ! Output:
        ! x         - value of x that decreases sufficiently FN (max step given direction)
        ! fvec      - value of FN(x,PARAM)
        ! fnorm     - norm of fvec (SSR/2)
        ! xconv     - used to check for spurious convergence
        ! check     - value is .false. if search was succesful
        ! Note:
        ! FN must have a dummy argument PARAM.
        ! Requires module -par_pass- defining the derived data type -parpass_nl-
        real(dp),intent(in)           :: xold(:),fvec_old(size(xold)),invD(size(xold),size(xold))
        real(dp),intent(in)           :: fnorm_old
        type(parpass_nl),intent(inout):: PARAM
        real(dp),intent(out)          :: x(:),fvec(:)
        real(dp),intent(out)          :: fnorm,xconv
        logical,intent(out)           :: check
        real(dp),parameter            :: ALF  = 1.0e-4_dp
        real(dp),parameter            :: TOLX = epsilon(x)
        real(dp)                      :: g(size(xold)),p(size(xold))
        real(dp)                      :: vabs,a,alam,alam2,alamin,b,disc,fnorm2,pabs,rhs1,rhs2,slope,tmplam,stpmax
        interface
            subroutine FN(outf,y,PARAM)
                use lib_kind
                use par_pass
                real(dp),intent(in)           :: y(:)
                type(parpass_nl),intent(inout):: PARAM
                real(dp),intent(out)          :: outf(:)
            end subroutine FN
        end interface
        g(:)   = matmul(fvec_old(:),invlap(invD(:,:)))              ! Compute ∇f for the line search
        p(:)   = matmul(invD(:,:),-fvec_old(:))
        vabs   = sqrt(dot_product(xold(:),xold(:)))
        stpmax = 100.0_dp*max(vabs,real(size(xold),dp))             ! Calculate stpmax for line searches
        check  = .false.
        pabs   = sqrt(dot_product(p(:),p(:)))
        if (pabs > stpmax) p(:) = p(:)*stpmax/pabs                  ! Scale if attempted step is too big
        slope  = dot_product(g,p)
        if (slope >= 0.0_dp) print *, 'Error (lnsrc): Roundoff problem'
        alamin = TOLX/maxval(abs(p(:))/max(abs(xold(:)),1.0_dp))    ! Compute λmin
        alam   = 1.0_dp                                             ! Always try full Newton step first
        do                                                          ! Start of iteration loop
            x(:) = xold(:) + alam*p(:)
            call FN(fvec,x,PARAM)
            fnorm = 0.5_dp*dot_product(fvec(:),fvec(:))
            if (alam < alamin) then                                 ! Convergence on Δx. For zero finding, the calling program should verify the convergence
                x(:)  = xold(:)
                check = .true.
                return
            elseif (fnorm <= fnorm_old+ALF*alam*slope) then         ! Sufficient function decrease
                return
            else                                                    ! Backtrack
                if (abs(alam-1.0_dp) < TOLX) then                   ! First time
                    tmplam = -slope/(2.0_dp*(fnorm-fnorm_old-slope))
                else                                                ! Subsequent backtracks
                    rhs1 = fnorm  - fnorm_old - alam*slope
                    rhs2 = fnorm2 - fnorm_old - alam2*slope
                    a = (rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
                    b = (-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
                    if (a == 0.0_dp) then
                        tmplam = -slope/(2.0_dp*b)
                    else
                        disc = b*b-3.0_dp*a*slope
                        if (disc < 0.0_dp) then
                            tmplam = 0.5_dp*alam
                        else if (b <= 0.0) then
                            tmplam = (-b+sqrt(disc))/(3.0_dp*a)
                        else
                            tmplam = -slope/(b+sqrt(disc))
                        end if
                    end if
                    if (tmplam > 0.5_dp*alam) tmplam = 0.5_dp*alam  ! λ ≤ 0.5λ1
                end if
            end if
            alam2  = alam
            fnorm2 = fnorm
            alam   = max(tmplam,0.1_dp*alam)                        ! λ ≥ 0.1λ1
        end do        
        xconv = maxval(abs(g(:))*max(abs(x(:)), 1.0_dp)/max(fnorm,0.5_dp*real(size(xold),dp)))
    end subroutine lnsrch
  
  
    subroutine hybroyden(x0,tol,FN,PARAM, x_sol,f_sol)
        ! Solves a system of non-linear equations using the Powell's hybrid
        ! method, i.e. rejects Newton step is new value makes no progress
        ! Usage:
        ! call hybroyden(x0,tol,FN,PARAM, x_sol,f_sol)
        ! Inputs:
        ! x0     - initial value for x
        ! tol    - stoping criteria
        ! FN     - name of the subroutine with the function to be solved
        ! PARAM  - type(PARPASS_NL) structure with the parameters of FN
        ! Output:
        ! x_sol  - solution to FN(x,PARAM)=0
        ! f_sol  - value of FN(x_opt,PARAM)
        ! Note:
        ! FN must have a dummy argument PARAM.
        ! Requires module -par_pass- defining the derived data type -parpass_nl-
        real(dp),intent(in)           :: x0(:),tol
        type(parpass_nl),intent(inout):: PARAM
        real(dp),intent(out)          :: x_sol(size(x0)),f_sol(size(x0))
        real(dp)                      :: Jk(size(x0,1),size(x0,1))
        real(dp)                      :: Jk1(size(x0,1),size(x0,1))
        real(dp)                      :: xk(size(x0,1),1),xk1(size(x0,1),1)
        real(dp)                      :: sk(size(x0,1),1)
        real(dp)                      :: fk(size(x0,1),1),fk1(size(x0,1),1)
        real(dp)                      :: yk(size(x0,1),1),den(1,1)
        real(dp)                      :: SSRk,SSRk1,lambda
        real(dp)                      :: diff
        integer                       :: n,k,j,iterssr,maxitssr=5
        interface
            subroutine FN(outf,y,PARAM)
                use lib_kind
                use par_pass
                real(dp),intent(in)           :: y(:)
                type(parpass_nl),intent(inout):: PARAM
                real(dp),intent(out)          :: outf(:)
            end subroutine FN
        end interface
        n = size(x0,1)
        ! Initializing
        k = 0
        xk1(:,1) = x0
        Jk1(:,:) = 0.0_dp
        forall(j = 1:n) Jk1(j,j) = 1.0_dp
        ! Broyden
        diff = 2.0_dp*tol
        do while ((diff > tol) .and. (k <= 1000))
            k  = k + 1
            xk = xk1
            Jk = Jk1
            ! Updating x
            call FN(fk(:,1),xk(:,1),PARAM)
            sk  = (-1.0_dp)*matmul(inv(Jk),fk)
            xk1 = xk + sk
            ! Updating the value of function
            call FN(fk1(:,1),xk1(:,1),PARAM)
            ! Quadratic sum of functions
            SSRk   = sum(fk(:,1)**2)
            SSRk1  = sum(fk1(:,1)**2)
            lambda = 0.5_dp
            ! No improvement with the Newton step?
            iterssr = 0
            do while ((SSRk1 >= SSRk) .and. (iterssr < maxitssr))
                iterssr = iterssr + 1
                sk  = lambda*sk
                xk1 = xk + sk
                call FN(fk1(:,1),xk1(:,1),PARAM)
                SSRk1  = sum(fk1(:,1)**2)
                lambda = lambda/2.0_dp
            end do
            ! Updating Jacobian guess
            yk  = fk1 - fk
            den = matmul(transpose(sk),sk)
            Jk1 = Jk + (matmul((yk-matmul(Jk,sk)),transpose(sk)))/den
            diff = maxval(abs(xk-xk1))
        end do
        x_sol = xk1(:,1)
        f_sol = fk1(:,1)
    end subroutine hybroyden
  
  
    subroutine fdjac(x0,FN,PARAM, D)
        ! Computes the numerical derivative of a system of equations
        ! Usage:
        ! call fdjac(x0,FN,PARAM, D)
        ! Inputs:
        ! x0     - initial value for x
        ! FN     - name of the subroutine with a function containing an NL system of equations
        ! PARAM  - type(PARPASS_NL) structure with the parameters of FN
        ! Output:
        ! D      - Jacobian matrix with numerical derivatives
        ! Note:
        ! FN must have a dummy argument PARAM.
        ! Requires module -par_pass- defining the derived data type -parpass_nl-
        real(dp),intent(in)           :: x0(:)
        type(parpass_nl),intent(inout):: PARAM
        real(dp),intent(out)          :: D(size(x0),size(x0))
        real(dp)                      :: xnow(size(x0,1))
        real(dp)                      :: f0(size(x0,1)),fnow(size(x0,1))
        real(dp)                      :: h0(size(x0,1)),hval
        real(dp)                      :: dh(size(x0,1)),dhval
        integer                       :: n,i,j
        interface
            subroutine FN(outf,y,PARAM)
                use lib_kind
                use par_pass
                real(dp),intent(in)           :: y(:)
                type(parpass_nl),intent(inout):: PARAM
                real(dp),intent(out)          :: outf(:)
            end subroutine FN
        end interface
        n    = size(x0,1)        ! Number of variables
        call FN(f0,x0,PARAM)     ! Vector with FN evaluated at x0
        h0   = 0.0_dp            ! Used to generate a basis vector
        hval = 1.0e-6_dp         ! Step used to compute derivatives
        D    = 0.0_dp            ! Where the approximated Jacobian will be stored
        ! Computing derivatives of equation i with respect to variable j
        do j = 1,n
            ! Perturbation on the j-th variable
            dh = h0
            dhval = hval*x0(j)
            dh(j) = dhval
            ! Pertubed vector of variables
            xnow  = x0 +  dh
            ! Trick to reduce finite precisione error
            dh    = xnow - x0
            ! Evaluating with perturbed x
            call FN(fnow,xnow,PARAM)
            do i = 1,n
                ! Numerical derivative: equation i, variable j
                D(i,j) = (fnow(i) - f0(i))/dh(j)
                print "(a,i3,5f14.8)", 'D', i, D(i,j), fnow(i) - f0(i), fnow(i), f0(i), dh(j)
            end do
        end do
    end subroutine fdjac
    
    
    subroutine fdjac_par(f0,x0,FN,PARAM, D)
        ! Computes the numerical derivative of a system of equations
        ! Using OpenMP to parallelize
        ! Usage:
        ! call fdjac_par(x0,FN,PARAM, D)
        ! Inputs:
        ! f0     - initial function value, i.e. FN(x0)
        ! x0     - initial value for x
        ! FN     - name of the subroutine with a function containing an NL system of equations
        ! PARAM  - type(PARPASS_NL) structure with the parameters of FN
        ! Output:
        ! D      - Jacobian matrix with numerical derivatives
        ! Note:
        ! FN must have a dummy argument PARAM.
        ! Requires module -par_pass- defining the derived data type -parpass_nl-
        real(dp),intent(in)           :: f0(:),x0(:)
        type(parpass_nl),intent(inout):: PARAM
        type(parpass_nl)              :: PARAM2
        real(dp),intent(out)          :: D(size(x0),size(x0))
        real(dp)                      :: xnow(size(x0,1))
        real(dp)                      :: fnow(size(x0,1))
        real(dp)                      :: h0(size(x0,1)),hval
        real(dp)                      :: dh(size(x0,1)),dhval
        integer                       :: n,j,i
        interface
            subroutine FN(outf,y,PARAM)
                use lib_kind
                use par_pass
                real(dp),intent(in)           :: y(:)
                type(parpass_nl),intent(inout):: PARAM
                real(dp),intent(out)          :: outf(:)
            end subroutine FN
        end interface
        n    = size(x0,1)        ! Number of variables
        h0   = 0.0_dp            ! Used to generate a basis vector
        hval = 1.0e-6_dp         ! Step used to compute derivatives
        D    = 0.0_dp            ! Where the approximated Jacobian will be stored
        ! Computing derivatives of equation i with respect to variable j
        !$omp parallel private(j,dh,dhval,xnow,fnow,i,PARAM2) 
        PARAM2 = PARAM
        !$omp do
        do j = 1,n
            ! Perturbation on the j-th variable
            dh = h0
            dhval = hval*x0(j)
            dh(j) = dhval
            ! Pertubed vector of variables
            xnow  = x0 +  dh
            ! Trick to reduce finite precisione error
            dh    = xnow - x0
            ! Evaluating with perturbed x
            !PARAM2%test_j = j
            call FN(fnow,xnow,PARAM2)
            ! Numerical derivative: all equations wrt variable j
            do i = 1,n
                D(i,j) = (fnow(i) - f0(i))/dh(j)
            end do
        end do
        !$omp end do
        !$omp end parallel
    end subroutine fdjac_par


    subroutine golden(xl,xr,tol,FN,PARAM, x_sol,f_sol)
        ! Implements the golden section search method to find a
        ! minimum (one-dimensional problems).
        ! Usage:
        ! call golden(xl,xr,tol,FN,PARAM, x_sol,f_sol)
        ! Inputs:
        ! xl,xr - left and right edges of the search interval
        ! tol   - stoping criteria
        ! FN    - name of the function to be solved, must be a (scalar) double precision function
        ! PARAM - type(parpass_g) structure with the parameters of FN
        ! Output:
        ! x_sol - solution to min{FN(X,PARAM)}
        ! f_sol - minimum of FN at x_sol
        ! Note:
        ! FN must have a dummy argument PARAM.
        ! Requires module -par_pass- defining the derived data type -parpass_g-
        real(dp),intent(inout)       :: xl,xr
        real(dp),intent(in)          :: tol
        type(parpass_g),intent(inout):: PARAM
        real(dp),intent(out)         :: x_sol,f_sol
        real(dp)                     :: a,b,c,d,rg,diff,fc,fd
        integer                      :: count
        interface
            real(dp) function FN(y,PARAM)
                use par_pass
                real(dp), intent(in)          :: y
                type(parpass_g), intent(inout):: PARAM
            end function FN
        end interface
        diff = 2.0_dp*tol
        a  = xl
        b  = xr
        rg = (3.0_dp-sqrt(5.0_dp))/2.0_dp
        c  = (1.0_dp-rg)*a + rg*b
        d  = rg*a + (1-rg)*b
        fc = FN(c,PARAM)
        fd = FN(d,PARAM)
        count = 0
        do while(diff > tol)
            if (fc >= fd) then
                a  = c
                c  = d
                d  = rg*a + (1.0_dp-rg)*b
                fc = fd 
                fd = FN(d,PARAM)
            else
                b  = d
                d  = c
                c  = (1.0_dp-rg)*a+rg*b
                fd = fc
                fc = FN(c,PARAM)
            end if
            diff = abs(b-a)
            count = count+1
        end do
        x_sol = (a+b)/2.0_dp
        f_sol = (fc+fd)/2.0_dp
    end subroutine golden

  
    subroutine amoeba(x0,lambda,FN,tol,max_feval,PARAM,savesx, x_opt,fn_opt)
        ! Performs the Nelder-Mead optimization (minimization) search in
        ! multi-dimensional minimization problems
        ! Usage:
        ! call amoeba(x0,lambda,FN,tol,max_feval,PARAM[,savesx], x_opt,fn_opt)
        ! Inputs:
        ! x0     - a (n+1,n) matrix containing the initial guesses for the solution
        ! lambda - initial step size for the simplex
        ! FN     - name of the subroutine with the function to be solved
        ! tol    - stoping criteria
        ! PARAM  - type(parpass_a) structure with the parameters of FN
        ! savesx - Save current vector in text file (1 or 0, optional)
        ! Output:
        ! x_opt  - solution to min {FN(x,PARAM)=0}
        ! fn_opt - value of FN(x_opt,PARAM)
        ! Note:
        ! FN must have a dummy argument PARAM.
        ! Requires module -par_pass- defining the derived data type -parpass_a-
        real(dp),intent(in)          :: x0(:),lambda(:),tol
        integer,intent(in)           :: max_feval
        type(parpass_nl),intent(inout):: PARAM
        integer,intent(in),optional  :: savesx
        real(dp),intent(out)         :: x_opt(size(x0)),fn_opt
        real(dp)                     :: x(size(x0)+1,size(x0)),f(size(x0)+1)
        real(dp)                     :: rho,xi,gam,sig
        real(dp)                     :: iden(size(x0),size(x0)),diff,f_r,f_e,f_c
        real(dp)                     :: x_bar(size(x0)),x_r(size(x0)),x_e(size(x0)),x_c(size(x0))
        integer                      :: i,n_dim,n_feval,i_low,i_high
        interface
            real(dp) function FN(y,PARAM)
                use par_pass
                real(dp),intent(in)          :: y(:)
                type(parpass_nl),intent(inout):: PARAM
            end function FN
        end interface
        ! Initial simplex
        n_dim  = size(x0)
        x      = 0.0_dp
        x(1,:) = x0
        iden   = 0.0_dp
        do i = 2,n_dim+1
            iden(i-1,i-1) = 1.0_dp*lambda(i-1)
            x(i,:) = x0 + iden(i-1,:)
        end do
        ! Define algorithm constants
        rho = 1.0_dp  ! rho > 0
        xi  = 2.0_dp  ! xi > max(rho, 1)
        gam = 0.5_dp  ! 0 < gam < 1
        sig = 0.5_dp  ! 0 < sig < 1
        ! Initialization
        f = 0.0_dp
        do i = 1,n_dim+1
            f(i) = FN(x(i,:),PARAM)
        end do
        n_feval = n_dim+1
        ! Sorting
        i_low  = minloc(f,1)
        i_high = maxloc(f,1)
        diff = 2.0_dp*tol
        do while ((diff > tol) .and. (n_feval <= max_feval))
            ! Compute the midpoint of the simplex opposite the worst point
            do i = 1,n_dim
                x_bar(i) = (sum(x(1:n_dim+1,i))-x(i_high,i))/n_dim
            end do
            ! Reflection though the centroid
            x_r = x_bar + rho*(x_bar-x(i_high,:))
            f_r = FN(x_r,PARAM)
            n_feval = n_feval + 1
            ! Accepting the point?
            if ((f(i_low) <= f_r) .and. (f_r <= f(n_dim))) then
                ! reflection
                x(i_high,:) = x_r
                f(i_high)   = f_r
            elseif (f_r  <  f(i_low)) then
                ! Test for possible expansion
                x_e = x_bar + rho*xi*(x_bar-x(i_high,:))
                f_e = FN(x_e,PARAM)
                n_feval = n_feval+1
                ! Can we accept the expanded point?
                if (f_e < f_r) then
                    ! expansion
                    x(i_high,:) = x_e
                    f(i_high)   = f_e
                else
                    ! eventual reflection
                    x(i_high,:) = x_r
                    f(i_high)   = f_r
                end if
            else if ((f(n_dim) <= f_r) .and. (f_r < f(i_high))) then
                ! Outside contraction
                x_c = x_bar + rho*gam*(x_bar-x(i_high,:))
                f_c = FN(x_c,PARAM)
                n_feval = n_feval + 1
                ! Accepting the contracted point
                if (f_c <= f_r) then
                    ! Outside contraction
                    x(i_high,:) = x_c
                    f(i_high)   = f_c
                else
                    ! Shrink
                    do i = 1,n_dim+1
                        if (i /= i_low) then
                            x(i,:) = x(i_low,:) + sig*(x(i,:)-x(i_low,:))
                            f(i)   = FN(x(i,:),PARAM)
                        end if
                    end do
                    n_feval = n_feval + n_dim
                end if
            else
                ! Since f(xr)>=f(xn+1), try inside contraction
                x_c = x_bar - gam*(x_bar-x(i_high,:))
                f_c = FN(x_c,PARAM)
                n_feval = n_feval + 1
                !  Can we accept the contracted point?
                if (f_c < f(i_high)) then
                    ! Inside contraction
                    x(i_high,:) = x_c
                    f(i_high)   = f_c
                else
                    ! Shrink
                    do i = 1,n_dim+1
                        if (i /= i_low) then
                            x(i,:) = x(i_low,:) + sig*(x(i,:)-x(i_low,:))
                            f(i)   = FN(x(i,:),PARAM)
                        end if
                    end do
                    n_feval = n_feval + n_dim
                end if
            end if
            ! Re-sort the points
            i_low  = minloc(f,1)
            i_high = maxloc(f,1)
            diff   = abs(f(i_high)-f(i_low))
            if (present(savesx)) then
                if (savesx == 1) call txt_write(x,'Simplex_Amoeba')
            end if
        end do
        ! Solution
        x_opt  = x(1,:)
        fn_opt = f(1)
        if (diff > tol) then
            print *, 'Warning (amoeba): Nelder-Mead Warning!'
            print *, 'Warning (amoeba): The maximum number of function evaluations was exceeded'
        end if
    end subroutine amoeba
  
  
end module lib_solver