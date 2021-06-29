module lib_basic

! Module to implement some (Basic) numerical routines.
! By: Christian Bustamante
! Email: cbustaam@gmail.com
! Last modified: 16 Oct, 2020, 17:32
! v 2.60
!
! To compile this module:
! $ ifort -c Lib_Basic.f90 Lib_Kind.f90
! $ gfortran -c Lib_Basic.f90 Lib_Kind.f90
! To use it put this just after the -program- statement:
! use lib_basic
!
! Subroutines:
! (Basic) 01. linspace      - Creates a linearly spaced grid
! (Basic) 02. logspace      - Creates a log-linearly spaced grid starting from a negative number
! (Basic) 03. invlogspace   - Creates a log-linearly spaced grid that for whih its density grows backwards 
! (Basic) 04. doublelogspace- Creates a log-linearly spaced grid with more density near an specified interior point   
! (Basic) 05. inv           - Returns the inverse of a matrix
! (Basic) 06. kron          - Computes the kroenecker product of two matrices
! (Basic) 07. findmaxeig    - Computes the dominant eigenvector and eigenvalue using the power method
! (Basic) 08. gridlookup    - Find nearest (at or below) gridpoint
! (Basic) 09. wgridlookup   - Find nearest (at or below) gridpoint of X to x0 and weight for linear interpolation of that point
! (Basic) 10. rescale       - Rescale any variable/parameter to or from any desired scale.
! (Basic) 11. tick          - Returns the initial reference time (used later by tock_sec or tock_hms)
! (Basic) 12. tick_wc       - Returns the initial reference time from wall clock (used later by tock_sec_wp or tock_hms_wp). Less precise than tick
! (Basic) 13. tock_sec      - Returns the time elapsed in seconds
! (Basic) 14. tock_sec_wc   - Returns the time elapsed in seconds from wall clock. Less precise than tock_sec
! (Basic) 15. tock_hms      - Returns the time elapsed in hours : minutes : seconds
! (Basic) 16. tock_hms_wc   - Returns the time elapsed in hours : minutes : seconds from wall clock. Less precise than tock_sec
! (Basic) 17. cumsum        - Computes the cummulative sum of a vector
! (Basic) 18. quick_sort    - Sorts a vector form low to high
! (Basic) 19. sort_rows     - Sorts ascendently matrix rows by the specified column
! (Basic) 20. qsort_rows    - Sorts ascendently matrix rows by the specified column (quicker)
! (Basic) 21. date_str      - Returns current date in string format as YYYYMMDD
! (Basic) 22. time_str      - Returns current clock time in string format as HHMMSS
! (Basic) 23. vec           - Vectorizes a matrix
! (Basic) 24. unvec         - Un-vectorizes a vector and returns a matrix
! (Basic) 25. arr2mat       - Reshapes an 3-dimensional array into a matrix by stacking each submatrix
! (Basic) 26. mat2arr       - Un-vectorizes a vector and returns a matrix
! (Basic) 27. eye           - Generates an identity matrix of order n
! (Basic) 28. progress      - Prints progress bar, plus the percentage progress
!
! Note: depends on Lib_Kind

use lib_kind
implicit none
private
public :: linspace, logspace, invlogspace, doublelogspace, inv, kron, findmaxeig, gridlookup, wgridlookup, &
          rescale, tick, tick_wc, tock_sec, tock_sec_wc, tock_hms, tock_hms_wc, cumsum, quick_sort,        &
          sort_rows, qsort_rows, date_str, time_str, vec, unvec, arr2mat, mat2arr, eye, progress,          &
          dlogspace_eq

contains


    subroutine linspace(x1,xN,N, X)
        ! Creates a linearly spaced grid for X.
        ! Usage: 
        ! call linspace(x1,xN,N, X)
        ! Inputs:
        ! x1    - first point in the grid
        ! xN    - last point in the grid
        ! N     - number of points
        ! Output:
        ! X     - equally space grid
        integer,intent(in)   :: N
        real(dp),intent(in)  :: x1,xN
        real(dp),intent(out) :: X(:)
        integer              :: i
        X(1) = x1
        X(N) = xN
        do i = 2,(N-1)
            X(i) = x1 + real((i-1),dp)*real((xN-x1)/(N-1.0_dp),dp)
        end do
    end subroutine linspace
  
  
    subroutine logspace(x1,xN,N, X)
        ! Creates a log-linearly spaced vector for X.
        ! It can star from a negaritive value, i.e. x1<0
        ! Usage:
        ! call logspace(x1,xN,N, X)
        ! Inputs:
        ! x1    - first point in the grid
        ! xN    - last point in the grid
        ! N     - number of points
        ! Output:
        ! X     - log-linearly spaced grid
        integer,intent(in)   :: N
        real(dp),intent(in)  :: x1,xN
        real(dp),intent(out) :: X(:)
        real(dp)             :: X0(N)
        call linspace(log(x1 - 1.0_dp*x1 + 1.0_dp)/log(10.0_dp),log(xN - 1.0_dp*x1 + 1.0_dp)/log(10.0_dp),N,X0)
        X0 = 10.0_dp**X0
        X  = X0 + (x1 - 1.0_dp)
    end subroutine logspace
  
  
    subroutine invlogspace(x1,xN,N, X)
        ! Creates an inverse-logspaced grid that goes backwards, i.e., 
        ! more points near xN  rather than near x1.
        ! Usage: 
        ! call invlogspace(x1,xN,N, X)
        ! Inputs:
        ! x1    - first point in the grid
        ! xN    - last point in the grid
        ! N     - number of points
        ! Output:
        ! X     - inversely log-linearly spaced grid
        integer,intent(in)   :: N
        real(dp),intent(in)  :: x1,xN
        real(dp),intent(out) :: X(:)
        integer              :: i
        real(dp)             :: X0(N),diffx0(N-1)
        ! Generating a log-spaced grid and obtaining differences between gridpoints
        call logspace(x1,xN,N, X0)
        diffx0 = X0(2:N) - X0(1:N-1)
        ! Adding those differences, but in a backwards manner
        X(1) = x1
        do i = 2,N
            X(i) = X(i-1) + diffx0(N-i+1)
        end do
    end subroutine invlogspace
  
  
    subroutine doublelogspace(x1,xmid,xN,N, X)
        ! Generates a double spaced grid, this is a grid that is proportionally
        ! inversely logspaced in [xlow,xmid] and proportionaly logspaced in [xmid,xhigh]. 
        ! This means that there is a major concentration of points near xmid.
        ! It is required that x1 < xmid < xN.
        ! Usage: 
        ! call doublelogspace(x1,xmid,xN,N, X)
        ! Inputs:
        ! x1    - first point in the grid
        ! xmid  - inflection point in the grid
        ! xN    - last point in the grid
        ! N     - number of points
        ! Output:
        ! X     - double log-linearly spaced grid
        integer,intent(in)   :: N
        real(dp),intent(in)  :: x1,xmid,xN
        real(dp),intent(out) :: X(:)
        integer              :: xnumL,xnumH
        real(dp)             :: propH
        real(dp),allocatable :: XL(:),XH(:)
        ! Computing proportions of points in [x1,xmid] and [xmid,xN]
        propH = (xN-xmid)/(xN-x1)
        xnumH = int(propH*real(N,dp))
        xnumL = N - xnumH + 1
        allocate(XL(xnumL),XH(xnumH))
        ! Generating a inversely log-spaced grid in [x1,xmin]
        call invlogspace(x1,xmid,xnumL, XL)
        ! Generating a log-spaced grid in [xmin,xN]
        call logspace(xmid,xN,xnumH, XH)
        ! Combining the two grids (note that there xmid is in both grid, 
        ! so we just take one of those)
        X = [XL(1:xnumL-1), XH]
    end subroutine doublelogspace
    
    
    subroutine dlogspace_eq(x1,xmid,xN,N, X)
        ! Generates a double spaced grid, this is a grid that is equally
        ! inversely logspaced in [xlow,xmid] and equally logspaced in [xmid,xhigh]. 
        ! This means that there is the same concentration of points near xmid.
        ! It is required that x1 < xmid < xN.
        ! Usage: 
        ! call dlogspace_eq(x1,xmid,xN,N, X)
        ! Inputs:
        ! x1    - first point in the grid
        ! xmid  - inflection point in the grid
        ! xN    - last point in the grid
        ! N     - number of points
        ! Output:
        ! X     - double log-linearly spaced grid
        integer,intent(in)   :: N
        real(dp),intent(in)  :: x1,xmid,xN
        real(dp),intent(out) :: X(:)
        integer              :: xnumL,xnumH
        real(dp),allocatable :: XL(:),XH(:)
        ! Computing number of points in [x1,xmid] and [xmid,xN]
        xnumH = int(0.5_dp*real(N,dp))
        xnumL = N - xnumH + 1
        allocate(XL(xnumL),XH(xnumH))
        ! Generating a inversely log-spaced grid in [x1,xmin]
        call invlogspace(x1,xmid,xnumL, XL)
        ! Generating a log-spaced grid in [xmin,xN]
        call logspace(xmid,xN,xnumH, XH)
        ! Combining the two grids (note that there xmid is in both grid, 
        ! so we just take one of those)
        X = [XL(1:xnumL-1), XH]
    end subroutine dlogspace_eq
  
  
    function inv(a) result(ainv)
        ! Inverse matrix (Based on Doolittle LU factorization for Ax=b)
        ! by Alex G. December 2009 (with slight modifications)
        ! Source: https://ww2.odu.edu/~agodunov/computing/programs/book2/Ch06/Inverse.f90
        ! Usage: 
        ! invA  = inv(A)
        ! Inputs:
        ! A     - Matrix to be inverted
        ! Output:
        ! invA  - Inverse of matrix A
        real(dp),dimension(:,:),intent(in)      :: a(:,:)
        real(dp),dimension(size(a,1),size(a,2)) :: ainv
        real(dp),dimension(size(a,1),size(a,2)) :: acopy,L,U
        real(dp),dimension(size(a,1))           :: b,d,x
        real(dp)                                :: coeff
        integer                                 :: n,i, j, k
        acopy = a
        n = size(a,1)
        ! Step 0: Initialization for matrices L and U and b
        L = 0.0_dp
        U = 0.0_dp
        b = 0.0_dp
        ! Step 1: Forward elimination
        do k = 1,n-1
            do i = k+1,n
                coeff  = acopy(i,k)/acopy(k,k)
                L(i,k) = coeff
                do j = k+1,n
                    acopy(i,j) = acopy(i,j)-coeff*acopy(k,j)
                end do
            end do
        end do
        ! Step 2: Prepare L and U matrices
        ! L matrix is a matrix of the elimination coefficient + the diagonal elements are 1.0
        do i = 1,n
            L(i,i) = 1.0_dp
        end do
        ! U matrix is the upper triangular part of A
        do j = 1,n
            do i = 1,j
                U(i,j) = acopy(i,j)
            end do
        end do
        ! Step 3: Compute columns of the inverse matrix C
        do k = 1,n
            b(k) = 1.0_dp
            d(1) = b(1)
            ! Step 3a: Solve Ld=b using the forward substitution
            do i = 2,n
                d(i) = b(i)
                do j = 1,i-1
                    d(i) = d(i) - L(i,j)*d(j)
                end do
            end do
            ! Step 3b: Solve Ux=d using the back substitution
            x(n) = d(n)/U(n,n)
            do i = n-1,1,-1
                x(i) = d(i)
                do j = n,i+1,-1
                    x(i) = x(i)-U(i,j)*x(j)
                end do
                x(i) = x(i)/u(i,i)
            end do
            ! Step 3c: Fill the solutions x(n) into column k of C
            do i = 1,n
                ainv(i,k) = x(i)
            end do
            b(k) = 0.0_dp
        end do
    end function inv
  
  
    subroutine kron(A,B, K)
        ! Computes the kroenecker product of two matrices
        ! Usage: 
        ! call kron(A,B, K)
        ! Inputs:
        ! A     - matrix of dimension MA x NA
        ! B     - matrix of dimension MB x NB
        ! Output:
        ! K    - kroenecker product of A and B
        real(dp),intent(in)  :: A(:,:),B(:,:)
        real(dp),intent(out) :: K(:,:)
        integer              :: i,j,MA,NA,MB,NB
        MA = ubound(A,1)
        NA = ubound(A,2)
        MB = ubound(B,1)
        NB = ubound(B,2)
        if (size(K,1) /= MA*MB .or. size(K,2) /= NA*NB) then
            print *, 'Error (kron): K has invalid size'
            call abort
        end if
        forall(i=1:MA, j=1:NA)
            K(MB*(i-1)+1:MB*I,NB*(j-1)+1:NB*j) = A(i,j)*B
        end forall
    end subroutine kron
  
  
    subroutine findmaxeig(Matrix,steps, Eigvec,Eigval)
        ! Computes the dominant eigenvector and eigen value of a matrix using the
        ! power method (with scaling) for approximating eigenvalues.
        ! Usage:
        ! call findmaxeig(Matrix,steps, Eigvec,Eigval)
        ! Inputs:
        ! Matrix - dominant eigenvector and eigenvalue will be computed for it
        ! steps  - number of iterations
        ! Output:
        ! Eigvec - dominant eigenvector
        ! Eigval - dominant eigenvalue
        real(dp),intent(in)  :: Matrix(:,:)
        integer,intent(in)   :: steps
        real(dp),intent(out) :: Eigvec(size(Matrix,1)),Eigval
        real(dp)             :: Eigvec1(size(Matrix,1))
        real(dp)             :: tol,diff
        integer              :: n,i
        tol = 1e-8_dp
        n   = size(Matrix,1)
        Eigvec = 1.0_dp                       ! Initialize eigen vector to any value.
        do i = 1,steps
            Eigvec1 = matmul(Matrix,Eigvec)   ! Multiply input matrix by eigenvector
            Eigval  = maxval(abs(Eigvec1))    ! Find eigenvalue
            if (abs(Eigval) < 1.0D-7) then
                exit
            end if
            Eigvec1 = Eigvec1/Eigval
            diff = maxval(abs(Eigvec1-Eigvec))
            if (diff < tol) then
                exit
            end if
            Eigvec = Eigvec1
        end do
    end subroutine findmaxeig


    integer function gridlookup(X,x0) result(index)
        ! Find nearest (at or below) gridpoint of X to x0
        ! Usage:
        ! index = gridlookup(X,x0)
        ! Inputs:
        ! X     - one-dimensional grid
        ! x0    - point
        ! Output:
        ! index - index of the nearest (at or below) gridpoint to x0
        real(dp),intent(in) :: X(:),x0
        integer             :: nX,ilow,ihigh,inow,diff
        real(dp)            :: xnow
        nX    = size(X)
        ilow  = 1
        ihigh = nX
        diff  = 2
        do while (diff > 1)
            inow = (ilow+ihigh)/2
            xnow = X(inow)
            if (xnow > x0) then
                ihigh = inow
            else
                ilow  = inow
            end if
            diff = ihigh - ilow
        end do
        index = ilow
    end function gridlookup
    
    
    subroutine wgridlookup(nX,X,x0, index,weight)
        ! Find nearest (at or below) gridpoint of X to x0 and weight for 
        ! linear interpolation of that point
        ! Usage:
        ! call wgridlookup(nX,X,x0, index,weight)
        ! Inputs:
        ! nX     - grid length
        ! X      - one-dimensional grid
        ! x0     - point
        ! Output:
        ! index  - index of the nearest (at or below) gridpoint to x0
        ! weight - weight for index when linearly interpolating
        integer,intent(in)  :: nX
        real(dp),intent(in) :: X(:),x0
        integer,intent(out) :: index
        real(dp),intent(out):: weight
        integer             :: ilow,ihigh,inow,diff
        real(dp)            :: xnow
        ilow  = 1
        ihigh = nX
        diff  = 2
        do while (diff > 1)
            inow = (ilow+ihigh)/2
            xnow = X(inow)
            if (xnow > x0) then
                ihigh = inow
            else
                ilow  = inow
            end if
            diff = ihigh - ilow
        end do
        index = ilow
        weight = X(index+1) - x0
        weight = weight/(X(index+1) -  X(index))
        weight = max(weight, 0.0_dp)
        weight = min(weight, 1.0_dp)
    end subroutine wgridlookup
  
  
    function rescale(x,lb,ub,caser) result(y)
        ! Rescale any variable/parameter to or from any desired scale:
        ! non-negative real numbers, or a closed interval.
        ! Usage:
        ! y = rescale(x,lb,ub,caser)
        ! Inputs:
        ! x      - any real number to be re-scaled
        ! lb     - lower bound
        ! ub     - lower bound (larger than 1e8 if infinity wanted)
        ! caser  - 1 if rescaling to [lb,up], -1 if rescaling from.
        ! Output:
        ! y      - rescaled variable to/from the interval [lb,ub]
        real(dp),intent(in) :: x,lb,ub
        integer,intent(in)  :: caser
        real(dp)            :: y
        real(dp),parameter  :: hugen = 1e6_dp
        if (lb >= ub) then
            print *, 'Error (rescale): Lower bound must me less than upper bound'
        end if
        if ((caser /= 1) .and. (caser /= -1)) then
            print *, 'Error (rescale): Case must be 1 or -1'
        end if
        if (ub > hugen) then
            ! Scale: Non-negative numbers
            if (caser == 1) then
                ! Rescaling to (0,+Inf)
                y = exp(x)
            else if (caser == -1) then
                ! Rescaling back from (0,+Inf)
                y = log(x)
            end if
        else if (ub <= hugen) then
            ! Scale: Intervale [lb,ub]
            if (caser == 1) then
                ! Rescaling to [lb,ub]
                y = lb + (ub-lb)*exp(x)/(1.0_dp+exp(x))
            else if (caser == -1) then
                ! Rescaling back from [lb,ub]
                y = log((x-lb)/(ub-x))
            end if
        else
            print *, 'Error (rescale): Please learn to use this function properly'
        end if
    end function rescale
  
  
    subroutine tick(t_begin)
        ! Returns the inital reference time (used later by tock)
        ! Usage: 
        ! call tick(t_begin)
        ! Output:
        ! t_begin - inital reference time (real_dp)
        real(dp),intent(out) :: t_begin
        call cpu_time(t_begin)
    end subroutine tick
  
  
    subroutine tick_wc(t_begin)
        ! Returns the inital reference time (used later by tock_wc)
        ! Time is from wall clock. Needed when using multi-threading
        ! Usage: 
        ! call tick_wc(t_begin)
        ! Output:
        ! t_begin - inital reference time (integer)
        integer,intent(out) :: t_begin
        call system_clock(t_begin)
    end subroutine tick_wc
  
  
    subroutine tock_sec(t_begin, secs)
        ! Returns the time elapsed in seconds
        ! Usage: 
        ! call tock_sec(t_begin, secs)
        ! Inputs:
        ! t_begin - initial reference time (from tick, real_dp)
        ! Outputs:
        ! secs    - seconds elapsed
        real(dp),intent(in) :: t_begin
        real(dp),intent(out):: secs
        real(dp)            :: t_end
        call cpu_time(t_end)
        t_end = t_end - t_begin
        secs  = t_end
    end subroutine tock_sec
  
  
    subroutine tock_sec_wc(t_begin, secs)
        ! Returns the time elapsed in seconds
        ! Time is from wall clock. Needed when using multi-threading
        ! Usage: 
        ! call tock_sec_wc(t_begin, secs)
        ! Inputs:
        ! t_begin - initial reference time (from tick_wc, integer)
        ! Outputs:
        ! secs    - seconds elapsed
        integer, intent(in) :: t_begin
        real(dp),intent(out):: secs
        integer             :: now, clock_rate, count_max
        call system_clock(now,clock_rate,count_max)
        secs = real(now - t_begin,dp)/real(clock_rate,dp)
    end subroutine tock_sec_wc
  
  
    subroutine tock_hms(t_begin, hours,mins,secs)
        ! Returns the time elapsed in the form of three scalars
        ! hours : minutes : seconds
        ! Usage: 
        ! call tock_hms(t_begin, hours,mins,secs)
        ! Inputs:
        ! t_begin - initial reference time (from tick, real_dp)
        ! Outputs:
        ! hours   - hours elapsed
        ! minutes - minutes elapsed
        ! seconds - seconds elapsed
        real(dp),intent(in) :: t_begin
        integer,intent(out) :: hours,mins,secs
        real(dp)            :: t_end
        call cpu_time(t_end)
        t_end = t_end - t_begin
        hours = floor((t_end/(60*60)))
        mins  = floor((t_end-60*60*hours)/60)
        secs  = floor((t_end-60*60*hours-60*mins))    
    end subroutine tock_hms
  
  
    subroutine tock_hms_wc(t_begin, hours,mins,secs)
        ! Returns the time elapsed in the form of three scalars
        ! hours : minutes : seconds
        ! Time is from wall clock. Needed when using multi-threading
        ! Usage: 
        ! call tock_hms_wc(t_begin, hours,mins,secs)
        ! Inputs:
        ! t_begin - initial reference time (from tick_wc, integer)
        ! Outputs:
        ! hours   - hours elapsed
        ! minutes - minutes elapsed
        ! seconds - seconds elapsed
        integer, intent(in) :: t_begin
        integer,intent(out) :: hours,mins,secs
        integer             :: now, clock_rate
        real(dp)            :: t_end
        call system_clock(now,clock_rate)
        t_end = real(now - t_begin,dp)/real(clock_rate,dp)
        hours = floor((t_end/(60*60)))
        mins  = floor((t_end-60*60*hours)/60)
        secs  = floor((t_end-60*60*hours-60*mins))    
    end subroutine tock_hms_wc
  
  
    subroutine cumsum(x, cumx)
        ! Computes the cummulative sum of a vector
        ! Usage: 
        ! call cumsum(x, cumx)
        ! Inputs:
        ! x      - input vector
        ! Outputs:
        ! cumx   - cummulative sum of x
        real(dp),intent(in) :: x(:)
        real(dp),intent(out):: cumx(size(x))
        integer             :: n,i
        n = size(x)
        cumx(1) = x(1)
        do i = 2,n
            cumx(i) = cumx(i-1) + x(i)
        end do  
    end subroutine cumsum
  
  
    recursive subroutine quick_sort(a,first,last)
        ! Sorts vector a
        ! By t-nissie
        ! Source: https://gist.github.com/t-nissie/479f0f16966925fa29ea
        ! Usage:
        ! call quick_sort(a,first,last)
        ! Inputs:
        ! a     - vector to be sorted
        ! first - integer with the lowest position in a
        ! last  - integer with the highest position in a
        ! Note: the vector a is, implicititly, inout.
        real(dp) :: a(:),x,t
        integer  :: first,last
        integer  :: i,j
        x = a((first+last)/2)
        i = first
        j = last
        do
            do while (a(i) < x)
                i = i+1
            end do
            do while (x < a(j))
                j = j-1
            end do
            if (i >= j) exit
            t    = a(i)
            a(i) = a(j)
            a(j) = t
            i = i+1
            j = j-1
        end do
        if (first < i-1) call quick_sort(a,first,i-1)
        if (j+1 < last)  call quick_sort(a,j+1,  last)
    end subroutine quick_sort
  
  
    subroutine sort_rows(Araw,col2ord, Asort)
        ! Sorts ascendently matrix rows by the specified column
        ! Usage:
        ! call sortrows(Araw,col2ord, Asort)
        ! Inputs:
        ! Araw    - matrix to be row-sorted
        ! col2ord - column that will be used to order matrix A
        ! Output:
        ! Asort   - sorted matrix
        real(dp),intent(in)  :: Araw(:,:)
        integer, intent(in)  :: col2ord
        real(dp),intent(out) :: Asort(size(Araw,1),size(Araw,2))
        integer              :: nrow,irow,krow
        real(dp)             :: buff(size(Araw,2))
        nrow  = size(Araw,1)
        Asort = Araw
        do irow = 1,nrow
            krow = minloc(Asort(irow:nrow, col2ord),DIM=1) + irow - 1
            buff(:)       = Asort(irow,:)
            Asort(irow,:) = Asort(krow,:)
            Asort(krow,:) = buff(:)
        end do
    end subroutine sort_rows
    
    
    recursive subroutine qsort_rows(Amat,first,last,col2ord)
        ! Sorts ascendenly a matrix by the specified column
        ! Usage:
        ! call qsort_rows(Amat,first,last,col2ord)
        ! Inputs:
        ! Amat    - matrix to be sorted
        ! first   - integer with the lowest position in a
        ! last    - integer with the highest position in a
        ! col2ord - column to order ascendenly
        ! Note: the vector a is, implicititly, inout.
        integer  :: col2ord
        real(dp) :: Amat(:,:),x,Row(size(Amat,2))
        integer  :: first,last
        integer  :: i,j
        x = Amat((first+last)/2, col2ord)
        i = first
        j = last
        do
            do while (Amat(i, col2ord) < x)
                i = i+1
            end do
            do while (x < Amat(j, col2ord))
                j = j-1
            end do
            if (i >= j) exit
            Row = Amat(i, :)
            Amat(i,:) = Amat(j,:)
            Amat(j,:) = Row
            i = i+1
            j = j-1
        end do
        if (first < i-1) call qsort_rows(Amat,first,i-1, col2ord)
        if (j+1 < last)  call qsort_rows(Amat,j+1,  last,col2ord)
    end subroutine qsort_rows
  
  
    subroutine date_str(dstr)
        ! Returns current date in string format as YYYYMMDD
        ! Usage:
        ! call date_str(dstr)
        ! Output:
        ! dstr  - string with date
        character(8),intent(out) :: dstr
        integer                  :: vdate(8)
        character(4)             :: yyyy
        character(2)             :: dd,mm
        call date_and_time(values=vdate)
        write (yyyy,'(i4)') vdate(1)
        write (mm,'(i2.2)') vdate(2)
        write (dd,'(i2.2)') vdate(3)
        dstr = yyyy//mm//dd
    end subroutine date_str
  
  
    subroutine time_str(tstr)
        ! Returns current clock time in string format as HHMMSS
        ! Usage:
        ! call time_str(tstr)
        ! Output:
        ! tstr  - string with clock time
        character(6),intent(out) :: tstr
        integer                  :: timec(3)
        character(2)             :: hh,mm,ss
        call itime(timec)
        write (hh,'(i2.2)') timec(1)
        write (mm,'(i2.2)') timec(2)
        write (ss,'(i2.2)') timec(3)
        tstr = hh//mm//ss
    end subroutine time_str
  
  
    function vec(X) result(xv)
        ! Vectorizes a matrix (2-dimensional array)
        ! Usage:
        ! xv = vec(X)
        ! Input:
        ! X      - 2-dimensional array
        ! Output:
        ! xv     - vectorized matrix
        real(dp),intent(in) :: X(:,:)
        real(dp)            :: xv(size(X,1)*size(X,2))
        integer             :: n1,n2,i,pos0
        n1 = size(X,1)
        n2 = size(X,2)
        do i = 1,n2
            pos0 = (i-1)*n1
            xv(pos0+1:pos0+n1) = X(:,i)
        end do
    end function vec
  
  
    function unvec(xv,n1,n2) result(X)
        ! Unvectorizes a vector, reshaping it as a 2-dimensional array
        ! Usage:
        ! X = vec(xv,n1,n2)
        ! Inputs:
        ! xv     - vectorized matrix
        ! n1     - length of row dimension
        ! n2     - length of column dimension
        ! Outputs:
        ! X      - 2-dimensional array
        real(dp),intent(in) :: xv(:)
        integer, intent(in) :: n1,n2
        real(dp)            :: X(n1,n2)
        integer             :: i,pos0
        do i = 1,n2
            pos0 = (i-1)*n1
            X(:,i) = xv(pos0+1:pos0+n1)
        end do
    end function unvec
    
    
    function arr2mat(x3) result(xM)
        ! Reshapes an 3-dimensional array into a matrix by stacking each 
        ! matrix (:,:,i) on top of each other with i=1,2,...,N3
        ! Usage:
        ! xM = vec(x3)
        ! Input:
        ! x3     - 3-dimensional array
        ! Output:
        ! xM     - stacked matrix
        real(dp),intent(in) :: x3(:,:,:)
        real(dp)            :: xM(size(x3,3)*size(x3,1),size(x3,2))
        integer             :: n1,n2,n3,i,idx1,idx2
        n1 = size(x3,1)
        n2 = size(x3,2)
        n3 = size(x3,3)
        do i = 1,n3
            idx1 = (i-1)*n1+1
            idx2 = i*n1
            xM(idx1:idx2, :) = x3(:,:,i)            
        end do    
    end function arr2mat
    
    
    function mat2arr(xM,n1,n2,n3) result(x3)
        ! Reshapes an matrix into a 3-dimensional array. 
        ! Undoes what -arr2mat- does
        ! Usage:
        ! x3 = vec(xM,n1,n2,n3)
        ! Input:
        ! xM     - stacked matrix
        ! Output:
        ! x3     - 3-dimensional array
        real(dp),intent(in) :: xM(:,:)
        integer,intent(in)  :: n1,n2,n3
        real(dp)            :: x3(n1,n2,n3)
        integer             :: i,idx1,idx2
        if (n3*n1 /= size(xM,1)) print "(a)", '[MAT2ARR] Error: your dimensions are wrong!'
        do i = 1,n3
            idx1 = (i-1)*n1+1
            idx2 = i*n1
            x3(:,:,i) = xM(idx1:idx2, :)
        end do    
    end function mat2arr
  
  
    function eye(n) result(Ident)
        ! Generates an identity matrix of order n
        ! Usage:
        ! X = eye(n)
        ! Inputs:
        ! n      - order of the identity matrix
        ! Output:
        ! Ident  - identity matrix of order n
        integer, intent(in) :: n
        real(dp)            :: Ident(n,n)
        integer             :: i  
        forall(i = 1:n) Ident(i,i) = 1.0
    end function eye
  
  
    subroutine progress(j,jmax,nsegi)
        ! Prints progress bar, plus the percentage progress
        ! Progress can only be measured when there is a deterministic total
        ! number of replications. E.g., when simulating period 23 out of 200
        ! Usage:
        ! call progress(j,jmax,nsegi)
        ! Inputs:
        ! j      - current iteration or period
        ! jmax   - total number of interations or periods
        ! nsegui - number of segments for the bar (optional, default is 20)
        integer,intent(in)          :: j,jmax
        integer,intent(in),optional :: nsegi
        integer,parameter           :: nobar = 8
        integer                     :: nseg,sseg,k,prog,nstar,rstar
        real                        :: pcent
        character(len=8)            :: myfmt
        character(len=:),allocatable:: bar
        ! Number of segments
        nseg = 20
        if (present(nsegi)) nseg = nsegi
        if (mod(100,nseg) /= 0) nseg = 20
        ! Length of segments
        sseg = 100/nseg
        ! Allocating progress bar
        allocate(character(len=nobar+nseg) :: bar)
        ! Writing empty bar
        write(unit=bar(1:nobar-1),fmt="(a)")  " ???% |"
        do k = nobar,nobar+nseg-1
            bar(k:k) = " "
        end do
        bar(nobar+nseg:nobar+nseg) = "|"
        ! Calculate progress (real and integer variables)
        pcent = real(j,dp)/real(jmax,dp)
        prog  = int(100*pcent)
        ! Write actual progress (number)
        write(unit=bar(2:4),fmt="(i3)") prog
        ! Write actual progress (stars)
        nstar = prog/sseg
        rstar = mod(prog,sseg)
        if (rstar >= 3) nstar = nstar+1
        do k = 1,nstar
            bar(nobar-1+k:nobar-1+k) = "*"
        end do
        ! Print the progress bar
        write(myfmt,"(a,i2,a)") "(a1,a",nobar+nseg,")"
        write(unit=6,fmt=myfmt,advance="no") char(13), bar
        if (prog /= 100) then
            flush(unit=6)
        else
            write(unit=6,fmt=*)
        end if
        return
    end subroutine progress
  
  
end module lib_basic