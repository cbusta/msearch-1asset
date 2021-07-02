module lib_interp

! Module to implement some numerical interpolation routines.
! By: Christian Bustamante
! Email: mail@cbustamante.co
! Last modified: 12 Feb 2021, 19:23
! v 1.61
!
! To compile this module:
! $ ifort -c Lib_Interp.f90 Lib_Basic.f90 Lib_Kind.f90 Lib_Smalapack.f90 -lblas -llapack
! $ gfortran -c Lib_Interp.f90 Lib_Basic.f90 Lib_Kind.f90 Lib_Smalapack.f90 -lblas -llapack
! To use it put this just after the -program- statement:
! use lib_interp
!
! Subroutines:
! (Interp) 01. linterp       - Linear interpolation (as subroutine)
! (Interp) 02. linterpx      - Linear interpolation (as function)
! (Interp) 03. bilinterp     - Bilinear interpolation
! (Interp) 04. bilinterp_vec - Bilinear interpolation when function is vectorized
! (Interp) 05. trilinterp    - Trilinear interpolation
! (Interp) 06. tetralinterp  - Tetralinear interpolation   
! (Interp) 07. cheb_grid     - Returns the grid of "zeros" for the Chebyshev polynomials
! (Interp) 08. cheb_coeff    - Returns the coefficients fo Chebyshev polynomial interpolation
! (Interp) 09. eval_cheb     - Evaluates Chebyshev polynomial to compute the approximate F(x)
! (Interp) 10. spline_nak    - Returns the coefficients of a piecewise polynomial cubic spline using the not-a-knot endpoint condition
! (Interp) 11. spline_naka   - First part of -spline_nak-. Returns matrix with coefficients for slope equations needed by -spline_nakb-
! (Interp) 12. spline_nakb   - Second part of -spline_nak-. Returns the coefficients of a piecewise polynomial cubic spline using the outcome of -spline_naka-
! (Interp) 13. eval_spline   - Evaluates the spline, this is, finds the right polynomial coefficients and use them to compute the approximate F(x)
! (Interp) 14. spline2d_nak  - Returns the coefficients of a piecewise polynomial cubic bivariate spline using the not-a-knot endpoint condition
! (Interp) 15. eval_spline2d - Evaluates the bivariate spline, this is, finds the right polynomial coefficients and use them to compute the approximate F(x,y)
!
! Note: depends on Lib_Kind and for some routines on Lib_Basic.f90 and Lib_Smalapack.f90

use lib_kind
use lib_basic
use lib_smalapack
implicit none
private
public :: linterp, linterpx, bilinterp, bilinterp_vec, trilinterp, tetralinterp,  &
          cheb_grid, cheb_coeff, eval_cheb, spline_nak, spline_naka, spline_nakb, &
          eval_spline, spline2d_nak, eval_spline2d

contains


    subroutine linterp(F,X,x0, f0)
        ! Returns the lineraly interpolated value of F at x0
        ! Usage:
        ! call linterp(F,X,x0, f0)
        ! Inputs:
        ! F      - values of the function at every point in X
        ! X      - points at which F(x) is known
        ! x0     - point at which F will be interpolated
        ! Output:
        ! f0     - linearly interpolated value of F at x0
        real(dp),intent(in)  :: F(:),X(:),x0
        real(dp),intent(out) :: f0
        real(dp)             :: w
        integer              :: n,ind
        n   = size(X)
        ind = gridlookup(X,x0)
        if (ind==0) ind = 1
        if (ind==n) ind = n-1
        w  = (X(ind+1)-x0)/(X(ind+1)-X(ind))
        f0 = w*F(ind) + (1.0_dp-w)*F(ind+1)
    end subroutine linterp
  
  
    real(dp) function linterpx(F,X,x0,extra) result(f0)
        ! Returns the lineraly interpolated value of F at x0
        ! Usage:
        ! f0 = linterpx(F,X,x0[,extra])
        ! Inputs:
        ! F      - values of the function at every point in X
        ! X      - points at which F(x) is known
        ! x0     - point at which F will be interpolated
        ! extra  - extrapolate if out of bounds? if yes, extra=1 (optional)
        ! Output:
        ! f0     - linearly interpolated value of F at x0
        real(dp),intent(in)  :: F(:),X(:),x0
        integer,intent(in),optional :: extra
        real(dp)             :: xd,w,D1,Delta
        integer              :: n,ind,nextra
        nextra = 1
        if (present(extra)) nextra = extra
        n   = size(X)
        ind = gridlookup(X,x0)
        if (ind==0) ind = 1
        if (ind==n) ind = n-1
        if (x0<X(1)) then
            xd = X(1)
        elseif (x0 > X(n)) then
            xd = X(n)
        else
            xd = x0
        end if
        w  = (X(ind+1)-x0)/(X(ind+1)-X(ind))
        if (nextra == 1) then
            Delta = x0 - xd
            D1 = (F(ind+1)-F(ind))/(X(ind+1)-X(ind))
            f0 = w*F(ind) + (1.0_dp-w)*F(ind+1) + D1*Delta
        else
            f0 = w*F(ind) + (1.0_dp-w)*F(ind+1)
        end if
    end function linterpx
  
  
    subroutine bilinterp(F,X,Y,x0,y0, f0)
        ! Returns the bilineraly interpolated value of F at (x0,y0)
        ! Usage:
        ! call bilinterp(F,X,Y,x0,y0, f0)
        ! Inputs:
        ! F      - values of the function at every point in X x Y
        ! X      - points at which F(x,:) is known
        ! Y      - points at which F(:,y) is known
        ! x0     - x-axis point at which F will be interpolated
        ! y0     - y-axis point at which F will be interpolated
        ! Output:
        ! f0     - bilinearly interpolated value of F at (x0,y0)
        real(dp),intent(in)  :: F(:,:),X(:),Y(:),x0,y0
        real(dp),intent(out) :: f0
        real(dp)             :: x1,y1,wx,wY,fx1,fx2,minX,maxX,minY,maxY
        integer              :: nX,nY
        integer              :: indX,indY
        real(dp),parameter   :: err = 1.0e-4_dp
100     format(1x,a,f10.6,a,f10.6,a,f10.6,a)
101     format(1x,a)
        nX   = size(X)
        nY   = size(Y)
        minX = X(1)
        maxX = X(nX)
        minY = Y(1)
        maxY = Y(nY)
        if (x0<minX-err .or. x0>maxX+err) then
            write(*,100) 'Error (bilinterp): x0 is not in X, i.e. ',x0,' \notin [',minX,',',maxX,']'
            write(*,101) 'Automatically bounded (but seriously, check your code)'
        end if
        if (y0<minY-err .or. y0>maxY+err) then
            write(*,100) 'Error (bilinterp): y0 is not in Y, i.e. ',y0,' \notin [',minY,',',maxY,']'
            write(*,101) 'Automatically bounded (but seriously, check your code)'
        end if
        x1   = min(max(x0,minX),maxX)
        y1   = min(max(y0,minY),maxY)
        indX = gridlookup(X,x1)
        indY = gridlookup(Y,y1)
        if (indX==0)  indX = 1
        if (indX==nX) indX = nX-1
        if (indY==0)  indY = 1
        if (indY==nY) indY = nY-1
        wX  = (X(indX+1)-x1)/(X(indX+1)-X(indX))
        wY  = (Y(indY+1)-y1)/(Y(indY+1)-Y(indY))
        fx1 = wX*F(indX,indY)   + (1.0_dp-wX)*F(indX+1,indY)
        fx2 = wX*F(indX,indY+1) + (1.0_dp-wX)*F(indX+1,indY+1)
        f0  = wY*fx1 + (1.0_dp-wY)*fx2
    end subroutine bilinterp

  
    subroutine bilinterp_vec(F,X,Y,x0,y0, f0)
        ! Returns the bilineraly interpolated value of F at (x0,y0)
        ! when F is vectorized. Note that dimension x changes faster than y!
        ! Usage:
        ! call bilinterp_vec(F,X,Y,x0,y0, f0)
        ! Inputs:
        ! F      - values of the function at every point in (X,Y) x 1
        ! X      - points at which F is known
        ! Y      - points at which F is known
        ! x0     - x-axis point at which F will be interpolated
        ! y0     - y-axis point at which F will be interpolated
        ! Output:
        ! f0     - bilinearly interpolated value of F at (x0,y0)
        real(dp),intent(in)  :: F(:),X(:),Y(:),x0,y0
        real(dp),intent(out) :: f0
        real(dp)             :: x1,y1,wx,wY,fx1,fx2,minX,maxX,minY,maxY
        integer              :: nX,nY
        integer              :: indX,indY,ind2
        real(dp),parameter   :: err = 1.0e-4_dp
100     format(1x,a,f10.6,a,f10.6,a,f10.6,a)
101     format(1x,a)
        nX   = size(X)
        nY   = size(Y)
        minX = X(1)
        maxX = X(nX)
        minY = Y(1)
        maxY = Y(nY)
        if (x0<minX-err .or. x0>maxX+err) then
            write(*,100) 'Error (bilinterp): x0 is not in X, i.e. ',x0,' \notin [',minX,',',maxX,']'
            write(*,101) 'Automatically bounded (but seriously, check your code)'
        end if
        if (y0<minY-err .or. y0>maxY+err) then
            write(*,100) 'Error (bilinterp): y0 is not in Y, i.e. ',y0,' \notin [',minY,',',maxY,']'
            write(*,101) 'Automatically bounded (but seriously, check your code)'
        end if
        x1   = min(max(x0,minX),maxX)
        y1   = min(max(y0,minY),maxY)
        indX = gridlookup(X,x1)
        indY = gridlookup(Y,y1)
        if (indX==0)  indX = 1
        if (indX==nX) indX = nX-1
        if (indY==0)  indY = 1
        if (indY==nY) indY = nY-1
        ind2 = (indY-1)*nX + indX
        wX   = (X(indX+1)-x1)/(X(indX+1)-X(indX))
        wY   = (Y(indY+1)-y1)/(Y(indY+1)-Y(indY))
        fx1  = wX*F(ind2)    + (1.0_dp-wX)*F(ind2+1)
        fx2  = wX*F(ind2+nX) + (1.0_dp-wX)*F(ind2+nX+1)
        f0   = wY*fx1 + (1.0_dp-wY)*fx2
    end subroutine bilinterp_vec
  
  
    subroutine trilinterp(F,X,Y,Z,x0,y0,z0, f0)
        ! Returns the tri-lineraly interpolated value of F at (x0,y0,z0)
        ! Usage:
        ! call trilinterp(F,X,Y,Z,x0,y0,z0, f0)
        ! Inputs:
        ! F      - values of the function at every point in X x Y x Z
        ! X      - points at which F(x,:,:) is known
        ! Y      - points at which F(:,y,:) is known
        ! Z      - points at which F(:,:,z) is known
        ! x0     - x-axis point at which F will be interpolated
        ! y0     - y-axis point at which F will be interpolated
        ! z0     - z-axis point at which F will be interpolated
        ! Output:
        ! f0     - bilinearly interpolated value of F at (x0,y0,z0)
        real(dp),intent(in)  :: F(:,:,:),X(:),Y(:),Z(:)
        real(dp),intent(in)  :: x0,y0,z0
        real(dp),intent(out) :: f0
        real(dp)             :: x1,y1,z1,wx,wY,wZ,fx(2,2),fy(2),minX,maxX,minY,maxY,minZ,maxZ
        integer              :: nX,nY,nZ
        integer              :: indX,indY,indZ
        real(dp),parameter   :: err = 1.0e-4_dp
100     format(1x,a,f10.6,a,f10.6,a,f10.6,a)
101     format(1x,a)
        nX   = size(X)
        nY   = size(Y)
        nZ   = size(Z)
        minX = X(1)
        maxX = X(nX)
        minY = Y(1)
        maxY = Y(nY)
        minZ = Z(1)
        maxZ = Z(nZ)
        if (x0<minx-err .or. x0>maxX+err) then
            write(*,100) 'Error (trilinterp): x0 is not in X, i.e. ',x0,' \notin [',minX,',',maxX,']'
            write(*,101) 'Automatically bounded (but seriously, check your code)'
        end if
        if (y0<minY-err .or. y0>maxY+err) then
            write(*,100) 'Error (trilinterp): y0 is not in Y, i.e. ',y0,' \notin [',minY,',',maxY,']'
            write(*,101) 'Automatically bounded (but seriously, check your code)'
        end if
        if (z0<minZ-err .or. z0>maxZ+err) then
            write(*,100) 'Error (trilinterp): z0 is not in Z, i.e. ',z0,' \notin [',minZ,',',maxZ,']'
            write(*,101) 'Automatically bounded (but seriously, check your code)'
        end if
        x1   = min(max(x0,minX),maxX)
        y1   = min(max(y0,minY),maxY)
        z1   = min(max(z0,minZ),maxZ)
        indX = gridlookup(X,x1)
        indY = gridlookup(Y,y1)
        indZ = gridlookup(Z,z1)
        if (indX==0)  indX = 1
        if (indX==nX) indX = nX-1
        if (indY==0)  indY = 1
        if (indY==nY) indY = nY-1
        if (indZ==0)  indZ = 1
        if (indZ==nZ) indZ = nZ-1
        wX  = (X(indX+1)-x1)/(X(indX+1)-X(indX))
        wY  = (Y(indY+1)-y1)/(Y(indY+1)-Y(indY))
        wZ  = (Z(indZ+1)-z1)/(Z(indZ+1)-Z(indZ))
        fx(1,1) = wX*F(indX,indY,indZ)     + (1.0_dp-wX)*F(indX+1,indY,indZ)
        fx(2,1) = wX*F(indX,indY+1,indZ)   + (1.0_dp-wX)*F(indX+1,indY+1,indZ)
        fx(1,2) = wX*F(indX,indY,indZ+1)   + (1.0_dp-wX)*F(indX+1,indY,indZ+1)
        fx(2,2) = wX*F(indX,indY+1,indZ+1) + (1.0_dp-wX)*F(indX+1,indY+1,indZ+1)
        fy(1) = wY*fx(1,1) + (1.0_dp-wY)*fx(2,1)
        fy(2) = wY*fx(1,2) + (1.0_dp-wY)*fx(2,2)
        f0  = wZ*fy(1) + (1.0_dp-wZ)*fy(2)
    end subroutine trilinterp   
  
  
    subroutine tetralinterp(F,X,Y,Z,W,x0,y0,z0,w0, f0)
        ! Returns the tetra-lineraly interpolated value of F at (x0,y0,z0,w0)
        ! Usage:
        ! call tetralinterp(F,X,Y,Z,W,x0,y0,z0,w0, f0)
        ! Inputs:
        ! F      - values of the function at every point in X x Y x Z x W
        ! X      - points at which F(x,:,:,:) is known
        ! Y      - points at which F(:,y,:,:) is known
        ! Z      - points at which F(:,:,z,:) is known
        ! Z      - points at which F(:,:,:,w) is known
        ! x0     - x-axis point at which F will be interpolated
        ! y0     - y-axis point at which F will be interpolated
        ! z0     - z-axis point at which F will be interpolated
        ! z0     - w-axis point at which F will be interpolated
        ! Output:
        ! f0     - bilinearly interpolated value of F at (x0,y0,z0,w0)
        real(dp),intent(in)  :: F(:,:,:,:),X(:),Y(:),Z(:),W(:)
        real(dp),intent(in)  :: x0,y0,z0,w0
        real(dp),intent(out) :: f0
        real(dp)             :: x1,y1,z1,w1,wx,wY,wZ,wW,fx(2,2,2),fy(2,2),fz(2),minX,maxX,minY,maxY,minZ,maxZ,minW,maxW
        integer              :: nX,nY,nZ,nW
        integer              :: indX,indY,indZ,indW        
        real(dp),parameter   :: err = 1.0e-4_dp
100     format(1x,a,f10.6,a,f10.6,a,f10.6,a)
101     format(1x,a)
        nX   = size(X)
        nY   = size(Y)
        nZ   = size(Z)
        nW   = size(W)
        minX = X(1)
        maxX = X(nX)
        minY = Y(1)
        maxY = Y(nY)
        minZ = Z(1)
        maxZ = Z(nZ)
        minW = W(1)
        maxW = W(nW)
        if (x0<minX-err .or. x0>maxX+err) then
            write(*,100) 'Error (tetralinterp): x0 is not in X, i.e. ',x0,' \notin [',minX,',',maxX,']'
            write(*,101) 'Automatically bounded (but seriously, check your code)'
        end if
        if (y0<minY-err .or. y0>maxY+err) then
            write(*,100) 'Error (tetralinterp): y0 is not in Y, i.e. ',y0,' \notin [',minY,',',maxY,']'
            write(*,101) 'Automatically bounded (but seriously, check your code)'
        end if
        if (z0<minZ-err .or. z0>maxZ+err) then
            write(*,100) 'Error (tetralinterp): z0 is not in Z, i.e. ',z0,' \notin [',minZ,',',maxZ,']'
            write(*,101) 'Automatically bounded (but seriously, check your code)'
        end if
        if (w0<minW-err .or. w0>maxW+err) then
            write(*,100) 'Error (tetralinterp): w0 is not in W, i.e. ',w0,' \notin [',minW,',',maxW,']'
            write(*,101) 'Automatically bounded (but seriously, check your code)'
        end if
        x1   = min(max(x0,minX),maxX)
        y1   = min(max(y0,minY),maxY)
        z1   = min(max(z0,minZ),maxZ)
        w1   = min(max(w0,minW),maxW)
        indX = gridlookup(X,x1)
        indY = gridlookup(Y,y1)
        indZ = gridlookup(Z,z1)
        indW = gridlookup(W,w1)
        if (indX==0)  indX = 1
        if (indX==nX) indX = nX-1
        if (indY==0)  indY = 1
        if (indY==nY) indY = nY-1
        if (indZ==0)  indZ = 1
        if (indZ==nZ) indZ = nZ-1
        if (indW==0)  indW = 1
        if (indW==nW) indW = nW-1
        wX  = (X(indX+1)-x1)/(X(indX+1)-X(indX))
        wY  = (Y(indY+1)-y1)/(Y(indY+1)-Y(indY))
        wZ  = (Z(indZ+1)-z1)/(Z(indZ+1)-Z(indZ))
        wW  = (W(indW+1)-w1)/(W(indW+1)-W(indW))
        fx(1,1,1) = wX*F(indX,indY,indZ,indW)       + (1.0_dp-wX)*F(indX+1,indY,indZ,indW)
        fx(2,1,1) = wX*F(indX,indY+1,indZ,indW)     + (1.0_dp-wX)*F(indX+1,indY+1,indZ,indW)
        fx(1,2,1) = wX*F(indX,indY,indZ+1,indW)     + (1.0_dp-wX)*F(indX+1,indY,indZ+1,indW)
        fx(2,2,1) = wX*F(indX,indY+1,indZ+1,indW)   + (1.0_dp-wX)*F(indX+1,indY+1,indZ+1,indW)
        fx(1,1,2) = wX*F(indX,indY,indZ,indW+1)     + (1.0_dp-wX)*F(indX+1,indY,indZ,indW+1)
        fx(2,1,2) = wX*F(indX,indY+1,indZ,indW+1)   + (1.0_dp-wX)*F(indX+1,indY+1,indZ,indW+1)
        fx(1,2,2) = wX*F(indX,indY,indZ+1,indW+1)   + (1.0_dp-wX)*F(indX+1,indY,indZ+1,indW+1)
        fx(2,2,2) = wX*F(indX,indY+1,indZ+1,indW+1) + (1.0_dp-wX)*F(indX+1,indY+1,indZ+1,indW+1)
        fy(1,1) = wY*fx(1,1,1) + (1.0_dp-wY)*fx(2,1,1)
        fy(1,2) = wY*fx(1,1,2) + (1.0_dp-wY)*fx(2,1,2)
        fy(2,1) = wY*fx(1,2,1) + (1.0_dp-wY)*fx(2,2,1)
        fy(2,2) = wY*fx(1,2,2) + (1.0_dp-wY)*fx(2,2,2)
        fz(1) = wZ*fy(1,1) + (1.0_dp-wZ)*fy(2,1)
        fz(2) = wZ*fy(1,2) + (1.0_dp-wZ)*fy(2,2)
        f0  = wW*fz(1) + (1.0_dp-wW)*fz(2)
    end subroutine tetralinterp


    subroutine cheb_grid(a,b,m, Z,X)
        ! Returns the grid that translates from the interval [-1,1] to [a,b]
        ! Usage:
        ! call cheb_grid(a,b,m, Z,X)
        ! Inputs:
        ! a,b    - upper and lower bounds of the interpolated space
        ! m      - order of the polynomial
        ! Y      - points to be interpolated
        ! Outputs:
        ! Z      - zeros of the Chebyshev polynomials on [-1,1]
        ! X      - grid on [a,b]
        real(dp)             :: pi_num = 4.0_dp*atan(1.0_dp)
        integer,intent(in)   :: m
        real(dp),intent(in)  :: a,b
        integer              :: k
        real(dp),intent(out) :: Z(m),X(m)
        ! Retrieving the zeros of the univariate polynomial
        do k = 1,m
            Z(k) = -cos((pi_num/(2*m))*(2*k-1))
            X(k) = (Z(k)+1)*(b-a)/2 + a
        end do
    end subroutine cheb_grid
  
  
    subroutine cheb_coeff(m,Z,Y, Coeffs)
        ! Returns the coefficients fo Chebyshev polynomial interpolation
        ! Usage:
        ! call cheb_coeff(m,Z,Y, Coeffs)
        ! Inputs:
        ! m      - order of the polynomial
        ! Z      - zeros of the univariate polynomial
        ! Y      - points to be interpolated
        ! Output:
        ! Coeffs - chebyshev polynomial coefficients
        integer,intent(in)   :: m
        real(dp),intent(in)  :: Z(m),Y(m)
        real(dp),intent(out) :: Coeffs(m,1)
        integer              :: j
        real(dp)             :: Y1(m,1)
        real(dp)             :: T(m,m),XTX(m,m),XTY(m,1)
        ! Arranging matrices
        T(:,1) = 1
        T(:,2) = Z(:)
        do j = 3,m
            T(:,j) = 2*Z(:)*T(:,j-1) - T(:,j-2)
        end do
        ! Finding the coefficients
        Y1(:,1)  = Y
        XTX = matmul(transpose(T),T)
        XTY = matmul(transpose(T),Y1)
        Coeffs = matmul(inv(XTX),XTY)
    end subroutine cheb_coeff


    subroutine eval_cheb(a,b,m,x,Coeffv, F_App)
        ! Evaluates Chebyshev polynomial to compute the approximate F(x)
        ! Usage:
        ! call eval_cheb(a,b,m,x,Coeffv, F_App)
        ! Inputs:
        ! a,b    - lower and upper bound of the interpolated space
        ! m      - order of the polynomial
        ! x      - point to be evaluated
        ! Coeffv - coefficients of the polynomial
        ! Output:
        ! F_App - approximate value of the function evaluated at x
        integer,intent(in)   :: m
        real(dp),intent(in)  :: a,b,x,Coeffv(:)
        real(dp),intent(out) :: F_App
        integer              :: j
        real(dp)             :: z,Tch(m)        
        z = 2.0_dp*((x-a))/(b-a) - 1.0_dp
        Tch(1) = 1.0_dp
        Tch(2) = z
        do j = 3,m
            Tch(j) = 2.0_dp*z*Tch(j-1) - Tch(j-2)
        end do
        F_App = dot_product(Coeffv,Tch)
    end subroutine eval_cheb


    subroutine spline_nak(X,Fknot,r, Coeff)
        ! Returns the coefficients of a piecewise polynomial cubic spline
        ! using the not-a-knot endpoint condition.
        ! Usage:
        ! call spline_nak(X,Fknot,r, Coeff)
        ! Inputs:
        ! X     - knots
        ! Fknot - function evaluated at the knots
        ! r     - number of interior points
        ! Output:
        ! Coeff - cubic splines coefficients
        integer,intent(in)   :: r
        real(dp),intent(in)  :: X(r+2),Fknot(r+2)
        real(dp),intent(out) :: Coeff(r+1,4)
        integer              :: i
        real(dp)             :: DeltaT(r+1),DivDiff(r+1)
        real(dp)             :: BigT(r,r),BigT2(r,r),Fvec(r,1),omega1,omegar
        real(dp)             :: SlopeSol(r,1),s(r+2)
        !Setting up
        DeltaT  = X(2:r+2) - X(1:r+1)
        DivDiff = (Fknot(2:r+2) - Fknot(1:r+1))/DeltaT
        omega1  = DeltaT(2) - ((DeltaT(1)**2)/DeltaT(2))
        omegar  = DeltaT(r) - ((DeltaT(r+1)**2)/DeltaT(r))
        ! Big T matrix
        BigT = 0.0_dp
        do i = 2,r-1
            BigT(i,i)   = 2.0_dp*(DeltaT(i)+DeltaT(i+1))
            BigT(i,i-1) = DeltaT(i+1)
            BigT(i,i+1) = DeltaT(i)
        end do
        BigT(1,1) = 2.0_dp*(DeltaT(1) + DeltaT(2)) - omega1
        BigT(r,r) = 2.0_dp*(DeltaT(r) + DeltaT(r+1)) - omegar
        BigT(1,2) = DeltaT(1)+((DeltaT(1)**2)/DeltaT(2))
        BigT(r,r-1) = DeltaT(r+1) + ((DeltaT(r+1)**2)/DeltaT(r))
        ! Vector f
        Fvec = 0.0_dp
        Fvec(1,1) = 3.0_dp*(DeltaT(2)*DivDiff(1) + DeltaT(1)*DivDiff(2))         &
                  - 2*(DeltaT(2)*DivDiff(1) - ((DeltaT(1)**2)/DeltaT(2))*DivDiff(2))
        Fvec(r,1) = 3.0_dp*(DeltaT(r+1)*DivDiff(r) + DeltaT(r)*DivDiff(r+1))     &
                  - 2.0_dp*(DeltaT(r)*DivDiff(r+1) - ((DeltaT(r+1)**2)/DeltaT(r))*DivDiff(r))
        do i = 2,r-1
            Fvec(i,1) = 3.0_dp*(DeltaT(i+1)*DivDiff(i) + DeltaT(i)*DivDiff(i+1))
        end do
        ! Solving for the slopes
        BigT2 = BigT
        SlopeSol = matmul(inv(BigT2),Fvec)
        s(2:r+1) = SlopeSol(:,1)
        s(1)   = 2.0_dp*(DivDiff(1) - ((DeltaT(1)/DeltaT(2))**2)*DivDiff(2))     &
               - (1.0_dp-(DeltaT(1)/DeltaT(2))**2)*s(2)                          &
               + ((DeltaT(1)/DeltaT(2))**2)*s(3)
        s(r+2) = 2.0_dp*(DivDiff(r+1) - ((DeltaT(r+1)/DeltaT(r))**2)*DivDiff(r)) &
               - (1.0_dp-(DeltaT(r+1)/DeltaT(r))**2)*s(r+1)                      &
               + ((DeltaT(r+1)/DeltaT(r))**2)*s(r)
        ! Solving for the coefficients
        Coeff(:,1) = Fknot(1:r+1)
        Coeff(:,2) = s(1:r+1)
        Coeff(:,3) = (3.0_dp*DivDiff - 2.0_dp*s(1:r+1) - s(2:r+2))/DeltaT
        Coeff(:,4) = (-2.0_dp*DivDiff + s(1:r+1) + s(2:r+2))/(DeltaT**2)
    end subroutine spline_nak
  
  
    subroutine spline_naka(X,r, Spmat)
        ! Returns the inverse of the matrix of coefficients at each knot
        ! needed for the piecewise polynomial cubic spline (not-a-knot
        ! endpoint condition). This only depends of the knots values, not
        ! on the function values. The outcome of this subroutine is used
        ! by -spline_nakb-.
        ! Usage:
        ! call spline_naka(X,r, Spmat)
        ! Inputs:
        ! X     - knots
        ! r     - number of interior points
        ! Output:
        ! Spmat - coefficients for the system of equations for slopes
        integer,intent(in)   :: r
        real(dp),intent(in)  :: X(r+2)
        real(dp),intent(out) :: Spmat(r,r)
        integer              :: i
        real(dp)             :: DeltaT(r+1)
        real(dp)             :: BigT(r,r),omega1,omegar
        ! Setting up
        DeltaT  = X(2:r+2) - X(1:r+1)
        omega1  = DeltaT(2) - ((DeltaT(1)**2)/DeltaT(2))
        omegar  = DeltaT(r) - ((DeltaT(r+1)**2)/DeltaT(r))
        ! Big T matrix
        BigT = 0.0_dp
        do i = 2,r-1
            BigT(i,i)   = 2.0_dp*(DeltaT(i) + DeltaT(i+1))
            BigT(i,i-1) = DeltaT(i+1)
            BigT(i,i+1) = DeltaT(i)
        end do
        BigT(1,1) = 2.0_dp*(DeltaT(1) + DeltaT(2)) - omega1
        BigT(r,r) = 2.0_dp*(DeltaT(r) + DeltaT(r+1)) - omegar
        BigT(1,2) = DeltaT(1) + ((DeltaT(1)**2)/DeltaT(2))
        BigT(r,r-1) = DeltaT(r+1) + ((DeltaT(r+1)**2)/DeltaT(r))
        ! Solving for the slopes
        Spmat = invlap(BigT)
    end subroutine spline_naka
  
  
    subroutine spline_nakb(Spmat,X,Fknot,r, Coeff)
        ! Returns the coefficients of a piecewise polynomial cubic spline
        ! using the not-a-knot endpoint condition.
        ! Requires the matrix Spmat obtained from -spline_naka-.
        ! Usage:
        ! call spline_nakb(Spmat,X,Fknot,r, Coeff)
        ! Inputs:
        ! Spmat - coefficients for the system of equations for slopes
        ! X     - knots
        ! Fknot - function evaluated at the knots
        ! r     - number of interior points
        ! Output:
        ! Coeff - cubic splines coefficients
        integer,intent(in)   :: r
        real(dp),intent(in)  :: Spmat(r,r),X(r+2),Fknot(r+2)
        real(dp),intent(out) :: Coeff(r+1,4)
        integer              :: i
        real(dp)             :: DeltaT(r+1),DivDiff(r+1)
        real(dp)             :: Fvec(r,1),SlopeSol(r,1),s(r+2)
        ! Setting up
        DeltaT  = X(2:r+2) - X(1:r+1)
        DivDiff = (Fknot(2:r+2) - Fknot(1:r+1))/DeltaT
        ! Vector f
        Fvec = 0.0_dp
        Fvec(1,1) = 3.0_dp*(DeltaT(2)*DivDiff(1) + DeltaT(1)*DivDiff(2))         &
                  - 2.0_dp*(DeltaT(2)*DivDiff(1) - ((DeltaT(1)**2)/DeltaT(2))*DivDiff(2))
        Fvec(r,1) = 3.0_dp*(DeltaT(r+1)*DivDiff(r) + DeltaT(r)*DivDiff(r+1))     &
                  - 2.0_dp*(DeltaT(r)*DivDiff(r+1) - ((DeltaT(r+1)**2)/DeltaT(r))*DivDiff(r))
        do i = 2,r-1
            Fvec(i,1) = 3.0_dp*(DeltaT(i+1)*DivDiff(i) + DeltaT(i)*DivDiff(i+1))
        end do
        ! Solving for the slopes
        SlopeSol = matmul(Spmat,Fvec)
        s(2:r+1) = SlopeSol(:,1)
        s(1)   = 2.0_dp*(DivDiff(1) - ((DeltaT(1)/DeltaT(2))**2)*DivDiff(2))     &
               - (1.0_dp-(DeltaT(1)/DeltaT(2))**2)*s(2)                          &
               + ((DeltaT(1)/DeltaT(2))**2)*s(3)
        s(r+2) = 2.0_dp*(DivDiff(r+1) - ((DeltaT(r+1)/DeltaT(r))**2)*DivDiff(r)) &
               - (1.0_dp-(DeltaT(r+1)/DeltaT(r))**2)*s(r+1)                      &
               + ((DeltaT(r+1)/DeltaT(r))**2)*s(r)
        ! Solving for the coefficients
        Coeff(:,1) = Fknot(1:r+1)
        Coeff(:,2) = s(1:r+1)
        Coeff(:,3) = (3.0_dp*DivDiff - 2.0_dp*s(1:r+1) - s(2:r+2))/DeltaT
        Coeff(:,4) = (-2.0_dp*DivDiff + s(1:r+1) + s(2:r+2))/(DeltaT**2)
    end subroutine spline_nakb

  
    subroutine eval_spline(Coeffm,Knots,r,x,order, F_App)
        ! Evaluates the cubic spline, this is, finds the right polynomial
        ! coefficients and use them to compute the approximate F(x).
        ! If the value to evaluete is outside the bounds of Knots, a
        ! first order (derivative) approximation is performed.
        ! Usage:
        ! call eval_spline(Coeff,Knots,r,x, F_App)
        ! Inputs:
        ! Coeff - r+1x4 matrix of the cubic spline coefficients
        ! Knots - r+2 vector of knots
        ! r     - number of interior points
        ! x     - point to be evaluated
        ! order - order of taylor expansion when x is out of bounds (def=1)
        ! Output:
        ! F_App - approximate value of the function evaluated at x
        integer,intent(in)          :: r
        real(dp),intent(in)         :: Coeffm(r+1,4),Knots(:),x
        integer,intent(in),optional :: order
        real(dp),intent(out)        :: F_App
        real(dp)                    :: D1,D2,Delta,xd
        integer                     :: Ind,j,norder
        norder = 1
        if (present(order)) norder = order
        Ind = gridlookup(Knots,x)
        if (Ind==0)   Ind = 1
        if (Ind==r+2) Ind = r+1
        if (x < Knots(1)) then
            xd = Knots(1)
        elseif (x > Knots(r+2)) then
            xd = Knots(r+2)
        else
            xd = x
        end if
        F_App = 0.0
        D1 = 0.0
        do j = 1,4
            F_App = Coeffm(Ind,j)*(xd-Knots(Ind))**(j-1) + F_App
            if (j >= 2) then
                D1 = (j-1)*Coeffm(Ind,j)*(xd-Knots(Ind))**(j-2) + D1
            end if
        end do
        D2 = 2.0_dp*Coeffm(Ind,3) + 6.0_dp*Coeffm(Ind,4)*(xd-Knots(Ind))
        Delta = x - xd
        if (norder == 1) then
            F_App = F_App + D1*Delta
        elseif (norder == 2) then
            F_App = F_App + D1*Delta + 0.5_dp*D2*Delta**2
        end if
    end subroutine eval_spline
    
  
    subroutine spline2d_nak(X,Y,F,rx,ry, Coeff)
        ! Returns the coefficients of a piecewise polynomial bivariate cubic
        ! spline using the not-a-knot endpoint condition.
        ! Usage:
        ! call spline2d_nak(X,Y,F,rx,ry, Coeff)
        ! Inputs:
        ! X     - knots, rx+2 colmun vector
        ! Y     - knots, ry+2 colmun vector
        ! F     - function evaluated at the knots, r+2 column vector
        ! rx    - number of interior points
        ! ry    - number of interior points
        ! Output:
        ! Coeff - cubic splines coefficients
        integer,intent(in)   :: rx,ry
        real(dp),intent(in)  :: X(rx+2),Y(ry+2),F(rx+2,ry+2)
        real(dp),intent(out) :: Coeff(ry+1,4,4,rx+1)
        integer              :: i,l
        real(dp)             :: Coeffx(rx+1,4,ry+2),Datay(ry+2,4,rx+1)
        ! Fix yknots, and fit spline for x at each yknot
        do i = 1,ry+2
            call spline_nak(X,F(:,i),rx, Coeffx(:,:,i))
        end do
        ! Rearranging step 1 data
        do i = 1,4
            do l = 1,rx+1
                Datay(:,i,l)  = Coeffx(l,i,:)
            end do
        end do
        ! Fitting univariate splines over y for the coefficients of x
        do i = 1,4
            do l = 1,rx+1
                call spline_nak(Y,Datay(:,i,l),ry, Coeff(:,:,i,l))
            end do
        end do
    end subroutine spline2d_nak
  
  
    subroutine eval_spline2d(Coeffm,Knotsx,Knotsy,rx,ry,x,y, F_App2)
        ! Evaluates the bivariate spline, this is, finds the right polynomial
        ! coefficients and use them to compute the approximate F(x,y).
        ! Usage:
        ! call eval_spline2d(Coeffm,Knotsx,Knotsy,rx,ry,x,y, F_App2)
        ! Inputs:
        ! Coeff  - (ry+1,4,4,rx+1) matrix of the cubic spline coefficients
        ! Knotsx - rx+2 vector of knots for x
        ! Knotsy - ry+2 vector of knots for y
        ! rx     - number of interior points of X
        ! ry     - number of interior points of Y
        ! x      - point to be evaluated
        ! Output:
        ! F_App2 - approximate value of the function evaluated at (x,y)
        integer,intent(in)          :: rx,ry
        real(dp),intent(in)         :: Coeffm(ry+1,4,4,rx+1),Knotsx(:),Knotsy(:)
        real(dp),intent(in)         :: x,y
        real(dp),intent(out)        :: F_App2
        real(dp)                    :: D(ry+1,4)
        integer                     :: i,l
        do i = 1,4
            do l = 1,rx+1
                call eval_spline(Coeffm(:,:,i,l),Knotsy,ry,y, F_App=D(l,i))
            end do
        end do
        call eval_spline(D,Knotsx,rx,x, F_App=F_App2)
    end subroutine eval_spline2d
  
  
end module lib_interp

