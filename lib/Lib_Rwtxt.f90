module lib_rwtxt
   
! Module to export and import data to and from a text file.
! By: Christian Bustamante
! Email: mail@cbustamante.co
! Last modified: 16 Oct 2020, 17:37
! v 2.21
!
! To compile this module:
! $ ifort -c Lib_Rwtxt.f90 Lib_Kind.f90
! $ gfortran -c Lib_Rwtxt.f90 Lib_Kind.f90
! To use it put this just after the -program- statement:
! use lib_rwtxt
!
!
! Subroutines:
! (Rwtxt) 01. txt_write      - Writes a n-dimensional array (n=1,2,3,4,5,6 including scalars) in the specified text file [Interface]
! (Rwtxt) 02. txt_read       - Reads a n-dimensional array (n=1,2,3,4,5,6 including scalars) from the specified text file [Interface]
! (Rwtxt) 03. txtw_data_scl  - Writes a scalar in the specified text file
! (Rwtxt) 04. txtw_data_vec  - Writes a vector in the specified text file
! (Rwtxt) 05. txtw_data_mat  - Writes a matrix in the specified text file
! (Rwtxt) 06. txtw_data_arr3 - Writes a 3-dimensional array in the specified text file
! (Rwtxt) 07. txtw_data_arr4 - Writes a 4-dimensional array in the specified text file
! (Rwtxt) 08. txtw_data_arr5 - Writes a 5-dimensional array in the specified text file
! (Rwtxt) 09. txtw_data_arr6 - Writes a 6-dimensional array in the specified text file
! (Rwtxt) 10. txtr_data_scl  - Reads an scalar from the specified text file
! (Rwtxt) 11. txtr_data_vec  - Reads an vector from the specified text file
! (Rwtxt) 12. txtr_data_mat  - Reads an matrix from the specified text file
! (Rwtxt) 13. txtr_data_arr3 - Reads an 3-dimensional array from the specified text file
! (Rwtxt) 14. txtr_data_arr4 - Reads an 4-dimensional array from the specified text file   
! (Rwtxt) 15. txtr_data_arr5 - Reads an 5-dimensional array from the specified text file   
! (Rwtxt) 16. txtr_data_arr6 - Reads an 6-dimensional array from the specified text file   
!
! Note: depends on Lib_Kind 
   
use lib_kind
implicit none
private
public :: txt_write, txt_read

interface txt_write
    ! Writes a n-dimensional array (n=1,2,3,4, including scalars) in the specified text file.
    ! Supports double precission arrays only - real(dp). Integers should be converted into
    ! dp with the function real(*,dp).
    ! Usage:
    ! call txtw_data(wdata,filename)
    ! Inputs:
    ! wdata    - data (n-dimensional array)
    ! filename - name of text file (with explicit file extension)
    module procedure txtw_data_scl
    module procedure txtw_data_vec
    module procedure txtw_data_mat
    module procedure txtw_data_arr3
    module procedure txtw_data_arr4
    module procedure txtw_data_arr5
    module procedure txtw_data_arr6
end interface

interface txt_read
    ! Reads a n-dimensional array (n=1,2,3,4, including scalars) from the specified text file.
    ! The text file is expected to be written by txt_write
    ! Dimensions are implicitly specified by the dimensions of the -inout- variable
    ! Supports double precission arrays only - real(dp). Integers should be converted into
    ! dp with the function real(*,dp).
    ! Usage:
    ! call txtw_data(wdata,filename)
    ! Inputs:
    ! filename - name of text file (with explicit file extension)
    ! rdata    - array that was read
    ! Note: note that when reading, you do not need to open and close the h5 file!
    module procedure txtr_data_scl
    module procedure txtr_data_vec
    module procedure txtr_data_mat
    module procedure txtr_data_arr3
    module procedure txtr_data_arr4
    module procedure txtr_data_arr5
    module procedure txtr_data_arr6
end interface
   
contains
   
   
    subroutine txtw_data_scl(wdata,filename)
        ! Writes a scalar in the specified text file.
        ! Belongs to the interface -txt_data-. See usage there.
        real(dp),intent(in)          :: wdata
        character (len=*),intent(in) :: filename
        character (60)               :: fnametxt
        fnametxt = trim(filename)
        open(11, file = fnametxt)
        write (11,*) wdata
        close(11)
    end subroutine txtw_data_scl
   
   
    subroutine txtw_data_vec(wdata,filename)
        ! Writes a 1-dimensional array in the specified text file.
        ! Belongs to the interface -txt_data-. See usage there.
        real(dp),intent(in)          :: wdata(:)
        character (len=*),intent(in) :: filename
        character (60)               :: fnametxt
        integer                      :: i
        fnametxt = trim(filename)
        open(11, file = fnametxt)
        do i = 1,size(wdata,1)
            write (11,*) wdata(i)
        end do
        close(11)
    end subroutine txtw_data_vec
   
   
    subroutine txtw_data_mat(wdata,filename)
        ! Writes a 2-dimensional array in the specified text file.
        ! Belongs to the interface -txt_data-. See usage there.
        real(dp),intent(in)          :: wdata(:,:)
        character (len=*),intent(in) :: filename
        character (60)               :: fnametxt
        integer                      :: i
        fnametxt = trim(filename)
        open(11, file = fnametxt)
        do i = 1,size(wdata,1)
            write (11,"(100000(f24.16,:))") wdata(i,:)
        end do
        close(11)
    end subroutine txtw_data_mat
   
   
    subroutine txtw_data_arr3(wdata,filename)
        ! Writes a 3-dimensional array in the specified text file.
        ! Belongs to the interface -txt_data-. See usage there.
        real(dp),intent(in)          :: wdata(:,:,:)
        character (len=*),intent(in) :: filename
        character (60)               :: fnametxt
        integer                      :: i,i2,i3,n1,n2,n3,icount
        integer, allocatable         :: dimvec(:)
        real(dp),allocatable         :: wvec(:)
        fnametxt = trim(filename)
        n1 = size(wdata,1)
        n2 = size(wdata,2)
        n3 = size(wdata,3)
        allocate(wvec(n1*n2*n3), dimvec(n1*n2*n3))
        ! Rearranging
        icount = 1
        do i3 = 1,n3
            do i2 = 1,n2
                wvec(icount:icount+n1-1) = wdata(:,i2,i3)
                icount = icount+n1
            end do
        end do
        dimvec = 0
        dimvec(1:3) = [n1,n2,n3]
        ! Writing
        open (11, file = fnametxt)
        do i = 1,n1*n2*n3
            write (11,"(f24.16,a,i10)") wvec(i), " ", dimvec(i)
        end do
        close(11)
    end subroutine txtw_data_arr3
   
   
    subroutine txtw_data_arr4(wdata,filename)
        ! Writes a 4-dimensional array in the specified text file.
        ! Belongs to the interface -txt_data-. See usage there.
        real(dp),intent(in)          :: wdata(:,:,:,:)
        character (len=*),intent(in) :: filename
        character (60)               :: fnametxt
        integer                      :: i,i2,i3,i4,n1,n2,n3,n4,icount
        integer, allocatable         :: dimvec(:)
        real(dp),allocatable         :: wvec(:)
        fnametxt = trim(filename)
        n1 = size(wdata,1)
        n2 = size(wdata,2)
        n3 = size(wdata,3)
        n4 = size(wdata,4)
        allocate(wvec(n1*n2*n3*n4), dimvec(n1*n2*n3*n4))
        ! Rearranging
        icount = 1
        do i4 = 1,n4
            do i3 = 1,n3
                do i2 = 1,n2
                    wvec(icount:icount+n1-1) = wdata(:,i2,i3,i4)
                    icount = icount+n1
                end do
            end do
        end do
        dimvec = 0
        dimvec(1:4) = [n1,n2,n3,n4]
        ! Writing
        open (11, file = fnametxt)
        do i = 1,n1*n2*n3*n4
            write (11,"(f24.16,i10)") wvec(i), dimvec(i)
        end do
        close(11)
    end subroutine txtw_data_arr4
  
  
    subroutine txtw_data_arr5(wdata,filename)
        ! Writes a 5-dimensional array in the specified text file.
        ! Belongs to the interface -txt_data-. See usage there.
        real(dp),intent(in)          :: wdata(:,:,:,:,:)
        character (len=*),intent(in) :: filename
        character (60)               :: fnametxt
        integer                      :: i,i2,i3,i4,i5,n1,n2,n3,n4,n5,icount
        integer, allocatable         :: dimvec(:)
        real(dp),allocatable         :: wvec(:)
        fnametxt = trim(filename)
        n1 = size(wdata,1)
        n2 = size(wdata,2)
        n3 = size(wdata,3)
        n4 = size(wdata,4)
        n5 = size(wdata,5)
        allocate(wvec(n1*n2*n3*n4*n5), dimvec(n1*n2*n3*n4*n5))
        ! Rearranging
        icount = 1
        Do i5 = 1,n5
            do i4 = 1,n4
                do i3 = 1,n3
                    do i2 = 1,n2
                        wvec(icount:icount+n1-1) = wdata(:,i2,i3,i4,i5)
                        icount = icount+n1
                    end do
                end do
            end do
        end do
        dimvec = 0
        dimvec(1:5) = [n1,n2,n3,n4,n5]
        ! Writing
        open (11, file = fnametxt)
        do i = 1,n1*n2*n3*n4*n5
            write (11,"(f24.16,a,i10)") wvec(i), " ", dimvec(i)
        end do
        close(11)
    end subroutine txtw_data_arr5
  
  
    subroutine txtw_data_arr6(wdata,filename)
        ! Writes a 6-dimensional array in the specified text file.
        ! Belongs to the interface -txt_data-. See usage there.
        real(dp),intent(in)          :: wdata(:,:,:,:,:,:)
        character (len=*),intent(in) :: filename
        character (60)               :: fnametxt
        integer                      :: i,i2,i3,i4,i5,i6,n1,n2,n3,n4,n5,n6,icount
        integer, allocatable         :: dimvec(:)
        real(dp),allocatable         :: wvec(:)
        fnametxt = trim(filename)
        n1 = size(wdata,1)
        n2 = size(wdata,2)
        n3 = size(wdata,3)
        n4 = size(wdata,4)
        n5 = size(wdata,5)
        n6 = size(wdata,6)
        allocate(wvec(n1*n2*n3*n4*n5*n6), dimvec(n1*n2*n3*n4*n5*n6))
        ! Rearranging
        icount = 1
        do i6 = 1,n6
            do i5 = 1,n5
                do i4 = 1,n4
                    do i3 = 1,n3
                        do i2 = 1,n2
                            wvec(icount:icount+n1-1) = wdata(:,i2,i3,i4,i5,i6)
                            icount = icount+n1
                        end do
                    end do
                end do
            end do
        end do
        dimvec = 0
        dimvec(1:6) = [n1,n2,n3,n4,n5,n6]
        ! Writing
        open (11, file = fnametxt)
        do i = 1,n1*n2*n3*n4*n5*n6
            write (11,"(f24.16,a,i10)") wvec(i), " ", dimvec(i)
        end do
        close(11)
    end subroutine txtw_data_arr6
  
  
    subroutine txtr_data_scl(filename,rdata)
        ! Reads an scalar from the specified text file. The dataset is
        ! read as real(dp).
        ! Belongs to the interface -txtr_data-. See usage there.
        character (len=*),intent(in) :: filename
        real(dp),intent(inout)       :: rdata
        character (60)               :: fnametxt
        fnametxt = trim(filename)
        open (12, file = fnametxt)
        read (12,*) rdata
        close(12)
    end subroutine txtr_data_scl
  
  
    subroutine txtr_data_vec(filename,rdata)
        ! Reads an 1-dimensional array from the specified text file. The
        ! dataset is read as real(dp).
        ! Belongs to the interface -txtr_data-. See usage there.
        character (len=*),intent(in) :: filename
        real(dp),intent(inout)       :: rdata(:)
        integer                      :: n1,i
        character (60)               :: fnametxt
        fnametxt = trim(filename)
        n1 = size(rdata,1)
        open (12, file = fnametxt)
        read (12,*) (rdata(i), i=1,n1)
        close(12)
    end subroutine txtr_data_vec
  
  
    subroutine txtr_data_mat(filename,rdata)
        ! Reads an 2-dimensional array from the specified text file. The
        ! dataset is read as real(dp).
        ! Belongs to the interface -txtr_data-. See usage there.
        character (len=*),intent(in) :: filename
        real(dp),intent(inout)       :: rdata(:,:)
        integer                      :: n1,n2,i1,i2
        character (60)               :: fnametxt
        fnametxt = trim(filename)
        n1 = size(rdata,1)
        n2 = size(rdata,2)
        open (12, file = fnametxt)
        read (12,*) ((rdata(i1,i2), i2=1,n2), i1=1,n1)
        close(12)
    end subroutine txtr_data_mat
  
  
    subroutine txtr_data_arr3(filename,rdata)
        ! Reads an 3-dimensional array from the specified text file. The
        ! dataset is read as real(dp).
        ! Belongs to the interface -txtr_data-. See usage there.
        character (len=*),intent(in) :: filename
        real(dp),intent(inout)       :: rdata(:,:,:)
        integer                      :: i,i2,i3,n1,n2,n3,icount,temp
        real(dp),allocatable         :: rvec(:)
        character (60)               :: fnametxt
        fnametxt = trim(filename)
        n1 = size(rdata,1)
        n2 = size(rdata,2)
        n3 = size(rdata,3)
        allocate(rvec(n1*n2*n3))
        ! reading vectorized array
        open (12, file = fnametxt)
        do i = 1,n1*n2*n3
            read (12,"(f24.16,i10)") rvec(i), temp
        end do
        close(12)
        ! rearranging
        icount = 1
        do i3 = 1,n3
            do i2 = 1,n2
                rdata(:,i2,i3) = rvec(icount:icount+n1-1)
                icount = icount+n1
            end do
        end do
    end subroutine txtr_data_arr3
  
  
    subroutine txtr_data_arr4(filename,rdata)
        ! Reads an 4-dimensional array from the specified text file. The
        ! dataset is read as real(dp).
        ! Belongs to the interface -txtr_data-. See usage there.
        character (len=*),intent(in) :: filename
        real(dp),intent(inout)       :: rdata(:,:,:,:)
        integer                      :: i,i2,i3,i4,n1,n2,n3,n4,icount,temp
        real(dp),allocatable         :: rvec(:)
        character (60)               :: fnametxt
        fnametxt = trim(filename)
        n1 = size(rdata,1)
        n2 = size(rdata,2)
        n3 = size(rdata,3)
        n4 = size(rdata,4)
        allocate(rvec(n1*n2*n3*n4))
        ! reading vectorized array
        open (12, file = fnametxt)
        do i = 1,n1*n2*n3*n4
            read (12,"(f24.16,i10)") rvec(i), temp
        end do
        close(12)
        ! rearranging
        icount = 1
        do i4 = 1,n4
            do i3 = 1,n3
                do i2 = 1,n2
                    rdata(:,i2,i3,i4) = rvec(icount:icount+n1-1)
                    icount = icount+n1
                end do
            end do
        end do
    end subroutine txtr_data_arr4
  
  
    subroutine txtr_data_arr5(filename,rdata)
        ! Reads an 5-dimensional array from the specified text file. The
        ! dataset is read as real(dp).
        ! Belongs to the interface -txtr_data-. See usage there.
        character (len=*),intent(in) :: filename
        real(dp),intent(inout)       :: rdata(:,:,:,:,:)
        integer                      :: i,i2,i3,i4,i5,n1,n2,n3,n4,n5,icount,temp
        real(dp),allocatable         :: rvec(:)
        character (60)               :: fnametxt
        fnametxt = trim(filename)
        n1 = size(rdata,1)
        n2 = size(rdata,2)
        n3 = size(rdata,3)
        n4 = size(rdata,4)
        n5 = size(rdata,5)
        allocate(rvec(n1*n2*n3*n4*n5))
        ! reading vectorized array
        open (12, file = fnametxt)
        do i = 1,n1*n2*n3*n4*n5
            read (12,"(f24.16,i10)") rvec(i), temp
        end do
        close(12)
        ! rearranging
        icount = 1
        do i5 = 1,n5
            do i4 = 1,n4
                do i3 = 1,n3
                    do i2 = 1,n2
                        rdata(:,i2,i3,i4,i5) = rvec(icount:icount+n1-1)
                        icount = icount+n1
                    end do
                end do
            end do
        end do
    end subroutine txtr_data_arr5
  
    subroutine txtr_data_arr6(filename,rdata)
        ! Reads an 6-dimensional array from the specified text file. The
        ! dataset is read as real(dp).
        ! Belongs to the interface -txtr_data-. See usage there.
        character (len=*),intent(in) :: filename
        real(dp),intent(inout)       :: rdata(:,:,:,:,:,:)
        integer                      :: i,i2,i3,i4,i5,i6,n1,n2,n3,n4,n5,n6,icount,temp
        real(dp),allocatable         :: rvec(:)
        character (60)               :: fnametxt
        fnametxt = trim(filename)
        n1 = size(rdata,1)
        n2 = size(rdata,2)
        n3 = size(rdata,3)
        n4 = size(rdata,4)
        n5 = size(rdata,5)
        n6 = size(rdata,6)
        allocate(rvec(n1*n2*n3*n4*n5*n6))
        ! reading vectorized array
        open (12, file = fnametxt)
        do i = 1,n1*n2*n3*n4*n5*n6
            read (12,"(f24.16,i10)") rvec(i), temp
        end do
        close(12)
        ! rearranging
        icount = 1
        do i6 = 1,n6
            do i5 = 1,n5
                do i4 = 1,n4
                    do i3 = 1,n3
                        do i2 = 1,n2
                            rdata(:,i2,i3,i4,i5,i6) = rvec(icount:icount+n1-1)
                            icount = icount+n1
                        end do
                    end do
                end do
            end do
        end do
    end subroutine txtr_data_arr6
   
   
end module lib_rwtxt
