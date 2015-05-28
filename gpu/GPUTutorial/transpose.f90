!-----------------------------------------------------------------------
! Fortran Library GPU Tutorial
! written by Viktor K. Decyk, UCLA
      module transpose
      implicit none
!
      contains
!
      subroutine transpose0(a,b)
! simple 2d transpose of length nx, ny
! a = transpose(b)
      implicit none
      real, dimension(:,:) :: a, b
! local data
      integer :: j, k, nx, ny
! nx/ny = size of arrays in x/y
      nx = min(size(a,2),size(b,1)); ny = min(size(a,1),size(b,2))
!
      do k = 1, ny
         do j = 1, nx
            a(k,j) = b(j,k)
         enddo
      enddo
!
      end subroutine
!
      subroutine transpose2(a,b,mx,my)
! segmented 2d transpose of length nx, ny, with block size mx, my
! a = transpose(b)
      implicit none
      integer :: mx, my
      real, dimension(:,:) :: a, b
! local data
      integer :: j, k, nx, ny, idx, idy, joff, koff, nbx, nby
! scratch fast, local array
      real, dimension(mx+1,my) :: s
! nx/ny = size of arrays in x/y
      nx = min(size(a,2),size(b,1)); ny = min(size(a,1),size(b,2))
! nbx/nby = number of blocks in x/y
      nbx = (nx - 1)/mx + 1; nby = (ny - 1)/my + 1
!
      do idy = 0, nby-1
         koff = my*idy
         do idx = 0, nbx-1
            joff = mx*idx
!
! copy in from slow memory with stride 1 into block of fast memory
            do k = 1, min(my,ny-koff)
               do j = 1, min(mx,nx-joff)
                  s(j,k) = b(j+joff,k+koff)
               enddo
            enddo
!
! copy out to slow memory with stride 1 from block of fast memory
            do j = 1, min(mx,nx-joff)
               do k = 1, min(my,ny-koff)
                  a(k+koff,j+joff) = s(j,k)
               enddo
            enddo
!
         enddo
      enddo
!
      end subroutine
!
      end module

