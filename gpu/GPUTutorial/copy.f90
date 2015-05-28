!-----------------------------------------------------------------------
! Fortran GPU Tutorial: Copy
! written by Viktor K. Decyk, UCLA
      module copy
      implicit none
!
      contains
!
      subroutine copy0(a,b)
! simple 1d copy of length nx
! a = b
      implicit none
      real, dimension(:) :: a, b
! local data
      integer :: j, nx
! nx = size of arrays in x
      nx = min(size(a,1),size(b,1))
!
      do j = 1, nx
         a(j) = b(j)
      enddo
!
      end subroutine
!
      subroutine copy1(a,b,mx)
! segmented 1d copy of length nx, with block size mx
! a = b
      implicit none
      integer :: mx
      real, dimension(:) :: a, b
! local data
      integer :: j, id, nx, nbx, joff
! nx = size of arrays in x
      nx = min(size(a,1),size(b,1))
! nbx = number of blocks
      nbx = (nx - 1)/mx + 1
!
      do id = 0, nbx-1
         joff = mx*id
         do j = 1, min(mx,nx-joff)
            a(j+joff) = b(j+joff)
         enddo
      enddo
!
      end subroutine
!
      subroutine copy2(a,b,mx)
! segmented 2d copy of length nx, ny, with block size mx
! a = b
      implicit none
      integer :: mx
      real, dimension(:,:) :: a, b
! local data
      integer :: j, k, id, nx, ny, nbx, joff
! nx/ny = size of arrays in x/y
      nx = min(size(a,1),size(b,1)); ny = min(size(a,2),size(b,2))
! nbx = number of blocks in x
      nbx = (nx - 1)/mx + 1
!
      do k = 1, ny
         do id = 0, nbx-1
            joff = mx*id
            do j = 1, min(mx,nx-joff)
               a(j+joff,k) = b(j+joff,k)
            enddo
         enddo
      enddo
!
      end subroutine
!
      subroutine saxpy2(a,b,s,mx)
! segmented 2d vector multiply of length nx, ny, with block size mx
! a = s*b + a
      implicit none
      integer :: mx
      real :: s
      real, dimension(:,:) :: a, b
! local data
      integer :: j, k, id, nx, ny, nbx, joff
! nx/ny = size of arrays in x/y
      nx = min(size(a,1),size(b,1)); ny = min(size(a,2),size(b,2))
! nbx = number of blocks in x
      nbx = (nx - 1)/mx + 1
!
      do k = 1, ny
         do id = 0, nbx-1
            joff = mx*id
            do j = 1, min(mx,nx-joff)
               a(j+joff,k) = s*b(j+joff,k) + a(j+joff,k)
            enddo
         enddo
      enddo
!
      end subroutine
!
      subroutine copy3(a,b,mx,my)
! segmented 2d copy of length nx, ny, with block size mx, my
! a = b
      implicit none
      integer :: mx, my
      real, dimension(:,:) :: a, b
! local data
      integer :: j, k, idx, idy, nx, ny, nbx, nby, joff, koff
! nx/ny = size of arrays in x/y
      nx = min(size(a,1),size(b,1)); ny = min(size(a,2),size(b,2))
! nbx/nby = number of blocks in x/y
      nbx = (nx - 1)/mx + 1; nby = (ny - 1)/my + 1
!
      do idy = 0, nby-1
         koff = my*idy
         do idx = 0, nbx-1
            joff = mx*idx
            do k = 1, min(my,ny-koff)
               do j = 1, min(mx,nx-joff)
                  a(j+joff,k+koff) = b(j+joff,k+koff)
               enddo
            enddo
         enddo
      enddo
!
      end subroutine
!
      end module

