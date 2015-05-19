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
      integer :: j, js, jb, nx, nbx
! nx = size of arrays in x
      nx = min(size(a,1),size(b,1))
! nbx = number of blocks
      nbx = (nx - 1)/mx + 1
!
      do jb = 1, nbx
         do js = 1, min(mx,nx-mx*(jb-1))
            j = js + mx*(jb - 1)
            a(j) = b(j)
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
      integer :: j, k, nx, ny, js, jb, nbx
! nx/ny = size of arrays in x/y
      nx = min(size(a,1),size(b,1)); ny = min(size(a,2),size(b,2))
! nbx = number of blocks in x
      nbx = (nx - 1)/mx + 1
!
      do k = 1, ny
         do jb = 1, nbx
            do js = 1, min(mx,nx-mx*(jb-1))
               j = js + mx*(jb - 1)
               a(j,k) = b(j,k)
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
      integer :: j, k, nx, ny, js, jb, nbx
! nx/ny = size of arrays in x/y
      nx = min(size(a,1),size(b,1)); ny = min(size(a,2),size(b,2))
! nbx = number of blocks in x
      nbx = (nx - 1)/mx + 1
!
      do k = 1, ny
         do jb = 1, nbx
            do js = 1, min(mx,nx-mx*(jb-1))
               j = js + mx*(jb - 1)
               a(j,k) = s*b(j,k) + a(j,k)
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
      integer :: j, k, nx, ny, js, ks, jb, kb, nbx, nby
! nx/ny = size of arrays in x/y
      nx = min(size(a,1),size(b,1)); ny = min(size(a,2),size(b,2))
! nbx/nby = number of blocks in x/y
      nbx = (nx - 1)/mx + 1; nby = (ny - 1)/my + 1
!
      do kb = 1, nby
         do jb = 1, nbx
            do ks = 1, min(my,ny-my*(kb-1))
               k = ks + my*(kb - 1)
               do js = 1, min(mx,nx-mx*(jb-1))
                  j = js + mx*(jb - 1)
                  a(j,k) = b(j,k)
               enddo
            enddo
         enddo
      enddo
!
      end subroutine
!
      end module

