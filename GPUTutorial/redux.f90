!-----------------------------------------------------------------------
! Fortran GPU Tutorial: Reduction
! written by Viktor K. Decyk, UCLA
      module redux
      implicit none
!
      contains
!
      subroutine sum0(a,sa)
! simple 1d sum reduction of length nx
! sa = sum(a)
      implicit none
      real :: sa
      real, dimension(:) :: a
! local data
      integer :: j
!
      sa = 0.0
      do j = 1, size(a,1)
         sa = sa + a(j)
      enddo
!
      end subroutine
!
      subroutine sum1(a,sa,mx)
! 1d sum reductions, each of length mx
! sa = sum(a)
      implicit none
      integer :: mx
      real :: sa
      real, dimension(:) :: a
! local data
      integer :: j, js, jb, nx, nbx
      real :: t
! nx = size of a array
      nx = size(a,1)
! nbx = number of blocks
      nbx = (nx - 1)/mx + 1
!
      sa = 0.0
      do jb = 1, nbx
      t = 0.0
      do js = 1, min(mx,nx-mx*(jb-1))
         j = js + mx*(jb - 1)
         t = t + a(j)
      enddo
      sa = sa + t

      enddo
!
      end subroutine
!
      subroutine sum2(a,d,mx)
! segmented 1d sum reductions, each of length mx
! forall (j = 1:nbx); d(j) = sum(a(1+mx*(j-1):min(nx,mx*j))); end forall
      implicit none
      integer :: mx
      real, dimension(:) :: a, d
! local data
      integer :: j, js, jb, nx, nbx
      real :: t
! nx = size of a array
      nx = size(a,1)
! nbx = number of blocks
      nbx = (nx - 1)/mx + 1
!
      do jb = 1, nbx
      t = 0.0
      do js = 1, min(mx,nx-mx*(jb-1))
         j = js + mx*(jb - 1)
         t = t + a(j)
      enddo
      d(jb) = t
      enddo
!
      end subroutine
!
      end module

