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
      integer :: j, k, nx, ny, js, ks, jb, kb, joff, koff, nbx, nby
! scratch fast, local array
      real, dimension(mx+1,my) :: s
! nx/ny = size of arrays in x/y
      nx = min(size(a,2),size(b,1)); ny = min(size(a,1),size(b,2))
! nbx/nby = number of blocks in x/y
      nbx = (nx - 1)/mx + 1; nby = (ny - 1)/my + 1
!
      do kb = 1, nby
         koff = my*(kb - 1)
         do jb = 1, nbx
            joff = mx*(jb - 1)
!
! copy in from slow memory with stride 1 into block of fast memory
            do ks = 1, min(my,ny-koff)
               k = ks + koff
               do js = 1, min(mx,nx-joff)
                  j = js + joff
                  s(js,ks) = b(j,k)
               enddo
            enddo
!
! copy out to slow memory with stride 1 from block of fast memory
            do js = 1, min(mx,nx-joff)
               j = js + joff
               do ks = 1, min(my,ny-koff)
                  k = ks + koff
                  a(k,j) = s(js,ks)
               enddo
            enddo
!
         enddo
      enddo
!
      end subroutine
!
      end module

