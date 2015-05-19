!-----------------------------------------------------------------------
! Interface file for omplib.f
      module omplib_h
      implicit none
!
      interface
         subroutine INIT_OMP(nth)
         implicit none
         integer :: nth
         end subroutine
      end interface
!
      interface
         subroutine SETNTHSIZE(nth)
         implicit none
         integer :: nth
         end subroutine
      end interface
!
      interface
         function GETNTHSIZE()
         implicit none
         integer GETNTHSIZE
         end function
      end interface
!
      end module