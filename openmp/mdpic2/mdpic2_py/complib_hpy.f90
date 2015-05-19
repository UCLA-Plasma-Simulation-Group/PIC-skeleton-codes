!-----------------------------------------------------------------------
! Interface file for omplib.c
!     module complib_h
!     implicit none
!
      interface
         subroutine cinit_omp(nth)
         implicit none
         integer, intent(in) :: nth
         end subroutine
      end interface
!
      interface
         subroutine csetnthsize(nth)
         implicit none
         integer, intent(in) :: nth
         end subroutine
      end interface
!
      interface
         function cgetnthsize()
         implicit none
         integer cgetnthsize
         end function
      end interface
!
!     end module
      end