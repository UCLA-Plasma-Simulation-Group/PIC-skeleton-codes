!-----------------------------------------------------------------------
! Interface file for dtimer.c
!     module dtimer_h
!     implicit none
!
      interface
         subroutine dtimer(dtime,itime,icntrl)
         implicit none
         double precision, intent(inout) :: dtime
         integer, dimension(4), intent(inout) :: itime
         integer, intent(in) :: icntrl
         end subroutine
      end interface
!
!     end module
      end
