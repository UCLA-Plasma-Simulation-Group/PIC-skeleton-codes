!-----------------------------------------------------------------------
! Fortran2003 interface to Unix gettimeofday timer
! written by Viktor K. Decyk, UCLA
      module dtimer_c
      use iso_c_binding
      implicit none
!
      type, bind(C) :: timeval
         integer(c_long) :: tv_sec
         integer(c_long) :: tv_usec
      end type
!
      interface
         subroutine dtimer(time,itime,icntrl) bind(C,name='dtimer')
         use iso_c_binding
         import :: timeval
         implicit none
         real(c_double) :: time
         type (timeval) :: itime
         integer(c_int), value :: icntrl
         end subroutine
      end interface
!
      end module


