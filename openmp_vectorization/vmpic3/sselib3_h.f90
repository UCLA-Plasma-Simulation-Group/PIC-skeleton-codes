!-----------------------------------------------------------------------
! Interface file for sselib3.c
      module sselib3_h
      implicit none
!
      interface
         integer function check_sse2()
         implicit none
         end function
      end interface
!
      interface
         integer function check_avx()
         implicit none
         end function
      end interface
!
   end module
