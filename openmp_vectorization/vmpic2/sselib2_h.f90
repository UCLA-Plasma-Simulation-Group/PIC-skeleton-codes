!-----------------------------------------------------------------------
! Interface file for sselib2.c
      module sselib2_h
      implicit none
!
      interface
         subroutine csse2iscan2(isdata,nths)
         implicit none
         integer :: nths
         integer, dimension(nths) :: isdata
         end subroutine
      end interface
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
