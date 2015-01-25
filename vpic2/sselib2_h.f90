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
      end module
