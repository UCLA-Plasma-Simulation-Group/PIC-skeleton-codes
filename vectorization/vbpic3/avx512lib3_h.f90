!-----------------------------------------------------------------------
! Interface file for avx512lib3.c
      module avx512lib3_h
      implicit none
!
      interface
         subroutine cknciscan2(isdata,nths)
         implicit none
         integer :: nths
         integer, dimension(nths) :: isdata
         end subroutine
      end interface
!
      end module
