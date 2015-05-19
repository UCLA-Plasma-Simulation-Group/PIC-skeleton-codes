!-----------------------------------------------------------------------
      integer function fprecision()
! function determines if default reals are actually doubles
      implicit none
      real :: prec
! ndprec = (0,1) = (no,yes) = (normal,autodouble) precision used
      if (digits(prec) > 24) then
         fprecision = 1
      else
         fprecision = 0
      endif
      end function fprecision
