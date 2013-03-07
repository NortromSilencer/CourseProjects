!     
! File:   intsqw.f
! Author: Peisi
!
! Created on March 7, 2013, 12:32 PM
!

!----------------------------------------------------------------------
      Subroutine IntSqw(X,X1,X2,Y1,Y2,S)
      Implicit Real(8) (A-H,O-Z)

      S=0.D0
      if (Y1 == 0.D0) then
        S=Y2
        return
      end if
      if (Y2 == 0.D0) then
        S=Y1
        return
      end if
      S=(Y2-Y1)/(X2-X1)*(X-X1)+Y1
      return
      end
