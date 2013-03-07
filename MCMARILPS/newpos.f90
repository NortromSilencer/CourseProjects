!     
! File:   newpos.f
! Author: Peisi
!
! Created on March 7, 2013, 12:41 PM
!

!--------------------------------------------------------------------
!     Calculate the new location of neutron from the previous
!     of the corresponding angles.
!--------------------------------------------------------------------
      Subroutine NewPos(R,XL,D)
      Implicit Real(8) (A-H,O-Z)
      Dimension R(2,3),D(2,3)
      do I=1,3
         R(2,I)=R(1,I)+XL*D(2,I)
      end do
      end

