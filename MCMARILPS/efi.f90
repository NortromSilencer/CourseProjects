!     
! File:   efi.f
! Author: Peisi
!
! Created on March 7, 2013, 12:33 PM
!

!----------------------------------------------------------------------
!     Approximate efficiency Function for 3He tube of an inch 
!     in diameter, normal incidence
!----------------------------------------------------------------------
      Real(8) Function EFI(X)
      Implicit Real(8) (A-H,O-Z)
      if (X <= -0.5D0) then
        EFI=1.D0/(1.D0+DEXP(2.D0*(X+0.5D0)))
      else
        EFI=1.D0/(1.D0+DEXP(1.3D0*(X+0.5D0)))
      end if
      return
      end

