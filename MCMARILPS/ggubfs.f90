!     
! File:   ggubfs.f
! Author: Peisi
!
! Created on March 7, 2013, 12:40 PM
!

!--------------------------------------------------------------------
!     Generates a random number between 0 and 1.
!--------------------------------------------------------------------
Subroutine GGUBFS(dSeed,randol)
    Implicit Real(8) (A-H,O-Z)
    Dimension randol(1000)
    Data D2P31M/2147483647.D0/
    Data D2P31/2147483711.D0/
    do I=1,1000
        dSeed=DMOD(16807.D0*dSeed,D2P31M)
        R=dSeed/D2P31
        if (R == 0.D0) R=.0000001D0
        randol(I)=R
    end do
    return
end

