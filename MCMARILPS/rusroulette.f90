!     
! File:   rusroulette.f
! Author: Peisi
!
! Created on March 7, 2013, 12:36 PM
!

!--------------------------------------------------------------------
!     Russian Roulette
!--------------------------------------------------------------------
      Subroutine RUSROULETTE(randol,WGT,WCO,ifLAG)
      Implicit Real(8) (A-H,O-Z)
      ifLAG=1
      if (wgt.eq.0.d0) then
	iflag=0
	return
      end if
      if (WGT < WCO) then
        if (randol < .5D0) then
           WGT=WGT*2.D0
          else
           ifLAG=0
        end if
      end if
      return
      end

