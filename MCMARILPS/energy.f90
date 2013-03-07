!     
! File:   energy.f
! Author: Peisi
!
! Created on March 7, 2013, 12:39 PM
!

!--------------------------------------------------------------------
!     Calculates the energy transfer When a collision occurs
!--------------------------------------------------------------------
      Subroutine ENERGY(IMED,U,ifLAG,RANdo)
      Implicit Real(8) (A-H,O-Z)
      Parameter (NQEXP=250,NWEXP=550,NQEXPC=250,NWEXPC=550)
      Common /EXPCOORD/QEXP(NQEXP),HWEXP(NWEXP)
      Common /EXPCOORDC/QEXPC(NQEXPC),HWEXPC(NWEXPC)
      Common /EXPERSAM/DATSAM(NQEXP,NWEXP)
      Common /EXPERCAN/DATCAN(NQEXPC,NWEXPC)
      Dimension endIST(NWEXP),endISTC(NWEXPC)
      ifLAG=1
      STTT=0.D0
      S1=0.D0
      E0=10.D0**U
      call EtKernel(IMED,E0,endIST,endISTC)
      if (IMED == 0) then
        N=NWEXP
        DX=HWEXP(2)-HWEXP(1)
      else
        N=NWEXPC
        DX=HWEXPC(2)-HWEXPC(1)
      endif
      do 100 I=1,N
      if (IMED == 0) then
         S2=endIST(I)
         else
         S2=endISTC(I)
      endif
      STTT=.5D0*(S1+S2)*DX+STTT
100   S1=S2
      if (STTT == 0.D0) then
        ifLAG=0
        return
      endif
      S1=0.D0
      S=0.D0
      do 11 I=1,N
      if (IMED == 0) then
         S2=endIST(I)
         else
         S2=endISTC(I)
      endif
      S2=S2/STTT
      DS=.5D0*(S1+S2)*DX
      if (S+DS-RANdo)33,22,22
33    S=S+DS
      S1=S2
11    CONTINUE
      return
22    if (IMED == 0) then
      A=(S2-S1)/DX
      if (a.eq.0.d0) then
       hw=rando*dx+hwexp(i-1)
       goto 44
      end if
      B=(S1*HWEXP(I)-S2*HWEXP(I-1))/DX
      BA=B/A
      RAD=BA*BA+2.D0*(RANdo-S)/A+(HWEXP(I-1))**2.D0+2*BA*(HWEXP(I-1))
      HW=-BA+DSQRT(RAD)
      if (.not.(HW > HWEXP(I-1) .and. HW < HWEXP(I))) HW=-BA-DSQRT(RAD)
44    E=E0-HW
      if(E <= 0.D0) E=E0-HWEXP(I-1)
      if(E <= 0.D0) E=1.D-4
      U=DLOG10(E)
      end if
      if (IMED == 1) then
      A=(S2-S1)/DX
       if (a.eq.0.d0) then
       hw=rando*dx+hwexp(i-1)
       goto 55
      end if
      B=(S1*HWEXPC(I)-S2*HWEXPC(I-1))/DX
      BA=B/A
      RAD=BA*BA+2.D0*(RANdo-S)/A+(HWEXPC(I-1))**2.D0+2*BA*(HWEXPC(I-1))
      HW=-BA+DSQRT(RAD)
      if(.not.(HW > HWEXPC(I-1) .and. HW < HWEXPC(I))) HW=-BA-DSQRT(RAD)
55    E=E0-HW
      if(E <= 0.D0) E=E0-HWEXPC(I-1)
      if(E <= 0.D0) E=1.D-4
      U=DLOG10(E)
      end if
      return
      end

