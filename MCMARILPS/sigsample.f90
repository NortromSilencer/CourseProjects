!     
! File:   sigsample.f
! Author: Peisi
!
! Created on March 7, 2013, 12:35 PM
!

!--------------------------------------------------------------------
!
!  Calculate the value of the Function
!
!
!                 P(E,E',O) * exp( -s (E')*D ) Ef(E')
!                                    T
!  P is the probability of scattering based on experimental Data. 
!  The normalized by the effective experimental Section out
!
!  where E 'is the output energy.
!
!  ifLAG=0 : Calculated without correction for attenuation and efficiency
!            detector
!  ifLAG other than 0: complete calculation.
!--------------------------------------------------------------------

      Real(8) Function SIGSAMPLE(ifLAG,N,NQ,X)
      Implicit Real(8) (A-H,O-Z)
      Implicit Integer(4) (I-N)
      Parameter (NQEXP=250,NWEXP=550)
      Common /DATOSN/EINC,U,theta
      Common /DISTS/DISCE1,DISAG1,DISCI1,DISVA,DISCI2,DISAG2,DISCE2
      Common /DATOSSAM/SE,SE0,SM,PS,UTAB(68),SETOT(68),SEM(68),PSCAT(68)
      Common /DATOSCAN/SEC,SE0C,SMC,PSC,SETOTC(68),SEMC(68),PSCATC(68)
      Common /EXPCOORD/QEXP(NQEXP),HWEXP(NWEXP)
      Common /EXPERSAM/DATSAM(NQEXP,NWEXP)
      Data HQS2M/2.07219093D-03/
      EACT=10.D0**U
      E=10.D0**X
      call VALINTERP(U)
      SEACT=SE
      call VALINTERP(X)
      HW=EACT-E
      
! Calculate the real Q for this final energy and angle
      if (N == 1) then
          Q=QEXP(NQ)
        else
          Q=DSQRT((E+EACT-2.D0*DSQRT(E*EACT)*DCOS(theta))/HQS2M)
      end if

      call Locate(QEXP,NQEXP,Q,JQ)
      call Locate(HWEXP,NWEXP,HW,JW)
      S2Dif=0.D0
      if(JQ >= 1 .and. JQ <= NQEXP-1 .and. JW >= 1 .and. JW <= NWEXP-1) then
          call IntSqw(Q,QEXP(JQ),QEXP(JQ+1),DATSAM(JQ,JW),DATSAM(JQ+1,JW),S1)
          call IntSqw(Q,QEXP(JQ),QEXP(JQ+1),DATSAM(JQ,JW+1),DATSAM(JQ+1,JW+1),S2)
          call IntSqw(HW,HWEXP(JW),HWEXP(JW+1),S1,S2,S2Dif)
      end if
      S2Dif=S2Dif*(SE0/SEACT)*DSQRT(EINC/EACT)
      if (ifLAG == 0) then
          SIGSAMPLE=S2Dif
          return
      end if

      TRANS=DEXP(-SMC*(DISCE1+DISCE2)-SM*(DISAG1+DISAG2)-SMC*(DISCI1+DISCI2))
      SIGSAMPLE=S2Dif*TRANS*EFI(X)

      return
      end

