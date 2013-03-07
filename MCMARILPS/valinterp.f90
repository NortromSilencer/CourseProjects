!     
! File:   valinterp.f
! Author: Peisi
!
! Created on March 7, 2013, 12:39 PM
!

!--------------------------------------------------------------------
!     Interpolated input Data table
!--------------------------------------------------------------------
      Subroutine VALINTERP(U)
      Implicit Real(8) (A-H,O-Z)
      Common /DATOSSAM/SE,SE0,SM,PS,UTAB(68),SETOT(68),SEM(68),PSCAT(68)
      Common /DATOSCAN/SEC,SE0C,SMC,PSC,SETOTC(68),SEMC(68),PSCATC(68)
      do 1 I=1,68
      if(UTAB(I)-U)1,2,2
1     CONTINUE
      I=68
2     CONTINUE
      if (I == 1) then
        SE=SETOT(1)
        SM=SEM(1)
        PS=PSCAT(1)
        SEC=SETOTC(1)
        SMC=SEMC(1)
        PSC=PSCATC(1)
        return
      end if
      X=(U-UTAB(I-1))/(UTAB(I)-UTAB(I-1))
      SE= X*(SETOT(I)-SETOT(I-1))+SETOT(I-1)
      SM= X*(SEM(I)-SEM(I-1))+SEM(I-1)
      PS= X*(PSCAT(I)-PSCAT(I-1))+PSCAT(I-1)
      SEC= X*(SETOTC(I)-SETOTC(I-1))+SETOTC(I-1)
      SMC= X*(SEMC(I)-SEMC(I-1))+SEMC(I-1)
      PSC= X*(PSCATC(I)-PSCATC(I-1))+PSCATC(I-1)
      return
      end

