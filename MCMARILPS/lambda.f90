!     
! File:   lambda.f
! Author: Peisi
!
! Created on March 7, 2013, 1:09 PM
!

!--------------------------------------------------------------------
!     Calculate the mean free path of the neutron and the distance to 
!     the next interaction. Modify the weight of the neutron.
!--------------------------------------------------------------------
Subroutine LAMBDA(XL, randol, WGT, ICAN, IMED, TRANS)
    Implicit Real(8) (A - H, O - Z)
    Common /DATOSSAM/SE, SE0, SM, PS, UTAB(68), SETOT(68), SEM(68), PSCAT(68)
    Common /DATOSCAN/SEC, SE0C, SMC, PSC, SETOTC(68), SEMC(68), PSCATC(68)
    Common /DISTS/DISCE1, DISAG1, DISCI1, DISVA, DISCI2, DISAG2, DISCE2


    TRANS = DEXP(-SMC * (DISCE1 + DISCE2) - SM * (DISAG1 + DISAG2)- SMC * (DISCI1 + DISCI2))
    if (TRANS == 1.D0) then
        WGT = 0.D0
        return
    end if

    A = 1.D0/(1.D0 - TRANS)
    X1 = (1.D0 - DEXP(-SMC * DISCE1)) * A
    X2 = (1.D0 - DEXP(-SMC * DISCE1 - SM * DISAG1)) * A
    X3 = (1.D0 - DEXP(-SMC * DISCE1 - SM * DISAG1 - SMC * DISCI1)) * A
    X4 = (1.D0 - DEXP(-SMC * DISCE1 - SM * DISAG1 - SMC * (DISCI1 + DISCI2))) * A
    X5 = (1.D0 - DEXP(-SMC * DISCE1 - SM * (DISAG1 + DISAG2)- SMC * (DISCI1 + DISCI2))) * A
    if (randol  <  X1) then
        XL = -DLOG(1.D0 - randol/A)/SMC
        ICAN = 1
        IMED = 1 ! SCATT. IN CAN EXT 1
    end if
    if (randol >= X1 .and. randol  <  X2) then
        XL = DISCE1 - (DLOG(1.D0 - randol/A) + SMC * DISCE1)/SM
        IMED = 0 ! SCATT. IN SAMPLE 1
    end if
    if (randol >= X2 .and. randol  <  X3) then
        XL = DISCE1 + DISAG1 - (DLOG(1.D0-(randol/A))+SMC*DISCE1+SM*DISAG1)/SMC
        ICAN = 1
        IMED = 1 ! SCATT. IN CAN INT 1
    end if
    if (randol >= X3 .and. randol  <  X4) then
        XL = DISCE1 + DISAG1 + DISCI1 - (DLOG(1.D0-(randol/A))+SMC*DISCE1+SM*DISAG1+SMC*DISCI1)/SMC
        ICAN = 1
        IMED = 1 ! SCATT. IN CAN INT 2
    end if
    if (randol >= X4 .and. randol  <  X5) then
        XL = DISCE1 + DISAG1 + DISCI1 + DISCI2 - &
            (DLOG(1.D0-(randol/A))+SMC*DISCE1+SM*DISAG1+ &
            SMC*(DISCI1+DISCI2))/SM
        IMED = 0 ! SCATT. IN SAMPLE 2
    end if
    if (randol >= X5 .and. randol  <  1.D0) then
        XL = DISCE1 + DISAG1 + DISCI1 + DISCI2 + DISAG2 - &
            (DLOG(1.D0-(randol/A))+SMC*DISCE1+SM*(DISAG1+DISAG2)+ &
            SMC*(DISCI1+DISCI2))/SMC
        ICAN = 1
        IMED = 1 ! SCATT. IN CAN EXT 2
    end if

    WGT = WGT * (1.D0 - TRANS)

    if (IMED == 0) then
        WGT = WGT * PS
    else
        WGT = WGT * PSC
    end if

    return
end
