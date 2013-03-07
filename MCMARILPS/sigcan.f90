!     
! File:   sigcan.f
! Author: Peisi
!
! Created on March 7, 2013, 12:33 PM
!

!--------------------------------------------------------------------
!     Same as above, for the can
!--------------------------------------------------------------------
Real(8) Function SIGCAN(ifLAG, N, NQ, X)
    Implicit Real(8) (A - H, O - Z)
    Implicit Integer(4) (I - N)
    Parameter (NQEXPC = 250, NWEXPC = 550)
    Common /DATOSN/EINC, U, theta
    Common /DISTS/DISCE1, DISAG1, DISCI1, DISVA, DISCI2, DISAG2, DISCE2
    Common /DATOSSAM/SE, SE0, SM, PS, UTAB(68), SETOT(68), SEM(68), PSCAT(68)
    Common /DATOSCAN/SEC, SE0C, SMC, PSC, SETOTC(68), SEMC(68), PSCATC(68)
    Common /EXPCOORDC/QEXPC(NQEXPC), HWEXPC(NWEXPC)
    Common /EXPERCAN/DATCAN(NQEXPC, NWEXPC)
    Data HQS2M/2.07219093D-03/
    EACT = 10**U
    E = 10.D0**X
    call VALINTERP(U)
    SEACTC = SEC
    call VALINTERP(X)
    HW = EACT - E
    ! Calculate the real Q for this final energy and angle
    if (N == 1) then
        Q = QEXPC(NQ)
    else
        Q = DSQRT((E + EACT - 2.D0 * DSQRT(E * EACT) * DCOS(theta))/HQS2M)
    end if

    call Locate(QEXPC, NQEXPC, Q, JQ)
    call Locate(HWEXPC, NWEXPC, HW, JW)
    S2Dif = 0.D0

    if (JQ >= 1 .and. JQ <= NQEXPC - 1 .and. JW >= 1 .and. JW <= NWEXPC - 1) then
        call IntSqw(Q, QEXPC(JQ), QEXPC(JQ + 1), DATCAN(JQ, JW), DATCAN(JQ + 1, JW), S1)
        call IntSqw(Q, QEXPC(JQ), QEXPC(JQ + 1), DATCAN(JQ, JW + 1), DATCAN(JQ + 1, JW + 1), S2)
        call IntSqw(HW, HWEXPC(JW), HWEXPC(JW + 1), S1, S2, S2Dif)
    end if
    S2Dif = S2Dif * (SE0C/SEACTC) * DSQRT(EINC/EACT)
    if (ifLAG == 0) then
        SIGCAN = S2Dif
        return
    end if

    !  Eye, here we need attenuated neutrons can, but not in sample

    TRANS = DEXP(-SMC * (DISCE1 + DISCE2) - SM * (DISAG1 + DISAG2)- SMC * (DISCI1 + DISCI2))
    SIGCAN = S2Dif * TRANS * EFI(X)

    return
end

