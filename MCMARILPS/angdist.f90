!     
! File:   angdist.f
! Author: Peisi
!
! Created on March 7, 2013, 12:37 PM
!

!--------------------------------------------------------------------
!     Angular distribution generated from experimental Data. 
!     Do not use the synthetic model.
!--------------------------------------------------------------------
Subroutine AngularDistribution(IMED, E0, E, DSDW)
    Implicit Real(8) (A - H, O - Z)
    !Implicit none
    Parameter (NQEXP = 250, NWEXP = 550, NQEXPC = 250, NWEXPC = 550)
    Common /EXPCOORD/QEXP(NQEXP), HWEXP(NWEXP)
    Common /EXPCOORDC/QEXPC(NQEXPC), HWEXPC(NWEXPC)
    Common /EXPERSAM/DATSAM(NQEXP, NWEXP)
    Common /EXPERCAN/DATCAN(NQEXPC, NWEXPC)
    Dimension DSDW(90)
    Data PI/3.1415926535897932D0/, HQS2M/2.07219093D - 03/
    HW = E0 - E
    I = 1
    do T = 2.D0, 180.D0, 2.D0
        theta = T/180.D0 * PI
        Q = DSQRT((E0 + E - 2.D0 * DSQRT(E0 * E) * DCOS(theta))/HQS2M)
        if (IMED == 0) then
            call Locate(QEXP, NQEXP, Q, JQ)
            call Locate(HWEXP, NWEXP, HW, JW)
            DSDW(I) = 0.D0
            if (JQ >= 1 .and. JQ <= NQEXP - 1 .and. JW >= 1 .and. JW <= NWEXP - 1) then
                call IntSqw(Q, QEXP(JQ), QEXP(JQ + 1), DATSAM(JQ, JW), DATSAM(JQ + 1, JW), S1)
                call IntSqw(Q, QEXP(JQ), QEXP(JQ + 1), DATSAM(JQ, JW + 1), DATSAM(JQ + 1, JW + 1), S2)
                call IntSqw(HW, HWEXP(JW), HWEXP(JW + 1), S1, S2, DSDW(I))
            end if
        end if

        if (IMED == 1) then
            call Locate(QEXPC, NQEXPC, Q, JQ)
            call Locate(HWEXPC, NWEXPC, HW, JW)
            DSDW(I) = 0.D0
            if (JQ >= 1 .and. JQ <= NQEXPC - 1 .and. JW >= 1 .and. JW <= NWEXPC - 1) then
                call IntSqw(Q, QEXPC(JQ), QEXPC(JQ + 1), DATCAN(JQ, JW), DATCAN(JQ + 1, JW), S1)
                call IntSqw(Q, QEXPC(JQ), QEXPC(JQ + 1), DATCAN(JQ, JW + 1), DATCAN(JQ + 1, JW + 1), S2)
                call IntSqw(HW, HWEXPC(JW), HWEXPC(JW + 1), S1, S2, DSDW(I))
            end if
        end if
        !   DSDW(I)=DSDW(I)*DSIN(theta)/Q     <----- Was abolished, NO VA
        !   Multiply by NO SIN (theta), because that is done in the 
        !   ROUTINE angles. As this program has already converted the experimental Data to 
        !   probability P(E0,E,theta) is not required Jacobian in Q.

        I = I + 1
    end do
    return
end

