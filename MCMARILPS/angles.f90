!     
! File:   angles.f
! Author: Peisi
!
! Created on March 7, 2013, 12:29 PM
!

!--------------------------------------------------------------------
!     Calculates scattering angles when an interaction. 
!     Change the values of position and angle, preserving 
!     Altima(?) arising from iteration.
!     Calculates new direction versor.
!--------------------------------------------------------------------
Subroutine Angles(IMED, A, R, D, U0, U, rand1, rand2)
!    Implicit Real(8) (A - H, O - Z)
    Implicit none
    Real(8) Dimension A(2, 2), R(2, 3), D(2, 3), DSDW(90)
    Real(8) Pi
    Data Pi/3.1415926535897932D0/
    Logical lFlag
    E0 = 10.D0**U0
    E = 10.D0**U
    call AngularDistribution(IMED, E0, E, DSDW)

    STTT = 0.D0
    S1 = 0.D0
    I = 1
    deltaT = 2.D0
    do T = 2.D0, 180.D0, deltaT
        theta = T/180.D0 * Pi
        S2 = DSDW(I) * E
        S2 = S2 * 2.D0 * Pi * DSIN(theta) * Pi/180.D0
        STTT = .5D0 * (S1 + S2) * deltaT + STTT
        I = I + 1
        S1 = S2
    end do
    if (STTT <= 0.D0) then
        return
    endif

    P = 0.D0
    S1 = 0.D0
    I = 1
    lFlag = .true.
    do T = 2.D0, 180.D0, deltaT
        theta = T/180.D0 * Pi
        S2 = DSDW(I) * E
        S2 = S2 * 2.D0 * Pi * DSIN(theta) * Pi/180.D0/STTT
        DS = .5D0 * (S1 + S2) * deltaT
        !/*
        ! *If the Value of expr is: 	Control Transfers To: 
        ! *Less than 0 	Statement label1 
        ! *Equal to 0 	Statement label2 
        ! *Greater than 0 	Statement label3 
        ! */
        if ((P + DS - rand1) < 0) then
            2 P = P + DS
            S1 = S2
            I = I + 1
        else
            3 AA = (S2 - S1)/deltaT

            !     if AA = 0 we must avoid division by zero
            if (aa .eq. 0.d0) then
            th = rand1 * deltaT + t - deltaT
            goto 4
            end if
            B = (S1 * T - S2 * (T - deltaT))/deltaT
            BA = B/AA
            RAD = BA * BA + 2.D0 * (rand1 - P)/AA + (T - deltaT)**2.D0 + 2 * BA * (T - deltaT)
            TH = -BA + DSQRT(RAD)
            if (.not.(TH > T - deltaT .and. TH < T)) TH = -BA - DSQRT(RAD)
            lFlag = .false.
            exit
        endif
    end do
    if (lFlag) TH = 180.D0
    A(2, 1) = TH
    A(2, 2) = 360.D0 * rand2

    do, J = 1, 3
        R(1, J) = R(2, J)
    end do
    A(1, 2) = 0.D0
    if (DABS(D(2, 3)) < 1.D0) then
        A(1, 2) = DACOS(D(2, 1)/DSQRT(1.D0 - D(2, 3) * D(2, 3))) * 180.D0/Pi
    end if
    A(1, 1) = DACOS(D(2, 3)) * 180.D0/Pi

    T1 = A(1, 1)/180.D0 * Pi
    T2 = A(2, 1)/180.D0 * Pi
    F1 = A(1, 2)/180.D0 * Pi
    F2 = A(2, 2)/180.D0 * Pi
    D(2, 1) = DCOS(T1) * DCOS(F1) * DSIN(T2) * DCOS(F2) - DSIN(F1) * DSIN(T2) * DSIN(F2) + DSIN(T1) * DCOS(F1) * DCOS(T2)
    D(2, 2) = DCOS(T1) * DSIN(F1) * DSIN(T2) * DCOS(F2) + DCOS(F1) * DSIN(T2) * DSIN(F2) + DSIN(T1) * DSIN(F1) * DCOS(T2)
    D(2, 3) = DCOS(T1) * DCOS(T2) - DSIN(T1) * DSIN(T2) * DCOS(F2)

    return
end