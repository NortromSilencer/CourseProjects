!     
! File:   etkernel.f
! Author: Peisi
!
! Created on March 7, 2013, 12:38 PM
!

!--------------------------------------------------------------------
!     Kernel-based energies transferring experimental Data.
!     EYE: when energy is larger than that of experiment
!     up-scattering for possible processes, missing Data in the
!     zone between the parabola of higher energy and the energy
!     original.
!--------------------------------------------------------------------
  Subroutine EtKernel(IMED, E0, endIST, endISTC)
    Implicit Real(8) (A - H, O - Z)
    Parameter (NQEXP = 250, NWEXP = 550, NQEXPC = 250, NWEXPC = 550)
    Common /EXPCOORD/QEXP(NQEXP), HWEXP(NWEXP)
    Common /EXPCOORDC/QEXPC(NQEXPC), HWEXPC(NWEXPC)
    Common /EXPERSAM/DATSAM(NQEXP, NWEXP)
    Common /EXPERCAN/DATCAN(NQEXPC, NWEXPC)
    Dimension endIST(NWEXP), endISTC(NWEXPC)
    Data HQS2M/2.07219093D-03/

    XK0 = DSQRT(E0/HQS2M)
    ! Fix the limits of integration in Q
    if (IMED == 0) then
        N = NWEXP
    else
        N = NWEXPC
    endif
    loop1:do JW = 1, N
        if (IMED == 0) then
            endIST(JW) = 0.D0
            P = XK0**2 - HWEXP(JW)/HQS2M
            H = HWEXP(JW)
        else
            endISTC(JW) = 0.D0
            P = XK0**2 - HWEXPC(JW)/HQS2M
            H = HWEXPC(JW)
        endif
        if (P < 0.D0) cycle loop1
        R = DSQRT(P)
        if (H < 0.D0) then
            QINF = -XK0 + R
            QSUP = XK0 + R
        else
            QINF = XK0 - R
            QSUP = XK0 + R
        end if
        ! Integral in Q to determine the kernel. The integrand is Q * S (Q, w)
        if (IMED == 0) then
            do JQ = 1, NQEXP
                if (QEXP(JQ) >= QINF .and. QEXP(JQ) < QSUP .and. E0 > HWEXP(JW)) then
                    endIST(JW) = endIST(JW) + QEXP(JQ) * DATSAM(JQ, JW)/DSQRT(E0 * (E0 - HWEXP(JW)))
                end if
            end do
        end if
        if (IMED == 1) then
            do JQ = 1, NQEXPC
                if (QEXPC(JQ) >= QINF .and. QEXPC(JQ) < QSUP .and. E0 > HWEXPC(JW)) then
                    endISTC(JW) = endISTC(JW) + QEXPC(JQ) * DATCAN(JQ, JW)/DSQRT(E0 * (E0 - HWEXPC(JW)))
                end if
            end do
        end if
        end do loop1
    return
  end

