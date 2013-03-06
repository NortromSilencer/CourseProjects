! REVISED VERSION 7 CALCULATED DISTANCES, AS A RESULT OF THE INTERSECTION WITH EACH CAN
! ESTIMATIED 4 DISTANCES BEFORE BECAUSE ONLY HAD MUCH IN MEDIA SPANNING BUT
! NO DETAILED THAT PORTION THEREOF. This is important because we WANT TO KNOW the actual distances 
! AND REAL TIME DELAY IN THE SAMPLE.
! Corrected error in SIGMUE And SIGCAN. It was wrong computation of Q, using a variable
! E0 was always 0, and had to use EACT
!   Multiple Scattering Argonne experiment 

!   Relationship spectra are recorded:
!   ESPS(I,J,K):     Spectrum accumulating counts: I th Energy, 
!                    J th scattering order and K th Q.
!                    The ESPS (I, 10, K) contains the scattering orders 10 or 
!                    greater, and ESPS (I, 11, K), the sum of K = 2 ... 10 (the 
!                    multiple total).
!   EPSPS2(I,J,K):   Accumulates squared elements of the previous
!   RESPS(I,J,K):    Mean value of ESPS.
!   REPSPS2(I,J,K):  Error on ESPS.
!
!   ESP1C(I,J):      Single scattering in the can; I th Energy and J th Q
!   EPSP1C2(I,J,K):  Accumulates squared elements of the previous
!   RESP1C(I,J):     Mean Value of ESP1C
!   REPSP1C2(I,J,K): Error on ESP1C.
!
!   ESP1CNAS(I,J):   Single scattering in the can; I th Energy and J th Q without attenuation.
!   ESP1CNAS2(I,J):  Accumulates squared elements of the previous
!   RESP1CNAS(I,J):  Mean Value of ESP1CNAS.
!   RESP1CNAS2(I,J): Error on ESP1CNAS.
!
!   RESPT(I,J):      Total Scattering with I th Energy and J th Q.
!   REPSPT2(I,J,K):  Accumulates squared elements of the previous
!   ESP1NAS(I,J):    Single Scattering without attenuation
!   RESP1NAS(I,J):   Mean value of the previous.
!----------------------------------------------------------------------
      PROGRAM MCMARI
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER 
     *(NQ=250,NW=550,NQEXP=250,NWEXP=550,NQEXPC=250,NWEXPC=550)
      REAL*8 KBT
      CHARACTER*11 title*45,FILMUE*45,FILCAN*45
      DIMENSION POS(2,3),D(2,3),A(2,2),RANDOL(1000),
     *ESPS(NW,11,NQ),ESPS2(NW,11,NQ),RESPS(NW,11,NQ),RESPS2(NW,11,NQ),
     *ESP1C(NW,NQ),ESP1C2(NW,NQ),RESP1C(NW,NQ),RESP1C2(NW,NQ),
     *ESP1CNAS(NW,NQ),ESP1CNAS2(NW,NQ),RESP1CNAS(NW,NQ),
     *RESP1CNAS2(NW,NQ),
     *RESPT(NW,NQ),RESPT2(NW,NQ),ESP1NAS(NW,NQ),RESP1NAS(NW,NQ),
     *NEVENT(10,2),QSAL(NQ),HWANALIZ(NW)
      COMMON /DATOSN/E0,U,THETA
      COMMON /DISTS/DISCE1,DISAG1,DISCI1,DISVA,DISCI2,DISAG2,DISCE2
      COMMON /DATOSSAM/SE,SE0,SM,PS,UTAB(68),SETOT(68),SEM(68),PSCAT(68)
      COMMON /DATOSCAN/SEC,SE0C,SMC,PSC,SETOTC(68),SEMC(68),PSCATC(68)
      COMMON /EXPCOORD/QEXP(NQEXP),HWEXP(NWEXP)
      COMMON /EXPCOORDC/QEXPC(NQEXPC),HWEXPC(NWEXPC)
      COMMON /EXPERSAM/DATSAM(NQEXP,NWEXP)
      COMMON /EXPERCAN/DATCAN(NQEXPC,NWEXPC)
      COMMON /MUESTRA/KBT
      COMMON /NEUTRON/VDIR(3),VPOS(3)
      COMMON /GEOMUEST/RC1,HC1,RC2,HC2,RC3,HC3,RC4,HC4,
     +C01(3),C02(3),C03(3),C04(3)
      DATA PI/3.1415926535897932D0/,HQS2M/2.07219093D-03/
!----------------------------------------------------------------------
!     Seed generator For Lahey Fortran Compiler (Lahey Computer System Inc.)
!11    DSEED=RRAND()*10000.
!--------------------------------
!     For Microsoft Fortran (Microsoft Fortran PowerStation, Microsoft)
!11    CALL GETTIM(I,J,K,L)
!      DSEED=L+K*10000./60.
!     For Linux on DEC Alpha or Intel x86/x64 PC
11    CALL RANDOM_SEED
      CALL RANDOM_NUMBER(DSEED)
      DSEED=DSEED*10000
!     For Linux with g77 (GNU Fortran compiler which is part of gcc)

!11    j=time()
!      call srand(j)
!      dseed=10000*rand()

      IF (DSEED.LE.0.D0) GOTO 11
      WRITE(6,*)'Initial Seed=',DSEED
!----------------------------------------------------------------------
! Several calls to ggubfs to get rid of the initial seed (for Intel Fortran compiler)

      call ggubfs(dseed,randol)
      call ggubfs(dseed,randol)
      call ggubfs(dseed,randol)
      call ggubfs(dseed,randol)
!----------------------------------------------------------------------
! Initializes Accumulates and Output spectra
      DO 10 J=1,NW
        DO 10 L=1,NQ
           ESP1C(J,L)=0.D0
           ESP1C2(J,L)=0.D0
           RESP1C(J,L)=0.D0
           RESP1C2(J,L)=0.D0
           ESP1CNAS(J,L)=0.D0
           ESP1CNAS2(J,L)=0.D0
           RESP1CNAS(J,L)=0.D0
           RESP1CNAS2(J,L)=0.D0
           RESPT(J,L)=0.D0
           RESPT2(J,L)=0.D0
           ESP1NAS(J,L)=0.D0
           RESP1NAS(J,L)=0.D0
           DO 10 K=1,11
             ESPS(J,K,L)=0.D0
             ESPS2(J,K,L)=0.D0
             RESPS(J,K,L)=0.D0
10           RESPS2(J,K,L)=0.D0
!----------------------------------------------------------------------
! GR  : Initializes NEVENT(j,k), i.e.the number of neutrons scattered j times in the sample (k=1) or in the can (k=2)
! LPS : Initializes NEVENT vector: number of events depending on the number of scattering J and type of neutron (from can or sample)             

      DO 12 J=1,10
        DO 12 K=1,2
12        NEVENT(J,K)=0
      NUMCIC0=1
      NTOT=0
      NOUT=0
      NHIS=0
      TPRO=0.D0
!----------------------------------------------------------------------
! GR  : Reading experimental total effective cross sections of shield and can (micro and macroscopic) and scattering probabilities (front and abs).
! LPS : Reading total cross-sections and can sample, microscopic and macroscopic, and scattering probabilities (vs. absorption).
! LPS : The total cross-sections are used experimental.

      OPEN(1,FILE='mcsminp.dat',STATUS='OLD')
      READ(1,*)
      READ(1,*)
      DO 2 I=1,68
2     READ(1,*)UTAB(I),SETOT(I),SEM(I),PSCAT(I),SETOTC(I),SEMC(I),
     #PSCATC(I)
      CLOSE(1)
!----------------------------------------------------------------------
! Reading Simulation Parameters
      CALL PARAMETROS(title,FILMUE,FILCAN,NCIC,NT,AH,WH,E0,XLS,WCO,
     +ICONT)
      CALL TIME (HORA1)
      HORA2=HORA1      

      
      
      
      CALL PANTALLA(title,RC1,RC3,AH,WH,XLS,FILMUE,NT*NCIC,NTOT,0.D0)

      U0=DLOG10(E0)
!.... reading experimental data for sample

      
      OPEN(1,FILE=FILMUE)
      READ(1,*)
      READ(1,*)
      READ(1,'(10(1PE12.5))')QEXP
      READ(1,*)
      READ(1,'(10(1PE12.5))')HWEXP
      READ(1,*)
      DO I=1,NQEXP
        READ(1,'(10(1PE12.5))')(DATSAM(I,J),J=1,NWEXP)
      END DO
      CLOSE(1)

!... and converting energy in eV
      DO I=1,NWEXP
          HWEXP(I)=HWEXP(I)/1000.D0
      END DO
c... reading experimental data for can
      OPEN(1,FILE=FILCAN)
      READ(1,*)
      READ(1,*)
      READ(1,'(10(1PE12.5))')QEXPC
      READ(1,*)
      READ(1,'(10(1PE12.5))')HWEXPC
      READ(1,*)
      DO I=1,NQEXPC
        READ(1,'(10(1PE12.5))')(DATCAN(I,J),J=1,NWEXPC)
      END DO
      CLOSE(1)

!... and converting energy in eV
      DO I=1,NWEXPC
          HWEXPC(I)=HWEXPC(I)/1000.D0
      END DO

!----------------------------------------------------------------------
! Initialize energy transfer vector (eV) for the analyzer
      DO 13 I=1,NW
13    HWANALIZ(I)=HWEXP(I)

!... and the momentum transfer vector (1/A).

      DO 14 I=1,NQ
14    QSAL(I)=QEXP(I)
!----------------------------------------------------------------------
! Transforms sample and can data to probabilities. alfa is the normalization. 
!... first the sample
      ALFA=0.D0
      DE=HWEXP(2)-HWEXP(1)
      DQ=QEXP(2)-QEXP(1)
      DO 500 I=1,NQEXP
        SUME=0.D0
        DO 501 J=1,NWEXP
501       SUME=SUME+DATSAM(I,J)*DE
        ALFA=ALFA+SUME*QEXP(I)*DQ
500   CONTINUE
      WRITE(6,*)'ALFA SAMPLE= ',ALFA

      DO 201 I=1,NQEXP
      DO 201 J=1,NWEXP
      DATSAM(I,J)=DATSAM(I,J)*(DSQRT(E0*(E0-HWEXP(J)))/HQS2M)
     #       /2.D0/PI/ALFA
201   CONTINUE
!...
      ALFA=0.D0
      DE=HWEXP(2)-HWEXP(1)
      DQ=QEXP(2)-QEXP(1)
      DO 900 I=1,NQEXP
        SUME=0.D0
        DO 901 K=1,NWEXP
901       SUME=SUME+2.D0*PI/(DSQRT(E0*(E0-HWEXP(K)))/HQS2M)*
     #         DATSAM(I,K)*DE
        ALFA=ALFA+SUME*QEXP(I)*DQ
900   CONTINUE
      WRITE(6,*)'ALFA SAMPLE= ',ALFA
!... then the can
      ALFA=0.D0
      DE=HWEXPC(2)-HWEXPC(1)
      DQ=QEXPC(2)-QEXPC(1)
      DO 510 I=1,NQEXPC
        SUME=0.D0
        DO 511 J=1,NWEXPC
511       SUME=SUME+DATCAN(I,J)*DE
        ALFA=ALFA+SUME*QEXP(I)*DQ
510   CONTINUE
      WRITE(6,*)'ALFA CAN= ',ALFA

      DO 211 I=1,NQEXPC
      DO 211 J=1,NWEXPC
      DATCAN(I,J)=DATCAN(I,J)*(DSQRT(E0*(E0-HWEXPC(J)))/HQS2M)
     #       /2.D0/PI/ALFA
211   CONTINUE
!...
      ALFA=0.D0
      DE=HWEXPC(2)-HWEXPC(1)
      DQ=QEXPC(2)-QEXPC(1)
      DO 910 I=1,NQEXPC
        SUME=0.D0
        DO 911 K=1,NWEXPC
911       SUME=SUME+2.D0*PI/(DSQRT(E0*(E0-HWEXPC(K)))/HQS2M)*
     #         DATCAN(I,K)*DE
        ALFA=ALFA+SUME*QEXPC(I)*DQ
910   CONTINUE
      WRITE(6,*)'ALFA CAN= ',ALFA

!----------------------------------------------------------------------
! reads the previous RUN if this is the continuation

       IF (ICONT.EQ.1) THEN
         OPEN(1,FILE='CONT.DAT')
         READ (1,*) NUMCIC0
         READ (1,*) NTOT
         READ (1,*) NOUT
         READ (1,*) NHIS
         READ (1,*) TPRO
         READ(1,'(10I7)')NEVENT
         DO 64 I=1,NQ
         DO 61 J=1,NW
61       READ(1,999)ESP1C(J,I),ESP1C2(J,I),ESP1NAS(J,I)
         DO 62 J=1,NW
62       READ(1,999)(ESPS(J,K,I),K=1,11)
         DO 63 J=1,NW
63       READ(1,999)(ESPS2(J,K,I),K=1,11)
64       CONTINUE
         CLOSE(1)
       END IF
!----------------------------------------------------------------------
! Initialize random number vector and call the random number generator routine

      KRAND=0
      CALL GGUBFS(DSEED,RANDOL)
!----------------------------------------------------------------------
! Starts the great loop of NCIC cycles and NT neutrons
c Generally NCIC=1 and the various RUNs are controlled by BATCHITO(.sh)
      DO 1000 NUMCIC=NUMCIC0,NCIC
      DO 200 N=1,NT
         NTOT=NTOT+1
         KRAND=KRAND+1
         IF (KRAND.GT.1000) THEN          ! IF JUST THE RANDOM GENERATED
          CALL GGUBFS(DSEED,RANDOL)       ! NEW
          KRAND=1
         END IF
         NS=1                             ! NUMBER neutron scattering
         U=U0                             ! INITIAL STUPOR
         WGT=1.D0                         ! INITIAL WEIGHT
         DELTE=0.D0                       ! CORR. AL HW effective
                                          ! DELAY DUE TO SCATT. MULT.
         ICAN=0

!... parameter interpolation for the incident energy 

         CALL VALINTERP(U)
         SE0=SE
         SE0C=SEC
         IRAND1=MOD(KRAND,1000)+1
         IRAND=MOD(IRAND1,1000)+1
C.. Neutron Entry position
!... I think that D is the direction, generally taken along Z

         CALL RCERO(RC1,AH,WH,POS,A,RANDOL(IRAND1),RANDOL(IRAND))
         D(1,1)=DSIN(A(2,1)/180.D0*PI)*DCOS(A(2,2)/180.D0*PI)
         D(1,2)=DSIN(A(2,1)/180.D0*PI)*DSIN(A(2,2)/180.D0*PI)
         D(1,3)=DCOS(A(2,1)/180.D0*PI)
         DO 20,I=1,3
20       D(2,I)=D(1,I)

!... Multiple scattering loop begins

100      CONTINUE
            IRAND=MOD(IRAND,1000)+1
! Initializes vector direction and position of the neutron (auxiliary)
            VDIR(1)=D(2,1)
            VDIR(2)=D(2,2)
            VDIR(3)=D(2,3)
            VPOS(1)=POS(1,1)
            VPOS(2)=POS(1,2)
            VPOS(3)=POS(1,3)
! GR  : calculate flight distances  in the various Media considering actual position and direction
! LPS : Flight distances calculated in different ways according to position and heading            
            CALL DISTAN(1)
! GR :      gives(?) flight distance and changes, lowering it, the weight
! LPS:      Lots flight distance and updates the weight
            CALL LAMBDA(XL,RANDOL(IRAND),WGT,ICAN,IMED,TRANS)
            IF (NS.EQ.1) THEN
                NHIS=NHIS+1
                TPRO=TPRO+TRANS
            END IF
            IRAND=MOD(IRAND,1000)+1
!... computes the ultimate energy correction due to the delay in the sample
            IF (NS.GT.1) DELTE=DELTE+XL/XLS/DSQRT(10.D0**U)
!... call "rulrus" (Russian Roulette), if the neutron dies, a new one is created
            CALL RULRUS(RANDOL(IRAND),WGT,WCO,IFLAG)
            IF (IFLAG.EQ.0) GOTO 200
!... update position of the neutron
            CALL POSNUE(POS,XL,D)
!.... calculate the contribution for each detector
! 1) Calculate DISTAN POS(2,...)
! 2) Calculate angle theta
! 3) Calculate the contribution 1 / s(E0) d2s/dOdE (Q, w) exp[-s(E) D] 
!    in the energies of the analyzer
! 4) sum for each detector
         VPOS(1)=POS(2,1)
         VPOS(2)=POS(2,2)
         VPOS(3)=POS(2,3)
         VDIR(1)=0.D0                     ! Detectors in the yz plane
         LL=NS
         IF (NS.GT.10) LL=10
 
!... Case neutron in the sample  -------------------------------------------
      IF (IMED.EQ.0) THEN                                                  |
         IF (ICAN.EQ.0) THEN                                               |
           IF (NS.LE.10) NEVENT(LL,1)=NEVENT(LL,1)+1                       |
         ELSE                                                              |
           IF (NS.LE.10) NEVENT(LL,2)=NEVENT(LL,2)+1                       |
         END IF                                                            |
!...                                                                       |
         DO 30 K=1,NW                                                      |
! Explore the final energies. Calculates EFIN final energies such that     |
! takes into account the delay in the sample for HWANALIZ (K) prefixed     |
                                                                           |
         EFIN=1.D0/(DELTE-1.D0/DSQRT(E0-HWANALIZ(K)))**2                   |
         UFIN=DLOG10(EFIN)                                                 |
                                                                           |
         DO 30 I=1,NQ                                                      |
         VAL=0.D0                                                          |
! Explore the actual Q channel QSAL(I), and based on that                  |
! delivers the output angle of the energy deduced effective HWANALIZ E0-(K)|
         T=(2*E0-HWANALIZ(K)-HQS2M*QSAL(I)**2)/                            |
     #   (2.D0*DSQRT(E0*(E0-HWANALIZ(K))))                                 |
         IF(DABS(T).GE.1.D0) GOTO 30                                       |
         T=DACOS(T)                                                        |
                                                                           |
         VDIR(2)=DSIN(T)                  ! Directions to the detector     |
         VDIR(3)=DCOS(T)                  ! given angle                    |
         THETA=DACOS(VDIR(2)*D(2,2)+VDIR(3)*D(2,3))  !actual interaction Ang. |
         CALL DISTAN(1)                   ! Distances to travel            |
                                                                           |
         VAL=SIGMUE(1,LL,I,UFIN)                                           |
         ESPS(K,LL,I)=ESPS(K,LL,I)+VAL*WGT                                 |
         ESPS2(K,LL,I)=ESPS2(K,LL,I)+(VAL*WGT)**2.D0                       |
         IF (LL.EQ.1) ESP1NAS(K,I)=ESP1NAS(K,I)+SIGMUE(0,LL,I,UFIN)*WGT    |
30       CONTINUE                                                          |
      END IF                                                               |
c End case of neutron within the sample     -------------------------------|      
! Case of the neutron in can                --------------------------------
                                                                           |
      IF (IMED.EQ.1) THEN                                                  |
! Generates scattering angular distribution for the given energy in the can|
         IF (NS.LE.10) NEVENT(LL,2)=NEVENT(LL,2)+1                         |
         DO 29 K=1,NW                                                      |
! Explore the final energies. Calculates EFIN final energies such that     |
! takes into account the delay in the sample for HWANALIZ (I) prefixed     |                                                         |
!
         EFIN=1.D0/(DELTE-1.D0/DSQRT(E0-HWANALIZ(K)))**2                   |
         UFIN=DLOG10(EFIN)                                                 |
                                                                           |
         DO 29 I=1,NQ                                                      |
         VAL=0.D0                                                          |
! Explore the Q channel QSAL cash given by (I), and based on that takes 
! the output angle deduced from the effective energy-HWANALIZ E0 (K).
         T=(2*E0-HWANALIZ(K)-HQS2M*QSAL(I)**2)/                            |
     #   (2.D0*DSQRT(E0*(E0-HWANALIZ(K))))                                 |
         IF(DABS(T).GE.1.D0) GOTO 29                                       |
         T=DACOS(T)                                                        |
         VDIR(2)=DSIN(T)                  ! Directions to the detector     |
         VDIR(3)=DCOS(T)                  ! given angle                    |
         THETA=DACOS(VDIR(2)*D(2,2)+VDIR(3)*D(2,3))                        |
                                                                           |
         CALL DISTAN(1)                   ! Distances to travel            |
                                                                           |
         VAL=SIGCAN(1,LL,I,UFIN)                                           |
         VAL0=SIGCAN(0,LL,I,UFIN)                                          |
         IF (LL.EQ.1) THEN                                                 |
           ESP1C(K,I)=ESP1C(K,I)+VAL*WGT                                   |
           ESP1C2(K,I)=ESP1C2(K,I)+(VAL*WGT)**2.D0                         |
           ESP1CNAS(K,I)=ESP1CNAS(K,I)+VAL0*WGT                            |
           ESP1CNAS2(K,I)=ESP1CNAS2(K,I)+(VAL0*WGT)**2.D0                  |
         ELSE                                                              |
           ESPS(K,LL,I)=ESPS(K,LL,I)+VAL*WGT                               |
           ESPS2(K,LL,I)=ESPS2(K,LL,I)+(VAL*WGT)**2.D0                     |
         END IF                                                            |
29       CONTINUE                                                          |
      END IF                                                               |
!    End of the neutron case in can-----------------------------------------
!... Start a new scattering analysis

         NS=NS+1
c
!... GR : the neutron has a new energy and direction
!    LPS: Lots new energy and angle         
         IRAND=MOD(IRAND,1000)+1
         U1=U
         CALL ENERGIA(IMED,U,IFLAG,RANDOL(IRAND))
	 IF (IFLAG.EQ.0) THEN
	   NOUT=NOUT+1
	   GOTO 200
	 END IF
         IRAND1=MOD(IRAND,1000)+1
         IRAND=MOD(IRAND1,1000)+1
         CALL ANGULOS(IMED,A,POS,D,U1,U,RANDOL(IRAND1),RANDOL(IRAND))
         CALL VALINTERP(U)
         GOTO 100
!... when the neutron history ends we can consider a new neutron
200    CONTINUE
!... end of the NT loop
!----------------------------------------------------------------------
!... Calculates mean values and errors
      DO 40 K=1,NW
      DO 40 I=1,NQ
         IF (NTOT.GT.0) RESP1C(K,I)=ESP1C(K,I)/DBLE(NTOT)
         IF (NTOT.GT.0) RESP1CNAS(K,I)=ESP1CNAS(K,I)/DBLE(NTOT)
         IF (NTOT.GT.0) RESP1NAS(K,I)=ESP1NAS(K,I)/DBLE(NTOT)
         IF (NTOT.GT.1)
     #RESP1C2(K,I)=DSQRT(DABS(ESP1C2(K,I)/DBLE(NTOT)-(ESP1C(K,I)/
     #DBLE(NTOT))**2.D0)/DBLE(NTOT-1))
         IF (NTOT.GT.1)
     #RESP1CNAS2(K,I)=DSQRT(DABS(ESP1CNAS2(K,I)/DBLE(NTOT)-
     #(ESP1CNAS(K,I)/DBLE(NTOT))**2.D0)/DBLE(NTOT-1))
      DO 40 LL=1,10
         NN=NEVENT(LL,1)+NEVENT(LL,2)
         IF (NTOT.GT.0) RESPS(K,LL,I)=ESPS(K,LL,I)/DBLE(NTOT)
40       IF (NN.GT.1)
     #RESPS2(K,LL,I)=DSQRT(DABS(ESPS2(K,LL,I)/DBLE(NN)-(ESPS(K,LL,I)/
     #DBLE(NN))**2.D0)/DBLE(NN-1))
      RTPRO=TPRO/DBLE(NHIS)
c
!... Sum total multiple scattering and its error

      DO 42 K=1,NW
      DO 42 I=1,NQ
      RESPS(K,11,I)=0.D0
      RESPS2(K,11,I)=0.D0
      DO 41 LL=2,10
      RESPS(K,11,I)=RESPS(K,11,I)+RESPS(K,LL,I)
41    RESPS2(K,11,I)=RESPS2(K,11,I)+RESPS2(K,LL,I)**2
42    RESPS2(K,11,I)=DSQRT(RESPS2(K,11,I))
c
!... Sum the total scattering and its error

      DO 51 K=1,NW
      DO 51 I=1,NQ
      RESPT(K,I)=RESPS(K,1,I)+RESP1C(K,I)+RESPS(K,11,I)
51    RESPT2(K,I)=DSQRT(RESPS2(K,1,I)**2+RESP1C2(K,I)**2+
     #            RESPS2(K,11,I)**2)
      

! Write outputs
!      CALL SYSTEM ('ECHO ON')
!      CALL SYSTEM ('IF EXIST MCSMR00.BAK DEL MCSMR00.BAK')
!      CALL SYSTEM ('RENAME MCSMR00.OUT MCSMR00.BAK')
! Write on the LOG file
      OPEN (1,FILE='MCSMR00.LOG')
      WRITE (1,*)'Total Neutron Number=     ',NTOT
      WRITE (1,*)'NEUTRONES FUERA DE RANGO= ',NOUT
      WRITE (1,*)'TRANSMISI?¡°N PROMEDIO ............ ',RTPRO
      WRITE(1,'(10I8)')NEVENT
      CLOSE(1)
c
      OPEN (1,FILE='MCSMR00.OUT')
      DO 31 I=1,NQ
      WRITE(1,*)' Q= ',QSAL(I)
      WRITE(1,998)
998   FORMAT('    HW          SINGLE       S/ATEN         CAN    
     #MULTIPLE        TOTAL   SING/TOT         ATEN    ATENC')

      DO 31 K=1,NW
      FMUL=0.D0
      AT=0.D0
      ATC=0.D0
      IF(RESPT(K,I).NE.0.D0) FMUL=RESPS(K,1,I)/RESPT(K,I)
      IF(RESP1NAS(K,I).NE.0.D0) AT=RESPS(K,1,I)/RESP1NAS(K,I)
      IF(RESP1CNAS(K,I).NE.0.D0) ATC=RESP1C(K,I)/RESP1CNAS(K,I)
31    WRITE(1,999)HWANALIZ(K),RESPS(K,1,I),RESP1NAS(K,I),RESP1C(K,I),
     #            RESPS(K,11,I),RESPT(K,I),FMUL,AT,ATC
      CLOSE(1)

!    ------------------------------------------------------------------
!      CALL SYSTEM ('IF EXIST MCSME00.BAK DEL MCSME00.BAK')
!      CALL SYSTEM ('RENAME MCSME00.OUT MCSME00.BAK')
!      OPEN (1,FILE='MCSME00.OUT')
!      WRITE(1,*)'ERRORES'
!      WRITE(1,*)
!      DO 32 I=1,NQ
!      WRITE(1,*)' Q= ',QSAL(I)
!      DO 32 K=1,NW
!32    WRITE(1,999)HWANALIZ(K),RESPS2(K,1,I),RESP1C2(K,I),RESPS2(K,11,I),
!     #            RESPT2(K,I)
!      CLOSE(1)

!      OPEN (1,FILE='CONT.DAT')
!      WRITE (1,*) NUMCIC0
!      WRITE (1,*) NTOT
!      WRITE (1,*) NOUT
!      WRITE (1,*) NHIS
!      WRITE (1,*) TPRO
!      WRITE(1,'(10I7)')NEVENT
!      DO 68 I=1,NQ
!      DO 65 J=1,NW
!65    WRITE(1,999)ESP1C(J,I),ESP1C2(J,I),ESP1NAS(J,I)
!      DO 66 J=1,NW
!66    WRITE(1,999)(ESPS(J,K,I),K=1,11)
!      DO 67 J=1,NW
!67    WRITE(1,999)(ESPS2(J,K,I),K=1,11)
!68    CONTINUE
!      CLOSE(1)      
!      ERAT-------------------------------------------------------------
999   FORMAT(1X,12(1PE12.5,1X))

      CALL TIME (HORA1)

      CALL PANTALLA(HORA1,HORA2,TITULO,RC1,RC3,AH,WH,XLS,FILMUE,
     #NT*NCIC,NTOT,RTPRO)

      HORA2=HORA1

1000  CONTINUE
      END
!---------------------------------------------------------------------
! Reading parameters from mcsmpar.dat
!---------------------------------------------------------------------
      SUBROUTINE PARAMETROS(title,FILMUE,FILCAN,NCIC,NT,AH,WH,E0,XLS,
     +WCO,ICONT)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOMUEST/RC1,HC1,RC2,HC2,RC3,HC3,RC4,HC4,
     +C01(3),C02(3),C03(3),C04(3)
      CHARACTER*45 title,FILMUE,FILCAN
      OPEN(1,FILE='mcsmpar.dat',STATUS='OLD')
      READ(1,'(A)')title
      WRITE(6,'(1x,A)')title
      READ(1,'(A)')FILMUE
      WRITE(6,'(1x,A)')FILMUE
      READ(1,'(A)')FILCAN
      WRITE(6,'(1x,A)')FILCAN
      READ(1,*)NCIC
      WRITE(6,'('' Number of cycles'',I4)')NCIC
      READ(1,*)NT
      WRITE(6,'('' Number of neutrons per Cycle '',I6)')NT
!
      READ(1,*)DECE
      WRITE(6,'('' External Can External Diameter '',F6.3)')DECE
      READ(1,*)DICE
      WRITE(6,'('' External Can Internal Diameter '',F6.3)')DICE
      READ(1,*)HCE
      WRITE(6,'('' External Can Height '',F6.3)')HCE
      READ(1,*)ETCE
      WRITE(6,'('' External Can Cover Thickness '',F6.3)')ETCE
!
      READ(1,*)DECI
      WRITE(6,'('' Internal Can External Diameter '',F6.3)')DECI
      READ(1,*)DICI
      WRITE(6,'('' Internal Can Internal Diameter '',F6.3)')DICI
      READ(1,*)HCI
      WRITE(6,'('' Internal Can Height '',F6.3)')HCI
      READ(1,*)ETCI
      WRITE(6,'('' Internal Can Cover Thickness '',F6.3)')ETCI
      
!--------------------------------------------------------------------------      
      
!... Derived measures are measures that define the 4 cylinder materials
! cylinder 1  
      rc1=dece/2.d0      ! cylinder radius 1
      hc1=hce+2.d0*etce  ! cylinder height 1
      c01(ix)=0.d0       ! cylinder coordinates 1
      c01(iy)=0.d0
      c01(iz)=0.d0
! cylinder 2
      rc2=dice/2.d0      ! cylinder radius 2
      hc2=hci+2.d0*etci  ! cylinder height 2
      c02(ix)=0.d0       ! cylinder coordinates 2
      c02(iy)=0.d0
      c02(iz)=0.d0
! cylinder 3
      rc3=deci/2.d0      ! cylinder radius 3
      hc3=hci+2.d0*etci  ! cylinder height 3
      c03(ix)=0.d0       ! cylinder coordinates 3
      c03(iy)=0.d0
      c03(iz)=0.d0
! cylinder 4 ... VOID
      rc4=dici/2.d0      ! cylinder radius 4
      hc4=hci+2.d0*etci  ! cylinder height 4
      c04(ix)=0.d0       ! cylinder coordinates 4
      c04(iy)=0.d0
      c04(iz)=0.d0

!-------------------------------------------

      READ(1,*)AH
      WRITE(6,'('' Beam Height '',F6.3,'' cm'')')AH
      READ(1,*)WH
      WRITE(6,'('' Beam Width '',F6.3,'' cm'')')WH
      READ(1,*)E0
      WRITE(6,'('' Incident Energy '',1PE12.5)')E0
      READ(1,*)XLS
      WRITE(6,'('' Outgoing Flight Length '',1PE12.5)')XLS
      READ(1,*)WCO
      WRITE(6,'('' Cutting Weight '',1PE10.3)')WCO
      READ(1,*)ICONT
      WRITE(6,'('' Continuing '',I3)')ICONT
      CLOSE(1)
!...

      RETURN
      END
      
!--------------------------------------------------------------------
!     GR  : write parameters on visual shell .sh
!     LPS : Screen output using dialog utility      
!--------------------------------------------------------------------
      SUBROUTINE PANTALLA(title,R1,R2,AH,WH,XLS,ARFIL,NPED,NTOT,TPRO)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER* 11 ARFIL*45,title*45

      open(7,file='info.dat')
     WRITE(7,*)'current Time ',HORA1, '  previous drop ',HORA2
      WRITE(7,*)title
      WRITE(7,*)'RC1= ',R1,'    RC3= ',R2
      WRITE(7,*)'AH=  ',AH, '   WH= ',WH
      WRITE(7,*)'XLS= ',XLS
      WRITE(7,*)'Sample file: ', ARFIL
      WRITE(7,*)'Number total order  ',NPED
      WRITE(7,*)'Neutron calendar  ',NTOT
      WRITE(7,*)'Average Transmission ',TPRO
      close(7)
      call system('./infobox4')

      END

!--------------------------------------------------------------------
!    Calculates the position on the anterior sample face, where
!    the neutron enters, and its entry angle.
!--------------------------------------------------------------------
      SUBROUTINE RCERO(RADIO,AH,WH,R,A,RAND1,RAND2)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(2,3),A(2,2)
      A(1,1)=0.D0
      A(1,2)=0.D0
!---- Aqua-position is calculated input
      R(1,1)=AH*(RAND2-.5D0)
      R(1,2)=WH*(RAND1-.5D0)
      R(1,3)=-DSQRT(RADIO**2.D0-R(1,2)**2.D0)
!---- Aqua-placed entry angle
      A(2,1)=0.D0
      A(2,2)=0.D0
      RETURN
      END

!----------------------------------------------------------------------
!  calculate intersection distances with limits of the 
!  different zones
!
!   +--------------------------------------+
!   |               CAN (1)                |
!   | +----+-+--------------------+-+----+ |
!   | |    | |                    | |    | |
!   | |    | |                    | |    | |
!   | |    | |                    | |    | |
!   | |    | |                    | |    | |
!   | |    | |                    | |    | |
!   | |    | |                    | |    | |
!   | |    | |                    | |    | |  ^ X
!   | |    | |                    | |    | |  |
!   | | W  | |                    | |  W | |  |
!   | | A  |C|                    |C|  A | |  |
!   | | T  |A|                    |A|  T | |  |
!   | | E  |N|       VOID (4)     |N|  E | |  |
!   | | R  | |                    | |  R | |  +------------> Z
!   | |    |3|                    |3|    | |
!   | | 2  | |                    | |  2 | |
!   | |    | |                    | |    | |
!   | |    | |                    | |    | |
!   | |    | |                    | |    | |
!   | |    | |                    | |    | |
!   | |    | |                    | |    | |
!   | |    | |                    | |    | |
!   | |    | |                    | |    | |
!   | |    | |                    | |    | |
!   | |    | |                    | |    | |
!   | +----+-+--------------------+-+----+ |
!   |               CAN (1)                |
!   +--------------------------------------+



      subroutine distan(iaxis)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOMUEST/RC1,HC1,RC2,HC2,RC3,HC3,RC4,HC4,
     +C01(3),C02(3),C03(3),C04(3)
      COMMON /NEUTRON/D(3),R0(3)
      COMMON /DISTS/DISCE1,DISAG1,DISCI1,DISVA,DISCI2,DISAG2,DISCE2

      ix=iaxis
      iy=modulo3(iaxis+1)
      iz=modulo3(iaxis+2)

!  Basic measurements in centimeters
!      hce=4.94d0        ! Total length external light water can
!      dece=4.5d0        ! Outer diameter of outer can
!      dice=4.4d0        ! Inner diameter of outer can
!      deci=4.3d0        ! Outer diameter of inner can
!      dici=4.2d0        ! Inner diameter of inner can
!      etce=1.0 d0       ! Thickness can lids external
!      etci=0.0 d0       ! Thickness can tops internal


      s=dsqrt(d(1)**2+d(2)**2+d(3)**2)
      d(1)=d(1)/s
      d(2)=d(2)/s
      d(3)=d(3)/s

! Calculates distances cylinder can

      call discil(iaxis,r0,d,c01,rc1,hc1,x11,x21,icompl)
      call discil(iaxis,r0,d,c02,rc2,hc2,x12,x22,icompl)
      call discil(iaxis,r0,d,c03,rc3,hc3,x13,x23,icompl)
      call discil(iaxis,r0,d,c04,rc4,hc4,x14,x24,icompl)
      
! When you enter the neutron from the outside, may be slightly "out of sample", 
! AND THE VALUE OF INTERSECTION WITH OUTSIDE RIGHT CAN
! is X2E , and not X1E.

      if(x21.gt.0.d0.and.x21.gt.x1e) x1e=x2e

      dc1= dabs(x21-x11)
      dc2= dabs(x22-x12)
      dc3= dabs(x23-x13)
      dc4= dabs(x24-x14)

      disce1=x12-x11
      disag1=x13-x12
      disci1=x14-x13
      disva =x24-x14
      disci2=x23-x24
      disag2=x22-x23
      disce2=x21-x22



      end

!-------------------------------------------------------------------------
! Given a cylinder, a point, a line and a flight direction,
! calculated intersection points of line and cylinder
!                                  x
!                                 |
!                              ___|___
!                             /   |   \
!                            |\___T___/|
!               Neutrons     |    |    |
!              ----------->  |    |____|__________
!                            |    /    |          z
!                            |   /     |
!                            | _/_____ |
!                            |//      \|
!                             \_______/
!                            /
!                           y
!   input:
!   iaxis: cylinder axis. 1: axis x; 2: axis y; 3:axis z
!   r0:  position vector of the point where the distance is calculated
!   dir: Flight direction vector
!   c0:  cylinder center coordinates
!   rc:  cylinder radius
!   hc:  cylinder height (partitioned h / 2 for x> 0 and h / 2 for x <0)
!   output:
!   x1: first intersection with the cylinder. if is zero, there is no intersection
!   x2: second intersection with the cylinder. if is zero, there is no intersection
c

      subroutine discil(iaxis,r0,d,c0,rc,hc,x1,x2,icompl)
      implicit real*8(a-h,o-z)
      dimension r0(3),d(3),c0(3),x(4)
c
      ix=iaxis
      iy=modulo3(iaxis+1)
      iz=modulo3(iaxis+2)
      
!     It raises the quadratic equation to be solved
!     is the intersection of the line with the circle.
      a=d(iy)**2+d(iz)**2
      b=2.d0*(d(iy)*(r0(iy)-c0(iy))+d(iz)*(r0(iz)-c0(iz)))
      c=(r0(iy)-c0(iy))**2+(r0(iz)-c0(iz))**2-rc**2

      do i=1,4
        x(i)=0.d0
      end do

      call bascara(a,b,c,x(1),x(2),dmin,icompl)

! intersection with the cylinder in x
! if x is greater than the limit, there is no intersection and 
! is put to zero
      xx=r0(ix)+x(1)*d(ix)
      if(xx.lt.c0(ix)-hc/2.d0.or.xx.gt.c0(ix)+hc/2.d0) x(1)=0.d0
      xx=r0(ix)+x(2)*d(ix)
      if(xx.lt.c0(ix)-hc/2.d0.or.xx.gt.c0(ix)+hc/2.d0) x(2)=0.d0
c
!...     Find the coordinates of intersection with the limits.
      if (d(ix).ne.0.d0) then
       x(3)=(-r0(ix)+c0(ix)+hc/2.d0)/d(ix)
       x(4)=(-r0(ix)+c0(ix)-hc/2.d0)/d(ix)
      end if
c
!... Examines the distance to the center of the cylinder, for which occurs
! the intersection with limits. If greater than the radius of the cylinder, 
! because the intersection is not, and setting the distance to zero.
      xx=(r0(iy)+x(3)*d(iy)-c0(iy))**2+(r0(iz)+x(3)*d(iz)-c0(iz))**2
      if (xx.ge.rc**2) x(3)=0.d0
      xx=(r0(iy)+x(4)*d(iy)-c0(iy))**2+(r0(iz)+x(4)*d(iz)-c0(iz))**2
      if (xx.ge.rc**2) x(4)=0.d0

!      
!... Order the x. Interested only positive x and belonging to the cylinder.  
      do j=1,3
      xmin=x(j)
      do i=j+1,4
        if(x(i).lt.xmin) then
        xmin=x(i)
        xx=x(j)
        x(j)=xmin
        x(i)=xx
        end if
      end do
      end do
!      
!
!... select the two shorter distances
      x1=0.d0
      x2=0.d0
      do i=1,4
        if(x(i).gt.0.d0) then
        x1=x(i)
        j=i
        goto 1
        end if
      end do
1     continue
      do i=j+1,4
        if(x(i).gt.0.d0) then
        x2=x(i)
        goto 2
        end if
      end do
2     continue  
!
!... Finish tidying strict Ascending
      if (x1.gt.x2) then
        xx=x2
        x2=x1
        x1=xx
      end if
!      

      return
      end
!----------------------------------------------------------------------      
      subroutine bascara(a,b,c,x1,x2,dmin,icompl)
      implicit real*8(a-h,o-z)
      x1=0.d0
      x2=0.d0
      icompl=0
      discr=b**2-4.d0*a*c
      if (discr.le.0.d0.or.a.eq.0.d0) then
        icompl=1
        return
      end if

      x1=(-b+dsqrt(discr))/2.d0/a
      x2=(-b-dsqrt(discr))/2.d0/a
      return
      end
!--------------------------------------------------------------------
      integer*4 function modulo3(i)
      j=mod(i,3)
      if(j.eq.0) j=3
      modulo3=j
      return
      end
      
!--------------------------------------------------------------------
!     Calculate the mean free path of the neutron and the distance to 
!     the next interaction. Modify the weight of the neutron.
!--------------------------------------------------------------------
      SUBROUTINE LAMBDA(XL,RANDOL,WGT,ICAN,IMED,TRANS)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /DATOSSAM/SE,SE0,SM,PS,UTAB(68),SETOT(68),SEM(68),PSCAT(68)
      COMMON /DATOSCAN/SEC,SE0C,SMC,PSC,SETOTC(68),SEMC(68),PSCATC(68)
      COMMON /DISTS/DISCE1,DISAG1,DISCI1,DISVA,DISCI2,DISAG2,DISCE2


      TRANS=DEXP(-SMC*(DISCE1+DISCE2)-SM*(DISAG1+DISAG2)-
     +       SMC*(DISCI1+DISCI2))
       IF (TRANS.EQ.1.D0) THEN
         WGT=0.D0
         RETURN
       END IF
       
       A=1.D0/(1.D0-TRANS)
       X1=(1.D0-DEXP(-SMC*DISCE1))*A
       X2=(1.D0-DEXP(-SMC*DISCE1-SM*DISAG1))*A
       X3=(1.D0-DEXP(-SMC*DISCE1-SM*DISAG1-SMC*DISCI1))*A
       X4=(1.D0-DEXP(-SMC*DISCE1-SM*DISAG1-SMC*(DISCI1+DISCI2)))*A
       X5=(1.D0-DEXP(-SMC*DISCE1-SM*(DISAG1+DISAG2)-
     +    SMC*(DISCI1+DISCI2)))*A
       IF (RANDOL.LT.X1) THEN
         XL=-DLOG(1.D0-RANDOL/A)/SMC
         ICAN=1
         IMED=1				! SCATT. IN CAN EXT 1
       END IF
       IF (RANDOL.GE.X1.AND.RANDOL.LT.X2) THEN
         XL=DISCE1-(DLOG(1.D0-RANDOL/A)+SMC*DISCE1)/SM
         IMED=0				! SCATT. IN SAMPLE 1
       END IF
       IF (RANDOL.GE.X2.AND.RANDOL.LT.X3) THEN
         XL=DISCE1+DISAG1-
     #   (DLOG(1.D0-(RANDOL/A))+SMC*DISCE1+SM*DISAG1)/SMC
         ICAN=1
         IMED=1				! SCATT. IN CAN INT 1
       END IF
       IF (RANDOL.GE.X3.AND.RANDOL.LT.X4) THEN
         XL=DISCE1+DISAG1+DISCI1-
     #   (DLOG(1.D0-(RANDOL/A))+SMC*DISCE1+SM*DISAG1+SMC*DISCI1)/SMC
         ICAN=1
         IMED=1				! SCATT. IN CAN INT 2
       END IF
       IF (RANDOL.GE.X4.AND.RANDOL.LT.X5) THEN
         XL=DISCE1+DISAG1+DISCI1+DISCI2-
     #   (DLOG(1.D0-(RANDOL/A))+SMC*DISCE1+SM*DISAG1+
     #   SMC*(DISCI1+DISCI2))/SM
         IMED=0				! SCATT. IN SAMPLE 2
       END IF
       IF (RANDOL.GE.X5.AND.RANDOL.LT.1.D0) THEN
         XL=DISCE1+DISAG1+DISCI1+DISCI2+DISAG2-
     #   (DLOG(1.D0-(RANDOL/A))+SMC*DISCE1+SM*(DISAG1+DISAG2)+
     #   SMC*(DISCI1+DISCI2))/SMC
         ICAN=1
         IMED=1				! SCATT. IN CAN EXT 2
       END IF

       WGT=WGT*(1.D0-TRANS)

      IF (IMED.EQ.0) THEN
        WGT=WGT*PS
      ELSE
        WGT=WGT*PSC
      END IF

      RETURN
      END

!--------------------------------------------------------------------
!     Calculate the new location of neutron from the previous
!     of the corresponding angles.
!--------------------------------------------------------------------
      SUBROUTINE POSNUE(R,XL,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(2,3),D(2,3)
      DO 10,I=1,3
10     R(2,I)=R(1,I)+XL*D(2,I)
      RETURN
      END

!--------------------------------------------------------------------
!     Generates a random number between 0 and 1.
!--------------------------------------------------------------------
      SUBROUTINE GGUBFS(DSEED,RANDOL)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION RANDOL(1000)
      DATA D2P31M/2147483647.D0/
      DATA D2P31/2147483711.D0/
      DO 1 I=1,1000
      DSEED=DMOD(16807.D0*DSEED,D2P31M)
      R=DSEED/D2P31
      IF (R.EQ.0.D0) R=.0000001D0
      RANDOL(I)=R
1     CONTINUE
      RETURN
      END
!--------------------------------------------------------------------
!     Interpolated input data table
!--------------------------------------------------------------------
      SUBROUTINE VALINTERP(U)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /DATOSSAM/SE,SE0,SM,PS,UTAB(68),SETOT(68),SEM(68),PSCAT(68)
      COMMON /DATOSCAN/SEC,SE0C,SMC,PSC,SETOTC(68),SEMC(68),PSCATC(68)
      DO 1 I=1,68
      IF(UTAB(I)-U)1,2,2
1     CONTINUE
      I=68
2     CONTINUE
      IF (I.EQ.1) THEN
        SE=SETOT(1)
        SM=SEM(1)
        PS=PSCAT(1)
        SEC=SETOTC(1)
        SMC=SEMC(1)
        PSC=PSCATC(1)
        RETURN
      END IF
      X=(U-UTAB(I-1))/(UTAB(I)-UTAB(I-1))
      SE= X*(SETOT(I)-SETOT(I-1))+SETOT(I-1)
      SM= X*(SEM(I)-SEM(I-1))+SEM(I-1)
      PS= X*(PSCAT(I)-PSCAT(I-1))+PSCAT(I-1)
      SEC= X*(SETOTC(I)-SETOTC(I-1))+SETOTC(I-1)
      SMC= X*(SEMC(I)-SEMC(I-1))+SEMC(I-1)
      PSC= X*(PSCATC(I)-PSCATC(I-1))+PSCATC(I-1)
      RETURN
      END

!--------------------------------------------------------------------
!     Calculates the energy transfer When a collision occurs
!--------------------------------------------------------------------
      SUBROUTINE ENERGIA(IMED,U,IFLAG,RANDO)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NQEXP=250,NWEXP=550,NQEXPC=250,NWEXPC=550)
      COMMON /EXPCOORD/QEXP(NQEXP),HWEXP(NWEXP)
      COMMON /EXPCOORDC/QEXPC(NQEXPC),HWEXPC(NWEXPC)
      COMMON /EXPERSAM/DATSAM(NQEXP,NWEXP)
      COMMON /EXPERCAN/DATCAN(NQEXPC,NWEXPC)
      DIMENSION ENDIST(NWEXP),ENDISTC(NWEXPC)
c
      IFLAG=1
      STTT=0.D0
      S1=0.D0
      E0=10.D0**U
      CALL ETKERNEL(IMED,E0,ENDIST,ENDISTC)
      IF (IMED.EQ.0) THEN
        N=NWEXP
        DX=HWEXP(2)-HWEXP(1)
      ELSE
        N=NWEXPC
        DX=HWEXPC(2)-HWEXPC(1)
      ENDIF
      DO 100 I=1,N
      IF (IMED.EQ.0) THEN
         S2=ENDIST(I)
         ELSE
         S2=ENDISTC(I)
      ENDIF
      STTT=.5D0*(S1+S2)*DX+STTT
100   S1=S2
      IF (STTT.EQ.0.D0) THEN
        IFLAG=0
        RETURN
      ENDIF
      S1=0.D0
      S=0.D0
      DO 11 I=1,N
      IF (IMED.EQ.0) THEN
         S2=ENDIST(I)
         ELSE
         S2=ENDISTC(I)
      ENDIF
      S2=S2/STTT
      DS=.5D0*(S1+S2)*DX
      IF (S+DS-RANDO)33,22,22
33    S=S+DS
      S1=S2
11    CONTINUE
      RETURN
22    IF (IMED.EQ.0) THEN
      A=(S2-S1)/DX
      if (a.eq.0.d0) then
       hw=rando*dx+hwexp(i-1)
       goto 44
      end if
      B=(S1*HWEXP(I)-S2*HWEXP(I-1))/DX
      BA=B/A
      RAD=BA*BA+2.D0*(RANDO-S)/A+(HWEXP(I-1))**2.D0+2*BA*(HWEXP(I-1))
      HW=-BA+DSQRT(RAD)
      IF (.NOT.(HW.GT.HWEXP(I-1).AND.HW.LT.HWEXP(I))) HW=-BA-DSQRT(RAD)
44    E=E0-HW
      IF(E.LE.0.D0) E=E0-HWEXP(I-1)
      IF(E.LE.0.D0) E=1.D-4
      U=DLOG10(E)
      END IF
      IF (IMED.EQ.1) THEN
      A=(S2-S1)/DX
       if (a.eq.0.d0) then
       hw=rando*dx+hwexp(i-1)
       goto 55
      end if
      B=(S1*HWEXPC(I)-S2*HWEXPC(I-1))/DX
      BA=B/A
      RAD=BA*BA+2.D0*(RANDO-S)/A+(HWEXPC(I-1))**2.D0+2*BA*(HWEXPC(I-1))
      HW=-BA+DSQRT(RAD)
      IF(.NOT.(HW.GT.HWEXPC(I-1).AND.HW.LT.HWEXPC(I))) HW=-BA-DSQRT(RAD)
55    E=E0-HW
      IF(E.LE.0.D0) E=E0-HWEXPC(I-1)
      IF(E.LE.0.D0) E=1.D-4
      U=DLOG10(E)
      END IF
      RETURN
      END

!--------------------------------------------------------------------
!     Kernel-based energies transferring experimental data.
!     EYE: when energy is larger than that of experiment
!     up-scattering for possible processes, missing data in the
!     zone between the parabola of higher energy and the energy
!     original.
!--------------------------------------------------------------------
      SUBROUTINE ETKERNEL(IMED,E0,ENDIST,ENDISTC)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NQEXP=250,NWEXP=550,NQEXPC=250,NWEXPC=550)
      COMMON /EXPCOORD/QEXP(NQEXP),HWEXP(NWEXP)
      COMMON /EXPCOORDC/QEXPC(NQEXPC),HWEXPC(NWEXPC)
      COMMON /EXPERSAM/DATSAM(NQEXP,NWEXP)
      COMMON /EXPERCAN/DATCAN(NQEXPC,NWEXPC)
      DIMENSION ENDIST(NWEXP),ENDISTC(NWEXPC)
      DATA HQS2M/2.07219093D-03/

      XK0=DSQRT(E0/HQS2M)
! Fix the limits of integration in Q
      IF (IMED.EQ.0) THEN
        N=NWEXP
      ELSE
        N=NWEXPC
      ENDIF
      DO JW=1,N
      IF (IMED.EQ.0) THEN
        ENDIST(JW)=0.D0
        P=XK0**2-HWEXP(JW)/HQS2M
        H=HWEXP(JW)
      ELSE
        ENDISTC(JW)=0.D0
        P=XK0**2-HWEXPC(JW)/HQS2M
        H=HWEXPC(JW)
      ENDIF
      IF(P.LT.0.D0) GOTO 1
      R=DSQRT(P)
          IF (H.LT.0.D0) THEN
              QINF=-XK0+R
              QSUP=XK0+R
            ELSE
              QINF=XK0-R
              QSUP=XK0+R
          END IF
! Integral in Q to determine the kernel. The integrand is Q * S (Q, w)
      IF (IMED.EQ.0) THEN
       DO  JQ=1,NQEXP
         IF(QEXP(JQ).GE.QINF.AND.QEXP(JQ).LT.QSUP.AND.E0.GT.HWEXP(JW))
     #     THEN
             ENDIST(JW)=ENDIST(JW)+QEXP(JQ)*DATSAM(JQ,JW)
     #                  /DSQRT(E0*(E0-HWEXP(JW)))
         END IF
       END DO
      END IF
      IF (IMED.EQ.1) THEN
       DO  JQ=1,NQEXPC
        IF(QEXPC(JQ).GE.QINF.AND.QEXPC(JQ).LT.QSUP.AND.E0.GT.HWEXPC(JW))
     #     THEN
             ENDISTC(JW)=ENDISTC(JW)+QEXPC(JQ)*DATCAN(JQ,JW)
     #                  /DSQRT(E0*(E0-HWEXPC(JW)))
        END IF
       END DO
      END IF

1     CONTINUE
      END DO

      RETURN
      END


!--------------------------------------------------------------------
!     Calculates scattering angles when an interaction. 
!     Change the values of position and angle, preserving 
!     Altima(?) arising from iteration.
!     Calculates new direction versor.
!--------------------------------------------------------------------
      SUBROUTINE ANGULOS(IMED,A,R,D,U0,U,RAND1,RAND2)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,2),R(2,3),D(2,3),DSDW(90)
      DATA PI/3.1415926535897932D0/
      E0=10.D0**U0
      E=10.D0**U
      CALL ANGDIST(IMED,E0,E,DSDW)

      STTT=0.D0
      S1=0.D0
      I=1
      DELTAT=2.D0
      DO 100 T=2.D0,180.D0,DELTAT
      THETA=T/180.D0*PI
      S2=DSDW(I)*E
      S2=S2*2.D0*PI*DSIN(THETA)*PI/180.D0
      STTT=.5D0*(S1+S2)*DELTAT+STTT
      I=I+1
100   S1=S2
      IF (STTT.LE.0.D0) THEN
	RETURN
      ENDIF

      P=0.D0
      S1=0.D0
      I=1
      DO 1 T=2.D0,180.D0,DELTAT
      THETA=T/180.D0*PI
      S2=DSDW(I)*E
      S2=S2*2.D0*PI*DSIN(THETA)*PI/180.D0/STTT
      DS=.5D0*(S1+S2)*DELTAT
      IF (P+DS-RAND1)2,3,3
2     P=P+DS
      S1=S2
      I=I+1
1     CONTINUE
      TH=180.D0
      GOTO 4

3     AA=(S2-S1)/DELTAT

!     If AA = 0 we must avoid division by zero
      if (aa.eq.0.d0) then
       th=rand1*DELTAT+t-DELTAT
       goto 4
      end if
      B=(S1*T-S2*(T-DELTAT))/DELTAT
      BA=B/AA
      RAD=BA*BA+2.D0*(RAND1-P)/AA+(T-DELTAT)**2.D0+2*BA*(T-DELTAT)
      TH=-BA+DSQRT(RAD)
      IF (.NOT.(TH.GT.T-DELTAT.AND.TH.LT.T)) TH=-BA-DSQRT(RAD)
4     A(2,1)=TH
      A(2,2)=360.D0*RAND2

      DO 20,J=1,3
20    R(1,J)=R(2,J)
      A(1,2)=0.D0
      IF (DABS(D(2,3)).LT.1.D0) THEN
      A(1,2)=DACOS(D(2,1)/DSQRT(1.D0-D(2,3)*D(2,3)))*180.D0/PI
      END IF
      A(1,1)=DACOS(D(2,3))*180.D0/PI

      T1=A(1,1)/180.D0*PI
      T2=A(2,1)/180.D0*PI
      F1=A(1,2)/180.D0*PI
      F2=A(2,2)/180.D0*PI
      D(2,1)=DCOS(T1)*DCOS(F1)*DSIN(T2)*DCOS(F2)-
     #     DSIN(F1)*DSIN(T2)*DSIN(F2)+DSIN(T1)*DCOS(F1)*DCOS(T2)
      D(2,2)=DCOS(T1)*DSIN(F1)*DSIN(T2)*DCOS(F2)+
     #     DCOS(F1)*DSIN(T2)*DSIN(F2)+DSIN(T1)*DSIN(F1)*DCOS(T2)
      D(2,3)=DCOS(T1)*DCOS(T2)-DSIN(T1)*DSIN(T2)*DCOS(F2)

      RETURN
      END
!--------------------------------------------------------------------
!     Angular distribution generated from experimental data. 
!     Do not use the synthetic model.
!--------------------------------------------------------------------
      SUBROUTINE ANGDIST(IMED,E0,E,DSDW)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NQEXP=250,NWEXP=550,NQEXPC=250,NWEXPC=550)
      COMMON /EXPCOORD/QEXP(NQEXP),HWEXP(NWEXP)
      COMMON /EXPCOORDC/QEXPC(NQEXPC),HWEXPC(NWEXPC)
      COMMON /EXPERSAM/DATSAM(NQEXP,NWEXP)
      COMMON /EXPERCAN/DATCAN(NQEXPC,NWEXPC)
      DIMENSION DSDW(90)
      DATA PI/3.1415926535897932D0/,HQS2M/2.07219093D-03/
      HW=E0-E
      I=1
      DO T=2.D0,180.D0,2.D0
       THETA=T/180.D0*PI
       Q=DSQRT((E0+E-2.D0*DSQRT(E0*E)*DCOS(THETA))/HQS2M)
					      
       IF (IMED.EQ.0) THEN
         CALL LOCATE(QEXP,NQEXP,Q,JQ)
         CALL LOCATE(HWEXP,NWEXP,HW,JW)
         DSDW(I)=0.D0
         IF(JQ.GE.1.AND.JQ.LE.NQEXP-1.AND.JW.GE.1.AND.JW.LE.NWEXP-1) 
     +     THEN
           CALL INTSQW
     #    (Q,QEXP(JQ),QEXP(JQ+1),DATSAM(JQ,JW),DATSAM(JQ+1,JW),S1)
           CALL INTSQW
     #    (Q,QEXP(JQ),QEXP(JQ+1),DATSAM(JQ,JW+1),DATSAM(JQ+1,JW+1),S2)
           CALL INTSQW
     #    (HW,HWEXP(JW),HWEXP(JW+1),S1,S2,DSDW(I))
         END IF
       END IF

       IF (IMED.EQ.1) THEN
       CALL LOCATE(QEXPC,NQEXPC,Q,JQ)
       CALL LOCATE(HWEXPC,NWEXPC,HW,JW)
       DSDW(I)=0.D0
        IF(JQ.GE.1.AND.JQ.LE.NQEXPC-1.AND.JW.GE.1.AND.JW.LE.NWEXPC-1)
     +    THEN
          CALL INTSQW
     #    (Q,QEXPC(JQ),QEXPC(JQ+1),DATCAN(JQ,JW),DATCAN(JQ+1,JW),S1)
          CALL INTSQW
     #    (Q,QEXPC(JQ),QEXPC(JQ+1),DATCAN(JQ,JW+1),DATCAN(JQ+1,JW+1),S2)
          CALL INTSQW
     #    (HW,HWEXPC(JW),HWEXPC(JW+1),S1,S2,DSDW(I))
        END IF
       END IF
!       DSDW(I)=DSDW(I)*DSIN(THETA)/Q     <----- Was abolished, NO VA
!   Multiply by NO SIN (THETA), because that is done in the 
!   ROUTINE angles. As this program has already converted the experimental data to 
!   probability P(E0,E,THETA) is not required Jacobian in Q.

       I=I+1
       END DO
       RETURN
       END

!--------------------------------------------------------------------
!     Russian Roulette
!--------------------------------------------------------------------
      SUBROUTINE RULRUS(RANDOL,WGT,WCO,IFLAG)
      IMPLICIT REAL * 8 (A-H,O-Z)
      IFLAG=1
      if (wgt.eq.0.d0) then
	iflag=0
	return
      end if
      IF (WGT.LT.WCO) THEN
        IF (RANDOL.LT..5D0) THEN
           WGT=WGT*2.D0
          ELSE
           IFLAG=0
        END IF
      END IF
      RETURN
      END
!--------------------------------------------------------------------
C
!  Calculate the value of the function
C
C
!                 P(E,E',O) * exp( -s (E')*D ) Ef(E')
!                                    T
!  P is the probability of scattering based on experimental data. 
!  The normalized by the effective experimental Section out
C
!  where E 'is the output energy.
C
!  IFLAG=0 : Calculated without correction for attenuation and efficiency
!            detector
!  IFLAG other than 0: complete calculation.
!--------------------------------------------------------------------

      REAL *8 FUNCTION SIGMUE(IFLAG,N,NQ,X)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      PARAMETER (NQEXP=250,NWEXP=550)
      COMMON /DATOSN/EINC,U,THETA
      COMMON /DISTS/DISCE1,DISAG1,DISCI1,DISVA,DISCI2,DISAG2,DISCE2
      COMMON /DATOSSAM/SE,SE0,SM,PS,UTAB(68),SETOT(68),SEM(68),PSCAT(68)
      COMMON /DATOSCAN/SEC,SE0C,SMC,PSC,SETOTC(68),SEMC(68),PSCATC(68)
      COMMON /EXPCOORD/QEXP(NQEXP),HWEXP(NWEXP)
      COMMON /EXPERSAM/DATSAM(NQEXP,NWEXP)
      DATA HQS2M/2.07219093D-03/
      EACT=10.D0**U
      E=10.D0**X
      CALL VALINTERP(U)
      SEACT=SE
      CALL VALINTERP(X)
      HW=EACT-E
      
! Calculate the real Q for this final energy and angle
      IF (N.EQ.1) THEN
          Q=QEXP(NQ)
        ELSE
          Q=DSQRT((E+EACT-2.D0*DSQRT(E*EACT)*DCOS(THETA))/HQS2M)
      END IF

      CALL LOCATE(QEXP,NQEXP,Q,JQ)
      CALL LOCATE(HWEXP,NWEXP,HW,JW)
      S2DIF=0.D0
      IF(JQ.GE.1.AND.JQ.LE.NQEXP-1.AND.JW.GE.1.AND.JW.LE.NWEXP-1) THEN
          CALL INTSQW
     #    (Q,QEXP(JQ),QEXP(JQ+1),DATSAM(JQ,JW),DATSAM(JQ+1,JW),S1)
          CALL INTSQW
     #    (Q,QEXP(JQ),QEXP(JQ+1),DATSAM(JQ,JW+1),DATSAM(JQ+1,JW+1),S2)
          CALL INTSQW
     #    (HW,HWEXP(JW),HWEXP(JW+1),S1,S2,S2DIF)
      END IF
      S2DIF=S2DIF*(SE0/SEACT)*DSQRT(EINC/EACT)
      IF (IFLAG.EQ.0) THEN
          SIGMUE=S2DIF
          RETURN
      END IF

      TRANS=DEXP(-SMC*(DISCE1+DISCE2)-SM*(DISAG1+DISAG2)-
     +      SMC*(DISCI1+DISCI2))
      SIGMUE=S2DIF*TRANS*EFI(X)

      RETURN
      END

!--------------------------------------------------------------------
!     Same as above, for the can
!--------------------------------------------------------------------
      REAL*8 FUNCTION SIGCAN(IFLAG,N,NQ,X)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      PARAMETER (NQEXPC=250,NWEXPC=550)
      COMMON /DATOSN/EINC,U,THETA
      COMMON /DISTS/DISCE1,DISAG1,DISCI1,DISVA,DISCI2,DISAG2,DISCE2
      COMMON /DATOSSAM/SE,SE0,SM,PS,UTAB(68),SETOT(68),SEM(68),PSCAT(68)
      COMMON /DATOSCAN/SEC,SE0C,SMC,PSC,SETOTC(68),SEMC(68),PSCATC(68)
      COMMON /EXPCOORDC/QEXPC(NQEXPC),HWEXPC(NWEXPC)
      COMMON /EXPERCAN/DATCAN(NQEXPC,NWEXPC)
      DATA HQS2M/2.07219093D-03/
      EACT=10**U
      E=10.D0**X
      CALL VALINTERP(U)
      SEACTC=SEC
      CALL VALINTERP(X)
      HW=EACT-E
! Calculate the real Q for this final energy and angle
      IF (N.EQ.1) THEN
          Q=QEXPC(NQ)
        ELSE
          Q=DSQRT((E+EACT-2.D0*DSQRT(E*EACT)*DCOS(THETA))/HQS2M)
      END IF

      CALL LOCATE(QEXPC,NQEXPC,Q,JQ)
      CALL LOCATE(HWEXPC,NWEXPC,HW,JW)
      S2DIF=0.D0

      IF(JQ.GE.1.AND.JQ.LE.NQEXPC-1.AND.JW.GE.1.AND.JW.LE.NWEXPC-1) THEN
        CALL INTSQW
     #  (Q,QEXPC(JQ),QEXPC(JQ+1),DATCAN(JQ,JW),DATCAN(JQ+1,JW),S1)
        CALL INTSQW
     #  (Q,QEXPC(JQ),QEXPC(JQ+1),DATCAN(JQ,JW+1),DATCAN(JQ+1,JW+1),S2)
        CALL INTSQW
     #  (HW,HWEXPC(JW),HWEXPC(JW+1),S1,S2,S2DIF)
      END IF
      S2DIF=S2DIF*(SE0C/SEACTC)*DSQRT(EINC/EACT)
      IF (IFLAG.EQ.0) THEN
          SIGCAN=S2DIF
          RETURN
      END IF

!  Eye, here we need attenuated neutrons can, but not in sample

      TRANS=DEXP(-SMC*(DISCE1+DISCE2)-SM*(DISAG1+DISAG2)-
     +      SMC*(DISCI1+DISCI2))
      SIGCAN=S2DIF*TRANS*EFI(X)

      RETURN
      END

!----------------------------------------------------------------------
!     Approximate efficiency function for 3He tube of an inch 
!     in diameter, normal incidence
!----------------------------------------------------------------------
      REAL*8 FUNCTION EFI(X)
      IMPLICIT REAL * 8 (A-H,O-Z)
      IF (X.LE.-0.5D0) THEN
        EFI=1.D0/(1.D0+DEXP(2.D0*(X+0.5D0)))
      ELSE
        EFI=1.D0/(1.D0+DEXP(1.3D0*(X+0.5D0)))
      END IF
      RETURN
      END
!----------------------------------------------------------------------
      SUBROUTINE LOCATE(XX,N,X,J)
      IMPLICIT REAL * 8 (A-H,O-Z)
      DIMENSION XX(N)
      JL=0
      JU=N+1
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GE.XX(1)).EQV.(X.GE.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GOTO 10
      ENDIF
      IF(X.EQ.XX(1))THEN
        J=1
      ELSE IF(X.EQ.XX(N))THEN
        J=N-1
      ELSE
        J=JL
      ENDIF
      RETURN
      END

!----------------------------------------------------------------------
      SUBROUTINE INTSQW(X,X1,X2,Y1,Y2,S)
      IMPLICIT REAL * 8 (A-H,O-Z)

      S=0.D0
      IF (Y1.EQ.0.D0) THEN
        S=Y2
        RETURN
      END IF
      IF (Y2.EQ.0.D0) THEN
        S=Y1
        RETURN
      END IF
      S=(Y2-Y1)/(X2-X1)*(X-X1)+Y1
      RETURN
      END

