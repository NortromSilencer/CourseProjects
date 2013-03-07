! REVISED VERSION 7 CALCULATED DISTANCES, AS A RESULT OF THE INTERSECTION WITH EACH CAN
! ESTIMATIED 4 DISTANCES BEFORE BECAUSE ONLY HAD MUCH IN MEDIA SPANNING BUT
! NO DETAILED THAT PORTION THEREOF. This is important because we WANT TO KNOW the actual distances 
! AND real TIME DELAY IN THE SAMPLE.
! Corrected error in SIGSAMPLE And SIGCAN. It was wrong computation of Q, using a variable
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
      Implicit Real(8) (A-H,O-Z)
      Parameter 
     *(NQ=250,NW=550,NQEXP=250,NWEXP=550,NQEXPC=250,NWEXPC=550)
      Real(8) KBT
      CHARACTER*11 title*45,FILMUE*45,FILCAN*45
      Dimension POS(2,3),D(2,3),A(2,2),randol(1000),
     *ESPS(NW,11,NQ),ESPS2(NW,11,NQ),RESPS(NW,11,NQ),RESPS2(NW,11,NQ),
     *ESP1C(NW,NQ),ESP1C2(NW,NQ),RESP1C(NW,NQ),RESP1C2(NW,NQ),
     *ESP1CNAS(NW,NQ),ESP1CNAS2(NW,NQ),RESP1CNAS(NW,NQ),
     *RESP1CNAS2(NW,NQ),
     *RESPT(NW,NQ),RESPT2(NW,NQ),ESP1NAS(NW,NQ),RESP1NAS(NW,NQ),
     *NEVENT(10,2),QSAL(NQ),HWANALIZ(NW)
      Common /DATOSN/E0,U,theta
      Common /DISTS/DISCE1,DISAG1,DISCI1,DISVA,DISCI2,DISAG2,DISCE2
      Common /DATOSSAM/SE,SE0,SM,PS,UTAB(68),SETOT(68),SEM(68),PSCAT(68)
      Common /DATOSCAN/SEC,SE0C,SMC,PSC,SETOTC(68),SEMC(68),PSCATC(68)
      Common /EXPCOORD/QEXP(NQEXP),HWEXP(NWEXP)
      Common /EXPCOORDC/QEXPC(NQEXPC),HWEXPC(NWEXPC)
      Common /EXPERSAM/DATSAM(NQEXP,NWEXP)
      Common /EXPERCAN/DATCAN(NQEXPC,NWEXPC)
      Common /MUESTRA/KBT
      Common /NEUTRON/VDIR(3),VPOS(3)
      Common /GEOMUEST/RC1,HC1,RC2,HC2,RC3,HC3,RC4,HC4,
     +C01(3),C02(3),C03(3),C04(3)
      Data PI/3.1415926535897932D0/,HQS2M/2.07219093D-03/
!----------------------------------------------------------------------
!     Seed generator For Lahey Fortran Compiler (Lahey Computer System Inc.)
!11    dSeed=Rrand()*10000.
!--------------------------------
!     For Microsoft Fortran (Microsoft Fortran PowerStation, Microsoft)
!11    call GETTIM(I,J,K,L)
!      dSeed=L+K*10000./60.
!     For Linux on DEC Alpha or Intel x86/x64 PC
!11    call RANdoM_SEED
!      call RANdoM_NUMBER(dSeed)
!      dSeed=dSeed*10000
!     For Linux with g77 (GNU Fortran compiler which is part of gcc)

11    j=time()
      call srand(j)
      dSeed=10000*rand()

      if (dSeed <= 0.D0) GOTO 11
      WRITE(6,*)'Initial Seed=',dSeed
!----------------------------------------------------------------------
! Several calls to ggubfs to get rid of the initial seed (for Intel Fortran compiler)

      call ggubfs(dSeed,randol)
      call ggubfs(dSeed,randol)
      call ggubfs(dSeed,randol)
      call ggubfs(dSeed,randol)
!----------------------------------------------------------------------
! Initializes Accumulates and Output spectra
      do 10 J=1,NW
        do 10 L=1,NQ
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
           do 10 K=1,11
             ESPS(J,K,L)=0.D0
             ESPS2(J,K,L)=0.D0
             RESPS(J,K,L)=0.D0
10           RESPS2(J,K,L)=0.D0
!----------------------------------------------------------------------
! GR  : Initializes NEVENT(j,k), i.e.the number of neutrons scattered j times in the sample (k=1) or in the can (k=2)
! LPS : Initializes NEVENT vector: number of events depending on the number of scattering J and type of neutron (from can or sample)             

      do 12 J=1,10
        do 12 K=1,2
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
      do 2 I=1,68
2     READ(1,*)UTAB(I),SETOT(I),SEM(I),PSCAT(I),SETOTC(I),SEMC(I),
     #PSCATC(I)
      CLOSE(1)
!----------------------------------------------------------------------
! Reading Simulation Parameters
      call PARAMETROS(title,FILMUE,FILCAN,NCIC,NT,AH,WH,E0,XLS,WCO,
     +ICONT)
!      call TIME (HORA1)
!      HORA2=HORA1      

      
      
      
      call PANTALLA(title,RC1,RC3,AH,WH,XLS,FILMUE,NT*NCIC,NTOT,0.D0)

      U0=DLOG10(E0)
!.... reading experimental Data for sample

      
      OPEN(1,FILE=FILMUE)
      READ(1,*)
      READ(1,*)
      READ(1,'(10(1PE12.5))')QEXP
      READ(1,*)
      READ(1,'(10(1PE12.5))')HWEXP
      READ(1,*)
      do I=1,NQEXP
        READ(1,'(10(1PE12.5))')(DATSAM(I,J),J=1,NWEXP)
      end do
      CLOSE(1)

!... and converting energy in eV
      do I=1,NWEXP
          HWEXP(I)=HWEXP(I)/1000.D0
      end do
c... reading experimental Data for can
      OPEN(1,FILE=FILCAN)
      READ(1,*)
      READ(1,*)
      READ(1,'(10(1PE12.5))')QEXPC
      READ(1,*)
      READ(1,'(10(1PE12.5))')HWEXPC
      READ(1,*)
      do I=1,NQEXPC
        READ(1,'(10(1PE12.5))')(DATCAN(I,J),J=1,NWEXPC)
      end do
      CLOSE(1)

!... and converting energy in eV
      do I=1,NWEXPC
          HWEXPC(I)=HWEXPC(I)/1000.D0
      end do

!----------------------------------------------------------------------
! Initialize energy transfer vector (eV) for the analyzer
      do 13 I=1,NW
13    HWANALIZ(I)=HWEXP(I)

!... and the momentum transfer vector (1/A).

      do 14 I=1,NQ
14    QSAL(I)=QEXP(I)
!----------------------------------------------------------------------
! Transforms sample and can Data to probabilities. alfa is the normalization. 
!... first the sample
      ALFA=0.D0
      DE=HWEXP(2)-HWEXP(1)
      DQ=QEXP(2)-QEXP(1)
      do 500 I=1,NQEXP
        SUME=0.D0
        do 501 J=1,NWEXP
501       SUME=SUME+DATSAM(I,J)*DE
        ALFA=ALFA+SUME*QEXP(I)*DQ
500   CONTINUE
      WRITE(6,*)'ALFA SAMPLE= ',ALFA

      do 201 I=1,NQEXP
      do 201 J=1,NWEXP
      DATSAM(I,J)=DATSAM(I,J)*(DSQRT(E0*(E0-HWEXP(J)))/HQS2M)
     #       /2.D0/PI/ALFA
201   CONTINUE
!...
      ALFA=0.D0
      DE=HWEXP(2)-HWEXP(1)
      DQ=QEXP(2)-QEXP(1)
      do 900 I=1,NQEXP
        SUME=0.D0
        do 901 K=1,NWEXP
901       SUME=SUME+2.D0*PI/(DSQRT(E0*(E0-HWEXP(K)))/HQS2M)*
     #         DATSAM(I,K)*DE
        ALFA=ALFA+SUME*QEXP(I)*DQ
900   CONTINUE
      WRITE(6,*)'ALFA SAMPLE= ',ALFA
!... then the can
      ALFA=0.D0
      DE=HWEXPC(2)-HWEXPC(1)
      DQ=QEXPC(2)-QEXPC(1)
      do 510 I=1,NQEXPC
        SUME=0.D0
        do 511 J=1,NWEXPC
511       SUME=SUME+DATCAN(I,J)*DE
        ALFA=ALFA+SUME*QEXP(I)*DQ
510   CONTINUE
      WRITE(6,*)'ALFA CAN= ',ALFA

      do 211 I=1,NQEXPC
      do 211 J=1,NWEXPC
      DATCAN(I,J)=DATCAN(I,J)*(DSQRT(E0*(E0-HWEXPC(J)))/HQS2M)
     #       /2.D0/PI/ALFA
211   CONTINUE
!...
      ALFA=0.D0
      DE=HWEXPC(2)-HWEXPC(1)
      DQ=QEXPC(2)-QEXPC(1)
      do 910 I=1,NQEXPC
        SUME=0.D0
        do 911 K=1,NWEXPC
911       SUME=SUME+2.D0*PI/(DSQRT(E0*(E0-HWEXPC(K)))/HQS2M)*
     #         DATCAN(I,K)*DE
        ALFA=ALFA+SUME*QEXPC(I)*DQ
910   CONTINUE
      WRITE(6,*)'ALFA CAN= ',ALFA

!----------------------------------------------------------------------
! reads the previous RUN if this is the continuation

       if (ICONT == 1) then
         OPEN(1,FILE='CONT.DAT')
         READ (1,*) NUMCIC0
         READ (1,*) NTOT
         READ (1,*) NOUT
         READ (1,*) NHIS
         READ (1,*) TPRO
         READ(1,'(10I7)')NEVENT
         do 64 I=1,NQ
         do 61 J=1,NW
61       READ(1,999)ESP1C(J,I),ESP1C2(J,I),ESP1NAS(J,I)
         do 62 J=1,NW
62       READ(1,999)(ESPS(J,K,I),K=1,11)
         do 63 J=1,NW
63       READ(1,999)(ESPS2(J,K,I),K=1,11)
64       CONTINUE
         CLOSE(1)
       end if
!----------------------------------------------------------------------
! Initialize random number vector and call the random number generator routine

      Krand=0
      call GGUBFS(dSeed,randol)
!----------------------------------------------------------------------
! Starts the great loop of NCIC cycles and NT neutrons
c Generally NCIC=1 and the various RUNs are controlled by BATCHITO(.sh)
      do 1000 NUMCIC=NUMCIC0,NCIC
      do 200 N=1,NT
         NTOT=NTOT+1
         Krand=Krand+1
         if (Krand > 1000) then          ! if JUST THE RANdoM GENERATED
          call GGUBFS(dSeed,randol)       ! NEW
          Krand=1
         end if
         NS=1                             ! NUMBER neutron scattering
         U=U0                             ! INITIAL STUPOR
         WGT=1.D0                         ! INITIAL WEIGHT
         DELTE=0.D0                       ! CORR. AL HW effective
                                          ! DELAY DUE TO SCATT. MULT.
         ICAN=0

!... Parameter interpolation for the incident energy 

         call VALINTERP(U)
         SE0=SE
         SE0C=SEC
         Irand1=MOD(Krand,1000)+1
         Irand=MOD(Irand1,1000)+1
C.. Neutron Entry position
!... I think that D is the direction, generally taken along Z

         call RCERO(RC1,AH,WH,POS,A,randol(Irand1),randol(Irand))
         D(1,1)=DSIN(A(2,1)/180.D0*PI)*DCOS(A(2,2)/180.D0*PI)
         D(1,2)=DSIN(A(2,1)/180.D0*PI)*DSIN(A(2,2)/180.D0*PI)
         D(1,3)=DCOS(A(2,1)/180.D0*PI)
         do 20,I=1,3
20       D(2,I)=D(1,I)

!... Multiple scattering loop begins

100      CONTINUE
            Irand=MOD(Irand,1000)+1
! Initializes vector direction and position of the neutron (auxiliary)
            VDIR(1)=D(2,1)
            VDIR(2)=D(2,2)
            VDIR(3)=D(2,3)
            VPOS(1)=POS(1,1)
            VPOS(2)=POS(1,2)
            VPOS(3)=POS(1,3)
! GR  : calculate flight distances  in the various Media considering actual position and direction
! LPS : Flight distances calculated in different ways according to position and heading            
            call DISTAN(1)
! GR :      gives(?) flight distance and changes, lowering it, the weight
! LPS:      Lots flight distance and updates the weight
            call LAMBDA(XL,randol(Irand),WGT,ICAN,IMED,TRANS)
            if (NS == 1) then
                NHIS=NHIS+1
                TPRO=TPRO+TRANS
            end if
            Irand=MOD(Irand,1000)+1
!... computes the ultimate energy correction due to the delay in the sample
            if (NS > 1) DELTE=DELTE+XL/XLS/DSQRT(10.D0**U)
!... call "rulrus" (Russian Roulette), if the neutron dies, a new one is created
            call RUSROULETTE(randol(Irand),WGT,WCO,ifLAG)
            if (ifLAG == 0) GOTO 200
!... update position of the neutron
            call NewPos(POS,XL,D)
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
         if (NS > 10) LL=10
 
!... Case neutron in the sample  -------------------------------------------
      if (IMED == 0) then                                                  |
         if (ICAN == 0) then                                               |
           if (NS <= 10) NEVENT(LL,1)=NEVENT(LL,1)+1                       |
         else                                                              |
           if (NS <= 10) NEVENT(LL,2)=NEVENT(LL,2)+1                       |
         end if                                                            |
!...                                                                       |
         do 30 K=1,NW                                                      |
! Explore the final energies. Calculates EFIN final energies such that     |
! takes into account the delay in the sample for HWANALIZ (K) prefixed     |
                                                                           |
         EFIN=1.D0/(DELTE-1.D0/DSQRT(E0-HWANALIZ(K)))**2                   |
         UFIN=DLOG10(EFIN)                                                 |
                                                                           |
         do 30 I=1,NQ                                                      |
         VAL=0.D0                                                          |
! Explore the actual Q channel QSAL(I), and based on that                  |
! delivers the output angle of the energy deduced effective HWANALIZ E0-(K)|
         T=(2*E0-HWANALIZ(K)-HQS2M*QSAL(I)**2)/                            |
     #   (2.D0*DSQRT(E0*(E0-HWANALIZ(K))))                                 |
         if(DABS(T) >= 1.D0) GOTO 30                                       |
         T=DACOS(T)                                                        |
                                                                           |
         VDIR(2)=DSIN(T)                  ! Directions to the detector     |
         VDIR(3)=DCOS(T)                  ! given angle                    |
         theta=DACOS(VDIR(2)*D(2,2)+VDIR(3)*D(2,3))  !actual interaction Ang. |
         call DISTAN(1)                   ! Distances to travel            |
                                                                           |
         VAL=SIGSAMPLE(1,LL,I,UFIN)                                           |
         ESPS(K,LL,I)=ESPS(K,LL,I)+VAL*WGT                                 |
         ESPS2(K,LL,I)=ESPS2(K,LL,I)+(VAL*WGT)**2.D0                       |
         if (LL == 1) ESP1NAS(K,I)=ESP1NAS(K,I)+SIGSAMPLE(0,LL,I,UFIN)*WGT    |
30       CONTINUE                                                          |
      end if                                                               |
c End case of neutron within the sample     -------------------------------|      
! Case of the neutron in can                --------------------------------
                                                                           |
      if (IMED == 1) then                                                  |
! Generates scattering angular distribution for the given energy in the can|
         if (NS <= 10) NEVENT(LL,2)=NEVENT(LL,2)+1                         |
         do 29 K=1,NW                                                      |
! Explore the final energies. Calculates EFIN final energies such that     |
! takes into account the delay in the sample for HWANALIZ (I) prefixed     |                                                         |
!
         EFIN=1.D0/(DELTE-1.D0/DSQRT(E0-HWANALIZ(K)))**2                   |
         UFIN=DLOG10(EFIN)                                                 |
                                                                           |
         do 29 I=1,NQ                                                      |
         VAL=0.D0                                                          |
! Explore the Q channel QSAL cash given by (I), and based on that takes 
! the output angle deduced from the effective energy-HWANALIZ E0 (K).
         T=(2*E0-HWANALIZ(K)-HQS2M*QSAL(I)**2)/                            |
     #   (2.D0*DSQRT(E0*(E0-HWANALIZ(K))))                                 |
         if(DABS(T) >= 1.D0) GOTO 29                                       |
         T=DACOS(T)                                                        |
         VDIR(2)=DSIN(T)                  ! Directions to the detector     |
         VDIR(3)=DCOS(T)                  ! given angle                    |
         theta=DACOS(VDIR(2)*D(2,2)+VDIR(3)*D(2,3))                        |
                                                                           |
         call DISTAN(1)                   ! Distances to travel            |
                                                                           |
         VAL=SIGCAN(1,LL,I,UFIN)                                           |
         VAL0=SIGCAN(0,LL,I,UFIN)                                          |
         if (LL == 1) then                                                 |
           ESP1C(K,I)=ESP1C(K,I)+VAL*WGT                                   |
           ESP1C2(K,I)=ESP1C2(K,I)+(VAL*WGT)**2.D0                         |
           ESP1CNAS(K,I)=ESP1CNAS(K,I)+VAL0*WGT                            |
           ESP1CNAS2(K,I)=ESP1CNAS2(K,I)+(VAL0*WGT)**2.D0                  |
         else                                                              |
           ESPS(K,LL,I)=ESPS(K,LL,I)+VAL*WGT                               |
           ESPS2(K,LL,I)=ESPS2(K,LL,I)+(VAL*WGT)**2.D0                     |
         end if                                                            |
29       CONTINUE                                                          |
      end if                                                               |
!    End of the neutron case in can-----------------------------------------
!... Start a new scattering analysis

         NS=NS+1
c
!... GR : the neutron has a new energy and direction
!    LPS: Lots new energy and angle         
         Irand=MOD(Irand,1000)+1
         U1=U
         call ENERGY(IMED,U,ifLAG,randol(Irand))
	 if (ifLAG == 0) then
	   NOUT=NOUT+1
	   GOTO 200
	 end if
         Irand1=MOD(Irand,1000)+1
         Irand=MOD(Irand1,1000)+1
         call Angles(IMED,A,POS,D,U1,U,randol(Irand1),randol(Irand))
         call VALINTERP(U)
         GOTO 100
!... when the neutron history ends we can consider a new neutron
200    CONTINUE
!... end of the NT loop
!----------------------------------------------------------------------
!... Calculates mean values and errors
      do 40 K=1,NW
      do 40 I=1,NQ
         if (NTOT > 0) RESP1C(K,I)=ESP1C(K,I)/DBLE(NTOT)
         if (NTOT > 0) RESP1CNAS(K,I)=ESP1CNAS(K,I)/DBLE(NTOT)
         if (NTOT > 0) RESP1NAS(K,I)=ESP1NAS(K,I)/DBLE(NTOT)
         if (NTOT > 1)
     #RESP1C2(K,I)=DSQRT(DABS(ESP1C2(K,I)/DBLE(NTOT)-(ESP1C(K,I)/
     #DBLE(NTOT))**2.D0)/DBLE(NTOT-1))
         if (NTOT > 1)
     #RESP1CNAS2(K,I)=DSQRT(DABS(ESP1CNAS2(K,I)/DBLE(NTOT)-
     #(ESP1CNAS(K,I)/DBLE(NTOT))**2.D0)/DBLE(NTOT-1))
      do 40 LL=1,10
         NN=NEVENT(LL,1)+NEVENT(LL,2)
         if (NTOT > 0) RESPS(K,LL,I)=ESPS(K,LL,I)/DBLE(NTOT)
40       if (NN > 1)
     #RESPS2(K,LL,I)=DSQRT(DABS(ESPS2(K,LL,I)/DBLE(NN)-(ESPS(K,LL,I)/
     #DBLE(NN))**2.D0)/DBLE(NN-1))
      RTPRO=TPRO/DBLE(NHIS)
c
!... Sum total multiple scattering and its error

      do 42 K=1,NW
      do 42 I=1,NQ
      RESPS(K,11,I)=0.D0
      RESPS2(K,11,I)=0.D0
      do 41 LL=2,10
      RESPS(K,11,I)=RESPS(K,11,I)+RESPS(K,LL,I)
41    RESPS2(K,11,I)=RESPS2(K,11,I)+RESPS2(K,LL,I)**2
42    RESPS2(K,11,I)=DSQRT(RESPS2(K,11,I))
c
!... Sum the total scattering and its error

      do 51 K=1,NW
      do 51 I=1,NQ
      RESPT(K,I)=RESPS(K,1,I)+RESP1C(K,I)+RESPS(K,11,I)
51    RESPT2(K,I)=DSQRT(RESPS2(K,1,I)**2+RESP1C2(K,I)**2+
     #            RESPS2(K,11,I)**2)
      

! Write outputs
!      call SYSTEM ('ECHO ON')
!      call SYSTEM ('if EXIST MCSMR00.BAK DEL MCSMR00.BAK')
!      call SYSTEM ('RENAME MCSMR00.OUT MCSMR00.BAK')
! Write on the LOG file
      OPEN (1,FILE='MCSMR00.LOG')
      WRITE (1,*)'Total Neutron Number=     ',NTOT
      WRITE (1,*)'NEUTRONES FUERA DE RANGO= ',NOUT
      WRITE (1,*)'TRANSMISI?“N PROMEDIO ............ ',RTPRO
      WRITE(1,'(10I8)')NEVENT
      CLOSE(1)
c
      OPEN (1,FILE='MCSMR00.OUT')
      do 31 I=1,NQ
      WRITE(1,*)' Q= ',QSAL(I)
      WRITE(1,998)
998   FORMAT('    HW          SINGLE       S/ATEN         CAN    
     #MULTIPLE        TOTAL   SING/TOT         ATEN    ATENC')

      do 31 K=1,NW
      FMUL=0.D0
      AT=0.D0
      ATC=0.D0
      if(RESPT(K,I).NE.0.D0) FMUL=RESPS(K,1,I)/RESPT(K,I)
      if(RESP1NAS(K,I).NE.0.D0) AT=RESPS(K,1,I)/RESP1NAS(K,I)
      if(RESP1CNAS(K,I).NE.0.D0) ATC=RESP1C(K,I)/RESP1CNAS(K,I)
31    WRITE(1,999)HWANALIZ(K),RESPS(K,1,I),RESP1NAS(K,I),RESP1C(K,I),
     #            RESPS(K,11,I),RESPT(K,I),FMUL,AT,ATC
      CLOSE(1)

!    ------------------------------------------------------------------
!      call SYSTEM ('if EXIST MCSME00.BAK DEL MCSME00.BAK')
!      call SYSTEM ('RENAME MCSME00.OUT MCSME00.BAK')
!      OPEN (1,FILE='MCSME00.OUT')
!      WRITE(1,*)'ERRORES'
!      WRITE(1,*)
!      do 32 I=1,NQ
!      WRITE(1,*)' Q= ',QSAL(I)
!      do 32 K=1,NW
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
!      do 68 I=1,NQ
!      do 65 J=1,NW
!65    WRITE(1,999)ESP1C(J,I),ESP1C2(J,I),ESP1NAS(J,I)
!      do 66 J=1,NW
!66    WRITE(1,999)(ESPS(J,K,I),K=1,11)
!      do 67 J=1,NW
!67    WRITE(1,999)(ESPS2(J,K,I),K=1,11)
!68    CONTINUE
!      CLOSE(1)      
!      ERAT-------------------------------------------------------------
999   FORMAT(1X,12(1PE12.5,1X))

!      call TIME (HORA1)

!      call PANTALLA(HORA1,HORA2,TITULO,RC1,RC3,AH,WH,XLS,FILMUE,
!     #NT*NCIC,NTOT,RTPRO)
!      call PANTALLA(TITULO,RC1,RC3,AH,WH,XLS,FILMUE,NT*NCIC,NTOT,RTPRO)
      call PANTALLA(title,RC1,RC3,AH,WH,XLS,FILMUE,NT*NCIC,NTOT,RTPRO)

!      HORA2=HORA1

1000  CONTINUE
      end
!---------------------------------------------------------------------
! Reading Parameters from mcsmpar.dat
!---------------------------------------------------------------------
      Subroutine PARAMETROS(title,FILMUE,FILCAN,NCIC,NT,AH,WH,E0,XLS,
     +WCO,ICONT)
      Implicit Real(8) (A-H,O-Z)
      Common /GEOMUEST/RC1,HC1,RC2,HC2,RC3,HC3,RC4,HC4,
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

      return
      end
      
!--------------------------------------------------------------------
!     GR  : write Parameters on visual shell .sh
!     LPS : Screen output using dialog utility      
!--------------------------------------------------------------------
      Subroutine PANTALLA(title,R1,R2,AH,WH,XLS,ARFIL,NPED,NTOT,TPRO)
      Implicit Real(8) (A-H,O-Z)
      CHARACTER* 11 ARFIL*45,title*45

      open(7,file='info.dat')
!     WRITE(7,*)'current Time ',HORA1, '  previous drop ',HORA2
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

      end

!--------------------------------------------------------------------
!    Calculates the position on the anterior sample face, where
!    the neutron enters, and its entry angle.
!--------------------------------------------------------------------
      Subroutine RCERO(RADIO,AH,WH,R,A,rand1,rand2)
      Implicit Real(8) (A-H,O-Z)
      Dimension R(2,3),A(2,2)
      A(1,1)=0.D0
      A(1,2)=0.D0
!---- Aqua-position is calculated input
      R(1,1)=AH*(rand2-.5D0)
      R(1,2)=WH*(rand1-.5D0)
      R(1,3)=-DSQRT(RADIO**2.D0-R(1,2)**2.D0)
!---- Aqua-placed entry angle
      A(2,1)=0.D0
      A(2,2)=0.D0
      return
      end

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



      Subroutine distan(iaxis)
      Implicit Real(8) (A-H,O-Z)
      Common /GEOMUEST/RC1,HC1,RC2,HC2,RC3,HC3,RC4,HC4,
     +C01(3),C02(3),C03(3),C04(3)
      Common /NEUTRON/D(3),R0(3)
      Common /DISTS/DISCE1,DISAG1,DISCI1,DISVA,DISCI2,DISAG2,DISCE2

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

      if(x21 > 0.d0 .and. x21 > x1e) x1e=x2e

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

      Subroutine discil(iaxis,r0,d,c0,rc,hc,x1,x2,icompl)
      Implicit Real(8)(a-h,o-z)
      Dimension r0(3),d(3),c0(3),x(4)
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
      if(xx.lt.c0(ix)-hc/2.d0.or.xx > c0(ix)+hc/2.d0) x(1)=0.d0
      xx=r0(ix)+x(2)*d(ix)
      if(xx.lt.c0(ix)-hc/2.d0.or.xx > c0(ix)+hc/2.d0) x(2)=0.d0
c
!...     Find the coordinates of intersection with the limits.
      if (d(ix).ne.0.d0) then
       x(3)=(-r0(ix)+c0(ix)+hc/2.d0)/d(ix)
       x(4)=(-r0(ix)+c0(ix)-hc/2.d0)/d(ix)
      end if
c
!... Examines the distance to the center of the cylinder, for which occurs
! the intersection with limits. if greater than the radius of the cylinder, 
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
        if(x(i) > 0.d0) then
        x1=x(i)
        j=i
        goto 1
        end if
      end do
1     continue
      do i=j+1,4
        if(x(i) > 0.d0) then
        x2=x(i)
        goto 2
        end if
      end do
2     continue  
!
!... Finish tidying strict Ascending
      if (x1 > x2) then
        xx=x2
        x2=x1
        x1=xx
      end if
!      

      return
      end
!----------------------------------------------------------------------      
      Subroutine bascara(a,b,c,x1,x2,dmin,icompl)
      Implicit Real(8)(a-h,o-z)
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
      Integer(4) Function modulo3(i)
      j=mod(i,3)
      if(j.eq.0) j=3
      modulo3=j
      return
      end
      


























