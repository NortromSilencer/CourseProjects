C VERSION REVISADA CALCULA 7 DISTANCIAS, COMO RESULTADO DE LA INTERSECCI脫N CON CADA CAN
c ANTES CALCULABA 4 DSTANCIAS PORQUE TENIA EN CUANTA SOLO LOS MEDIOS QUE ATRAVESABA, PERO
C NO QUE PORCION DETALLADA DE ELLOS. ESTO ES IMPORTANTE PORQUE QUEREMOS SABER DISTANCIAS REALES
C Y TIEMPOS DE DEMORA REALES EN LA MUESTRA.
C Corrigio error en SIGMUE Y SIGCAN. Estaba mal el c髆puto de Q, usaba una variable
C E0 que era siempre 0, y hab韆 que usar EACT
C   Scattering multiple un experimento de Argonne

C   Relaci贸n de los espectros que se registran:
C   ESPS(I,J,K):     Espectro acumulante saliente para la energ铆a I-esima
C                    el orden de scattering J-simo, y el Q K-simo. El
C                    ESPS(I,10,K) contiene los 贸rdenes de scattering 10 o
C                    mayor, y el ESPS(I,11,K), la suma de K=2...10 (el
C                    m煤ltiple total).
C   EPSPS2(I,J,K):   Acumula los cuadrados del anterior.
C   RESPS(I,J,K):    Valor medio del ESPS.
C   REPSPS2(I,J,K):  Error de ESPS.
C   ESP1C(I,J):      Scattering simple en el can para la energia I-esima
C                    y el Q K-simo.
C   EPSP1C2(I,J,K):  Acumula los cuadrados del anterior.
C   RESP1C(I,J):     Valor medio de ESP1C
C   REPSP1C2(I,J,K): Error de ESP1C.
C   ESP1CNAS(I,J):   Scattering simple en el can para la energia I-esima
C                    y el Q K-simo, sin atenuaci贸n.
C   ESP1CNAS2(I,J):  Acumula los cuadrados del anterior.
C   RESP1CNAS(I,J):  Valor medio de ESP1CNAS.
C   RESP1CNAS2(I,J): Error de ESP1CNAS.
C   RESPT(I,J):      Scattering total para la energ铆a I-esima y el Q K-simo.
C   REPSPT2(I,J,K):  Acumula los cuadrados del anterior.
C   ESP1NAS(I,J):    Scattering simple sin atenuacion.
C   RESP1NAS(I,J):   Valor medio del anterior.


      PROGRAM MCMARI
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER 
     *(NQ=250,NW=550,NQEXP=250,NWEXP=550,NQEXPC=250,NWEXPC=550)
      REAL*8 KBT
      CHARACTER*11 HORA1,HORA2,TITULO*45,FILMUE*45,FILCAN*45
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

C     Generador de semilla Lahey Fortran
C11    DSEED=RRAND()*10000.
C---------------------------------------------------------------------
C     Para Microsoft Fortran
C11    CALL GETTIM(I,J,K,L)
C      DSEED=L+K*10000./60.
C     Para LINUX ALPHA o INTEL
11    CALL RANDOM_SEED
      CALL RANDOM_NUMBER(DSEED)
      DSEED=DSEED*10000
c     Para g77 Linux

c11    j=time()
c      call srand(j)
c      dseed=10000*rand()

      IF (DSEED.LE.0.D0) GOTO 11
      WRITE(6,*)'SEMILLA INICIAL=',DSEED

c  Varias llamadas a ggubfs para deshacerse de la semilla inicial (para compilador INTEL)

      call ggubfs(dseed,randol)
      call ggubfs(dseed,randol)
      call ggubfs(dseed,randol)
      call ggubfs(dseed,randol)

C Inicializa espectros acumulantes y de salida
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

C Inicializa vector NEVENT: numero de eventos en funcion del numero
C de scatterings J y del tipo de neutron (de can o de muestra)

      DO 12 J=1,10
        DO 12 K=1,2
12        NEVENT(J,K)=0
      NUMCIC0=1
      NTOT=0
      NOUT=0
      NHIS=0
      TPRO=0.D0

C Lectura de secciones eficaces totales de muestra y can, microscopicas
C y macroscopicas, y probabilidades de scattering (frente a absorcion).
C Las secciones eficaces totales que se usan son las experimentales.

      OPEN(1,FILE='mcsminp.dat',STATUS='OLD')
      READ(1,*)
      READ(1,*)
      DO 2 I=1,68
2     READ(1,*)UTAB(I),SETOT(I),SEM(I),PSCAT(I),SETOTC(I),SEMC(I),
     #PSCATC(I)
      CLOSE(1)

C Lee parametros generales

      CALL PARAMETROS(TITULO,FILMUE,FILCAN,NCIC,NT,AH,WH,E0,XLS,WCO,
     +ICONT)
      CALL TIME (HORA1)
      HORA2=HORA1
      


      CALL PANTALLA(HORA1,HORA2,TITULO,RC1,RC3,AH,WH,XLS,FILMUE,
     #NT*NCIC,NTOT,0.D0)
      U0=DLOG10(E0)
      
C Lee datos experimentales de muestra y del can

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

C Pasa las energias a eV
      DO I=1,NWEXP
          HWEXP(I)=HWEXP(I)/1000.D0
      END DO

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

C Pasa las energias a eV
      DO I=1,NWEXPC
          HWEXPC(I)=HWEXPC(I)/1000.D0
      END DO


C Inicializa el vector de transferencias de energia (eV) del analizador
      DO 13 I=1,NW
13    HWANALIZ(I)=HWEXP(I)

C Inicializa el vector de transferencias de impulso (1/A) del analizador

      DO 14 I=1,NQ
14    QSAL(I)=QEXP(I)

C Transforma los files de muestra y can a probabilidad 

      ALFA=0.D0
      DE=HWEXP(2)-HWEXP(1)
      DQ=QEXP(2)-QEXP(1)
      DO 500 I=1,NQEXP
        SUME=0.D0
        DO 501 J=1,NWEXP
501       SUME=SUME+DATSAM(I,J)*DE
        ALFA=ALFA+SUME*QEXP(I)*DQ
500   CONTINUE
      WRITE(6,*)'ALFA MUESTRA= ',ALFA

      DO 201 I=1,NQEXP
      DO 201 J=1,NWEXP
      DATSAM(I,J)=DATSAM(I,J)*(DSQRT(E0*(E0-HWEXP(J)))/HQS2M)
     #       /2.D0/PI/ALFA
201   CONTINUE

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
      WRITE(6,*)'ALFA MUESTRA= ',ALFA

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


C Lee corrida anterior si se define esta como de continuacion

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

C Inicializa indice del vector de numeros random, y llama a la rutina random

      KRAND=0
      CALL GGUBFS(DSEED,RANDOL)

C Comienza el gran loop de NCIC ciclos de NT neutrones

      DO 1000 NUMCIC=NUMCIC0,NCIC
      DO 200 N=1,NT
         NTOT=NTOT+1
         KRAND=KRAND+1
         IF (KRAND.GT.1000) THEN          ! SI SE ACABAN LOS RANDOM GENERA
          CALL GGUBFS(DSEED,RANDOL)       ! NUEVOS
          KRAND=1
         END IF
         NS=1                             ! NUMERO DE SCATTERINGS DEL NEUTRON
         U=U0                             ! LETARGIA INICIAL
         WGT=1.D0                         ! PESO INICIAL
         DELTE=0.D0                       ! CORR. AL HW EFECTIVO
                                          ! POR DEMORA DEBIDO A SCATT. MULT.
         ICAN=0

C Interpola parametros para la energia incidente

         CALL VALINTERP(U)
         SE0=SE
         SE0C=SEC
         IRAND1=MOD(KRAND,1000)+1
         IRAND=MOD(IRAND1,1000)+1

C Posicion de entrada del neutron

         CALL RCERO(RC1,AH,WH,POS,A,RANDOL(IRAND1),RANDOL(IRAND))
         D(1,1)=DSIN(A(2,1)/180.D0*PI)*DCOS(A(2,2)/180.D0*PI)
         D(1,2)=DSIN(A(2,1)/180.D0*PI)*DSIN(A(2,2)/180.D0*PI)
         D(1,3)=DCOS(A(2,1)/180.D0*PI)
         DO 20,I=1,3
20       D(2,I)=D(1,I)

C Comienza loop de scattering multiple

100      CONTINUE
            IRAND=MOD(IRAND,1000)+1
C Inicializa vectores direccion y posicion del neutron (auxiliares)
            VDIR(1)=D(2,1)
            VDIR(2)=D(2,2)
            VDIR(3)=D(2,3)
            VPOS(1)=POS(1,1)
            VPOS(2)=POS(1,2)
            VPOS(3)=POS(1,3)
C Calcula distancias de vuelo en distintos medios segun posicion y direccion
C actual
            CALL DISTAN(1)     
C Sortea la distancia de vuelo y actualiza el peso
            CALL LAMBDA(XL,RANDOL(IRAND),WGT,ICAN,IMED,TRANS)
            IF (NS.EQ.1) THEN
                NHIS=NHIS+1
                TPRO=TPRO+TRANS
            END IF
            IRAND=MOD(IRAND,1000)+1
C Computa correccion a la energia final debido a la demora dentro de la muestra
            IF (NS.GT.1) DELTE=DELTE+XL/XLS/DSQRT(10.D0**U)
C Llama a la ruleta rusa. Si el neutron muere llama uno nuevo
            CALL RULRUS(RANDOL(IRAND),WGT,WCO,IFLAG)
            IF (IFLAG.EQ.0) GOTO 200
C Actualiza posicion del neutron
            CALL POSNUE(POS,XL,D)

C Calcula la contribucion a cada detector
C 1) Calcula DISTAN POS(2,...)
C 2) Calcula angulo theta 
C 3) Calcula la contribucion 1/s(E0) d2s/dOdE (Q,w) exp[-s(E) D]
C    en las energias del analizador
C 4) suma para cada detector
         VPOS(1)=POS(2,1)
         VPOS(2)=POS(2,2)
         VPOS(3)=POS(2,3)
         VDIR(1)=0.D0                     ! Detectores en el plano yz
         LL=NS
         IF (NS.GT.10) LL=10

C Caso del neutron dentro de la muestra  -----------------------------------
      IF (IMED.EQ.0) THEN                                                  |
         IF (ICAN.EQ.0) THEN                                               |
           IF (NS.LE.10) NEVENT(LL,1)=NEVENT(LL,1)+1                       |
         ELSE                                                              |
           IF (NS.LE.10) NEVENT(LL,2)=NEVENT(LL,2)+1                       |
         END IF                                                            |
                                                                           |
         DO 30 K=1,NW                                                      |
C Explora las energias finales. Calcula las energias finales EFIN, tales   |
C que se tenga en cuenta la demora en la muestra para los HWANALIZ(K)      |
C prefijados                                                               |
         EFIN=1.D0/(DELTE-1.D0/DSQRT(E0-HWANALIZ(K)))**2                   |
         UFIN=DLOG10(EFIN)                                                 |
                                                                           |
         DO 30 I=1,NQ                                                      |
         VAL=0.D0                                                          |
C Explora el canal de Q efectivo dado por QSAL(I), y en base a eso         |
C saca el  谩ngulo de salida deducido de la energ铆a efectiva E0-HWANALIZ(K) |
         T=(2*E0-HWANALIZ(K)-HQS2M*QSAL(I)**2)/                            |
     #   (2.D0*DSQRT(E0*(E0-HWANALIZ(K))))                                 |
         IF(DABS(T).GE.1.D0) GOTO 30                                       |
         T=DACOS(T)                                                        |
                                                                           |
         VDIR(2)=DSIN(T)                  ! Direcciones al detector del    |
         VDIR(3)=DCOS(T)                  ! angulo dado                    |
         THETA=DACOS(VDIR(2)*D(2,2)+VDIR(3)*D(2,3))  !Ang. real de interacc|
         CALL DISTAN(1)                   ! Distancias a recorrer          |
                                                                           |
         VAL=SIGMUE(1,LL,I,UFIN)                                           |
         ESPS(K,LL,I)=ESPS(K,LL,I)+VAL*WGT                                 |
         ESPS2(K,LL,I)=ESPS2(K,LL,I)+(VAL*WGT)**2.D0                       |
         IF (LL.EQ.1) ESP1NAS(K,I)=ESP1NAS(K,I)+SIGMUE(0,LL,I,UFIN)*WGT    |
30       CONTINUE                                                          |
      END IF                                                               |
C Fin del caso del neutron dentro de la muestra ____________________________
C Caso del neutron en el can             -----------------------------------
                                                                           |
      IF (IMED.EQ.1) THEN                                                  |
C Genera distribucion angular de scattering del can para la energia dada   |
         IF (NS.LE.10) NEVENT(LL,2)=NEVENT(LL,2)+1                         |
         DO 29 K=1,NW                                                      |
C Explora las energias finales. Calcula las energias finales EFIN, tales   |
C que se tenga en cuenta la demora en la muestra para los HWANALIZ(I)      |
C prefijados                                                               |
         EFIN=1.D0/(DELTE-1.D0/DSQRT(E0-HWANALIZ(K)))**2                   |
         UFIN=DLOG10(EFIN)                                                 |
                                                                           |
         DO 29 I=1,NQ                                                      |
         VAL=0.D0                                                          |
C Explora el canal de Q efectivo dado por QSAL(I), y en base a eso         |
C saca el 谩ngulo de salida deducido de la energ铆a efectiva E0-HWANALIZ(K). |
         T=(2*E0-HWANALIZ(K)-HQS2M*QSAL(I)**2)/                            |
     #   (2.D0*DSQRT(E0*(E0-HWANALIZ(K))))                                 |
         IF(DABS(T).GE.1.D0) GOTO 29                                       |
         T=DACOS(T)                                                        |
         VDIR(2)=DSIN(T)                  ! Direcciones al detector del    |
         VDIR(3)=DCOS(T)                  ! angulo dado                    |
         THETA=DACOS(VDIR(2)*D(2,2)+VDIR(3)*D(2,3))                        |
                                                                           |
         CALL DISTAN(1)                   ! Distancias a recorrer          |
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
C Fin del caso del neutron en el can     -----------------------------------
C Comienza a analizar un nuevo scattering

         NS=NS+1

C Sortea nueva energia y angulo

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

200    CONTINUE

C Calcula valores medios y errores

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

C Suma el total del scattering multiple y su error

      DO 42 K=1,NW
      DO 42 I=1,NQ
      RESPS(K,11,I)=0.D0
      RESPS2(K,11,I)=0.D0
      DO 41 LL=2,10
      RESPS(K,11,I)=RESPS(K,11,I)+RESPS(K,LL,I)
41    RESPS2(K,11,I)=RESPS2(K,11,I)+RESPS2(K,LL,I)**2
42    RESPS2(K,11,I)=DSQRT(RESPS2(K,11,I))

C Suma el scattering total y su error

      DO 51 K=1,NW
      DO 51 I=1,NQ
      RESPT(K,I)=RESPS(K,1,I)+RESP1C(K,I)+RESPS(K,11,I)
51    RESPT2(K,I)=DSQRT(RESPS2(K,1,I)**2+RESP1C2(K,I)**2+
     #            RESPS2(K,11,I)**2)


C Escribe salidas

C      CALL SYSTEM ('ECHO ON')
C      CALL SYSTEM ('IF EXIST MCSMR00.BAK DEL MCSMR00.BAK')
C      CALL SYSTEM ('RENAME MCSMR00.OUT MCSMR00.BAK')
      OPEN (1,FILE='MCSMR00.LOG')
      WRITE (1,*)'NEUTRONES TOTALES=        ',NTOT
      WRITE (1,*)'NEUTRONES FUERA DE RANGO= ',NOUT
      WRITE (1,*)'TRANSMISI脫N PROMEDIO ............ ',RTPRO
      WRITE(1,'(10I8)')NEVENT
      CLOSE(1)

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


C      CALL SYSTEM ('IF EXIST MCSME00.BAK DEL MCSME00.BAK')
C      CALL SYSTEM ('RENAME MCSME00.OUT MCSME00.BAK')
C      OPEN (1,FILE='MCSME00.OUT')
C      WRITE(1,*)'ERRORES'
C      WRITE(1,*)
C      DO 32 I=1,NQ
C      WRITE(1,*)' Q= ',QSAL(I)
C      DO 32 K=1,NW
C32    WRITE(1,999)HWANALIZ(K),RESPS2(K,1,I),RESP1C2(K,I),RESPS2(K,11,I),
C     #            RESPT2(K,I)
C      CLOSE(1)

C      OPEN (1,FILE='CONT.DAT')
C      WRITE (1,*) NUMCIC0
C      WRITE (1,*) NTOT
C      WRITE (1,*) NOUT
C      WRITE (1,*) NHIS
C      WRITE (1,*) TPRO
C      WRITE(1,'(10I7)')NEVENT
C      DO 68 I=1,NQ
C      DO 65 J=1,NW
C65    WRITE(1,999)ESP1C(J,I),ESP1C2(J,I),ESP1NAS(J,I)
C      DO 66 J=1,NW
C66    WRITE(1,999)(ESPS(J,K,I),K=1,11)
C      DO 67 J=1,NW
C67    WRITE(1,999)(ESPS2(J,K,I),K=1,11)
C68    CONTINUE
C      CLOSE(1)

999   FORMAT(1X,12(1PE12.5,1X))


      CALL TIME (HORA1)

      CALL PANTALLA(HORA1,HORA2,TITULO,RC1,RC3,AH,WH,XLS,FILMUE,
     #NT*NCIC,NTOT,RTPRO)

      HORA2=HORA1

1000  CONTINUE
      END

C---------------------------------------------------------------------
      SUBROUTINE PARAMETROS(TITULO,FILMUE,FILCAN,NCIC,NT,AH,WH,E0,XLS,
     +WCO,ICONT)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOMUEST/RC1,HC1,RC2,HC2,RC3,HC3,RC4,HC4,
     +C01(3),C02(3),C03(3),C04(3)
      CHARACTER*45 TITULO,FILMUE,FILCAN
      OPEN(1,FILE='mcsmpar.dat',STATUS='OLD')
      READ(1,'(A)')TITULO
      WRITE(6,'(1x,A)')TITULO
      READ(1,'(A)')FILMUE
      WRITE(6,'(1x,A)')FILMUE
      READ(1,'(A)')FILCAN
      WRITE(6,'(1x,A)')FILCAN
      READ(1,*)NCIC
      WRITE(6,'('' N潞 de ciclos '',I4)')NCIC
      READ(1,*)NT
      WRITE(6,'('' N潞 de neutrones por ciclo '',I6)')NT

      READ(1,*)DECE
      WRITE(6,'('' Diametro externo can externo '',F6.3)')DECE
      READ(1,*)DICE
      WRITE(6,'('' Diametro interno can externo '',F6.3)')DICE
      READ(1,*)HCE
      WRITE(6,'('' Altura can externo '',F6.3)')HCE
      READ(1,*)ETCE
      WRITE(6,'('' Espesor tapa can externo '',F6.3)')ETCE

      READ(1,*)DECI
      WRITE(6,'('' Diametro externo can interno '',F6.3)')DECI
      READ(1,*)DICI
      WRITE(6,'('' Diametro interno can interno '',F6.3)')DICI
      READ(1,*)HCI
      WRITE(6,'('' Altura can interno '',F6.3)')HCI
      READ(1,*)ETCI
      WRITE(6,'('' Espesor tapa can interno '',F6.3)')ETCI

c--------------------------------------------------------------------------

c Medidas derivadas, son las medidas de los 4 cilindros que delimitan materiales
c Cilindro 1
      rc1=dece/2.d0      ! radio cilindro 1
      hc1=hce+2.d0*etce  ! altura cinidro 1
      c01(ix)=0.d0       ! coordenadas del cilindro 1
      c01(iy)=0.d0
      c01(iz)=0.d0
c Cilindro 2
      rc2=dice/2.d0      ! radio cilindro 2
      hc2=hci+2.d0*etci  ! altura cinidro 2
      c02(ix)=0.d0       ! coordenadas del cilindro 2
      c02(iy)=0.d0
      c02(iz)=0.d0
c Cilindro 3
      rc3=deci/2.d0      ! radio cilindro 3
      hc3=hci+2.d0*etci  ! altura cinidro 3
      c03(ix)=0.d0       ! coordenadas del cilindro 3
      c03(iy)=0.d0
      c03(iz)=0.d0
c Cilindro 4
      rc4=dici/2.d0      ! radio cilindro 4
      hc4=hci+2.d0*etci  ! altura cinidro 4
      c04(ix)=0.d0       ! coordenadas del cilindro 4
      c04(iy)=0.d0
      c04(iz)=0.d0

c-------------------------------------------

      READ(1,*)AH
      WRITE(6,'('' Altura del haz '',F6.3,'' cm'')')AH
      READ(1,*)WH
      WRITE(6,'('' Ancho del haz '',F6.3,'' cm'')')WH
      READ(1,*)E0
      WRITE(6,'('' Energ铆a incidente '',1PE12.5)')E0
      READ(1,*)XLS
      WRITE(6,'('' Longitud de vuelo saliente '',1PE12.5)')XLS
      READ(1,*)WCO
      WRITE(6,'('' Peso de corte '',1PE10.3)')WCO
      READ(1,*)ICONT
      WRITE(6,'('' Continuaci贸n '',I3)')ICONT
      CLOSE(1)


      RETURN
      END

C--------------------------------------------------------------------
C     Salida por pantalla usando la utilidad dialog 
C--------------------------------------------------------------------
      SUBROUTINE PANTALLA(HORA1,HORA2,TITULO,R1,R2,AH,WH,XLS,ARFIL,
     #NPED,NTOT,TPRO)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER* 11 HORA1,HORA2,ARFIL*45,TITULO*45

      open(7,file='info.dat')
      WRITE(7,*)'Hora actual ',HORA1, '  Bajada anterior ',HORA2
      WRITE(7,*)titulo
      WRITE(7,*)'RC1= ',R1,'    RC3= ',R2
      WRITE(7,*)'AH=  ',AH, '   WH= ',WH
      WRITE(7,*)'XLS= ',XLS
      WRITE(7,*)'Archivo de muestra: ', ARFIL
      WRITE(7,*)'N潞 total pedido     ',NPED
      WRITE(7,*)'Neutrones corridos  ',NTOT
      WRITE(7,*)'Trasmisi贸n promedio ',TPRO
      close(7)
      call system('./infobox4')

      END

C--------------------------------------------------------------------
C     Calcula la posici贸n sobre la cara anterior de la muestra, por
C     donde ingresa el neutr贸n, y el 谩ngulo de entrada del mismo.
C--------------------------------------------------------------------
      SUBROUTINE RCERO(RADIO,AH,WH,R,A,RAND1,RAND2)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(2,3),A(2,2)
      A(1,1)=0.D0
      A(1,2)=0.D0
C---- Aqu铆 se calcula la posici贸n de entrada
      R(1,1)=AH*(RAND2-.5D0)
      R(1,2)=WH*(RAND1-.5D0)
      R(1,3)=-DSQRT(RADIO**2.D0-R(1,2)**2.D0)
C---- Aqu铆 se coloca el 谩ngulo de entrada
      A(2,1)=0.D0
      A(2,2)=0.D0
      RETURN
      END

c
c  Calcula las distancias de intersecci贸n con los l铆mites de las
C  diferentes zonas.
c
c   +--------------------------------------+
c   |               CAN (1)                |
c   | +----+-+--------------------+-+----+ |
c   | |    | |                    | |    | |
c   | |    | |                    | |    | |
c   | |    | |                    | |    | |
c   | |    | |                    | |    | |
c   | |    | |                    | |    | |
c   | |    | |                    | |    | |
c   | |    | |                    | |    | |  ^ X
c   | |    | |                    | |    | |  |
c   | | A  | |                    | |  A | |  |
c   | | G  |C|                    |C|  G | |  |
c   | | U  |A|                    |A|  U | |  |
c   | | A  |N|      VACIO (4)     |N|  A | |  |
c   | |    | |                    | |    | |  +------------> Z
c   | | 2  |3|                    |3|  2 | |
c   | |    | |                    | |    | |
c   | |    | |                    | |    | |
c   | |    | |                    | |    | |
c   | |    | |                    | |    | |
c   | |    | |                    | |    | |
c   | |    | |                    | |    | |
c   | |    | |                    | |    | |
c   | |    | |                    | |    | |
c   | |    | |                    | |    | |
c   | |    | |                    | |    | |
c   | +----+-+--------------------+-+----+ |
c   |               CAN (1)                |
c   +--------------------------------------+



      subroutine distan(ieje)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOMUEST/RC1,HC1,RC2,HC2,RC3,HC3,RC4,HC4,
     +C01(3),C02(3),C03(3),C04(3)
      COMMON /NEUTRON/D(3),R0(3)
      COMMON /DISTS/DISCE1,DISAG1,DISCI1,DISVA,DISCI2,DISAG2,DISCE2

      ix=ieje
      iy=modulo3(ieje+1)
      iz=modulo3(ieje+2)

c  Medidas basicas en centimetros
c      hce=4.94d0        ! Longitud total externa can agua liviana
c      dece=4.5d0        ! Diametro externo del can externo
c      dice=4.4d0        ! Diametro interno del can externo
c      deci=4.3d0        ! Diametro externo del can interno
c      dici=4.2d0        ! Diametro interno del can interno
c      etce=1.0 d0       ! Espesor de las tapas del can externo
c      etci=0.0 d0       ! Espesor de las tapas del can interno


      s=dsqrt(d(1)**2+d(2)**2+d(3)**2)
      d(1)=d(1)/s
      d(2)=d(2)/s
      d(3)=d(3)/s

C Calcula distancias al cilindro del can

      call discil(ieje,r0,d,c01,rc1,hc1,x11,x21,icompl)
      call discil(ieje,r0,d,c02,rc2,hc2,x12,x22,icompl)
      call discil(ieje,r0,d,c03,rc3,hc3,x13,x23,icompl)
      call discil(ieje,r0,d,c04,rc4,hc4,x14,x24,icompl)

c CUANDO INGRESA EL NEUTRON DESDE AFUERA, PUEDE ESTAR LIGERAMENTE "FUERA DE
C LA MUESTRA", Y EL VALOR DE INTERSECCION CORRECTO CON EL EXTERIOR DEL CAN
C ES X2E, Y NO X1E.

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

C-------------------------------------------------------------------------
C Dado un cilindro, un punto y un recta y una direcci贸n de vuelo,
C calcula los puntos de intersecci贸n de la recta con el cilindro.
C                                  x
C                                 |
C                              ___|___
C                             /   |   \
C                            |\___T___/|
C               Neutrones    |    |    |
C              ----------->  |    |____|__________
C                            |    /    |          z
C                            |   /     |
C                            | _/_____ |
C                            |//      \|
C                             \_______/
C                            /
C                           y
C   Entradas:
C   ieje: eje del cilindro. 1: eje x; 2: eje y; 3:eje z
C   r0:  vector posici贸n del punto desde donde se calcula la distancia
C   dir: vector direcci贸n de vuelo
C   c0:  coordenadas del centro del cilindro
C   rc:  radio del cilindro
C   hc:  altura del cilindro (se reparte h/2 para x>0 y h/2 para x<0)
C   Salidas:
C   x1: primera intersecci贸n con el cilindro. Si es nula, es que no hay 1潞 intersecci贸n
C   x2: segunda intersecci贸n con el cilindro. Si es nula, es que no hay 2潞 intersecci贸n


      subroutine discil(ieje,r0,d,c0,rc,hc,x1,x2,icompl)
      implicit real*8(a-h,o-z)
      dimension r0(3),d(3),c0(3),x(4)
      
      ix=ieje
      iy=modulo3(ieje+1)
      iz=modulo3(ieje+2)
      
C     Plantea la ecuaci贸n cuadr谩tica que hay que resolver      
C     Es la intersecci贸n de la recta con la circunferencia.
      a=d(iy)**2+d(iz)**2
      b=2.d0*(d(iy)*(r0(iy)-c0(iy))+d(iz)*(r0(iz)-c0(iz)))
      c=(r0(iy)-c0(iy))**2+(r0(iz)-c0(iz))**2-rc**2

      do i=1,4
        x(i)=0.d0
      end do

      call bascara(a,b,c,x(1),x(2),dmin,icompl)

c Examina la x para la que ocurre la intersecci贸n con el cilindro.
c Si x es mayor que la x de la tapa, no hay intersecci贸n con el cilindro
c finito, y la pone a cero.
      xx=r0(ix)+x(1)*d(ix)
      if(xx.lt.c0(ix)-hc/2.d0.or.xx.gt.c0(ix)+hc/2.d0) x(1)=0.d0
      xx=r0(ix)+x(2)*d(ix)
      if(xx.lt.c0(ix)-hc/2.d0.or.xx.gt.c0(ix)+hc/2.d0) x(2)=0.d0

c     Calcula las coordenadas de intersecci贸n con las tapas.
      if (d(ix).ne.0.d0) then
       x(3)=(-r0(ix)+c0(ix)+hc/2.d0)/d(ix)
       x(4)=(-r0(ix)+c0(ix)-hc/2.d0)/d(ix)
      end if

c Examina la distancia al centro del cilindro, para la que ocurre
c la intersecci贸n con las tapas. Si es mayor que el radio del cilindro,
c es porque la intersecci贸n no ocurre, y pone la distancia a cero.
      xx=(r0(iy)+x(3)*d(iy)-c0(iy))**2+(r0(iz)+x(3)*d(iz)-c0(iz))**2
      if (xx.ge.rc**2) x(3)=0.d0
      xx=(r0(iy)+x(4)*d(iy)-c0(iy))**2+(r0(iz)+x(4)*d(iz)-c0(iz))**2
      if (xx.ge.rc**2) x(4)=0.d0
      

c Ordeno las x. Interesan solo las x positivas y que pertenecen al cilindro.      
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

      
c selecciona las dos distancias menores
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

c Termina poniendo en orden estricto de menor a mayor
      if (x1.gt.x2) then
        xx=x2
        x2=x1
        x1=xx
      end if

      
      return
      end
      
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

      integer*4 function modulo3(i)
      j=mod(i,3)
      if(j.eq.0) j=3
      modulo3=j
      return
      end

C--------------------------------------------------------------------
C     Calcula el camino libre medio del neutr髇 y la distancia a la
C     siguiente interacci髇. Modifica el peso del neutr髇.
C--------------------------------------------------------------------
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
         IMED=1				! SCATT. EN CAN EXT 1
       END IF
       IF (RANDOL.GE.X1.AND.RANDOL.LT.X2) THEN
         XL=DISCE1-(DLOG(1.D0-RANDOL/A)+SMC*DISCE1)/SM
         IMED=0				! SCATT. EN MUESTRA 1
       END IF
       IF (RANDOL.GE.X2.AND.RANDOL.LT.X3) THEN
         XL=DISCE1+DISAG1-
     #   (DLOG(1.D0-(RANDOL/A))+SMC*DISCE1+SM*DISAG1)/SMC
         ICAN=1
         IMED=1				! SCATT. EN CAN INT 1
       END IF
       IF (RANDOL.GE.X3.AND.RANDOL.LT.X4) THEN
         XL=DISCE1+DISAG1+DISCI1-
     #   (DLOG(1.D0-(RANDOL/A))+SMC*DISCE1+SM*DISAG1+SMC*DISCI1)/SMC
         ICAN=1
         IMED=1				! SCATT. EN CAN INT 2
       END IF
       IF (RANDOL.GE.X4.AND.RANDOL.LT.X5) THEN
         XL=DISCE1+DISAG1+DISCI1+DISCI2-
     #   (DLOG(1.D0-(RANDOL/A))+SMC*DISCE1+SM*DISAG1+
     #   SMC*(DISCI1+DISCI2))/SM
         IMED=0				! SCATT. EN MUESTRA 2
       END IF
       IF (RANDOL.GE.X5.AND.RANDOL.LT.1.D0) THEN
         XL=DISCE1+DISAG1+DISCI1+DISCI2+DISAG2-
     #   (DLOG(1.D0-(RANDOL/A))+SMC*DISCE1+SM*(DISAG1+DISAG2)+
     #   SMC*(DISCI1+DISCI2))/SMC
         ICAN=1
         IMED=1				! SCATT. EN CAN EXT 2
       END IF

       WGT=WGT*(1.D0-TRANS)

      IF (IMED.EQ.0) THEN
        WGT=WGT*PS
      ELSE
        WGT=WGT*PSC
      END IF

      RETURN
      END

C--------------------------------------------------------------------
C     Calcula la nueva ubicaci髇 de neutr髇 a partir de la anterior y
C     de los 醤gulos correspondientes.
C--------------------------------------------------------------------
      SUBROUTINE POSNUE(R,XL,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(2,3),D(2,3)
      DO 10,I=1,3
10     R(2,I)=R(1,I)+XL*D(2,I)
      RETURN
      END

C--------------------------------------------------------------------
C     Genera un n鷐ero al azar entre 0 y 1.
C--------------------------------------------------------------------
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
C--------------------------------------------------------------------
C     Interpola en tabla de datos de input
C--------------------------------------------------------------------
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

C--------------------------------------------------------------------
C     Calcula la energ铆a de dispersi贸n cuando ocurre una interacci贸n.
C--------------------------------------------------------------------
      SUBROUTINE ENERGIA(IMED,U,IFLAG,RANDO)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NQEXP=250,NWEXP=550,NQEXPC=250,NWEXPC=550)
      COMMON /EXPCOORD/QEXP(NQEXP),HWEXP(NWEXP)
      COMMON /EXPCOORDC/QEXPC(NQEXPC),HWEXPC(NWEXPC)
      COMMON /EXPERSAM/DATSAM(NQEXP,NWEXP)
      COMMON /EXPERCAN/DATCAN(NQEXPC,NWEXPC)
      DIMENSION ENDIST(NWEXP),ENDISTC(NWEXPC)

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

C--------------------------------------------------------------------
C     Kernel de transferencia de energias basado en los datos experi-
C     mentales. OJO: cuando la energia es mayor que la del experimento
C     por posibles procesos de upscattering, faltaran datos en la
C     zonita entre la parabola de mayor energia y la de la energia
C     original.
C--------------------------------------------------------------------
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
C Fija los l铆mites de integraci贸n en Q
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
C Integra en Q para determinar el kernel. El integrando es Q*S(Q,w)
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


C--------------------------------------------------------------------
C     Calcula los 谩ngulos de dispersi贸n cuando ocurre una interacci贸n.
C     Cambia los valores de posici贸n y 谩ngulo, conservando los que
C     surgen de la 煤ltima iteraci贸n.
C     Calcula el nuevo versor direcci贸n.
C--------------------------------------------------------------------
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

C     Si AA=0 hay que evitar la division por cero
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
C--------------------------------------------------------------------
C     Genera la distribuci髇 angular a partir de los datos experimentales.
C     No utiliza el modelo sint閠ico.
C--------------------------------------------------------------------
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
C       DSDW(I)=DSDW(I)*DSIN(THETA)/Q     <----- FUE SUPRIMIDO, NO VA
C   NO HAY QUE MULTIPLICAR POR SIN(THETA), PORQUE ESO SE HACE EN LA
C   RUTINA ANGULOS. COMO ESTE PROGRAMA YA CONVIRTIO LOS DATOS EXPERIMENTALES
C   A PROBABILIDAD P(E0,E,THETA) NO SE REQUIERE NINGUN JACOBIANO EN Q.

       I=I+1
       END DO
       RETURN
       END

C--------------------------------------------------------------------
C     Ruleta rusa
C--------------------------------------------------------------------
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
C--------------------------------------------------------------------
C
C  Calcula el valor de la funci贸n 
C
C
C                 P(E,E',O) * exp( -s (E')*D ) Ef(E')
C                                    T
C  P es la probabilidad de scattering basada en los datos experimentales.
C  La normalizamos por la seccione eficaz experimental
C
C  donde E' es la energia de salida.
C
C  IFLAG=0 : Calcula sin atenuacion ni correccion por eficiencia del
C            detector
C  IFLAG distinto de 0: Calculo completo.
C--------------------------------------------------------------------

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
      
C Calcula el Q verdadero para esta energ final y 爊gulo
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

C--------------------------------------------------------------------
C     Idem anterior, para el can
C--------------------------------------------------------------------
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
C Calcula el Q verdadero para esta energ final y 爊gulo
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

c  Ojo, aca podemos necesitar los neutrones atenuados en can, pero no en muestra

      TRANS=DEXP(-SMC*(DISCE1+DISCE2)-SM*(DISAG1+DISAG2)-
     +      SMC*(DISCI1+DISCI2))
      SIGCAN=S2DIF*TRANS*EFI(X)

      RETURN
      END

C----------------------------------------------------------------------
C     Funcion eficiencia aproximada para un tubo de 3He de una pulgada
C     de diametro, a incidencia normal
C----------------------------------------------------------------------
      REAL*8 FUNCTION EFI(X)
      IMPLICIT REAL * 8 (A-H,O-Z)
      IF (X.LE.-0.5D0) THEN
        EFI=1.D0/(1.D0+DEXP(2.D0*(X+0.5D0)))
      ELSE
        EFI=1.D0/(1.D0+DEXP(1.3D0*(X+0.5D0)))
      END IF
      RETURN
      END
C----------------------------------------------------------------------
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

C----------------------------------------------------------------------
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

