C Corrige mediante factores. Versión del 2011.
C Lee el input de una corrida de Monte Carlo, y le aplica los factores
C de atenuación y scattering múltiple, generando al final un nuevo input
C para la siguiente corrida de Monte Carlo, (o el resultado final).
C La salida del Monte Carlo tiene el mismo número de puntos que el ínput.
C Se invoca como ./COMARI02.EXE N1 N2, donde N1 es el numero de corrida a
C procesar, y N2 en el dato a procesar segun la lista que aparece abajo.

      PROGRAM COMARI
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NQ=250,NW=550)
      DIMENSION QEXP(NQ),HWEXP(NW),DATSAM(NQ,NW),FMUL(NQ,NW),AT(NQ,NW),
     #TOTAL(NQ,NW),SQWG(NQ,NW)
      CHARACTER*50 FILNAME,TIT*15,RUNS(8)*5,SUF*1,NC*5
      EXTERNAL GETARG
c      DATA RUNS
c     +/'14645','14647','14656','14662','14669','14678','14686','14694'/
      DATA PI/3.1415926535897932D0/,HQS2M/2.07219093/
C Lee datos experimentales de muestra y del can
      CALL GETARG(1,SUF)
      CALL GETARG(2,NC)
      READ(NC,*)II


      FILNAME='S'//NC//'.SQE0'
      WRITE(6,*)'[33m[1m File elegido ',FILNAME,'[0m'

      E0=250.0D0

      WRITE(6,*)'[33m[1m Numero de corrida [0m',SUF

      OPEN(1,FILE=FILNAME)
      READ(1,*)
      READ(1,*)
      READ(1,'(10(1PE12.5))')QEXP
      READ(1,*)
      READ(1,'(10(1PE12.5))')HWEXP
      READ(1,*)
      DO I=1,NQ
        READ(1,'(10(1PE12.5))')(DATSAM(I,J),J=1,NW)
      END DO
      CLOSE(1)

C Busca su normalización
      ALFA=0.D0
      DE=HWEXP(2)-HWEXP(1)
      DQ=QEXP(2)-QEXP(1)
      DO 500 I=1,NQ
        SUME=0.D0
        DO 501 J=1,NW
501       SUME=SUME+DATSAM(I,J)*DE
        ALFA=ALFA+SUME*QEXP(I)*DQ
500   CONTINUE
      WRITE(6,*)'ALFA MUESTRA= ',ALFA
      ALFAVIEJO=ALFA

C Lee datos de output de la corrida anterior, de single/total y atenuación
      TIT='MCSMR00.OUT'//SUF
      OPEN (1,FILE=TIT,STATUS='OLD')
      DO 31 I=1,NQ
      READ (1,*)
      READ (1,*)
      DO 31 J=1,NW
31    READ(1,*)X,X,X,X,X,TOTAL(I,J),FMUL(I,J),AT(I,J)
      CLOSE(1)

C Corrige DATSAM por scattering múltiple y atenuación
      DO 40 I=1,NQ
      DO 40 J=1,NW
      IF(AT(I,J).NE.0.D0) DATSAM(I,J)=DATSAM(I,J)*FMUL(I,J)/AT(I,J)
40    CONTINUE


C Nuevo alfa
      ALFA=0.D0
      DE=HWEXP(2)-HWEXP(1)
      DQ=QEXP(2)-QEXP(1)
      DO 510 I=1,NQ
        SUME=0.D0
        DO 511 J=1,NW
511       SUME=SUME+DATSAM(I,J)*DE
        ALFA=ALFA+SUME*QEXP(I)*DQ
510   CONTINUE
      WRITE(6,*)'ALFA NUEVO= ',ALFA

c Reescalea la SQW para que su normalización siga siendo la vieja

      DO I=1,NQ
       DO J=1,NW
          DATSAM(I,J)=DATSAM(I,J)*ALFAVIEJO/ALFA
       END DO
      END DO

C Crea el input para la nueva corrida Monte Carlo

      TIT='S'//NC//'.SQE'//SUF
      OPEN(1,FILE=TIT)
      WRITE(1,*)NQ,NW
      WRITE(1,*)' ESCALA Q '
      WRITE(1,100)QEXP
      WRITE(1,*)' ESCALA HW (eV) '
      WRITE(1,100)HWEXP
      WRITE(1,*)' DATOS MUESTRA '
      DO 1 I=1,NQ
1     WRITE(1,100)(DATSAM(I,J),J=1,NW)
      CLOSE(1)

c Crea la salida para el scilab
      
C HALLA EL MAXIMO DE SQW
c      VSQW=DATSAM(1,1)
c      DO I=1,NQ
c        DO J=1,NW
c          IF (DATSAM(I,J).GT.VSQW)  VSQW=DATSAM(I,J)
c        END DO
c      END DO
C DONDE HAY CEROS EN SQW PONE UN VALOR ASOCIADO AL MAXIMO, PARA
C UNA BUENA REPRESENTACION GRAFICA
c      DO I=1,NQ
c        DO J=1,NW
c          SQWG(I,J)=DATSAM(I,J)
c          IF (DATSAM(I,J).EQ.0.D0)  SQWG(I,J)=1.D0*VSQW
c        END DO
c      END DO
c ---------------------------------------------------------------
c
c          FILNAME='Q'//RUNS(II)
c          OPEN(1,FILE=FILNAME)
c          WRITE(1,100)(QEXP(I),I=1,NQ)
c          CLOSE(1)
c          FILNAME='E'//RUNS(II)
c          OPEN(1,FILE=FILNAME)
c          WRITE(1,100)((HWEXP(I)),I=1,NW)
c          CLOSE(1)
c          FILNAME='S'//RUNS(II)
c          OPEN(1,FILE=FILNAME)
c          DO  I=1,NQ
c          WRITE(1,100)(SQWG(I,J),J=1,NW)
c          END DO
c          CLOSE(1)


100   FORMAT(10(1PE12.5))

      END

