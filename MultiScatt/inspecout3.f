C  Inspecciona los files MCSMR00.OUT del directorio en el que se encuentra.
C  Para compilar ejecutar el comando
C  g77 -O3 inspecout.f -o inspecout.exe -lgrace_np
C  porque usa como interfaz al GRACE.


      IMPLICIT REAL*8(A-H,O-Z)
c Muestras 2 -
c      PARAMETER (NN=144,MM=541,NQ=500)
c Muestra1
      PARAMETER (NE=550,NQ=250)
      REAL*8 EX(NE),QX(NQ)
      REAL*8 SING(NQ,NE),SSAT(NQ,NE),CAN(NQ,NE),XMUL(NQ,NE),TOTAL(NQ,NE)
     *,FMUL(NQ,NE),AT(NQ,NE)

      character*1 SUF,IND,IND2,Z*5,TIT*15,COM*35

      WRITE(6,*)'[33m[1mTipo de inspección[0m'
      WRITE(6,*)'[37;42m[1m1- Una corrida, varios Q         [0m'
      WRITE(6,*)'[37;46m[1m2- Varias corridas, un Q         [0m'
      WRITE(6,*)'[37;42m[1m3- Corrida inicial y final, un Q [0m'
      READ(5,*)ITIPO
      GOTO(1111,2222,3333)ITIPO

1111  WRITE(6,*)'[33m[1mNúmero de corrida a inspeccionar[0m'
      READ(5,'(A)')SUF
C Lee datos de output de la corrida pedida.
      TIT='MCSMR00.OUT'//SUF
      OPEN (1,FILE=TIT,STATUS='OLD')
      DO 31 I=1,NQ
      READ (1,'(7x,f7.4)')QX(I)
      READ (1,*)
      DO 31 J=1,NE
31    READ(1,*)EX(J),SING(I,J),SSAT(I,J),CAN(I,J),XMUL(I,J),TOTAL(I,J),
     *FMUL(I,J),AT(I,J)
      CLOSE(1)
      I=GraceOpenf(2048)
      I=GraceCommandf('load "Default.agr"')

      DO 2000 II=1,4
      WRITE(6,*)'[33m[1mValor aproximado de Q para el que mira[0m'
      READ(5,*)QPED
      CALL LOCATE(QX,NQ,QPED,IPED)
      WRITE(6,*)
     *'[33m[1m Elija Q a visualizar, por el número de orden[0m'
c      IF (IPED-2.GE.1) 
c     * WRITE(6,*)'[37;42m[1m',IPED-2,' - ',QX(IPED-2),'[0m'
c      IF (IPED-1.GE.1) 
c     * WRITE(6,*)'[37;46m[1m',IPED-1,' - ',QX(IPED-1),'[0m'  
      IF (IPED.GE.1) 
     * WRITE(6,*)'[37;42m[1m',IPED,' - ',QX(IPED),'[0m'  
c      IF (IPED+1.GE.1) 
c     * WRITE(6,*)'[37;46m[1m',IPED+1,' - ',QX(IPED+1),'[0m'  
c      IF (IPED+2.GE.1) 
c     * WRITE(6,*)'[37;42m[1m',IPED+2,' - ',QX(IPED+2),'[0m'  
c      READ(5,*)IPED

      write(ind,'(i1)')ii-1
      write(ind2,'(i1)')ii+3
      com='focus g'//ind
      I=GraceCommandf(com)
      write(z,'(f5.2)')QX(IPED)
      com='subtitle "Q= '//z//'"' 
      write(6,*)com
      I=GraceCommandf(com)

      OPEN(1,FILE='AUX.DAT')
      DO J=1,NE
      WRITE(1,999)EX(J),SING(IPED,J),SSAT(IPED,J),CAN(IPED,J),
     *XMUL(IPED,J),TOTAL(IPED,J),FMUL(IPED,J),AT(IPED,J)
      END DO
      CLOSE (1)
      I=GraceCommandf('read nxy "AUX.DAT" ')
      I=GraceCommandf('s1 hidden true')
      I=GraceCommandf('s0 color 1')
      I=GraceCommandf('s2 color 2')
      I=GraceCommandf('s3 color 3')
      I=GraceCommandf('s4 color 4')
      I=GraceCommandf('s0 legend "single"')
      I=GraceCommandf('s2 legend "can"')
      I=GraceCommandf('s3 legend "multiple"')
      I=GraceCommandf('s4 legend "total"')
      I=GraceCommandf('legend on')
      I=GraceCommandf('legend loctype view')
      I=GraceCommandf('legend 0.520229007634, 0.679389312977')

      com='move g'//ind//'.s5 to g'//ind2//'.s0 '
      I=GraceCommandf(com)
      com='move g'//ind//'.s6 to g'//ind2//'.s4 '
      I=GraceCommandf(com)

      I=GraceCommandf('autoscale')
      I=GraceCommandf('redraw')
      com='focus g'//ind2
      I=GraceCommandf(com)
      I=GraceCommandf('autoscale')
      I=GraceCommandf('redraw')

2000  continue      
      GOTO 1000

2222  I=GraceOpenf(2048)
c      I=GraceCommandf('autoscale')
      I=GraceCommandf('load "Default2.agr"')
      I=GraceCommandf('with string')
      I=GraceCommandf('string on')
      I=GraceCommandf('string loctype view')
      I=GraceCommandf('string 0.536585365853, 0.963414634147')
      I=GraceCommandf('string color 1')
      I=GraceCommandf('string rot 0')
      I=GraceCommandf('string font 0')
      I=GraceCommandf('string just 8')
      I=GraceCommandf('string char size 1.350000')

      do 3000 ii=1,8
      write(ind,'(i1)')ii-1
      write(suf,'(i1)')ii
      TIT='MCSMR00.OUT'//SUF
      OPEN (1,FILE=TIT,STATUS='OLD')
      DO 32 I=1,NQ
      READ (1,'(7x,f7.4)')QX(I)
      READ (1,*)
      DO 32 J=1,NE
32    READ(1,*)EX(J),SING(I,J),SSAT(I,J),CAN(I,J),XMUL(I,J),TOTAL(I,J),
     *FMUL(I,J),AT(I,J)
      CLOSE(1)
      if (ii.gt.1) goto 400
      WRITE(6,*)'[33m[1mValor aproximado de Q para el que mira[0m'
      READ(5,*)QPED
      CALL LOCATE(QX,NQ,QPED,IPED)
      WRITE(6,*)
     *'[33m[1m Elija Q a visualizar, por el número de orden[0m'
      IF (IPED-2.GE.1) 
     * WRITE(6,*)'[37;42m[1m',IPED-2,' - ',QX(IPED-2),'[0m'
      IF (IPED-1.GE.1) 
     * WRITE(6,*)'[37;46m[1m',IPED-1,' - ',QX(IPED-1),'[0m'  
      IF (IPED.GE.1) 
     * WRITE(6,*)'[37;42m[1m',IPED,' - ',QX(IPED),'[0m'  
      IF (IPED+1.GE.1) 
     * WRITE(6,*)'[37;46m[1m',IPED+1,' - ',QX(IPED+1),'[0m'  
      IF (IPED+2.GE.1) 
     * WRITE(6,*)'[37;42m[1m',IPED+2,' - ',QX(IPED+2),'[0m'  
      READ(5,*)IPED
      write(z,'(f5.2)')QX(IPED)
      com='string def "p=xxx bar    Q='//z//'"'       
      I=GraceCommandf(com)
400   OPEN(1,FILE='AUX.DAT')
      DO J=1,NE
      WRITE(1,999)EX(J),SING(IPED,J)
      END DO
      CLOSE (1)
      I=GraceCommandf('focus g0')
      I=GraceCommandf('legend on')
      I=GraceCommandf('legend loctype view')
      I=GraceCommandf('legend 0.147709923664, 0.868702290076')
      I=GraceCommandf('read "AUX.DAT"')
      com='s'//ind//' legend "'//suf//'"'
      I=GraceCommandf(com)
      I=GraceCommandf('redraw')
      read(5,*)
      OPEN(1,FILE='AUX.DAT')
      DO J=1,NE
      WRITE(1,999)EX(J),SSAT(IPED,J)
      END DO
      CLOSE (1)
      I=GraceCommandf('focus g1')
      I=GraceCommandf('read "AUX.DAT"')
      I=GraceCommandf('redraw')
      read(5,*)
      OPEN(1,FILE='AUX.DAT')
      DO J=1,NE
      WRITE(1,999)EX(J),XMUL(IPED,J)
      END DO
      CLOSE (1)
      I=GraceCommandf('focus g2')
      I=GraceCommandf('read "AUX.DAT"')
      I=GraceCommandf('redraw')
      read(5,*)
      OPEN(1,FILE='AUX.DAT')
      DO J=1,NE
      WRITE(1,999)EX(J),TOTAL(IPED,J)
      END DO
      CLOSE (1)
      I=GraceCommandf('focus g3')
      I=GraceCommandf('read "AUX.DAT"')
      I=GraceCommandf('redraw')
      read(5,*)
      OPEN(1,FILE='AUX.DAT')
      DO J=1,NE
      WRITE(1,999)EX(J),AT(IPED,J)
      END DO
      CLOSE (1)
      I=GraceCommandf('focus g4')
      I=GraceCommandf('read "AUX.DAT"')
      I=GraceCommandf('redraw')
      read(5,*)
      OPEN(1,FILE='AUX.DAT')
      DO J=1,NE
      WRITE(1,999)EX(J),FMUL(IPED,J)
      END DO
      CLOSE (1)
      I=GraceCommandf('focus g5')
      I=GraceCommandf('read "AUX.DAT"')

      I=GraceCommandf('redraw')


      I=GraceCommandf('legend on')
      I=GraceCommandf('legend loctype view')
      I=GraceCommandf('legend 0.520229007634, 0.679389312977')


3000  continue
      goto 1000

3333  I=GraceOpenf(2048)
c      I=GraceCommandf('autoscale')
      I=GraceCommandf('load "Default0.agr"')
      I=GraceCommandf('with string')
      I=GraceCommandf('string on')
      I=GraceCommandf('string loctype view')
      I=GraceCommandf('string 0.536585365853, 0.963414634147')
      I=GraceCommandf('string color 1')
      I=GraceCommandf('string rot 0')
      I=GraceCommandf('string font 0')
      I=GraceCommandf('string just 8')
      I=GraceCommandf('string char size 1.350000')

      WRITE(6,*)'[33m[1mCorrida final[0m'
      READ(5,*)ii
      write(suf,'(i1)')ii
      TIT='MCSMR00.OUT1'
      OPEN (1,FILE=TIT,STATUS='OLD')
      DO 33 I=1,NQ
      READ (1,'(7x,f7.4)')QX(I)
      READ (1,*)
      DO 33 J=1,NE
33    READ(1,*)EX(J),SING(I,J),SSAT(I,J),CAN(I,J),XMUL(I,J),TOTAL(I,J),
     *FMUL(I,J),AT(I,J)
      CLOSE(1)
      WRITE(6,*)'[33m[1mValor aproximado de Q para el que mira[0m'
      READ(5,*)QPED
      CALL LOCATE(QX,NQ,QPED,IPED)
      WRITE(6,*)
     *'[33m[1m Elija Q a visualizar, por el número de orden[0m'
      IF (IPED-2.GE.1) 
     * WRITE(6,*)'[37;42m[1m',IPED-2,' - ',QX(IPED-2),'[0m'
      IF (IPED-1.GE.1) 
     * WRITE(6,*)'[37;46m[1m',IPED-1,' - ',QX(IPED-1),'[0m'  
      IF (IPED.GE.1) 
     * WRITE(6,*)'[37;42m[1m',IPED,' - ',QX(IPED),'[0m'  
      IF (IPED+1.GE.1) 
     * WRITE(6,*)'[37;46m[1m',IPED+1,' - ',QX(IPED+1),'[0m'  
      IF (IPED+2.GE.1) 
     * WRITE(6,*)'[37;42m[1m',IPED+2,' - ',QX(IPED+2),'[0m'  
      READ(5,*)IPED
      write(z,'(f5.2)')QX(IPED)
      com='string def "Q='//z//'"' 
      I=GraceCommandf(com)
C Lee el single sin atenuación de la primera corrida, que son los datos crudos experimentales
      OPEN(1,FILE='AUX.DAT')
      DO J=1,NE
      WRITE(1,999)EX(J),SSAT(IPED,J)
      END DO
      CLOSE (1)
      I=GraceCommandf('focus g0')
      I=GraceCommandf('legend on')
      I=GraceCommandf('legend loctype view')
      I=GraceCommandf('legend 0.147709923664, 0.868702290076')
      I=GraceCommandf('read "AUX.DAT"')
      com='s0 legend "Raw"'
      I=GraceCommandf(com)
      read(5,*)

C Lee las componentes de la corrida solicitada, para comparar con los datos crudos experimentales
      TIT='MCSMR00.OUT'//SUF
      OPEN (1,FILE=TIT,STATUS='OLD')
      DO 34 I=1,NQ
      READ (1,'(7x,f7.4)')QX(I)
      READ (1,*)
      DO 34 J=1,NE
34    READ(1,*)EX(J),SING(I,J),SSAT(I,J),CAN(I,J),XMUL(I,J),TOTAL(I,J),
     *FMUL(I,J),AT(I,J)
      CLOSE(1)

      OPEN(1,FILE='AUX.DAT')
      DO J=1,NE
      WRITE(1,999)EX(J),SING(IPED,J)
      END DO
      CLOSE (1)
      I=GraceCommandf('legend on')
      I=GraceCommandf('legend loctype view')
      I=GraceCommandf('legend 0.147709923664, 0.868702290076')
      I=GraceCommandf('read "AUX.DAT"')
      com='s1 legend "Single"'
      I=GraceCommandf(com)
      I=GraceCommandf('redraw')
      read(5,*)

      OPEN(1,FILE='AUX.DAT')
      DO J=1,NE
      WRITE(1,999)EX(J),CAN(IPED,J)
      END DO
      CLOSE (1)
      I=GraceCommandf('legend on')
      I=GraceCommandf('legend loctype view')
      I=GraceCommandf('legend 0.147709923664, 0.868702290076')
      I=GraceCommandf('read "AUX.DAT"')
      com='s2 legend "Can"'
      I=GraceCommandf(com)
      I=GraceCommandf('redraw')
      read(5,*)

      OPEN(1,FILE='AUX.DAT')
      DO J=1,NE
      WRITE(1,999)EX(J),XMUL(IPED,J)
      END DO
      CLOSE (1)
      I=GraceCommandf('legend on')
      I=GraceCommandf('legend loctype view')
      I=GraceCommandf('legend 0.147709923664, 0.868702290076')
      I=GraceCommandf('read "AUX.DAT"')
      com='s3 legend "Multiple"'
      I=GraceCommandf(com)
      I=GraceCommandf('redraw')
      read(5,*)

      OPEN(1,FILE='AUX.DAT')
      DO J=1,NE
      WRITE(1,999)EX(J),TOTAL(IPED,J)
      END DO
      CLOSE (1)
      I=GraceCommandf('legend on')
      I=GraceCommandf('legend loctype view')
      I=GraceCommandf('legend 0.147709923664, 0.868702290076')
      I=GraceCommandf('read "AUX.DAT"')
      com='s4 legend "Total"'
      I=GraceCommandf(com)
      I=GraceCommandf('redraw')

c      I=GraceCommandf('legend on')
c      I=GraceCommandf('legend loctype view')
c      I=GraceCommandf('legend 0.520229007634, 0.679389312977')




      I=GraceClosePipeF()

100   FORMAT(10(1PE12.5))
999   FORMAT(1X,12(1PE12.5,1X))
1000  END

C----------------------------------------------------------------------
C     Dado un vector ordenado, y un numero, encuentra entre cuales
C     valores del vector se halla ese numero
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


