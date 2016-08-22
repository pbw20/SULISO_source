C     PROGRAM UJISO
C
C     Calculates isotopic partition-function ratios for subsets 
C     containing up to 1000 atoms.
C
C     IAN H. WILLIAMS, Department of Chemistry, University of Bath
C
C     version 7: April 2015
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*4 UJI,STOP,PROG
      LOGICAL TSTATE 
      COMMON /USE/IR,IW,NAT,NC,IFF,NT,T(7),QE(2,7),QZ(2,7),QT(2,7),
     *SLM1,SLF1,RCF1,SLM2,SLF2,RCF2,SCALE
      DATA UJI,STOP /'UJI ','STOP'/
C
      IR=5
      IW=6
      WRITE(IW,200)
C
C     PROG = UJI  to run calculations for NISO isotopologues (up to 9)
C     PROG = STOP to finish
C     NISO = number of isotopologues (format I1) excluding parent 
C     IFF  = zero, unless NF = IFF (format I5) frequencies are to be
C            input directly
C
    1 READ(IR,100) PROG,NISO,IFFREQ,SCALE
      IF(PROG.EQ.UJI) THEN
C
C     NT = number of temperatures (degrees K, up to 7 allowed)
C
        READ(IR,101) NT,(T(I),I=1,NT)
C
C     Take each isotopologue in turn
C
        DO 4 ICALC=1,NISO+1
        IF(ICALC.EQ.1) THEN
          SLF1    = 0.0D0
          DO 2  I = 1,NT
          QE(1,I) = 0.0D0
    2     QZ(1,I) = 0.0D0
        ELSE
          SLF2    = 0.0D0
          DO 3  I = 1,NT
          QE(2,I) = 0.0D0
    3     QZ(2,I) = 0.0D0
        ENDIF
        CALL PF(ICALC,TSTATE)
        IF(ICALC.GT.1) CALL RATIO(ICALC,TSTATE)
    4   CONTINUE
      ELSE IF(PROG.EQ.STOP) THEN
        WRITE(IW,201)
        STOP 
      ELSE
        WRITE(IW,202) PROG
        STOP
      ENDIF
      GO TO 1
C
  100 FORMAT(A4,I1,I5,F6.3)
  101 FORMAT(I2,3X,7F10.2)
  200 FORMAT(20('*'),' UJISO: June 2016 ',20('*'))
  201 FORMAT(/30('*'),' The End ',30('*'))
  202 FORMAT('Job fails because keyword ',A,' is unrecognised')
      END
C
      SUBROUTINE PF(ICALC,TSTATE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SKIP,TSTATE,OMIT
      INTEGER IX(9)
      COMMON /USE/IR,IW,NAT,NC,IFF,NT,T(7),QE(2,7),QZ(2,7),QT(2,7),
     *SLM1,SLF1,RCF1,SLM2,SLF2,RCF2,SCALE
      COMMON /MOL/F(4501500),FREQ(3000),X(3,1000),W(1000)
C
      DO 1 L=1,6
    1 IX(L)=0
      WRITE(IW,203) ICALC,SCALE
C
C     Obtain geometry and isotopic masses for this isotopologue
C
      CALL XMI(ICALC)
C
C     Check whether force constants or frequencies are to be read in
C
      NFREQ = NC
      IF (IFF.NE.0) THEN
        NF = IFF
        READ(IR,102) (FREQ(I),I=1,NF)
        WRITE(IW,205)
        WRITE(IW,206) (FREQ(I),I=1,NF)
      ELSE
C
C     Obtain Cartesian force constants for this species
C
        IF (ICALC.EQ.1) CALL FCI
      ENDIF
C
C     Calculate vibrational frequencies (FREQ)
C
      CALL VIBI
C
      TSTATE = .FALSE.
      OMIT   = .FALSE.
      READ(IR,103) ITS,NX,(IX(L),L=1,NX)
      IF (NX.GT.0) THEN
C
C     These are the frequencies to be excluded from evaluation of
C     VP, EX and ZP (but not the transition frequency)
C
        WRITE(IW,207) NX
        WRITE(IW,206) (FREQ(IX(L)),L=1,NX)
      ELSEIF (NX.LT.0) THEN
        OMIT = .TRUE.
      ENDIF
      IF(ITS.NE.0) THEN
C
C     This is a transition structure whose transition frequency is
C     number ITS in the array FREQ
C
        TSTATE = .TRUE.
        RCF = -FREQ(ITS)
        IF (ICALC.EQ.1) THEN
           RCF1  = RCF
        ELSE
           RCF2  = RCF
        ENDIF
        WRITE(IW,204) RCF 
      ENDIF
      IF(OMIT) THEN
        IF(TSTATE) THEN
          WRITE(IW,221)
        ELSE
          WRITE(IW,222)
        ENDIF
      ELSE
        IF(TSTATE) THEN
          WRITE(IW,208)
        ELSE
          WRITE(IW,209)
        ENDIF
      ENDIF
C
C     Loop over frequencies
C
      DO 2  I=1,NFREQ
      SKIP = .FALSE.
      DO 3  L=1,NX
    3 IF(I.EQ.IX(L)) SKIP = .TRUE.
C
C     Exclude this frequency from evaluations of VP, EX and ZP
C
      IF(SKIP) GO TO 2
      FRI = FREQ(I)*SCALE
      IF(FRI.LT.0.0D0) THEN
        IF(OMIT.AND..NOT.TSTATE) GO TO 2
        IF(OMIT.AND.I.NE.ITS)   GO TO 2
        FRI = -FRI
        IF(I.NE.ITS) WRITE(IW,219) FRI
C
C     Imaginary frequencies other than the transition frequency are
C     flagged with a warning but are nonetheless treated as real and
C     included in evaluations of VP, EX and ZP
C
      ENDIF
      IF (ICALC.EQ.1) SLF1 = SLF1 + DLOG(FRI)
      IF (ICALC.GT.1) SLF2 = SLF2 + DLOG(FRI)
C
C     The transition frequency is included in VP but not in EX and ZP
C
C     Loop over temperatures
C
      DO 4  J=1,NT
      U   = 1.438918006D0*FRI/T(J)
      UO2 = U/2.0D0
      IF(I.EQ.ITS) THEN
C
C     Quantum correction for transition vector (Bell tunnelling)
C
        GAMMA = UO2/DSIN(UO2)
        IF (ICALC.EQ.1) THEN
          QT(1,J) = GAMMA
        ELSE
          QT(2,J) = GAMMA
        ENDIF
      ELSE
C
C     ZP is the zero-point energy divided by kT
C     EXCG is the natural log of the 'excitational' partition function
C     Note: the full vibrational partition function (counting from the 
C     classical potential energy minimum) is exp(-U/2)/(1-exp(-U))
C
        IF (ICALC.EQ.1) THEN
          QE(1,J) = QE(1,J) - DLOG(1.0D0 - DEXP(-U))
          QZ(1,J) = QZ(1,J) - U/2.0D0
        ELSE
          QE(2,J) = QE(2,J) - DLOG(1.0D0 - DEXP(-U))
          QZ(2,J) = QZ(2,J) - U/2.0D0
        ENDIF
      ENDIF
    4 CONTINUE
    2 CONTINUE
      RETURN
C
C 102 FORMAT(6F12.4)
  102 FORMAT(6F14.7)
  103 FORMAT(20I4)
  203 FORMAT(/20('*'),' Calculations for isotopologue',I2,1X,20('*')/
     * 'Frequencies scaled by ',F5.3)
  204 FORMAT(/'Transition frequency =',F10.4,'i cm-1')
  205 FORMAT(/'Input frequencies:')
  206 FORMAT(6F12.4)
  207 FORMAT(/'The',I2,' rejected frequencies are:')
  208 FORMAT(/'EX and ZP evaluated over 3N-1 frequencies excluding the',
     * ' transition frequency')
  209 FORMAT(/'EX and ZP evaluated over 3N frequencies')
  219 FORMAT(/'Warning: imaginary frequency ',F8.4,'i cm-1 treated as',
     *' real for evaluation of EX, ZP and the frequency term')
  221 FORMAT(/'Imaginary frequencies (except transition frequency)',
     *' omitted from evaluation of EX, ZP and the frequency term')
  222 FORMAT(/'Imaginary frequencies omitted from evaluation of EX, ZP',
     *' and the frequency term')
      END
C
      SUBROUTINE RATIO(ICALC,TSTATE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL TSTATE
      COMMON /USE/IR,IW,NAT,NC,IFF,NT,T(7),QE(2,7),QZ(2,7),QT(2,7),
     *SLM1,SLF1,RCF1,SLM2,SLF2,RCF2,SCALE
      COMMON /MOL/F(4501500),FREQ(3000),X(3,1000),W(1000)
C
C     Check for internal consistency of calculations:
C     PRODMR & PRODFR should be equal.
C     The ratio of the products of atomic masses is 'light/heavy' but
C     the ratio of the products of frequencies is 'heavy/light'.
C     (N.B. for a TS, the vibrational product VP = PRODFR includes
C           the transition frequency)
C
      PRODMR = DEXP(1.5D0*(SLM1-SLM2))
      PRODFR = DEXP(SLF2-SLF1)
      WRITE(IW,210) ICALC,PRODMR,PRODFR
C
C     Calculate contributions to isotopic partition-function ratio
C     at each temperature
C
      IF(TSTATE) THEN
C
C     The ratio of isotopic transition frequencies is “light/heavy”
C
        RRCF   = RCF1/RCF2
        WRITE(IW,218) ICALC,RRCF
      ENDIF
        DO 1  I=1,NT
        WRITE(IW,211) ICALC,T(I)
        REXC   = DEXP(QE(2,I) - QE(1,I))
        RZPE   = DEXP(QZ(2,I) - QZ(1,I))
        PFR    = PRODMR*REXC*RZPE
        WRITE(IW,212) PFR,PRODMR,REXC,RZPE
        IF(TSTATE) THEN
C
C     Quantum correction to transition-frequency factor
C
          QCRRCF   = QT(1,I)/QT(2,I)
          QPFR     = PFR/QCRRCF
          WRITE(IW,220) QCRRCF,QPFR
        ENDIF
    1   CONTINUE
        RETURN
C
  210 FORMAT(/'Product Rule for isotopologue',I3/3X,
     * 'Isotopic mass term (IM) = ',F20.10/7X,'Frequency term (VP) = ',
     * F20.10)
  211 FORMAT(/'Partition Function Ratio for isotopologues 1 and',
     * I2,' at ',F10.2,' K')
  212 FORMAT(/'Q(heavy)/Q(light) =',8X,'IM',7X,'*',8X,'EX',7X,'*',
     * 8X,'ZP'/4E18.10)
  218 FORMAT(/'Transition frequency ratio TF1/TF',I1,' =',6X,F15.10)
  220 FORMAT(/'Quantum correction to transition frequency ratio =',
     * E18.10/'Quantum corrected isotopic partition ratio =',6X,E18.10)
      END
C
      SUBROUTINE XMI(ICALC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*80 TITLE
      CHARACTER*4 ATOM(50),IAT(1000),IA,LAST
      CHARACTER*1 STAR
C
C     Input isotopic masses and molecular geometry
C
      COMMON /USE/IR,IW,NAT,NC,IFF,NT,T(7),QE(2,7),QZ(2,7),QT(2,7),
     *SLM1,SLF1,RCF1,SLM2,SLF2,RCF2,SCALE
      COMMON /MOL/F(4501500),FREQ(3000),X(3,1000),W(1000)
      DIMENSION Z(50)
      DATA Z/1.007825D0,2.014102D0,3.016050D0,12.0D0,13.003354826D0,
     * 14.003242D0,14.003074002D0,15.00010897D0,15.99491463D0,
     * 16.9991312D0,17.991603D0,18.998405D0,31.97207070D0,33.96786665D0,
     * 34.96885D0,36.9659D0,6.01513D0,7.01601D0,9.01219D0,10.01294D0,
     * 11.009305D0,28.085611D0,30.97376D0,78.9184D0,80.9163D0,
     * 22.98977D0,9999.0D0,24.309846D0,26.98153D0,35.458127D0,
     * 39.09778D0,40.004845D0,65.386887D0,79.906561D0,85.468D0,
     * 126.90435D0,132.9051D0,10.814017D0,1.00790D0,12.01100D0,
     * 14.00670D0,15.99940D0,32.060D0,4.002602D0,20.1797D0,39.948D0,
     * 83.798D0,32.97145843D0,35.96708062D0,11.0114336D0/
      DATA ATOM/'  1H','  2H','  3H',' 12C',' 13C',' 14C',' 14N',
     * ' 15N',' 16O',' 17O',' 18O','   F',' 32S',' 34S','35CL','37CL',
     * ' 6LI',' 7LI',' 9BE',' 10B',' 11B','  SI','   P','79BR','81BR',
     * '  NA','  SP','  MG','  AL','  CL','   K','  CA','  ZN','  BR',
     * '  RB','   I','  CS','   B','   H','   C','   N','   O','   S',
     * '  HE','  NE','  AR','  KR',' 33S',' 36S',' 11C'/
      DATA LAST /'*   '/
C
      IF(ICALC.GT.1) GO TO 4
      READ(IR,104) TITLE
      WRITE(IW,213) TITLE
      READ(IR,108) STAR
C
C   Data for isotopologue 1: symbol and cartesian coordinates in Ångström
C
      N=1
    1 IZ=0
      READ(IR,105) IA,(X(J,N),J=1,3)
      IF(IA.EQ.LAST) GO TO 3
      DO 2 M=1,50
    2 IF(ATOM(M).EQ.IA) IZ=M
      IF(IZ.EQ.0) GO TO 99
      W(N)=Z(IZ)
      WRITE(IW,214) N,IA,W(N),(X(J,N),J=1,3)
      N=N+1
      GO TO 1
C
    3 CONTINUE
C     READ(IR,105) ONEMORE
      NAT=N-1
      NC=NAT*3
      GO TO 7
C
C     Data for subsequent isotopologues: only isotopic symbols required
C     (Cartesian coordinates are retrieved from array X)
C
    4 READ(IR,104) TITLE
      WRITE(IW,213) TITLE
      READ(IR,106) (IAT(I),I=1,NAT)
      DO 6  I=1,NAT
      IZ=0
      IA=IAT(I)
      DO 5  M=1,50
    5 IF(ATOM(M).EQ.IA) IZ=M
      IF(IZ.EQ.0) GO TO 99
      W(I)=Z(IZ)
    6 WRITE(IW,214) I,IA,W(I),(X(J,I),J=1,3)
C
C     Calculate log of product of atomic molecular masses, TLM
C
    7 TLM=0.0D0
      DO 8  I=1,NAT
    8 TLM = TLM + DLOG(W(I))
      IF(ICALC.EQ.1) THEN
        SLM1  = TLM
      ELSE
        SLM2  = TLM
      ENDIF
      RETURN
C
   99 WRITE(IW,215) IA
      STOP
C
  104 FORMAT(64A)
  105 FORMAT(A4,1X,3F12.7)
C 105 FORMAT(A4,1X,3F10.6)
C 105 FORMAT(A4,3F14.7)
  106 FORMAT(16(A4,1X))
  108 FORMAT(A1)
  213 FORMAT('Structure ',A//4X,'Atom',6X,'Mass',8X,'X',13X,'Y',13X,'Z',2X,
     *'(Angstrom)')
  220 FORMAT(A//4X,'Atom',6X,'Mass',8X,'X',13X,'Y',13X,'Z',2X,
     *'(Angstrom)')
  214 FORMAT(I3,2X,A4,2X,F10.6,3(2X,F12.8))
  215 FORMAT(/'Atomic species ',A4,' is not recognised')
      END
C
      SUBROUTINE FCI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /USE/IR,IW,NAT,NC,IFF,NT,T(7),QE(2,7),QZ(2,7),QT(2,7),
     *SLM1,SLF1,RCF1,SLM2,SLF2,RCF2,SCALE
      COMMON /MOL/F(4501500),FREQ(3000),X(3,1000),W(1000)
      NNC= NC*(NC+1)/2
      READ(IR,107) (F(I),I=1,NNC)
      RETURN
  107 FORMAT(1X,3E20.12)
C 107 FORMAT(6F12.8)
      END
C
      SUBROUTINE VIBI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION V(4501500),G(3000)
      COMMON /USE/IR,IW,NAT,NC,IFF,NT,T(7),QE(2,7),QZ(2,7),QT(2,7),
     *SLM1,SLF1,RCF1,SLM2,SLF2,RCF2,SCALE
      COMMON /MOL/F(4501500),FREQ(3000),X(3,1000),W(1000)
      DATA TFACT/2.6436411815482D07/,ZERO/0.0D0/
      L=0
      DO 1  N=1,NAT
      ROOT=DSQRT(W(N))
      DO 1  M=1,3
      L=L+1
    1 G(L)=ROOT
      IJ=0
      DO 2  I=1,NC
      GI=G(I)
      DO 2  J=1,I
      IJ=IJ+1
    2 V(IJ)=F(IJ)/(GI*G(J))
      DO 4  I=1,NC
    4 FREQ(I)=ZERO
      CALL EIGEN(NC,V)
      DO 3  I=1,NC
      TEMP=FREQ(I)
      FREQ(I)=DSQRT(DABS(TFACT*TEMP))
    3 FREQ(I)=DSIGN(FREQ(I),TEMP)
      WRITE(IW,216)
      WRITE(IW,217) (FREQ(I),I=1,NC)
      RETURN
  216 FORMAT(/'Frequencies of normal modes (cm-1)')
  217 FORMAT(6F12.4)
      END
C
      SUBROUTINE EIGEN(N,A)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /USE/IR,IW,NAT,NC,IFF,NT,T(7),QE(2,7),QZ(2,7),QT(2,7),
     *SLM1,SLF1,RCF1,SLM2,SLF2,RCF2,SCALE
      COMMON /MOL/F(4501500),E(3000),X(3,1000),WT(1000)
      DIMENSION W(5,3000),A(45015000)
      EPS3=1.D-30
      ZERO=0.D0
      LL=(N*(N+1))/2+1
      EPS=1.D-8
      IORD=-1
      NM1=N-1
      IF(N.EQ.2) GOTO 80
      NM2=N-2
      DO 70 K=1,NM2
         KP1=K+1
         W(2,K)=A((K*(K+1))/2)
         SUM=0.
         DO 10 J=KP1,N
            W(2,J)=A((J*(J-1))/2+K)
   10    SUM=W(2,J)**2+SUM
         S=SIGN(SQRT(SUM),W(2,KP1))
         W(1,K)=-S
         W(2,KP1)=W(2,KP1)+S
         A(K+(KP1*(KP1-1))/2)=W(2,KP1)
         H=W(2,KP1)*S
         IF(ABS(H).LT.1.D-35) GOTO 70
         SUMM=0.D0
         DO 50 I=KP1,N
            SUM=0.D0
            DO 20 J=KP1,I
   20       SUM=SUM+A(J+(I*(I-1))/2)*W(2,J)
            IF(I.GE.N) GOTO 40
            IP1=I+1
            DO 30 J=IP1,N
   30       SUM=SUM+A(I+(J*(J-1))/2)*W(2,J)
   40       W(1,I)=SUM/H
   50    SUMM=W(1,I)*W(2,I)+SUMM
         U=SUMM*0.5D0/H
         DO 60 J=KP1,N
            W(1,J)=W(2,J)*U-W(1,J)
            DO 60 I=KP1,J
   60    A(I+(J*(J-1))/2)=W(1,I)*W(2,J)+W(1,J)*W(2,I)+A(I+(J*(J-1))/2)
   70 A((K*(K+1))/2)=H
   80 W(2,NM1)=A((NM1*(NM1+1))/2)
      W(2,N)=A((N*(N+1))/2)
      W(1,NM1)=A(NM1+(N*(N-1))/2)
      W(1,N)=0.D0
      GERSCH=ABS(W(2,1))+ABS(W(1,1))
      DO 90 I=1,NM1
   90 GERSCH=MAX(ABS(W(2,I+1))+ABS(W(1,I))+ABS(W(1,I+1)),GERSCH)
      DEL=EPS*GERSCH
      DO 100 I=1,N
         W(3,I)=W(1,I)
  100    E(I)=W(2,I)
      IF(DEL.EQ.ZERO)  GOTO  210
      K=N
  110 L=K
  120 IF(ABS(W(3,L-1)).LT.DEL) GOTO 130
      L=L-1
      IF(L.GT.1)  GOTO 120
  130 IF(L.EQ.K)  GOTO 160
      WW=(E(K-1)+E(K))*0.5D0
      R=E(K)-WW
      Z=SIGN(SQRT(W(3,K-1)**2+R*R),R)+WW
      EE=E(L)-Z
      E(L)=EE
      FF=W(3,L)
      R=SQRT(EE*EE+FF*FF)
      J=L
      GOTO 150
  140 R=SQRT(E(J)**2+W(3,J)**2)
      W(3,J-1)=S*R
      EE=E(J)*C
      FF=W(3,J)*C
  150 C=E(J)/R
      S=W(3,J)/R
      WW=E(J+1)-Z
      E(J)=(FF*C+WW*S)*S+EE+Z
      E(J+1)=C*WW-S*FF
      J=J+1
      IF(J.LT.K) GOTO 140
      W(3,K-1)=E(K)*S
      E(K)=E(K)*C+Z
      GOTO 110
  160 K=K-1
      IF(K.GT.1) GOTO 110
      SORTER=1.D0
      IF(IORD.LT.0) SORTER=-1.D0
      J=N
  170 L=1
      II=1
      LL=1
      DO 190 I=2,J
         IF((E(I)-E(L))*SORTER .GT. 0.D0) GOTO 180
         L=I
         GOTO 190
  180    II=I
         LL=L
  190 CONTINUE
      IF(II.EQ.LL) GOTO 200
      WW=E(LL)
      E(LL)=E(II)
      E(II)=WW
  200 J=II-1
      IF(J.GE.2) GOTO 170
  210 RETURN
      END


