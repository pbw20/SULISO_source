C     PROGRAM LIPFR (force constants input as 3E20.12)
C
C     Calculates isotopic partition function ratios for species 
C     with up to a maximum number of atoms set in the header file
C     CAM_SIZE.
C
C     IAN H. WILLIAMS, Department of Chemistry, University of Bath
C
C     version 1.9: 16 October 2016
C
C     This version uses "CAM_SIZE" to determine the array sizes.
C     Number of atoms in system is set at (for example) MAXNAT = 1000
C     -> number of Cartesian coordinates, MAXNC = MAXNAT*3
C     -> No. of valence coords, MAXNI = MAXNAT*3 + MAXNAT/2
C     -> no. of independent internal coords, MAXNINT = MAXNC-6
C     -> max. length linear array MAXDIAG = MAXNC*(MAXNC+1)/2
C     -> parameter(MAXAPO =MAXNAT + 1)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'CAM_SIZE'
      CHARACTER*4 IPFR,STOP,PROG,IFLIB,LIBR
      LOGICAL LIBRATIONS,BIGELEISEN,TSTATE
      COMMON /USE/IR,IW,NAT,NC,IFFREQ,NT,T(5),TLN(5),WM,TLM,SYM,GEL
      DATA IPFR,STOP,LIBR /'IPFR','STOP','LIBR'/
C
      IR=5
      IW=6
      WRITE(IW,200)
C
C     PROG= IPFR to run calculations for NISO isotopologues (up to 9)
C     PROG= STOP to finish
C
    1 READ(IR,100) PROG,NISO,IFLIB
      IF(PROG.EQ.IPFR) GO TO 11
      IF(PROG.EQ.STOP) GO TO 12
      WRITE(IW,202) PROG
      STOP
C
C     NT = number of temperatures (degrees K, up to 5 allowed)
C
   11 READ(IR,101) NT,(T(I),I=1,NT)
      DO 2  I=1,NT
    2 TLN(I)=DLOG(T(I))
      LIBRATIONS = .FALSE.
      IF(IFLIB.EQ.LIBR) LIBRATIONS = .TRUE.
C
C     Take each isotopologue in turn
C
      DO 3 ICALC=1,NISO
      CALL PF(ICALC,LIBRATIONS,BIGELEISEN,TSTATE)
    3 IF(ICALC.GT.1) CALL RATIO(ICALC,BIGELEISEN,TSTATE)
C
      GO TO 1
   12 WRITE(IW,201)
      STOP 
C
  100 FORMAT(A4,I2,1X,A4)
  101 FORMAT(I2,3X,5F10.4)
  200 FORMAT(20('*'),' LIPFR: 16 October 2016 ',20('*'))
  201 FORMAT(/30('*'),' The End ',30('*'))
  202 FORMAT('JOB FAILS BECAUSE KEYWORD ',A,' IS UNRECOGNISED')
      END
C
      SUBROUTINE PF(ICALC,LIBRATIONS,BIGELEISEN,TSTATE)
      INCLUDE 'CAM_SIZE'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SKIP,TSTATE,LIBRATIONS,BIGELEISEN
      INTEGER IX(9)
      DOUBLE PRECISION PM(3),ZP(5),EX(5)
      COMMON /USE/IR,IW,NAT,NC,IFFREQ,NT,T(5),TLN(5),WM,TLM,SYM,GEL
      COMMON /MOL/F(MAXDIAG),FREQ(MAXNC),X(3,MAXNAT),W(MAXNAT)
      COMMON /BOTH/QT(2,5),QR(2,5),QE(2,5),QZ(2,5),QTUN(2,5),
     *SLM1,TMLN1,SLI1,SLF1,QT1,QR1,QE1,QZ1,FTSL1,RCF1,
     *SLM2,TMLN2,SLI2,SLF2,QT2,QR2,QE2,QZ2,FTSL2,RCF2
      DATA CTR,CROT,R,HCOK,PILN/-3.665071893D0,-3.18871834D0,
     *1.985836D-03,1.438918006D0,1.144729886D0/
      DATA ZERO, ONE,  TWO,  ZP5  ,ONEP5,TWOP5
     *    /0.0D0,1.0D0,2.0D0,0.5D0,1.5D0,2.5D0/
C
      WRITE(IW,203) ICALC
C
C     Obtain geometry and isotopic masses for this isotopologue
C
      CALL XMI(ICALC)
C
C     WM is the total molecular mass
C     TLM is the natural log of the product of the atomic masses
C
      IF(LIBRATIONS) GO TO 9
      IF(ICALC.EQ.1) THEN
        SLM1  = TLM
        TMLN1 = DLOG(WM)
      ELSE
        SLM2  = TLM
        TMLN2 = DLOG(WM)
      ENDIF
C
      NR=0
      IF(NAT.EQ.1) GO TO 8
C
C     CALCULATE PRINCIPAL MOMENTS OF INERTIA, PM
C
      CALL MOMIN(NR,PM)
C
C     PMLN is the natural log of the product of the moments of inertia
C
      PMLN = ZERO
      DO 1 I=1,NR
    1 PMLN = PMLN + DLOG(PM(I))
C
      IF(ICALC.EQ.1) SLI1 = PMLN
      IF(ICALC.GT.1) SLI2 = PMLN
C
C     Check whether force constants or frequencies are to be read in
C
    9 NFREQ = NC
      IF (IFFREQ.NE.0) THEN
        NFREQ = IFFREQ
        READ(IR,102) (FREQ(I),I=1,NFREQ)
        WRITE(IW,205)
        WRITE(IW,206) (FREQ(I),I=1,NFREQ)
		GO TO 10
      ENDIF
C
C     Obtain Cartesian force constants for this species
C
      IF (ICALC.EQ.1) CALL FCI
C
C     Calculate vibrational frequencies (FREQ)
C
      IF (IFFREQ.EQ.0) CALL VIBI
C
   10 TSTATE = .FALSE.
      BIGELEISEN = .FALSE.
      IF (ICALC.EQ.1) THEN
          SLF1  = ZERO
          FTSL1 = ZERO
      ELSE
          SLF2  = ZERO
          FTSL2 = ZERO
      ENDIF
      DO 2 L=1,6
    2 IX(L)=ZERO
      READ(IR,103) ITS,NX,(IX(L),L=1,NX)
      IF (NX.NE.0) THEN
        WRITE(IW,207) NX
        WRITE(IW,206) (FREQ(IX(L)),L=1,NX)
      ENDIF
      IF(ITS.NE.0) THEN
         TSTATE = .TRUE.
         IF(LIBRATIONS) BIGELEISEN = .TRUE.
         RCF = -FREQ(ITS)
        IF (ICALC.EQ.1) THEN
           RCF1  = RCF
        ELSE
           RCF2  = RCF
        ENDIF
        WRITE(IW,204) RCF 
      ENDIF
      DO 3 I=1,NT
      ZP(I)=ZERO
    3 EX(I)=ZERO
    8 CONTINUE
      IF (LIBRATIONS) THEN
         WRITE(IW,210)
         IF (TSTATE) THEN
           WRITE(IW,211)
           IF (ICALC.EQ.1) THEN
              FTSL1 = DLOG(RCF)
           ELSE
              FTSL2 = DLOG(RCF)
           ENDIF
         ELSE
           WRITE(IW,212)
         ENDIF
      ELSE
         IF (TSTATE) THEN
           WRITE(IW,213)
           IF (ICALC.EQ.1) THEN
              FTSL1 = DLOG(RCF)
           ELSE
              FTSL2 = DLOG(RCF)
           ENDIF
         ELSE
           WRITE(IW,214)
         ENDIF
      ENDIF
C
C     Loop over frequencies
C
      DO 6  I=1,NFREQ
      SKIP = .FALSE.
      DO 4  L=1,NX
    4 IF(I.EQ.IX(L)) SKIP = .TRUE.
      IF(SKIP) GO TO 6
      FRI = FREQ(I)
      IF(FRI.LT.ZERO) FRI = -FRI
      IF (ICALC.EQ.1) SLF1 = SLF1 + DLOG(FRI)
      IF (ICALC.GT.1) SLF2 = SLF2 + DLOG(FRI)
C
C     Loop over temperatures
C
      DO 5  J=1,NT
      U=HCOK*FRI/T(J)
      UO2 = U/2.0D0
      IF(I.EQ.ITS) THEN
C
C     Quantum correction for transition vector (Bell tunnelling)
C
        GAMMA = UO2/DSIN(UO2)
        IF (ICALC.EQ.1) THEN
          QTUN(1,J) = GAMMA
        ELSE
          QTUN(2,J) = GAMMA
        ENDIF
      ELSE
C
C     ZP is the zero-point energy divided by kT
C     EXCG is the natural log of the “excitational” partition function
C     Note: the full vibrational partition function (counting from the 
C     classical potential energy minimum) is exp(-U/2)/(1-exp(-U))
C
        ZP(J)=ZP(J)+UO2
        XU=DEXP(-U)
        EX(J)=EX(J)-DLOG(ONE-XU)
      ENDIF
    5 CONTINUE
    6 CONTINUE
C
C     Calculate partition functions at each temperature
C
      DO 7  I=1,NT
      WRITE(IW,208) T(I)
      IF(LIBRATIONS) THEN
C
C     Translational and rotational degrees of freedom are treated as
C     external vibrations of the system, and so make no separate
C     contribution to the free energy.
C
        TR = ZERO
        RT = ZERO
      ELSE
C
C     TR is the natural log of the translational partition function
C
        TR = CTR  +  TWOP5*TLN(I)  +  ONEP5*DLOG(WM)
C
C     RT is the natural log of the rotational partition function
C
        RT = - SYMM
     *       + ZP5*DFLOAT(NR-2)*PILN
     *       + ZP5*DFLOAT(NR)*(CROT+TLN(I))
     *       + ZP5*PMLN
      ENDIF
C
C     G0E0 is the natural log of the partition function for the whole 
C          system, relative to the zero-point vibrational energy level,
C          at a constant pressure of 1 atmosphere.
C
C     GEL is the natural log of the electronic degeneracy
C
      G0E0=TR+RT+EX(I)+GEL
C
C     G0EP is the natural log of the partition function for the whole
C          system, relative to the classical potential energy minimum,
C          at a constant pressure of 1 atmosphere.
C
      G0EP=G0E0-ZP(I)
      WRITE(IW,209) G0EP,G0E0,ZP(I),TR,RT,EX(I),GEL
C
C     Now prepare quantities for the top and bottom of the IPF ratio
C
      IF (ICALC.EQ.1) THEN
        QT(1,I) = TR
        QR(1,I) = RT
        QE(1,I) = EX(I)
        QZ(1,I) = -ZP(I)
      ELSE
        QT(2,I) = TR
        QR(2,I) = RT
        QE(2,I) = EX(I)
        QZ(2,I) = -ZP(I)
      ENDIF
    7 CONTINUE
C
      RETURN
C
  102 FORMAT(6F12.4)
  103 FORMAT(20I4)
  203 FORMAT(/20('*'),' Calculations for isotopologue',I2,1X,
     * 20('*')/)
  204 FORMAT(/'Transition frequency =',F10.4,'i cm-1')
  205 FORMAT(/'Input frequencies:')
C 206 FORMAT(6F13.7)
  206 FORMAT(6F12.4)
  207 FORMAT(/'The',I2,' rejected frequencies are:')
  208 FORMAT(/'Temperature =',F8.2,' K, Std. State = 1 Atmosphere')
  209 FORMAT(/'-(G-PE)/RT=',F12.6,' = ','-(G-E0)/RT=',F11.6,
     *'  -(ZPE/RT=',F12.6,')',/10X,'TRANS',9X,'ROT',12X,'EXC',11X,
     *'ELEC'/5X,4(2X,F12.6))
  210 FORMAT(/'Translational and rotational degrees of freedom treated a
     *s librations:')
  211 FORMAT(1X,'EXC and ZPE evaluated over 3N-1 real TS frequencies')
  212 FORMAT(1X,'EXC and ZPE evaluated over 3N frequencies')
  213 FORMAT(1X,'EXC and ZPE evaluated over 3N-7 TS frequencies')
  214 FORMAT(1X,'EXC and ZPE evaluated over 3N-6 frequencies')
      END
C
      SUBROUTINE RATIO(ICALC,BIGELEISEN,TSTATE)
      INCLUDE 'CAM_SIZE'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL TSTATE,LIBRATIONS,BIGELEISEN
      COMMON /USE/IR,IW,NAT,NC,IFFREQ,NT,T(5),TLN(5),WM,TLM,SYM,GEL
      COMMON /MOL/F(MAXDIAG),FREQ(MAXNC),X(3,MAXNAT),W(MAXNAT)
      COMMON /BOTH/QT(2,5),QR(2,5),QE(2,5),QZ(2,5),QTUN(2,5),
     *SLM1,TMLN1,SLI1,SLF1,QT1,QR1,QE1,QZ1,FTSL1,RCF1,
     *SLM2,TMLN2,SLI2,SLF2,QT2,QR2,QE2,QZ2,FTSL2,RCF2
      DATA ZP5,ONEP5 /0.5D0,1.5D0/
C     Check for internal consistency of calculations:
C     the Teller-Redlich Product Rule requires TMMI & TF to be equal.
C
      TMMI = DEXP(ONEP5*(TMLN1-SLM1-TMLN2+SLM2)+ZP5*(SLI1-SLI2))
      TF   = DEXP(SLF1-SLF2)
      WRITE(IW,213) ICALC,TMMI,TF
C
C     CALCULATE CONTRIBUTIONS TO ISOTOPIC PARTITION FUNCTION RATIO 
C     AT EACH TEMPERATURE
C
      IF(.NOT. BIGELEISEN) THEN
        IF(.NOT. TSTATE) THEN
          DO 1  I=1,NT
          WRITE(IW,214) ICALC,T(I)
          RMMILN = (QT(2,I) + QR(2,I)) - (QT(1,I) + QR(1,I))
          RMMI   = DEXP(RMMILN)
          REXCLN = QE(2,I) - QE(1,I)
          REXC   = DEXP(REXCLN)
          RZPELN = QZ(2,I) - QZ(1,I)
          RZPE   = DEXP(RZPELN)
          PFRLN  = RMMILN + REXCLN + RZPELN
          PFR    = DEXP(PFRLN)
          WRITE(IW,215) ICALC,PFRLN,ICALC,PFR,
     *                 RMMILN,RMMI,REXCLN,REXC,RZPELN,RZPE
    1     CONTINUE
        ELSE
C
C     The ratio of isotopic transition frequencies is “light/heavy”
C
        RRCF   = RCF1/RCF2
        WRITE(IW,218) ICALC,RRCF
          DO 2  I=1,NT
          WRITE(IW,214) ICALC,T(I)
          RMMILN = (QT(2,I) + QR(2,I)) - (QT(1,I) + QR(1,I))
          RMMI   = DEXP(RMMILN)
          REXCLN = QE(2,I) - QE(1,I)
          REXC   = DEXP(REXCLN)
          RZPELN = QZ(2,I) - QZ(1,I)
          RZPE   = DEXP(RZPELN)
          PFRLN  = RMMILN + REXCLN + RZPELN
          PFR    = DEXP(PFRLN)
          WRITE(IW,215) ICALC,PFRLN,ICALC,PFR,
     *                 RMMILN,RMMI,REXCLN,REXC,RZPELN,RZPE
C
C     Quantum correction to transition-frequency factor
C
          QCRRCF   = QTUN(2,I)/QTUN(1,I)
          QPFR     = PFR*QCRRCF
          WRITE(IW,220) QCRRCF,QPFR
    2     CONTINUE
        ENDIF
      ELSE
        DO 3  I=1,NT
        WRITE(IW,214) ICALC,T(I)
        RRCFLN = FTSL2 - FTSL1
        RRCF   = DEXP(RRCF)
        REXCLN = QE(2,I) - QE(1,I)
        REXC   = DEXP(REXCLN)
        RZPELN = QZ(2,I) - QZ(1,I)
        RZPE   = DEXP(RZPELN)
        PFRLN  = RRCFLN + REXCLN + RZPELN
        PFR    = DEXP(PFRLN)
        WRITE(IW,216) ICALC,PFRLN,ICALC,PFR,
     *               RRCFLN,RRCF,REXCLN,REXC,RZPELN,RZPE
C
C     Quantum correction to transition-frequency factor
C
          QCRRCF   = QTUN(2,I)/QTUN(1,I)
          QPFR     = PFR*QCRRCF
          WRITE(IW,220) QCRRCF,QPFR
    3   CONTINUE
      ENDIF
      RETURN
C
  213 FORMAT(/'Teller-Redlich Product Rule for isotopologue',I3/3X,
     * 'Mass and Moment of Inertia term = ',F16.8/20X,
     * 'Frequency term = ',F16.8)
  214 FORMAT(/'Partition Function Ratio (IPFR) for isotopologues 1 and',
     * I2,' at ',F8.2,' K')
  215 FORMAT(/5X,'ln(Q',I1,'/Q1)  = ',F12.6,8X,'Q',I1,'/Q1  = ',
     * F15.6//5X,'ln(MMI)    = ',F12.6,8X,'MMI    = ',F15.6/
     * 5X,'ln(EXC)    = ',F12.6,8X,'EXC    = ',F15.6/
     * 5X,'ln(ZPE)    = ',F12.6,8X,'ZPE    = ',F15.6)
  216 FORMAT(/5X,'ln(Q',I1,'/Q1)  = ',F12.6,8X,'Q',I1,'/Q1  = ',
     * F15.6//5X,'ln(RRCF)   = ',F12.6,8X,'RRCF   = ',F15.6/
     * 5X,'ln(EXC)    = ',F12.6,8X,'EXC    = ',F15.6/
     * 5X,'ln(ZPE)    = ',F12.6,8X,'ZPE    = ',F15.6)
  218 FORMAT(/'Transition frequency ratio TF1/TF',I1,' =',4X,F13.8)
  220 FORMAT(/'Quantum correction to transition frequency ratio =',
     * F12.8/'Quantum corrected isotopic partition ratio =',2X,F16.8)
      END
C
      SUBROUTINE XMI(ICALC)
      INCLUDE 'CAM_SIZE'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*76 TITLE
      CHARACTER*4 ATOM(50),LAST,LSTA,IAT(1000),IA
C
C     Input isotopic masses and molecular geometry
C
      COMMON /USE/IR,IW,NAT,NC,IFFREQ,NT,T(5),TLN(5),WM,TLM,SYM,GEL
      COMMON /MOL/F(MAXDIAG),FREQ(MAXNC),X(3,MAXNAT),W(MAXNAT)
      COMMON /BOTH/QT(2,5),QR(2,5),QE(2,5),QZ(2,5),QTUN(2,5),
     *SLM1,TMLN1,SLI1,SLF1,QT1,QR1,QE1,QZ1,FTSL1,RCF1,
     *SLM2,TMLN2,SLI2,SLF2,QT2,QR2,QE2,QZ2,FTSL2,RCF2
      DIMENSION Z(50)
      DATA Z/1.007825D0,2.014102D0,3.016050D0,12.0D0,13.003354826D0,
     * 14.003242D0,14.003074002D0,15.00010897D0,15.99491463D0,
     * 16.9991312D0,17.991603D0,18.998405D0,31.97207070D0,33.96786665D0,
     * 34.96885D0,36.9659D0,6.01513D0,7.01601D0,9.01219D0,10.01294D0,
     * 11.009305D0,28.085611D0,30.97376D0,78.9184D0,80.9163D0,
     * 22.98977D0,999.0D0,24.309846D0,26.98153D0,35.458127D0,
     * 39.09778D0,40.004845D0,65.386887D0,79.906561D0,101.07D0,
     * 126.90435D0,132.9051D0,10.814017D0,1.00790D0,12.01100D0,
     * 14.00670D0,15.99940D0,32.060D0,4.002602D0,20.1797D0,39.948D0,
     * 83.798D0,32.97145843D0,35.96708062D0,11.0114336D0/
      DATA ATOM/'  1H','  2H','  3H',' 12C',' 13C',' 14C',' 14N',
     * ' 15N',' 16O',' 17O',' 18O','   F',' 32S',' 34S','35CL','37CL',
     * ' 6LI',' 7LI',' 9BE',' 10B',' 11B','  SI','   P','79BR','81BR',
     * '  NA','  SP','  MG','  AL','  CL','   K','  CA','  ZN','  BR',
     * '  RU','   I','  CS','   B','   H','   C','   N','   O','   S',
     * '  HE','  NE','  AR','  KR',' 33S',' 36S',' 11C'/
      DATA LAST,LSTA /'END ','*   '/
C
C     Symmetry number (IS), ground-state electronic degeneracy (IG), 
C     "Are frequencies to be read in?" (IIFREQ), and TITLE.
C
      READ(IR,100) IS,IG,IFFREQ,TITLE
      WRITE(IW,200) TITLE,IS,IG
      SYM=DLOG(DFLOAT(IS))
      GEL=DLOG(DFLOAT(IG))
      IF(ICALC.GT.1) GO TO 4
C
C     Data for isotopolog 1: symbol and cartesian coordinates in Ångström
C
      N=1
    1 IZ=0
      READ(IR,101) IA,(X(J,N),J=1,3)
      IF(IA.EQ.LAST.OR.IA.EQ.LSTA) GO TO 3
      DO 2 M=1,50
    2 IF(ATOM(M).EQ.IA) IZ=M
      IF(IZ.EQ.0) GO TO 99
      W(N)=Z(IZ)
      WRITE(IW,201) N,IA,W(N),(X(J,N),J=1,3)
      N=N+1
      GO TO 1
C
    3 NAT=N-1
      NC=NAT*3
      GO TO 7
C
C     Data for subsequent isotopologs: only isotopic symbols required
C     (Cartesian coordinates are retrieved from array X)
C
    4 READ(IR,102) (IAT(I),I=1,NAT)
      DO 6  I=1,NAT
      IZ=0
      IA=IAT(I)
      DO 5  M=1,50
    5 IF(ATOM(M).EQ.IA) IZ=M
      IF(IZ.EQ.0) GO TO 99
      W(I)=Z(IZ)
    6 WRITE(IW,201) I,IA,W(I),(X(J,I),J=1,3)
C
C     Calculate molecular mass, WM
C
    7 WM=ZERO
      TLM=ZERO
      DO 8  I=1,NAT
      TLM = TLM + DLOG(W(I))
    8 WM  = WM  +      W(I)
      RETURN
C
   99 WRITE(IW,209) IA
      STOP
C
  100 FORMAT(3I2,74A)
  101 FORMAT(A4,1X,3F12.7)
  102 FORMAT(16(A4,1X))
  200 FORMAT(1X,A//'Symmetry number = ',I2,6X,
     *'Electronic degeneracy = ',I2//4X,'Atom',6X,
     *'Mass',8X,'X',13X,'Y',13X,'Z',2X,'(Ångström)'/)
  201 FORMAT(I3,2X,A4,2X,F10.6,3(2X,F12.8))
  209 FORMAT(/'ATOMIC SPECIES ',A4,' IS NOT RECOGNISED')
      END
C
      SUBROUTINE MOMIN (NR,PM)
      INCLUDE 'CAM_SIZE'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /USE/IR,IW,NAT,NC,IFFREQ,NT,T(5),TLN(5),WM,TLM,SYM,GEL
      COMMON /MOL/F(MAXDIAG),FREQ(MAXNC),X(3,MAXNAT),W(MAXNAT)
      DIMENSION PM(3),CM(3),Q(1000),D(6),A(3,3),EIG(3)
      DATA ZERO,TM3/0.0D0,0.001D0/
	  WRITE(IW,*) NAT
      DO M=1,3
        CM(M) = ZERO
        DO   K=1,NAT
          CM(M)=CM(M) + W(K)*X(M,K)
        END DO
      CM(M)=CM(M)/WM
      END DO
      DO I=1,NAT
      Q(I)=ZERO
         DO M=1,3
            XM=X(M,I) - CM(M)
            Q(I) = Q(I) + XM*XM
         END DO
      END DO
      DO I=1,3
        DO J=1,I
          A(I,J)=ZERO
          DO K=1,NAT
            XIK=X(I,K) - CM(I)
            XJK=X(J,K) - CM(J)
            IF(I-J .EQ. 0) THEN
               A(I,J) = A(I,J) + W(K) * ( Q(K) - XIK*XIK )
            ELSE
               A(I,J) = A(I,J) - W(K) * XIK*XJK
            ENDIF
          END DO
        END DO
      END DO
      A(1,2) = A(2,1)
      A(1,3) = A(3,1)
      A(2,3) = A(3,2)
      CALL DSYEVC3(A,EIG) 
      DO I=1,3
        IF(EIG(I)-TM3 .GT. 0) THEN
          NR=NR+1
          PM(NR)=EIG(I)
        ENDIF
      END DO
      WRITE(IW,200) NR,(PM(I),I=1,NR)
      RETURN
  200 FORMAT(/'The',I2,' Principal Moments (Amu Å**2) are:'/3F20.10)
      END
C      
      SUBROUTINE FCI
      INCLUDE 'CAM_SIZE'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /USE/IR,IW,NAT,NC,IFFREQ,NT,T(5),TLN(5),WM,TLM,SYM,GEL
      COMMON /MOL/F(MAXDIAG),FREQ(MAXNC),X(3,MAXNAT),W(MAXNAT)
      NNC= NC*(NC+1)/2
      READ(IR,101) (F(I),I=1,NNC)
      RETURN
  101 FORMAT(1X,3E20.16)
C  101 FORMAT(1X,3E20.12)
C 101 FORMAT(6F12.8)
      END
C
      SUBROUTINE VIBI
      INCLUDE 'CAM_SIZE'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION V(MAXDIAG),G(MAXNC)
      COMMON /USE/IR,IW,NAT,NC,IFFREQ,NT,T(5),TLN(5),WM,TLM,SYM,GEL
      COMMON /MOL/F(MAXDIAG),FREQ(MAXNC),X(3,MAXNAT),W(MAXNAT)
      DATA TFACT/2.6436411815482D07/,ZERO/0.0D0/
      L=0
	  PRINT*, NAT
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
      WRITE(IW,200)
      WRITE(IW,201) (FREQ(I),I=1,NC)
      RETURN
  200 FORMAT(/'Frequencies of normal modes (cm-1)')
  201 FORMAT(6F13.7)
C  201 FORMAT(6F12.4)
      END
C
      SUBROUTINE EIGEN(N,A)
      INCLUDE 'CAM_SIZE'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /USE/IR,IW,NAT,NC,IFFREQ,NT,T(5),TLN(5),WM,TLM,SYM,GEL
      COMMON /MOL/F(MAXDIAG),E(MAXNC),X(3,MAXNAT),WT(MAXNAT)
      DIMENSION W(5,MAXNC),A(MAXDIAG)
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

      SUBROUTINE DSYEVC3(A,W)
C
C     Calculates the eigenvalues (W) of a symmetric 3x3 matrix A 
C     using Cardano's analytical algorithm.
C     Joachim Kopp, Numerical diagonalization of hermitian 3x3 matrices
C     arXiv.org preprint: physics/0610206
C
C     Arguments
      DOUBLE PRECISION A(3,3),W(3)
C     Parameters
      DOUBLE PRECISION SQRT3
      PARAMETER ( SQRT3 = 1.73205080756887729352744634151D0 )
C     Local Variables
      DOUBLE PRECISION M,C1,C0,DE,DD,EE,FF,P,SQRTP,Q,C,S,PHI
C
C     Determine coefficients of characteristic polynomial. We write
C           | A   D   F  |
C      A =  | D*  B   E  |
C           | F*  E*  C  |
C
      DE    = A(1,2) * A(2,3)
      DD    = A(1,2)**2
      EE    = A(2,3)**2
      FF    = A(1,3)**2
      M     = A(1,1) + A(2,2) + A(3,3)
      C1    = ( A(1,1)*A(2,2) + A(1,1)*A(3,3) + A(2,2)*A(3,3) )
     *         - (DD + EE + FF)
      C0    = A(3,3)*DD + A(1,1)*EE + A(2,2)*FF - A(1,1)*A(2,2)*A(3,3)
     *         - 2.0D0 * A(1,3)*DE
C
      P     = M**2 - 3.0D0 * C1
      Q     = M*(P - (3.0D0/2.0D0)*C1) - (27.0D0/2.0D0)*C0
      SQRTP = SQRT(ABS(P))
C
      PHI   = 27.0D0 * ( 0.25D0 * C1**2 * (P - C1)
     *          + C0 * (Q + (27.0D0/4.0D0)*C0) )
      PHI   = (1.0D0/3.0D0) * ATAN2(SQRT(ABS(PHI)), Q)
C
      C     = SQRTP * COS(PHI)
      S     = (1.0D0/SQRT3) * SQRTP * SIN(PHI)
C
      W(2) = (1.0D0/3.0D0) * (M - C)
      W(3) = W(2) + S
      W(1) = W(2) + C
      W(2) = W(2) - S
C
	  PRINT*,W(1),W(2),W(3)
      END


