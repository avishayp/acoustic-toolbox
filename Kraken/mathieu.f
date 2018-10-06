      SUBROUTINE MTHU (Q, IORD, IRA, IEO, INORM, MAXORD)

C     *********************************************************************
C     *                                                                   *
C     *     This set of routines is based on N. Toyama and K. Shogen,     *
C     *     "Computation of the value of the even and odd Mathieu         *
C     *     functions of order N for a given parameter S and an           *
C     *     argument X," IRE Transactions on Antennas and Propagation     *
C     *     Vol. Ap-32(5), 537-539 (1984).  However, my implementation    *
C     *     differs in several significant aspects, particularily in      *
C     *     calling method, the normalization options provided, the       *
C     *     calculation of the series coefficients, D, the use of         *
C     *     expansions for large and small Q, the method of calculating   *
C     *     the continued fraction poynomial for finding the roots and    *
C     *     other points were attempts have been made to improve          *
C     *     stability.                                                    *
C     *                               Stephen C. Wales, November 1984     *
C     *                                                                   *
C     *********************************************************************

C     This routine initializes the various parameters and calls the zero and 
C     series calculating subroutines if needed.  The arguments are:
C     Q the parameter (must be geater than 0)
C     IORD the order (must be geater than 0)
C     IRA = 1 Radial functions of the third type (= first + i x second)
C         = 2 Angular functions
C     IEO = Even-odd switch, even if even, odd if odd
C     INORM is the normalization switch.  The normalization may be specified
C           by type if negative or zero, or by author if positive.  See
C           subroutine DNORM for a detailed explanation of the possible options.
C     MAXORD causes the zeros to be located up to MAXORD the first time,
C         thereby saving computation time.  MAXORD will be updated if IORD
C         is greater, but this is a terrible waste of computation time since
C         finding the zeros is 90 % of the battle and the routine must find
C         each one up to the one you're using. Eg. If the user doesn't set
C         this parameter and works up to N zeros, the program will have to do
C         N*(N+1)/2 zero searches instead of N.

C     Two subroutines are called:
C     ABCALC to calculate the zeros of the continued fraction.
C     DCALC to calculate the series coefficients.

      IMPLICIT REAL*8 (A-G,O-Z)

      COMMON /MPARMS/ QP, IORDP, IRAP, IEOP, INORMP, MAXORDP
      COMMON /MCOEFS/ ICOEF, MAXZ, AB(0:200,4), QOLD(4), 
     1                                     D(0:300), ND, NDMAX, MZ(4)

C   * TEST FOR LEGAL VALUES *
                                                  ! LEGAL VALUES OF ORDER
      IF (IORD.LT.MOD(IEO,2) .OR. IORD.GT.300)  THEN
          WRITE (6,111)
111       FORMAT (' ILLEGAL ORDER FOR MATHIEU FUNCTION')
          IF (IORD.LT.0)  THEN
              WRITE (6,112) IORD
112       FORMAT (' ORDER IS LESS THAN ZERO, IORD =',I5)
          ELSE IF (IORD.EQ.0)  THEN
              WRITE (6,113)
113       FORMAT (' ORDER = 0 FOR AN ODD FUNCTION')
          ELSE IF (IORD.GT.300)  THEN
              WRITE (6,114)
114       FORMAT (' ORDER > 300 !')
          END IF
          STOP '* ERROR IN MATHIEU FUNCTION CALL - ORDER *'

      ELSE IF (Q.LT.0)  THEN                      ! LEGAL VALUES OF Q
          WRITE (6,115)
115       FORMAT (' ILLEGAL VALUE OF Q FOR MATHIEU FUNCTION')
          WRITE (6,116) Q
116       FORMAT (' Q = ',1PE17.10,' (< 0)')
          STOP '* ERROR IN MATHIEU FUNCTION CALL - Q *'

      ELSE IF (X.LT.0)  THEN                      ! LEGAL VALUES OF X
          WRITE (6,117)
117       FORMAT (' ILLEGAL VALUE OF X FOR MATHIEU FUNCTION')
          WRITE (6,118) X
118       FORMAT (' X = ',1PE17.10, '(< 0)')
          STOP '* ERROR IN MATHIEU FUNCTION CALL - X *'
      END IF

C   * LOAD PARAMETERS IN MPARMS COMMON *
C     (AN ADDITIONAL P IS ADDED TO PARAMETERS SO THEY MAY PASSED THROUGH
C      THROUGH THE MPARMS COMMON)
                                                  ! * MPARMS COMMON *
      QP = Q
      IORDP = IORD
      IRAP = IRA
      IEOP = MOD (IEO,2)
      INORMP = INORM
      MAXORDP = MAXORD
C                         ! *> use of maxord gives sysaccviol err <*
      IF (MAXORDP.LT.IORD)  MAXORDP = IORD 
                                                  ! * ABPARM COMMON *
      MAXZ = MAXORDP/2 + 1                        ! MAX. NUM. OF ZEROS TO FIND
      ICOEF = 2*IEOP + 1 + MOD(IORD,2)            ! 1 = EVEN FUNC, EVEN ORDER
                                                  ! 2 = EVEN FUNC, ODD ORDER
                                                  ! 3 = ODD FUNC, EVEN ORDER
                                                  ! 4 = ODD FUNC, ODD ORDER
      ITYPE = (1-2*IEOP) * IORD                   ! = IORD IF EVEN FUNCTION
                                                  ! = -IORD IF ODD FUNCTION

      IF (IFIRST.EQ.0)  THEN                      ! INITIALIZE OLD Q ARRAY
          DO 10 I = 1,4
   10     QOLD(I) = -1.D0
          QDOLD = -1.D0
          IFIRST = 1
      END IF

C   * * * CALCULATE ZEROS AND THE COEFFICIENTS FOR THE SERIES EXPANSION * * *

C     (RECALCULATE ZEROS AND COEFFICIENTS IF Q HAS CHANGED)
C     (RECALCULATE COEFFICIENTS IF ZERO, A/B, BEING USED IS NO LONGER
C      THE SAME FOR ANY REASON, DIFF. Q, CHANGE OF ORDER OR FUNCTION TYPE)

      IF (ABS(Q-QOLD(ICOEF)).GT.1.D-10 .OR. MAXZ.GT.MZ(ICOEF))  THEN
          CALL ABCALC
      END IF

      IF (ABS(Q-QDOLD).GT.1.D-10 .OR. ITYPE.NE.ITOLD)  THEN
          CALL DCALC
          QDOLD = Q
          ITOLD = ITYPE
      END IF
      RETURN
      END


C**********************************************************************C
      FUNCTION RMTHU( X, Q, IORD, IEO, IDER, INORM, MAXORD )

C     This is essentially a dummy function for convenience purposes.
C     Most of the parameters are defined in Mathieu above, the ones that
C     aren't are X the argument and Q the parameter of the Matheiu function.
C     and IDER the derivative to be calculated, eg. IDER = 0, no derivatve,
C     IDER = 1, first derivative, etc.  In the acoustic case,
C     Q = K**2 (A**2-B**2) / 4, and X = LN ((A+B)/(A-B)) / 2, where A
C     and B are the semi-major and semi-minor axes of the ellipse and K the
C     acoustic wavenumber, 2*pi/wavelength.  The quantity 2 sqrt(Q) cosh(X)
C     is merely KA.

      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 RMTHU, RMTHUV

      U = SQRT(Q) * EXP ( X)
      V = SQRT(Q) * EXP (-X)
      RMTHU = RMTHUV (U, V, IORD, IEO, IDER, INORM, MAXORD)
      RETURN
      END



      FUNCTION RMTHUV( U, V, IORD, IEO, IDER, INORM, MAXORD )

C     Calculate the radial Mathieu functions. All the parameters are
C     described in MTHU or RMTHU above except U and V.  U and V
C     are an optional set of parameters/arguments to be used instead of
C     X and Q normally used, they are related to X and Q as defined in
C     RMTHU above.  In the acoustic case, U = K * (A+B) / 2 and
C     V = K * (A-B) / 2.  Conversely Q = U * V and X = LN (U/V) / 2.
C     The last formula presents the poblem with the X,Q formulation
C     since X approaches infinity as A approaches B for the circular
C     cylinder.

      IMPLICIT REAL*8 (A-G,O-Z)
      IMPLICIT COMPLEX*16 (H)
      COMPLEX*16 AI
      PARAMETER (PI=3.1415926535897933D0, PIO2=PI/2.D0,
     &           AI=(0.0D0,1.D0) )
      REAL*8 BJU(0:500), BJV(0:500), BYU(0:500), D1BJV(0:500)

      COMPLEX*16 RMTHUV, HJU(0:500), D1HJU(0:500), D1HJU1, DRMTHU

      COMMON /MCOEFS/ ICOEF, MAXZ, AB(0:200,4), QOLD(4), 
     1                                     D(0:300), ND, NDMAX, MZ(4)

C   * INTIALIZE FOR CALCULATIONS *

      Q = U * V
      CALL MTHU (Q, IORD, 1, IEO, INORM, MAXORD)

C   * GET BESSEL FUNCTIONS, THESE CALLS ARE TO SLATEC LIBRARY *

      NDX = ND + NDMAX + 3                        ! MAX. ORDER OF BESSEL FUNC.

      IF (V.NE.VOLD .OR. NDX.GT.NDVOLD)  THEN     ! BESSEL FUNCS. FOR V TERMS
          IF (V.NE.VOLD)  THEN
              CALL DBESJ (V, 0.D0, NDX, BJV(0), NZJ)
          ELSE
              NS = NDVOLD                         ! LAST POINT IS RECALCULATED
              NT = NDX - NDVOLD + 1
              ALPHA = DBLE (NS)
              CALL DBESJ (V, ALPHA, NT, BJV(NS), NZJ)
          END IF
          D1BJV(0) = V * BJV(1)                   ! CALC. FIRST DERIVATIVE
          DO 5 I = 1,NDX-1
    5     D1BJV(I) = -V * (BJV(I-1) - BJV(I+1)) / 2.D0
          VOLD = V
          NDVOLD = NDX
      END IF

      IF (U.NE.UOLD .OR. NDX.GT.NDUOLD)  THEN     ! BESSEL FUNCS. FOR U TERMS
          IF (U.NE.UOLD)  THEN
              CALL DBESJ (U, 0.D0, NDX, BJU(0), NZJ)
              CALL DBESY (U, 0.D0, NDX, BYU(0))
          ELSE
              NS = NDUOLD                         ! LAST POINT IS RECALCULATED
              NT = NDX - NDUOLD + 1
              ALPHA = DBLE (NS)
              CALL DBESJ (U, ALPHA, NT, BJU(NS), NZJ)
              CALL DBESY (U, ALPHA, NT, BYU(NS))
          END IF
          DO 10 I = 0,NDX                         ! FORM HANKEL FUNCTION
   10     HJU(I) = DCMPLX (BJU(I), BYU(I))
          D1HJU(0) = -U * HJU(1)                  ! CALC. FIRST DERIVATIVE
          DO 15 I = 1,NDX-1
   15     D1HJU(I) = U * (HJU(I-1) - HJU(I+1)) / 2.D0
CD         WRITE (6,119) NT, NS
C119       FORMAT (I5,' BESSEL FUNCTIONS CALCULATED BEGINING AT',I5)
CD         IF (NDZ.GT.0)  WRITE (6,1110) NDZ
C1110        FORMAT ('   LAST',I5,
C     1                   ' ELEMENTS OF BESSEL FUNCTION WERE ZERO')
          UOLD = U
          NDUOLD = NDX
          DO 20 II = 0,NDX-2                      ! WRONSKIAN TEST
          WTEST = DBLE(HJU(II)) * DIMAG(D1HJU(II))
     1                              - DIMAG(HJU(II)) * DBLE(D1HJU(II))
          WTEST = WTEST*PIO2 - 1.D0
          IF (ABS(WTEST).GT.1.D-8 .AND. IHERCT.LT.20)  THEN
              IHERCT = IHERCT + 1
              WRITE (6,1111) II, WTEST, HJU(II),DIMAG(D1HJU(II)),
     1                       DIMAG(HJU(II)), D1HJU(II)
1111          FORMAT (' ERROR IN HANKEL FUNCTIONS, ',
     1                                 'WTEST(',I5,') =',1PE17.10,/
     2            ' = ', E17.10,' * ',E17.10,' - ',E17.10,' * ',E17.10)
          END IF
   20     CONTINUE
      END IF


C   * * * COMPLEX RADIAL FUNCTIONS * * *

      RMTHUV = (0.0D0, 0.0D0)
      MIEO2 = MOD (IEO,2)
      DO 100 K = 0,ND
      K1 = K - NDMAX
      K2 = K + NDMAX + MOD(IORD,2)
      IF (K1.LT.0)  THEN
          K1 = IABS (K1)
          BJV1 = (-1.D0)**K1 * BJV(K1)
          HJU1 = (-1.D0)**K1 * HJU(K1)
      ELSE
          BJV1 = BJV(K1)
          HJU1 = HJU(K1)
      END IF
      IF (MIEO2.EQ.1)  HJU1 = - HJU1
      HFACT = (-1.D0)**K * D(K) * (BJV1*HJU(K2) + BJV(K2)*HJU1)
  100 RMTHUV = RMTHUV + HFACT

C   * FIRST DERIVATIVE OF THE COMPLEX RADIAL FUNCTIONS *

      DRMTHU = (0.0D0, 0.0D0)
      DO 120 K = 0,ND
      K1 = K - NDMAX
      K2 = K + NDMAX + MOD(IORD,2)
      IF (K1.LT.0)  THEN
          K1 = IABS (K1)
          BJV1 = (-1.D0)**K1 * BJV(K1)
          HJU1 = (-1.D0)**K1 * HJU(K1)
          D1BJV1 = (-1.D0)**K1 * D1BJV(K1)
          D1HJU1 = (-1.D0)**K1 * D1HJU(K1)
      ELSE
          BJV1 = BJV(K1)
          HJU1 = HJU(K1)
          D1BJV1 = D1BJV(K1)
          D1HJU1 = D1HJU(K1)
      END IF
      IF (MIEO2.EQ.1)  THEN
          HJU1 = -HJU1
          D1HJU1 = -D1HJU1
      END IF
  120 DRMTHU = DRMTHU + (-1)**K * D(K) *
     1                       (BJV1*D1HJU(K2) + D1BJV1*HJU(K2)
     2                               + D1BJV(K2)*HJU1 + BJV(K2)*D1HJU1)

C   * SECOND DERIVATIVE OF THE COMPLEX RADIAL FUNCTIONS *

      IF (IDER.GE.2)  THEN
        WRITE (6,1114)
1114    FORMAT (' **************************************************',/
     1          ' *  SECOND AND HIGHER DERIVATIVES OF THE COMPLEX  *',/
     2          ' * RADIAL FUNCTIONS HAVE NOT YET BEEN IMPLEMENTED *',/
     3          ' **************************************************')
          STOP 'ERROR IN RMTHUV FUNCTION ROUTINE'
      END IF

C   * FINAL NORMALIZATION AND SIGN *

      IF (NDMAX+MOD(IORD,2).EQ.0)  THEN
          RMTHUV = RMTHUV / 2.D0
          DRMTHU = DRMTHU / 2.D0
      END IF
      RMTHUV = (-1.D0)**(IORD/2) * RMTHUV
      DRMTHU = (-1.D0)**(IORD/2) * DRMTHU
      CALL DNORM (RMTHUV)
      CALL DNORM (DRMTHU)

C   * CHECK ANSWER FOR VALIDITY BY FORMING WRONSKIAN

      IF (INORM.NE.3)  THEN
          WTEST = DBLE(RMTHUV) * DIMAG(DRMTHU)
     1                   - DIMAG(RMTHUV) * DBLE(DRMTHU)
          WTEST = WTEST * PIO2 - 1.D0
          IF (ABS(WTEST).GT.1.D-8 .AND. IMERCT.LT.20) THEN
              IMERCT = IMERCT + 1
              WRITE (6,1118) WTEST, RMTHUV, DRMTHUV, IORD,
     1                       IEO, Q, U, V
1118          FORMAT (' RMTHUV ERROR - ERROR IN WRONSKIAN ',
     1                                   'CALCULATION IS',1PE17.10,/
     2                ' RM, DRM:', 1PE17.10,', ', E17.10,/
     3                ' IORD, IEO, Q:', 0PI5, ', ', I5, 1PE17.10,/
     4                ' U, V:', E17.10, ', ',E17.10)
          END IF
          IF (ABS(WTEST).GT.1.D-2)  STOP 'RMTHUV FUNCTION ERROR'
      END IF

C   * IF DERIVATIVE REQUIRED SET RMTHUV TO DRMTHU *

      IF (IDER.EQ.1)  RMTHUV = DRMTHU

      RETURN
      END


C**********************************************************************C
      FUNCTION AMTHU( X, Q, IORD, IEO, IDER, INORM, MAXORD )

C     All of the parameters have been expalined elsewhere.  The argument
C     X is now interpreted as the angle in radians.

      IMPLICIT REAL*8 (A-G,O-Z)
      COMPLEX*16 FMTHU

      COMMON /MCOEFS/ ICOEF, MAXZ, AB(0:200,4), QOLD(4), 
     1                                     D(0:300), ND, NDMAX, MZ(4)

      CALL MTHU (Q, IORD, 2, IEO, INORM, MAXORD)
      AMTHU = 0.0D0

C   * ANGULAR FUNCTIONS *

      IF (IDER.EQ.0)  THEN
          IF (ICOEF.EQ.1)  THEN                   ! EVEN FUNC, EVEN ORDER
              DO 210 K = 0,ND
  210         AMTHU = AMTHU + D(K) * COS ((2*K)*X)

          ELSE IF (ICOEF.EQ.2)  THEN              ! EVEN FUNC, ODD ORDER
              DO 220 K = 0,ND
  220         AMTHU = AMTHU + D(K) * COS ((2*K+1)*X)

          ELSE IF (ICOEF.EQ.3)  THEN              ! ODD FUNC, EVEN ORDER
              DO 230 K = 1,ND
  230         AMTHU = AMTHU + D(K) * SIN ((2*K)*X)

          ELSE IF (ICOEF.EQ.4)  THEN              ! ODD FUNC, ODD ORDER
              DO 240 K = 0,ND
  240         AMTHU = AMTHU + D(K) * SIN ((2*K+1)*X)
          END IF

C   * FIRST DERIVATIVE OF THE ANGULAR FUNCTIONS *

      ELSE IF (IDER.EQ.1)  THEN
          IF (ICOEF.EQ.1)  THEN                   ! EVEN FUNC, EVEN ORDER
              DO 410 K = 0,ND
  410         AMTHU = AMTHU - (2*K) * D(K) * SIN ((2*K)*X)

          ELSE IF (ICOEF.EQ.2)  THEN              ! EVEN FUNC, ODD ORDER
              DO 420 K = 0,ND
  420         AMTHU = AMTHU - (2*K+1) * D(K) * SIN ((2*K+1)*X)

          ELSE IF (ICOEF.EQ.3)  THEN              ! ODD FUNC, EVEN ORDER
              DO 430 K = 1,ND
  430         AMTHU = AMTHU + (2*K) * D(K) * COS ((2*K)*X)

          ELSE IF (ICOEF.EQ.4)  THEN              ! ODD FUNC, ODD ORDER
              DO 440 K = 0,ND
  440         AMTHU = AMTHU + (2*K+1) * D(K) * COS ((2*K+1)*X)
          END IF

C   * SECOND DERIVATIVE OF THE ANGULAR FUNCTIONS *

      ELSE IF (IDER.EQ.2)  THEN
          IF (ICOEF.EQ.1)  THEN                   ! EVEN FUNC, EVEN ORDER
              DO 610 K = 0,ND
  610         AMTHU = AMTHU - (2*K)**2 * D(K) * COS ((2*K)*X)

          ELSE IF (ICOEF.EQ.2)  THEN              ! EVEN FUNC, ODD ORDER
              DO 620 K = 0,ND
  620         AMTHU = AMTHU - (2*K+1)**2 * D(K) * COS ((2*K+1)*X)

          ELSE IF (ICOEF.EQ.3)  THEN              ! ODD FUNC, EVEN ORDER
              DO 630 K = 1,ND
  630         AMTHU = AMTHU - (2*K)**2 * D(K) * SIN ((2*K)*X)

          ELSE IF (ICOEF.EQ.4)  THEN              ! ODD FUNC, ODD ORDER
              DO 640 K = 0,ND
  640         AMTHU = AMTHU - (2*K+1)**2 * D(K) * SIN ((2*K+1)*X)
          END IF

C   * THIRD DERIVATIVE OF THE ANGULAR FUNCTIONS *

      ELSE IF (IDER.GE.3)  THEN
          WRITE (6,1122)
1122      FORMAT (' *********************************************',/
     1            ' * THIRD & HIGHER DERIVATIVES OF THE ANGULAR *',/
     2            ' * FUNCTIONS HAVE NOT YET BEEN IMPLEMENTED.  *',/
     3            ' *********************************************')
          STOP 'ERROR IN AMTHU FUNCTION ROUTINE'
      END IF

C   * DO FINAL NORMALIZATION *

      FMTHU = (1.D0, 0.0D0)
      CALL DNORM (FMTHU)
      AMTHU = AMTHU * FMTHU
      RETURN
      END


C**********************************************************************C
      SUBROUTINE DCALC

C     This routine calculates the coefficients for the expansions of the
C     Mathieu functions.  The method used differs substantially from Toyama 
C     and Shogen.  The continued fraction representations for the ratios of 
C     the amplitudes of the coefficients found in Abramowitz and Stegun,
C     National Bureau of Standards, AMS-55, p. 723, 20.2.5 to 20.2.20.
C     are used.  Both the up and down representaitons are used.  Beginning 
C     in the center with a value of 1, the down representation, GDN, 
C     eq. 20.2.20, is used to obtain the lower half, and the up representation,
C     eq. 20.2.19,  is used to calculate the upper half.
C     THE RESULTING COEFFICIENTS ARE NOT NORMALIZED, THIS IS DONE BY DNORM.

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (PI=3.141592653589793D0, PRECR=1.D-16)
      REAL*8  GUP(400), GDN(300)

      COMMON /MPARMS/ Q, IORD, IRA, IEO, INORM, MAXORD      
      COMMON /MCOEFS/ ICOEF, MAXZ, AB(0:200,4), QOLD(4), 
     1                                     D(0:300), ND, NDMAX, MZ(4)

      IR = IORD / 2
      A = AB(IR,ICOEF)
CD     WRITE (6,1) A, IR, ICOEF
CD   1 FORMAT (' CHAR.VALUE =',1PE17.10,'(IR=',0PI5,', ICOEF=',I5,')')


C   * * * SPECIAL CASE FOR Q VERY SMALL, < 1.D-8 * * *

      IF (Q.LT.1.D-8)  THEN
          ND = IR + 4
          NDMAX = IR
          D(IR) = 1.D0
          Q4 = Q / 4.D0
          J = 0
          DO 5 I = IR-1,0,-1                      ! BELOW D(IR)
          J = J + 1
    5     D(I) = (Q4 / DBLE(J*(IORD-1))) * D(I+1)
          J = 0
          DO 7 I = IR+1,IR+4                      ! ABOVE D(IR)
          J = J + 1
    7     D(I) = -(Q4 / DBLE(J*(IORD+J))) * D(I-1)
          RETURN
      END IF


C   * * * NORMAL ROUTINE FOR CALCULATING THE SERIES COEFICIENTS D * * *

C   * SET UP CONSTANTS *

      IF (ICOEF.EQ.1)  THEN
          PHI0 =  2.0D0
          PHI1 =  0.0D0
      ELSE IF (ICOEF.EQ.2)  THEN
          PHI0 =  1.0D0
          PHI1 = -1.0D0
      ELSE IF (ICOEF.EQ.3)  THEN
          PHI0 =  0.0D0
          PHI1 =  1.0D37
      ELSE IF (ICOEF.EQ.4)  THEN
          PHI0 =  1.0D0
          PHI1 =  1.0D0
      END IF
      DAS = MOD (ICOEF-1,2)

C   * CALCULATE REASONABLE STARTING AND STOPPING POINTS *

      IF (A.GT.201.D0)  THEN                      ! NEAREST N**2
          M = NINT((SQRT(A)-DAS)/2.D0)
      ELSE
          M = 10
      END IF
      M2 = M**2
      ACC = 0.0D0
      
      DO 10 I = M,300
      ACC = ACC + LOG10(((2*I)**2-M2)/Q)
      IF (ACC.LE.35.D0) MSTOP = I                 ! STOPPING PT. FOR SERIES
      IF (ACC.GT.75.D0)  THEN
          MCF = I                                 ! ST.PT. FOR CONT. FRAC. CALC.
          GO TO 20
      END IF
   10 CONTINUE
      MCF = 300
 
C   * CALCULATE G'S *
C     COMMENT: THIS FORMULATION AND ALL OTHERS GETS INTO TROUBLE WHEN
C              A = (2*K+DAS)**2

   20 GDN(1) = (A - DAS**2)/Q + PHI1              ! RUNNING UP
      V = (A - (2.D0+DAS)**2) / Q
      GDN(2) = V - PHI0/GDN(1)

      DO 30 K = 3,MSTOP
      V = (A - (2.D0*DBLE(K-1)+DAS)**2) / Q
      GDN(K) = V - 1.D0/GDN(K-1)
   30 CONTINUE

      IF (ICOEF.EQ.3)  THEN
          KSIGN = 1
          KSTAR = 1
      ELSE
          KSIGN = 0
          KSTAR = 0
      END IF
      GUP(MCF+1) = 0.0D0                          ! RUNNING DOWN
      GMINER = 1.D0
      DO 40 K = MCF,2,-1
      VQ = A - (2.D0*DBLE(K)+DAS)**2
      GUP(K) = Q / (VQ - Q*GUP(K+1))
      IF (GUP(K).GT.0.0D0 .AND. KSIGN.EQ.KSTAR)  KSIGN = K
      GREL = ABS (GUP(K)-GDN(K)) / (ABS (GUP(K)) + PRECR)
      IF (GREL.LT.GMINER)  THEN
          GMINER = GREL
          KBEST = K
      END IF
   40 CONTINUE
      VQ = A - (2.D0+DAS)**2
      IF (ICOEF.NE.3) THEN
          GUP(1) = Q * PHI0 / (VQ - Q*GUP(2))
          IF (GUP(1).GT.0.0D0 .AND. KSIGN.EQ.KSTAR)  KSIGN = 1
          GREL = ABS (GUP(1)-GDN(1)) / (ABS (GUP(1)) + PRECR)
          IF (GREL.LT.GMINER)  THEN
              GMINER = GREL
              KBEST = 1
          END IF
      END IF

C   * CALCULATE UNNORMALIZED D'S AND LOCATE MAXIMUM VALUE *

      D(KBEST) = 1.0D0
      DO 50 K = KBEST,1,-1                        ! * LOW END *
      K1 = K - 1
   50 D(K1) = D(K) / GDN(K)
      DO 60 K = KBEST+1,MCF                       ! * HIGH END *
      D(K) = GUP(K) * D(K-1)
      IF (ABS(D(K)).LT.1.D-25)  THEN
          ND = K
          GO TO 100
      END IF
   60 CONTINUE
     
  100 WM = 0.D0
      D2 = 0.0D0
      DO 110 K = 0,ND
      WM = WM + K * D(K)**2
  110 D2 = D2 + D(K)**2
      WM = WM / D2
      NDMAX = NINT (WM)
CD     WRITE (6,120) KSIGN, KBEST, GMINER, NDMAX, D(NDMAX)
CD 120 FORMAT (' KSIGN =', I5,',  KBEST =', I5,/
CD    1        ' GMINER (QUALITY OF FIT) =',1PE17.10,/
CD    2        ' NDMAX, D:', 0PI5,', ', 1PE17.10 )
CD     WRITE (6,130) 0, REAL( D(0)/D(NDMAX) )
CD 130 FORMAT (' K, D          : ', I5,', ', 1PE17.10 )
CD     DO 140 K = 1,ND
CD 140 WRITE (6,150) K, REAL(D(K)/D(NDMAX)), REAL(GDN(K)), REAL(GUP(K))
CD 150 FORMAT (' K, D, GDN, GUP: ', I5, 3(', ',1PE17.10) )

C   * TEST FOR ACCEPTABLE MATCH BETWEEN THE TWO SERIES *

      IF (GMINER.GT.1.D-6)  THEN 
          WRITE (6,160) KBEST, GMINER
  160     FORMAT (' UNABLE TO OBTAIN GOOD FIT BETWEEN UP & DOWN G''S',
     1           ' IN DCALC',/
     2            ' KBEST =', I5,',    GMINER =', 1PE17.10 )
          STOP 'STOP - ERROR IN DCALC, BAD FIT'
      END IF
      RETURN
      END


C**********************************************************************C
      SUBROUTINE DNORM (FMTHU)

C     This routine normalizes the result to the requested normalization.
C     FMTHU is the unnormalized Mathieu function.
C     Several normalizations are possible as specified by the INORM parameter.
C     The JNORM parameter is the effective parameter internal to this 
C     subroutine.  If INORM is negative or zero it is set equal to JNORM,  
C     if positive then a lookup table is provided in JNORMS to convert INORM
C     to the proper representation.   The JNORMS table is as follows:
C                  JNORM   JNORM
C                 Radial  Angular      Authors
C     INORM =  1    -3      -2   Abramowitz and Stegun
C           =  2    -3      -1   Morse and Feshbach, Toyama and Shogen
C           =  3    -2      -2   McLachlin (Radial not implemented)
C           =  4    -3      -2   Abramowitz and Stegun
C           =  5    -3      -2   Normalized for plane wave expansions,
C                                A&S * sqrt(2), M&F * sqrt(2/M)
C           In addition certain of these need factors like SQRT (PI/2) etc.,
C           see SPECIAL CASES below.
C                 Normalized by ...
C     JNORM = -4, principle coefficient, IORD/2, 
C           = -3, the central coefficient, NDMAX = sum (n*Dn**2) / sum Dn**2
C           = -2, sum squared equal to one,
C           = -1, sum equal to one for even functions, unity value at X = 0,
C                 sum of order of coeff. times coeff., (2*K+IP)*D(K), equal to
C                 one for odd functions, unity slope at X = 0,
C           = 0,  lead coefficient equal to one

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (PI=3.141592653589793D0)

      INTEGER JNORMS(2,5)
      COMPLEX*16  FMTHU

      COMMON /MPARMS/ Q, IORD, IRA, IEO, INORM, MAXORD      
      COMMON /MCOEFS/ ICOEF, MAXZ, AB(0:200,4), QOLD(4), 
     1                                     D(0:300), ND, NDMAX, MZ(4)

      DATA JNORMS/-3,-2,-3,-1,-2,-2,-3,-2,-3,-2/

      IR = IORD / 2

C   * DETERMINE TYPE OF NORMALIZATION TO BE USED *

      IF (INORM.LE.0)  THEN                       ! NORMALIZATION SPECIFIED
          JNORM = INORM
      ELSE IF (INORM.LE.5)  THEN                  ! TABLE LOOKUP FOR
          JNORM = JNORMS(IRA,INORM)               !      VARIOUS AUTHORS
      ELSE
          JNORM = 0
      END IF
CD     WRITE (6,1136) INORM, JNORM
C1136   FORMAT (' INORM, JNORM:', I5, I5 )

C   * CALCULATE NORMALIZATION FACTOR *

      IF (JNORM.EQ.-4)  THEN                      ! * NORM. BY PRINC. COEFF.*
          DNRML = D(IR)

      ELSE IF (JNORM.EQ.-3)  THEN                 ! * NORM. BY LARGEST COEFF.*
          DNRML = D(NDMAX)

      ELSE IF (JNORM.EQ.-2)  THEN                 ! * NORM. BY SUM SQUARED *
          IF (ICOEF.EQ.1)  THEN
              DNRML = D(0)**2
          ELSE
              DNRML = 0.0D0
          END IF
          DO 110 K = 0,ND
  110     DNRML = DNRML + D(K)**2
          DNRML = SQRT (DNRML)
          IF (ICOEF.EQ.3)  THEN                   ! CORRECT FOR POSS. SIGN ERR.
              DNRML = DNRML * SIGN (1.D0,D(1))
          ELSE
              DNRML = DNRML * SIGN (1.D0,D(0))
          END IF

      ELSE IF (JNORM.EQ.-1)  THEN                 ! * NORMALIZE TO UNITY SUM *
          DNRML = 0.0D0
          IF (IEO.EQ.0)  THEN                     ! EVEN SOLUTIONS
              DO 120 K = 0,ND
  120         DNRML = DNRML + D(K)
          ELSE                                    ! ODD SOLUTIONS
              IP = MOD (IORD,2)
              DO 130 K = 0,ND
  130         DNRML = DNRML + DBLE (2*K+IP) * D(K)
          END IF

      ELSE                                        ! * NORMALIZE BY LEAD COEFF.*
          IF (ICOEF.EQ.3)  THEN
              DNRML = D(1)
          ELSE
              DNRML = D(0)
          END IF
      END IF

      IF (ABS(DNRML).LE.1.D-35)  THEN
          WRITE (6,1137) DNRML, INORM, JNORM
1137      FORMAT (' *** WARNING - NORMALIZATION ERROR ***',/
     1            '     DNRML =', 1PE17.10,/
     2            '     INORM =', E17.10,/
     3            '     JNORM =', E17.10,/
     4            '     USING DNRML = 1.D0')
          DNRML = 1.0D0
          WRITE (6,1142) IORD, NDMAX, REAL(D(NDMAX))
1142      FORMAT ('  IORD, NDMAX, D(NDMAX):', 2(I5,', '), 1PE17.10 )
      END IF

C   * SPECIAL CASES *

      IF (INORM.EQ.2)  THEN                       ! MORSE & FESH., TOY.& SHOGEN
          IF (IRA.EQ.1)  DNRML = SQRT (2.D0/PI) * DNRML        
      ELSE IF (INORM.EQ.3)  THEN                  ! MCLACHLIN
          IF (IRA.EQ.1)  STOP ' MCLACHLIN RADIAL NOT IMPLEMENTED'
      ELSE IF (INORM.EQ.5)  THEN
          IF (IRA.EQ.2)  DNRML = DNRML / SQRT(2.D0)
      END IF
CD     WRITE (6,1143) DNRML
C1143   FORMAT (' DNRML =', 1PE17.10 )

C   * DO NORMALIZATION *

      FMTHU = FMTHU / DNRML

      RETURN
      END


C**********************************************************************C
      SUBROUTINE ABCALC

C     This routine finds the Nth zero of the continued fraction defining the
C     characteristic values for the parameter Q.  It is similar to Toyama and
C     Shogen in that it uses a step search to find all of the zeros up to the 
C     one desired, but that is about it.  It calculates the first zero using
C     a call to AINIT and then improves on it by searching if needed.  Other
C     zeros are found by stepping through at a step size that depends on the
C     zero being located.  Large and small Q approximations are used wherever
C     possible and improved on by seach if needed.

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (PRECR=1.D-16, PRECA=1.D-32)

      COMMON /MPARMS/ Q, IORD, IRA, IEO, INORM, MAXORD
      COMMON /MCOEFS/ ICOEF, MAXZ, AB(0:200,4), QOLD(4), 
     1                                     D(0:300), ND, NDMAX, MZ(4)

      EXTERNAL ABFUNC

C   * FIND INITIAL VALUE OF A/B *

      AR = AINIT (Q, ICOEF)
      AZ = ABFUNC (AR, Q, ICOEF)
CD     WRITE (6,10) AR, AZ
CD  10 FORMAT (' INITIAL A/B ESTIMATE & ACCURACY: ',1PE17.10,',',E18.10)
      IF (ABS(AZ).GT.PRECR*ABS(AR)+PRECA)  THEN
          XDEL = 0.1D0
          CALL SEARCH (ABFUNC, AR, Q, ICOEF, XDEL, 1)
      END IF

C   * PLACE IN PROPER POSITION & SET INITIAL LOOP PARAMETER *

      IF (ICOEF.NE.3)  THEN
          AB(0,ICOEF) = AR
CD         WRITE (6,20) AR
CD  20     FORMAT (' AB(0,ICOEF) = ', 1PE17.10)
          IR1 = 1
      ELSE
          AB(0,ICOEF) = 0.0D0
          AB(1,ICOEF) = AR
CD         WRITE (6,30) 0.0D0, AR
CD  30     FORMAT (' AB(0,ICOEF) = ', 1PE17.10,/,
CD    1            ' AB(1,ICOEF) = ', 1PE17.10)
          IR1 = 2
      END IF

C   * SET EVEN-ODD SWITCH PARAMETER *

      IC1 = MOD (ICOEF-1,2)                       ! 0 = EVEN, 1 = ODD

C   * LOOP OVER K FINDING THE FIRST MAXZ ZEROS OF A (OR B) *

CD     WRITE (6,40) IR1, MAXZ
CD  40 FORMAT (' IR LOOP FROM', I5,' TO', I5 )
      DO 100 IR = IR1,MAXZ
      IOR = 2*IR + IC1
      XDEL = (IOR**2 - (IOR-2)**2) / 2.D0

      IF (IOR.LE.7)  THEN                         ! DETERMINE LIMIT FOR SMALL Q
          QTEST = 5.D-2 + 8.D-5 * IOR**4
      ELSE
          QTEST = 10.D0** ((-15.D0+14.D0*LOG10(DBLE(IOR)))/8.D0)
      END IF
CD     WRITE (6,50) QTEST, IOR
CD  50 FORMAT (' QTEST, IOR:', 1PE17.10,', ',0PI5 )

      IF (Q.GT.10.D0*QTEST)  THEN                 ! MEDIUM AND LARGE Q'S
CD         WRITE (6,60)
CD  60     FORMAT (' A FOUND BY SEARCH')
          AR = AR + XDEL
          CALL SEARCH (ABFUNC, AR, Q, ICOEF, XDEL, 0)
      ELSE                                        ! SMALL Q'S
          IF (ICOEF.LE.2)  THEN                   ! FOR A'S        
              AR = ASMALLQ (Q, IOR)
          ELSE
              AR = ASMALLQ (Q,-IOR)               ! FOR B'S
          END IF
          AZ = ABFUNC (AR, Q, ICOEF)
CD         WRITE (6,70) AZ
CD  70     FORMAT (' A FOUND BY SMALL A FORMULA ACC. = ',1PE17.10)
          IF (ABS(AZ).GT.PRECR*ABS(AR)+PRECA .AND.
     1                                    Q.GT.QTEST*1.D-3)  THEN
              XDEL = XDEL / 100.D0
              CALL SEARCH (ABFUNC, AR, Q, ICOEF, XDEL, 1)
CD             WRITE (6,80) ABFUNC (AR, Q, ICOEF)
CD  80         FORMAT (' A IMPROVED BY SEARCH, ACC. = ',1PE17.10)
          END IF
      END IF
CD     WRITE (6,90) IR, AR
CD  90 FORMAT (' AB(', I5, ',ICOEF) =', E17.10)
  100 AB(IR,ICOEF) = AR

      QOLD(ICOEF) = Q
      MZ(ICOEF) = MAXZ
CD     WRITE (6,110) ICOEF, QOLD(ICOEF), ICOEF, MZ(ICOEF)
CD 110 FORMAT ('QOLD(',I5,') =', 1PE17.10,',  MZ(',I5,') = ',I5)
      RETURN
      END



C**********************************************************************C
      FUNCTION ABFUNC (A, Q, ICOEF)

C     THIS ROUTINE CALCULATES THE CONTINUED FRACTION FOR GIVEN VALUES
C     OF A AND Q. THE CONTINUED FRACTION USED IS DETERMINED BY ICOEF.
C     NTERM IS THE NUMBER OF TERMS CALCULATED.
C     This routine claculates the continued fraction for the given values of
C     A and Q according to the parameter ICOEF.  It follows Toyama and Shogen
C     in that it finds the zeros by a polymomial expansion of the numerator
C     (AN here).  It differs in that it figures out how many terms are needed,
C     normalizes the polynomial to prevent overflow, and uses an iterative
C     method to calculate the polynomial to any order.

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (PRECR=1.D-16,PRECA=1.D-32)

C   * BEGIN BY CALCULATING NUMBER OF TERMS NEEDED *

      IF (A.LE.0.0D0)  THEN
          AA = 0
      ELSE
          AA = A
      END IF
      NTERM = NINT (SQRT (Q+AA)) + 10

C   * ABFUNC CALCULATIONS *

      IF (ICOEF.EQ.1)  THEN                       ! EVEN FUNCTION, EVEN ORDER
          AN = A
          BN = 2.D0
          DIP = 0.D0

      ELSE IF (ICOEF.EQ.2)  THEN                  ! EVEN FUNCTION, ODD ORDER
          AN = A - 1.D0 - Q
          BN = 1.D0
          DIP = 1.D0

      ELSE IF (ICOEF.EQ.3)  THEN                  ! ODD FUNCTION, EVEN ORDER
          AN = A - 4.D0
          BN = 1.D0
          DIP = 2.D0

      ELSE IF (ICOEF.EQ.4)  THEN                  ! ODD FUNCTION, ODD ORDER
          AN = A - 1.D0 + Q
          BN = 1.D0
          DIP = 1.D0

      ELSE
          WRITE (6,1155) ICOEF
1155      FORMAT (' * * * ILLEGAL VALUE FOR ICOEF (',I5, ') * * *' )
          STOP '   ICOEF ERROR'
      END IF

C   * CALCULATE THE CONTINUED FRACTION FOR MATHIEU FUNCTIONS *
C     CURRENTLY ONLY CALCULATES THE NUMERATOR TO FIND THE ZEROS,
C     REMOVE THE COMMENTS TO CACULATE THE CONTINUED FRACTION.

C     CN = 1.D0
C     DN = 0.D0
      ABSA = ABS (A)
      QSQ = Q**2

      DO 100 I = 1,NTERM
      SQ = (2.D0*DBLE(I) + DIP)**2
      FN = A - SQ
      DIV = - MAX (1.D0, Q, SQ)
      AOLD = AN
      AN = (FN*AN - QSQ*BN) / DIV
      BN = AOLD / DIV
C     COLD = CN
C     CN = (FN*CN - QSQ*DN) / DIV
C     DN = CN / DIV
      CFOLD = ABFUNC
      ABFUNC = AN                                              ! / CN 
C      IF (ABS(ABFUNC-CFOLD).LT.(PRECR*MAX(ABS(ABFUNC),ABSA)+PRECA)
C     1                                             .AND. I.GT.2)  THEN
C          NTERM = I
CD         WRITE (6,1157) A, ABFUNC, NTERM
C          RETURN
C      END IF
  100 CONTINUE
CD     WRITE (6,1157) A, ABFUNC, NTERM
C1157   FORMAT ( ' A, ABFUNC, MAX NTERM:',1PE17.10,', ',E17.10,', ',0PI5)

      RETURN
      END


C**********************************************************************C
      FUNCTION ASMALLQ (Q, IR)

C     This function calculates a zero in the continued frqction equations
C     for small Q.  IR is the zero number, even solutions are indicated by 
C     IR >= 0 and odd solutions by IR < 0.  The equations are from Abramowitz
C     and Stegun, N.B.S. AMS-55, eqs. 20.2.25

      IMPLICIT REAL*8 (A-H,O-Z)

      Q2 = Q**2
      Q4 = Q2**2
      IF (IR.GE.0)  THEN
          QP = -Q
      ELSE
          QP = Q
      END IF

      IF (IR.EQ.0)  THEN                          ! A0(Q)
          ASMALLQ = -Q2 * (1.D0/2.D0 - Q2 * (7.D0/128.D0 - Q2 
     1                    * (29.D0/2304.D0 - Q2*68687.D0/18874368.D0)))

      ELSE IF (ABS(IR).EQ.1)  THEN                ! A1(-Q), B1(Q)
          ASMALLQ = 1.D0 - QP - Q2/8.D0 + Q2*QP/64.D0 - Q4/1536.D0 
     1              - 11.D0*Q4*QP/36864.D0 + 49.D0*Q4*Q2/589824.D0
     2            - 55.D0*Q4*Q2*QP/9437184.D0 - 83.D0*Q4*Q4/35389440.D0

      ELSE IF (IR.EQ.2)  THEN                     ! A2(Q)    
          ASMALLQ = 4.D0 + Q2 * (5.D0/12.D0 - Q2 * (763.D0/13824.D0
     1                - Q2 * (1002401.D0/79626240.D0
     2                - Q2 * 1669068401.D0/458647142400.D0)))

      ELSE IF (IR.EQ.-2)  THEN                    ! B2(Q)
          ASMALLQ = 4.D0 - Q2 * (1.D0/12.D0 - Q2 * (5.D0/13824.D0
     1                - Q2 * (289.D0/79626240.D0
     2                - Q2 * 21391.D0/458647142400.D0)))
    
      ELSE IF (ABS(IR).EQ.3)  THEN                ! A3(-Q), B3(Q)
          ASMALLQ = 9.D0 + Q2/16.D0 - Q2*QP/64.D0 + 13.D0*Q4/20480.D0
     1                + 5.D0*Q4*QP/16384.D0 - 1961.D0*Q4*Q2/23592960.D0
     2                + 609.D0*Q4*Q2*QP/104857600.D0

      ELSE IF (IR.EQ.4)  THEN                     ! A4(Q)
          ASMALLQ = 16.D0 + Q2/30.D0 - 317.D0*Q4/864000.D0
     1                + 5701.D0*Q4*Q2/2.7216D9

      ELSE IF (IR.EQ.-4)  THEN                    ! B4(Q)
          ASMALLQ = 16.D0 + Q2/30.D0 + 433.D0*Q4/864000.D0
     1                + 10049.D0*Q4*Q2/2.7216D9

      ELSE IF (ABS(IR).EQ.5)  THEN                ! A5(-Q), B5(Q)
          ASMALLQ = 25.D0 + Q2/48.D0 + 11.D0*Q4/774144.D0
     1                - Q4*QP/1474456.D0 + 37.D0*Q4*Q2/891813888.D0

      ELSE IF (IR.EQ.6)  THEN                     ! A6(Q)
          ASMALLQ = 36.D0 + Q2/70.D0 + 187.D0*Q4/43904000.D0
     1                + 6743617.D0*Q4*Q2/9.29359872D13

      ELSE IF (IR.EQ.-6)  THEN                    ! B6(Q)
          ASMALLQ = 36.D0 + Q2/70.D0 + 187.D0*Q4/43904000.D0
     1                + 5861633.D0*Q4*Q2/9.29359872D13


      ELSE IF (ABS(IR).GE.7)  THEN                ! ALL OTHER AR(Q), BR(Q)
          R2 = DBLE(IR**2)
          R212 = (R2-1.D0)**2
          D1 = 2.D0 * (R2-1.D0)
          D2 = 16.D0 * (R2-4.D0) * R212 * D1
          D3 = 2.D0 * (R2-9.D0) * R212 * D2
          ASMALLQ = R2 + Q2/D1 + (5.D0*R2+7.D0)*Q4/D2
     1                         + (9.D0*R2**2+58.D0*R2+29.D0)*Q4*Q2/D3

      END IF
      RETURN
      END


C**********************************************************************C
      FUNCTION ALARGEQ (Q, IR)

C     This function calculates a zero in the continued fraction equations
C     for large Q.  IR is the zero number for even solutions and the zero
C     number minus one for odd solutions.  The equations are from Abramowitz
C     and Stegun, N.B.S. AMS-55, eqs. 20.3.30.

C     THIS FUNCTION GIVES AN APPROXIMATION FOR LARGE VALUES OF Q
C     FROM ABRAMOWITZ AND STEGUN, N.B.S. AMS-55 EQ. 20.3.30

      REAL*8 ALARGEQ

      W = 2.D0*IR + 1.D0
      W2 = W**2
      W4 = W2*W2
      W6 = W4*W2
      D0 = (W2 + 3.D0) / (4.D0*W)
      D1 = (5.D0 + 34.D0/W2 + 9.D0/W4) / 4.D0
      D2 = (33.D0 + 410.D0/W2 + 405.D0/W4) / (4.D0*W)
      D3 = (63.D0 + 1260.D0/W2 + 2943.D0/W4 + 486.D0/W6) / W2
      D4 = (527.D0 + 15617.D0/W2 + 69001.D0/W4 + 41607.D0/W6) / (W*W2)
      T5SQRTP = 32.D0 * SQRT(Q) / W2

      ALARGEQ = -2.D0*Q + 2.D0*W*SQRT(Q) - (W2+1.D0)/8.D0
     1         - D0/T5SQRTP - D1/T5SQRTP**2
     2         - D2/T5SQRTP**3 - D3/T5SQRTP**4 - D4/T5SQRTP**5

      RETURN
      END


C**********************************************************************C
      FUNCTION AINIT (Q, ICOEF)

C     This function calculates the first zero in the continued fraction
C     equations.  These are A0, A1, B2, and B1 for ICOEF = 1 to 4 respectively.
C     It calls ASMALLQ or ALARGEQ according to the size of Q.  Its principal
C     contribution is to figure out which one to call.

      IMPLICIT REAL*8 (A-H,O-Z)

      Q2 = Q**2
      IF (ICOEF.EQ.1)  THEN                       ! CALCULATE A0
          IF (Q.LT.1.75D0)  THEN
              AINIT = ASMALLQ (Q, 0)
          ELSE
              AINIT = ALARGEQ (Q, 0)
          END IF

      ELSE IF (ICOEF.EQ.2)  THEN                  ! CALCULATE A1
          IF (Q.LT.4.D0)  THEN
              AINIT = ASMALLQ (Q, 1)
          ELSE
              AINIT = ALARGEQ (Q, 1)
          END IF

      ELSE IF (ICOEF.EQ.3)  THEN                  ! CALCULATE B2
          IF (Q.LT.6.D0)  THEN
              AINIT = ASMALLQ (Q,-2)
          ELSE
              AINIT = ALARGEQ (Q, 1)
          END IF

      ELSE IF (ICOEF.EQ.4)  THEN                  ! CALCULATE B1
          IF (Q.LT.4.D0)  THEN
              AINIT = ASMALLQ (Q ,-1)
          ELSE
              AINIT = ALARGEQ (Q,0)
          END IF
      END IF
      RETURN
      END


C**********************************************************************C
      SUBROUTINE SEARCH (Y, X, Q, IC, XDEL, IFLIP)

C     ZERO CROSSING ROUTINE, SEARCH FOR A ZERO CROSSING FIRST BY STEPPING
C     UNTIL THE ZERO IS BRACKETED AND THEN USE A COMBINATION OF STEPPING 
C     AND REGULA-FALSI METHOD TO FIND THE ZERO.   
C     Modified zero crossing routine.  It begins by seaching for the next
C     zero using a simple step search.  When the zero has been bracketed it
C     uses a regula-falsi followed by an "accelerator" in alteration.  The
C     acclerator consists of dividing the remaining interval by a 10 from
C     the point just found by regula-falsi, it should take care of functions
C     which are ill conditioned for regula-falsi.
C     Parameters in the call:
C         Y is the external function to be evalueated.
C         X is the initial estimate of the zero location, returned as the
C           location of the zero.
C         Q and IC are parameters of the function Y.
C         XDEL is the step size, it may be changed by this routine.
C         IFLIP if equal to one means the direction of the search may be 
C           in the step search.

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (PRECR=1.D-15, PRECA=1.D-35)
      EXTERNAL Y

      X0 = X
      X  = X + XDEL
      Y0 = Y (X0, Q, IC)
      Y1 = Y (X, Q, IC)
      IF (ABS(Y1).GT.ABS(Y0) .AND. IFLIP.EQ.1)  THEN
          XT = X0
          X0 = X
          X  = XT
          YT = Y0
          Y0 = Y1
          Y1 = YT
          XDEL = -XDEL
      END IF

C   * SEARCH FOR ZERO CROSSING LOOP *

      DO 10 I = 1,1000
      IF (Y1*Y0.LE.0.0D0)  GO TO 15
      X0 = X
      Y0 = Y1
      X = X0 + XDEL
      Y1 = Y (X, Q, IC)
   10 CONTINUE

C   * SEARCH FAILED *

      WRITE (6,12) X0, Y0, X, Y1, XDEL
   12   FORMAT ( ' SEARCH ROUTINE FAILED AFTER 1000 ITERATIONS', /
     1           ' X0, Y0 =', 1PE17.10, ', ',E17.10, /
     2           ' X1, Y1 =', E17.10, ', ',E17.10, /
     3           ' XDEL =', E17.10)
      STOP 'SEARCH FAILURE'

C   * CROSSED A ZERO - USE REGULA FALSI AND AN ACCELERATOR 
C                                                 THE REST OF THE WAY *

   15 YZ = Y1
      DO 20 ICOUNT = 1, 100
      IF ((ABS(X-X0).LE.ABS(X0)*PRECR) .OR.
     1                         (ABS(YZ).LE.PRECA))  GO TO 30
      XZ = X0 - (X-X0) * Y0 / (Y1-Y0)             ! REGULA FALSI
      YZ = Y (XZ, Q, IC)
      IF (SIGN(1.D0,YZ).EQ.SIGN(1.D0,Y1))  THEN
          X = XZ
          Y1 = YZ
      ELSE
          X0 = XZ
          Y0 = YZ
      END IF

      XDEL = (X - X0) / 2.D0                      ! DIVIDE INTERVAL IN HALF
      X2 = (X + X0) / 2.D0
      IF (X2.EQ.X0 .OR. X2.EQ.X)  GO TO 30
      Y2 = Y (X2, Q, IC)

      IF (SIGN(1.D0,Y2).EQ.SIGN(1.D0,Y1))  THEN
          X  = X2
          Y1 = Y2
      ELSE
          X0 = X2
          Y0 = Y2
      END IF
   20 CONTINUE
      IF ((ABS(X-X0).LE.ABS(X0)*PRECR) .OR.
     1                         (ABS(YZ).LE.PRECA))  GO TO 30

      WRITE (6,25) X0, X, Y0, Y1, ICOUNT
   25    FORMAT ('    X0, X1:', 1PE17.10, ', ',E17.10, /
     1           '    Y0, Y1:', E17.10, ', ',E17.10, /
     1           '    ICOUNT =', I5 )
      STOP 'SEARCH FAILURE'

   30 CONTINUE
CD     WRITE (6,35) X, X0, Y1, Y0, ICOUNT
CD  35    FORMAT ('       X, X0:', 1PE17.10, ', ',E17.10, /
CD    1           '       Y, Y0:', E17.10, ', ',E17.10, /
CD    2           '       ICOUNT =', I5 )
      IF (ABS(Y0).LT.ABS(Y1))  THEN
          X = X0
          Y1 = Y0
      END IF
      RETURN
      END
