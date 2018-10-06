      FUNCTION TWERSK( OPT, OMEGA, BUMDEN, XI, ETA, KX, RHO0, C0 )

C     SURFACE IMPEDANCE MODEL COURTESY STEVE WALES WITH
C     MINOR MODIFICATIONS (MBP 9/23/85)

C     CALCULATES THE SURFACE IMPEDANCE USING TWERSKY SCATTER MODEL

C         OPT IS A MEMBER OF THE SET {'S','H','T','I'}
C            'S' MEANS SOFT PROTUBERANCES
C            'H' MEANS HARD PROTUBERANCES
C            'T' MEANS SOFT PROTUBERANCES, AMPLITUDE ONLY
C            'I' MEANS HARD PROTUBERANCES, AMPLITUDE ONLY

C         ( THE OPTIONS WITH A HARD BACKING HAVE BEEN REMOVED
C           SINCE THEY ARE IMPROPER FOR ICE )

C         BUMDEN = N, THE NUMBER OF PROTUBERANCES PER METER.
C                  ( BUMP DENSITY )
C         XI     = HEIGHT OF PROTUBERANCE FROM THE SURFACE
C         ETA    = HALF-WIDTH OF THE PROTUBERANCE.
C                  (ETA & XI ARE THE LENGTHS OF THE SEMI-MAJOR AND
C                   SEMI-MINOR AXES OF THE ELLIPSE)

      IMPLICIT REAL*8  (A-H, O-Z)
      PARAMETER ( AI=(0.0D0,1.0D0), PI=3.141592653589793D0,
     &          RADDEG=180.D0/PI )

      COMPLEX*16 FS, FM, TWERSK, KX, KZ, COSPHI, SINPHI
      REAL*8 K
      CHARACTER*1 OPT

C     ------ SURFACE IS ALWAYS SOFT:
      ISHS = 0
C     ------ SET PROTUBERANCE HARDNESS FLAG:
      IF ( OPT(1:1) .EQ. 'S' .OR. OPT .EQ. 'T' ) THEN
         ISHP = 0
      ELSE
         ISHP = 1
      ENDIF

      K  = OMEGA / C0
      KZ = SQRT( K**2 - KX**2 )
      SINPHI = KX / K
      COSPHI = KZ / K
      PHI = RADDEG * DATAN2( DBLE(SINPHI), DBLE(COSPHI) )
      FS = FM ( PHI, 180.D0-PHI, K, XI, ETA, ISHS, ISHP )

      IF ( ISHS .EQ. 1 )  THEN                     ! HARD SURFACE
         TWERSK = -RHO0 * C0 * K / (BUMDEN * FS)
      ELSE                                         ! SOFT SURFACE
         TWERSK = -RHO0 * C0 * BUMDEN * FS / ( COSPHI**2 * K)
      END IF

C     ------ FOR AMPLITUDE ONLY THROW-OUT IMAGINARY PART OF IMPEDANCE
      IF ( OPT .EQ. 'T' .OR. OPT .EQ. 'I' ) THEN
         TWERSK = DBLE( TWERSK )
      ENDIF

      RETURN
      END
C**********************************************************************C
      FUNCTION FM (PHI, OEPHI0, K, XI, ETA, ISHS, ISHP)

C     Burke and Twersky's F Function - based on direct Mathieu function
C     evaluation except for high frequencies.
C     This function calculates the F function for low and high frequencies
C     scattering from elliptical cylinders imbedded in a hard or 
C     or soft plane.  It first determines the appropriate routine to call
C     for the ellipse, high or low frequency, and then adds or subtracts
C     the result according to hard or soft surface.

C     Input parameters used are:
C         PHI - Scattered energy angle
C         OEPHI0 - Incident energy angle at 180 - PHI0
C         ISHS - Soft or hard surface [0 or 1]
C         ISHP - Soft or hard protuberance [0 or 1]
C         XI - Protuberance height (semi-axis) from plane
C         ETA - Protuberance half-width (semi-axis) in plane
C         K - Wavenumber

      IMPLICIT REAL*8 (A-F, O-Z)
      REAL*8 K
      COMPLEX*16 FM, GM, GBTM

      PHI0 = 180.D0 - OEPHI0		! RECALCULATE PHI0
      ISAS = 2*ISHS - 1			! -1, SOFT, CALC. ANTISYMETRIC MODES
					! +1, HARD, CALC. SYMETRIC MODES

C   * LOW FREQUENCY APPROXIMATION *

      IF (K*MAX(XI,ETA).LT.50.D0)  THEN
          FM = 2.D0 * GM (PHI, PHI0, K, XI, ETA, ISHP, ISAS)

C   * HIGH FREQUENCY APPROXIMATION *

      ELSE
          FM = GBTM (PHI, PHI0, PHI0, K, XI, ETA, ISHP)
     1             + ISAS * GBTM (PHI, OEPHI0, PHI0, K, XI, ETA, ISHP)
      END IF

      RETURN
      END
C**********************************************************************C

      FUNCTION GM (PHI, PHI0, K, XI, ETA, ISHP, ISAS)

C     Input parameters used are:
C         PHI - Scattered energy angle
C         PHI0 - Incident energy angle at PI - PHI0
C         K - Wavenumber of the waves 
C         XI - Semi-axis of the ellipse along the axis from which the
C              angles are measured.
C         ETA - Semi- axis of the ellipse perpendicular to the axis from
C               which the angles are measured.  The in-plane axis if the
C               ellipses are in a surface.
C         ISHP - Type of protuberance, 0 = soft, 1 = hard
C         ISAS - < 0 Antisymetric modes
C              - = 0 All modes
C                > 0 Symetric modes

      IMPLICIT REAL*8 (A-F,O-Z)
      PARAMETER ( DEGRAD=3.141592653589793D0/180.D0, PREC=1.D-8 )

      REAL*8 K
      COMPLEX*16 GM, RMTHUV, RM, FACTOR

      RMIN = 0.0    ! MBP 5/20/92 (previously uninitialized)

C   * ROTATE SO THAT COORDINATE SYSTEM IS ALONG MAJOR AXIS *
C     PICK OFF MODES WHICH ARE SYM/ANTISYM ABOUT 0 OR 90 DEGS.

      IF (XI.GE.ETA)  THEN		! WANT SYM/ANTISYM MODES ABOUT 90 DEGS
          PHIN = DEGRAD * PHI0
          POUT = DEGRAD * PHI
          A = XI
          B = ETA
          IF (ISAS.LT.0)  THEN		! ANTISYMETRIC MODES 
              IOP = 1
              IOD = 2
          ELSE IF (ISAS.EQ.0)  THEN	! ALL MODES
              IOP = 0
              IOD = 1
          ELSE IF (ISAS.GT.0)  THEN	! SYMETRIC MODES
              IOP = 0
              IOD = 2
          END IF
          IEO1 = 0
          IEO2 = 1

      ELSE				! WANT SYM/ANTISYM MODES ABOUT 0 DEGS
          PHIN = DEGRAD * (90.D0 + PHI0)
          POUT = DEGRAD * (90.D0 + PHI)
          A = ETA
          B = XI
          IF (ISAS.LT.0)  THEN		! ANTISYMETRIC MODES 
              IEO1 = 1
              IEO2 = 1
          ELSE IF (ISAS.EQ.0)  THEN	! ALL MODES
              IEO1 = 0
              IEO2 = 1
          ELSE IF (ISAS.GT.0)  THEN	! SYMETRIC MODES
              IEO1 = 0
              IEO2 = 0
          END IF
          IOP = 0
          IOD = 1
      END IF

C   * COORDINATES ETC *

      NORD = NINT (K*A + 10*(A+B)/A)
      U = K * (A+B) / 2.D0
      V = K * (A-B) / 2.D0
      Q = U * V

C   * FORM THE SUM *

      GM = (0.0D0, 0.0D0)
      ARMAX = 0.0D0

      DO 20 IEO = IEO1,IEO2			! * EVEN/ODD LOOP *

      DO 10 IORD = IEO+IOP,NORD,IOD		! * ORDER LOOP *
      RM  = RMTHUV (U, V, IORD, IEO, ISHP, 5, NORD)
      ABSR = ABS (RM)
      IF (ABSR.LT.RMIN .OR. IORD.EQ.IEO1+IOP)  RMIN = ABSR
      IF (ABSR*PREC.GE.RMIN)  GO TO 15		! NORMALLY PASSED WHEN BESSEL 
						! FUNC. OF 2ND TYPE (B2) BLOWS UP.
      SM1 = AMTHU (PHIN, Q, IORD, IEO, 0, 5, NORD)
      SM2 = AMTHU (POUT, Q, IORD, IEO, 0, 5, NORD)
      FACTOR = SM1 * SM2 * (-DBLE (RM) / RM)
      ABSRX = ABS (-DBLE (RM) / RM)
      IF (ABSRX.GT.ARMAX)  ARMAX = ABSR
C      IF (ABSRX.LE.PREC*ARMAX)  GO TO 15	! LESS SENSITIVE TEST OF RATIO,
						! CAN BE PASSED IF B1 IS SMALL.
   10 GM = GM + FACTOR

      IF (ABSRX/ARMAX.GT.1.D-12)  THEN
          WRITE (6,*) '     * * * ERROR IN GM * * *'
          WRITE (6,*) '       CONVERGENCE FAILURE'
          WRITE (6,*) '     ORDER OF APPROX.=', NORD
          WRITE (6,*) '     ABSRX/ARMAX =',ABSRX/ARMAX
      END IF

   15 CONTINUE
   20 CONTINUE

      RETURN
      END
C**********************************************************************C

      FUNCTION GBTM (PHI, PHI0, PHIZ, K, XI, ETA, ISHP)

      IMPLICIT INTEGER (I,J)
      IMPLICIT REAL*8  (A-E, H, K-Z)
      IMPLICIT COMPLEX*16 (F,G)
      COMPLEX*16 AI
      PARAMETER ( AI=(0.0D0,1.0D0), PI=3.141592653589793D0, PISQ=PI*PI,
     1          DEGRAD=PI/180.D0, RADDEG=180.D0/PI )

      REAL*8   FRAC, GAMMA


C   * * * * *  HIGH FREQUENCY APPROXIMATIONS  * * * * *
C             (ACCORDING TO BURKE & TWERSKY)

C   * PARAMETERS ALWAYS NEEDED *

C     REALLY SHOULD BE USING PHIZ FOR CALCULATION OF L, BUT WE CAN
C     USE PHI0 SINCE IF PHI0 = 180 - PHIZ, COT**2(PHI0) = COT**2(PHIZ)

      RPHI0 = DEGRAD * PHI0					! PHI0 IN RADIANS
      RHO = ETA / XI					! WIDTH TO HEIGHT RATIO
      RHOSQ = RHO**2					! RATIO SQUARED
      SINPHI0 = SIN (RPHI0)				! SIN OF INC. ANGLE
      COSPHI0 = COS (RPHI0)				! COS OF INC. ANGLE
      IF (ABS(SINPHI0).GT.1.D-10)  THEN			! COT OF INC. ANGLE
          COTSQPHI0 = (COSPHI0 / SINPHI0)**2
      ELSE
          COTSQPHI0 = 1.0D20
      END IF
      L = (1.0D0 + RHOSQ**2 * COTSQPHI0) /
     1                          (1.0D0 + RHOSQ * COTSQPHI0)
      L = XI * SQRT (L)					! A MESSY CONSTANT
      GAMMA = TAN (DEGRAD*PHIZ) / RHOSQ
      GAMMA = ATAN (GAMMA)				! ANOTHER ONE

C   * DIFFRACTION (SHADOW-FORMING) TERM CALCULATIONS *

      CPHI0 = (PHI-PHI0) * DEGRAD / 2.D0			! (SCAT-INC)/2
      CPHI1 = (PHI+PHI0) * DEGRAD / 2.D0			! (SCAT+INC)/2
      SINCPHI0 = SIN (CPHI0)				! SIN OF (SCAT-INC)/2
      COSCPHI0 = COS (CPHI0)				! COS OF (SCAT-INC)/2
      COSPHIGAM = COS (CPHI1-GAMMA)			! ANOTHER MESSY CONST.
      IF (ABS(SINCPHI0).GT.1.D-10)  THEN		! * DIFFRACTION TERM *
          COTCPHI0 = COSCPHI0 / SINCPHI0		! COTAN OF (SCAT-INC)/2
          GD = -COTCPHI0 * SIN (2.D0*K*L*COSPHIGAM*SINCPHI0) / 2.D0
      ELSE
          GD = -K * L * COSPHIGAM			! SPEC. CASE FOR SIN 0
      END IF

C   * SET UP SIGN FOR HARD OR SOFT PROTUBERANCES *

      ISIGNP = -1					! SOFT 
      IF (ISHP.EQ.1)  ISIGNP = 1			! HARD 

C   * REFLECTION TERM FOR SEPERATED ANGLES *

      SINCPHI1 = SIN (CPHI1)
      COSCPHI1 = COS (CPHI1)
      P1 = RHOSQ * COSCPHI1**2 + SINCPHI1**2
      EXPL = -2.D0 * K * XI * ABS(SINCPHI0) * SQRT(P1)
      GL = (RHO/2.D0) * SQRT (AI*PI*K*XI) * SQRT (ABS(SINCPHI0))
     1                                    * EXP (AI*EXPL) / P1**(0.75D0)

C   * ANGLES CLOSE TOGETHER - FORWARD DIRECTION *
C     NOTE: CPHI0 = B & T'S BETA

      FRAC  = ABS (K * CPHI0 * XI)
      IF (FRAC.LE.PI .AND. FRAC.GT.1.D-10)  THEN
          SINPG = SIN (RPHI0-GAMMA)			! MORE MESSY CONSTS.
          COSPG = COS (RPHI0-GAMMA)
          GLSA = K*L*CPHI0 * SINPG
     1            + K*(CPHI0**2) * (L*COSPG - AI*K*PI*XI*ETA/2.D0)
          FRAC = (1.D0 + COS(FRAC)) / 2.D0		! TRANS. PT. FRAC*PI/2
          GL = FRAC * GLSA + (1.D0-FRAC) * GL
      END IF

C   * FINAL CALCULATION *

      GBTM = GD + ISIGNP * GL

      RETURN
      END
