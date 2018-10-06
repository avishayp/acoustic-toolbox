PROGRAM BOUNCE

  ! Program for computing the reflection coefficient, R(theta)
  ! Michael B. Porter

  USE KrakencMod
  USE RefCoMod
  USE sspMod
  IMPLICIT NONE
  INTEGER            :: IAllocStat
  REAL (KIND=8)      :: kMin, kmax, freq0
  COMPLEX   (KIND=8) :: CRCI
  CHARACTER (LEN=80) :: FileRoot

  ! get the file root for naming all input and output files
  ! should add some checks here ...

  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )

  ! *** Read in environmental info ***

  TITLE = 'BOUNCE-'
  iProf = 1

  CALL ReadEnvironment( FileRoot, TITLE, freq, MaxMedium, &
       TopOpt, HSTop, NG,    &
       BotOpt, HSBot, ENVFile, PRTFile )

  READ(  ENVFile, *    ) clow, chigh        ! Spectral limits
  WRITE( PRTFile, "( /, ' cLow = ', G12.5, 'm/s      cHigh = ', G12.5, 'm/s' )" ) cLow, cHigh

  IF ( clow == 0.0 ) CALL ERROUT( PRTFile, 'F', 'BOUNCE', 'clow must be greater than zero' )

  READ(  ENVFile, * ) RMax                  ! Maximum range for calculations
  WRITE( PRTFile, * ) 'RMax = ', RMax
  CLOSE( ENVFile )

  omega2 = ( 2.0 * pi * freq ) ** 2

  CALL ReadReflectionCoefficient( FileRoot, BotOpt( 1 : 1 ), TopOpt( 1 : 1 ), PRTFile ) ! Optionally read in bottom reflection coefficient data

  freq0 = freq   ! save the reference frequency for scaling the grid

  CALL UpdateSSPLoss( freq, freq0 )

  ! update loss in the halfspaces based on the frequency
  ! depth of 1e20 is a large number to ensure bio loss is not included
  IF ( HSTop%BC == 'A' ) THEN
     HSTop%cP  = CRCI( 1D20, HSTop%alphaR, HSTop%alphaI, freq, freq0, SSP%AttenUnit, HSTop%beta, HSTop%fT, bio, NBioLayers )
     HSTop%cS  = CRCI( 1D20, HSTop%betaR,  HSTop%betaI,  freq, freq0, SSP%AttenUnit, HSTop%beta, HSTop%fT, bio, NBioLayers )
  END IF

  IF ( HSBot%BC == 'A' ) THEN
     HSBot%cP  = CRCI( 1D20, HSBot%alphaR, HSBot%alphaI, freq, freq0, SSP%AttenUnit, HSBot%beta, HSBot%fT, bio, NBioLayers )
     HSBot%cS  = CRCI( 1D20, HSBot%betaR,  HSBot%betaI,  freq, freq0, SSP%AttenUnit, HSBot%beta, HSBot%fT, bio, NBioLayers )
  END IF

  ! *** Write internal reflection coefficient data ***

  WRITE( PRTFile, * ) 'Writing internal refl. coeff. table'
  ! Compute NkTab

  kmin = SQRT( omega2 ) / chigh
  kmax = SQRT( omega2 ) / clow
  IF ( chigh > 1.0E6 ) kmin = 0.0

  NkTab = INT( 1000.0 * RMax * ( kmax - kmin ) / ( 2.0 * pi ) )
  WRITE( PRTFile, * ) 'NkTab = ', NkTab

  ALLOCATE( xTab( NkTab ), fTab( NkTab ), gTab( NkTab), ITab( NkTab ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'BOUNCE', 'Too many points in reflection coefficient' )

  CALL ComputeReflCoef( FileRoot, kmin, kmax )

END PROGRAM BOUNCE

!**********************************************************************C

SUBROUTINE Initialize

  ! Initializes arrays defining difference equations

  USE KrakencMod
  USE sspMod
  IMPLICIT NONE
  LOGICAL           :: ElasticFlag = .FALSE.
  INTEGER           :: iAllocStat, ii, j, Med, N1, NPts
  REAL (KIND=8)     :: Twoh
  COMPLEX  (KIND=8) :: cP2, cS2
  COMPLEX (KIND=8), ALLOCATABLE :: cP( : ), cS( : )
  CHARACTER (LEN=8) :: TASK

  cmin          = 1.0E6
  FirstAcoustic = 0
  Loc( 1 )      = 0
  NPTS          = SUM( N( 1 : SSP%NMedia ) ) + SSP%NMedia

  ALLOCATE ( B1( NPTS ), B2( NPTS ), B3( NPTS ), B4( NPTS ), rho( NPTS ), CP( NPTS ), cS( NPTS ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) &
       CALL ERROUT( PRTFile, 'F', 'BOUNCE - INIT', 'Insufficient memory: Reduce mesh.' )

  ! *** Loop over media ***

  MediumLoop: DO Med = 1, SSP%NMedia
     IF ( Med /= 1 ) Loc( Med ) = Loc( Med - 1 ) + N( Med - 1 ) + 1
     N1 = N(   Med ) + 1
     ii = Loc( Med ) + 1

     ! *** EvaluateSSP reads in the data for a medium ***

     TASK = 'TAB'
     CALL EvaluateSSP( cP( ii ), cS( ii ), rho( ii ), Med, N1, freq, TASK, ENVFile, PRTFile )

     ! *** Load diagonals of the finite-difference equations ***

     IF ( cS( ii ) == ( 0.0, 0.0 ) ) THEN ! Case of an acoustic medium ---

        SSP%Material( Med ) = 'ACOUSTIC'
        IF ( FirstAcoustic == 0 ) FirstAcoustic = Med
        LastAcoustic = Med

        cmin = MIN( MINVAL( DBLE( cP( ii : ii + N( Med ) ) ) ), cmin )
        B1( ii : ii + N( Med ) ) = -2.0 +  h( Med ) ** 2 * omega2 / cP( ii : ii + N( Med ) ) ** 2

     ELSE                                ! Case of an elastic medium ---

        IF ( SSP%sigma( Med ) /= 0.0 ) THEN
           WRITE( PRTFile, * ) 'Rough elastic interface not allowed'
           WRITE( PRTFile, * ) 'PROGRAM ABORTING'
           STOP 'ERROR IN BOUNCE: Rough elastic interface not allowed'
        ENDIF

        SSP%Material( Med ) = 'ELASTIC'
        ElasticFlag = .TRUE.
        Twoh   = 2.0 * h( Med )

        DO J = ii, ii + N( Med )
           cmin = MIN( DBLE( cS( J ) ), cmin )

           cP2 = cP( J ) ** 2
           cS2 = cS( J ) ** 2

           B1( J )  = Twoh / ( rho( J ) * cS2 )
           B2( J )  = Twoh / ( rho( J ) * cP2 )
           B3( J )  = 4.0 * Twoh * rho( J ) * cS2 * ( cP2 - cS2 ) / cP2
           B4( J )  = Twoh * ( cP2 - 2.0 * cS2 ) / cP2
           rho( J ) = Twoh * omega2 * rho( J )
        END DO
     ENDIF
  END DO MediumLoop

  ! *** Bottom properties ***

  IF ( HSBot%BC( 1 : 1 ) == 'A' ) THEN
     IF ( HSBot%cS /= ( 0.0, 0.0 ) ) THEN   ! Elastic bottom:
        ElasticFlag = .TRUE.
        cmin = MIN( cmin, DBLE(  HSBot%cS ) )
     ELSE                                   ! Acoustic bottom:
        cmin = MIN( cmin, DBLE(  HSBot%cP ) )
     ENDIF
  ENDIF

  ! *** Top properties ***

  IF ( HSTop%BC( 1 : 1 ) == 'A' ) THEN
     IF (  HSTop%cS /= ( 0.0, 0.0 ) ) THEN   ! Elastic top:
        ElasticFlag = .TRUE.
        cmin = MIN( cmin, DBLE(  HSTop%cS ) )
     ELSE                                    ! Acoustic top:
        cmin = MIN( cmin, DBLE(  HSTop%cP ) )
     ENDIF
  ENDIF

  IF ( ElasticFlag ) cmin = 0.9 * cmin
  clow = MAX( clow, 0.99 * cmin )

END SUBROUTINE Initialize
!**********************************************************************C
SUBROUTINE ComputeReflCoef( FileRoot, kmin, kmax )

  ! Computes the reflection coefficient for k in [kmin, kmax]

  USE KrakencMod
  USE RefCoMod
  USE sspMod
  IMPLICIT NONE
  INTEGER                  :: ik, IPow, itheta, Loops
  REAL      (KIND=8)       :: k0, c0, Deltak, k1, kMin, kMax
  COMPLEX   (KIND=8)       :: x, f, g
  REAL      (KIND=8)       :: kx( NkTab ), kz( NkTab ), theta( NkTab ), R( NkTab ), phase( NkTab )
  COMPLEX   (KIND=8)       :: RCmplx( NkTab ), R1, R2
  CHARACTER (LEN=80)       :: FileRoot

  N( 1 : SSP%NMedia ) = NG( 1 : SSP%NMedia )
  h( 1 : SSP%NMedia ) = ( SSP%Depth( 2 : SSP%NMedia + 1 ) - SSP%Depth( 1 : SSP%NMedia ) ) / N( 1 : SSP%NMedia )

  HV( 1 ) = h( 1 )
  CALL Initialize

  Deltak = ( kmax - kmin ) / ( NkTab - 1 )

  Wavenumber: DO ik = 1, NkTab
     k1 = kmin + ( ik - 1 ) * Deltak
     x  = k1 ** 2

     CALL BCImpedance( x, 'BOT', HSBot, f, g, IPow )  ! Bottom impedance
     CALL AcousticLayers( x, f, g, IPow  )  ! Shoot through acoustic layers
     xTab( ik ) = DBLE( x )
     fTab( ik ) = f
     gTab( ik ) = g
     ITab( ik ) = IPow
  END DO Wavenumber

  IF ( HSTop%BC( 1 : 1 ) == 'A' ) THEN
     c0 = DBLE( HSTop%cP )                    ! use upper halfspace speed for reference if a halfspace was specified
  ELSE
     c0 = 1500
  END IF

  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) 'Reference sound speed = ', c0

  k0    = SQRT( omega2 ) / c0 ! free-space wavenumber
  kx    = SQRT( xTab )        ! horizontal wavenumber

  WHERE( k0 > kx ) ! propagating waves
     kz     = SQRT( k0 ** 2 - kx ** 2 )   ! vertical   wavenumber
     theta  = RadDeg * ATAN2( kz, kx )    ! angle of incidence
     RCmplx =  - ( fTab - i * kz * gTab ) / ( fTab + i * kz * gTab )   ! complex reflection coef.
     R      = ABS( RCmplx )
     phase  = RadDeg * ATAN2( AIMAG( RCmplx ), REAL( RCmplx ) )
  ELSEWHERE        ! evanescent waves
     kz     = 0
     theta  = 0
     RCmplx = 1.0   ! complex reflection coef.
     R      = 1.0
     phase  = 180.0
  END WHERE

  ! unwrap the phase by counting loops in the complex plane
  Loops = 0
  Angle: DO itheta = NkTab - 1, 1, -1
     R1 = RCmplx( itheta )
     R2 = RCmplx( itheta + 1 )
     IF ( AIMAG( R1 ) > 0 .AND. AIMAG( R1 ) < 0 .AND. REAL( R2 ) < 0 ) Loops = Loops + 1
     IF ( AIMAG( R1 ) < 0 .AND. AIMAG( R2 ) > 0 .AND. REAL( R2 ) < 0 ) Loops = Loops - 1
     phase( itheta ) = phase( itheta ) - Loops * 360
  END DO Angle

  OPEN( FILE = TRIM( FileRoot ) // '.irc', UNIT = IRCFile, STATUS = 'UNKNOWN' )   ! Internal Reflection Coef. format
  WRITE( IRCFile, * ) '''', TITLE, '''', freq
  WRITE( IRCFile, * ) NkTab
  WRITE( IRCFile, FMT = "( 5G15.7, I5 )" ) ( xTab( ik ), fTab( ik ), gTab( ik ), ITab( ik ), ik = 1, NkTab )

  OPEN( FILE = TRIM( FileRoot ) // '.brc', UNIT = BRCFile, STATUS = 'UNKNOWN' )   ! Bottom Reflection Coef. format
  !WRITE( BRCFile, * ) TITLE
  !WRITE( BRCFile, * ) freq
  WRITE( BRCFile, * ) NkTab
  DO ik = NkTab, 1, -1
     WRITE( BRCFile, * ) theta( ik ), R( ik ), phase( ik )
  END DO

END SUBROUTINE ComputeReflCoef
!**********************************************************************C
SUBROUTINE AcousticLayers( x, F, G, IPow )

  ! Shoot through acoustic layers

  USE KrakencMod
  USE sspMod
  IMPLICIT NONE
  INTEGER,       PARAMETER :: IPowF = -5
  REAL (KIND=8), PARAMETER :: Roof = 1.0E5, Floor = 1.0E-5
  INTEGER                  :: ii, IPow, Med
  REAL (KIND=8)            :: rhoM
  COMPLEX (KIND=8)         :: x, f, g, P0, P1, P2, h2k2

  IF ( FirstAcoustic == 0 ) RETURN

  ! *** Loop over successive acoustic media ***

  MediumLoop: DO Med = LastAcoustic, FirstAcoustic, -1
     h2k2 = h( Med ) ** 2 * x
     ii   = Loc( Med ) + N( Med ) + 1
     rhoM = rho( ii )
     P1   = -2.0 * g
     P2   = ( B1( ii ) - h2k2 ) * g - 2.0 * h( Med ) * f * rhoM

     ! *** Shoot through a single medium ***

     DO ii = Loc( Med ) + N( Med ), Loc( Med ) + 1, -1

        P0 = P1
        P1 = P2
        P2 = ( h2k2 - B1( ii ) ) * P1 - P0

        DO WHILE ( ABS( DBLE( P2 ) ) > Roof )
           P0   = Floor * P0
           P1   = Floor * P1
           P2   = Floor * P2
           IPow = IPow - IPowF
        END DO
     END DO

     ! g = P'/rho and g = -P since fP + gP'/rho = 0
     f = -( P2 - P0 ) / ( 2.0 * h( Med ) ) / rhoM
     g = -P1
  END DO MediumLoop

END SUBROUTINE AcousticLayers
