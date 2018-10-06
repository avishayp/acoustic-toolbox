PROGRAM FIELD

  ! Generate file of replica vectors
  ! Useage: Field FileRoot
  ! where
  ! FileRoot.mod contains the modes and
  ! FileRoot.shd contains the output shade file

  USE MathConstantsMod
  USE SdRdRMod
  USE beampatternMod
  USE interpMod
  USE SubTabulate
  USE ReadModes
  USE EvaluateCMMod
  USE EvaluateADMod

  IMPLICIT NONE
  INTEGER,   PARAMETER :: FLPFile = 10, PRTFile = 6, SHDFile = 25, MaxM = 20000, MaxNfreq = 1000 ! MaxM also in Evaluate, EvaluateAD, EvaluateCM
  INTEGER              :: iird, IR, iProf, NProf, Nrr, IAllocStat, iRec, IS, M, MLimit, MSrc, IOStat
  REAL                 :: zMin, zMax, atten = 0
  COMPLEX              :: k( MaxM )
  CHARACTER   (LEN=50) :: Opt
  CHARACTER   (LEN=80) :: SHDTitle, Title, FileRoot
  CHARACTER   (LEN=1 ) :: Comp
  CHARACTER   (LEN=10) :: PlotType = '          '
  REAL,    ALLOCATABLE :: rr( : ), rProf( : )
  COMPLEX, ALLOCATABLE :: phiS( :, : ), phiR( :, : ), P( :, : ), C( : ), Ptemp( : )
  REAL (KIND=8), ALLOCATABLE :: kz2( : ), thetaT( : ), S( : )
  REAL        (KIND=8) :: omega, c0
  INTEGER              :: ifreq
  
  ! get the file root for naming all input and output files
  ! should add some checks here ...

  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )

  ! Open the print file
  OPEN( UNIT = PRTFile, FILE = 'field.prt', STATUS = 'UNKNOWN', IOSTAT = iostat )

  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) '__________________________________________________________________________'
  WRITE( PRTFile, * ) 'Running FIELD'
  WRITE( PRTFile, * ) 'Sums modes, producing pressure'
  WRITE( PRTFile, * )

  ! open the field paramaters file (FLPFile)
  OPEN( FILE = TRIM( FileRoot ) // '.flp', UNIT = FLPFile, STATUS = 'OLD', FORM = 'FORMATTED', IOSTAT = IOStat, ACTION = 'READ' )
  IF ( IOStat /= 0 ) THEN
     WRITE( PRTFile, * ) 'FLPFile = ', TRIM( FileRoot ) // '.flp'
     CALL ERROUT( PRTFile, 'F', 'FIELD', 'Unable to open FLPFile' )
  END IF

  SHDTitle( 1 : 1 ) = '$'

  READ( FLPFile, * ) SHDTitle
  READ( FLPFile, * ) Opt
  READ( FLPFile, * ) MLimit

  SELECT CASE ( Opt( 1 : 1 ) )
  CASE ( 'X' )
     WRITE( PRTFile, * ) 'Line source'
  CASE ( 'R' )
     WRITE( PRTFile, * ) 'Point source'
  CASE ( 'S' )
     WRITE( PRTFile, * ) 'Scaled cylindrical (Point source with cylindrical scaling removed)'
  CASE DEFAULT
     CALL ERROUT( PRTFile, 'F', 'FIELD', 'Unknown value for Option( 1 : 1 ); should be X or R' )
  END SELECT

  Comp   = Opt( 3 : 3 )
  MLimit = MIN( MaxM, MLimit )

  ! *** Read profile ranges ***

  READ( FLPFile, * ) NProf

  IF ( NProf == 1      ) THEN
     WRITE( PRTFile, * ) 'Range-independent calculation'
  ELSE
     WRITE( PRTFile, * ) 'Range-dependent calculation'
     SELECT CASE ( Opt( 2 : 2 ) )
     CASE ( 'C' )
        WRITE( PRTFile, * ) 'Coupled modes'
     CASE ( 'A' )
        WRITE( PRTFile, * ) 'Adiabatic modes'
     CASE DEFAULT
        CALL ERROUT( PRTFile, 'F', 'FIELD', 'Unknown option for type of calc; should be Adiabatic or Coupled modes' )
     END SELECT
  END IF

  ! optionally read in a source beam pattern
  SBPFlag = Opt( 3 : 3 )
  CALL ReadPAT( FileRoot, PRTFile )  ! Source Beam Pattern

  ! profile ranges ...
  ALLOCATE( rProf( MAX( 3, NProf + 1 ) ), Stat = IAllocStat )   ! NProf + 1 profiles (one added at `infinite' range)

  IF ( IAllocStat /= 0 ) THEN
     WRITE( PRTFile, * ) 'NProf = ', NProf
     CALL ERROUT( PRTFile, 'F', 'FIELD', 'Dynamic memory allocation failed: Too many profiles' )
  ENDIF

  IF ( NProf > 2 ) rProf( 3 ) = -999.9
  READ( FLPFile, * ) rProf( 1 : NProf )
  CALL SubTab( rProf, NProf )

  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) 'Number of profiles   = ', NProf
  IF ( NProf >= 1  ) WRITE( PRTFile, "( 5G14.6 )" ) ( rProf( iProf ), iProf = 1, MIN( NProf, Number_to_Echo ) )
  IF ( NProf > Number_to_Echo  ) WRITE( PRTFile, * ) ' ... ', rProf( NProf )
  
  ! EvaluateAD/EvaluateCM need a profile at zero range
  IF ( rProf( 1 ) /= 0.0 ) CALL ERROUT( PRTFile, 'F', 'FIELD', 'The first profile must be at 0 km' )

  CALL ReadRcvrRanges( FLPFile, PRTFile )         ! Read receiver ranges
  zMin = -HUGE( zMin )
  zMax = +HUGE( zMax )
  CALL ReadSdRd( FLPFile, PRTFile, zMin, zMax )   ! Read source/receiver depths

  ALLOCATE( phiS( MaxM, Pos%Nsd ), phiR( MaxM, Pos%Nrd ), C( MaxM ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'FIELD', 'Dynamic memory allocation failed: Too many receiver depths' )

  ! *** Read receiver ranges (offsets from vertical) ***

  READ( FLPFile, * ) Nrr

  IF ( Nrr /= Pos%Nrd ) THEN
     WRITE( PRTFile, * ) 'Nrr, Nrd = ', Nrr, Pos%Nrd
     CALL ERROUT( PRTFile, 'W', 'FIELD', 'Nrr being set to Nrd' )
     Nrr = Pos%Nrd
  ENDIF

  ALLOCATE( rr( Nrr ) )
  IF ( Nrr > 1 ) rr( 2 ) = -999.9
  IF ( Nrr > 2 ) rr( 3 ) = -999.9
  READ( FLPFile, * ) rr( 1 : Nrr )

  CALL SubTab( rr, Pos%Nrd )
  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) 'Number of receiver range offsets = ', Pos%Nrd
  WRITE( PRTFile, * ) 'Receiver range offsets (km)'

  IF ( Pos%Nrd >= 1 ) WRITE( PRTFile,  "( 5G14.6 )" ) ( rr( IR ), IR = 1, MIN( Pos%Nrd, Number_to_Echo ) )
  IF ( Pos%Nrd > Number_to_Echo ) WRITE( PRTFile, * ) ' ... ', rr( Pos%Nrd )

  WRITE( PRTFile, * )
  ALLOCATE ( P( Pos%Nrd, Pos%Nr ), Ptemp( Pos%Nr ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'FIELD', ' Insufficient memory to allocate P( Nrd, Nr )' )

  !  *** Read in modes ***

  IProf = 1
  allocate( freqVec( MaxNfreq ) )   !!! must match dimension in GetModes.f90 contained in ReadModes.f90
  !!! the treatment of multiple frequencies is ugly here because we don't know the number of frequencies until we
  !!! read the mode header
  Nfreq = MaxNfreq

  FreqLoop: DO ifreq = 1, MaxNfreq
     IF ( ifreq > Nfreq ) EXIT FreqLoop  !!! Nfreq is updated afer calling GetModes below
     CALL GetModes( FileRoot, iProf, ifreq, MaxM, Pos%sd, Pos%Nsd, 'N' , k, phiS, MSrc, freqVec, Nfreq, Title )
     CALL GetModes( FileRoot, iProf, ifreq, MaxM, Pos%rd, Pos%Nrd, Comp, k, phiR, MSrc, freqVec, Nfreq, Title )

     ! Generate header

     IF ( SHDTitle( 1 : 1 ) == '$' ) SHDTitle = Title

     IF ( ifreq == 1 ) THEN
        CALL WriteHeader( TRIM( FileRoot ) // '.shd', SHDTitle, atten, PlotType )
        iRec = 10
     END IF

     ! *** MAIN LOOP: For each source evaluate and write the field ***

     SourceDepths: DO IS = 1, Pos%Nsd
        M = MIN( MLimit, MSrc )   ! Set number of propagating modes

        C( 1 : MSrc ) = phiS( 1 : MSrc, is )

        ! apply the source beam pattern
        IF ( SBPFlag == '*' .AND. IS == 1 ) THEN
           ALLOCATE( kz2( MSrc ), thetaT( MSrc ), S( MSrc ) )
           c0    = 1500   !!! reference sound speed, should be speed at the source depth
           omega = 2 * pi * freqVec( ifreq )
           kz2   = REAL( omega ** 2 / c0 ** 2 - k( 1 : MSrc ) ** 2 )      ! vertical wavenumber squared
           WHERE ( kz2 < 0 ) kz2 = 0                                      ! remove negative values

           thetaT = RadDeg * ATAN( SQRT( kz2 ) / REAL( k( 1 : MSrc ) ) )  ! calculate the angle in degrees
           CALL interp1( SrcBmPat( :, 1 ), SrcBmPat( :, 2 ), thetaT, S )
           C( 1 : Msrc ) = C( 1 : Msrc ) * REAL( S )         ! apply the shading
        END IF

        IF ( NProf == 1 ) THEN   ! Range-independent case
           CALL Evaluate(                              C, phiR,         Pos%Nrd, Pos%R, Pos%Nr, rr, k, M, Opt, P )
        ELSE
           SELECT CASE ( Opt( 2 : 2 ) )
           CASE ( 'C' )   ! Coupled mode theory
              CALL EvaluateCM( FileRoot, rProf, NProf, C, phiR, Pos%rd, Pos%Nrd, Pos%R, Pos%Nr,     k, M, Opt, P )
           CASE ( 'A' )   ! Adiabatic mode theory
              CALL EvaluateAD( FileRoot, rProf, NProf, ifreq, C, phiR, Pos%rd, Pos%Nrd, Pos%R, Pos%Nr,        M, Opt, P )
           END SELECT
        ENDIF

        ! write out the field
        RcvrDepths: DO iird = 1, Pos%Nrd
           IRec  = IRec + 1
           Ptemp = P( iird, 1 : Pos%Nr )   ! temporary variable, avoids compiler warning
           WRITE( SHDFile, REC = IRec ) Ptemp
        END DO RcvrDepths

     END DO SourceDepths

  END DO FreqLoop

  ! do some clean up (that the compiler should really take care of automatically)
  CLOSE( SHDFile )
  IF ( ALLOCATED( kz2 ) ) DEALLOCATE( kz2, thetaT, S )

  WRITE( PRTfile, * ) 'Field completed successfully'
END PROGRAM
