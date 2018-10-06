PROGRAM BELLHOP

  ! BELLHOP Beam tracing for ocean acoustics

  ! Copyright (C) 2009 Michael B. Porter

  ! This program is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU General Public License as published by
  ! the Free Software Foundation, either version 3 of the License, or
  ! (at your option) any later version.

  ! This program is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ! GNU General Public License for more details.

  ! You should have received a copy of the GNU General Public License
  ! along with this program.  If not, see <http://www.gnu.org/licenses/>.

  ! First version (1983) originally developed with Homer Bucker, Naval Ocean Systems Center
  
  USE MathConstantsMod
  USE bellhopMod
  USE BeamPatternMod
  USE bdryMod
  USE RefCoMod
  USE SdRdRMod
  USE angleMod
  USE sspMod

  IMPLICIT NONE
  LOGICAL, PARAMETER   :: ThreeD = .FALSE., Inline = .FALSE.
  INTEGER              :: jj
  CHARACTER ( LEN=2 )  :: AttenUnit
  CHARACTER ( LEN=80 ) :: FileRoot
  COMPLEX   (KIND= 8 ) :: CRCI   ! not needed if READIN handles the conversion of attenuation units

  ! get the file root for naming all input and output files
  ! should add some checks here ...

  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )
  ! Open the print file
  OPEN( UNIT = PRTFile, FILE = TRIM( FileRoot ) // '.prt', STATUS = 'UNKNOWN', IOSTAT = iostat )

  ! Read in or otherwise initialize inline all the variables used by BELLHOP 

  IF ( Inline ) THEN
     ! NPts, Sigma not used by BELLHOP
     Title = 'BELLHOP- Calibration case with envfil passed as parameters'
     freq = 250
     ! NMedia variable is not used by BELLHOP

     ! *** Boundary information (type of boundary condition and, if a halfspace, then halfspace info)

     AttenUnit         = 'W'
     Bdry%Top%HS%BC    = 'V'
     Bdry%Top%HS%Depth = 0
     Bdry%Bot%HS%Depth = 100
     Bdry%Bot%HS%Opt   = 'A_'
     Bdry%Bot%HS%BC    = 'A'
     Bdry%Bot%HS%cp    = CRCI( 1e20, 1590D0,    0.5D0, freq, AttenUnit, bio )   ! compressional wave speed
     Bdry%Bot%HS%cs    = CRCI( 1e20,    0D0,      0D0, freq, AttenUnit, bio )   ! shear         wave speed
     Bdry%Bot%HS%rho   = 1.2                               ! density

     ! *** sound speed in the water column ***

     SSP%Type = 'C'   ! interpolation method for SSP
     SSP%NPts = 2     ! number of SSP points
     SSP%z(  1 : 2 ) = [    0,  100 ]
     SSP%c(  1 : 2 ) = [ 1500, 1500 ]
     SSP%cz( 1 : 2 ) = [    0,    0 ]   ! user should really not have to supply this ...

     ! *** source and receiver positions ***

     Pos%Nsd = 1
     Pos%Nrd = 100
     Pos%Nr  = 500

     ALLOCATE( Pos%sd( Pos%Nsd ), Pos%ws( Pos%Nsd ), Pos%isd( Pos%Nsd ) )
     ALLOCATE( Pos%rd( Pos%Nrd ), Pos%wr( Pos%Nrd ), Pos%ird( Pos%Nrd ) )
     ALLOCATE( Pos%r(  Pos%Nr  ) )

     Pos%sd( 1 ) = 50.
     !Pos%rd     = [ 0, 50, 100 ]
     !Pos%r      = 1000. * [ 1, 2, 3, 4, 5 ]   ! meters !!!
     Pos%rd      = [ ( jj, jj = 1, Pos%Nrd ) ]
     Pos%r       = 10. * [ ( jj, jj = 1 , Pos%Nr ) ]   ! meters !!!

     Beam%RunType = 'C'
     Beam%Type    = 'G   '
     Beam%deltas  = 0
     Beam%Box%z   = 101.
     Beam%Box%r   = 5100 ! meters

     Angles%Nalpha = 1789
     ! Angles%alpha  = [ -80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80 ] ! -89 89
     Angles%alpha  = ( 180. / Angles%Nalpha ) * [ ( jj, jj = 1, Angles%Nalpha ) ] - 90.

     ! *** altimetry ***

     ALLOCATE( Top( 2 ), Stat = IAllocStat )
     IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'BELLHOP', 'Insufficient memory for altimetry data'  )
     Top( 1 )%x = [ -sqrt( huge( Top( 1 )%x( 1 ) ) ) / 1.0d5, 0.d0 ]
     Top( 2 )%x = [  sqrt( huge( Top( 1 )%x( 1 ) ) ) / 1.0d5, 0.d0 ]

     CALL ComputeBdryTangentNormal( Top, 'Top' )

     ! *** bathymetry ***

     ALLOCATE( Bot( 2 ), Stat = IAllocStat )
     IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'BELLHOP', 'Insufficient memory for bathymetry data'  )
     Bot( 1 )%x = [ -sqrt( huge( Bot( 1 )%x( 1 ) ) ) / 1.0d5, 5000.d0 ]
     Bot( 2 )%x = [  sqrt( huge( Bot( 1 )%x( 1 ) ) ) / 1.0d5, 5000.d0 ]

     CALL ComputeBdryTangentNormal( Bot, 'Bot' )

     ALLOCATE( RBot( 1 ), Stat = IAllocStat )   ! bottom reflection coefficient
     ALLOCATE( RTop( 1 ), Stat = iAllocStat )   ! top    reflection coefficient

     ! *** Source Beam Pattern ***
     NSBPPts = 2
     ALLOCATE( SrcBmPat( 2, 2 ), Stat = IAllocStat )
     IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'BELLHOP-ReadPat', 'Insufficient memory'  )
     SrcBmPat( 1, : ) = [ -180.0, 0.0 ]
     SrcBmPat( 2, : ) = [  180.0, 0.0 ]
     SrcBmPat( :, 2 ) = 10 ** ( SrcBmPat( :, 2 ) / 20 )  ! convert dB to linear scale !!!
  ELSE
     CALL ReadEnvironment(  FileRoot, ThreeD )
     CALL ReadATI( FileRoot, Bdry%Top%HS%Opt( 5 : 5 ), Bdry%Top%HS%Depth, PRTFile )                          ! AlTImetry
     CALL ReadBTY( FileRoot, Bdry%Bot%HS%Opt( 2 : 2 ), Bdry%Bot%HS%Depth, PRTFile )                          ! BaThYmetry
     CALL ReadReflectionCoefficient( FileRoot, Bdry%Bot%HS%Opt( 1 : 1 ), Bdry%Top%HS%Opt( 2 : 2 ), PRTFile ) ! (top and bottom)
     SBPFlag = Beam%RunType( 3 : 3 )
     CALL ReadPat( FileRoot, PRTFile )                                                                       ! Source Beam Pattern
  END IF

  CALL OpenOutputFiles( FileRoot, ThreeD )
  CALL BellhopCore

END PROGRAM BELLHOP

! **********************************************************************!

SUBROUTINE BellhopCore

  USE bellhopMod
  USE RefCoMod
  USE bdryMod
  USE angleMod
  USE SdRdRMod
  USE ArrMod
  USE BeamPatternMod
  USE sspMod

  IMPLICIT NONE
  LOGICAL, PARAMETER   :: ThreeD = .FALSE.
  INTEGER, PARAMETER   :: SHDFile = 25, RAYFile = 21, ArrivalsStorage = 20000000
  INTEGER              :: IBPvec( 1 ), ibp, is, iBeamWindow2, Ird1, Irec, NalphaOpt, iSeg
  REAL                 :: Tstart, Tstop
  REAL        (KIND=8) :: Amp0, alpha0, DalphaOpt, xs( 2 ), RadMax, s, &
                          c, cimag, gradc( 2 ), crr, crz, czz, rho
  COMPLEX, ALLOCATABLE :: U( :, : )
  COMPLEX     (KIND=8) :: epsilon, PickEpsilon
  COMPLEX     (KIND=8) :: CRCI   ! function to convert soundspeed from user units to program units

  CALL CPU_TIME( Tstart )

  omega = 2.0 * pi * freq

  IF ( Beam%deltas == 0.0 ) THEN
     Beam%deltas = ( Bdry%Bot%HS%Depth - Bdry%Top%HS%Depth ) / 10.0   ! Automatic step size selection
     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) 'Step length,      deltas = ', Beam%deltas, 'm (automatically selected)'
  END IF

  Angles%alpha  = DegRad * Angles%alpha   ! convert to radians
  Angles%Dalpha = 0.0
  IF ( Angles%Nalpha /= 1 ) &
       Angles%Dalpha = ( Angles%alpha( Angles%Nalpha ) - Angles%alpha( 1 ) ) / ( Angles%Nalpha - 1 )  ! angular spacing between beams

  ! convert range-dependent geoacoustic parameters from user to program units
  IF ( atiType( 2 : 2 ) == 'L' ) THEN
     DO iSeg = 1, NatiPts
        Top( iSeg )%HS%cp = CRCI( 1e20, Top( iSeg )%HS%alphaR, Top( iSeg )%HS%alphaI, freq, freq, 'W ', &
             betaPowerLaw, ft, bio, NBioLayers )   ! compressional wave speed
        Top( iSeg )%HS%cs = CRCI( 1e20, Top( iSeg )%HS%betaR,  Top( iSeg )%HS%betaI,  freq, freq, 'W ', &
             betaPowerLaw, ft, bio, NBioLayers )   ! shear         wave speed
     END DO
  END IF
   
  IF ( btyType( 2 : 2 ) == 'L' ) THEN
     DO iSeg = 1, NbtyPts
        Bot( iSeg )%HS%cp = CRCI( 1e20, Bot( iSeg )%HS%alphaR, Bot( iSeg )%HS%alphaI, freq, 'W ', &
             betaPowerLaw, ft, bio, NBioLayers )   ! compressional wave speed
        Bot( iSeg )%HS%cs = CRCI( 1e20, Bot( iSeg )%HS%betaR,  Bot( iSeg )%HS%betaI,  freq, 'W ', &
             betaPowerLaw, ft, bio, NBioLayers )   ! shear         wave speed
     END DO
  END IF

  SELECT CASE ( Beam%RunType( 5 : 5 ) )
  CASE ( 'I' )
     Nrd_per_range = 1         ! irregular grid
  CASE DEFAULT
     Nrd_per_range = Pos%Nrd   ! rectilinear grid
  END SELECT

  ! for a TL calculation, allocate space for the pressure matrix
  SELECT CASE ( Beam%RunType( 1 : 1 ) )
  CASE ( 'C', 'S', 'I' )        ! TL calculation
     ALLOCATE ( U( Nrd_per_range, Pos%Nr ), Stat = IAllocStat )
     IF ( IAllocStat /= 0 ) &
          CALL ERROUT( PRTFile, 'F', 'BELLHOP', 'Insufficient memory for TL matrix: reduce Nr * Nrd'  )
  CASE ( 'A', 'a', 'R', 'E' )   ! Arrivals calculation
     ALLOCATE ( U( 1, 1 ), Stat = IAllocStat )   ! open a dummy variable
  END SELECT

  ! for an arrivals run, allocate space for arrivals matrices
  SELECT CASE ( Beam%RunType( 1 : 1 ) )
  CASE ( 'A', 'a' )
     MaxNArr = MAX( ArrivalsStorage / ( Nrd_per_range * Pos%Nr ), 10 )   ! allow space for at least 10 arrivals
     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) '( Maximum # of arrivals = ', MaxNArr, ')'

     ALLOCATE ( Arr( Nrd_per_range, Pos%Nr, MaxNArr ), NArr( Nrd_per_range, Pos%Nr ), Stat = IAllocStat )
     IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'BELLHOP', &
          'Insufficient memory to allocate arrivals matrix; reduce parameter ArrivalsStorage' )
  CASE DEFAULT
     MaxNArr = 1
     ALLOCATE ( Arr( Nrd_per_range, Pos%Nr, 1 ), NArr( Nrd_per_range, Pos%Nr ), Stat = IAllocStat )
  END SELECT

  NArr( 1:Nrd_per_range, 1:Pos%Nr ) = 0

  WRITE( PRTFile, * )

  SourceDepth: DO is = 1, Pos%Nsd
     xs = [ 0.0, Pos%sd( is ) ]   ! source coordinate

     SELECT CASE ( Beam%RunType( 1 : 1 ) )
     CASE ( 'C', 'S', 'I' ) ! TL calculation, zero out pressure matrix
        U = 0.0
     CASE ( 'A', 'a' )      ! Arrivals calculation, zero out arrival matrix
        NArr = 0
     END SELECT

     CALL EvaluateSSP( xs, c, cimag, gradc, crr, crz, czz, rho, Freq, 'TAB' )
     RadMax = 5 * c / freq  ! 5 wavelength max radius

     ! Are there enough beams?
     DalphaOpt = SQRT( c / ( 6.0 * freq * Pos%r( Pos%Nr ) ) )
     NalphaOpt = 2 + INT( ( Angles%alpha( Angles%Nalpha ) - Angles%alpha( 1 ) ) / DalphaOpt )

     IF ( Beam%RunType( 1 : 1 ) == 'C' .AND. Angles%Nalpha < NalphaOpt ) THEN
        CALL ERROUT( PRTFile, 'W', 'BELLHOP', 'Too few beams' )
        WRITE( PRTFile, * ) 'Nalpha should be at least = ', NalphaOpt
     ENDIF

     ! Trace successive beams

     ElevationAngle: DO ialpha = 1, Angles%Nalpha

        IF ( Angles%iSingle == 0 .OR. ialpha == Angles%iSingle ) THEN    ! Single beam run?
        !IF ( ialpha == 900 .or. ialpha == 2222 ) THEN    ! Single beam run?

           alpha0 = Angles%alpha( ialpha ) * RadDeg   ! take-off angle in degrees

           IBPvec = maxloc( SrcBmPat( :, 1 ), mask = SrcBmPat( :, 1 ) < alpha0 )  ! index of ray angle in beam pattern
           IBP    = IBPvec( 1 )
           IBP    = MAX( IBP, 1 )               ! don't go before beginning of table
           IBP    = MIN( IBP, NSBPPts - 1 )     ! don't go past end of table

           ! linear interpolation to get amplitude
           s    = ( alpha0  - SrcBmPat( IBP, 1 ) ) / ( SrcBmPat( IBP + 1, 1 ) - SrcBmPat( IBP, 1 ) )
           Amp0 = ( 1 - s ) * SrcBmPat( IBP, 2 ) + s * SrcBmPat( IBP + 1, 2 )

           ! Lloyd mirror pattern for semi-coherent option
           IF ( Beam%RunType( 1 : 1 ) == 'S' ) &
              Amp0 = Amp0 * SQRT( 2.0 ) * ABS( SIN( omega / c * xs( 2 ) * SIN( Angles%alpha( ialpha ) ) ) )

           ! show progress ...
           IF ( MOD( ialpha - 1, max( Angles%Nalpha / 50, 1 ) ) == 0 ) THEN
              WRITE( PRTFile, FMT = "( 'Tracing beam ', I7, F10.2 )" ) ialpha, alpha0
           END IF

           CALL TraceRay2D( xs, Angles%alpha( ialpha ), Amp0 )   ! Trace a ray

           IF ( Beam%RunType( 1 : 1 ) == 'R' ) THEN   ! Write the ray trajectory to RAYFile
              CALL WriteRay( alpha0, Beam%Nsteps )
           ELSE                                       ! Compute the contribution to the field

              epsilon = PickEpsilon( Beam%Type( 1 : 2 ), omega, c, gradc, Angles%alpha( ialpha ), &
                   Angles%Dalpha, Beam%rLoop, Beam%epsMultiplier ) ! 'optimal' beam constant

              SELECT CASE ( Beam%Type( 1 : 1 ) )
              CASE ( 'R' )
                 iBeamWindow2 = Beam%iBeamWindow **2
                 RadMax       = 50 * c / freq  ! 50 wavelength max radius
                 CALL InfluenceCervenyRayCen(   U, epsilon, Angles%alpha( ialpha ), iBeamWindow2, RadMax )
              CASE ( 'C' )
                 iBeamWindow2 = Beam%iBeamWindow **2
                 RadMax       = 50 * c / freq  ! 50 wavelength max radius
                 CALL InfluenceCervenyCart(     U, epsilon, Angles%alpha( ialpha ), iBeamWindow2, RadMax )
              CASE ( 'g' )
                 CALL InfluenceGeoHatRayCen(    U, Angles%alpha( ialpha ), Angles%Dalpha )
              CASE ( 'S' )
                 CALL InfluenceSGB(             U, Angles%alpha( ialpha ), Angles%Dalpha )
              CASE ( 'B' )
                 CALL InfluenceGeoGaussianCart( U, Angles%alpha( ialpha ), Angles%Dalpha )
              CASE DEFAULT
                 CALL InfluenceGeoHatCart(      U, Angles%alpha( ialpha ), Angles%Dalpha )
              END SELECT

           END IF
        END IF
     END DO ElevationAngle

     ! write results to disk

     SELECT CASE ( Beam%RunType( 1 : 1 ) )
     CASE ( 'C', 'S', 'I' )   ! TL calculation
        CALL ScalePressure( Angles%Dalpha, ray2D( 1 )%c, Pos%r, U, Nrd_per_range, Pos%Nr, Beam%RunType, freq )
        IRec = 10 + Nrd_per_range * ( is - 1 )
        RcvrDepth: DO Ird1 = 1, Nrd_per_range
           IRec = IRec + 1
           WRITE( SHDFile, REC = IRec ) U( Ird1, 1 : Pos%Nr )
        END DO RcvrDepth
     CASE ( 'A' )             ! arrivals calculation, ascii
        CALL WriteArrivalsASCII(  Pos%r, Nrd_per_range, Pos%Nr, Beam%RunType( 4 : 4 ) )
     CASE ( 'a' )             ! arrivals calculation, binary
        CALL WriteArrivalsBinary( Pos%r, Nrd_per_range, Pos%Nr, Beam%RunType( 4 : 4 ) )
     END SELECT

  END DO SourceDepth

  ! close all files
  SELECT CASE ( Beam%RunType( 1 : 1 ) )
  CASE ( 'C', 'S', 'I' )      ! TL calculation
     CLOSE( SHDFile )
  CASE ( 'A', 'a' )           ! arrivals calculation
     CLOSE( ARRFile )
  CASE ( 'R' )                ! ray trace
     CLOSE( RAYFile )
  END SELECT

  ! Display run time
  CALL CPU_TIME( Tstop )
  WRITE( PRTFile, "( /, ' CPU Time = ', G15.3, 's' )" ) Tstop - Tstart
END SUBROUTINE BellhopCore

! **********************************************************************!

FUNCTION PickEpsilon( BeamType, omega, c, gradc, alpha, Dalpha, rLoop, EpsMultiplier )

  ! Picks the optimum value for epsilon

  IMPLICIT NONE
  INTEGER,            PARAMETER     :: PRTFile = 6
  COMPLEX,            PARAMETER     :: i = ( 0.0, 1.0 )
  REAL      (KIND=8), INTENT( IN  ) :: omega, c, gradc( 2 ) ! angular frequency, sound speed and gradient
  REAL      (KIND=8), INTENT( IN  ) :: alpha, Dalpha        ! angular spacing for ray fan
  REAL      (KIND=8), INTENT( IN  ) :: epsMultiplier, Rloop ! multiplier, loop range
  CHARACTER (LEN= 2), INTENT( IN  ) :: BeamType
  LOGICAL, SAVE      :: INIFlag = .TRUE.
  REAL      (KIND=8) :: HalfWidth
  REAL      (KIND=8) :: cz
  COMPLEX   (KIND=8) :: epsilonOpt
  COMPLEX   (KIND=8) :: PickEpsilon
  CHARACTER (LEN=40) :: TAG

  SELECT CASE ( BeamType( 1 : 1 ) )
  CASE ( 'C', 'R' )   ! Cerveny beams
     TAG    = 'Cerveny style beam'
     SELECT CASE ( BeamType( 2 : 2 ) )
     CASE ( 'F' )
        TAG       = 'Space filling beams'
        halfwidth = 2.0 / ( ( omega / c ) * Dalpha )
        epsilonOpt    = i * 0.5 * omega * halfwidth ** 2
     CASE ( 'M' )
        TAG       = 'Minimum width beams'
        halfwidth = SQRT( 2.0 * c * 1000.0 * rLoop / omega )
        epsilonOpt    = i * 0.5 * omega * halfwidth ** 2
     CASE ( 'W' )
        TAG       = 'WKB beams'
        halfwidth = HUGE( halfwidth )
        cz        = gradc( 2 )
        IF ( cz == 0.0 ) THEN
           epsilonOpt = 1.0D10
        ELSE
           epsilonOpt = ( -SIN( alpha ) / COS( alpha ** 2 ) ) * c * c / cz
        ENDIF
     END SELECT

  CASE ( 'G', 'g' )   ! geometric hat beams
     TAG        = 'Geometic hat beams'
     halfwidth  = 2.0 / ( ( omega / c ) * Dalpha )
     epsilonOpt = i * 0.5 * omega * halfwidth ** 2

  CASE ( 'B' )   ! geometric Gaussian beams
     TAG        = 'Geometric Gaussian beams'
     halfwidth  = 2.0 / ( ( omega / c ) * Dalpha )
     epsilonOpt = i * 0.5 * omega * halfwidth ** 2

  CASE ( 'S' )   ! simple Gaussian beams
     TAG        = 'Simple Gaussian beams'
     halfwidth  = 2.0 / ( ( omega / c ) * Dalpha )
     epsilonOpt = i * 0.5 * omega * halfwidth ** 2
  END SELECT

  PickEpsilon = EpsMultiplier * epsilonOpt

  ! On first call write info to prt file
  IF ( INIFlag ) THEN
     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) TAG
     WRITE( PRTFile, * ) 'halfwidth  = ', halfwidth
     WRITE( PRTFile, * ) 'epsilonOpt = ', epsilonOpt
     WRITE( PRTFile, * ) 'EpsMult    = ', EpsMultiplier
     WRITE( PRTFile, * )
     INIFlag = .FALSE.
  END IF

END FUNCTION PickEpsilon

! **********************************************************************!

SUBROUTINE TraceRay2D( xs, alpha, Amp0 )

  ! Traces the beam corresponding to a particular take-off angle

  USE bellhopMod
  USE bdryMod
  USE RefCoMod
  USE sspMod
  IMPLICIT NONE
  REAL     (KIND=8), INTENT( IN ) :: alpha, Amp0  ! initial angle, amplitude
  REAL     (KIND=8), INTENT( IN ) :: xs( 2 )      ! x-y coordinate of the source
  INTEGER           :: is, is1                    ! index for a step along the ray
  REAL     (KIND=8) :: c, cimag, gradc( 2 ), crr, crz, czz, rho
  REAL     (KIND=8) :: dEndTop( 2 ), dEndBot( 2 ), TopnInt( 2 ), BotnInt( 2 ), ToptInt( 2 ), BottInt( 2 )
  REAL     (KIND=8) :: DistBegTop, DistEndTop, DistBegBot, DistEndBot ! Distances from ray beginning, end to top and bottom
  REAL     (KIND=8) :: sss

  ! Initial conditions

  iSmallStepCtr = 0
  CALL EvaluateSSP( xs, c, cimag, gradc, crr, crz, czz, rho, Freq, 'TAB' )
  ray2D( 1 )%c         = c
  ray2D( 1 )%x         = xs
  ray2D( 1 )%t         = [ COS( alpha ), SIN( alpha ) ] / c
  ray2D( 1 )%p         = [ 1.0, 0.0 ]
  ray2D( 1 )%q         = [ 0.0, 1.0 ]
  ray2D( 1 )%tau       = 0.0
  ray2D( 1 )%Amp       = Amp0
  ray2D( 1 )%Phase     = 0.0
  ray2D( 1 )%NumTopBnc = 0
  ray2D( 1 )%NumBotBnc = 0

  ! second component of qv is not used in geometric beam tracing
  ! set I.C. to 0 in hopes of saving run time
  IF ( Beam%RunType( 2 : 2 ) == 'G' ) ray2D( 1 )%q = [ 0.0, 0.0 ]

  CALL GetTopSeg( xs( 1 ) )   ! identify the top    segment above the source
  CALL GetBotSeg( xs( 1 ) )   ! identify the bottom segment below the source

  ! convert range-dependent geoacoustic parameters from user to program units
  ! compiler is not accepting the copy of the whole structure at once ...
  IF ( atiType( 2 : 2 ) == 'L' ) THEN
     Bdry%Top%HS%cp  = Top( IsegTop )%HS%cp   ! grab the geoacoustic info for the new segment
     Bdry%Top%HS%cs  = Top( IsegTop )%HS%cs
     Bdry%Top%HS%rho = Top( IsegTop )%HS%rho
  END IF

  IF ( btyType( 2 : 2 ) == 'L' ) THEN
     Bdry%Bot%HS%cp  = Bot( IsegBot )%HS%cp
     Bdry%Bot%HS%cs  = Bot( IsegBot )%HS%cs
     Bdry%Bot%HS%rho = Bot( IsegBot )%HS%rho
  END IF

  ! Trace the beam (note that Reflect alters the step index is)
  is = 0
  CALL Distances2D( ray2D( 1 )%x, Top( IsegTop )%x, Bot( IsegBot )%x, dEndTop,    dEndBot,  &
       Top( IsegTop )%n, Bot( IsegBot )%n, DistBegTop, DistBegBot )

  IF ( DistBegTop <= 0 .OR. DistBegBot <= 0 ) THEN
     Beam%Nsteps = 1
     WRITE( PRTFile, * ) 'Terminating the ray trace because the source is on or outside the boundaries'
     RETURN       ! source must be within the medium
  END IF

  Stepping: DO istep = 1, MaxN - 1
     is  = is + 1
     is1 = is + 1

     CALL Step2D( ray2D( is ), ray2D( is1 ),  &
          Top( IsegTop )%x, Top( IsegTop )%n, &
          Bot( IsegBot )%x, Bot( IsegBot )%n )

     ! New altimetry segment?
     IF ( ray2D( is1 )%x( 1 ) < rTopSeg( 1 ) .OR. &
          ray2D( is1 )%x( 1 ) > rTopSeg( 2 ) ) THEN
        CALL GetTopSeg( ray2D( is1 )%x( 1 ) )
        IF ( atiType( 2 : 2 ) == 'L' ) THEN
           Bdry%Top%HS%cp  = Top( IsegTop )%HS%cp   ! grab the geoacoustic info for the new segment
           Bdry%Top%HS%cs  = Top( IsegTop )%HS%cs
           Bdry%Top%HS%rho = Top( IsegTop )%HS%rho
        END IF
     END IF

     ! New bathymetry segment?
     IF ( ray2D( is1 )%x( 1 ) < rBotSeg( 1 ) .OR. &
          ray2D( is1 )%x( 1 ) > rBotSeg( 2 ) ) THEN
        CALL GetBotSeg( ray2D( is1 )%x( 1 ) )
        IF ( btyType( 2 : 2 ) == 'L' ) THEN
           Bdry%Bot%HS%cp  = Bot( IsegBot )%HS%cp   ! grab the geoacoustic info for the new segment
           Bdry%Bot%HS%cs  = Bot( IsegBot )%HS%cs
           Bdry%Bot%HS%rho = Bot( IsegBot )%HS%rho
        END IF
     END IF

     ! Reflections?
     ! Tests that ray at step is is inside, and ray at step is+1 is outside
     ! to detect only a crossing from inside to outside
     ! DistBeg is the distance at step is,   which is saved
     ! DistEnd is the distance at step is+1, which needs to be calculated

     CALL Distances2D( ray2D( is1 )%x, Top( IsegTop )%x, Bot( IsegBot )%x, dEndTop,    dEndBot,  &
          Top( IsegTop )%n, Bot( IsegBot )%n, DistEndTop, DistEndBot )

     IF      ( DistBegTop > 0.0d0 .AND. DistEndTop <= 0.0d0 ) THEN  ! test top reflection

        IF ( atiType == 'C' ) THEN
           sss     = DOT_PRODUCT( dEndTop, Top( IsegTop )%t ) / Top( IsegTop )%Len   ! proportional distance along segment
           TopnInt = ( 1 - sss ) * Top( IsegTop )%Noden + sss * Top( 1 + IsegTop )%Noden
           ToptInt = ( 1 - sss ) * Top( IsegTop )%Nodet + sss * Top( 1 + IsegTop )%Nodet
        ELSE
           TopnInt = Top( IsegTop )%n   ! normal is constant in a segment
           ToptInt = Top( IsegTop )%t
        END IF

        CALL Reflect2D( is, Bdry%Top%HS, 'TOP', ToptInt, TopnInt, Top( IsegTop )%kappa, RTop, NTopPTS )
        ray2D( is + 1 )%NumTopBnc = ray2D( is )%NumTopBnc + 1

        CALL Distances2D( ray2D( is + 1 )%x, Top( IsegTop )%x, Bot( IsegBot )%x, dEndTop,    dEndBot,  &
             Top( IsegTop )%n, Bot( IsegBot )%n, DistEndTop, DistEndBot )

     ELSE IF ( DistBegBot > 0.0d0 .AND. DistEndBot <= 0.0d0 ) THEN  ! test bottom reflection

        IF ( btyType == 'C' ) THEN
           sss     = DOT_PRODUCT( dEndBot, Bot( IsegBot )%t ) / Bot( IsegBot )%Len   ! proportional distance along segment
           BotnInt = ( 1 - sss ) * Bot( IsegBot )%Noden + sss * Bot( 1 + IsegBot )%Noden
           BottInt = ( 1 - sss ) * Bot( IsegBot )%Nodet + sss * Bot( 1 + IsegBot )%Nodet
        ELSE
           BotnInt = Bot( IsegBot )%n   ! normal is constant in a segment
           BottInt = Bot( IsegBot )%t
        END IF

        CALL Reflect2D( is, Bdry%Bot%HS, 'BOT', BottInt, BotnInt, Bot( IsegBot )%kappa, RBot, NBotPTS )
        ray2D( is + 1 )%NumBotBnc = ray2D( is )%NumBotBnc + 1
        CALL Distances2D( ray2D( is + 1 )%x, Top( IsegTop )%x, Bot( IsegBot )%x, dEndTop,    dEndBot, &
             Top( IsegTop )%n, Bot( IsegBot )%n, DistEndTop, DistEndBot )

     END IF

     ! Has the ray left the box, lost its energy, escaped the boundaries, or exceeded storage limit?
     IF ( ABS( ray2D( is + 1 )%x( 1 ) ) > Beam%Box%r .OR. &
          ABS( ray2D( is + 1 )%x( 2 ) ) > Beam%Box%z .OR. ray2D( is + 1 )%Amp < 0.005 .OR. &
          ( DistBegTop < 0.0 .AND. DistEndTop < 0.0 ) .OR. &
          ( DistBegBot < 0.0 .AND. DistEndBot < 0.0 ) ) THEN
          ! ray2D( is + 1 )%t( 1 ) < 0 ) THEN ! this last test kills off a backward traveling ray
        Beam%Nsteps = is + 1
        EXIT Stepping
     ELSE IF ( is >= MaxN - 3 ) THEN
        CALL ERROUT( PRTFile, 'W', 'TraceRay2D', 'Insufficient storage for ray trajectory' )
        Beam%Nsteps = is
        EXIT Stepping
     END IF

     DistBegTop = DistEndTop
     DistBegBot = DistEndBot

  END DO Stepping
END SUBROUTINE TraceRay2D

! **********************************************************************!

SUBROUTINE Distances2D( rayx, Topx, Botx, dTop, dBot, Topn, Botn, DistTop, DistBot )

  ! Calculates the distances to the boundaries
  ! Formula differs from JKPS because code uses outward pointing normals

  IMPLICIT NONE
  REAL (KIND=8), INTENT( IN  ) :: rayx( 2 )              ! ray coordinate
  REAL (KIND=8), INTENT( IN  ) :: Topx( 2 ), Botx( 2 )   ! top, bottom coordinate
  REAL (KIND=8), INTENT( IN  ) :: Topn( 2 ), Botn( 2 )   ! top, bottom normal vector (outward)
  REAL (KIND=8), INTENT( OUT ) :: dTop( 2 ), dBot( 2 )   ! vector pointing from top, bottom bdry to ray
  REAL (KIND=8), INTENT( OUT ) :: DistTop, DistBot       ! distance (normal to bdry) from the ray to top, bottom boundary

  dTop    = rayx - Topx  ! vector pointing from top    to ray
  dBot    = rayx - Botx  ! vector pointing from bottom to ray
  DistTop = -DOT_PRODUCT( Topn, dTop )
  DistBot = -DOT_PRODUCT( Botn, dBot )

END SUBROUTINE Distances2D

! **********************************************************************!

SUBROUTINE Step2D( ray0, ray2, Topx, Topn, Botx, Botn )

  ! Does a single step along the ray
  ! x denotes the ray coordinate, (r,z)
  ! t denotes the scaled tangent to the ray (previously (rho, zeta))
  ! c * t would be the unit tangent

  USE bellhopMod
  USE BdryMod
  USE sspMod
  IMPLICIT NONE
  TYPE( ray2DPt )    :: ray0, ray1, ray2
  REAL (KIND=8 ), INTENT( IN ) :: Topx( 2 ), Topn( 2 ), Botx( 2 ), Botn( 2 )
  INTEGER            :: iSegz0, iSegr0
  REAL     (KIND=8 ) :: gradc0( 2 ), gradc1( 2 ), gradc2( 2 ), &
       c0, cimag0, crr0, crz0, czz0, csq0, cnn0_csq0, &
       c1, cimag1, crr1, crz1, czz1, csq1, cnn1_csq1, &
       c2, cimag2, crr2, crz2, czz2, urayt0( 2 ), urayt1( 2 ), &
       h, halfh, hw0, hw1, ray2n( 2 ), RM, RN, gradcjump( 2 ), cnjump, csjump, w0, w1, rho 

  ! The numerical integrator used here is a version of the polygon (a.k.a. midpoint, leapfrog, or Box method), and similar
  ! to the Heun (second order Runge-Kutta method).
  ! However, it's modified to allow for a dynamic step change, while preserving the second-order accuracy).

  ! *** Phase 1 (an Euler step)

  CALL EvaluateSSP( ray0%x, c0, cimag0, gradc0, crr0, crz0, czz0, rho, Freq, 'TAB' )

  csq0      = c0 * c0
  cnn0_csq0 = crr0 * ray0%t( 2 )**2 - 2.0 * crz0 * ray0%t( 1 ) * ray0%t( 2 ) + czz0 * ray0%t( 1 )**2
  iSegz0    = iSegz     ! make note of current layer
  iSegr0    = iSegr

  h = Beam%deltas       ! initially set the step h, to the basic one, deltas
  urayt0 = c0 * ray0%t  ! unit tangent

  CALL ReduceStep2D( ray0%x, urayt0, iSegz0, iSegr0, Topx, Topn, Botx, Botn, Beam%deltas, h ) ! reduce h to land on boundary
  halfh = 0.5 * h   ! first step of the modified polygon method is a half step

  ray1%x = ray0%x + halfh * urayt0
  ray1%t = ray0%t - halfh * gradc0 / csq0
  ray1%p = ray0%p - halfh * cnn0_csq0 * ray0%q
  ray1%q = ray0%q + halfh * c0        * ray0%p

  ! *** Phase 2

  CALL EvaluateSSP( ray1%x, c1, cimag1, gradc1, crr1, crz1, czz1, rho, Freq, 'TAB' )
  csq1      = c1 * c1
  cnn1_csq1 = crr1 * ray1%t( 2 )**2 - 2.0 * crz1 * ray1%t( 1 ) * ray1%t( 2 ) + czz1 * ray1%t( 1 )**2

  ! The Munk test case with a horizontally launched ray caused problems.
  ! The ray vertexes on an interface and can ping-pong around that interface.
  ! Have to be careful in that case about big changes to the stepsize (that invalidate the leap-frog scheme) in phase II.
  ! A modified Heun or Box method could also work.

  urayt1 = c1 * ray1%t   ! unit tangent

  CALL ReduceStep2D( ray0%x, urayt1, iSegz0, iSegr0, Topx, Topn, Botx, Botn, Beam%deltas, h ) ! reduce h to land on boundary

  ! use blend of f' based on proportion of a full step used.
  w1  = h / ( 2.0d0 * halfh )
  w0  = 1.0d0 - w1
  hw0 = h * w0
  hw1 = h * w1

  ray2%x   = ray0%x   + hw0 * urayt0              + hw1 * urayt1
  ray2%t   = ray0%t   - hw0 * gradc0 / csq0       - hw1 * gradc1 / csq1
  ray2%p   = ray0%p   - hw0 * cnn0_csq0 * ray0%q  - hw1 * cnn1_csq1 * ray1%q
  ray2%q   = ray0%q   + hw0 * c0        * ray0%p  + hw1 * c1        * ray1%p
  ray2%tau = ray0%tau + hw0 / CMPLX( c0, cimag0, KIND=8 ) + hw1 / CMPLX( c1, cimag1, KIND=8 )

  ray2%Amp       = ray0%Amp
  ray2%Phase     = ray0%Phase
  ray2%NumTopBnc = ray0%NumTopBnc
  ray2%NumBotBnc = ray0%NumBotBnc

  ! If we crossed an interface, apply jump condition

  CALL EvaluateSSP( ray2%x, c2, cimag2, gradc2, crr2, crz2, czz2, rho, Freq, 'TAB' )
  ray2%c = c2

  IF ( iSegz /= iSegz0 .OR. iSegr /= iSegr0 ) THEN
     gradcjump =  gradc2 - gradc0
     ray2n     = [ -ray2%t( 2 ), ray2%t( 1 ) ]   ! ray normal

     cnjump    = DOT_PRODUCT( gradcjump, ray2n  )
     csjump    = DOT_PRODUCT( gradcjump, ray2%t )

     IF ( iSegz /= iSegz0 ) THEN         ! crossing in depth
        RM = +ray2%t( 1 ) / ray2%t( 2 )  ! this is tan( alpha ) where alpha is the angle of incidence
     ELSE                                ! crossing in range
        RM = -ray2%t( 2 ) / ray2%t( 1 )  ! this is tan( alpha ) where alpha is the angle of incidence
     END IF

     RN     = RM * ( 2 * cnjump - RM * csjump ) / c2
     ray2%p = ray2%p - ray2%q * RN

  END IF

END SUBROUTINE Step2D

! **********************************************************************!

SUBROUTINE ReduceStep2D( x0, urayt, iSegz0, iSegr0, Topx, Topn, Botx, Botn, deltas, h )

  ! calculate a reduced step size, h, that lands on any points where the environment changes

  USE BdryMod
  USE sspMod
  IMPLICIT NONE
  INTEGER,       INTENT( IN  ) :: iSegz0, iSegr0                             ! SSP layer the ray is in
  REAL (KIND=8), INTENT( IN  ) :: x0( 2 ), urayt( 2 )                        ! ray coordinate and tangent
  REAL (KIND=8), INTENT( IN  ) :: Topx( 2 ), Topn( 2 ), Botx( 2 ), Botn( 2 ) ! Top, bottom coordinate and normal
  REAL (KIND=8), INTENT( IN  ) :: deltas                                     ! default step size
  REAL (KIND=8), INTENT( OUT ) :: h                                          ! reduced step size 
  REAL (KIND=8)                :: x( 2 ), d( 2 ), d0( 2 ), h1, h2, h3, h4, h5, h10

  ! Detect interface or boundary crossing and reduce step, if necessary, to land on that crossing.
  ! Keep in mind possibility that user put source right on an interface
  ! and that multiple events can occur (crossing interface, top, and bottom in a single step).

  x = x0 + h * urayt ! make a trial step

  ! interface crossing in depth
  h1 = huge( h1 )
  IF ( ABS( urayt( 2 ) ) > EPSILON( h1 ) ) THEN
     IF      ( SSP%z( iSegz0     ) > x(  2 ) ) THEN
        h1 = ( SSP%z( iSegz0     ) - x0( 2 ) ) / urayt( 2 )
     ELSE IF ( SSP%z( iSegz0 + 1 ) < x(  2 ) ) THEN
        h1 = ( SSP%z( iSegz0 + 1 ) - x0( 2 ) ) / urayt( 2 )
     END IF
  END IF

  ! top crossing
  h2 = huge( h2 )
  d  = x - Topx              ! vector from top to ray
  IF ( DOT_PRODUCT( Topn, d ) > EPSILON( h2 ) ) THEN
     d0  = x0 - Topx         ! vector from top    node to ray origin
     h2 = -DOT_PRODUCT( d0, Topn ) / DOT_PRODUCT( urayt, Topn )
  END IF

  ! bottom crossing
  h3 = huge( h3 )
  d  = x - Botx              ! vector from bottom to ray
  IF ( DOT_PRODUCT( Botn, d ) > EPSILON( h2 ) ) THEN
     d0  = x0 - Botx         ! vector from bottom node to ray origin
     h3 = -DOT_PRODUCT( d0, Botn ) / DOT_PRODUCT( urayt, Botn )
  END IF

  !!!! following tests could all be merged based on a single min or max range
  
  ! top segment crossing in range
  h4 = huge( h4 )
  IF ( ABS( urayt( 1 ) )  > EPSILON( h4 ) ) THEN
     IF       ( x(  1 ) < rTopSeg( 1 ) ) THEN
        h4 = -( x0( 1 ) - rTopSeg( 1 ) ) / urayt( 1 )
     ELSE IF  ( x(  1 ) > rTopSeg( 2 ) ) THEN
        h4 = -( x0( 1 ) - rTopSeg( 2 ) ) / urayt( 1 )
     END IF
  END IF

  ! bottom segment crossing in range
  h5 = huge( h5 )
  IF ( ABS( urayt( 1 ) )  > EPSILON( h5 ) ) THEN
     IF       ( x(  1 ) < rBotSeg( 1 ) ) THEN
        h5 = -( x0( 1 ) - rBotSeg( 1 ) ) / urayt( 1 )
     ELSE IF  ( x(  1 ) > rBotSeg( 2 ) ) THEN
        h5 = -( x0( 1 ) - rBotSeg( 2 ) ) / urayt( 1 )
     END IF
  END IF

  ! ocean segment crossing in r
  h10 = huge( h10 )

  IF ( SSP%Type == 'Q' ) THEN
     IF ( ABS( urayt( 1 ) ) > EPSILON( h1 ) ) THEN
        IF        ( x(  1 ) < SSP%Seg%r( iSegr0     ) ) THEN
           h10 = -( x0( 1 ) - SSP%Seg%r( iSegr0     ) ) / urayt( 1 )
           ! write( *, * ) 'ocean segment crossing in r'

        ELSE IF   ( x(  1 ) > SSP%Seg%r( iSegr0 + 1 ) ) THEN
           h10 = -( x0( 1 ) - SSP%Seg%r( iSegr0 + 1 ) ) / urayt( 1 )
           ! write( *, * ) 'ocean segment crossing in r'

        END IF
     END IF

  END IF

  h = MIN( h, h1, h2, h3, h4, h5, h10 )  ! take limit set by shortest distance to a crossing

  IF ( h < 1.0d-4 * deltas ) THEN   ! is it taking an infinitesimal step?
     h = 1.0d-5 * deltas            ! make sure we make some motion
     iSmallStepCtr = iSmallStepCtr + 1   ! keep a count of the number of sequential small steps
  ELSE
     iSmallStepCtr = 0   ! didn't do a small step so reset the counter
  END IF
  
END SUBROUTINE ReduceStep2D

! **********************************************************************!

SUBROUTINE Reflect2D( is, HS, BotTop, tBdry, nBdry, kappa, RefC, Npts )

  USE bellhopMod
  USE sspMod
  USE RefCoMod
  IMPLICIT NONE
  INTEGER,              INTENT( IN ) :: Npts
  REAL     (KIND=8),    INTENT( IN ) :: tBdry( 2 ), nBdry( 2 )  ! Tangent and normal to the boundary
  REAL     (KIND=8),    INTENT( IN ) :: kappa                   ! Boundary curvature
  CHARACTER (LEN=3),    INTENT( IN ) :: BotTop                  ! Flag indicating bottom or top reflection
  TYPE( HSInfo ),       INTENT( IN ) :: HS                      ! half-space properties
  TYPE(ReflectionCoef), INTENT( IN ) :: RefC( NPts )            ! reflection coefficient
  INTEGER,              INTENT( INOUT ) :: is
  INTEGER           :: is1
  REAL     (KIND=8) :: c, cimag, gradc( 2 ), crr, crz, czz, rho                   ! derivatives of sound speed
  REAL     (KIND=8) :: RM, RN, Tg, Th, rayt( 2 ), rayn( 2 ), rayt_tilde( 2 ), rayn_tilde( 2 ), cnjump, csjump  ! for curvature change
  REAL     (KIND=8) :: ck, co, si, cco, ssi, pdelta, rddelta, sddelta, theta_bot ! for beam shift
  COMPLEX  (KIND=8) :: gamma1, gamma2, gamma1Sq, gamma2Sq, GK, Refl   ! for tabulated reflection coef.
  COMPLEX  (KIND=8) :: ch, a, b, d, sb, delta, ddelta                 ! for beam shift
  TYPE(ReflectionCoef) :: RInt

  is  = is + 1
  is1 = is + 1

  Tg = DOT_PRODUCT( ray2D( is )%t, tBdry )  ! component of ray tangent, along boundary
  Th = DOT_PRODUCT( ray2D( is )%t, nBdry )  ! component of ray tangent, normal to boundary

  ray2D( is1 )%NumTopBnc = ray2D( is )%NumTopBnc
  ray2D( is1 )%NumBotBnc = ray2D( is )%NumBotBnc
  ray2D( is1 )%x         = ray2D( is )%x
  ray2D( is1 )%t         = ray2D( is )%t - 2.0 * Th * nBdry  ! changing the ray direction

  ! Calculate the change in curvature
  ! Based on formulas given by Muller, Geoph. J. R.A.S., 79 (1984).

  CALL EvaluateSSP( ray2D( is )%x, c, cimag, gradc, crr, crz, czz, rho, Freq, 'TAB' )   ! just to get c

  ! incident unit ray tangent and normal
  rayt = c * ray2D( is )%t           ! unit tangent to ray
  rayn = [ -rayt( 2 ), rayt( 1 ) ]   ! unit normal to ray

  ! reflected unit ray tangent and normal (the reflected tangent, normal system has a different orientation)
  rayt_tilde = c * ray2D( is1 )%t                       ! unit tangent to ray
  rayn_tilde = -[ -rayt_tilde( 2 ), rayt_tilde( 1 ) ]   ! unit normal to ray

  RN = -2 * kappa / c ** 2 / Th    ! boundary curvature correction

  ! get the jumps (this could be simplified, e.g. jump in rayt is roughly 2 * Th * nbdry
  cnjump = -DOT_PRODUCT( gradc, rayn_tilde - rayn  )
  csjump = -DOT_PRODUCT( gradc, rayt_tilde - rayt )

  IF ( BotTop == 'TOP' ) THEN
     cnjump = -cnjump   ! this is because the (t,n) system of the top boundary has a different sense to the bottom boundary
     RN     = -RN
  END IF

  RM = Tg / Th   ! this is tan( alpha ) where alpha is the angle of incidence
  RN = RN + RM * ( 2 * cnjump - RM * csjump ) / c ** 2

  SELECT CASE ( Beam%Type( 3 : 3 ) )
  CASE ( 'D' )
     RN = 2.0 * RN
  CASE ( 'Z' )
     RN = 0.0
  END SELECT

  ray2D( is1 )%c   = c
  ray2D( is1 )%tau = ray2D( is )%tau
  ray2D( is1 )%p   = ray2D( is )%p + ray2D( is )%q * RN
  ray2D( is1 )%q   = ray2D( is )%q

  ! account for phase change

  SELECT CASE ( HS%BC )
  CASE ( 'R' )                 ! rigid
     ray2D( is1 )%Amp   = ray2D( is )%Amp
     ray2D( is1 )%Phase = ray2D( is )%Phase
  CASE ( 'V' )                 ! vacuum
     ray2D( is1 )%Amp   = ray2D( is )%Amp
     ray2D( is1 )%Phase = ray2D( is )%Phase + pi
  CASE ( 'F' )                 ! file
     RInt%theta = RadDeg * ABS( ATAN2( Th, Tg ) )           ! angle of incidence (relative to normal to bathymetry)
     IF ( RInt%theta > 90 ) RInt%theta = 180. - RInt%theta  ! reflection coefficient is symmetric about 90 degrees
     CALL InterpolateReflectionCoefficient( RInt, RefC, Npts, PRTFile )
     ray2D( is1 )%Amp   = ray2D( is )%Amp * RInt%R
     ray2D( is1 )%Phase = ray2D( is )%Phase + RInt%phi
  CASE ( 'A', 'G' )            ! half-space
     GK       = omega * Tg     ! wavenumber in direction parallel to bathymetry
     gamma1Sq = ( omega / c     ) ** 2 - GK ** 2 - i * tiny( omega )   ! tiny prevents g95 giving -zero, and wrong branch cut
     gamma2Sq = ( omega / HS%cP ) ** 2 - GK ** 2 - i * tiny( omega )
     gamma1   = SQRT( -gamma1Sq )
     gamma2   = SQRT( -gamma2Sq )

     Refl = ( HS%rho * gamma1 - rho * gamma2 ) / ( HS%rho * gamma1 + rho * gamma2 )

     IF ( ABS( Refl ) < 1.0E-5 ) THEN   ! kill a ray that has lost its energy in reflection
        ray2D( is1 )%Amp   = 0.0
        ray2D( is1 )%Phase = ray2D( is )%Phase
     ELSE
        ray2D( is1 )%Amp   = ABS( Refl ) * ray2D(  is )%Amp
        ray2D( is1 )%Phase = ray2D( is )%Phase + ATAN2( AIMAG( Refl ), REAL( Refl ) )
        ! compute beam-displacement Tindle, Eq. (14)
        ! needs a correction to beam-width as well ...
        !  IF ( REAL( gamma2Sq ) < 0.0 ) THEN
        !     rhoW   = 1.0   ! density of water
        !     rhoWSq  = rhoW  * rhoW
        !     rhoHSSq = rhoHS * rhoHS
        !     DELTA = 2 * GK * rhoW * rhoHS * ( gamma1Sq - gamma2Sq ) /
        ! &( gamma1 * i * gamma2 *
        ! &( -rhoWSq * gamma2Sq + rhoHSSq * gamma1Sq ) )
        !     RV( is + 1 ) = RV( is + 1 ) + DELTA
        !  END IF

        if ( Beam%Type( 4 : 4 ) == 'S' ) then   ! beam displacement & width change (Seongil's version)

           ch = ray2D( is )%c / conjg( HS%cP )
           co = ray2D( is )%t( 1 ) * ray2D( is )%c
           si = ray2D( is )%t( 2 ) * ray2D( is )%c
           ck = omega / ray2D( is )%c

           a   = 2 * HS%rho * ( 1 - ch * ch )
           b   = co * co - ch * ch
           d   = HS%rho * HS%rho * si * si + b
           sb  = sqrt( b )
           cco = co * co
           ssi = si * si

           delta   = a * co / si / ( ck * sb * d )    
           pdelta  = real( delta ) / ( ray2D( is )%c / co)

           ddelta  = -a / ( ck*sb*d ) - a*cco / ssi / (ck*sb*d) + a*cco / (ck*b*sb*d) &
                -a*co / si / (ck*sb*d*d) * (2* HS%rho * HS%rho *si*co-2*co*si)
           rddelta = -real( ddelta )
           sddelta = rddelta / abs( rddelta )        

           ! next 3 lines have an update by Diana McCammon to allow a sloping bottom
           theta_bot = datan( tBdry( 2 ) / tBdry( 1 ))  ! bottom angle
           ray2D( is1 )%x( 1 ) = ray2D( is1 )%x( 1 ) + real( delta ) *dcos( theta_bot )   ! range displacement
           ray2D( is1 )%x( 2 ) = ray2D( is1 )%x( 2 ) + real( delta ) *dsin( theta_bot )   ! depth displacement
           ray2D( is1 )%tau    = ray2D( is1 )%tau + pdelta             ! phase change
           ray2D( is1 )%q      = ray2D( is1 )%q + sddelta * rddelta * si * c * ray2D( is )%p   ! beam-width change
        endif

     ENDIF
  CASE DEFAULT
     WRITE( PRTFile, * ) 'HS%BC = ', HS%BC
     CALL ERROUT( PRTFile, 'F', 'Reflect2D', 'Unknown boundary condition type' )
  END SELECT

END SUBROUTINE Reflect2D

!*******************************************************************!

SUBROUTINE WriteRay( alpha0, Nsteps1 )

  ! Compress the ray data keeping every iSkip point, points near surface or bottom, and last point.
  ! Write to RAYFile.
  ! During an eigenray calculation, subsets of the full ray may be passed
  ! These have lengths Nsteps1 vs. Nsteps for the entire ray

  USE bellhopMod
  IMPLICIT NONE
  INTEGER,       PARAMETER    :: RAYFile = 21
  INTEGER,       INTENT( IN ) :: Nsteps1
  REAL (KIND=8), INTENT( IN ) :: alpha0   ! take-off angle of this ray
  INTEGER                     :: N2, iSkip

  ! compression

  N2    = 1
  iSkip = MAX( Nsteps1 / 5000, 1 )   ! max #pts written is about 5000

  Stepping: DO iStep = 2, Nsteps1
     ! ensure that we always write ray pts. near bdry reflections (works only for flat bdry)
     IF ( MIN( Bdry%Bot%HS%Depth - ray2D( iStep )%x( 2 ),  ray2D( iStep )%x( 2 ) - Bdry%Top%HS%Depth ) < 0.2 .OR. &
          MOD( iStep, iSkip ) == 0 .OR. iStep == Nsteps1 ) THEN
        N2 = N2 + 1
        ray2D( N2 )%x = ray2D( iStep )%x
     END IF
  END DO Stepping

  ! write to ray file

  WRITE( RAYFile, * ) alpha0
  WRITE( RAYFile, * ) N2, ray2D( Nsteps1 )%NumTopBnc, ray2D( Nsteps1 )%NumBotBnc
  DO iStep = 1, N2
     WRITE( RAYFile, * ) SNGL( ray2D( iStep )%x )
  END DO

END SUBROUTINE WriteRay

