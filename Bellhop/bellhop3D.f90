PROGRAM BELLHOP3D

  ! Gaussian beam tracing in three dimensions
  ! Michael B. Porter

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

  ! Original version done at NRL in 1986
  ! Published at SACLANT Undersea Research Center 1989
  !
  ! Converted to Fortran 2003 and greatly enhanced for outside users in 2010
  ! with the support of the Agency for Defense Development, Republic of Korea
  ! Supported also by the U.S. Office of Naval Research
  !
  ! Changes included:
  !    Use of standard BELLHOP input files
  !    Inclusion of multiple sources and receivers
  !    Integrated option for Nx2D or full 3D
  !    Implementation of several standard beam options for both the Nx2D and 3D cases
  !    Reading in:
  !       oceanography from a SSP file
  !       bathymetry or altimetry  from a BTY or ATI file
  !       a reflection coefficient from a BRC or TRC file
  !       a source beam pattern from a SBP file

  ! Loose ends:
  !    Nx2D version does not handle jumps in cx, cy (routine step2d)
  !    cannot specify isingle( 2 ) for alpha and beta
  !    Trilinear interpolation (hex) ignores cxy values.
  !    If the water depth for the lower half space is much larger than that for the SSP
  !    in the water column, it will use a step that is too large and exit the ray box.
 
  !    Cerveny beams (rarely used):
  !       influenceC calls SSP; need to select EvaluateSSP2D or EvaluateSSP3D for that to work in BELLHOP3D
  !       efficiency changes for Cerveny beams as well (and in BELLHOP)
  !       space filling or minimum width options are applied to both alpha and beta--- often only want space filling in azimuth

  !    Influend3DGeoHat writes no eigenray info if number of receiver ranges NR=1
  !    r( 1 ) = 1 m in BELLHOP plus logic for reversing
  !    geogaussian should pre-calculate dtauds and dqds like geohat?
  !    fix automatic deltas selection
  !    detect sd below bottom with bty file

  !    Desired additional features:
  !       Terrain following option for RD
  !       Variable bottom type vs.lat/long

  USE bellhopMod
  USE RefCoMod
  USE bdry3DMod
  USE BeamPatternMod

  IMPLICIT NONE
  LOGICAL,   PARAMETER :: ThreeD = .TRUE.
  CHARACTER ( LEN=80 ) :: FileRoot

  ! get the file root for naming all input and output files
  ! should add some checks here ...

  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )
  ! Open the print file
  OPEN( UNIT = PRTFile, FILE = TRIM( FileRoot ) // '.prt', STATUS = 'UNKNOWN', IOSTAT = iostat )

  ! Read in control data

  CALL ReadEnvironment(    FileRoot, ThreeD )  
  CALL ReadATI3D( FileRoot, Bdry%Top%HS%Opt( 5 : 5 ), Bdry%Top%HS%Depth, PRTFile )    ! AlTImetry
  CALL ReadBTY3D( FileRoot, Bdry%Bot%HS%Opt( 2 : 2 ), Bdry%Bot%HS%Depth, PRTFile )    ! BaThYmetry

  CALL ReadReflectionCoefficient(  FileRoot, Bdry%Bot%HS%Opt( 1 : 1 ), Bdry%Top%HS%Opt( 2 : 2 ), PRTFile )    ! (top and bottom)
  SBPFlag = Beam%RunType( 3 : 3 )
  CALL ReadPAT( FileRoot,                                 PRTFile )    ! Source Beam Pattern
  CALL OpenOutputFiles( FileRoot, ThreeD )

  CALL BellhopCore

END PROGRAM BELLHOP3D

! **********************************************************************!

SUBROUTINE BellhopCore

  USE bellhopMod
  USE RefCoMod
  USE bdry3DMod
  USE angleMod
  USE SdRdRMod
  USE ArrMod
  USE BeamPatternMod
  USE sspMod

  IMPLICIT NONE
  LOGICAL,   PARAMETER :: ThreeD = .TRUE.
  INTEGER,   PARAMETER :: SHDFile = 25, RAYFile = 21, ArrivalsStorage = 2000000
  INTEGER              :: IBPvec( 1 ), ibp, iBeamWindow2, ird, isd, itheta, isx, isy, iRec, ir
  REAL      ( KIND=8 ) :: Tstart, Tstop
  REAL      ( KIND=8 ) :: Amp0, alpha0, RadMax, S
  REAL      ( KIND=8 ) :: c0, cimag0, gradc( 3 ), cxx, cyy, czz, cxy, cxz, cyz, rho, xs( 3 )
  REAL      ( KIND=8 ), ALLOCATABLE :: x_rcvrMat( :, :, : ), t_rcvr( :, : )
  COMPLEX   ( KIND=8 ) :: epsilon( 2 )
  COMPLEX, ALLOCATABLE :: P( :, :, : ), U( :, : )

  CALL CPU_TIME( Tstart )

  omega = 2.0 * pi * freq

  IF ( Beam%deltas == 0.0 ) Beam%deltas = ( Bdry%Bot%HS%Depth - Bdry%Top%HS%Depth ) / 10.0   ! Automatic step size selection

  Angles%alpha  = DegRad * Angles%alpha   ! convert to radians
  Angles%Dalpha = 0.0
  IF ( Angles%Nalpha /= 1 ) &
     Angles%Dalpha = ( Angles%alpha( Angles%Nalpha ) - Angles%alpha( 1 ) ) / ( Angles%Nalpha - 1 )  ! angular spacing between beams
  SELECT CASE ( Beam%RunType( 5 : 5 ) )
  CASE ( 'I' )
     Nrd_per_range = 1         ! irregular grid
  CASE ( 'R' )
     Nrd_per_range = Pos%Nrd   ! rectilinear grid
  END SELECT

  ! for a TL calculation, allocate space for the pressure matrix
  SELECT CASE ( Beam%RunType( 1 : 1 ) )
  CASE ( 'C', 'S', 'I' )        ! TL calculation
     ALLOCATE ( P( Pos%Ntheta, Nrd_per_range, Pos%Nr ), Stat = IAllocStat )
     ALLOCATE ( U( Nrd_per_range, Pos%Nr ), Stat = IAllocStat )    ! used for 2D option
     IF ( IAllocStat /= 0 ) &
          CALL ERROUT( PRTFile, 'F', 'BELLHOP', 'Insufficient memory for TL matrix: reduce Nr * Nrd'  )
  CASE ( 'A', 'a', 'R', 'E' )   ! Arrivals calculation
     ALLOCATE ( P( 1, 1, 1 ), Stat = IAllocStat )   ! open a dummy variable
     ALLOCATE ( U( 1, 1 ),    Stat = IAllocStat )   ! open a dummy variable
  END SELECT

  ! for an arrivals run, allocate space for arrivals matrices
  SELECT CASE ( Beam%RunType( 1 : 1 ) )
  CASE ( 'A', 'a' )
     MaxNArr = MAX( ArrivalsStorage / ( Pos%Ntheta * Nrd_per_range * Pos%Nr ), 10 )   ! allow space for at least 10 arrivals
     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) '( Maximum # of arrivals = ', MaxNArr, ')'

     ALLOCATE ( Arr3D( Pos%Ntheta, Nrd_per_range, Pos%Nr, MaxNArr ), &
               NArr3D( Pos%Ntheta, Nrd_per_range, Pos%Nr ), Stat = IAllocStat )
     IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'BELLHOP', &
          'Insufficient memory to allocate arrivals matrix; reduce parameter ArrivalsStorage' )

     ! For a 2D Arrivals run, also need to allocate a structure for that
     IF ( Beam%RunType( 6 : 6 ) == '2' ) THEN
         ALLOCATE ( Arr( Nrd_per_range, Pos%Nr, MaxNArr ), &
                   NArr( Nrd_per_range, Pos%Nr ), Stat = IAllocStat )
         IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'BELLHOP', &
            'Insufficient memory to allocate arrivals matrix; reduce parameter ArrivalsStorage' )
         NArr( 1 : Nrd_per_range, 1 : Pos%Nr ) = 0
     END IF
  CASE DEFAULT
     MaxNArr = 1
     ALLOCATE ( Arr3D( Pos%Ntheta, Nrd_per_range, Pos%Nr, 1 ), NArr3D( Pos%Ntheta, Nrd_per_range, Pos%Nr ), Stat = IAllocStat )
  END SELECT

  NArr3D( 1 : Pos%Ntheta, 1 : Nrd_per_range, 1 : Pos%Nr ) = 0
  
  ALLOCATE( x_rcvrMat( 2, Pos%Ntheta, Pos%Nr ), t_rcvr( 2, Pos%Ntheta ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'BELLHOP', &
       'Insufficient memory to x_rcvrMat; reduce Nr or Ntheta' )

  ! tangent along receiver bearing line
  t_rcvr( 1, : ) = COS( DegRad * Pos%theta( 1 : Pos%Ntheta ) )
  t_rcvr( 2, : ) = SIN( DegRad * Pos%theta( 1 : Pos%Ntheta ) )

  WRITE( PRTFile, * )

  SourceDepth: DO isd = 1, Pos%Nsd   ! Loop over source depth

     Source_x: DO isx = 1, Pos%Nsx    ! loop over source x-coordinate

        Source_y: DO isy = 1, Pos%Nsy   ! loop over source y-coordinate
           P = 0.0 !  Zero out field matrix
           U = 0.0
           
           ! IF ( r( 1 ) == 0.0 ) r( 1 ) = 1.0
           xs = [ Pos%sx( isx ), Pos%sy( isy ), Pos%sd( isd ) ]
           WRITE( PRTFile, * )
           WRITE( PRTFile, * ) 'xs = ', xs

           ! positions of rcvrs in the x-y plane; this is pre-calculated for InfluenceGeoHatCart
           ! It is not clear that the pre-calculation saves time ...
           DO ir = 1, Pos%Nr
              DO itheta = 1, Pos%Ntheta
                 x_rcvrMat( 1 : 2, itheta, ir ) = xs( 1 : 2 ) + Pos%r( ir ) * t_rcvr( :, itheta )  ! x-y coordinate of the receiver
              END DO
           END DO
  
           ! *** Compute 'optimal' beam constant ***

           CALL EvaluateSSP3D( xs, c0, cimag0, gradc, cxx, cyy, czz, cxy, cxz, cyz, rho, freq, 'TAB' )
           ray2D( 1 )%c = c0
           CALL PickEpsilon( Beam%Type( 1 : 2 ), omega, c0, Angles%Dalpha, Angles%Dbeta, Beam%rLoop, Beam%epsMultiplier, epsilon ) ! beam constant

           ! *** Trace successive beams ***

           AzimuthalAngle: DO ibeta = 1, Angles%Nbeta   ! this is also the receiver bearing angle for a 2D run
              IF ( Angles%iSingle == 0 .OR. ibeta == 2 ) THEN    ! Single beam run?
              !IF ( ibeta == 2 ) THEN    ! Single beam run?
              !IF ( mod( ibeta+1, 2 ) == 0 ) THEN    ! Single beam run?
                 WRITE( PRTFile, FMT = "( 'Tracing azimuthal beam ', I4, F10.2 )" ) ibeta, RadDeg * Angles%beta( ibeta )
                 ! WRITE( *,       FMT = "( 'Tracing beam ', I4, F10.2 )" ) ibeta, RadDeg * Angles%beta( ibeta )

                 ElevationAngle: DO ialpha = 1, Angles%Nalpha

                    IF ( Angles%iSingle == 0 .OR. ialpha == Angles%iSingle ) THEN    ! Single beam run?
                    !IF ( ialpha  == 11 .or. ialpha == 12 ) THEN    ! Single beam run?

                       alpha0 = RadDeg * Angles%alpha( ialpha )    ! take-off angle in degrees
                       ! WRITE( PRTFile, FMT = "( '   Tracing declination beam ', I4, F10.2 )" ) ialpha, alpha0
                       !write( *, * ) '    ialpha', ialpha, alpha0

                       IBPvec = maxloc( SrcBmPat( :, 1 ), mask = SrcBmPat( :, 1 ) < alpha0 )  ! index of ray angle in beam pattern
                       IBP    = IBPvec( 1 )
                       IBP    = MAX( IBP, 1 )               ! don't go before beginning of table
                       IBP    = MIN( IBP, NSBPPts - 1 )     ! don't go past end of table

                       ! linear interpolation to get amplitude
                       s    = ( alpha0  - SrcBmPat( IBP, 1 ) ) / ( SrcBmPat( IBP + 1, 1 ) - SrcBmPat( IBP, 1 ) )
                       Amp0 = ( 1 - s ) * SrcBmPat( IBP, 2 ) + s * SrcBmPat( IBP + 1, 2 )

                       SELECT CASE ( Beam%RunType( 6 : 6 ) )   ! flag for 2D or 3D calculation
                       CASE ( '2' )   ! Nx2D calculation, neglecting horizontal refraction
                          CALL TraceRay2D(  xs, Angles%alpha( ialpha ), Angles%beta( ibeta ), Amp0 )
                          IF ( Beam%RunType( 1 : 1 ) /= 'R' ) THEN     ! If not a ray trace run, calculate the field
                             SELECT CASE ( Beam%Type( 1 : 1 ) )
                             CASE ( 'R' )
                                iBeamWindow2 = Beam%iBeamWindow ** 2
                                RadMax       = 50 * c0 / freq  ! 50 wavelength max radius
                                CALL InfluenceCervenyRayCen(   U, epsilon( 1 ), Angles%alpha( ialpha ), IBeamWindow2, RadMax )
                             CASE ( 'C' )
                                CALL ERROUT( PRTFile, 'F', 'BELLHOP3D', 'Run Type ''C'' not supported at this time' )
                                iBeamWindow2 = Beam%iBeamWindow ** 2
                                RadMax       = 50 * c0 / freq  ! 50 wavelength max radius
                                CALL InfluenceCervenyCart(     U, epsilon( 1 ), Angles%alpha( ialpha ), IBeamWindow2, RadMax )
                             CASE ( 'g' )
                                CALL InfluenceGeoHatRayCen(    U,       Angles%alpha( ialpha ), Angles%Dalpha )
                             CASE ( 'S' )
                                CALL InfluenceSGB(             U,       Angles%alpha( ialpha ), Angles%Dalpha )
                             CASE ( 'B' )
                                CALL InfluenceGeoGaussianCart( U,       Angles%alpha( ialpha ), Angles%Dalpha )
                             CASE DEFAULT
                                CALL InfluenceGeoHatCart(      U,       Angles%alpha( ialpha ), Angles%Dalpha )
                             END SELECT
                          END IF

                       CASE ( '3' )   ! full 3D calculation
                          CALL TraceRay3D( xs, Angles%alpha( ialpha ), Angles%beta( ibeta ), epsilon, Amp0 )

                          IF ( Beam%RunType( 1 : 1 ) /= 'R' ) THEN     ! If not a ray trace run, calculate the field

                             SELECT CASE ( Beam%RunType( 2 : 2 ) )
                             CASE ( 'C' )   ! Cerveny style beams
                                CALL ERROUT( PRTFile, 'F', 'BELLHOP3D', 'Run Type ''C'' not supported at this time' )

                                ! option: assemble f, g, h from p-q
                                !ray3D( 1 : Beam%Nsteps )%f    =    ray3D( 1 : Beam%Nsteps )%p_tilde( 1 ) * &
                                !                                   ray3D( 1 : Beam%Nsteps )%q_hat( 2 ) - &
                                !                                   ray3D( 1 : Beam%Nsteps )%q_tilde( 2 ) * &
                                !                                   ray3D( 1 : Beam%Nsteps )%p_hat( 1 )
                                !ray3D( 1 : Beam%Nsteps )%g    = -( ray3D( 1 : Beam%Nsteps )%p_tilde( 2 ) * &
                                !                                   ray3D( 1 : Beam%Nsteps )%q_hat( 1 ) - &
                                !                                   ray3D( 1 : Beam%Nsteps )%q_tilde( 1 ) * &
                                !                                   ray3D( 1 : Beam%Nsteps )%p_hat( 2 ) )
                                !ray3D( 1 : Beam%Nsteps )%h    =    ray3D( 1 : Beam%Nsteps )%p_tilde( 2 ) * &
                                !                                   ray3D( 1 : Beam%Nsteps )%q_hat( 2 ) - &
                                !                                   ray3D( 1 : Beam%Nsteps )%q_tilde( 2 ) * &
                                !                                   ray3D( 1 : Beam%Nsteps )%p_hat( 2 )
                                !ray3D( 1 : Beam%Nsteps )%DetP =    ray3D( 1 : Beam%Nsteps )%p_tilde( 1 ) * &
                                !                                   ray3D( 1 : Beam%Nsteps )%p_hat( 2 ) - &
                                !                                   ray3D( 1 : Beam%Nsteps )%p_tilde( 2 ) * &
                                !                                   ray3D( 1 : Beam%Nsteps )%p_hat( 1 )
                                !ray3D( 1 : Beam%Nsteps )%DetQ =    ray3D( 1 : Beam%Nsteps )%q_tilde( 1 ) * &
                                !                                   ray3D( 1 : Beam%Nsteps )%q_hat( 2 ) - &
                                !                                   ray3D( 1 : Beam%Nsteps )%q_tilde( 2 ) * &
                                !                                   ray3D( 1 : Beam%Nsteps )%q_hat( 1 )

                                !CALL Influence3D( xs, Angles%alpha( ialpha ), iBeamWindow, P )
                             CASE ( 'g' )        ! Geometric-beams with hat-shape
                                CALL Influence3DGeoHatRayCen(       xs, Angles%alpha( ialpha ), Angles%beta( ibeta ), &
                                                              Angles%Dalpha, Angles%Dbeta, P )
                             CASE ( 'G', '^' )   ! Geometric-beams with hat-shape in Cartesian coordinates
                                CALL Influence3DGeoHatCart(         xs, Angles%alpha( ialpha ), Angles%beta( ibeta ), &
                                                              Angles%Dalpha, Angles%Dbeta, P, x_rcvrMat, t_rcvr )
                             CASE ( 'b' )        ! Geometric-beams with Gaussian-shape
                                CALL Influence3DGeoGaussianRayCen(  xs, Angles%alpha( ialpha ), Angles%beta( ibeta ), &
                                                              Angles%Dalpha, Angles%Dbeta, P )
                             CASE ( 'B' )        ! Geometric-beams with Gaussian-shape
                                CALL Influence3DGeoGaussianCart(  xs, Angles%alpha( ialpha ), Angles%beta( ibeta ), &
                                                              Angles%Dalpha, Angles%Dbeta, P, x_rcvrMat, t_rcvr )
                            CASE DEFAULT
                                CALL ERROUT( PRTFile, 'F', 'BELLHOP3D', 'Invalid Run Type' )
                             END SELECT
                          END IF
                       END SELECT

                       ! Optionally dump rays to a disk file
                       IF ( Beam%RunType( 1 : 1 ) == 'R' ) THEN
                          CALL WriteRay3D( Angles%alpha( ialpha ), Angles%beta( ibeta ), Beam%Nsteps, xs )
                       ENDIF
                    ENDIF   ! closes iSingle test
                 END DO ElevationAngle

                 ! for a 2D TL run, scale the pressure and copy the 2D slice into the radial of the 3D field
                 IF ( Beam%RunType( 6 : 6 ) == '2' ) THEN  ! 2D calculation
                    SELECT CASE ( Beam%RunType( 1 : 1 ) )
                    CASE ( 'C', 'S', 'I' )   ! TL calculation
                       CALL ScalePressure( Angles%Dalpha, ray2D( 1 )%c, Pos%r, U, Nrd_per_range, Pos%Nr, Beam%RunType, freq )
                       P( ibeta, :, : ) = U
                       U = 0   ! clear out the pressure field on the radial in prep for next radial
                    CASE ( 'A' )             ! arrivals calculation, ascii
                       NArr3D( ibeta, :, :    ) = NArr( :, : )
                       Arr3D(  ibeta, :, :, : ) = Arr(  :, :, : )
                       Narr = 0   ! this clears out the 2D arrival structure
                    CASE ( 'a' )             ! arrivals calculation, binary
                       NArr3D( ibeta, :, :    ) = NArr( :, : )
                       Arr3D(  ibeta, :, :, : ) = Arr(  :, :, : )
                       Narr = 0   ! this clears out the 2D arrival structure
                    END SELECT
                 END IF
              ENDIF   ! closes iSingle test
           END DO AzimuthalAngle

           ! *** Scale the complex pressure field ***

           IF ( Beam%RunType( 6 : 6 ) == '3' ) THEN  ! 3D calculation
              SELECT CASE ( Beam%RunType( 1 : 1 ) )
              CASE ( 'C', 'S', 'I' )   ! TL calculation
                 CALL ScalePressure3D( Angles%Dalpha, Angles%Dbeta, ray2D( 1 )%c, epsilon, P, &
                                       Pos%Ntheta, Nrd_per_range, Pos%Nr, Beam%RunType, freq )
              END SELECT
           END IF

           ! Write out the field

           SELECT CASE ( Beam%RunType( 1 : 1 ) )
           CASE ( 'C', 'S', 'I' )   ! TL calculation
              DO ird = 1, Pos%Nrd
                 RcvrBearing: DO itheta = 1, Pos%Ntheta
                    iRec = 10 + ( isx    - 1 ) * Pos%Nsy * Pos%Ntheta * Pos%Nsd * Pos%Nrd + &
                                ( isy    - 1 ) *           Pos%Ntheta * Pos%Nsd * Pos%Nrd + &
                                ( itheta - 1 ) *                        Pos%Nsd * Pos%Nrd + &
                                ( isd    - 1 ) *                                  Pos%Nrd + ird
                    WRITE( SHDFile, REC = IRec ) P( itheta, ird, 1 : Pos%Nr )

                 END DO RcvrBearing
              END DO
              CASE ( 'A' )             ! arrivals calculation, ascii
                 CALL WriteArrivalsASCII3D(  Pos%r, Pos%Ntheta, Nrd_per_range, Pos%Nr, Beam%RunType( 4:4 ) )
              CASE ( 'a' )             ! arrivals calculation, binary
                 CALL WriteArrivalsBinary3D( Pos%r, Pos%Ntheta, Nrd_per_range, Pos%Nr, Beam%RunType( 4:4 ) )
           END SELECT
        END DO Source_y
     END DO Source_x
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

SUBROUTINE PickEpsilon( BeamType, omega, c, Dalpha, Dbeta, rLoop, EpsMultiplier, epsilon )

  ! Picks the optimum value for epsilon

  IMPLICIT NONE
  INTEGER,            PARAMETER     :: PRTFile = 6
  COMPLEX,            PARAMETER     :: i = ( 0.0, 1.0 )
  REAL      (KIND=8), INTENT( IN  ) :: omega, c             ! angular frequency, sound speed
  REAL      (KIND=8), INTENT( IN  ) :: Dalpha, Dbeta        ! angular spacing for ray fan
  REAL      (KIND=8), INTENT( IN  ) :: epsMultiplier, Rloop ! multiplier, loop range
  COMPLEX   (KIND=8), INTENT( OUT ) :: epsilon( 2 )         ! beam initial conditions
  CHARACTER (LEN= 2), INTENT( IN  ) :: BeamType
  LOGICAL, SAVE      :: INIFlag = .TRUE.
  REAL      (KIND=8) :: HalfWidth( 2 )
  COMPLEX   (KIND=8) :: epsilonOpt( 2 )
  CHARACTER (LEN=80) :: TAG

  SELECT CASE ( BeamType( 1 : 1 ) )
  CASE ( 'C', 'R' )   ! Cerveny beams
     TAG    = 'Cerveny style beam'

     SELECT CASE ( BeamType( 2 : 2 ) )
     CASE ( 'F' )
        TAG            = 'Space filling beams'
        HalfWidth( 1 ) = 2.0 / ( ( omega / c ) * Dalpha )
        HalfWidth( 2 ) = 0.0
        IF ( Dbeta /= 0.0 ) HalfWidth( 2 ) = 2.0 / ( ( omega / c ) * Dbeta )
        epsilonOpt     = i * 0.5 * omega * HalfWidth ** 2
     CASE ( 'M' )
        TAG             = 'Minimum width beams'
        HalfWidth(  1 ) = SQRT( 2.0 * c * 1000.0 * rLoop / omega )
        HalfWidth(  2 ) = HalfWidth( 1 )
        epsilonOPT      = i * 0.5 * omega * HalfWidth **2
     CASE ( 'C' )
        TAG    = 'Cerveny style beam'
     END SELECT

  CASE ( 'g' )
     TAG             = 'Geometric beam, hat-shaped, Ray coord.'
     HalfWidth(  1 ) = 0.0
     HalfWidth(  2 ) = 0.0
     epsilonOPT      = 0.0

  CASE ( 'G', '^' )
     TAG             = 'Geometric beam, hat-shaped, Cart. coord.'
     HalfWidth(  1 ) = 0.0
     HalfWidth(  2 ) = 0.0
     epsilonOPT      = 0.0

  CASE ( 'b' )
     TAG             = 'Geometric beam, Gaussian-shaped, Ray coord.'
     HalfWidth(  1 ) = 0.0
     HalfWidth(  2 ) = 0.0
     epsilonOPT      = 0.0

  CASE ( 'B' )
     TAG             = 'Geometric beam, Gaussian-shaped, Cart. coord.'
     HalfWidth(  1 ) = 0.0
     HalfWidth(  2 ) = 0.0
     epsilonOPT      = 0.0

  CASE ( 'S' )
     TAG        = 'Simple Gaussian beams'
     halfwidth  = 2.0 / ( ( omega / c ) * Dalpha )
     epsilonOpt = i * 0.5 * omega * halfwidth ** 2
  END SELECT

  epsilon = epsMultiplier * epsilonOPT

  ! On first call write info to prt file
  IF ( INIFlag ) THEN
     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) TAG
     WRITE( PRTFile, * ) 'HalfWidth1  = ', HalfWidth( 1 )
     WRITE( PRTFile, * ) 'HalfWidth2  = ', HalfWidth( 2 )
     WRITE( PRTFile, * ) 'epsilonOPT1 = ', epsilonOPT( 1 )
     WRITE( PRTFile, * ) 'epsilonOPT2 = ', epsilonOPT( 2 )
     WRITE( PRTFile, * ) 'EpsMult     = ', EpsMultiplier
     WRITE( PRTFile, * )
     INIFlag = .FALSE.
  END IF

END SUBROUTINE PickEpsilon


!**********************************************************************!

SUBROUTINE TraceRay2D( xs, alpha, beta, Amp0 )

  ! Traces the beam corresponding to a particular take-off angle

  USE bellhopMod
  USE bdry3DMod
  USE RefCoMod
  IMPLICIT NONE
  REAL     (KIND=8), INTENT( IN ) :: alpha, beta, Amp0 ! initial angles, amplitude
  REAL     (KIND=8), INTENT( IN ) :: xs( 3 )     ! x-y-z coordinate of the source
  INTEGER           :: is, is1                   ! index for a step along the ray
  REAL     (KIND=8) :: x( 3 )                    ! ray coordinate
  REAL     (KIND=8) :: c, cimag, gradc( 2 ), crr, crz, czz, rho
  REAL     (KIND=8) :: DistBegTop, DistEndTop, DistBegBot, DistEndBot ! Distances from ray beginning, end to top and bottom
  REAL     (KIND=8) :: tradial( 2 ), BotnInt( 3 ), TopnInt( 3 ), s1, s2

  ! *** Initial conditions ***

  iSmallStepCtr = 0
  tradial = [ COS( beta ), SIN( beta ) ]
  ray2D( 1 )%x = [ 0.0D0, xs( 3 ) ]

  CALL EvaluateSSP2D( ray2D( 1 )%x, c, cimag, gradc, crr, crz, czz, rho, xs, tradial, freq )
  ray2D( 1 )%t         = [ COS( alpha ), SIN( alpha ) ] / c
  ray2D( 1 )%p         = [ 1.0, 0.0 ]
  ray2D( 1 )%q         = [ 0.0, 1.0 ]
  ray2D( 1 )%tau       = 0.0
  ray2D( 1 )%Amp       = Amp0
  ray2D( 1 )%Phase     = 0.0
  ray2D( 1 )%NumTopBnc = 0
  ray2D( 1 )%NumBotBnc = 0

  ! *** Trace the beam ***

  CALL GetTopSeg3D( xs )   ! identify the top    segment above the source
  CALL GetBotSeg3D( xs )   ! identify the bottom segment below the source

  ! Trace the beam (note that Reflect alters the step index is)
  is = 0
  CALL Distances3D( xs, Topx, Botx, Topn, Botn, DistBegTop, DistBegBot )

  IF ( DistBegTop <= 0 .OR. DistBegBot <= 0 ) THEN
     Beam%Nsteps = 1
     RETURN       ! source must be within the medium
  END IF

  Stepping: DO istep = 1, MaxN - 1
     is  = is + 1
     is1 = is + 1
     CALL Step2D( ray2D( is ), ray2D( is1 ), xs, tradial )

     ! convert polar coordinate of ray to x-y coordinate
     x( 1 ) = xs( 1 ) + ray2D( is1 )%x( 1 ) * tradial( 1 )
     x( 2 ) = xs( 2 ) + ray2D( is1 )%x( 1 ) * tradial( 2 )
     x( 3 ) = ray2D( is1 )%x( 2 )

     CALL GetTopSeg3D( x )    ! identify the top    segment above the source
     CALL GetBotSeg3D( x )    ! identify the bottom segment below the source

     IF ( IsegTopx == 0 .OR. IsegTopy == 0 .OR. IsegBotx == 0 .OR. IsegBoty == 0 ) THEN ! we escaped the box
        Beam%Nsteps = is
        EXIT Stepping
     END IF

     ! Reflections?
     ! Tests that ray at step is is inside, and ray at step is+1 is outside
     ! to detect only a crossing from inside to outside
     ! DistBeg is the distance at step is,   which is saved
     ! DistEnd is the distance at step is+1, which needs to be calculated
  
     CALL Distances3D( x, Topx, Botx, Topn, Botn, DistEndTop, DistEndBot )

     IF      ( DistBegTop > 0.0d0 .AND. DistEndTop <= 0.0d0 ) THEN  ! test top reflection
        IF ( atiType == 'C' ) THEN

           x = [ xs( 1 ) + ray2D( is + 1 )%x( 1 ) * tradial( 1 ),   &
                 xs( 2 ) + ray2D( is + 1 )%x( 1 ) * tradial( 2 ),   &
                           ray2D( is + 1 )%x( 2 ) ]

           s1     = ( x( 1 ) - Topx( 1 ) ) / Top_deltax   ! proportional distance along segment
           s2     = ( x( 2 ) - Topx( 2 ) ) / Top_deltay   ! proportional distance along segment

           TopnInt = Top( IsegTopx,     IsegTopy     )%Noden * ( 1 - s1 ) * ( 1 - s2 ) +  &
                     Top( IsegTopx + 1, IsegTopy     )%Noden * ( s1     ) * ( 1 - s2 ) +  &
                     Top( IsegTopx + 1, IsegTopy + 1 )%Noden * ( s1     ) * ( s2     ) +  &
                     Top( IsegTopx,     IsegTopy + 1 )%Noden * ( 1 - s1 ) * ( s2     )
        ELSE
           TopnInt = Topn   ! normal is constant in a segment
        END IF
        
        CALL Reflect2D( is, Bdry%Top%HS, 'TOP', TopnInt, RTop, NTopPTS, xs, tradial )
        ray2D( is + 1 )%NumTopBnc = ray2D( is )%NumTopBnc + 1

        x = [ xs( 1 ) + ray2D( is + 1 )%x( 1 ) * tradial( 1 ),   &
              xs( 2 ) + ray2D( is + 1 )%x( 1 ) * tradial( 2 ),   &
                        ray2D( is + 1 )%x( 2 ) ]

        CALL Distances3D( x, Topx, Botx, Topn, Botn, DistEndTop, DistEndBot )

     ELSE IF ( DistBegBot > 0.0d0 .AND. DistEndBot <= 0.0d0 ) THEN  ! test bottom reflection 
        IF ( btyType == 'C' ) THEN

           x = [ xs( 1 ) + ray2D( is + 1 )%x( 1 ) * tradial( 1 ),   &
                 xs( 2 ) + ray2D( is + 1 )%x( 1 ) * tradial( 2 ),   &
                           ray2D( is + 1 )%x( 2 ) ]

           s1     = ( x( 1 ) - Botx( 1 ) ) / Bot_deltax   ! proportional distance along segment
           s2     = ( x( 2 ) - Botx( 2 ) ) / Bot_deltay   ! proportional distance along segment

           BotnInt = Bot( IsegBotx,     IsegBoty     )%Noden * ( 1 - s1 ) * ( 1 - s2 ) +  &
                     Bot( IsegBotx + 1, IsegBoty     )%Noden * ( s1     ) * ( 1 - s2 ) +  &
                     Bot( IsegBotx + 1, IsegBoty + 1 )%Noden * ( s1     ) * ( s2     ) +  &
                     Bot( IsegBotx,     IsegBoty + 1 )%Noden * ( 1 - s1 ) * ( s2     )
        ELSE
           BotnInt = Botn   ! normal is constant in a segment
        END IF

        CALL Reflect2D( is, Bdry%Bot%HS, 'BOT', BotnInt, RBot, NBotPTS, xs, tradial )
        ray2D( is + 1 )%NumBotBnc = ray2D( is )%NumBotBnc + 1

        x = [ xs( 1 ) + ray2D( is + 1 )%x( 1 ) * tradial( 1 ),   &
              xs( 2 ) + ray2D( is + 1 )%x( 1 ) * tradial( 2 ),   &
                        ray2D( is + 1 )%x( 2 ) ]

        CALL Distances3D( x, Topx, Botx, Topn, Botn, DistEndTop, DistEndBot )

     END IF

     ! Has the ray left the box, lost its energy, escaped the boundaries, or exceeded storage limit?
     !!!! this should be modified to have a single box
     !!!! no need to test x( 1 ), for instance, against several limits; calculate one limit in advance
     IF ( ABS( x( 1 ) - xs( 1 ) ) > Beam%Box%x .OR. &
          ABS( x( 2 ) - xs( 2 ) ) > Beam%Box%y .OR. &
          ABS( x( 3 ) - xs( 3 ) ) > Beam%Box%z .OR. &
          x( 1 ) < MAX( BotGlobalx( 1            ), TopGlobalx( 1            ) ) .OR. &
          x( 2 ) < MAX( BotGlobaly( 1            ), TopGlobaly( 1            ) ) .OR. &
          x( 1 ) > MIN( BotGlobalx( NBTYPts( 1 ) ), TopGlobalx( NATIPts( 1 ) ) ) .OR. &
          x( 2 ) > MIN( BotGlobaly( NBTYPts( 2 ) ), TopGlobaly( NATIPts( 2 ) ) ) .OR. &
          ray2D( is + 1 )%Amp < 0.005 .OR. &
          ! ray2D( is + 1 )%t( 1 ) < 0  .OR. & ! kills off a backward traveling ray
          iSmallStepCtr > 50 ) THEN
        Beam%Nsteps = is + 1
        EXIT Stepping
     ELSE IF ( is >= MaxN - 3 ) THEN
        CALL ERROUT( PRTFile, 'W', 'TraceRay2D', 'Insufficient storage for ray trajectory' )
        WRITE( PRTFile, * ) 'Angles are  alpha = ', alpha * RadDeg, '    beta = ', beta * RadDeg
        Beam%Nsteps = is
        EXIT Stepping
     END IF

     DistBegTop = DistEndTop
     DistBegBot = DistEndBot

  END DO Stepping

END SUBROUTINE TraceRay2D

! **********************************************************************!

SUBROUTINE Step2D( ray0, ray2, xs, tradial )

  ! Does a single step along the ray
  ! x denotes the ray coordinate, (r,z)
  ! t denotes the scaled tangent to the ray (previously (rho, zeta))
  ! c * t would be the unit tangent

  USE bellhopMod
  USE Bdry3DMod
  USE sspMod
  IMPLICIT NONE
  REAL (KIND=8), INTENT( IN ) :: xs( 3 ), tradial( 2 )   ! coordinate of source and ray bearing angle
  TYPE( ray2DPt )    :: ray0, ray1, ray2
  INTEGER            :: iSegx0, iSegy0, iSegz0
  REAL     (KIND=8 ) :: gradc0( 2 ), gradc1( 2 ), gradc2( 2 ), rho, &
                        c0, cimag0, crr0, crz0, czz0, csq0, cnn0_csq0, &
                        c1, cimag1, crr1, crz1, czz1, csq1, cnn1_csq1, &
                        c2, cimag2, crr2, crz2, czz2, urayt0( 2 ), urayt1( 2 ), &
                        h, halfh, hw0, hw1, ray2n( 2 ), RM, RN, gradcjump( 2 ), cnjump, csjump, w0, w1, &
                        rayx3D( 3 ), rayt3D( 3 ) 

  ! The numerical integrator used here is a version of the polygon (a.k.a. midpoint, leapfrog, or Box method), and similar
  ! to the Heun (second order Runge-Kutta method).
  ! However, it's modified to allow for a dynamic step change, while preserving the second-order accuracy).

  ! *** Phase 1 (an Euler step)

  CALL EvaluateSSP2D( ray0%x, c0, cimag0, gradc0, crr0, crz0, czz0, rho, xs, tradial, freq )

  csq0      = c0 * c0
  cnn0_csq0 = crr0 * ray0%t( 2 )**2 - 2.0 * crz0 * ray0%t( 1 ) * ray0%t( 2 ) + czz0 * ray0%t( 1 )**2

  iSegx0 = iSegx     ! make note of current layer
  iSegy0 = iSegy
  iSegz0 = iSegz

  h = Beam%deltas            ! initially set the step h, to the basic one, deltas

  urayt0 = c0 * ray0%t  ! unit tangent
  rayx3D = [ xs( 1 ) + ray0%x( 1 ) * tradial( 1 ), xs( 2 ) + ray0%x( 1 ) * tradial( 2 ), ray0%x( 2 ) ]
  rayt3D = [           urayt0( 1 ) * tradial( 1 ),           urayt0( 1 ) * tradial( 2 ), urayt0( 2 ) ]

  CALL ReduceStep3D( rayx3D, rayt3D, iSegx0, iSegy0, iSegz0, Beam%deltas, h ) ! reduce h to land on boundary

  halfh = 0.5 * h   ! first step of the modified polygon method is a half step

  ray1%x = ray0%x + halfh * urayt0
  ray1%t = ray0%t - halfh * gradc0 / csq0
  ray1%p = ray0%p - halfh * cnn0_csq0 * ray0%q
  ray1%q = ray0%q + halfh * c0        * ray0%p

  ! *** Phase 2

  CALL EvaluateSSP2D( ray1%x, c1, cimag1, gradc1, crr1, crz1, czz1, rho, xs, tradial, freq )
  csq1      = c1 * c1
  cnn1_csq1 = crr1 * ray1%t( 2 )**2 - 2.0 * crz1 * ray1%t( 1 ) * ray1%t( 2 ) + czz1 * ray1%t( 1 )**2

  ! The Munk test case with a horizontally launched ray caused problems.
  ! The ray vertexes on an interface and can ping-pong around that interface.
  ! Have to be careful in that case about big changes to the stepsize (that invalidate the leap-frog scheme) in phase II.
  ! A modified Heun or Box method could also work.

  urayt1 = c1 * ray1%t   ! unit tangent
  rayt3D = [           urayt1( 1 ) * tradial( 1 ),           urayt1( 1 ) * tradial( 2 ), urayt1( 2 ) ]

  CALL ReduceStep3D( rayx3D, rayt3D, iSegx0, iSegy0, iSegz0, Beam%deltas, h ) ! reduce h to land on boundary

  ! use blend of f' based on proportion of a full step used.
  w1  = h / ( 2.0d0 * halfh )
  w0  = 1.0d0 - w1
  hw0 = h * w0
  hw1 = h * w1

  ray2%x   = ray0%x   + hw0 * c0 * ray0%t        + hw1 * c1 * ray1%t
  ray2%t   = ray0%t   - hw0 * gradc0 / csq0      - hw1 * gradc1 / csq1
  ray2%p   = ray0%p   - hw0 * cnn0_csq0 * ray0%q - hw1 * cnn1_csq1 * ray1%q
  ray2%q   = ray0%q   + hw0 * c0        * ray0%p + hw1 * c1        * ray1%p
  ray2%tau = ray0%tau + hw0 / CMPLX( c0, cimag0, KIND=8 ) + hw1 / CMPLX( c1, cimag1, KIND=8 )

  ray2%Amp       = ray0%Amp
  ray2%Phase     = ray0%Phase
  ray2%NumTopBnc = ray0%NumTopBnc
  ray2%NumBotBnc = ray0%NumBotBnc

  ! If we crossed an interface, apply jump condition

  CALL EvaluateSSP2D( ray2%x, c2, cimag2, gradc2, crr2, crz2, czz2, rho, xs, tradial, freq )
  ray2%c = c2

  !!! this needs modifying like the full 3D version to handle jumps in the x-y direction
  IF ( iSegz /= iSegz0 ) THEN
     gradcjump =  gradc2 - gradc0  ! this is precise only for c-linear layers
     ray2n     = [ -ray2%t( 2 ), ray2%t( 1 ) ]

     cnjump    = DOT_PRODUCT( gradcjump, ray2n  )
     csjump    = DOT_PRODUCT( gradcjump, ray2%t )

     RM        = ray2%t( 1 ) / ray2%t( 2 )
     RN        = RM * ( 2 * cnjump - RM * csjump ) / c2
     RN        = -RN
     ray2%p    = ray2%p + ray2%q * RN

  END IF

END SUBROUTINE Step2D

! **********************************************************************!

SUBROUTINE Reflect2D( is, HS, BotTop, nBdry3d, RefC, Npts, xs, tradial )

  USE bellhopMod
  USE RefCoMod
  USE norms
  IMPLICIT NONE
  INTEGER,              INTENT( IN ) :: Npts
  REAL     (KIND=8),    INTENT( IN ) :: xs( 3 ), tradial( 2 )
  REAL     (KIND=8),    INTENT( IN ) :: nBdry3d( 3 )            ! Normal to the boundary
  CHARACTER (LEN=3),    INTENT( IN ) :: BotTop                  ! Flag indicating bottom or top reflection
  TYPE( HSInfo ),       INTENT( IN ) :: HS                      ! half-space properties
  TYPE(ReflectionCoef), INTENT( IN ) :: RefC( NPts )            ! reflection coefficient
  INTEGER,              INTENT( INOUT ) :: is
  INTEGER           :: is1
  REAL     (KIND=8) :: c, cimag, gradc( 2 ), crr, crz, czz, cn, cs, rho             ! derivatives of sound speed
  REAL     (KIND=8) :: RM, RN, Tg, Th, rayt( 2 ), rayn( 2 ), rayt_tilde( 2 ), rayn_tilde( 2 ), cnjump, csjump  ! for curvature change
  REAL     (KIND=8) :: ck, co, si, cco, ssi, pdelta, rddelta, sddelta   ! for beam shift
  COMPLEX  (KIND=8) :: gamma1, gamma2, gamma1Sq, gamma2Sq, GK, Refl, ch, a, b, d, sb, delta, ddelta
  REAL     (KIND=8) :: kappa                                            ! Boundary curvature
  REAL     (KIND=8) :: tBdry( 2 ), nBdry( 2 )
  REAL     (KIND=8) :: nBdry3dCone( 3 ), theta, phi, Rray   ! for cone reflection
  TYPE(ReflectionCoef) :: RInt

  is  = is + 1
  is1 = is + 1

  ! following should perhaps be pre-calculated in Bdry3DMod
  nBdry( 1 ) = DOT_PRODUCT( nBdry3d( 1 : 2 ), tradial )
  nBdry( 2 ) = nBdry3d( 3 )

!!$  ! special case of a conical seamount
!!$
!!$  IF ( BotTop == 'BOT' ) THEN
!!$     phi   = 15 * DegRad   ! 15 degree angle of seamount
!!$     Rray  = ray2D( is )%x( 1 )   ! cylindrical range of ray from source
!!$     ! get bearing from cone to ray using (x,y) coordinate of the ray
!!$     theta = atan2( xs( 2 ) + Rray * tradial( 2 ), xs( 1 ) + Rray * tradial( 1 ) )
!!$
!!$     nBdry3dCone =  [ -cos( theta ) * sin( phi ), -sin( theta ) * sin( phi ), cos( phi ) ]
!!$
!!$     ! nBdry3dCone = -[ -cos( theta ) * tan( phi ), -sin( theta ) * tan( phi ), 1.0D0 ]
!!$     ! nBdry3dCone = -[ ray2D( is )%x( 1 ), ray2D( is )%x( 2 ), -1.0D0  ]
!!$     ! nBdry3dCone = NBdry3dCone / NORM2b( NBdry3dCone )
!!$
!!$     nBdry( 1 ) = DOT_PRODUCT( nBdry3dCone( 1 : 2 ), tradial )
!!$     nBdry( 2 ) = nBdry3dCone( 3 )
!!$  END IF

  Th = DOT_PRODUCT( ray2D( is )%t, nBdry )  ! component of ray tangent, normal to boundary

  tBdry = ray2D( is )%t - Th * nBdry  ! calculate the component of the tangent along the boundary
  tBdry = tBdry / NORM2b( tBdry )
  ! could also calculate tbdry as +/- of [ nbdy( 2), -nbdry( 1 ) ], but need sign

  Tg = DOT_PRODUCT( ray2D( is )%t, tBdry )  ! component of ray tangent, along boundary

  ray2D( is1 )%NumTopBnc = ray2D( is )%NumTopBnc
  ray2D( is1 )%NumBotBnc = ray2D( is )%NumBotBnc
  ray2D( is1 )%x         = ray2D( is )%x
  ray2D( is1 )%t         = ray2D( is )%t - 2.0 * Th * nBdry  ! changing the ray direction
  ! Calculate the change in curvature
  ! Based on formulas given by Muller, Geoph. J. R.A.S., 79 (1984).

  CALL EvaluateSSP2D( ray2D( is1 )%x, c, cimag, gradc, crr, crz, czz, rho, xs, tradial, freq )

  ! incident unit ray tangent and normal
  rayt = c * ray2D( is )%t           ! unit tangent to ray
  rayn = [ -rayt( 2 ), rayt( 1 ) ]   ! unit normal to ray

  ! reflected unit ray tangent and normal (the reflected tangent, normal system has a different orientation)
  rayt_tilde = c * ray2D( is1 )%t                       ! unit tangent to ray
  rayn_tilde = -[ -rayt_tilde( 2 ), rayt_tilde( 1 ) ]   ! unit normal to ray

  kappa = 0   ! no boundary curvature in BELLHOP3D
  RN = 2 * kappa / c**2 / Th    ! boundary curvature correction

  ! get the jumps (this could be simplified, e.g. jump in rayt is roughly 2 * Th * nbdry
  cnjump = -DOT_PRODUCT( gradc, rayn_tilde - rayn  )
  csjump = -DOT_PRODUCT( gradc, rayt_tilde - rayt )

  IF ( BotTop == 'TOP' ) THEN
     cnjump = -cnjump    ! flip sign for top reflection
     RN = -RN
  END IF

  RM = Tg / Th
  RN = RN + RM * ( 2 * cnjump - RM * csjump ) / c ** 2

  SELECT CASE ( Beam%Type( 2 : 2 ) )
  CASE ( 'D' )
     RN = 2.0 * RN
  CASE ( 'Z' )
     RN = 0.0
  END SELECT

  ray2D( is1 )%c   = c
  ray2D( is1 )%p   = ray2D( is )%p + ray2D( is )%q * RN
  ray2D( is1 )%q   = ray2D( is )%q
  ray2D( is1 )%tau = ray2D( is )%tau

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
  CASE ( 'A', 'G' )                 ! half-space
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

           print *, 'beam displacing ...'
           ray2D( is1 )%x( 1 ) = ray2D( is1 )%x( 1 ) + real( delta )   ! displacement
           ray2D( is1 )%tau    = ray2D( is1 )%tau + pdelta             ! phase change
           ray2D( is1 )%q      = ray2D( is1 )%q + sddelta * rddelta * si * c * ray2D( is )%p   ! beam-width change
        endif

     ENDIF
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
  INTEGER                     :: is, N2, iSkip

  ! compression

  N2    = 1
  iSkip = MAX( Nsteps1 / 5000, 1 )   ! max #pts written is about 5000

  Stepping: DO is = 2, Nsteps1
     ! ensure that we always write ray pts. near bdry reflections (works only for flat bdry)
     IF ( MIN( Bdry%Bot%HS%Depth - ray2D( is )%x( 2 ),  ray2D( is )%x( 2 ) - Bdry%Top%HS%Depth ) < 0.2 .OR. &
          MOD( is, iSkip ) == 0 .OR. is == Nsteps1 ) THEN
        N2 = N2 + 1
        ray2D( N2 )%x = ray2D( is )%x
     END IF
  END DO Stepping

  ! write to ray file

  WRITE( RAYFile, * ) alpha0
  WRITE( RAYFile, * ) N2, ray2D( Nsteps1 )%NumTopBnc, ray2D( Nsteps1 )%NumBotBnc
  DO is = 1, N2
     WRITE( RAYFile, * ) SNGL( ray2D( is )%x )
  END DO

END SUBROUTINE WriteRay

!**********************************************************************!

SUBROUTINE EvaluateSSP2D( x2D, c, cimag, gradc, crr, crz, czz, rho, xs, tradial, freq )

  ! Converts cartesian gradients to polar
  ! should this be inside SSPMod?

  USE sspMod
  IMPLICIT NONE
  REAL (KIND=8), INTENT( IN  ) :: x2D( 2 ), xs( 3 ), tradial( 2 ), freq
  REAL (KIND=8), INTENT( OUT ) :: c, cimag, gradc( 2 ), czz, crz, crr, rho
  REAL (KIND=8)                :: x( 3 ), gradc3D(3 ), cxx, cyy, cxy, cxz, cyz

  ! convert polar coordinate to cartesian
  x = [ xs( 1 ) + x2D( 1 ) * tradial( 1 ), xs( 2 ) + x2D( 1 ) * tradial( 2 ), x2D( 2 ) ]

  CALL EvaluateSSP3D( x, c, cimag, gradc3D, cxx, cyy, czz, cxy, cxz, cyz, rho, freq, 'TAB' )

  gradc( 1 )  = DOT_PRODUCT( tradial, gradc3D( 1 : 2 ) )  ! r derivative
  gradc( 2 )  = gradc3D( 3 )                              ! z derivative

  crz = tradial( 1 ) * cxz + tradial( 2 ) * cyz
  crr = cxx * ( tradial( 1 ) )**2 + 2.0 * cxy * tradial( 1 ) * tradial( 2 ) + cyy * ( tradial( 2 ) )**2

  RETURN
END SUBROUTINE EvaluateSSP2D

! **********************************************************************!

SUBROUTINE WriteRay3D( alpha0, beta0, Nsteps1, xs )

  ! Compress the ray data keeping every iSkip point, points near surface or bottom, and last point.
  ! Write to RAYFile.
  ! During an eigenray calculation, subsets of the full ray may be passed
  ! These have lengths Nsteps1 vs. Nsteps for the entire ray

  USE bellhopMod
  IMPLICIT NONE
  INTEGER, PARAMETER :: RAYFile = 21
  INTEGER,       INTENT( IN ) :: Nsteps1
  REAL (KIND=8), INTENT( IN ) :: alpha0, beta0   ! take-off angle of this ray
  REAL (KIND=8), INTENT( IN ) :: xs( 3 )              ! source location
  INTEGER :: is, N2, iSkip

  ! if Nx2D run, copy r-z rays to x-y-z rays

  IF ( Beam%RunType( 6 : 6 ) == '2' ) THEN
     ray3D%x( 1 )    = xs( 1 ) + ray2D%x( 1 ) * COS( beta0 )
     ray3D%x( 2 )    = xs( 2 ) + ray2D%x( 1 ) * SIN( beta0 )
     ray3D%x( 3 )    = ray2D%x( 2 )
     ray3D%NumTopBnc = ray2D%NumTopBnc
     ray3D%NumBotBnc = ray2D%NumBotBnc
  END IF

  ! compression

  N2    = 1
  iSkip = MAX( Nsteps1 / 5000, 1 )   ! max #points written is about 5000
  iSkip = 1
  
  Stepping: DO is = 2, Nsteps1
     ! ensure that we always write ray points near boundary reflections
     IF ( MIN( Bdry%Bot%HS%Depth - ray3D( is )%x( 3 ),  ray3D( is )%x( 3 ) - Bdry%Top%HS%Depth ) < 0.2 .OR. &
          MOD( is, iSkip ) == 0 .OR. is == Nsteps1 ) THEN
        N2 = N2 + 1
        ray3D( N2 )%x = ray3D( is )%x
     END IF
  END DO Stepping

  ! write to ray file

  WRITE( RAYFile, * ) alpha0
  WRITE( RAYFile, * ) N2, ray3D( Nsteps1 )%NumTopBnc, ray3D( Nsteps1 )%NumBotBnc
  DO is = 1, N2
     WRITE( RAYFile, * ) SNGL( ray3D( is )%x )
  END DO

END SUBROUTINE WriteRay3D

!**********************************************************************!

SUBROUTINE TraceRay3D( xs, alpha, beta, epsilon, Amp0 )

  ! Traces the beam corresponding to a particular take off angle

  USE bellhopMod
  USE RefCoMod
  USE bdry3DMod
  USE norms
  USE sspMod
  IMPLICIT NONE
  REAL     ( KIND=8 ), INTENT( IN ) :: xs( 3 ), Amp0  ! source coordinate, initial amplitude
  REAL     ( KIND=8 ), INTENT( IN ) :: alpha, beta    ! take-off angles of the ray
  COMPLEX  ( KIND=8 ), INTENT( IN ) :: epsilon( 2 )   ! beam initial conditions
  INTEGER             :: is, is1
  REAL     ( KIND=8 ) :: DistBegTop, DistEndTop, DistBegBot, DistEndBot, &
                         c, cimag, gradc( 3 ), cxx, cyy, czz, cxy, cxz, cyz, rho   ! soundspeed derivatives
  REAL     (KIND=8) :: TopnInt( 3 ), BotnInt( 3 )
  REAL     (KIND=8) :: s1, s2
  ! *** Initial conditions ***

  iSmallStepCtr = 0
  ray3D( 1 )%x    = xs
  CALL EvaluateSSP3D(  ray3D( 1 )%x, c, cimag, gradc, cxx, cyy, czz, cxy, cxz, cyz, rho, freq, 'TAB' )

  ray3D( 1 )%t    = [ COS( alpha ) * COS( beta ) / c, COS( alpha ) * SIN( beta ) / c, SIN( alpha ) / c ]
  !ray3D( 1 )%f    = epsilon( 2 )
  !ray3D( 1 )%g    = epsilon( 1 )
  !ray3D( 1 )%h    = 0.0
  !ray3D( 1 )%DetP = 1.0
  !ray3D( 1 )%DetQ = epsilon( 1 ) * epsilon( 2 )
  ray3D( 1 )%c    = c
  ray3D( 1 )%phi  = 0.0

  ray3D( 1 )%tau       = 0.0
  ray3D( 1 )%Amp       = Amp0
  ray3D( 1 )%Phase     = 0.0
  ray3D( 1 )%NumTopBnc = 0
  ray3D( 1 )%NumBotBnc = 0

  !ray3D( 1 )%p_tilde = [ 1.0, 0.0 ]   use these for complex Cerveny beams
  !ray3D( 1 )%q_tilde = [ epsilon( 1 ), CMPLX( 0.0, KIND=8 ) ]
  !ray3D( 1 )%p_hat   = [ 0.0, 1.0 ]
  !ray3D( 1 )%q_hat   = [ CMPLX( 0.0, KIND=8 ), epsilon( 2 ) ]

  ray3D( 1 )%p_tilde = [ 1.0, 0.0 ]
  ray3D( 1 )%q_tilde = [ 0.0, 0.0 ]
  ray3D( 1 )%p_hat   = [ 0.0, 1.0 ]
  ray3D( 1 )%q_hat   = [ 0.0, 0.0 ]

  ! dummy BotSeg info to force GetBotSeg to search for the active segment on first call
  xTopSeg = [ +big, -big ]
  yTopSeg = [ +big, -big ]
  xBotSeg = [ +big, -big ]
  yBotSeg = [ +big, -big ]
    
  CALL GetTopSeg3D( xs )   ! identify the top    segment above the source
  CALL GetBotSeg3D( xs )   ! identify the bottom segment below the source
  
  ! Trace the beam (note that Reflect alters the step index is)
  is = 0

  CALL Distances3D( ray3D( 1 )%x, Topx,  Botx, Topn, Botn, DistBegTop, DistBegBot )

  IF ( DistBegTop <= 0 .OR. DistBegBot <= 0 ) THEN
     Beam%Nsteps = 1
     RETURN       ! source must be within the medium
  END IF

  Stepping: DO istep = 1, MaxN - 1
     is  = is + 1
     is1 = is + 1

     CALL Step3D( ray3D( is ), ray3D( is1 ) )
     CALL GetTopSeg3D( ray3D( is1 )%x )   ! identify the top    segment above the source
     CALL GetBotSeg3D( ray3D( is1 )%x )   ! identify the bottom segment below the source

     IF ( IsegTopx == 0 .OR. IsegTopy == 0 .OR. IsegBotx == 0 .OR. IsegBoty == 0 ) THEN ! we escaped the box
        Beam%Nsteps = is
        EXIT Stepping
     END IF

     ! Reflections?
     ! Tests that ray at step is is inside, and ray at step is+1 is outside
     ! to detect only a crossing from inside to outside
     ! DistBeg is the distance at step is, which is saved
     ! DistEnd is the distance at step is+1, which needs to be calculated
  
     CALL Distances3D( ray3D( is1 )%x, Topx, Botx, Topn, Botn, DistEndTop, DistEndBot )

     IF      ( DistBegTop > 0.0d0 .AND. DistEndTop <= 0.0d0 ) THEN  ! test top reflection
        IF ( atiType == 'C' ) THEN
           s1 = ( ray3D( is1 )%x( 1 ) - Topx( 1 ) ) / Top_deltax   ! proportional distance along segment
           s2 = ( ray3D( is1 )%x( 2 ) - Topx( 2 ) ) / Top_deltay   ! proportional distance along segment

           TopnInt = Top( IsegTopx,     IsegTopy     )%Noden * ( 1 - s1 ) * ( 1 - s2 ) +  &
                     Top( IsegTopx + 1, IsegTopy     )%Noden * ( s1     ) * ( 1 - s2 ) +  &
                     Top( IsegTopx + 1, IsegTopy + 1 )%Noden * ( s1     ) * ( s2     ) +  &
                     Top( IsegTopx,     IsegTopy + 1 )%Noden * ( 1 - s1 ) * ( s2     )
        ELSE
           TopnInt = Topn   ! normal is constant in a segment
        END IF

        CALL Reflect3D( is, Bdry%Top%HS, 'TOP', TopnInt, RTop, NTopPTS )
        ray3D( is + 1 )%NumTopBnc = ray3D( is )%NumTopBnc + 1
        CALL Distances3D( ray3D( is + 1 )%x, Topx,  Botx, Topn, Botn, DistEndTop, DistEndBot )

     ELSE IF ( DistBegBot > 0.0d0 .AND. DistEndBot <= 0.0d0 ) THEN  ! test bottom reflection
                !write( *, * ) 'bottom reflection', is, DistBegBot, DistEndBot

        IF ( btyType == 'C' ) THEN
           s1 = ( ray3D( is1 )%x( 1 ) - Botx( 1 ) ) / Bot_deltax   ! proportional distance along segment
           s2 = ( ray3D( is1 )%x( 2 ) - Botx( 2 ) ) / Bot_deltay   ! proportional distance along segment

           BotnInt = Bot( IsegBotx,     IsegBoty     )%Noden * ( 1 - s1 ) * ( 1 - s2 ) +  &
                     Bot( IsegBotx + 1, IsegBoty     )%Noden * ( s1     ) * ( 1 - s2 ) +  &
                     Bot( IsegBotx + 1, IsegBoty + 1 )%Noden * ( s1     ) * ( s2     ) +  &
                     Bot( IsegBotx,     IsegBoty + 1 )%Noden * ( 1 - s1 ) * ( s2     )
        ELSE
           BotnInt = Botn   ! normal is constant in a segment
        END IF

        CALL Reflect3D( is, Bdry%Bot%HS, 'BOT', BotnInt, RBot, NBotPTS )
        !CALL Reflect3D( is, Bdry%Bot%HS, 'BOT', ( cos( theta ) * sin( phi ), , RBot, NBotPTS )
        ray3D( is + 1 )%NumBotBnc = ray3D( is )%NumBotBnc + 1
        CALL Distances3D( ray3D( is + 1 )%x, Topx, Botx, Topn, Botn, DistEndTop, DistEndBot )

     END IF

     ! Has the ray left the box, lost its energy, escaped the boundaries, or exceeded storage limit?
     IF ( ABS( ray3D( is + 1 )%x( 1 ) - xs( 1 ) ) > Beam%Box%x .OR. &
          ABS( ray3D( is + 1 )%x( 2 ) - xs( 2 ) ) > Beam%Box%y .OR. &
          ABS( ray3D( is + 1 )%x( 3 ) - xs( 3 ) ) > Beam%Box%z .OR. &
          ray3D( is + 1 )%x( 1 ) < MAX( BotGlobalx( 1            ), TopGlobalx( 1            ) ) .OR. &
          ray3D( is + 1 )%x( 2 ) < MAX( BotGlobaly( 1            ), TopGlobaly( 1            ) ) .OR. &
          ray3D( is + 1 )%x( 1 ) > MIN( BotGlobalx( NBTYPts( 1 ) ), TopGlobalx( NATIPts( 1 ) ) ) .OR. &
          ray3D( is + 1 )%x( 2 ) > MIN( BotGlobaly( NBTYPts( 2 ) ), TopGlobaly( NATIPts( 2 ) ) ) .OR. &
          ray3D( is + 1 )%Amp < 0.005 .OR. &
          iSmallStepCtr > 50 ) THEN
        Beam%Nsteps = is + 1
        EXIT Stepping
     ELSE IF ( is >= MaxN - 3 ) THEN
        Beam%Nsteps = is
        CALL ERROUT( PRTFile, 'W', 'TraceRay3D', 'Insufficient storage for ray trajectory' )
        WRITE( PRTFile, * ) 'Angles are  alpha = ', alpha * RadDeg, '    beta = ', beta * RadDeg
        EXIT Stepping
     END IF

     DistBegTop = DistEndTop
     DistBegBot = DistEndBot

  END DO Stepping
END SUBROUTINE TraceRay3D

! **********************************************************************!

SUBROUTINE Distances3D( rayx, Topx, Botx, Topn, Botn, DistTop, DistBot )

  ! Computes distances from ray to boundaries
  ! Formula differs from JKPS because code uses outward pointing normals

  IMPLICIT NONE
  REAL (KIND=8), INTENT( IN  ) :: rayx( 3 )             ! ray coordinate
  REAL (KIND=8), INTENT( IN  ) :: Topx( 3 ), Botx( 3 )  ! top, bottom boundary coordinate for node
  REAL (KIND=8), INTENT( IN  ) :: Topn( 3 ), Botn( 3 )  ! top, bottom boundary normal
  REAL (KIND=8), INTENT( OUT ) :: DistTop, DistBot      ! distance from the ray to top, bottom boundaries 
  REAL (KIND=8)                :: dTop( 3 ), dBot( 3 )

  dTop    = rayx - Topx  ! vector pointing from top    to ray
  dBot    = rayx - Botx  ! vector pointing from bottom to ray
  DistTop = -DOT_PRODUCT( Topn, dTop )
  DistBot = -DOT_PRODUCT( Botn, dBot )

END SUBROUTINE Distances3D

!************************************************************************

SUBROUTINE Step3D( ray0, ray2 )

  ! Does a single step along the ray
  ! x denotes the ray coordinate, ( x, y, z )
  ! t denotes the scaled tangent to the ray (previously (xi, eta, zeta) )
  ! c * t would be the unit tangent

  USE bellhopMod   ! needed only to get Type( ray3DPt )
  USE Bdry3DMod
  USE sspMod
  USE cross_products
  USE norms
  IMPLICIT NONE

  ! rays
  TYPE( ray3DPt ) :: ray0, ray1, ray2
  INTEGER         :: iSegx0, iSegy0, iSegz0
  REAL  (KIND=8 ) :: gradc0( 3 ), gradc1( 3 ), gradc2( 3 ), &
                     c0, cimag0, csq0, cxx0, cyy0, czz0, cxy0, cxz0, cyz0, cnn0, cmn0, cmm0, &
                     c1, cimag1, csq1, cxx1, cyy1, czz1, cxy1, cxz1, cyz1, cnn1, cmn1, cmm1, &
                     c2, cimag2,       cxx2, cyy2, czz2, cxy2, cxz2, cyz2, c_mat0( 2, 2 ), c_mat1( 2, 2 ), &
                     urayt0( 3 ), urayt1( 3 ), h, halfh, hw0, hw1, w0, w1, rho
  REAL  (KIND=8 ) :: gradcjump( 3 ), csjump, cn1jump, cn2jump, tBdry( 3 ), nBdry( 3 ), RM, R1, R2, Tg, Th, &
                     rayt( 3 ), rayn1( 3 ), rayn2( 3 )
  REAL  (KIND=8 ) :: e1( 3 ), e2( 3 )                              ! ray normals for ray-centered coordinates
  REAL  (KIND=8 ) :: p_tilde_in(  2 ), p_hat_in(  2 ), q_tilde_in(  2 ), q_hat_in(  2 ), &
                     p_tilde_out( 2 ), p_hat_out( 2 )

  ! The numerical integrator used here is a version of the polygon (a.k.a. midpoint, leapfrog, or Box method), and similar
  ! to the Heun (second order Runge-Kutta method).
  ! However, it's modified to allow for a dynamic step change, while preserving the second-order accuracy).

  ! *** Phase 1 (an Euler step)

  CALL EvaluateSSP3D( ray0%x, c0, cimag0, gradc0, cxx0, cyy0, czz0, cxy0, cxz0, cyz0, rho, freq, 'TAB' )
  CALL RayNormal( ray0%t, ray0%phi, c0, e1, e2 ) ! Compute ray normals e1 and e2
  CALL Get_c_partials( cxx0, cxy0, cxz0, cyy0, cyz0, czz0, e1, e2, cnn0, cmn0, cmm0 ) ! Compute second partials of c along ray normals

  csq0   = c0 * c0

  iSegx0 = iSegx     ! make note of current layer
  iSegy0 = iSegy
  iSegz0 = iSegz

  h = Beam%deltas            ! initially set the step h, to the basic one, deltas
  urayt0 = c0 * ray0%t  ! unit tangent

  CALL ReduceStep3D( ray0%x, urayt0, iSegx0, iSegy0, iSegz0, Beam%deltas, h ) ! reduce h to land on boundary

  halfh = 0.5 * h       ! first step of the modified polygon method is a half step

  ray1%x = ray0%x    + halfh * urayt0
  ray1%t = ray0%t    - halfh * gradc0 / csq0

  ! ray1%f    = ray0%f    + halfh * ( c0 * ray0%DetP - cnn0 / csq0 * ray0%DetQ )
  ! ray1%g    = ray0%g    + halfh * ( c0 * ray0%DetP - cmm0 / csq0 * ray0%DetQ )
  ! ray1%h    = ray0%h    + halfh * (                - cmn0 / csq0 * ray0%DetQ )
  ! ray1%DetP = ray0%DetP + halfh / csq0 * ( -cmm0 * ray0%f - cnn0 * ray0%g + 2.0 * cmn0 * ray0%h )
  ! ray1%DetQ = ray0%DetQ + halfh * c0 * ( ray0%f + ray0%g )

  ray1%phi  = ray0%phi  + halfh * ( 1.0 / c0 ) * ray0%t( 3 ) * &
              ( ray0%t( 2 ) * gradc0( 1 ) - ray0%t( 1 ) * gradc0( 2 ) ) / &
              ( ray0%t( 1 ) ** 2 + ray0%t( 2 ) ** 2 )

  c_mat0( 1, : ) = -[ cnn0, cmn0 ] / csq0
  c_mat0( 2, : ) = -[ cmn0, cmm0 ] / csq0

  ray1%p_tilde = ray0%p_tilde + halfh * MATMUL( c_mat0, ray0%q_tilde )
  ray1%q_tilde = ray0%q_tilde + halfh *         c0    * ray0%p_tilde 

  ray1%p_hat   = ray0%p_hat   + halfh * MATMUL( c_mat0, ray0%q_hat )
  ray1%q_hat   = ray0%q_hat   + halfh *         c0    * ray0%p_hat 

  ! *** Phase 2

  CALL EvaluateSSP3D( ray1%x, c1, cimag1, gradc1, cxx1, cyy1, czz1, cxy1, cxz1, cyz1, rho, freq, 'TAB' )
  rayt = c1 * ray1%t           ! unit tangent to ray
  CALL RayNormal_unit( rayt, ray2%phi, e1, e2 )
  CALL Get_c_partials( cxx1, cxy1, cxz1, cyy1, cyz1, czz1, e1, e2, cnn1, cmn1, cmm1 ) ! Compute second partials of c along ray normals

  csq1   = c1 * c1
  urayt1 = c1 * ray1%t   ! unit tangent

  CALL ReduceStep3D( ray0%x, urayt1, iSegx0, iSegy0, iSegz0, Beam%deltas, h ) ! reduce h to land on boundary

  ! use blend of f' based on proportion of a full step used.
  w1  = h / ( 2.0d0 * halfh )
  w0  = 1.0d0 - w1
  hw0 = h * w0
  hw1 = h * w1

  ray2%x = ray0%x    + hw0 * urayt0        + hw1 * urayt1
  ray2%t = ray0%t    - hw0 * gradc0 / csq0 - hw1 * gradc1 / csq1

  ! ERROR: need to do hw0 and hw1 blend here as well !!!
  ! ray2%f    = ray0%f    + h * ( c1 * ray1%DetP - cnn1 / csq1 * ray1%DetQ )
  ! ray2%g    = ray0%g    + h * ( c1 * ray1%DetP - cmm1 / csq1 * ray1%DetQ )
  ! ray2%h    = ray0%h    + h * (                - cmn1 / csq1 * ray1%DetQ )
  ! ray2%DetP = ray0%DetP + h / csq1 * ( -cmm1 * ray1%f - cnn1 * ray1%g + 2.0 * cmn1 * ray1%h )
  ! ray2%DetQ = ray0%DetQ + h * c1 * ( ray1%f + ray1%g )

  ray2%phi  = ray0%phi  + h * ( 1.0 / c1 ) * ray1%t( 3 ) * &
              ( ray1%t( 2 ) * gradc1( 1 ) - ray1%t( 1 ) * gradc1( 2 ) ) / &
              ( ray1%t( 1 ) ** 2 + ray1%t( 2 ) ** 2 )
  ray2%tau  = ray0%tau + hw0 / CMPLX( c0, cimag0, KIND=8 ) + hw1 / CMPLX( c1, cimag1, KIND=8 )

  ray2%Amp       = ray0%Amp
  ray2%Phase     = ray0%Phase
  ray2%NumTopBnc = ray0%NumTopBnc
  ray2%NumBotBnc = ray0%NumBotBnc

  c_mat1( 1, : ) = -[ cnn1, cmn1 ] / csq1
  c_mat1( 2, : ) = -[ cmn1, cmm1 ] / csq1

  ray2%p_tilde = ray0%p_tilde + hw0 * MATMUL( c_mat0, ray0%q_tilde ) + hw1 * MATMUL( c_mat1, ray1%q_tilde )
  ray2%q_tilde = ray0%q_tilde + hw0 *         c0    * ray0%p_tilde   + hw1 *         c1    * ray1%p_tilde

  ray2%p_hat   = ray0%p_hat   + hw0 * MATMUL( c_mat0, ray0%q_hat   ) + hw1 * MATMUL( c_mat1, ray1%q_hat )
  ray2%q_hat   = ray0%q_hat   + hw0 *         c0    * ray0%p_hat     + hw1 *         c1    * ray1%p_hat 

  ! *** If we crossed an interface, apply jump condition ***

  CALL EvaluateSSP3D( ray2%x, c2, cimag2, gradc2, cxx2, cyy2, czz2, cxy2, cxz2, cyz2, rho, freq, 'TAB' )

  ray2%c = c2

  IF ( iSegx /= iSegx0 .OR. &
       iSegy /= iSegy0 .OR. &
       iSegz /= iSegz0 ) THEN

     gradcjump =  gradc2 - gradc0  ! this is precise only for c-linear layers

     !!! what if we cross isegx, isegy, or isegz at the same time?
     IF ( iSegz /= iSegz0 ) THEN
        nBdry = [ 0D0, 0D0, -SIGN( 1.0D0, ray2%t( 3 ) ) ]   ! inward normal to layer
     ELSE IF ( iSegx /= iSegx0 ) THEN
        nBdry = [ -SIGN( 1.0D0, ray2%t( 1 ) ), 0D0, 0D0 ]   ! inward normal to x-sgement
     ELSE
        nBdry = [ 0D0, -SIGN( 1.0D0, ray2%t( 2 ) ), 0D0 ]   ! inward normal to y-sgement
     END IF

     Th    = DOT_PRODUCT( ray2%t, nBdry )   ! component of ray tangent, normal to boundary
     tBdry = ray2%t - Th * nBdry            ! tangent, along the boundary, in the reflection plane
     tBdry = tBdry / NORM2b( tBdry )
     Tg    = DOT_PRODUCT( ray2%t, tBdry )   ! component of ray tangent, along the boundary

     rayt = c2 * ray2%t                     ! unit tangent to ray

     rayn2 = cross_product( rayt, nBdry )   ! ray tangent x boundary normal gives refl. plane normal
     rayn2 = rayn2 / NORM2b( rayn2 )        !!! why does this need normalizing?
     rayn1 = -cross_product( rayt, rayn2 )  ! ray tangent x refl. plane normal is first ray normal

     ! normal and tangential derivatives of the sound speed
     cn1jump = DOT_PRODUCT( gradcjump, rayn1 )
     cn2jump = DOT_PRODUCT( gradcjump, rayn2 )
     csjump  = DOT_PRODUCT( gradcjump, rayt  )

     RM = Tg / Th   ! this is tan( alpha ) where alpha is the angle of incidence
     R1 = RM * ( 2 * cn1jump - RM * csjump ) / c2 ** 2
     R2 = RM * cn2jump / c2 ** 2

     ! *** curvature correction ***

     CALL RayNormal_unit( rayt, ray2%phi, e1, e2 )   ! Compute ray normals e1 and e2

     !!! is the compiler smart enough to pre-calculate all the following dot products?
     ! project p-q values in e1, e2 system, onto rayn1, rayn2 system
     p_tilde_in = DOT_PRODUCT( rayn1, e1 ) * ray2%p_tilde + DOT_PRODUCT( rayn1, e2 ) * ray2%p_hat
     p_hat_in   = DOT_PRODUCT( rayn2, e1 ) * ray2%p_tilde + DOT_PRODUCT( rayn2, e2 ) * ray2%p_hat
     q_tilde_in = DOT_PRODUCT( rayn1, e1 ) * ray2%q_tilde + DOT_PRODUCT( rayn1, e2 ) * ray2%q_hat
     q_hat_in   = DOT_PRODUCT( rayn2, e1 ) * ray2%q_tilde + DOT_PRODUCT( rayn2, e2 ) * ray2%q_hat

     ! here's the actual curvature change
     p_tilde_out = p_tilde_in - q_tilde_in * R1 - q_hat_in * R2  
     p_hat_out   = p_hat_in   - q_tilde_in * R2

     ! rotate back to e1, e2 system
     ray2%p_tilde = DOT_PRODUCT( rayn1, e1 ) * p_tilde_out + DOT_PRODUCT( rayn2, e1 ) * p_hat_out
     ray2%p_hat   = DOT_PRODUCT( rayn1, e2 ) * p_tilde_out + DOT_PRODUCT( rayn2, e2 ) * p_hat_out

     ray2%q_tilde = ray2%q_tilde
     ray2%q_hat   = ray2%q_hat

  END IF

END SUBROUTINE Step3D

! **********************************************************************!

SUBROUTINE ReduceStep3D( x0, urayt, iSegx0, iSegy0, iSegz0, deltas, h )

  USE Bdry3DMod
  USE sspMod
  IMPLICIT NONE
  INTEGER,       INTENT( IN  ) :: iSegx0, iSegy0, iSegz0
  REAL (KIND=8), INTENT( IN  ) :: x0( 3 ), urayt( 3 )  ! ray coordinate and tangent
  REAL (KIND=8), INTENT( IN  ) :: deltas               ! default step size
  REAL (KIND=8), INTENT( OUT ) :: h                    ! reduced step size
  REAL (KIND=8) :: d( 3 ), d0( 3 ), tri_n( 3 )
  REAL (KIND=8) :: x( 3 )                              ! ray coordinate after full trial step
  REAL (KIND=8) :: h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11      ! step sizes

  ! Detect interface or boundary crossing and reduce step, if necessary, to land on that crossing.
  ! Keep in mind possibility that user put source right on an interface
  ! and that multiple events can occur (crossing interface, top, and bottom in a single step).

  x = x0 + h * urayt ! make a trial step
  ! write( *, * ) 'x0, x', x0, x
  ! interface crossing in depth
  ! step reduction is not done for the top or bottom layer
  ! instead the SSP is extrapolated
  ! This prevents problems when the boundaries are outside the domain of the SSP
  h1 = huge( h1 )
  IF ( ABS( urayt( 3 ) ) > EPSILON( h1 ) ) THEN
     IF      ( SSP%z( iSegz0     ) > x(  3 ) .AND. iSegz0     > 1  ) THEN
        h1 = ( SSP%z( iSegz0     ) - x0( 3 ) ) / urayt( 3 )
        ! write( *, * ) 'layer crossing', iSegz0, h1
     ELSE IF ( SSP%z( iSegz0 + 1 ) < x(  3 ) .AND. iSegz0 + 1 < SSP%Nz ) THEN
        h1 = ( SSP%z( iSegz0 + 1 ) - x0( 3 ) ) / urayt( 3 )
        ! write( *, * ) 'layer crossing', iSegz0, h1
     END IF
  END IF

  ! top crossing
  h2 = huge( h2 )
  d  = x - Topx              ! vector from top to ray
  IF ( DOT_PRODUCT( Topn, d )  > EPSILON( h1 ) ) THEN
     d0 = x0 - Topx   ! vector from top    node to ray origin
     h2 = -DOT_PRODUCT( d0, Topn ) / DOT_PRODUCT( urayt, Topn )
     ! write( *, * ) 'top crossing'
  END IF

  ! bottom crossing
  h3 = huge( h3 )
  d  = x - Botx              ! vector from bottom to ray
  IF ( DOT_PRODUCT( Botn, d ) > EPSILON( h1 ) ) THEN
     d0 = x0 - Botx   ! vector from bottom node to ray origin
     h3 = -DOT_PRODUCT( d0, Botn ) / DOT_PRODUCT( urayt, Botn )
     ! write( *, * ) 'bottom crossing'
  END IF

  ! top segment crossing in x
  h4 = huge( h4 )
  IF ( ABS( urayt( 1 ) ) > EPSILON( h1 ) ) THEN
     IF       ( x(  1 ) < xTopSeg( 1 ) ) THEN
        h4 = -( x0( 1 ) - xTopSeg( 1 ) ) / urayt( 1 )
        ! write( *, * ) 'top segment crossing in x'

     ELSE IF  ( x(  1 ) > xTopSeg( 2 ) ) THEN
        h4 = -( x0( 1 ) - xTopSeg( 2 ) ) / urayt( 1 )
        ! write( *, * ) 'top segment crossing in x'
     END IF
  END IF

  ! bottom segment crossing in x
  h5 = huge( h5 )
  IF ( ABS( urayt( 1 ) ) > 1000.0 * EPSILON( h1 ) ) THEN
     IF      ( ( x(  1 ) - xBotSeg( 1 ) ) < 0.0 ) THEN
        h5 =  -( x0( 1 ) - xBotSeg( 1 ) ) / urayt( 1 )
        ! write( *, * ) 'bot segment crossing in x left', h5

     ELSE IF ( ( x(  1 ) - xBotSeg( 2 ) ) > 0.0 ) THEN
        h5 =  -( x0( 1 ) - xBotSeg( 2 ) ) / urayt( 1 )
        ! write( *, * ) 'bot segment crossing in x right', xBotSeg( 2 ), urayt( 1 ), x0( 1 ), x( 1 ), h5

     END IF
  END IF

  ! top segment crossing in y
  h6 = huge( h6 )
  IF ( ABS( urayt( 2 ) ) > 1000.0 * EPSILON( h1 ) ) THEN
     IF       ( x(  2 ) < yTopSeg( 1 ) ) THEN
        h6 = -( x0( 2 ) - yTopSeg( 1 ) ) / urayt( 2 )
        ! write( *, * ) 'top segment crossing in y'
     ELSE IF  ( x(  2 ) > yTopSeg( 2 ) ) THEN
        h6 = -( x0( 2 ) - yTopSeg( 2 ) ) / urayt( 2 )
        ! write( *, * ) 'top segment crossing in y'

     END IF
  END IF

  ! bottom segment crossing in y
  h7 = huge( h7 )
  IF ( ABS( urayt( 2 ) ) > 1000.0 * EPSILON( h1 ) ) THEN
     IF      ( ( x(  2 ) - yBotSeg( 1 ) ) < 0.0 ) THEN
        h7 =  -( x0( 2 ) - yBotSeg( 1 ) ) / urayt( 2 )
        ! write( *, * ) 'bot segment crossing in y down', h7

     ELSE IF ( ( x(  2 ) - yBotSeg( 2 ) ) > 0.0 ) THEN
        h7 =  -( x0( 2 ) - yBotSeg( 2 ) ) / urayt( 2 )
        ! write( *, * ) 'bot segment crossing in y up', h7

     END IF
  END IF

  ! triangle crossing within a top segment
  h8    = huge( h8 )
  d     = x  - Topx   ! vector from bottom node to ray end
  d0    = x0 - Topx   ! vector from bottom node to ray origin
  tri_n = [ -Top_deltay, Top_deltax, 0.0d0 ]

  IF ( ( DOT_PRODUCT( tri_n, d0 ) > 0.0d0 .AND. DOT_PRODUCT( tri_n, d ) <= 0.0d0 ) .OR. &
       ( DOT_PRODUCT( tri_n, d0 ) < 0.0d0 .AND. DOT_PRODUCT( tri_n, d ) >= 0.0d0 )  ) THEN
     h8 = -DOT_PRODUCT( d0, tri_n ) / DOT_PRODUCT( urayt, tri_n )
     ! write( *, * ) 'diagonal crossing'

  END IF

  ! triangle crossing within a bottom segment
  h9    = huge( h9 )
  d     = x  - Botx   ! vector from bottom node to ray end
  d0    = x0 - Botx   ! vector from bottom node to ray origin
  tri_n = [ -Bot_deltay, Bot_deltax, 0.0d0 ]

  IF ( ( DOT_PRODUCT( tri_n, d0 ) > 0.0d0 .AND. DOT_PRODUCT( tri_n, d ) <= 0.0d0 ) .OR. &
       ( DOT_PRODUCT( tri_n, d0 ) < 0.0d0 .AND. DOT_PRODUCT( tri_n, d ) >= 0.0d0 )  ) THEN
     h9 = -DOT_PRODUCT( d0, tri_n ) / DOT_PRODUCT( urayt, tri_n )
     ! write( *, * ) 'diagonal crossing'

  END IF

  ! ocean segment crossing in x
  h10 = huge( h10 )
  h11 = huge( h11 )

  IF ( SSP%Type == 'H' ) THEN
     IF ( ABS( urayt( 1 ) ) > EPSILON( h1 ) ) THEN
        IF        ( x(  1 ) < SSP%Seg%x( iSegx0     ) ) THEN
           h10 = -( x0( 1 ) - SSP%Seg%x( iSegx0     ) ) / urayt( 1 )
           ! write( *, * ) 'ocean segment crossing in x'

        ELSE IF   ( x(  1 ) > SSP%Seg%x( iSegx0 + 1 ) ) THEN
           h10 = -( x0( 1 ) - SSP%Seg%x( iSegx0 + 1 ) ) / urayt( 1 )
           ! write( *, * ) 'ocean segment crossing in x'

        END IF
     END IF

     ! ocean segment crossing in y
     IF ( ABS( urayt( 2 ) ) > EPSILON( h1 ) ) THEN
        IF        ( x(  2 ) < SSP%Seg%y( iSegy0     ) ) THEN
           h11 = -( x0( 2 ) - SSP%Seg%y( iSegy0     ) ) / urayt( 2 )
           ! write( *, * ) 'ocean segment crossing in y'
        ELSE IF   ( x(  2 ) > SSP%Seg%y( iSegy0 + 1 ) ) THEN
           h11 = -( x0( 2 ) - SSP%Seg%y( iSegy0 + 1 ) ) / urayt( 2 )
           ! write( *, * ) 'ocean segment crossing in y'

        END IF
     END IF
  END IF

  h = MIN( h, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11 )  ! take limit set by shortest distance to a crossing

  IF ( h < 1.0d-4 * deltas ) THEN   ! is it taking an infinitesimal step?
     h = 1.0d-5 * deltas     ! make sure we make some motion
     iSmallStepCtr = iSmallStepCtr + 1   ! keep a count of the number of sequential small steps
  ELSE
     iSmallStepCtr = 0   ! didn't do a small step so reset the counter
  END IF
END SUBROUTINE ReduceStep3D

! **********************************************************************!

SUBROUTINE Reflect3D( is, HS, BotTop, nBdry, RefC, Npts )

  ! *** If reflecting off surface need jump conditions ***

  USE bellhopMod
  USE RefCoMod
  USE cross_products
  USE norms
  USE sspMod
  IMPLICIT NONE
  INTEGER,           INTENT( IN    ) :: Npts                     ! Number of points in the reflection coefficient
  REAL    (KIND= 8), INTENT( INOUT ) :: nBdry( 3 )               ! normal to the boundary (changes if cone reflection)
  CHARACTER (LEN=3), INTENT( IN    ) :: BotTop                   ! bottom or top flag
  TYPE( HSInfo ),    INTENT( IN    ) :: HS                       ! halfspace parameters
  TYPE(ReflectionCoef), INTENT( IN ) :: RefC( NPts )             ! reflection coefficient
  INTEGER,           INTENT( INOUT ) :: is                       ! index of the ray step
  INTEGER          :: is1
  REAL    (KIND=8) :: rayt( 3 ), rayn1( 3 ), rayn2( 3 )             ! unit ray tangent and normals
  REAL    (KIND=8) :: rayt_tilde( 3 ), rayn1_tilde( 3 ), rayn2_tilde( 3 ), cn1jump, cn2jump, csjump 
  REAL    (KIND=8) :: c, cimag, gradc( 3 ), cxx, cyy, czz, cxy, cxz, cyz, rho   ! derivatives of sound speed in cartesian coordinates
  REAL    (KIND=8) :: cn1, cn2, cs                                  ! derivatives of sound speed in ray-centered coorindates
  REAL    (KIND=8) :: RM, R1, R2, R3, Tg, Th                            ! curvature corrections on reflection
  COMPLEX (KIND=8) :: gamma1, gamma2, gamma1Sq, gamma2Sq, GK, Refl
  TYPE(ReflectionCoef) :: RInt
  REAL    (KIND=8) :: tBdry( 3 )                                    ! tangent to the boundary
  REAL    (KIND=8) :: e1( 3 ), e2( 3 )                              ! ray normals for ray-centered coordinates
  REAL    (KIND=8) :: p_tilde_in(  2 ), p_hat_in(  2 ), q_tilde_in(  2 ), q_hat_in(  2 ), &
                      p_tilde_out( 2 ), p_hat_out( 2 )
  REAL    (KIND=8) :: Phi, theta, Radius, kappa,z_xx, z_xy, z_yy, x, y, &
       t_rot( 2 ), n_rot( 2 ), RotMat( 2, 2 ), kappaMat( 2, 2 ), DMat( 2, 2 )   ! for cone reflection

  is  = is + 1
  is1 = is + 1

  ! special case of a conical seamount

  kappa = 0
  z_xx = 0
  z_xy = 0
  z_yy = 0
  
!!$  IF ( BotTop == 'BOT' ) THEN
!!$     phi   = 15 * DegRad   ! 15 degree angle of seamount
!!$     theta = atan2(ray3D( is )%x( 2 ), ray3D( is )%x( 1 ) )   ! bearing from origin (cone location) to ray
!!$
!!$     nBdry =  [ -cos( theta ) * sin( phi ), -sin( theta ) * sin( phi ), cos( phi ) ]
!!$     Radius = NORM2B( ray3D( is )%x( 1 : 2 ) )   ! radius of seamount at the bounce point
!!$     kappa = -1 / Radius   ! curvature of seamount is the reciprocal of the radius at the bounce point
!!$
!!$     x = ray3D( is )%x( 1 )
!!$     y = ray3D( is )%x( 2 )
!!$     z_xx = ( 1 / Radius - x**2 / Radius**3 ) * tan( phi )
!!$     z_xy = - x * y / Radius**3 * tan( phi )
!!$     z_yy = ( 1 / Radius - y**2 / Radius**3 ) * tan( phi )
!!$  END IF

  ! z = z_xx / 2 * x^2 + z_xy * xy + z_yy * y^2; coefs are the kappa matrix
  kappaMat( 1, 1 ) = z_xx / 2
  kappaMat( 1, 2 ) = z_xy
  kappaMat( 2, 1 ) = z_xy
  kappaMat( 2, 2 ) = z_yy / 2
     
  Th = DOT_PRODUCT( ray3D( is )%t, nBdry )  ! component of ray tangent, normal to boundary

  ray3D( is1 )%NumTopBnc = ray3D( is )%NumTopBnc
  ray3D( is1 )%NumBotBnc = ray3D( is )%NumBotBnc
  ray3D( is1 )%x         = ray3D( is )%x
  ray3D( is1 )%t         = ray3D( is )%t - 2.0 * Th * nBdry   ! reflect the ray

  tBdry = ray3D( is )%t - Th * nBdry        ! tangent, along the boundary, in the reflection plane
  tBdry = tBdry / NORM2b( tBdry )
  Tg = DOT_PRODUCT( ray3D( is )%t, tBdry )  ! component of ray tangent, along the boundary

  ! Calculate the ray normals, rayn1, rayn2

  CALL EvaluateSSP3D( ray3D( is1 )%x, c, cimag, gradc, cxx, cyy, czz, cxy, cxz, cyz, rho, freq, 'TAB' )

  ! incident unit ray tangent and normal
  rayt  = c * ray3D( is )%t               ! unit tangent to ray
  rayn2 = -cross_product( rayt, nBdry )   ! ray tangent x boundary normal gives refl. plane normal
  rayn2 = rayn2 / NORM2b( rayn2 )         !!! does rayn2 need this normalization?
  rayn1 = -cross_product( rayt, rayn2 )   ! ray tangent x refl. plane normal is first ray normal

  !!! There is no jump in the normal to the reflection plane ???
  ! reflected unit ray tangent and normal (the reflected tangent, normal system has a different orientation)
  rayt_tilde  = c * ray3D( is1 )%t                         ! unit tangent to ray
  rayn2_tilde = -cross_product( rayt_tilde, nBdry )        ! ray tangent x boundary normal gives refl. plane normal
  rayn2_tilde = rayn2_tilde / NORM2b( rayn2_tilde )        !!! does rayn2 need this normalization?
  rayn1_tilde = -cross_product( rayt_tilde, rayn2_tilde )  ! ray tangent x refl. plane normal is first ray normal

  ! normal and tangential derivatives of the sound speed
  cn1jump = DOT_PRODUCT( gradc, rayn1_tilde - rayn1 )
  cn2jump = DOT_PRODUCT( gradc, rayn2_tilde - rayn2 )
  csjump  = DOT_PRODUCT( gradc, rayt_tilde  - rayt )

  ! rotation matrix to get surface curvature in and perpendicular to the reflection plane
  ! we use only the first two elements of the vectors because we want the projection in the x-y plane
  t_rot = rayn1( 1 : 2 ) / NORM2b( rayn1( 1 : 2 ) )
  n_rot = rayn2( 1 : 2 ) / NORM2b( rayn2( 1 : 2 ) )
  RotMat( 1 : 2, 1 ) = t_rot
  RotMat( 1 : 2, 2 ) = n_rot

  ! apply the rotation to get the matrix D of curvatures (see Popov 1977 for definition of DMat)
  DMat = MATMUL( kappaMat, RotMat )
  DMat = MATMUL( RotMat, DMat )

  !!! not sure if cn2 needs a sign flip also
  IF ( BotTop == 'TOP' ) THEN
     cn1jump = -cn1jump    ! flip sign for top reflection
     cn2jump = -cn2jump    ! flip sign for top reflection
  END IF

  ! Note that Tg, Th need to be multiplied by c to normalize tangent; hence, c^2 below
  RM = Tg / Th   ! this is tan( alpha ) where alpha is the angle of incidence
  !R1 = RM * ( 2 * cn1jump - RM * csjump ) / c ** 2
  !R2 = RM * cn2jump / c ** 2
  R1 = -2 / c ** 2 * DMat( 1, 1 ) / Th + RM * ( 2 * cn1jump - RM * csjump ) / c ** 2
  R2 = +2 / c *      DMat( 1, 2 )
  R3 = +2 *          DMat( 2, 2 ) * Th
 
  SELECT CASE ( Beam%Type( 2 : 2 ) )
  CASE ( 'D' )
     R1 = 2.0 * R1
     R2 = 2.0 * R2
     R3 = 2.0 * R3
  CASE ( 'Z' )
     R1 = 0.0
     R2 = 0.0
     R3 = 0.0
  END SELECT

  ! f, g, h continuation; needs curvature corrections
  ! ray3D( is1 )%f    = ray3D( is )%f ! + ray3D( is )%DetQ * R1
  ! ray3D( is1 )%g    = ray3D( is )%g ! + ray3D( is )%DetQ * R1
  ! ray3D( is1 )%h    = ray3D( is )%h ! + ray3D( is )%DetQ * R2
  ! ray3D( is1 )%DetP = ray3D( is )%DetP
  ! ray3D( is1 )%DetQ = ray3D( is )%DetQ

  ray3D( is1 )%c   = c
  ray3D( is1 )%tau = ray3D( is )%tau

  ! *** curvature correction ***

  CALL RayNormal( ray3D( is )%t, ray3D( is )%phi, ray3D( is )%c, e1, e2 )  ! Compute ray normals e1 and e2

  ! rotate p-q from e1, e2 system, onto rayn1, rayn2 system
 
  p_tilde_in = DOT_PRODUCT( rayn1, e1 ) * ray3D( is )%p_tilde + DOT_PRODUCT( rayn1, e2 ) * ray3D( is )%p_hat
  p_hat_in   = DOT_PRODUCT( rayn2, e1 ) * ray3D( is )%p_tilde + DOT_PRODUCT( rayn2, e2 ) * ray3D( is )%p_hat
  q_tilde_in = DOT_PRODUCT( rayn1, e1 ) * ray3D( is )%q_tilde + DOT_PRODUCT( rayn1, e2 ) * ray3D( is )%q_hat
  q_hat_in   = DOT_PRODUCT( rayn2, e1 ) * ray3D( is )%q_tilde + DOT_PRODUCT( rayn2, e2 ) * ray3D( is )%q_hat

  ! here's the actual curvature change
  p_tilde_out = p_tilde_in + q_tilde_in * R1 + q_hat_in * R2
  p_hat_out   = p_hat_in   - q_tilde_in * R2 + q_hat_in * R3

  ! rotate p-q back to e1, e2 system
  ray3d( is1 )%p_tilde = DOT_PRODUCT( rayn1, e1 ) * p_tilde_out + DOT_PRODUCT( rayn2, e1 ) * p_hat_out
  ray3d( is1 )%p_hat   = DOT_PRODUCT( rayn1, e2 ) * p_tilde_out + DOT_PRODUCT( rayn2, e2 ) * p_hat_out

  ray3D( is1 )%q_tilde = ray3D( is )%q_tilde
  ray3D( is1 )%q_hat   = ray3D( is )%q_hat

  ! Logic below fixes a bug when the dot product is infinitessimally greater than 1 (then ACos is complex)
  ray3D( is1 )%phi = ray3D( is )%phi + 2 * ACOS( MAX( MIN( DOT_PRODUCT( rayn1, e1 ), 1.0D0 ), -1.0D0 ) ) !!!What happens to torsion?

  ! account for phase change

  SELECT CASE ( HS%BC )
  CASE ( 'R' )                 ! rigid
     ray3D( is1 )%Amp   = ray3D( is )%Amp
     ray3D( is1 )%Phase = ray3D( is )%Phase
  CASE ( 'V' )                 ! vacuum
     ray3D( is1 )%Amp   = ray3D( is )%Amp
     ray3D( is1 )%Phase = ray3D( is )%Phase + pi
  CASE ( 'F' )                 ! file
     RInt%theta = RadDeg * ABS( ATAN2( Th, Tg ) )           ! angle of incidence (relative to normal to bathymetry)
     IF ( RInt%theta > 90 ) RInt%theta = 180. - RInt%theta  ! reflection coefficient is symmetric about 90 degrees
     CALL InterpolateReflectionCoefficient( RInt, RefC, Npts, PRTFile )
     ray3D( is1 )%Amp   = ray3D( is )%Amp * RInt%R
     ray3D( is1 )%Phase = ray3D( is )%Phase + RInt%phi
  CASE ( 'A', 'G' )            ! half-space
     GK       = omega * Tg     ! wavenumber in direction parallel to bathymetry
     gamma1Sq = ( omega / c     ) ** 2 - GK ** 2 - i * tiny( omega )   ! tiny prevents g95 giving -zero, and wrong branch cut
     gamma2Sq = ( omega / HS%cP ) ** 2 - GK ** 2 - i * tiny( omega )
     gamma1   = SQRT( -gamma1Sq )
     gamma2   = SQRT( -gamma2Sq )

     Refl = ( HS%rho * gamma1 - gamma2 ) / ( HS%rho * gamma1 + gamma2 )

     IF ( ABS( Refl ) < 1.0E-5 ) THEN   ! kill a ray that has lost its energy in reflection
        ray3D( is1 )%Amp   = 0.0
        ray3D( is1 )%Phase = ray3D( is )%Phase
     ELSE
        ray3D( is1 )%Amp   = ABS( Refl ) * ray3D(  is )%Amp
        ray3D( is1 )%Phase = ray3D( is )%Phase + ATAN2( AIMAG( Refl ), REAL( Refl ) )
     ENDIF
  END SELECT

END SUBROUTINE Reflect3D

!**********************************************************************!

SUBROUTINE RayNormal( t, phi, c, e1, e2 )

  ! Computes the ray normals

  USE norms
  IMPLICIT NONE
  REAL (KIND=8), INTENT( IN  ) :: t( 3 )             ! tangent vector (NOT normalized)
  REAL (KIND=8), INTENT( IN  ) :: phi                ! torsion
  REAL (KIND=8), INTENT( IN  ) :: c                  ! sound speed
  REAL (KIND=8), INTENT( OUT ) :: e1( 3 ), e2( 3 )   ! ray normals
  REAL (KIND=8)                :: RL                 ! length of part of the tangent vector

  RL = NORM2b( t( 1 : 2 ) )

  !  e1
  e1( 1 ) = ( c * t( 1 ) * t( 3 ) * COS( phi ) + t( 2 ) * SIN( phi ) ) / RL
  e1( 2 ) = ( c * t( 2 ) * t( 3 ) * COS( phi ) - t( 1 ) * SIN( phi ) ) / RL
  e1( 3 ) = -c * RL * COS( phi )

  !  e2
  e2( 1 ) = ( c * t( 1 ) * t( 3 ) * SIN( phi ) - t( 2 ) * COS( phi ) ) / RL
  e2( 2 ) = ( c * t( 2 ) * t( 3 ) * SIN( phi ) + t( 1 ) * COS( phi ) ) / RL
  e2( 3 ) = -c * RL * SIN( phi )

  RETURN
END SUBROUTINE RayNormal

!**********************************************************************!

SUBROUTINE RayNormal_unit( t, phi, e1, e2 )

  ! Computes the ray normals
  ! This version assumes t is a unit normal

  USE norms
  IMPLICIT NONE
  REAL (KIND=8), INTENT( IN  ) :: t( 3 )             ! tangent vector (normalized)
  REAL (KIND=8), INTENT( IN  ) :: phi                ! torsion
  REAL (KIND=8), INTENT( OUT ) :: e1( 3 ), e2( 3 )   ! ray normals
  REAL (KIND=8)                :: RL                 ! length of part of the tangent vector

  RL = NORM2b( t( 1 : 2 ) )

  !  e1
  e1( 1 ) = ( t( 1 ) * t( 3 ) * COS( phi ) + t( 2 ) * SIN( phi ) ) / RL
  e1( 2 ) = ( t( 2 ) * t( 3 ) * COS( phi ) - t( 1 ) * SIN( phi ) ) / RL
  e1( 3 ) = -RL * COS( phi )

  !  e2
  e2( 1 ) = ( t( 1 ) * t( 3 ) * SIN( phi ) - t( 2 ) * COS( phi ) ) / RL
  e2( 2 ) = ( t( 2 ) * t( 3 ) * SIN( phi ) + t( 1 ) * COS( phi ) ) / RL
  e2( 3 ) = -RL * SIN( phi )

  RETURN
END SUBROUTINE RayNormal_unit
!**********************************************************************!

SUBROUTINE Get_c_partials( cxx, cxy, cxz, cyy, cyz, czz, e1, e2, cnn, cmn, cmm )

  ! Computes the second partials of c along ray normals

  IMPLICIT NONE
  REAL (KIND=8), INTENT( IN  ) :: cxx, cxy, cxz, cyy, cyz, czz  ! curvature of sound speed (cartesian)
  REAL (KIND=8), INTENT( IN  ) :: e1( 3 ), e2( 3 )              ! principal normals
  REAL (KIND=8), INTENT( OUT ) :: cnn, cmn, cmm                 ! curvature of sound speed (ray-centered)

  cnn = cxx * e1( 1 )**2 + cyy * e1( 2 )**2 + czz * e1( 3 )**2 + 2.0 * cxy * e1( 1 ) * e1( 2 ) + &
        2.0 * cxz * e1( 1 ) * e1( 3 ) + 2.0 * cyz * e1( 2 ) * e1( 3 )

  cmn = cxx * e1( 1 ) * e2( 1 ) + cyy * e1( 2 ) * e2( 2 ) + czz * e1( 3 ) * e2( 3 ) + &
        cxy * ( e1( 1 ) * e2( 2 ) + e2( 1 ) * e1( 2 ) ) + cxz * ( e1( 1 ) * e2( 3 ) + e2( 1 ) * e1( 3 ) ) +  &
        cyz * ( e1( 2 ) * e2( 3 ) + e2( 2 ) * e1( 3 ) )

  cmm = cxx * e2( 1 )**2 + cyy * e2( 2 )**2 + czz * e2( 3 )**2 + 2.0 * cxy * e2( 1 ) * e2( 2 ) + &
        2.0 * cxz * e2( 1 ) * e2( 3 ) + 2.0 * cyz * e2( 2 ) * e2( 3 )

  RETURN
END SUBROUTINE Get_c_partials
