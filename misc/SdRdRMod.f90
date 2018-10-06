MODULE SdRdRMod
  ! Reads in source depths, receiver depths, receiver ranges, and receiver bearings

  USE monotonicMod
  USE SortMod
  USE SubTabulate

  IMPLICIT NONE
  SAVE
  INTEGER, PARAMETER         :: Number_to_Echo = 21
  INTEGER                    :: Nfreq          ! number of frequencies
  REAL (KIND=8), ALLOCATABLE :: freqVec( : )   ! frequency vector for braodband runs

  TYPE Position
     INTEGER              :: Nsx, Nsy, Nsd, Nrd, Nr, Ntheta   ! number of x, y, z, r, theta coordinates
     REAL                 :: Delta_r, Delta_theta
     INTEGER, ALLOCATABLE :: isd( : ), ird( : )
     REAL,    ALLOCATABLE :: sx( : ), sy( : )   ! source x, y coordinates
     REAL,    ALLOCATABLE :: sd( : ), rd( : ), ws( : ), wr( : ), r( : )
     REAL,    ALLOCATABLE :: theta( : )         ! receiver bearings
  END TYPE Position

  TYPE (Position ) :: Pos   ! structure containing source and receiver positions

CONTAINS

  SUBROUTINE ReadfreqVec( ENVFile, PRTFile, freq, BroadbandOption )

    ! Optionally reads a vector of source frequencies for a broadband run
    ! allocating and creating a frequency vector
    !
    ! If the broadband option is not selected, then the input freq is stored in the frequency vector

    IMPLICIT NONE
    INTEGER,       INTENT( IN ) :: ENVFile, PRTFile
    REAL (KIND=8), INTENT( IN ) :: freq             ! default frequency
    CHARACTER,     INTENT( IN ) :: BroadbandOption*( 1 )
    INTEGER                     :: IAllocStat, ifreq

    Nfreq = 1

    IF ( BroadbandOption == 'B' ) THEN
       READ( ENVFile, * ) Nfreq   ! Broadband run
       WRITE( PRTFile, * ) '__________________________________________________________________________'
       WRITE( PRTFile, * )
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Number of frequencies =', Nfreq
    END IF

    IF ( Nfreq <= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadEnvironment', 'Number of frequencies must be positive'  )

    IF ( ALLOCATED( FreqVec ) ) DEALLOCATE( FreqVec )
    ALLOCATE( FreqVec( MAX( 3, Nfreq ) ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadEnvironment', 'Too many frequencies'  )

    IF ( BroadbandOption == 'B' ) THEN
       WRITE( PRTFile, * ) 'Frequencies (Hz)'
       FreqVec( 3 ) = -999.9
       READ(  ENVFile, * ) FreqVec( 1 : Nfreq )
       CALL SubTab( FreqVec, Nfreq )

       WRITE( PRTFile, "( 5G14.6 )" ) ( freqVec( ifreq ), ifreq = 1, MIN( Nfreq, Number_to_Echo ) )
       IF ( Nfreq > Number_to_Echo ) WRITE( PRTFile, * ) ' ... ', freqVec( Nfreq )
    ELSE
       freqVec( 1 ) = freq
    END IF

    RETURN

  END SUBROUTINE ReadfreqVec

  !********************************************************************!

  SUBROUTINE Readsxsy( ENVFile, PRTFile, ThreeD )

    IMPLICIT NONE
    LOGICAL, INTENT( IN ) :: ThreeD   ! flag indicating whether this is a 3D run
    INTEGER, INTENT( IN ) :: ENVFile, PRTFile
    INTEGER               :: is, IAllocStat

    IF ( ThreeD ) THEN

       ! *** Read source x coordinates ***

       READ(  ENVFile, * ) Pos%Nsx
       WRITE( PRTFile, * ) '__________________________________________________________________________'
       WRITE( PRTFile, * )
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Number of source x coordinates = ', Pos%Nsx

       IF ( Pos%Nsx <= 0 ) CALL ERROUT( PRTFile, 'F', 'Readsxsy', 'Number of source x coordinates must be positive' )

       IF ( ALLOCATED( Pos%sx ) ) DEALLOCATE( Pos%sx )
       ALLOCATE( Pos%sx( MAX( 3, Pos%Nsx ) ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'Readsxsy', 'Too many sources'  )

       WRITE( PRTFile, * ) 'Source x coordinate (km)'
       Pos%sx( 3 ) = -999.9
       READ( ENVFile, * ) Pos%sx( 1 : Pos%Nsx )

       CALL SubTab( Pos%sx, Pos%Nsx )
       !CALL SORT(   Pos%sx, Pos%Nsx )

       WRITE( PRTFile, "( 5G14.6 )" ) ( Pos%sx( is ), is = 1, MIN( Pos%Nsx, Number_to_Echo ) )
       IF ( Pos%Nsx > Number_to_Echo ) WRITE( PRTFile, * ) ' ... ', Pos%sx( Pos%Nsx )

       ! *** Read source y coordinates ***

       READ(  ENVFile, * ) Pos%Nsy
       WRITE( PRTFile, * )
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Number of source y coordinates = ', Pos%Nsy

       IF ( Pos%Nsy <= 0 ) CALL ERROUT( PRTFile, 'F', 'Readsxsy', 'Number of source y coordinates must be positive' )

       IF ( ALLOCATED( Pos%sy ) ) DEALLOCATE( Pos%sy )
       ALLOCATE( Pos%sy( MAX( 3, Pos%Nsy ) ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'Readsxsy', 'Too many sources'  )

       WRITE( PRTFile, * ) 'Source y coordinate (km)'
       Pos%sy( 3 ) = -999.9
       READ( ENVFile, * ) Pos%sy( 1 : Pos%Nsy )

       CALL SubTab( Pos%sy, Pos%Nsy )
       !CALL SORT(   Pos%sy, Pos%Nsy )

       WRITE( PRTFile, "( 5G14.6 )" ) ( Pos%sy( is ), is = 1, MIN( Pos%Nsy, Number_to_Echo ) )
       IF ( Pos%Nsy > Number_to_Echo ) WRITE( PRTFile, * ) ' ... ', Pos%sy( Pos%Nsy )

       Pos%sx = 1000.0 * Pos%sx   ! convert km to m
       Pos%sy = 1000.0 * Pos%sy

!!$       IF ( .NOT. monotonic( Pos%sx, Pos%Nsx ) ) THEN
!!$          CALL ERROUT( PRTFile, 'F', 'SdRdRMod', 'Source x-coordinates are not monotonically increasing' )
!!$       END IF 
!!$ 
!!$       IF ( .NOT. monotonic( Pos%sy, Pos%Nsy ) ) THEN
!!$          CALL ERROUT( PRTFile, 'F', 'SdRdRMod', 'Source y-coordinates are not monotonically increasing' )
!!$       END IF 
    ELSE
       Pos%Nsx = 1
       Pos%Nsy = 1
       ALLOCATE( Pos%sx( 1 ), Pos%sy( 1 ) )
       Pos%sx( 1 ) = 0.
       Pos%sy( 1 ) = 0.
    END IF

    RETURN
  END SUBROUTINE Readsxsy

  !********************************************************************!

  SUBROUTINE ReadSdRd( ENVFile, PRTFile, zMin, zMax )

    ! Reads source and receiver depths
    ! zMin and zMax are limits for those depths; sources and receivers are shifted to be within those limits

    IMPLICIT NONE
    INTEGER, INTENT( IN ) :: ENVFile, PRTFile
    REAL,    INTENT( IN ) :: zMin, zMax
    !LOGICAL               :: monotonic
    INTEGER               :: is, ir, IAllocStat

    ! *** Read source depths ***

    READ(  ENVFile, * ) Pos%Nsd
    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) '__________________________________________________________________________'
    WRITE( PRTFile, * )
    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) 'Number of source   depths = ', Pos%Nsd

    IF ( Pos%Nsd <= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadSdRd', 'Number of sources must be positive'  )

    IF ( ALLOCATED( Pos%sd ) ) DEALLOCATE( Pos%sd, Pos%ws, Pos%isd )
    ALLOCATE( Pos%sd( MAX( 3, Pos%Nsd ) ), Pos%ws( Pos%Nsd ), Pos%isd( Pos%Nsd ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadSdRd', 'Too many sources'  )

    WRITE( PRTFile, * ) 'Source depths (m)'
    Pos%sd( 3 ) = -999.9
    READ( ENVFile, * ) Pos%sd( 1 : Pos%Nsd )

    CALL SubTab( Pos%sd, Pos%Nsd )
    !CALL SORT(   Pos%sd, Pos%Nsd )

    WRITE( PRTFile, "( 5G14.6 )" ) ( Pos%sd( is ), is = 1, MIN( Pos%Nsd, Number_to_Echo ) )
    IF ( Pos%Nsd > Number_to_Echo ) WRITE( PRTFile, * ) ' ... ', Pos%sd( Pos%Nsd )

    ! *** Read receiver depths ***

    READ(  ENVFile, * ) Pos%Nrd
    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) 'Number of receiver depths = ', Pos%Nrd

    IF ( Pos%Nrd <= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadSdRd', 'Number of receivers must be positive'  )

    IF ( ALLOCATED( Pos%rd ) ) DEALLOCATE( Pos%rd, Pos%wr, Pos%ird )
    ALLOCATE( Pos%rd( MAX( 3, Pos%Nrd ) ), Pos%wr( Pos%Nrd ), Pos%ird( Pos%Nrd ), Stat = IAllocStat  )
    IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadSdRd', 'Too many receivers'  )

    WRITE( PRTFile, * ) 'Receiver depths (m)'
    Pos%rd( 3 ) = -999.9
    READ( ENVFile, * ) Pos%rd( 1 : Pos%Nrd )

    CALL SubTab( Pos%rd, Pos%Nrd )
    !CALL SORT(   Pos%rd, Pos%Nrd )

    WRITE( PRTFile, "( 5G14.6 )" ) ( Pos%rd( ir ), ir = 1, MIN( Pos%Nrd, Number_to_Echo ) )
    IF ( Pos%Nrd > Number_to_Echo ) WRITE( PRTFile, * ) ' ... ', Pos%rd( Pos%Nrd )

    ! *** Check for sd/rd in upper or lower halfspace ***

    IF ( ANY( Pos%sd( 1 : Pos%Nsd ) < zMin ) ) THEN
       WHERE ( Pos%sd < zMin ) Pos%sd = zMin
       CALL ERROUT( PRTFile, 'W', 'SdRdRMod', 'Source above or too near the top bdry has been moved down' )
    END IF

    IF ( ANY( Pos%sd( 1 : Pos%Nsd ) > zMax ) ) THEN
       WHERE( Pos%sd > zMax ) Pos%sd = zMax
       CALL ERROUT( PRTFile, 'W', 'SdRdRMod', 'Source below or too near the bottom bdry has been moved up' ) 
    END IF

    IF ( ANY( Pos%rd( 1 : Pos%Nrd ) < zMin ) ) THEN
       WHERE( Pos%rd < zMin ) Pos%rd = zMin
       CALL ERROUT( PRTFile, 'W', 'SdRdRMod', 'Receiver above or too near the top bdry has been moved down' ) 
    END IF

    IF ( ANY( Pos%rd( 1 : Pos%Nrd ) > zMax ) ) THEN
       WHERE( Pos%rd > zMax ) Pos%rd = zMax
       CALL ERROUT( PRTFile, 'W', 'SdRdRMod', 'Receiver below or too near the bottom bdry has been moved up' ) 
    END IF

!!$    IF ( .NOT. monotonic( Pos%sd, Pos%Nsd ) ) THEN
!!$       CALL ERROUT( PRTFile, 'F', 'SdRdRMod', 'Source depths are not monotonically increasing' )
!!$    END IF 
!!$ 
!!$    IF ( .NOT. monotonic( Pos%rd, Pos%Nrd ) ) THEN
!!$       CALL ERROUT( PRTFile, 'F', 'SdRdRMod', 'Receiver depths are not monotonically increasing' )
!!$    END IF 

    RETURN
  END SUBROUTINE ReadSdRd

  !********************************************************************!

  SUBROUTINE ReadRcvrRanges( ENVFile, PRTFile )

    ! Read receiver ranges

    IMPLICIT NONE
    INTEGER, INTENT( IN ) :: ENVFile, PRTFile
    INTEGER               :: ir, IAllocStat

    READ(  ENVFile, * ) Pos%Nr
    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) 'Number of receiver ranges = ', Pos%Nr

    IF ( ALLOCATED( Pos%r ) ) DEALLOCATE( Pos%r )
    ALLOCATE( Pos%r( MAX( 3, Pos%Nr ) ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadRcvrRanges', 'Too many range points' )

    WRITE( PRTFile, * ) 'Receiver ranges (km)'
    Pos%r( 3 ) = -999.9
    READ( ENVFile, * ) Pos%r( 1 : Pos%Nr )

    CALL SubTab( Pos%r, Pos%Nr )
    CALL Sort(   Pos%r, Pos%Nr )

    WRITE( PRTFile, "( 5G14.6 )" ) ( Pos%r( ir ), ir = 1, MIN( Pos%Nr, Number_to_Echo ) )
    IF ( Pos%Nr > Number_to_Echo ) WRITE( PRTFile, * ) ' ... ', Pos%r( Pos%Nr )

    Pos%r( 1 : Pos%Nr ) = 1000.0 * Pos%r( 1 : Pos%Nr )   ! Convert ranges to meters

    ! calculate range spacing
    Pos%delta_r = 0.0
    IF ( Pos%Nr /= 1 ) Pos%delta_r = Pos%r( Pos%Nr ) - Pos%r( Pos%Nr - 1 )

    ! For a point source can't have receiver at origin
    ! IF ( OPT( 1 : 1 ) == 'R' .AND. Pos%r( 1 ) <= 0.0 ) 
    ! IF ( Pos%r( 1 ) <= 0.0 ) Pos%r( 1 ) = MIN( 1.0, Pos%r( 2 ) )

    IF ( .NOT. monotonic( Pos%r, Pos%Nr ) ) THEN
       CALL ERROUT( PRTFile, 'F', 'SdRdRMod', 'Receiver ranges are not monotonically increasing' )
    END IF 
 
    RETURN
  END SUBROUTINE ReadRcvrRanges

  !********************************************************************!

  SUBROUTINE ReadRcvrBearings( ENVFile, PRTFile )

    ! Read receiver bearings

    IMPLICIT NONE
    INTEGER, INTENT( IN ) :: ENVFile, PRTFile
    INTEGER               :: itheta, IAllocStat

    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) '__________________________________________________________________________'
    WRITE( PRTFile, * )

    READ(  ENVFile, * ) Pos%Ntheta
    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) 'Number of receiver bearings = ', Pos%Ntheta
    WRITE( PRTFile, * ) 'Receiver bearings (degrees)'

    IF ( ALLOCATED( Pos%theta ) ) DEALLOCATE( Pos%theta )
    ALLOCATE( Pos%theta( MAX( 3, Pos%Ntheta ) ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadBearings', 'Too many bearing angles' )

    Pos%theta( 3 ) = -999.9
    READ( ENVFile, * ) Pos%theta( 1 : Pos%Ntheta )

    CALL SubTab( Pos%theta, Pos%Ntheta )
    CALL Sort(   Pos%theta, Pos%Ntheta )

    WRITE( PRTFile, "( 5G14.6 )" ) ( Pos%theta( itheta ), itheta = 1, MIN( Pos%Ntheta, Number_to_Echo ) )
    IF ( Pos%Ntheta > Number_to_Echo ) WRITE( PRTFile, * ) ' ... ', Pos%theta( Pos%Ntheta )

    ! full 360-degree sweep? remove duplicate angle
    IF ( Pos%theta( Pos%Ntheta ) == Pos%theta( 1 ) + 360.0D0 ) Pos%Ntheta = Pos%Ntheta - 1

    ! calculate angular spacing
    Pos%Delta_theta = 0.0
    IF ( Pos%Ntheta /= 1 ) Pos%Delta_theta = Pos%theta( Pos%Ntheta ) - Pos%theta( Pos%Ntheta - 1 )

    IF ( .NOT. monotonic( Pos%theta, Pos%Ntheta ) ) THEN
       CALL ERROUT( PRTFile, 'F', 'SdRdRMod', 'Receiver bearings are not monotonically increasing' )
    END IF 
 
    RETURN
  END SUBROUTINE ReadRcvrBearings

END MODULE SdRdRMod
