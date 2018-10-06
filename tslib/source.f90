SUBROUTINE SOURCE( T, S, SD, Nsd, Ntout, omega, fmin, fmax, Pulse, PulseTitle, IniFlag )

  ! Evaluate the source time series

  ! stemp is a temporary workspace

  ! IniFlag, used to control filtering, has to be controlled on the outside
  ! since SPARC changes the filtering as a function of wavenumber.

  IMPLICIT NONE
  INTEGER, PARAMETER :: MaxNt = 1000000, MaxNST = 1000000
  INTEGER            :: it, Nsd, Ntout, Nt
  LOGICAL   IniFlag
  REAL      T( * ), TF( MaxNt ), SD( Nsd ), stempR( MaxNt ), deltat, omega, fmin, fmax, TSTART
  COMPLEX   S( Nsd, * ), stemp( MaxNt )
  CHARACTER (LEN=4) :: Pulse
  CHARACTER PulseTitle*( * )
  EQUIVALENCE ( stempR, stemp )
  COMPLEX, ALLOCATABLE :: SF( :, : )

  SAVE TF, SF, Nt

  ! If first call, tabulate the time series

  IF ( IniFlag ) THEN
     IF ( .NOT. ALLOCATED( SF ) ) ALLOCATE( SF( Nsd, MaxNST / Nsd ) )
     IF ( Pulse(1:1) == 'F' .OR. Pulse(1:1) == 'B' ) THEN
        ! From a file
        CALL STSHDR( PulseTitle, SD, Nsd )
        CALL SFILE( Pulse, PulseTitle, Nsd, TF, stempR, SF, Nt, MaxNt, MaxNST )
        deltat = TF( 2 ) - TF( 1 )
     ELSE
        ! Use one of the canned wavelets
        Nsd    = 1
        deltat = 0.2 / omega
        TSTART = -200.0 * deltat
        Nt     = 1024
        DO it = 1, Nt
           TF( it ) = TSTART + ( it - 1 ) * deltat
           CALL CANS( TF( it ), omega, Pulse, SF, Nsd, it, PulseTitle )
        END DO
     ENDIF

     ! Filter the time series
     CALL FILTER( Pulse, deltat, stemp, SF, Nsd, Nt, fmin, fmax )

     IniFlag = .FALSE.
  ENDIF

  CALL EVALU( T, Pulse, TF, SF, Nsd, Nt, S, Ntout )

END SUBROUTINE SOURCE
!**********************************************************************C
SUBROUTINE STSHDR( PulseTitle, SD, Nsd )

  !     Read header from time series file

  IMPLICIT NONE
  INTEGER, PARAMETER :: STSFIL = 10
  INTEGER, INTENT( INOUT ) :: Nsd
  REAL,    INTENT( OUT   ) :: SD( Nsd )
  CHARACTER PulseTitle*( * )

  ! *** Read in time series data ***

  OPEN( FILE = 'STSFIL', UNIT = STSFIL, STATUS = 'OLD', FORM = 'FORMATTED' )
  READ( STSFIL, * ) PulseTitle

  ! SD is used by PLOTTS
  READ( STSFIL, * ) Nsd, SD

END SUBROUTINE STSHDR

!**********************************************************************C

SUBROUTINE SFILE( Pulse, PulseTitle, Nsd, TF, stemp, SF, Nt, MaxNt, MaxNST )

  ! Time series from file
  ! stemp is used to hold a temporary real vector

  IMPLICIT NONE
  INTEGER, PARAMETER :: PRTFIL = 6, STSFIL = 10
  INTEGER            :: it, Nsd, Nt, MaxNt, MaxNST
  REAL      stemp( Nsd ), TF( * ), tempV( Nsd ), temp, TMax
  COMPLEX   SF( Nsd, * )
  CHARACTER (LEN=4) :: Pulse
  CHARACTER PulseTitle*( * )

  ! Read in time series data

  Nt = 0
  DO it = 1, MaxNt

     READ( STSFIL, *, END = 2000 ) TF( it ), stemp

     IF ( it > 1 ) THEN
        IF ( TF( it ) < TF( it - 1 ) ) CALL ERROUT( PRTFIL, 'F', 'SOURCE:SFILE', 'Time series not ordered in time' )
     ENDIF

     SF( :, it ) = stemp
     Nt = Nt + 1

     IF ( Nt * Nsd > MaxNST ) CALL ERROUT( PRTFIL, 'F', 'SOURCE:SFILE', 'Too many time series points' )
  END DO

  ! normal exit is by EOF in file, not by falling through loop
  CALL ERROUT( PRTFIL, 'F', 'SOURCE:SFILE', 'Too many time series points' )

  ! Time reversal
  ! note if Nt is odd, middle point doesn't move

2000 IF ( Pulse(1:1) == 'B' ) THEN
     TMax = TF( Nt )
     DO it = 1, Nt / 2
        tempV                = REAL( SF( :, it ) )
        SF( :, it )          = SF( :, Nt - it + 1 )
        SF( :, Nt - it + 1 ) = tempV

        temp                  = TMax - TF( it )
        TF( it )              = TMax - TF( Nt - it + 1 )
        TF( Nt - it + 1 )     = temp
     END DO
  ENDIF

  CLOSE( STSFIL )

END SUBROUTINE SFILE
!**********************************************************************C

SUBROUTINE EVALU( T, Pulse, TF, SF, Nsd, Nt, S, Ntout )

  ! *** Evaluate the source function ***

  REAL      T( * ), TF( * )
  COMPLEX   S( Nsd, * ), SF( Nsd, * )
  CHARACTER (LEN=4) :: Pulse

  SAVE it
  DATA it /1/

  DO itout = 1, Ntout

     ! Necessary to advance time window?

     IF ( T( itout ) < TF( it ) ) it = 1 ! Case of T reset to origin for new march
     DO WHILE ( T( itout ) > TF( it + 1 ) .AND. it < Nt - 1 )
        it = it + 1
     END DO

     ! Linear interpolation in time

     S( :, itout ) = 0.0   ! in case T outside time window

     IF ( T( itout ) >= TF( 1  ) .AND. T( itout ) <= TF( Nt ) ) THEN
        W = ( T( itout ) - TF( it ) ) / ( TF( it + 1 ) - TF( it ) )
        S( :, itout ) = SF( :, it ) + W * ( SF( :, it+1 ) - SF( :, it ) )
        IF ( Pulse( 3 : 3 ) == '-' ) S( :, itout ) = -S( :, itout )
     ENDIF

  END DO   ! next itout

END SUBROUTINE EVALU

!**********************************************************************C

SUBROUTINE FILTER( Pulse, deltat, stemp, SF, Nsd, Nt, fmin, fmax )

  ! Filter the source function
  ! Nt must be a power or 2

  IMPLICIT NONE
  INTEGER   is, Nsd, Nt
  REAL      deltat, fmin, fmax
  COMPLEX   stemp( Nt ), SF( Nsd, Nt )
  CHARACTER (LEN=4) :: Pulse

  DO is = 1, Nsd
     stemp = SF( is, : ) !  Copy time series into stemp         
     IF ( Pulse( 4 : 4 ) /= 'N' ) CALL BNDPAS(  stemp, Nt, deltat, fmin, fmax )! Band-pass filter
     IF ( Pulse( 2 : 2 ) == 'H' ) CALL PREENV(  stemp, Nt ) ! Form the pre-envelope
     IF ( Pulse( 2 : 2 ) == 'Q' ) CALL HILBERT( stemp, Nt ) ! Hilbert transform
     SF( is, 1:Nt ) = stemp( 1:Nt ) ! Copy time series back
  END DO

END SUBROUTINE FILTER
