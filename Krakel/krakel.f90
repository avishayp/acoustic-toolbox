PROGRAM KRAKEL

  ! Program for solving for ocean acoustic normal modes

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

  ! Originally developed as part of the author's dissertation under the supervision
  ! of Prof. Edward L. Reiss, Northwestern University

  USE krakelmod
  IMPLICIT NONE
  INTEGER            :: IFirst, ILast, IREC
  REAL (KIND=8)      :: omega, Error
  CHARACTER (LEN=80) :: FileRoot

  ! get the file root for naming all input and output files
  ! should add some checks here ...

  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )

  NV( 1 : 5 ) = [ 1, 2, 4, 8, 16 ]
  CALL GETPAR( FileRoot )
  M = MaxM

  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) 'Mesh multiplier   CPU seconds'

  ! Main loop: solve the problem for a sequence of meshes

  DO iSet = 1, NSETS

     N( 1 : NMedia ) = NG( 1 : NMedia ) * NV( iSet )
     h( 1 : NMedia ) = ( Depth( 2 : NMedia + 1 ) - Depth( 1 : NMedia ) ) / N( 1 : NMedia )
     hV( iSet ) = h( 1 )

     CALL SOLVE( FileRoot, Error )
     IF ( Error * 1000.0 * RMax < 1.0 ) GOTO 3000
  END DO

  ! Fall through indicates failure to converge
  WRITE( PRTFile, * ) 'WARNING: Too many meshes needed, check convergence'

  ! Solution complete: discard modes with phase velocity above cHigh

3000 omega = SQRT( omega2 )
  DO Mode = M, 1, -1
     IF ( omega / DBLE( SQRT( Extrap( 1, Mode ) ) ) < cHigh ) GOTO 5000
  END DO
5000 M = Mode

  ! Write eigenvalues to PRTFile and MODFile
  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) '   I          k             alpha          Phase Speed'

  k( 1 : M ) = SQRT( Extrap( 1, 1 : M ) )
  DO Mode = 1, M
     WRITE( PRTFile, "( I5, 3G18.10 )" ) Mode, k( Mode ), omega / DBLE( k( Mode ) )
  END DO

  WRITE( MODFile, REC=5 ) M, LRECL
  IFirst = 1
  DO IREC = 1, 1 + ( 2 * M - 1 ) / LRECL

     ILast  = MIN( M, IFirst + LRECL / 2 - 1 )
     WRITE( MODFile, REC = 6 + M + IREC ) CMPLX( k( IFirst : ILast ) )                            

     IFirst = ILast + 1
  END DO

END PROGRAM KRAKEL
!**********************************************************************C
SUBROUTINE GETPAR( FileRoot )

  ! Read in the ENVFile data

  USE krakelmod
  USE SdRdRMod
  USE RefCoMod
  IMPLICIT NONE
  REAL               :: zMin, zMax
  CHARACTER (LEN=80) :: FileRoot

  ! Read in environmental info

  TITLE = 'KRAKEL-  '

  CALL READIN( FileRoot, Title, Freq, MaxMedium, NMedia, TopOpt, HSTop, NG, Sigma, Depth, BotOpt, HSBot, ENVFile, PRTFile )

  READ(  ENVFile, *    ) cLow, cHigh                 ! Spectral limits
  WRITE( PRTFile, "( /, ' cLow = ', G12.5, 'm/s      cHigh = ', G12.5, 'm/s' )" ) cLow, cHigh
  IF ( cLow <= 0.0 .OR. cHigh <= 0.0 .OR. cLow >= cHigh ) &
       CALL ERROUT( PRTFile, 'F', 'KRAKEL-GETPAR', 'Need phase speeds cLow, cHigh > 0 and cLow < cHigh'  )

  READ(  ENVFile, * ) RMax                           ! Maximum range for calculations
  WRITE( PRTFile, * ) 'RMax = ', RMax

  ! Read source/receiver depths
  zMin = SNGL( Depth( 1 ) )
  zMax = SNGL( Depth( NMedia + 1 ) )
  CALL ReadSdRd( ENVFile, PRTFile, zMin, zMax )

  CLOSE( ENVFile )
  omega2 = ( 2.0 * pi * Freq ) ** 2

  CALL ReadReflectionCoefficient( FileRoot, HSBot%BC, HSTop%BC, PRTFile )   ! Optionally read in Bot, Top reflection coefficients
  M = MaxM

END SUBROUTINE GETPAR
!**********************************************************************C
SUBROUTINE INIT

  ! Initializes arrays defining difference equations
  ! Identical to KRAKENC apart for B1, B2, B3, B4 and rho

  USE krakelmod
  IMPLICIT NONE
  LOGICAL          :: ElasticFlag = .FALSE.
  INTEGER          :: IAllocStat, ii, j, Medium, N1, NPoints
  COMPLEX (KIND=8), ALLOCATABLE :: cP( : ), cS( : )
  COMPLEX (KIND=8) :: lambda, mu

  cMin     = HUGE( cMin )
  NFACT    = 0
  Loc( 1 ) = 0

  ! Allocate storage for finite-difference coefficients

  NPoints = SUM( N( 1 : NMedia ) ) + NMedia

  IF ( ALLOCATED( B1 ) ) DEALLOCATE( B1, B1C, B2, B3, B4, rho )
  ALLOCATE ( B1( NPoints ), B1C( NPoints ), B2( NPoints ), B3( NPoints ), B4( NPoints ), rho( NPoints ), &
       cP( NPoints ), cS( NPoints ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) &
       CALL ERROUT( PRTFile, 'F', 'KRAKEN - INIT', 'Insufficient memory: Reduce mesh.' )

  j = 0

  DO Medium = 1, NMedia   ! Loop over media
     IF ( Medium /= 1 ) Loc( Medium ) = Loc( Medium - 1 ) + N( Medium - 1 ) + 1
     N1 = N( Medium ) + 1

     IF ( Loc( Medium ) + N1 > MaxN ) THEN
        WRITE( PRTFile, * ) 'Insufficient storage for next mesh'
        WRITE( PRTFile, * ) 'PROGRAM ABORTING'
        STOP 'ERROR IN KRAKENC: Insufficient storage for next mesh'
     ENDIF

     CALL PROFIL( Depth, cP, cS, rho( Loc( Medium ) + 1 ), Medium, N1, FREQ, TopOpt( 1 : 1 ), TopOpt( 3 : 3 ), 'TAB' )

     IF ( cS( 1 ) == ( 0.0, 0.0 ) ) THEN

        ! Case of an acoustic medium

        Material( Medium ) = 'ACOUSTIC'
        IF ( NFACT == 0 ) NFACT = Medium
        NLACT = Medium

        DO ii = 1, N1
           cMin = MIN( DBLE( cP( ii ) ), cMin )
           j = j + 1
           B1( j ) = -2.0 +  h( Medium ) ** 2 * omega2 / cP( ii ) ** 2
        END DO
     ELSE

        ! Case of an elastic medium

        IF ( sigma( Medium ) /= 0.0 ) THEN
           WRITE( PRTFile, * ) 'Rough elastic interface not allowed'
           WRITE( PRTFile, * ) 'PROGRAM ABORTING'
           STOP 'ERROR IN KRAKENC: Rough elastic interface not allowed'
        ENDIF
        Material( Medium ) = 'ELASTIC'
        ElasticFlag        = .TRUE.

        DO ii = 1, N1
           cMin   = MIN( DBLE( cS( ii ) ), cMin )
           j      = j + 1
           lambda = rho( j ) * ( cP( ii ) ** 2 - 2.0 * cS( ii ) ** 2 )
           mu     = rho( j ) * cS( ii ) ** 2

           B1( j  ) = 1.0 / mu
           B2( j  ) = 1.0 / ( lambda + 2.0 * mu )
           B3( j  ) = 4.0 * mu * ( lambda + mu ) / ( lambda + 2.0 * mu )
           B4( j  ) = lambda / ( lambda + 2.0 * mu )
           rho( j ) = omega2 * rho( j )

        END DO

     ENDIF
  END DO ! Next Medium:

  ! Bottom properties

  SELECT CASE ( HSBot%BC )
  CASE ( 'A' )
     IF ( HSBot%cS /= ( 0.0, 0.0 ) ) THEN    ! Elastic  bottom half-space
        Material( NMedia + 1 ) = 'ELASTIC'
        ElasticFlag = .TRUE.
        cMin = MIN( cMin, DBLE( HSBot%cS ) )
     ELSE                                    ! Acoustic bottom half-space
        Material( NMedia + 1 ) = 'ACOUSTIC'
        cMin = MIN( cMin, DBLE( HSBot%cP ) )
     ENDIF
  END SELECT

  ! Top properties

  SELECT CASE ( HSTop%BC )
  CASE ( 'A' )
     IF ( HSTop%cS /= ( 0.0, 0.0 ) ) THEN    ! Elastic  top half-space
        ElasticFlag = .TRUE.
        cMin = MIN( cMin, DBLE( HSTop%cS ) )
     ELSE                                    ! Acoustic top half-space
        cMin = MIN( cMin, DBLE( HSTop%cP ) )
     ENDIF
  END SELECT

  ! If elastic medium then reduce cMin for Scholte wave
  IF ( ElasticFlag ) cMin = 0.85 * cMin
  cLow = MAX( cLow, 0.99 * cMin )

END SUBROUTINE INIT
!**********************************************************************C
SUBROUTINE SOLVE( FileRoot, Error )

  ! Solves the eigenvalue problem at the current mesh
  ! and produces a new extrapolation

  USE krakelmod
  IMPLICIT NONE
  INTEGER       :: j, Key
  REAL      (KIND=8)                :: T1, T2, F1, F2, TStart, TEnd, x1, x2
  REAL      (KIND=8), INTENT( OUT ) :: Error
  CHARACTER (LEN=80), INTENT(  IN ) :: FileRoot

  CALL CPU_TIME( Tstart )
  CALL INIT     ! set up the finite difference mesh                      
  CALL SOLVE2   ! solve for the eigenvalues

  Extrap( iSet, 1 : M ) = EVMat( iSet, 1 : M ) 
  IF ( iSet == 1 ) CALL VECTOR( FileRoot )   ! If this is the first mesh, compute the eigenvectors

  ! Richardson extrapolation to improve the accuracy          

  Error = HUGE( Error )
  KEY   = 2 * M / 3 + 1     ! index of element used to check convergence

  IF ( iSet > 1 ) THEN
     T1 = Extrap( 1, KEY )

     DO j = iSet, 2, -1
        DO Mode = 1, M
           F2 = Extrap( j,     Mode )
           F1 = Extrap( j - 1, Mode )
           x2 = NV( iSet  ) ** 2
           x1 = NV( j - 1 ) ** 2
           Extrap( j - 1, Mode ) = F2 - ( F1 - F2 ) / ( x2 / x1 - 1.0 )
        END DO
     END DO

     T2    = Extrap( 1, KEY )
     Error = ABS( T2 - T1 )
  ENDIF

  CALL CPU_TIME( Tend)   ! check elapsed time
  ET( iSet ) = Tend - Tstart
  WRITE( PRTFile, '( 1X, I8, 6X, G15.3 )' ) NV( iSet ), ET( iSet )

END SUBROUTINE SOLVE

!**********************************************************************C
SUBROUTINE SOLVE2

  ! Provides initial guess to root finder for each EVMat(I)

  USE krakelmod
  IMPLICIT NONE
  INTEGER            :: MaxIteration, NTotal, ii, j, IAllocStat, Iteration
  REAL      (KIND=8) :: Tolerance
  REAL      (KIND=8) :: x, P( 10 )
  CHARACTER (LEN=80) :: ErrorMessage

  x     = omega2 / cLow ** 2
  MaxIteration = 500

  ! solve1 has already allocated space for the following unless the problem has shear
  IF ( .NOT. ALLOCATED( k ) ) THEN
     M = 3000   ! this sets the upper limit on how many modes can be calculated
     ALLOCATE( EVMat( NSets, M ), Extrap( NSets, M ), k( M ), VG( M ), Stat = IAllocStat )
     IF ( IAllocStat /= 0 ) &
          CALL ERROUT( PRTFile, 'F', 'KRAKEN - SOLVE2', 'Insufficient memory (too many modes).' )
  END IF

  NTotal = SUM( N( NFact : NLact ) )

  DO Mode = 1, M

     ! For first or second meshes, use a high guess
     ! Otherwise use extrapolation to produce an initial guess

     x = 1.00001 * x
     IF ( TopOpt( 4 : 4 ) == 'R' ) x = k( Mode ) ** 2
     IF ( iSet >= 2 ) THEN

        DO ii = 1, iSet-1
           P( ii ) = EVMat( ii, Mode)
        END DO

        IF ( iSet >= 3 ) THEN
           DO ii = 1, iSet-2
              DO j = 1, iSet - ii - 1
                 P( j ) = ( ( hV( iSet ) ** 2 - hV( j + ii ) ** 2 ) * P( j ) - ( hV( iSet ) ** 2 &
                      -hV( j ) ** 2 ) * P( j+1 ) ) / ( hV( j ) ** 2 - hV( j + ii ) ** 2 )
              END DO
           END DO
           x = P( 1 )
        ENDIF

     ENDIF

     ! Use the secant method to refine the eigenvalue
     Tolerance = ABS( x ) * 10.0 ** ( 3.0 - PRECISION( x ) )
     CALL ZSecantX( x, Tolerance, Iteration, MaxIteration, ErrorMessage )

     IF ( ErrorMessage /= ' ' ) THEN
        WRITE( PRTFile, * ) ErrorMessage
        WRITE( PRTFile, * ) 'for iSet, Mode = ',iSet, Mode
     ENDIF

     EVMat( iSet, Mode ) = x

     ! Toss out modes outside user specified spectral limits
     IF ( SQRT( omega2 ) / cHigh > REAL( SQRT( x ) ) ) THEN
        M = Mode - 1
        RETURN
     ENDIF

     IF ( Mode > NTotal / 5 ) THEN
        WRITE( PRTFile, * ) 'Caution: mesh is too coarse to compute all modes in spectrum'
        M = Mode
        RETURN
     ENDIF
  END DO   ! next mode

  ! Check whether we have hit the ceiling on max # modes

  IF ( M == MaxM ) THEN
     WRITE( PRTFile, * ) 'Too many modes in phase speed interval--'
     WRITE( PRTFile, * ) 'program will compute as many as it can'
  ENDIF

END SUBROUTINE SOLVE2
!**********************************************************************C
SUBROUTINE FUNCT( x, Delta, IPOW )

  ! FUNCT( x ) = 0 is the dispersion relation

  USE krakelmod
  IMPLICIT NONE
  INTEGER,       INTENT( OUT ) :: IPOW
  INTEGER                      :: IPVT( MaxN ), NMat, lda, Ml, Mu, j, Info
  REAL (KIND=8)                :: A( 13, MaxN ), VDET( 2 )
  REAL (KIND=8), INTENT( IN  ) :: x
  REAL (KIND=8), INTENT( OUT ) :: Delta

  ! If we get to a k-value below the user's spectrum limit
  ! then force a zero by returning a new function value
  ! that is ten orders of magnitude less than the last value

  ! IF ( REAL( x ) <= omega2 / cHigh ** 2 ) THEN
  !    IPOW = IPOW - 10
  !    RETURN
  ! ENDIF

  CALL SETA( A, NMat, x ) ! Set up the matrix

  lda = 13
  ML  = 4
  MU  = 4

  CALL DGBTRF( NMat, NMat, ML, MU, A, lda, IPVT, INFO )  ! Factor the matrix

  IF ( INFO /= 0 ) THEN
     WRITE( PRTFile, * ) 'INFO = ', INFO
     STOP
  ENDIF

  CALL DGBDI( A, lda, NMat, ML, MU, IPVT, VDET )  ! Calculate the determinant

  Delta = VDET( 1 )
  IPOW  = INT( VDET( 2 ) )

  ! Deflate previous roots

  IF ( Mode > 1 ) THEN
     Delta = Delta / ( X - EVMat( iSet, Mode - 1 ) )
     IF ( Mode > 2 ) THEN
        DO j = 1, Mode-2
           Delta = ( EVMat( iSet, Mode - 1 ) - EVMat( iSet, j ) ) * Delta / &
                   ( X - EVMat( iSet, j ) )
        END DO
     ENDIF
  ENDIF

END SUBROUTINE FUNCT
!**********************************************************************C
SUBROUTINE VECTOR( FileRoot )

  ! Do inverse iteration to compute each of the eigenvectors
  ! and write these to the disk file

  USE krakelmod
  USE SdRdRMod
  IMPLICIT NONE
  INTEGER          :: IPVT( MaxN ), ii, Info, iz, j, Job, lda, ldb, Medium, ML, MU, NMat, NRHS, NTotal
  REAL             :: z( MaxN )
  REAL    (KIND=8) :: A( 13, MaxN ), EVector( MaxN ), x
  COMPLEX          :: EVectorS( MaxN )
  COMPLEX (KIND=8) :: cx
  CHARACTER (LEN=80), INTENT(  IN ) :: FileRoot

  ! Tabulate Z-coordinates & compute NMat, NTotal
  NMat   = 0   ! size of the matrix
  NTotal = 0   ! number of z-points
  j      = 0

  DO Medium = 1, NMedia
     j      = j + 1
     z( j ) = SNGL( Depth( Medium ) )

     DO ii = 1, N( Medium )
        j      = j + 1
        z( j ) = z( j - 1 ) + SNGL( h( Medium ) )
     END DO

     IF ( Material( Medium ) == 'ACOUSTIC' ) THEN
        NMat = NMat +       N( Medium ) + 1
     ELSE
        NMat = NMat + 4 * ( N( Medium ) + 1 )
     ENDIF

     NTotal = NTotal + N( Medium ) + 1
  END DO

  ! Open MODFile and write header

  LRECL = MAX( 2 * NMat, 32, 3 * NMedia ) ! Logical record length in `longwords' (4 bytes)
  OPEN ( FILE = TRIM( FileRoot ) //'.mod', UNIT = MODFile, ACCESS = 'DIRECT', RECL = 4 * LRECL, FORM = 'UNFORMATTED' )

  WRITE( MODFile, REC = 1 ) LRECL, TITLE, REAL( FREQ ),  NMedia, NTotal, NMat
  WRITE( MODFile, REC = 2 ) ( N( Medium ), Material( Medium ), Medium = 1, NMedia )
  WRITE( MODFile, REC = 3 ) HSTop%BC, CMPLX( HSTop%cP ), CMPLX( HSTop%cS ), REAL( HSTop%rho ), &
                            REAL( Depth( 1          ) ), &
                            HSBot%BC, CMPLX( HSBot%cP ), CMPLX( HSBot%cS ), REAL( HSBot%rho ), &
                            REAL( Depth( NMedia + 1 ) )
  WRITE( MODFile, REC = 4 ) ( REAL( Depth( Medium ) ), REAL( rho( Loc( Medium ) + 1 ) ), Medium = NFACT, NLACT )
  WRITE( MODFile, REC = 6 ) ( z( Iz ), Iz = 1, NTotal )

  ! MAIN LOOP: For each eigenvalue call SINVIT to get eigenvector

  DO Mode = 1, M
     x = 1.00000001 * EVMat( 1, Mode )

     CALL SETA( A, NMat, x ) ! Set up the matrix

     lda = 13
     ML  = 4
     MU  = 4
     CALL DGBTRF( NMat, NMat, ML, MU, A, lda, IPVT, INFO )  ! Factor the matrix

     IF ( INFO /= 0 ) THEN
        WRITE( PRTFile, * ) 'INFO = ', INFO
        STOP
     ENDIF

     JOB                 = 0
     EVector( 1 : NMat ) = 1.0
     NRHS                = 1      ! number of right hand sides
     ldb                 = NMat   ! first dimension of b
     CALL DGBTRS( 'N', NMat, ML, MU, NRHS, A, lda, IPVT, EVector, ldb, JOB )  ! Solve

     cx = x   ! convert to complex
     CALL Normalize( EVector, EVectorS, NMat, cx ) ! Normalize

     WRITE( MODFile, REC = 6 + Mode ) ( EVectorS( j ), j = 1, NMat )
  END DO   ! next mode

END SUBROUTINE VECTOR
!**********************************************************************C
SUBROUTINE SETA( A, NMat, x )

  ! ROUTINE FOR SETTING UP THE MATRIX A

  USE krakelmod
  IMPLICIT NONE
  INTEGER,          INTENT( OUT ) :: NMat
  INTEGER                         :: ii, j, l, Medium, IPow
  REAL    (KIND=8), INTENT( OUT ) :: A( 13, MaxN )
  REAL    (KIND=8)                :: gammaP, gammaS, gammaP2, gammaS2, rmu, factor, h_rho, two_by_h
  REAL    (KIND=8), INTENT( IN  ) :: x
  COMPLEX (KIND=8)                :: F, G

  ! Compute size of matrix
  NMat = 0
  DO Medium = 1, NMedia
     IF ( Material( Medium ) == 'ACOUSTIC' ) THEN
        NMat = NMat +       N( Medium ) + 1
     ELSE
        NMat = NMat + 4 * ( N( Medium ) + 1 )
     ENDIF
  END DO

  A = 0.0 ! Zero out the matrix

  ! MAIN LOOP ( over each Medium )

  L = 0   ! check this (L was uninitialized in previous version)
  j = 0
  DO Medium = 1, NMedia
     IF ( Medium /= 1 ) THEN

        ! Acousto-elastic interface

        IF ( Material( Medium - 1 ) == 'ACOUSTIC' .AND. Material( Medium ) == 'ELASTIC' ) THEN
           ! R2 Condition
           L = L + 1
           j = j + 1

           h_rho = h( Medium - 1 ) * rho( Loc( Medium - 1 ) + 1 )

           A( 10, j - 1 ) = 1.0
           A(  9, j     ) = 0.5 * ( B1( L ) - h( Medium - 1 ) ** 2 * x )
           A(  7, j + 2 ) = omega2 * h_rho

           ! R3 Condition
           A( 10, j     ) = 0.0
           A(  7, j + 3 ) = 1.0
           ! R4 Condition
           A( 11, j     ) = 1.0
           A(  7, j + 4 ) = 1.0
        ENDIF

        ! Elasto-acoustic interface

        IF ( Material( Medium - 1 ) == 'ELASTIC' .AND. Material( Medium ) == 'ACOUSTIC' ) THEN
           L = Loc( Medium ) + 1
           j = j + 5
           ! R3 Condition
           A(  9, j - 2 ) = 1.0
           ! R4 Condition
           A(  9, j - 1 ) = 1.0
           A(  8, j     ) = 1.0
           ! R2 Condition
           h_rho = h( Medium ) * rho( Loc( Medium ) + 1 )
           A( 12, j - 3 ) = -h_rho * omega2
           A(  9, j     ) = 0.5 * ( B1( L ) - h( Medium ) ** 2 * x )
           A(  8, j + 1 ) = 1.0
        ENDIF

        ! Acoustic/acoustic interface

        IF ( Material( Medium - 1 ) == 'ACOUSTIC' .AND. Material( Medium ) == 'ACOUSTIC' ) THEN
           ! Continuity of normal displacement
           L = L + 1
           j = j + 1
           h_rho = h( Medium - 1 ) * rho( Loc( Medium - 1 ) + 1 )
           A( 10, j - 1 ) = 1.0 / h_rho
           A(  9, j     ) = 0.5 * ( B1( L ) - h( Medium - 1 ) ** 2 * x ) / h_rho

           L = L + 1
           h_rho = h( Medium ) * rho( Loc( Medium ) + 1 )
           A(  8, j + 1 ) = 0.5 * ( B1( L ) - h( Medium     ) ** 2 * x ) / h_rho
           A(  7, j + 2 ) = 1.0 / h_rho
           ! Continuity of pressure
           j = j + 1
           A( 10, j - 1 ) =  1.0
           A(  9, j     ) = -1.0
        ENDIF

        ! Elastic/elastic interface

        IF ( Material( Medium - 1 ) == 'ELASTIC' .AND. Material(Medium) == 'ELASTIC' ) THEN
           ! Continuity of everything
           j = j + 1
           A( 11, j     ) =  1.0
           A(  7, j + 4 ) = -1.0
           j = j + 1
           A( 11, j     ) =  1.0
           A(  7, j + 4 ) = -1.0
           j = j + 1
           A( 11, j     ) =  1.0
           A(  7, j + 4 ) = -1.0
           j = j + 1
           A( 11, j     ) =  1.0
           A(  7, j + 4 ) = -1.0
        ENDIF
     ENDIF

     IF ( Material( Medium ) == 'ACOUSTIC' ) THEN

        ! Acoustic section

        L = Loc( Medium ) + 1
        IF ( Medium == 1 ) j = 1
        DO ii = 1, N( Medium ) - 1
           j = j + 1
           L = L + 1
           A( 10, j - 1 ) = 1.0
           A(  9, j     ) = B1( L ) - h( Medium ) * h( Medium ) * x
           A(  8, j + 1 ) = 1.0
        END DO
     ELSE

        ! Elastic section

        L = Loc( Medium )
        two_by_h = 2.0 / h( Medium )
        DO ii = 1, N( Medium )
           L = L + 1
           j = j + 1
           A( 11, j ) = two_by_h
           A( 12, j ) = x * B4( L )
           A( 13, j ) = x * B3( L ) - rho( L )
           j = j + 1
           A( 10, j ) = -1.0
           A( 11, j ) = two_by_h
           A( 13, j ) = -rho( L )
           j = j + 1
           A(  9, j ) = B1( L )
           A( 11, j ) = two_by_h
           A( 12, j ) = x
           j = j + 1
           A(  9, j ) = B2( L )
           A( 10, j ) = -B4( L )
           A( 11, j ) = two_by_h

           L = L + 1
           j = j + 1
           A( 7, j ) = -two_by_h
           A( 8, j ) = x * B4( L )
           A( 9, j ) = x * B3( L ) - rho( L )
           j = j + 1
           A( 6, j ) = -1.0
           A( 7, j ) = -two_by_h
           A( 9, j ) = -rho( L )
           j = j + 1
           A( 5, j ) = B1( L )
           A( 7, j ) = -two_by_h
           A( 8, j ) = x
           j = j + 1
           A( 5, j ) =  B2( L )
           A( 6, j ) = -B4( L )
           A( 7, j ) = -two_by_h
           j = j - 4
           L = L - 1
        END DO
     ENDIF
  END DO   ! next medium

  ! Top BC

  SELECT CASE ( Material( 1 ) )
  CASE ( 'ELASTIC' )      ! Elastic medium
     SELECT CASE ( HSTop%BC )
     CASE ( 'R' )         ! Rigid top
        A( 9, 1 ) = 1.0
        A( 9, 2 ) = 1.0
     CASE ( 'V' )         ! Free (vacuum) top
        A( 7, 3 ) = 1.0
        A( 7, 4 ) = 1.0
     END SELECT
  CASE( 'ACOUSTIC' )      ! Acoustic medium
     CALL BCImpedance( x, 'TOP', HSTop, F, G, IPOW )   ! Top impedance
     IF ( G == 0.0 ) THEN ! Free (vacuum) top
        A( 9, 1 ) = 1.0
        A( 8, 2 ) = 0.0
     ELSE
        h_rho     = rho( 1 ) * h( 1 )
        A( 9, 1 ) = 0.5 * ( B1( 1 ) - h( 1 ) ** 2 * x ) / h_rho + F / G
        A( 8, 2 ) = 1.0 / h_rho
     ENDIF
  END SELECT

  ! Bottom BC

  SELECT CASE ( Material( NMedia ) )
  CASE ( 'ELASTIC' )     ! Elastic medium
     SELECT CASE ( HSBot%BC )
     CASE ( 'R' )        ! Rigid bottom
        A( 11, NMat - 3 ) = 1.0
        A( 11, NMat - 2 ) = 1.0
     CASE ( 'V' )        ! Free (vacuum) bottom
        A(  9, NMat - 1 ) = 1.0
        A(  9, NMat     ) = 1.0
     CASE ( 'A' )        ! Elastic bottom
        gammaP2 = x - omega2 / HSBot%cP ** 2
        gammaS2 = x - omega2 / HSBot%cS ** 2
        gammaP  = SQRT( gammaP2 )
        gammaS  = SQRT( gammaS2 )
        Rmu     = HSBot%rho * HSBot%cS ** 2
        factor  = -Rmu / ( gammaS * gammaP - x )

        A( 11, NMat - 3 ) = factor * gammaP * ( x - gammaS2 )
        A( 10, NMat - 2 ) = factor * ( 2.0 * gammaS * gammaP - ( gammaS2 + x ) )
        A( 12, NMat - 3 ) = x * A( 10, NMat - 2 )
        A( 11, NMat - 2 ) = factor * gammaS * ( x - gammaS2 )
        A(  9, NMat - 1 ) = 1.0
        A(  9, NMat     ) = 1.0
     END SELECT
  CASE( 'ACOUSTIC' )    ! Acoustic medium
     CALL BCImpedance( x, 'BOT', HSBot, F, G, IPOW )   ! Bottom impedance

     IF ( G == 0.0 ) THEN ! Free (vacuum) bottom
        A( 10, NMat - 1 ) = 0.0
        A(  9, NMat     ) = 1.0
     ELSE
        L = L + 1
        h_rho             = h( NMedia ) * rho( Loc( NMedia ) + 1 )
        A( 10, NMat - 1 ) = 1.0 / h_rho
        A(  9, NMat     ) = 0.5 * ( B1( L ) - h( NMedia ) ** 2 * x ) / h_rho - F / G
     ENDIF
  END SELECT

END SUBROUTINE SETA
!**********************************************************************C
SUBROUTINE BCImpedance( x, BotTop, HS, F, G, IPow )                                                   

  ! Compute Boundary Condition IMPedance
  ! Same as KRAKENC except for ELASUP, ELASDN removed

  USE krakelmod
  use RefCoMod
  IMPLICIT NONE
  INTEGER, INTENT( OUT ) :: IPow
  INTEGER                :: IBot, ITop, iPower
  REAL    (KIND=8)       :: x, omega, RadDeg, rhoInside
  COMPLEX (KIND=8)       :: cx, kx, kz, Twersky, gammaP, RCmplx, CInside
  COMPLEX (KIND=8), INTENT( OUT ) :: F, G
  CHARACTER      (LEN=3) :: BotTop
  TYPE( HSInfo )         :: HS
  TYPE(ReflectionCoef)   :: RInt

  cx   = x   ! complex version of x for routines that expect a complex argument
  IPow = 0 

  ! Get rho, c just Inside the boundary               
  ! There is at least one acoustic layer in the problem, except
  ! in the case where BOUNCE is used to get the refl. coef. for       
  ! a purely elastic stack of layers.                                 

  SELECT CASE ( BotTop )
  CASE ( 'TOP' ) 
     IF ( NFact > 0 ) THEN
        Itop       = Loc( NFAct ) + N( NFAct ) + 1
        rhoInside  = rho( Itop )
        cInside    = SQRT( omega2 * h( NFAct ) ** 2 / ( 2.0 + B1( Itop ) ) )
     ENDIF
  CASE ( 'BOT' )
     IF ( NLACT > 0 ) THEN 
        Ibot       = Loc( NLACT ) + N( NLACT ) + 1 
        rhoInside  = rho( Ibot )
        cInside    = SQRT( omega2 * h( NLACT ) ** 2 / ( 2.0 + B1( Ibot ) ) )
     ENDIF
  END SELECT

  SELECT CASE ( HS%BC )
  CASE ( 'V' )                   ! Vacuum
     F = 1.0
     ! G = -CI * PEKRT( omega2 / cInside ** 2 - x ) * sigma(1) ** 2
     G = 0.0
  CASE (  'S', 'H', 'T', 'I'  )  ! Vacuum with Twersky scatter model
     omega = SQRT( omega2 )
     kx = SQRT( x )
     F = 1.0
     G = Twersky( omega, HS, kx, rhoInside, cInside )
     G = G / ( i * omega * rhoInside )
  CASE ( 'R' )                    ! Rigid
     F = 0.0
     G = 1.0
  CASE ( 'A' )                    ! Acousto-elastic half-space
     gammaP = SQRT( x - omega2 / HS%cP ** 2 )
     F = 1.0
     G = HS%rho / gammaP
  CASE ( 'F' )                    ! Tabulated reflection coefficient
     ! Compute the grazing angle THETA
     kx         = SQRT( x ) 
     kz         = SQRT( omega2 / CInside ** 2 - x ) 
     RadDeg     = 180.0 / pi
     RInt%theta = RadDeg * DATAN2( DBLE( kz ), DBLE( kx ) )

     ! Evaluate R( ThetaInt )
     IF ( BotTop == 'TOP' ) THEN
        CALL InterpolateReflectionCoefficient( RInt, RTop, NTopPts, PRTFile )
     ELSE
        CALL InterpolateReflectionCoefficient( RInt, RBot, NBotPts, PRTFile )
     ENDIF

     ! Convert R(THETA) to (f,g) in Robin BC
     RCmplx = RInt%R * EXP( i * RInt%phi )
     F      = 1.0
     G      = ( 1.0 + RCmplx ) / ( i * kz * ( 1.0 - RCmplx ) )
  CASE ( 'P' )                    ! Precalculated reflection coef
     CALL InterpolateIRC( cx, F, G, IPower, xTab, FTab, GTab, ITab, NkTab )
  END SELECT

  IF ( BotTop == 'TOP' ) G = -G    ! A top BC has the sign flipped relative to a bottom BC

END SUBROUTINE BCImpedance
!**********************************************************************C
SUBROUTINE Normalize( EVector, EVectorS, NMat, x )

  ! Convert the EVector to ( u, w, tau_zx, tau_zz ) and normalize

  USE krakelmod
  IMPLICIT NONE
  INTEGER          :: ii, Medium, NMat, j, kk
  REAL    (KIND=8) :: EVector( MaxN ), SqNorm2, rhoMedium
  COMPLEX          :: EVectorS( MaxN ), CIk
  COMPLEX (KIND=8) :: R1, R2, R3, R4, x, SqNorm, Contrib, RN, ScaleFac, gammaP

  CIk = i * SQRT( x )

  ! Scale down the mode

  SqNorm2 = 0.0

  DO ii = 1, NMat
     SqNorm2 = MAX( SqNorm2, ABS( EVector( ii ) ) )
  END DO

  EVector( 1 : NMat ) = EVector( 1 : NMat ) / SqNorm2

  ! Loop to compute norm of eigenvector

  SqNorm = 0.0
  j      = 1
  kk     = 1

  DO Medium = 1, NMedia
     rhoMedium = rho( Loc( Medium ) + 1 )

     DO ii = 1, N( Medium ) + 1
        IF ( Material( Medium ) == 'ELASTIC' ) THEN
           R1 = EVector( kk     )
           R2 = EVector( kk + 1 )
           R3 = EVector( kk + 2 )
           R4 = EVector( kk + 3 )

           Contrib = h( Medium ) * ( -x * B3( j ) * R1 ** 2 + B4( j ) * R1 * R4 - R2 * R3 )

           EVectorS( kk     ) = CIk * EVector( kk     )
           EVectorS( kk + 1 ) =       EVector( kk + 1 )
           EVectorS( kk + 2 ) = CIk * EVector( kk + 2 )
           EVectorS( kk + 3 ) =       EVector( kk + 3 )

           kk = kk + 4
        ELSE
           Contrib        = -h( Medium ) * EVector( kk ) ** 2 / ( rhoMedium * omega2 )
           EVectorS( kk ) = -EVector( kk )   ! sign flip implies this is tau_zz not pressure
           kk             = kk + 1
        ENDIF

        IF ( ii == 1 .OR. ii == N( Medium ) + 1 ) Contrib = 0.5 * Contrib
        SqNorm = SqNorm + Contrib
        j      = j + 1
     END DO
  END DO

  ! Bottom half-space contribution

  IF ( HSBot%BC == 'A' ) THEN
     IF ( Material( NMedia + 1 ) == 'ELASTIC' ) THEN
        WRITE( PRTFile, * ) 'Elastic halfspace normalization not implemented'
     ELSE
        gammaP  = SQRT( x - omega2 / HSBot%cP ** 2 )
        Contrib = -EVector( NMat ) ** 2 / ( 2.0 * gammaP * HSBot%rho * omega2 )
     ENDIF
  ENDIF

  SqNorm = SqNorm + Contrib

  ! Scale the mode

  RN       = -omega2 * SqNorm
  ScaleFac = 1.0 / SQRT( RN )
  EVectorS( 1 : NMat ) = ScaleFac * EVectorS( 1 : NMat )

END SUBROUTINE Normalize
