SUBROUTINE InfluenceR( U, epsilon, alpha, iBeamWindow2, RadMax )

  ! Computes the beam influence, i.e. the contribution of a single beam to the complex pressure
  ! This routine is for the Cerveny-style beams
  ! The 'R' in InfluenceR is for beams in Ray-centered coordinates

  USE bellhopMod
  USE SdRdRMod
  IMPLICIT NONE
  INTEGER,       INTENT( IN    ) :: IBeamWindow2
  REAL (KIND=8), INTENT( IN    ) :: alpha, RadMax                ! take-off angle
  COMPLEX,       INTENT( INOUT ) :: U( Nrd_per_range, Pos%Nr )   ! complex pressure field
  INTEGER          :: is, id, ir, ir1, ir2, KMAHV( MaxN ), KMAH, image
  REAL    (KIND=8) :: nA, nB, nSq, rA, rB, c, W, zr, znv, rnv, ratio1, n, Hermite
  COMPLEX (KIND=8) :: pVB( MaxN ), qVB( MaxN ), q, epsV( MaxN ), contri, gammaV( MaxN ), gamma, P_n, P_s, epsilon
  COMPLEX (KIND=8) :: tau

  ! need to add logic related to Nrd_per_range

  ! During reflection imag(q) is constant and adjacent normals cannot bracket a segment of the TL
  ! line, so no special treatment is necessary

  IF ( Beam%Type( 2 : 2 ) == 'C' ) THEN
     epsV( 1 : Beam%Nsteps ) = i * ABS( ray2D( 1 : Beam%Nsteps )%q( 1 ) / ray2D( 1 : Beam%Nsteps )%q( 2 ) )
  ELSE
     epsV( 1 : Beam%Nsteps ) = epsilon
  END IF

  pVB(    1 : Beam%Nsteps ) = ray2D( 1 : Beam%Nsteps )%p( 1 ) + epsV( 1 : Beam%Nsteps ) * ray2D( 1 : Beam%Nsteps )%p( 2 )
  qVB(    1 : Beam%Nsteps ) = ray2D( 1 : Beam%Nsteps )%q( 1 ) + epsV( 1 : Beam%Nsteps ) * ray2D( 1 : Beam%Nsteps )%q( 2 )
  gammaV( 1 : Beam%Nsteps ) = pVB(   1 : Beam%Nsteps ) / qVB( 1 : Beam%Nsteps )

  IF ( Beam%RunType( 4 : 4 ) == 'R' ) THEN
     Ratio1 = SQRT( ABS( COS( alpha ) ) )  ! point source
  ELSE
     Ratio1 = 1                            ! line  source
  END IF

  ! compute KMAH index
  ! Following is incorrect for 'Cerveny'-style beamwidth (narrow as possible)
  KMAHV(  1 ) = 1

  DO is = 2, Beam%Nsteps
     KMAHV(  is ) = KMAHV( is - 1 )
     CALL BranchCut( qVB( is - 1 ), qVB( is ), Beam%Type, KMAHV( is ) )
  END DO

  RcvrDepths: DO id = 1, Nrd_per_range
     zR = Pos%rd( id )

     Images: DO image = 1, Beam%Nimage
        ir1 = HUGE( ir1 )

        Stepping: DO is = 2, Beam%Nsteps

           ! Compute ray-centered coordinates, (znV, rnV)

           znV = -ray2D( is )%t( 1 ) * ray2D( is )%c
           rnV =  ray2D( is )%t( 2 ) * ray2D( is )%c

           IF ( ABS( znV ) < tiny( znV ) ) THEN   ! Check for normal parallel to TL-line
              CYCLE Stepping   ! skip to next step on ray
           END IF

           SELECT CASE ( image )     ! Images of beams
           CASE ( 1 )                ! True beam
              nB  = ( zR -                             ray2D( is )%x( 2 )   ) / znV
           CASE ( 2 )                ! Surface-reflected beam
              rnV = -rnV
              nB  = ( zR - ( 2.0 * Bdry%Top%HS%Depth - ray2D( is )%x( 2 ) ) ) / znV
           CASE ( 3 )                ! Bottom-reflected beam
              rnV = -rnV
              nB  = ( zR - ( 2.0 * Bdry%Bot%HS%Depth - ray2D( is )%x( 2 ) ) ) / znV
           END SELECT

           rB  = ray2D( is )%x( 1 ) + nB * rnV
           ir2 = MAX( MIN( INT( ( rB - Pos%r( 1 ) ) / Pos%Delta_r ) + 1, Pos%Nr ), 1 ) ! index of receiver

           ! detect and skip duplicate points (happens at boundary reflection)
           IF ( ir1 >= ir2 .OR. ABS( ray2D( is )%x( 1 ) - ray2D( is - 1 )%x( 1 ) ) < 1.0D3 * SPACING( ray2D( is )%x( 1 ) ) ) THEN
              rA  = rB
              nA  = nB
              ir1 = ir2
              CYCLE Stepping
           END IF

           Ranges: DO ir = ir1 + 1, ir2    ! Compute influence for each rcvr
              W     = ( Pos%r( ir ) - rA ) / ( rB - rA )
              q     =    qVB( is - 1 ) + W * (    qVB( is ) -    qVB( is - 1 ) )
              gamma = gammaV( is - 1 ) + W * ( gammaV( is ) - gammaV( is - 1 ) )
              n     = nA + W * ( nB - nA )
              nSq   = n * n
              IF ( AIMAG( gamma ) > 0 ) THEN
                 WRITE( PRTFile, * ) 'Unbounded beam'
                 CYCLE Ranges   ! next receiver range
              END IF

              IF ( -0.5 * omega * AIMAG( gamma ) * nSq < iBeamWindow2 ) THEN   ! Within beam window?
                 c      = ray2D( is - 1 )%c
                 tau    = ray2D( is - 1 )%tau + W * ( ray2D( is )%tau - ray2D( is - 1 )%tau )
                 contri = ratio1 * ray2D( is )%Amp * SQRT( c * ABS( epsV( is ) ) / q ) * &
                      EXP( -i * ( omega * ( tau + 0.5 * gamma * nSq ) - ray2D( is )%phase ) )

                 SELECT CASE ( Beam%Component )
                 CASE ( 'P' )   ! pressure
                 CASE ( 'V' )   ! vertical component
                    P_n    = -i * omega * gamma * n * contri
                    P_s    = -i * omega / c         * contri
                    contri = c * DOT_PRODUCT( [ P_n, P_s ], ray2D( is )%t ) 
                 CASE ( 'H' )   ! horizontal component
                    P_n    = -i * omega * gamma * n * contri
                    P_s    = -i * omega / c         * contri
                    contri = c * ( -P_n * ray2D( is )%t( 2 ) + P_s * ray2D( is )%t( 1 ) ) 
                 END SELECT

                 KMAH = KMAHV( is - 1 )
                 CALL BranchCut( qVB( is - 1 ), q, Beam%Type, KMAH ) ! Get correct branch of SQRT

                 IF ( KMAH  < 0  ) contri = -contri
                 IF ( image == 2 ) contri = -contri

                 SELECT CASE ( Beam%RunType( 1 : 1 ) )
                 CASE ( 'I', 'S' )   ! Incoherent or Semi-coherent TL
                    contri = ABS(contri)
                 END SELECT

                 U( id, ir ) = U( id, ir ) + CMPLX( Hermite( n, RadMax, 2 * RadMax ) * contri )
              END IF
           END DO Ranges
           rA  = rB
           nA  = nB
           ir1 = ir2
        END DO Stepping
     END DO Images
  END DO RcvrDepths

END SUBROUTINE InfluenceR

! **********************************************************************!

SUBROUTINE InfluenceC( U, epsilon, alpha, iBeamWindow2, RadMax )

  ! Computes the beam influence, i.e. the contribution of a single beam to the complex pressure
  ! This version uses a beam representation in Cartesian coordinates

  USE bellhopMod
  USE SdRdRMod
  USE sspMod
  IMPLICIT NONE
  INTEGER,       INTENT( IN    ) :: IBeamWindow2
  REAL (KIND=8), INTENT( IN    ) :: alpha, RadMax                ! take-off angle
  COMPLEX,       INTENT( INOUT ) :: U( Nrd_per_range, Pos%Nr )   ! complex pressure field
  INTEGER          :: KMAHV( MaxN ), KMAH, is, id, ir, ir1, ir2
  REAL    (KIND=8) :: x( 2 ), rayt( 2 ), rayn( 2 ), Tr, Tz, rA, rB, zr, W, ratio1, &
                      c, cimag, cs, cn, csq, gradc( 2 ), crr, crz, czz, Hermite, deltaz
  COMPLEX (KIND=8) :: pVB( MaxN ), qVB( MaxN ), q, epsV( MaxN ), contri, gammaV( MaxN ), gamma, const, epsilon
  COMPLEX (KIND=8) :: tau

  ! need to add logic related to Nrd_per_range

  ! During reflection imag(q) is constant and adjacent normals cannot bracket a segment of the TL
  ! line, so no special treatment is necessary
  
  IF ( Beam%Type( 2 : 2 ) == 'C' ) THEN
     epsV( 1 : Beam%Nsteps ) = i * ABS( ray2D( 1 : Beam%Nsteps )%q( 1 ) / ray2D( 1 : Beam%Nsteps )%q( 2 ) )
  ELSE
     epsV( 1 : Beam%Nsteps ) = epsilon
  END IF

  pVB( 1 : Beam%Nsteps ) = ray2D( 1 : Beam%Nsteps )%p( 1 ) + epsV( 1 : Beam%Nsteps ) * ray2D( 1 : Beam%Nsteps )%p( 2 )
  qVB( 1 : Beam%Nsteps ) = ray2D( 1 : Beam%Nsteps )%q( 1 ) + epsV( 1 : Beam%Nsteps ) * ray2D( 1 : Beam%Nsteps )%q( 2 )

  IF ( Beam%RunType( 4 : 4 ) == 'R' ) THEN
     Ratio1 = SQRT( ABS( COS( alpha ) ) )  ! point source
  ELSE
     Ratio1 = 1                            ! line source
  END IF

  ! Form gamma and KMAH index
  ! Treatment of KMAH index is incorrect for 'Cerveny' style beam width BeamType

  Stepping0: DO is = 1, Beam%Nsteps

     rayt = ray2D( is )%c * ray2D( is )%t   ! unit tangent
     rayn = [ rayt( 2 ), -rayt( 1 ) ]       ! unit normal

     CALL EvaluateSSP( ray2D( is )%x, c, cimag, gradc, crr, crz, czz, Freq, 'TAB' )

     csq = c * c
     cS  = DOT_PRODUCT( gradc, rayt )
     cN  = DOT_PRODUCT( gradc, rayn )

     Tr  = rayt(  1 )
     Tz  = rayt(  2 )

     gammaV( is ) = 0.0
     IF ( qVB( is ) /= 0.0 ) gammaV( is ) = 0.5 * ( pVB( is ) / qVB( is ) * Tr**2 + 2.0 * cN / csq * Tz * Tr - cS / csq * Tz**2 )

     IF ( is == 1 ) THEN
        KMAHV( 1 ) = 1
     ELSE
        KMAHV( is ) = KMAHV( is - 1 )
        CALL BranchCut( qVB( is - 1 ), qVB( is ), Beam%Type, KMAHV( is ) )
     END IF

     !RLTEMP = SQRT( -2.0 / ( omega * AIMAG( pVB( is ) / qVB( is ) ) ) )
     !RKTEMP = -cV( is ) * REAL( pVB( is ) / qVB( is ) )
     !write( *, * ) is, rltemp, rktemp

  END DO Stepping0

  Stepping: DO is = 3, Beam%Nsteps
     IF ( ray2D( is     )%x( 1 ) > Pos%r( Pos%Nr ) ) RETURN
     rA = ray2D( is - 1 )%x( 1 )
     rB = ray2D( is     )%x( 1 )
     IF ( ABS( rB - rA ) < 1.0D3 * SPACING( rB ) ) CYCLE Stepping   ! don't process duplicate points

     ! Compute upper index on rcvr line
     ! Assumes r is a vector of equally spaced points
     ir1 = MAX( MIN( INT( ( rA - Pos%r( 1 ) ) / Pos%Delta_r ) + 1, Pos%Nr ), 1 ) ! should be ", 0 )" ?
     ir2 = MAX( MIN( INT( ( rB - Pos%r( 1 ) ) / Pos%Delta_r ) + 1, Pos%Nr ), 1 )

     IF ( ir1 >= ir2 ) CYCLE Stepping

     Ranges: DO ir = ir1 + 1, ir2

        W     = ( Pos%r( ir ) - rA ) / ( rB - rA )

        x     = ray2D(    is - 1 )%x    + W * ( ray2D(       is )%x    -  ray2D(    is - 1 )%x   )
        rayt  = ray2D(    is - 1 )%t    + W * ( ray2D(       is )%t    -  ray2D(    is - 1 )%t   )
        c     = ray2D(    is - 1 )%c    + W * ( ray2D(       is )%c    -  ray2D(    is - 1 )%c   )
        q     = qVB(      is - 1 )      + W * ( qVB(         is )      -  qVB(      is - 1 )     )
        tau   = ray2D(    is - 1 )%tau  + W * ( ray2D(       is )%tau  -  ray2D(    is - 1 )%tau )
        gamma = gammaV(   is - 1 )      + W * ( gammaV(      is )      -  gammaV(   is - 1 )     )

        IF ( AIMAG( gamma ) > 0 ) THEN
           WRITE( PRTFile, * ) 'Unbounded beam'
           WRITE( PRTFile, * ) gammaV( is - 1 ), gammaV( is ), gamma
           CYCLE Ranges   ! next receiver range
        END IF

        const = ratio1 * SQRT( c * ABS( epsV( is - 1 ) ) / q )

        ! Get correct branch of SQRT
        KMAH = KMAHV( is - 1 )
        CALL BranchCut( qVB( is - 1 ), q, Beam%Type, KMAH )
        IF ( KMAH < 0 ) const = -const

        RcvrDepths: DO id = 1, Nrd_per_range
           zR = Pos%rd( id )
                                   ! True beam
              deltaz = zR - x( 2 )
              IF ( omega * AIMAG( gamma ) * deltaz**2 < iBeamWindow2 ) &
                   contri =           ray2D( is )%Amp * Hermite( deltaz, RadMax, 2.0 * RadMax ) * &
                   EXP( -i * ( omega * ( tau + rayt( 2 ) * deltaz + gamma * deltaz**2 ) - ray2D( is )%Phase ) )

           IF ( Beam%Nimage >= 2 ) THEN ! Surface reflected beam
              deltaz = -zR + 2.0 * Bdry%Top%HS%Depth - x( 2 )
              IF ( omega * AIMAG( gamma ) * deltaz**2 < iBeamWindow2 ) &
                   contri =  contri - ray2D( is )%Amp * Hermite( deltaz, RadMax, 2.0 * RadMax ) * &
                   EXP( -i * ( omega * ( tau + rayt( 2 ) * deltaz + gamma * deltaz**2 ) - ray2D( is )%Phase ) )
           END IF

           IF ( Beam%Nimage >= 3 ) THEN ! Bottom reflected beam
              deltaz = -zR + 2.0 * Bdry%Bot%HS%Depth - x( 2 )
              IF ( omega * AIMAG( gamma ) * deltaz**2 < iBeamWindow2 ) &
                   contri =  contri + ray2D( is )%Amp * Hermite( deltaz, RadMax, 2.0 * RadMax ) * &
                   EXP( -i * ( omega * ( tau + rayt( 2 ) * deltaz + gamma * deltaz**2 ) - ray2D( is )%Phase ) )
           END IF

           ! contribution to field
           SELECT CASE( Beam%RunType( 1 : 1 ) )
           CASE ( 'C' )        ! coherent
              contri = const * contri
           CASE ( 'I', 'S' )   ! incoherent or semi-coherent
              contri = ABS( const * contri )
           END SELECT
           U( id, ir ) = U( id, ir ) + CMPLX( contri )
        END DO RcvrDepths
     END DO Ranges
  END DO Stepping

END SUBROUTINE InfluenceC

! **********************************************************************!

SUBROUTINE InfluenceGeoHat( U, alpha, Dalpha )

  ! Computes the beam influence, i.e. the contribution of a single beam to the complex pressure
  ! This version uses geometric, hat-shaped beams

  USE bellhopMod
  USE SdRdRMod
  USE ArrMod

  IMPLICIT NONE
  REAL (KIND=8), INTENT( IN    ) :: alpha, dalpha                ! take-off angle, angular spacing
  COMPLEX,       INTENT( INOUT ) :: U( Nrd_per_range, Pos%Nr )   ! complex pressure field
  INTEGER          :: is, id, ir, irT( 1 ), irTT
  REAL    (KIND=8) :: x_ray( 2 ), rayt( 2 ), rayn( 2 ), x_rcvr( 2 ), s, n, phaseInt, &
                      rA, rB, rLen, RadMax, zMin, zMax, q0, qold, A, const, ratio1, W, &
                      Amp, phase, dqds, q, SrcAngle, RcvrAngle
  COMPLEX (KIND=8) :: delay, dtauds

  q0       = ray2D( 1 )%c / Dalpha   ! Reference for J = q0 / q
  SrcAngle = RadDeg * alpha          ! take-off angle in degrees
  phase    = 0.0
  qOld     = ray2D( 1 )%q( 1 )       ! used to track KMAH index
  rA       = ray2D( 1 )%x( 1 )       ! range at start of ray

  ! what if never satistified?
  ! what if there is a single receiver (ir = 0 possible)
  irT = MINLOC( Pos%r( 1 : Pos%Nr ), MASK = Pos%r( 1 : Pos%Nr ) > rA )   ! find index of first receiver to the right of rA
  ir  = irT( 1 )
  IF ( ray2D( 1 )%t( 1 ) < 0.0d0 ) ir = ir - 1  ! if ray is left-traveling, then we want the first receiver to the left of rA

  IF ( Beam%RunType( 4 : 4 ) == 'R' ) THEN
     Ratio1 = SQRT( ABS( COS( alpha ) ) )  ! point source
  ELSE
     Ratio1 = 1                            ! line  source
  END IF

  Stepping: DO is = 2, Beam%Nsteps
     rB     = ray2D( is     )%x( 1 )
     x_ray  = ray2D( is - 1 )%x

     ! compute normalized tangent (compute it because we need to measure the step length)
     rayt = ray2D( is )%x - ray2D( is - 1 )%x
     rlen = NORM2( rayt )
     IF ( rlen < 1.0D3 * SPACING( ray2D( is )%x( 1 ) ) ) CYCLE Stepping  ! if duplicate point in ray, skip to next step along the ray
     rayt = rayt / rlen                    ! unit tangent to ray
     rayn = [ -rayt( 2 ), rayt( 1 ) ]      ! unit normal  to ray

     dqds   = ray2D( is )%q( 1 ) - ray2D( is - 1 )%q( 1 )
     dtauds = ray2D( is )%tau    - ray2D( is - 1 )%tau

     ! phase shifts at caustics
     q  = ray2D( is - 1 )%q( 1 )
     IF ( q <= 0.0d0 .AND. qOld > 0.0d0 .OR. &
          q >= 0.0d0 .AND. qOld < 0.0d0 ) phase = phase + pi / 2.  ! phase shifts at caustics
     qold = q

     RadMax = MAX( ABS( ray2D( is - 1 )%q( 1 ) ), ABS( ray2D( is )%q( 1 ) ) ) / q0 / ABS( rayt( 1 ) ) ! beam radius projected onto vertical line
     zmin   = MIN( ray2D( is - 1 )%x( 2 ), ray2D( is )%x( 2 ) ) - RadMax  ! min depth of ray segment
     zmax   = MAX( ray2D( is - 1 )%x( 2 ), ray2D( is )%x( 2 ) ) + RadMax  ! max depth of ray segment

     ! is this a steep ray? Then don't try to get depth limits: it's too complicated
     IF ( ABS( rayt( 1 ) ) < 0.5 ) THEN
        zmin = -HUGE( zmin )
        zmax = +HUGE( zmax )
     END IF

     ! compute beam influence for this segment of the ray
     RcvrRanges: DO
        ! is r( ir ) contained in [ rA, rB ]? Then compute beam influence
        IF ( Pos%r( ir ) >= MIN( rA, rB ) .AND. Pos%r( ir ) < MAX( rA, rB ) ) THEN
           RcvrDepths: DO id = 1, Nrd_per_range
              IF ( Beam%RunType( 5 : 5 ) == 'I' ) THEN
                 IF ( Pos%rd( ir ) < zmin .OR. Pos%rd( ir ) > zmax ) CYCLE RcvrDepths
                 x_rcvr = [ Pos%r( ir ), Pos%rd( ir ) ]   ! irregular   grid
              ELSE
                 IF ( Pos%rd( id ) < zmin .OR. Pos%rd( id ) > zmax ) CYCLE RcvrDepths
                 x_rcvr = [ Pos%r( ir ), Pos%rd( id ) ]   ! rectilinear grid
              END IF

              s      =      DOT_PRODUCT( x_rcvr - x_ray, rayt ) / rlen ! proportional distance along ray
              n      = ABS( DOT_PRODUCT( x_rcvr - x_ray, rayn ) )      ! normal distance to ray
              q      = ray2D( is - 1 )%q( 1 ) + s * dqds               ! interpolated amplitude
              RadMax = ABS( q / q0 )                                   ! beam radius

              IF ( n < RadMax ) THEN
                 A     = 1 / RadMax
                 delay = ray2D( is - 1 )%tau + s * dtauds           ! interpolated delay
                 const = Ratio1 * SQRT( ray2D( is )%c / ABS( q ) ) * A * ray2D( is )%Amp
                 Amp   = const * ( RadMax - n )

                 ! phase shifts at caustics
                 phaseInt = phase
                 IF ( q <= 0.0d0 .AND. qOld > 0.0d0 .OR. &
                      q >= 0.0d0 .AND. qOld < 0.0d0 ) phaseInt = phase + pi / 2.

                 SELECT CASE( Beam%RunType( 1 : 1 ) )
                 CASE ( 'E' )      ! eigenrays
                    CALL WriteRay( SrcAngle, is )
                 CASE ( 'A', 'a' ) ! arrivals
                    RcvrAngle  = RadDeg * ATAN2( rayt( 2 ), rayt( 1 ) )
                    CALL AddArr( omega, id, ir, Amp, ray2D( is - 1 )%Phase + phaseInt, delay, &
                         SrcAngle, RcvrAngle, ray2D( is )%NumTopBnc, ray2D( is )%NumBotBnc )
                 CASE ( 'C'  )     ! coherent TL
                    U( id, ir ) = U( id, ir ) + CMPLX( Amp * EXP( -i * ( omega * delay - ray2D( is - 1 )%Phase - phaseInt ) ) )
                 CASE DEFAULT      ! incoherent/semi-coherent TL
                    W           = ( RadMax - n ) / RadMax   ! hat function: 1 on center, 0 on edge
                    U( id, ir ) = U( id, ir ) + SNGL( ( Amp / W ) ** 2 * W )
                 END SELECT
              END IF
           END DO RcvrDepths
        END IF

        ! bump receiver index, ir, towards rB
        IF ( Pos%r( ir ) < rB ) THEN
           IF ( ir >= Pos%Nr        ) EXIT   ! go to next step on ray
           irTT = ir + 1                     ! bump right
           IF ( Pos%r( irTT ) >= rB ) EXIT
        ELSE
           IF ( ir <= 1             ) EXIT   ! go to next step on ray
           irTT = ir - 1                     ! bump left
           IF ( Pos%r( irTT ) <= rB ) EXIT
        END IF
        ir = irTT
     END DO RcvrRanges

     rA = rB
  END DO Stepping

END SUBROUTINE InfluenceGeoHat

! **********************************************************************!

SUBROUTINE InfluenceGeoGaussian( U, alpha, Dalpha )

  ! Computes the beam influence, i.e. the contribution of a single beam to the complex pressure
  ! This version uses geometric, Gaussian beams

  USE bellhopMod
  USE SdRdRMod
  USE ArrMod

  IMPLICIT NONE
  INTEGER,       PARAMETER       :: BeamWindow = 4               ! beam window: kills beams outside e**(-0.5 * ibwin**2 )
  REAL (KIND=8), INTENT( IN    ) :: alpha, dalpha                ! take-off angle, angular spacing
  COMPLEX,       INTENT( INOUT ) :: U( Nrd_per_range, Pos%Nr )   ! complex pressure field
  INTEGER          :: is, id, ir, irT( 1 ), irTT
  REAL    (KIND=8) :: x_ray( 2 ), rayt( 2 ), rayn( 2 ), &
                      x_rcvr( Nrd_per_range, 2 ), s( Nrd_per_range ), n( Nrd_per_range ), phaseInt, &
                      rA, rB, rLen, RadMax( Nrd_per_range ), zMin, zMax, &
                      q( Nrd_per_range ), q0, qold, qnew, sigma( Nrd_per_range ), lambda, &
                      A( Nrd_per_range ), const( Nrd_per_range ), ratio1, W( Nrd_per_range ), &
                      Amp( Nrd_per_range ), phase, dqds, SrcAngle, RcvrAngle
  COMPLEX (KIND=8) :: delay( Nrd_per_range ), dtauds

  q0       = ray2D( 1 )%c / Dalpha   ! Reference for J = q0 / q
  SrcAngle = RadDeg * alpha          ! take-off angle in degrees
  phase    = 0
  qOld     = ray2D( 1 )%q( 1 )       ! used to track KMAH index
  rA       = ray2D( 1 )%x( 1 )       ! range at start of ray

  ! what if never satistified?
  ! what if there is a single receiver (ir = 0 possible)
  irT       = MINLOC( Pos%r( 1 : Pos%Nr ), MASK = Pos%r( 1 : Pos%Nr ) > rA )      ! find index of first receiver to the right of rA
  ir        = irT( 1 )

  IF ( ray2D( 1 )%t( 1 ) < 0.0d0 ) ir = ir - 1  ! if ray is left-traveling, get the first receiver to the left of rA

  ! sqrt( 2 pi ) represents a sum of Gaussians in free space
  IF ( Beam%RunType( 4 : 4 ) == 'R' ) THEN
     Ratio1 = SQRT( ABS( COS( alpha ) ) ) / SQRT( 2. * pi )   ! point source
  ELSE
     Ratio1 = 1 / SQRT( 2. * pi )                             ! line  source
  END IF

  Stepping: DO is = 2, Beam%Nsteps

     rB    = ray2D( is     )%x( 1 )
     x_ray = ray2D( is - 1 )%x

     ! compute normalized tangent (compute it because we need to measure the step length)
     rayt = ray2D( is )%x - ray2D( is - 1 )%x
     rlen = NORM2( rayt )
     IF ( rlen < 1.0D3 * SPACING( ray2D( is )%x( 1 ) ) ) CYCLE Stepping  ! if duplicate point in ray, skip to next step along the ray
     rayt = rayt / rlen
     rayn = [ -rayt( 2 ), rayt( 1 ) ]      ! unit normal to ray

     dqds   = ray2D( is )%q( 1 ) - ray2D( is - 1 )%q( 1 )
     dtauds = ray2D( is )%tau    - ray2D( is - 1 )%tau

     ! phase shifts at caustics
     qnew  = ray2D( is - 1 )%q( 1 )
     IF ( qnew <= 0.0 .AND. qOld > 0.0 .OR. &
          qnew >= 0.0 .AND. qOld < 0.0 ) phase = phase + pi / 2.    ! phase shifts at caustics
     qold = qnew

     ! calculate beam width
     lambda = ray2D( is - 1 )%c / ( omega / ( 2 * pi ) )
     sigma  = MAX( ABS( ray2D( is - 1 )%q( 1 ) ), ABS( ray2D( is )%q( 1 ) ) ) / q0 / ABS( rayt( 1 ) ) ! beam radius projected onto vertical line
     sigma  = MAX( sigma, MIN( 0.2 * is * Beam%deltas / lambda, pi * lambda ) )
     RadMax = BeamWindow * sigma

     zmin   = min( ray2D( is - 1 )%x( 2 ), ray2D( is )%x( 2 ) ) - RadMax   ! min depth of ray segment
     zmax   = max( ray2D( is - 1 )%x( 2 ), ray2D( is )%x( 2 ) ) + RadMax   ! max depth of ray segment

     ! is this a steep ray?
     ! If so, don't try to get depth limits: it's too complicated
     IF ( ABS( rayt( 1 ) ) < 0.5 ) THEN
        zmin = -HUGE( zmin )
        zmax = +HUGE( zmax )
     END IF

     ! compute beam influence for this segment of the ray
     RcvrRanges: DO
        ! is r( ir ) contained in [ rA, rB )? Then compute beam influence
        IF ( Pos%r( ir ) >= MIN( rA, rB ) .AND. Pos%r( ir ) < MAX( rA, rB ) ) THEN

           IF ( Beam%RunType( 5 : 5 ) == 'I' ) THEN
              x_rcvr( 1, : ) = [ Pos%r( ir ), Pos%rd( ir ) ]   ! irregular   grid
           ELSE
              do id = 1, Nrd_per_range
                 x_rcvr( id, : ) = [ Pos%r( ir ), Pos%rd( id ) ]   ! rectilinear grid
              end do
           END IF
           
           !IF ( x_rcvr( 2 ) < zmin .OR. x_rcvr( 2 ) > zmax ) CYCLE RcvrDepths

           do concurrent ( id = 1 : Nrd_per_range )
              s( id )      =      DOT_PRODUCT( x_rcvr( id, : ) - x_ray, rayt ) / rlen  ! proportional distance along ray
              n( id )      = ABS( DOT_PRODUCT( x_rcvr( id, : ) - x_ray, rayn ) )       ! normal distance to ray
           end do
        
           q      = ray2D( is - 1 )%q( 1 ) + s * dqds                ! interpolated amplitude
           sigma  = ABS( q / q0 )                                    ! beam radius

           ! calculate the beamwidth (must be at least pi * lambda, except in the nearfield
           sigma  = MAX( sigma, MIN( 0.2 * is * Beam%deltas / lambda, pi * lambda ) )

           WHERE ( n( : ) < BeamWindow * sigma( : ) )   ! Within beam window?
              A     = ABS( q0 / q )
              delay = ray2D( is - 1 )%tau + s * dtauds     ! interpolated delay
              const = Ratio1 * SQRT( ray2D( is )%c / ABS( q ) ) * ray2D( is )%Amp
              W     = EXP( -0.5 * ( n / sigma ) ** 2 ) / ( sigma * A )   ! Gaussian decay
              Amp   = const * W
           END WHERE

        RcvrDepths: DO id = 1, Nrd_per_range

           IF ( n( id ) < BeamWindow * sigma( id ) ) THEN   ! Within beam window?

              phaseInt = phase
              IF ( q( id ) <= 0.0d0 .AND. qOld > 0.0d0 .OR. &
                   q( id ) >= 0.0d0 .AND. qOld < 0.0d0 ) phaseInt = phase + pi / 2.  ! phase shifts at caustics

              SELECT CASE( Beam%RunType( 1 : 1 ) )
              CASE ( 'E' )                ! eigenrays
                 CALL WriteRay( SrcAngle, is )
              CASE ( 'A', 'a' )           ! arrivals
                 RcvrAngle  = RadDeg * ATAN2( rayt( 2 ), rayt( 1 ) )
                 CALL AddArr( omega, id, ir, Amp( id ), ray2D( is )%Phase + phaseInt, delay( id ), &
                      SrcAngle, RcvrAngle, ray2D( is )%NumTopBnc, ray2D( is )%NumBotBnc )
              CASE ( 'C' )                ! coherent TL
                 U( id, ir ) = U( id, ir ) + CMPLX( Amp( id ) * EXP( -i * ( omega * delay( id ) - ray2D( is )%Phase - phaseInt ) ) )
              CASE DEFAULT                ! incoherent/semicoherent TL
                 U( id, ir ) = U( id, ir ) + SNGL( 2. * const( id ) ** 2 * W( id ) )
              END SELECT
           END IF
        END DO RcvrDepths
        END IF

        ! receiver not bracketted; bump receiver index, ir, towards rB
        IF ( rB > Pos%r( ir ) ) THEN
           IF ( ir >= Pos%Nr        ) EXIT   ! go to next step on ray
           irTT = ir + 1                     ! bump right
           IF ( Pos%r( irTT ) >= rB ) EXIT   ! go to next step on ray
        ELSE
           IF ( ir <= 1             ) EXIT   ! go to next step on ray
           irTT = ir - 1                     ! bump left
           IF ( Pos%r( irTT ) <= rB ) EXIT   ! go to next step on ray
        END IF
        ir = irTT

     END DO RcvrRanges

     rA = rB
  END DO Stepping

END SUBROUTINE InfluenceGeoGaussian

! **********************************************************************!

SUBROUTINE InfluenceSGB( U, alpha, Dalpha )

  ! Computes the beam influence, i.e. the contribution of a single beam to the complex pressure
  ! This version uses a beam representation in Cartesian coordinates
  ! This is an implementation of Bucker's Simple Gaussian Beams

  USE bellhopMod
  USE SdRdRMod
  USE ArrMod

  IMPLICIT NONE
  INTEGER           :: id, ir, is
  REAL (KIND=8)     :: x( 2 ), rayt( 2 ), A, beta, cn, CPA, deltaz, DS, phase, q, qOld, rA,rB, Ratio1, &
                       sint, SX1, thet, W, SrcAngle, alpha, Dalpha
  COMPLEX           :: U( Nrd_per_range, Pos%Nr )
  COMPLEX (KIND=8)  :: contri, delay, tau

  Ratio1 = SQRT(  COS( alpha ) )
  phase  = 0
  qOld   = 1.0
  BETA   = 0.98  ! Beam Factor
  A      = -4.0 * LOG( BETA ) / Dalpha**2
  CN     = Dalpha * SQRT( A / pi )
  rA     = ray2D( 1 )%x( 1 )
  ir     = 1

  Stepping: DO is = 2, Beam%Nsteps

     rB = ray2D( is )%x( 1 )

     ! phase shifts at caustics
     q  = ray2D( is - 1 )%q( 1 )
     IF ( q < 0.0d0 .AND. qOld >= 0.0d0 .OR. &
          q > 0.0d0 .AND. qOld <= 0.0d0 ) phase = phase + pi / 2.  ! phase shifts at caustics
     qold = q

     RcvrRanges: DO WHILE ( ABS( rB - rA ) > 1.0D3 * SPACING( rA ) .AND. rB > Pos%r( ir ) )   ! Loop over bracketted receiver ranges

        W     = ( Pos%r( ir ) - rA ) / ( rB - rA )
        x     = ray2D( is - 1 )%x      + W * ( ray2D( is )%x      - ray2D( is - 1 )%x )
        rayt  = ray2D( is - 1 )%t      + W * ( ray2D( is )%t      - ray2D( is - 1 )%t )
        q     = ray2D( is - 1 )%q( 1 ) + W * ( ray2D( is )%q( 1 ) - ray2D( is - 1 )%q( 1 ) )
        tau   = ray2D( is - 1 )%tau    + W * ( ray2D( is )%tau    - ray2D( is - 1 )%tau )

        ! following is incorrect because ray doesn't always use a step of deltas
        SINT  =  ( is - 1 ) * Beam%deltas + W * Beam%deltas

        IF ( q < 0.0d0 .AND. qOld >= 0.0d0 .OR. q > 0.0d0 .AND. qOld <= 0.0d0 ) phase = phase + pi / 2. ! phase shifts at caustics

        RcvrDepths: DO id = 1, Nrd_per_range
           deltaz =  Pos%rd( id ) - x( 2 )   ! ray to rcvr distance
           ! Adeltaz    = ABS( deltaz )
           ! IF ( Adeltaz < RadMax ) THEN
           SELECT CASE( Beam%RunType( 1 : 1 ) )
           CASE ( 'E' )         ! eigenrays
              SrcAngle = RadDeg * alpha   ! take-off angle in degrees
              CALL WriteRay( SrcAngle, is )
           CASE DEFAULT         ! coherent TL
              CPA    = ABS( deltaz * ( rB - rA ) ) / SQRT( ( rB - rA )**2 + ( ray2D( is )%x( 2 ) - ray2D( is - 1 )%x( 2 ) )**2  )
              DS     = SQRT( deltaz **2 - CPA **2 )
              SX1    = SINT + DS
              thet   = ATAN( CPA / SX1 )
              delay  = tau + rayt( 2 ) * deltaz
              contri = Ratio1 * CN * ray2D( is )%Amp * EXP(-A * thet ** 2 - &
                       i * ( omega * delay - ray2D( is )%Phase - phase ) ) / SQRT( SX1 )
              U( id, ir ) = U( id, ir ) + CMPLX( contri )
           END SELECT
           ! END IF
        END DO RcvrDepths

        qOld = q
        ir   = ir + 1
        IF ( ir > Pos%Nr ) RETURN
     END DO RcvrRanges

     rA = rB
  END DO Stepping

END SUBROUTINE InfluenceSGB

! **********************************************************************!

SUBROUTINE BranchCut( q1C, q2C, BeamType, KMAH )

  ! Checks for a branch cut crossing and updates KMAH accordingly

  IMPLICIT NONE
  COMPLEX  (KIND=8), INTENT( IN )    :: q1C, q2C
  CHARACTER (LEN=4), INTENT( IN )    :: BeamType
  INTEGER,           INTENT( INOUT ) :: KMAH
  REAL     (KIND=8)                  :: q1, q2

  SELECT CASE ( BeamType( 2 : 2 ) )
  CASE ( 'W' )   ! WKBeams
     q1 = REAL( q1C )
     q2 = REAL( q2C )
     IF ( ( q1 < 0.0 .AND. q2 >= 0.0 ) .OR. &
          ( q1 > 0.0 .AND. q2 <= 0.0 ) ) KMAH = -KMAH
  CASE DEFAULT
     IF ( REAL( q2C ) < 0.0 ) THEN
        q1 = AIMAG( q1C )
        q2 = AIMAG( q2C )
        IF ( ( q1 < 0.0 .AND. q2 >= 0.0 ) .OR. &
             ( q1 > 0.0 .AND. q2 <= 0.0 ) ) KMAH = -KMAH
     END IF
  END SELECT

END SUBROUTINE BranchCut

! **********************************************************************!

FUNCTION Hermite( x, x1, x2 )

  ! Calculates a smoothing function based on the h0 hermite cubic
  ! x is the point where the function is to be evaluated
  ! returns:
  ! [  0, x1  ] = 1
  ! [ x1, x2  ] = cubic taper from 1 to 0
  ! [ x2, inf ] = 0

  IMPLICIT NONE
  REAL (KIND=8 ), INTENT( IN  ) :: x, x1, x2
  REAL (KIND=8 )                :: Hermite   ! returned function value (Hermite interpolate)
  REAL (KIND=8 )                :: Ax, u

  Ax  = ABS( x  )

  IF ( Ax <= x1 ) THEN
     Hermite = 1.0d0
  ELSE IF ( Ax >= x2 ) THEN
     Hermite = 0.0d0
  ELSE
     u       = ( Ax - x1 ) / ( x2 - x1 )
     Hermite = ( 1.0d0 + 2.0d0 * u ) * ( 1.0d0 - u ) ** 2
  END IF

  !hermit = hermit / ( 0.5 * ( x1 + x2 ) )

END FUNCTION Hermite

! **********************************************************************!

SUBROUTINE ScalePressure( Dalpha, c, r, U, Nrd, Nr, RunType, freq )

  ! Scale the pressure field

  IMPLICIT NONE
  REAL,              PARAMETER       :: pi = 3.14159265
  INTEGER,           INTENT( IN    ) :: Nrd, Nr
  REAL,              INTENT( IN    ) :: r( Nr )      ! ranges
  REAL     (KIND=8), INTENT( IN    ) :: Dalpha, freq, c ! angular spacing between rays, source frequency, nominal sound speed
  COMPLEX,           INTENT( INOUT ) :: U( Nrd, Nr )    ! Pressure field
  CHARACTER (LEN=5), INTENT( IN    ) :: RunType
  INTEGER                            :: ir
  REAL     (KIND=8)                  :: const, factor

  ! Compute scale factor for field
  SELECT CASE ( RunType( 2 : 2 ) )
  CASE ( 'C' )   ! Cerveny Gaussian beams in Cartesian coordinates
     const = -Dalpha * SQRT( freq ) / c
  CASE ( 'R' )   ! Cerveny Gaussian beams in Ray-centered coordinates
     const = -Dalpha * SQRT( freq ) / c
  CASE DEFAULT
     const = -1.0
  END SELECT

  IF ( RunType( 1 : 1 ) /= 'C' ) U = SQRT( REAL( U ) ) ! For incoherent run, convert intensity to pressure

  ! scale and/or incorporate cylindrical spreading
  Ranges: DO ir = 1, Nr
     IF ( RunType( 4 : 4 ) == 'X' ) THEN   ! line source
        factor = -4.0 * SQRT( pi ) * const
     ELSE                                  ! point source
        IF ( r ( ir ) == 0 ) THEN
           factor = 0.0D0                  ! avoid /0 at origin, return pressure = 0
        ELSE
           factor = const / SQRT( ABS( r( ir ) ) )
        END IF
     END IF
     U( :, ir ) = SNGL( factor ) * U( :, ir )
  END DO Ranges

END SUBROUTINE ScalePressure
