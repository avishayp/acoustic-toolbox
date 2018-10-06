SUBROUTINE Evaluate( C, phi, Nz, R, Nr, rr, k, M, Option, P )

  ! Given modes and wavenumbers, compute pressure field
  ! Normalized to pressure of point source at 1 meter
  ! Option:
  ! X     Cartesian   (x, z) coordinates
  ! S scaled cylindrical (sqrt(r) removed)
  ! R     Cylindrical (r, z) coordinates

  IMPLICIT NONE
  INTEGER, PARAMETER    :: MaxM = 20000, MinExp = -100
  REAL,    PARAMETER    :: pi = 3.1415926
  COMPLEX, PARAMETER    :: i = ( 0.0, 1.0 )
  INTEGER, INTENT( IN ) :: M, Nr, Nz
  REAL,    INTENT( IN ) :: rr( Nz ), r( Nr )
  COMPLEX, INTENT( IN ) :: C( M ), phi( MaxM, Nz ), k( MaxM )
  COMPLEX, INTENT( OUT) :: P( Nz, Nr )
  CHARACTER (LEN=50), INTENT( IN ) :: Option
  INTEGER               :: ir, iz
  COMPLEX               :: Hank( M ), ik( M ), const( M ), Cmat( M, Nz ), factor

  ! If no modes, return vanishing pressure
  IF ( M <= 0 ) THEN
     P = 0.0
     RETURN
  END IF

  ! Initialization
  factor = i * SQRT( 2.0 * pi ) * EXP( i * pi / 4.0 )

  IF ( Option( 1 : 1 ) == 'X' ) THEN   ! line source
     const( 1 : M ) = factor * C( 1 : M ) / k( 1 : M )
  ELSE                                 ! point source or scaled cylindrical
     const( 1 : M ) = factor * C( 1 : M ) / SQRT( k( 1 : M ) )
  END IF

  ik( 1 : M ) = -i * k( 1 : M )   ! use e{ i( w t - k r ) } form
  IF ( Option( 4 : 4 ) == 'I' ) ik = REAL( ik )   ! Incoherent case

  DO iz = 1, Nz
     Cmat( :, iz ) = const( : ) * phi( 1 : M, iz ) * EXP( ik( : ) * rr( iz ) )
  END DO

  Ranges: DO ir = 1, Nr
     ! eliminate underflows (can raise CPU time)
     !WHERE (  REAL( ik * r( ir ) ) > MinExp )
     Hank = EXP( ik * r( ir ) )
     !ELSEWHERE
     !   Hank = 0.0
     !END WHERE

!!$     Depths: DO iz = 1, Nz
!!$        IF ( Option( 4 : 4 ) /= 'I' )  THEN     ! coherent   case
!!$           T =       SUM(    Cmat( :, iz ) * Hank( : ) )
!!$        ELSE                                    ! incoherent case
!!$           T = SQRT( SUM(  ( Cmat( :, iz ) * Hank( : ) ) ** 2 ) )
!!$        END IF
!!$        ! Cylindrical or cartesian coordinates?
!!$        !IF ( Option( 1 : 1 ) == 'R' .AND. ABS( r( ir ) + rr( iz ) ) > TINY( R( 1 ) ) ) T = T / SQRT( r( ir ) + rr( iz ) )
!!$        P( iz, ir ) = T
!!$     END DO Depths

     IF ( Option( 4 : 4 ) /= 'I' )  THEN     ! coherent   case
        Depths: DO iz = 1, Nz
            P( iz, ir ) =      SUM(    Cmat( :, iz ) * Hank( : ) )
        END DO Depths
     ELSE                                    ! incoherent case
        DepthsInc: DO iz = 1, Nz
           P( iz, ir ) = SQRT( SUM(  ( Cmat( :, iz ) * Hank( : ) ) ** 2 ) )
        END DO DepthsInc
     END IF

     ! Cylindrical or cartesian coordinates?
     IF ( Option( 1 : 1 ) == 'R' ) THEN
        WHERE ( ABS( r( ir ) + rr( : ) ) > TINY( R( 1 ) ) )
           P( :, ir ) = P( :, ir ) / SQRT( r( ir ) + rr( : ) )
        END WHERE
     END IF
     
  END DO Ranges

END SUBROUTINE Evaluate
