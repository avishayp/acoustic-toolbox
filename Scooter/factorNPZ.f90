SUBROUTINE Factor( N, d, e, Rv1, Rv2, Rv4 )

! MICHAEL B. PORTER 7/1/85
!
! Based on TINVIT in EISPACK
! Gaussian elimination to factor a symmetric tridiagonal linear system
!
! On input
!    N is the order of the matrix
!    d contains the diagonal elements of the input matrix
!    e              subdiagonal
!      in its last N-1 positions.  e(1) is arbitrary
! On output
!    Rv1, Rv2, Rv3 and Rv4 contain the factored matrix

IMPLICIT NONE
INTEGER          :: i, N
COMPLEX (KIND=8) :: U, V, XU, d( N ), e( N ), Rv1( N ), Rv2( N ), Rv4( N )

! LU decomposition without interchanges

U = d( 1 )
V = e( 2 )

DO I = 2, N-1
   XU         = e( I ) / U
   Rv4( I   ) = XU
   Rv1( I-1 ) = 1.0 / U
   Rv2( I-1 ) = V
   U          = d( I ) - XU * V
   V          = e( I+1 )
END DO

! following is same as above with I = N except V=e(N+1) is not used
XU         = e( N ) / U
Rv4( N   ) = XU
Rv1( N-1 ) = 1.0 / U
Rv2( N-1 ) = V
U          = d(N) - XU * V

IF ( U == 0.0 ) WRITE( *, * ) 'Singular matrix'
Rv1( N ) = 1.0 / U
Rv2( N ) = 0.0

END SUBROUTINE Factor

SUBROUTINE BackSub( N, Rv1, Rv2, Rv4, b )

! MICHAEL B. PORTER 7/1/85

! Based on TINVIT in EISPACK, once upon a time ...
! performs back-substitution for a symmetric tridiagonal linear system
!
! On input
!    N is the order of the matrix
!    RV1,2,3,4 contain the LU factorizationi of A
!    B contains the right hand side (Ax=b)
! On output
!    b contains the solution

IMPLICIT NONE
INTEGER          :: i, N
COMPLEX (KIND=8), INTENT( IN    ) :: Rv1( N ), Rv2( N ), Rv4( N )
COMPLEX (KIND=8), INTENT( INOUT ) :: b( N )

! Forward elimination

DO I = 2, N
   b( I ) = b( I ) - Rv4( I ) * b( I - 1 )
   ! prevent underflows from propagating:
   !IF ( ABS(  REAL( b( I ) ) ) < SMALL .AND. ABS( AIMAG( b( I ) ) ) < SMALL ) b( I ) = 0.0
END DO

! Back-substitution (result in b)
b( N ) = b( N ) * Rv1( N )

IF ( N >= 2 ) THEN
   DO I = N - 1, 1, -1
      b( I ) = ( b( I ) - b( I + 1 ) * Rv2( I ) ) * Rv1( I )
   END DO
END IF

END SUBROUTINE BackSub
