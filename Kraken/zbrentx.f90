SUBROUTINE ZBRENTX( X, A, B, T, ErrorMessage ) 

  ! Michael B. Porter  8/84                                               

  ! FORTRAN conversion of an ALGOL program published in                      
  !    The Computer Journal 14(4):422-425 (1971)                   
  !    BY R. P. Brent

  ! Returns a zero X of the function F in the given interval [ A, B ],
  ! to within a tolerance of 6 * MACHEP * ABS( X ) + 2 * T, where         
  !    MACHEP is the relative machine precision and T is a positive tolerance.
  !    The procedure assumes (and tests) that FUNCT(A) AND FUNCT(B) have different signs.
  !    This is an extended range version that expects the function to have the form                                                 
  !         SUBROUTINE FUNCT( X, G, IPOW )
  !    where G * 10 ** IPOW gives FUNCT( X )
  !    The IPow variable gives a power of 10.

  IMPLICIT NONE
  INTEGER            :: IExpA, IExpB, IExpC
  REAL      (KIND=8) :: X, A, B, C, D, E, T, fa, fb, fc, F1, F2, MACHEP, M, P, Q, R, S, TEN, TOL
  CHARACTER (LEN=80) :: ErrorMessage

  ErrorMessage = ' ' 
  MACHEP = 1.0E-16 
  TEN    = 10.0 

  CALL FUNCT( A, fa, IExpA ) 
  CALL FUNCT( B, fb, IExpB ) 

  IF ( ( (fa > 0.0) .AND. (fb > 0.0) ) .OR.                   &
       ( (fa < 0.0) .AND. (fb < 0.0) ) ) THEN                 
     ErrorMessage = 'Function sign is the same at the interval endpoints' 
     RETURN 
  ENDIF

  ! INTERNAL ROOT                                                     

2000 C  = A 
  fc    = fa 
  IExpC = IExpA 
  E     = B - A 
  D     = E 

  ! EXTERNAL ROOT                                                     

  IF ( IExpA < IExpB ) THEN 
     F1 = fc * TEN ** ( IExpC - IExpB ) 
     F2 = fb 
  ELSE 
     F1 = fc 
     F2 = fb * TEN ** ( IExpB - IExpC ) 
  ENDIF

3000 IF ( ABS( F1 ) < ABS( F2 ) ) THEN 
     A     = B 
     B     = C 
     C     = A 
     fa    = fb 
     IExpA = IExpB 
     fb    = fc 
     IExpB = IExpC 
     fc    = fa 
     IExpC = IExpA 
  ENDIF

  TOL = 2.0 * MACHEP * ABS( B ) + T 
  M   = 0.5 * ( C - B ) 
  IF ( ( ABS( M ) > TOL) .AND. ( fb /= 0.0 ) ) THEN 

     ! SEE IF A BISECTION IS FORCED                                  
     IF ( IExpA < IExpB ) THEN 
        F1 = fa * TEN ** ( IExpA - IExpB ) 
        F2 = fb 
     ELSE 
        F1 = fa 
        F2 = fb * TEN ** ( IExpB - IExpA ) 
     ENDIF

     IF ( (ABS(E) < TOL) .OR. ( ABS( F1 ) <= ABS( F2 ) ) ) THEN                            
        E = M 
        D = E 
     ELSE 
        S = fb / fa * TEN ** ( IExpB - IExpA ) 
        IF ( A == C ) THEN 
           ! LINEAR INTERPOLATION                                 
           P = 2.0 * M * S 
           Q = 1.0 - S 
        ELSE 
           ! INVERSE QUADRATIC INTERPOLATION                      
           Q = fa / fc * TEN ** ( IExpA - IExpC ) 
           R = fb / fc * TEN ** ( IExpB - IExpC ) 
           P = S * ( 2.0 * M * Q * ( Q - R ) - ( B - A ) * ( R - 1.0 ) ) 
           Q = ( Q - 1.0 ) * ( R - 1.0 ) * ( S - 1.0 ) 
        ENDIF
        IF ( P > 0.0 ) THEN 
           Q = -Q 
        ELSE 
           P = -P 
        ENDIF
        S = E 
        E = D 
        IF ( ( 2.0 * P < 3.0 * M * Q - ABS( TOL * Q ) ) .AND. ( P < ABS( 0.5 * S * Q ) ) ) THEN                             
           D = P / Q 
        ELSE 
           E = M 
           D = E 
        ENDIF
     ENDIF

     A     = B 
     fa    = fb 
     IExpA = IExpB 

     IF ( ABS( D ) > TOL) THEN 
        B = B + D 
     ELSE 
        IF ( M > 0.0 ) THEN 
           B = B + TOL 
        ELSE 
           B = B - TOL 
        ENDIF
     ENDIF

     CALL FUNCT( B, fb, IExpB ) 
     IF ( ( fb > 0.0 ) .EQV. ( fc > 0.0 ) ) GOTO 2000 
     GOTO 3000 
  ENDIF

  X = B 

END SUBROUTINE ZBRENTX
