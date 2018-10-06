SUBROUTINE Influence3D( xs, alpha, iBeamWindow, P )

  ! Computes the beam influence, i.e. 
  ! the contribution of a single beam to the complex pressure

  USE bellhopMod
  USE SdRdRMod
  USE norms

  IMPLICIT NONE
  INTEGER,         INTENT( IN  ) :: iBeamWindow
  REAL ( KIND=8 ), INTENT( IN  ) :: alpha                  ! declination angle
  REAL ( KIND=8 ), INTENT( IN  ) :: xs( 3 )                ! source coordinate
  COMPLEX        , INTENT( OUT ) :: P( Ntheta, Nrd, Nr )   ! complex pressure
  INTEGER            :: KMAHV( MaxN ), Itheta, ir, is, id, ir1, ir2, KMAH, image, iBeamWindow2
  REAL    ( KIND=8 ) :: NA, NB, NSquared, LSquared1, mA, mB, MSquared, LSquared2, nm, &
       rA, rB, ST, CT, tau, W, xt( 3 ), zT, c, delta, e1( 3 ), e2( 3 ), QL, QR
  COMPLEX ( KIND=8 ) :: F, G, H, DetQ, contri

  iBeamWindow2 = iBeamWindow

  ! During reflection imag( q ) is constant and adjacent normals cannot bracket a segment of the TL
  ! line, so no special treatment is necessary

  ! *** Begin by forming KMAH index ***

  KMAHV( 1 ) = 1
  DO is = 2, Beam%Nsteps
     KMAHV( is ) = KMAHV( is - 1 )

     IF ( REAL( ray3D( is )%DetQ ) < 0.0 ) THEN
        QL = AIMAG( ray3D( is - 1 )%DetQ )
        QR = AIMAG( ray3D( is     )%DetQ )
        IF ( ( QL < 0.0 .AND. QR >= 0.0 ) .OR. &
             ( QL > 0.0 .AND. QR <= 0.0 ) ) KMAHV( is ) = -KMAHV( is )
     ENDIF

  END DO

  ReceiverDepths: DO id = 1, Nrd

     Radials: DO Itheta = 1, Ntheta
        CT = COS( DegRad * theta( Itheta ) )
        ST = SIN( DegRad * theta( Itheta ) )

        ! image: DO image = 1, 1
        image = 1
        ir1   = HUGE( ir1 )

        Stepping: DO is = 2, Beam%Nsteps

           ! Compute ray normals e1 and e2
           CALL RayNormal( ray3D( is )%t, ray3D( is )%phi, ray3D( is )%c, e1, e2 )

           ! Compute coordinates of intercept: nB, mB, rB
           xt = ray3D( is )%x - xs
           zT = ray3D( is )%x( 3 ) - rd( id )

           delta = -CT * ( e1( 2 ) * e2( 3 ) - e1( 3 ) * e2( 2 ) ) + ST * ( e1( 1 ) * e2( 3 ) - e1( 3 ) * e2( 1 ) )

           ! Check for ray normal || radial of rcvr line
           IF ( ABS( delta ) < 1.0E-10 ) EXIT Stepping
           nB = ( -CT * ( xt( 2 ) * e2( 3 ) - zT * e2( 2 ) ) + ST * ( xt( 1 ) * e2( 3 ) - zT * e2( 1 ) ) ) / delta
           mB = ( -CT * ( e1( 2 ) * zT - e1( 3 ) * xt( 2 ) ) + ST * ( e1( 1 ) * zT - e1( 3 ) * xt( 1 ) ) ) / delta
           rB = ( -xt( 1 ) * ( e1( 2 ) * e2( 3 ) - e1( 3 ) * e2( 2 ) ) + &
                   xt( 2 ) * ( e1( 1 ) * e2( 3 ) - e1( 3 ) * e2( 1 ) ) - &
                        zT * ( e1( 1 ) * e2( 2 ) - e1( 2 ) * e2( 1 ) ) ) / delta
           ir2 = MAX( MIN( INT( ( rB - r( 1 ) ) / Delta_r ) + 1, Nr ), 1 )   ! index of nearest rcvr before normal

           ! *** Compute contributions to bracketted receivers ***

           IF ( ( ir1 < ir2 ) .AND. ( NORM2( ray3D( is )%x - ray3D( is - 1 )%x ) > 1.0D3 * SPACING( ray3D( is )%x( 1 ) ) ) ) THEN

              DO ir = ir1 + 1, ir2
                 W   = ( r( ir ) - rA ) / ( rB - rA )
                 nSquared = ( nA + W * ( nB - nA ) )**2
                 mSquared = ( mA + W * ( mB - mA ) )**2

                 ! Within beam window?
                 LSquared1 = -2.0 / ( omega * AIMAG( ray3D( is )%f / ray3D( is )%DetQ ) )
                 LSquared2 = -2.0 / ( omega * AIMAG( ray3D( is )%g / ray3D( is )%DetQ ) )

                 IF ( nSquared < iBeamWindow2 * LSquared1 .AND. mSquared < iBeamWindow2 * LSquared2  ) THEN

                    f    = ray3D( is - 1 )%f    + W * ( ray3D( is )%f    - ray3D( is - 1 )%f    )
                    g    = ray3D( is - 1 )%g    + W * ( ray3D( is )%g    - ray3D( is - 1 )%g    )
                    h    = ray3D( is - 1 )%h    + W * ( ray3D( is )%h    - ray3D( is - 1 )%h    )
                    DetQ = ray3D( is - 1 )%DetQ + W * ( ray3D( is )%DetQ - ray3D( is - 1 )%DetQ )
                    c    = ray3D( is - 1 )%c    + W * ( ray3D( is )%c    - ray3D( is - 1 )%c    )
                    tau  = ray3D( is - 1 )%tau  + W * ( ray3D( is )%tau  - ray3D( is - 1 )%tau  )
                    KMAH = KMAHV( is - 1 )

                    IF ( REAL( DetQ ) < 0.0 ) THEN
                       QL = AIMAG( ray3D( is - 1 )%DetQ )
                       QR = AIMAG(                 DetQ )
                       IF ( ( QL < 0.0 .AND. QR >= 0.0 ) .OR. &
                            ( QL > 0.0 .AND. QR <= 0.0 ) ) KMAH = -KMAH
                    ENDIF

                    nm     = SQRT( nSquared * mSquared )
                    contri = ABS( COS( alpha ) ) * ray3D( is )%Amp * SQRT( c / DetQ ) * &
                         EXP( -i * omega *( tau + 0.5 * ( f * nSquared + 2.0 * h * nm + g * mSquared ) / DetQ - ray3D( is )%phase) )

                    IF ( KMAH  < 0  ) contri = -contri
                    IF ( image == 2 ) contri = -contri
                    P( Itheta, id, ir ) = P( Itheta, id, ir ) + contri
                 ENDIF
              END DO
           ENDIF

           nA  = nB
           mA  = mB
           rA  = rB
           ir1 = ir2
        END DO Stepping
        ! END DO image
     END DO Radials
  END DO ReceiverDepths

END SUBROUTINE Influence3D
