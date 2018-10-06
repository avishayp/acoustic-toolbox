FUNCTION CRCI( z, c, alpha, freq, freq0, AttenUnit, beta, fT, bio, NBioLayers )

  ! Converts real wave speed and attenuation to a single
  !  complex wave speed (with positive imaginary part)

  ! AttenUnit
  ! 6 Cases:    N Nepers/meter
  !             M dB/meter      (M for Meters)
  !             m dB/meter with a power law
  !             F dB/m-kHz      (F for Frequency dependent)
  !             W dB/wavelength (W for Wavelength)
  !             Q Q
  !             L Loss parameter
  !
  ! second letter adds volume attenuation according to standard laws:
  !             T for Thorp
  !             B for biological
  !
  ! freq is the current frequency
  ! freq0 is the reference frequency for which the dB/meter was specified (used only for 'm'

  ! Returns
  ! c     real      part of sound speed
  ! alpha imaginary part of sound speed

  USE MathConstantsMod
  IMPLICIT NONE
  INTEGER, PARAMETER               :: MaxBioLayers = 20
  INTEGER                          :: iBio, NBioLayers
  REAL     (KIND=8), INTENT( IN )  :: freq, freq0, alpha, c, z, beta, fT
  CHARACTER (LEN=2), INTENT( IN )  :: AttenUnit
  REAL     (KIND=8)                :: f2, omega, alphaT, Thorp, a
  COMPLEX  (KIND=8)                :: CRCI

  TYPE bioStructure
     REAL (KIND=8) :: Z1, Z2, f0, Q, a0
  END TYPE bioStructure

  TYPE( bioStructure ) :: bio( MaxBioLayers )

  omega = 2.0 * pi * freq

  !  Convert to Nepers/m 
  alphaT = 0.0
  SELECT CASE ( AttenUnit( 1 : 1 ) )
  CASE ( 'N' )
     alphaT = alpha
  CASE ( 'M' )   ! dB/m
     alphaT = alpha / 8.6858896D0
  CASE ( 'm' )
     alphaT = alpha / 8.6858896D0
     IF ( freq < fT ) THEN   ! frequency raised to the power beta
        alphaT = alphaT * ( freq / freq0 ) ** beta
     ELSE                              ! linear in frequency
        alphaT = alphaT * ( freq / freq0 ) * ( fT / freq0 ) ** ( beta - 1 )
     END IF     
  CASE ( 'F' )   ! dB/(m kHz)
     alphaT = alpha * freq / 8685.8896D0
  CASE ( 'W' )   ! dB/wavelength
     IF ( c /= 0.0         ) alphaT = alpha * freq / ( 8.6858896D0 * c )
     !        The following lines give f^1.25 Frequency dependence
     !        FAC = SQRT( SQRT( Freq / 50.0 ) )
     !        IF ( c /= 0.0 ) alphaT = FAC * alpha * Freq / ( 8.6858896D0 * c )
  CASE ( 'Q' )   ! Quality factor
     IF ( c * alpha /= 0.0 ) alphaT = omega / ( 2.0 * c * alpha )
  CASE ( 'L' )   ! loss parameter
     IF ( c /= 0.0         ) alphaT = alpha * omega / c
  END SELECT

  ! added volume attenuation
  SELECT CASE ( AttenUnit( 2 : 2 ) )
  CASE ( 'T' )
     f2 = ( Freq / 1000.0 ) **2

     ! Original formula from Thorp 1967
     ! Thorp = 40.0 * f2 / ( 4100.0 + f2 ) + 0.1 * f2 / ( 1.0 + f2 )   ! dB/kyard
     ! Thorp = Thorp / 914.4D0                 ! dB / m
     ! Thorp = Thorp / 8.6858896D0             ! Nepers / m

     ! Updated formula from JKPS Eq. 1.34
     Thorp = 3.3d-3 + 0.11 * f2 / ( 1.0 + f2 ) + 44.0 * f2 / ( 4100.0 + f2 ) + 3d-4 * f2   ! dB/km
     Thorp = Thorp / 8685.8896 ! Nepers / m

     alphaT = alphaT + Thorp
  CASE ( 'B' )   ! biological attenuation per Orest Diachok
     DO iBio = 1, NBioLayers
        IF ( z >= bio( iBio )%Z1 .AND. z <= bio( iBio )%Z2 ) THEN
           a = bio( iBio )%a0 / ( ( 1.0 - bio( iBio )%f0 ** 2 / freq ** 2  ) ** 2 + 1.0 / bio( iBio )%Q ** 2 )   ! dB/km
           a = a / 8685.8896   ! Nepers / m
           alphaT = alphaT + a
        END IF
     END DO
  END SELECT

  ! Convert Nepers/m to equivalent imaginary sound speed 
  alphaT = alphaT * c * c / omega
  CRCI   = CMPLX( c, alphaT, KIND=8 )

  RETURN
END FUNCTION CRCI
