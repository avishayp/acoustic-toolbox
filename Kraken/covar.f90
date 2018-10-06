PROGRAM COVAR 

!!! this version has not been tested since converting to using the module SRDRD
  !     Takes replica shade files and computes synthesized                
  !     data vectors with Gaussian random noise added                     
  !     Averages data vector to from covariance matrix                    

  !     Data vectors are produced for the phone positions in              
  !     the shade file (not those of the mode file).                      

  USE SdRdRMod                                              
  INTEGER, PARAMETER :: INPFIL = 5, PRTFIL = 6, SHDFIL = 20, DATFIL = 21, COVFIL = 22, &
       MaxM = 1500, MaxNR = 5000, MaxPhones = 101
  REAL,    PARAMETER :: pi = 3.141592                                                         
  INTEGER   iSeed( 4 ), MIN_LOC( 1 )
  REAL      CovarianceISO( MaxPhones, MaxPhones )
  COMPLEX   P( MaxNR ), phiS( MaxM ), phiR( MaxM, MaxPhones ), k( MaxM )
  COMPLEX, ALLOCATABLE :: d( : ), dNoisy( : ), Noise( : ), Covariance( :, : )
  CHARACTER Title*80, Title2*80, PlotType*10, NoiseType*1 

  SAVE iSeed 
  DATA iSeed /2301, 4562, 3219, 1823/ 

  ! *** Read data file header ***                                     

  CALL ReadHeader( SHDFIL, ' ', Title, Freq, Atten, deltaR, PlotType, xs, ys, theta )

  OPEN ( FILE = 'DATFIL', UNIT = DATFIL, STATUS = 'NEW' ) 
  OPEN ( FILE = 'COVFIL', UNIT = COVFIL, STATUS = 'NEW' ) 

  ! Read user specification for source position and compute indices   

  READ( INPFIL, * ) Title 
  READ( INPFIL, * ) SDepth, SRange 
  READ( INPFIL, * ) NoiseType, SNRDB 
  READ( INPFIL, * ) NTimes, iSeed( 4 ) 
  !     NTimes = 30000  ! **********************************              
  !     SNRDB  = 200.0  ! **********************************             
  ! Convert km to meters          
  SRange = 1000.0 * SRange 
  SNR    = 10.0 ** ( SNRDB / 10.0 ) 

  ! *** Calculate source depth index ***                              

  MIN_LOC = MINLOC( ABS( SD( 1 : NSD ) - SDepth ) )
  iSDepth = MIN_LOC( 1 )

  ! Check to see if within tolerance                              

  IF ( Freq == 0.0 ) THEN 
     Tolerance = 0.001 
  ELSE
     Tolerance = 1500.0 / ( 6.0 * Freq )   ! one wavelength
  END IF

  IF ( ABS( SD( iSDepth ) - SDepth ) > Tolerance ) THEN 
     WRITE( PRTFIL, * ) 'Sources exist at the following depths:' 
     WRITE( PRTFIL, * ) ( SD( I ), I = 1, NSD ) 
     CALL ERROUT( PRTFIL, 'F', 'COVAR', 'No source at specified depth' )                  
  END IF

  ! *** Calculate source range index ***                              

  MIN_LOC = MINLOC( ABS( R( 1 : NR ) - SRange ) )
  iSRange = MIN_LOC( 1 )

  ! Check to see if within tolerance                              

  IF ( ABS( R( iSRange ) - SRange ) > Tolerance ) THEN 
     WRITE( PRTFIL, * ) 'Sources exist at the following ranges:' 
     WRITE( PRTFIL, * ) R( 1 : NR ) 
     CALL ERROUT( PRTFIL, 'F', 'COVAR', 'No source at specified range' )                  
  END IF

  ! extract data vector                                           

  WRITE( PRTFIL, * ) 'Data vector being extracted for SDEP, SRAN', SD( iSDepth), R( iSRange )
  ALLOCATE( d( NRD ), dNoisy( NRD ) )

  IREC = ( iSDepth - 1 ) * NRD + 6 
  DO IRD1 = 1, NRD
     IREC = IREC + 1 
     READ ( SHDFIL, REC = IREC ) P( 1 : NR ) 
     d( IRD1 ) = P( iSRange )
  END DO

  ! Normalize data vector so that array power = 1                 

  d = d / SQRT( DOT_PRODUCT( d, d ) )

  PhonePower = 1.0 / NRD 

  ! Zero out covariance matrix                                    

  ALLOCATE( Noise( NRD ), Covariance( NRD, NRD ) )
  Noise      = 0.0
  Covariance = 0.0

  ! *** If making colored noise, need modal data ***                  

  IF ( NoiseType == 'C' ) THEN 
     IPROF   = 1
     SD( 1 ) = 1.0  ! noise sources are at 1 m
     NSD     = 1

     CALL GETMOD( IPROF, ' ', MaxM, SD, NSD, 'N', k, phiS, M, Freq, Title2 )                                 
     CALL GETMOD( IPROF, ' ', MaxM, RD, NRD, 'N', k, phiR, M, Freq, Title2 )                                 

     IF ( ANY( AIMAG( k( 1:M ) ) >= 0.0 ) ) THEN 
        CALL ERROUT( PRTFIL, 'F', 'COVAR', 'Modes present with no loss: Noise field is infinite ' )   
     END IF

     ! Construct phi and Power
     DO mode = 1, M 
        phiR( mode, 1 : NRD ) = phiS( mode ) * phiR( mode, 1 : NRD ) / &
             SQRT( -REAL( k( mode ) ) * AIMAG( k( mode ) ) )
        Power = SUM( phiR( mode, 1 : NRD ) * CONJG( phiR( mode, 1 : NRD ) ) )
     END DO
  END IF

  ! header for data vector file                                   

  WRITE( DATFIL, * ) '''', Title(1:60), '''' 
  WRITE( DATFIL, * ) Freq 
  WRITE( DATFIL, * ) NTimes, NRD 

  DO IRD1 = 1, NRD
     WRITE( DATFIL, * ) RD( IRD1 )
  END DO

  ! *** Main loop: generate realizations of noise field ***           

  AverageSignal = 0.0
  AverageNoise  = 0.0
  sigma         = SQRT( PhonePower / SNR )   ! Noise level on a phone

  DO iTime = 1, NTimes 

     IF (      NoiseType == 'W' ) THEN 
        CALL White( iSeed, sigma, Noise, NRD ) 
     ELSE IF ( NoiseType == 'C' ) THEN 
        CALL Color( iSeed, phiR, M, NRD, sigma, Power, Noise ) 
     ELSE 
        STOP 'Unknown noise type' 
     END IF

     dNoisy = d + Noise

     AverageSignal = SUM( ABS( d     ) ** 2 )
     AverageNoise  = SUM( ABS( Noise ) ** 2 )

     ! write noise vector to DATFIL                               
     !         WRITE( DATFIL, 1000 ) ( iTime, IRD, dNoisy( IRD ) , IRD = 1, NRD )               
1000 FORMAT( I3, I3, ' (', F11.7, ',', F11.7, ' )' ) 

     ! Compute cross-sensor correlation matrix R                  

     DO I = 1, NRD
        Covariance( I, : ) = Covariance( I, : ) + dNoisy( I ) * CONJG( dNoisy( : ) )
     END DO

  END DO

  AverageSignal = AverageSignal / ( NTimes * NRD ) 
  AverageNoise  = AverageNoise  / ( NTimes * NRD ) 

  WRITE( PRTFIL, * ) 'AverageSignal = ', AverageSignal 
  WRITE( PRTFIL, * ) 'AverageNoise  = ', AverageNoise 

  WRITE( PRTFIL, * ) 'Average signal (dB) = ', 10.0 * LOG10( AverageSignal )                                      
  WRITE( PRTFIL, * ) 'Average noise  (dB) = ', 10.0 * LOG10( AverageNoise  )                                      
  WRITE( PRTFIL, * ) 'Average SNR    (dB) = ', 10.0 * LOG10( AverageSignal / AverageNoise )                           

  ! *** Analytic formula for covariance matrix in isovelocity ocean

  CovarianceISO( 1:NRD, 1:NRD ) = 0.0

  DO I = 1, NRD
     DO J = 1, NRD
        Depth = 100.0 
        DO mode = 1, M 
           gamma = ( mode - .5 ) * pi / Depth 
           CovarianceISO( I, J ) = CovarianceISO( I, J ) + SIN( gamma * SD( 1 ) ) ** 2 * &
                SIN( gamma * RD( I ) ) * SIN( gamma * RD( J ) )
        END DO
     END DO
  END DO

  ! *** Write normalized cross-sensor correlation matrix column-wise  

  WRITE( COVFIL, * ) '''', Title(1:60), '''' 
  WRITE( COVFIL, * ) Freq 
  WRITE( COVFIL, * ) NRD 
  DO IRD1 = 1, NRD
     WRITE( COVFIL, * ) RD( IRD1 )
  END DO

  DO I = 1, NRD 
     DO J = 1, NRD 
        Covariance( I, J ) = NRD * Covariance( I,J ) / ( ( 1.0 + sigma**2 ) * NTimes )                     
        WRITE( COVFIL, 1000 ) I, J, Covariance( I, J ) 
        !               ABS( CovarianceISO( I, J ) * Covariance( 1, 1 ) / CovarianceISO( 1, 1 ) ) 
     END DO
  END DO

  STOP 
END PROGRAM COVAR

!----------------------------------------------------------------------C

SUBROUTINE White( iSeed, sigma, Noise, NRD ) 

  ! Generates Gaussian random white noise                             

  COMPLEX i 
  PARAMETER ( i = (0.0, 1.0), pi = 3.1415926 ) 

  INTEGER   iSeed( 4 ) 
  COMPLEX   Noise( * ) 

  DO IRD = 1, NRD
     CALL RANDOM_NUMBER( X )
     CALL RANDOM_NUMBER( Y ) 
     X = X + 0.000001
     Noise( IRD ) = sigma * SQRT( -LOG( X ) ) * EXP( i * 2.0 * pi * Y )
  END DO

  ! Note power = NRD * sigma ** 2                                 

  RETURN 
END SUBROUTINE White

!----------------------------------------------------------------------C

SUBROUTINE Color( iSeed, phi, M, NRD, sigma, Power, Noise ) 

  ! Generates colored noise vector                                    

  COMPLEX i
  PARAMETER ( MaxM = 1500, i = (0.0, 1.0), pi = 3.1415926 ) 

  INTEGER iSeed( 4 ) 
  COMPLEX phi( MaxM, NRD ), Noise( NRD ), Coef

  Noise( 1 : NRD ) = 0.0  ! Zero out noise vector

  ! loop over modes                                               
  DO mode = 1, M 

     ! compute the random coefficient                             
     CALL RANDOM_NUMBER( X )
     CALL RANDOM_NUMBER( Y ) 
     X = X + 0.000001

     Coef = SQRT( -LOG( X ) ) * EXP( i * 2.0 * pi * Y ) 
     Noise = Noise + Coef * phi( mode, : )
  END DO

  ! Normalize so that power = NRD * sigma**2                      

  Noise = sigma * Noise * SQRT( NRD / Power )

  RETURN 
END SUBROUTINE Color
