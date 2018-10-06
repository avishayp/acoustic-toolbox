MODULE KrakencMod

  USE MathConstantsMod
  SAVE

  INTEGER,          PARAMETER   :: ENVFile = 5, PRTFile = 6, MODFile = 20, EVMFile = 22, MaxM = 20000, MaxMedium = 500, NSets = 5
  INTEGER                       :: FirstAcoustic, LastAcoustic, NV( NSets ), ISet, M, &
                                   LRecordLength, IRecProfile = 1, ModeCount, Mode, IProf, ifreq
  REAL      (KIND=8)            :: ET( NSets ), hV( NSets ), VG( MaxM ), cMin, cLow, cHigh, Freq, omega2, RMax
  COMPLEX   (KIND=8)            :: k( MaxM ), EVMat( NSets, MaxM ), Extrap( NSets, MaxM )
  CHARACTER (LEN= 8)            :: TopOpt, BotOpt
  CHARACTER (LEN=80)            :: Title

  ! Halfspace properties
  TYPE HSInfo
     CHARACTER (LEN=1)          :: BC                            ! Boundary condition type
     REAL    (KIND=8)           :: alphaR, alphaI, betaR, betaI  ! P-wave, S-wave speeds (user units)
     REAL    (KIND=8)           :: beta, fT                      ! power law and transition frequency
     COMPLEX (KIND=8)           :: cP, cS                        ! P-wave, S-wave speeds (neper/m loss)
     REAL    (KIND=8)           :: rho, BumpDensity, eta, xi     ! density, boss parameters
  END TYPE HSInfo

  TYPE( HSInfo )                :: HSTop, HSBot

  ! finite-difference grid
  INTEGER                       :: Loc( MaxMedium ), NG( MaxMedium ), N( MaxMedium )
  REAL      (KIND=8)            :: H(   MaxMedium )

  ! storage for finite-difference equations
  REAL    (KIND=8), ALLOCATABLE :: rho( : )
  COMPLEX (KIND=8), ALLOCATABLE :: B1( : ), B2( : ), B3( : ), B4( : )

END MODULE KrakencMod
