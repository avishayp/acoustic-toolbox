MODULE krakelmod

  SAVE

  INTEGER,          PARAMETER   :: ENVFile = 5, PRTFile = 6, MODFile = 20, EVMFile = 22, &
       MaxM = 5000, MaxMedium = 500, NSets = 5, MaxN = 200000
  REAL    (KIND=8), PARAMETER   :: pi = 3.1415926535898D0
  COMPLEX (KIND=8), PARAMETER   :: i = ( 0.0, 1.0 ) 

  INTEGER                       :: NFact, NLact, NMedia, NV( NSets ), ISet, M, LRECL, ModeCount, Mode, IProf
  REAL    (KIND=8)              :: ET( NSets ), hV( NSets ), cMin, cLow, cHigh, Freq, omega2, RMax
  REAL    (KIND=8), ALLOCATABLE :: EVMat( :, : ), Extrap( :, : ), VG( : )
  COMPLEX (KIND=8), ALLOCATABLE :: k( : )

  ! Halfspace properties
  TYPE HSInfo
     CHARACTER (LEN=1)          :: BC                          ! Boundary condition type
     COMPLEX (KIND=8)           :: cP, cS                      ! P-wave, S-wave speeds
     REAL    (KIND=8)           :: rho, BumpDensity, eta, xi   ! density, boss parameters
  END TYPE HSInfo

  TYPE( HSInfo )                :: HSTop, HSBot

  ! media properties
  INTEGER                       :: Loc(   MaxMedium ), NG( MaxMedium ), N(     MaxMedium )
  REAL      (KIND=8)            :: Depth( MaxMedium ), H(  MaxMedium ), sigma( MaxMedium )
  CHARACTER (LEN= 8)            :: Material( MaxMedium ), TopOpt, BotOpt

  ! storage for finite-difference equations
  CHARACTER (LEN=80)            :: Title
  REAL (KIND=8), ALLOCATABLE    :: B1( : ), B1C( : ), B2( : ), B3( : ), B4( : ), rho( : )

END MODULE krakelmod

