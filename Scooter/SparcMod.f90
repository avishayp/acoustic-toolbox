MODULE SparcMod

  USE MathConstantsMod
  SAVE
  INTEGER,       PARAMETER :: ENVFile = 5, PRTFile = 6, GRNFile = 25, RTSFile = 35, MaxN = 17000, MaxMedium = 500, MaxIt = 1000000
  INTEGER              :: Loc( MaxMedium), N( MaxMedium ), NSig, Nk, NTot1, Nrr, NTout
  REAL                 :: C2R( MaxN ), C2I( MaxN ), Z( MaxN ), &
                          H( MaxMedium ), CMin, CLow, CHigh, omega2, DeltaK, &
                          Deltat, CrossT, CMax, TStart, V, TMult, alpha, beta, FMin, FMax
  REAL        (KIND=8) :: rho( MaxN )
  CHARACTER (LEN=8 )   :: TopOpt, BotOpt
  CHARACTER (LEN=4 )   :: Pulse
  CHARACTER (LEN=80)   :: Title
  REAL,    ALLOCATABLE :: k( : ), Tout( : ), RTSrd( :, : ), RTSrr( :, : )
  COMPLEX, ALLOCATABLE :: Green( :, :, : )

END MODULE SparcMod
