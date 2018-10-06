SUBROUTINE ReadHeader( SHDFile, FileName, Title, Atten, PlotType )

  ! Read header from disk file
  ! This routine is not currently used anywhere

  ! FileName is a SHDFIL for complex pressure or a GRNFIL for a Green's function
  ! Title   arbitrary title

  ! variables taken from SdRdRMod:
  ! FreqVec vector of frequencies
  ! theta   vector of bearing lines,   theta( 1 : Ntheta )
  ! sd      vector of source   depths, sd(    1 : Nsd    )
  ! rd      vector of receiver depths, rd(    1 : Nrd    )
  ! r       vector of receiver ranges, r(     1 : Nr     )

  USE SdRdRMod
  IMPLICIT NONE
  INTEGER, PARAMETER                :: PRTFile = 6
  REAL,               INTENT( OUT ) :: Atten           ! stabilizing attenuation for SCOOTER FFP runs
  INTEGER                           :: SHDFile         ! unit number of SHDFile
  CHARACTER (LEN=80), INTENT( OUT ) :: Title, FileName
  CHARACTER (LEN=10), INTENT( OUT ) :: PlotType
  INTEGER                           :: IAllocStat, IOStat, LRecL

  ! Open file, read header
  IF ( SHDFile == 0 ) SHDFile = 25
  IF ( FileName( 1 : 1 ) == ' ' ) FileName = 'SHDFIL'

  ! INQUIRE( FILE = FileName, RECL = IRECL )
  OPEN( UNIT = SHDFile,   FILE = FileName, STATUS = 'OLD', ACCESS = 'DIRECT', FORM = 'UNFORMATTED', RECL = 4, &
        IOSTAT = IOStaT, ACTION = 'READ' )
  IF ( IOStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadHeader', 'Unable to open shade file' )
  
  READ( SHDFile, REC = 1 ) LRecl
  CLOSE( UNIT = SHDFile )
  OPEN(  UNIT = SHDFile,   FILE = FileName, STATUS = 'OLD', ACCESS = 'DIRECT', FORM = 'UNFORMATTED', RECL = 4 * LRecl )

  READ( SHDFile, REC = 1 ) LRecl, Title
  READ( SHDFile, REC = 2 ) PlotType
  READ( SHDFile, REC = 3 ) Nfreq, Pos%Ntheta, Pos%Nsx, Pos%Nsy, Pos%Nsd, Pos%Nrd, Pos%Nr, atten

  ALLOCATE( FreqVec( Nfreq ), Pos%sd( Pos%Nsd ), Pos%rd( Pos%Nrd ), Pos%r( Pos%Nr ), Pos%theta( Pos%Ntheta ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadHeader', 'Too many source/receiver combinations' )

  READ( SHDFile, REC = 4 ) FreqVec
  READ( SHDFile, REC = 5 ) Pos%theta
  READ( SHDFile, REC = 6 ) Pos%sx
  READ( SHDFile, REC = 7 ) Pos%sy
  READ( SHDFile, REC = 8 ) Pos%sd
  READ( SHDFile, REC = 9 ) Pos%rd
  READ( SHDFile, REC = 10 ) Pos%r

  ! Pos%deltaR = Pos%r( Pos%Nr ) - Pos%r( Pos%Nr - 1 )

END SUBROUTINE ReadHeader
