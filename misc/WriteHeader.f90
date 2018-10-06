SUBROUTINE WriteHeader( FileName, Title, Atten, PlotType )

  ! Write header to disk file

  USE SdRdRMod
  IMPLICIT NONE
  INTEGER, PARAMETER      :: SHDFile = 25
  REAL,      INTENT( IN ) :: Atten             ! stabilizing attenuation (for wavenumber integration only)
  CHARACTER, INTENT( IN ) :: FileName*( * )    ! Name of the file (could be a shade file or a Green's function file)
  CHARACTER, INTENT( IN ) :: Title*( * )       ! Arbitrary title
  CHARACTER, INTENT( IN ) :: PlotType*( 10 )   ! 

  INTEGER LRecl

  ! receiver bearing angles
  IF ( .NOT. ALLOCATED( Pos%theta ) ) THEN
     ALLOCATE( Pos%theta( 1 ) )
     Pos%theta( 1 ) = 0   ! dummy bearing angle
     Pos%Ntheta     = 1
  END IF

  ! source x-coordinates
  IF ( .NOT. ALLOCATED( Pos%sx ) ) THEN
     ALLOCATE( Pos%sx( 1 ) )
     Pos%sx( 1 ) = 0   ! dummy x-coordinate
     Pos%Nsx     = 1
  END IF

  ! source y-coordinates
  IF ( .NOT. ALLOCATED( Pos%sy ) ) THEN
     ALLOCATE( Pos%sy( 1 ) )
     Pos%sy( 1 ) = 0   ! dummy y-coordinate
     Pos%Nsy     = 1
  END IF

  IF ( PlotType( 1 : 2 ) /= 'TL' ) THEN
     ! MAX( 41, ... ) below because Title is already 40 words (or 80 bytes)
     LRecl = MAX( 41, 2 * Nfreq, Pos%Ntheta, Pos%Nsx, Pos%Nsy, Pos%Nsd, Pos%Nrd, 2 * Pos%Nr )   ! words/record (Nr doubled for complex pressure storage)

     OPEN ( FILE = FileName, UNIT = SHDFile, STATUS = 'UNKNOWN', ACCESS = 'DIRECT', RECL = 4 * LRecl, FORM = 'UNFORMATTED')
     WRITE( SHDFile, REC = 1 ) LRecl, Title( 1 : 80 )
     WRITE( SHDFile, REC = 2 ) PlotType
     WRITE( SHDFile, REC = 3 ) Nfreq, Pos%Ntheta, Pos%Nsx, Pos%Nsy, Pos%Nsd, Pos%Nrd, Pos%Nr, atten
     WRITE( SHDFile, REC = 4 ) freqVec( 1 : Nfreq )
     WRITE( SHDFile, REC = 5 ) Pos%theta( 1 : Pos%Ntheta )
     WRITE( SHDFile, REC = 6 ) Pos%sx( 1 : Pos%Nsx )
     WRITE( SHDFile, REC = 7 ) Pos%sy( 1 : Pos%Nsy )
     WRITE( SHDFile, REC = 8 ) Pos%sd( 1 : Pos%Nsd )
     WRITE( SHDFile, REC = 9 ) Pos%rd( 1 : Pos%Nrd )
     WRITE( SHDFile, REC = 10 ) Pos%r( 1 : Pos%Nr  )
  ELSE   ! compressed format for TL from FIELD3D
     LRecl = MAX( 41, 2 * Nfreq, Pos%Ntheta, Pos%Nsd, Pos%Nrd, 2 * Pos%Nr )   ! words/record (Nr doubled for complex pressure storage)

     OPEN ( FILE = FileName, UNIT = SHDFile, STATUS = 'UNKNOWN', ACCESS = 'DIRECT', RECL = 4 * LRecl, FORM = 'UNFORMATTED')
     WRITE( SHDFile, REC = 1 ) LRecl, Title( 1 : 80 )
     WRITE( SHDFile, REC = 2 ) PlotType
     WRITE( SHDFile, REC = 3 ) Nfreq, Pos%Ntheta, Pos%Nsx, Pos%Nsy, Pos%Nsd, Pos%Nrd, Pos%Nr, atten
     WRITE( SHDFile, REC = 4 ) freqVec( 1 : Nfreq )
     WRITE( SHDFile, REC = 5 ) Pos%theta( 1 : Pos%Ntheta )
     WRITE( SHDFile, REC = 6 ) Pos%sx( 1 ), Pos%sx( Pos%Nsx )
     WRITE( SHDFile, REC = 7 ) Pos%sy( 1 ), Pos%sy( Pos%Nsy )
     WRITE( SHDFile, REC = 8 ) Pos%sd( 1 : Pos%Nsd )
     WRITE( SHDFile, REC = 9 ) Pos%rd( 1 : Pos%Nrd )
     WRITE( SHDFile, REC = 10 ) Pos%r( 1 : Pos%Nr  )
  END IF

END SUBROUTINE WriteHeader

!**********************************************************************!

SUBROUTINE WriteField( P, Nrd, Nr, IRec )

  ! Write the field to disk

  IMPLICIT NONE
  INTEGER, PARAMETER       :: SHDFile = 25
  INTEGER, INTENT( IN )    :: Nrd, Nr        ! Number of receiver depths, ranges
  COMPLEX, INTENT( IN )    :: P( Nrd, Nr )   ! Pressure field
  INTEGER, INTENT( INOUT ) :: IRec           ! last record read
  INTEGER                  :: ird

  DO ird = 1, Nrd
     IRec = IRec + 1
     WRITE( SHDFile, REC = IRec ) P( ird, : )
  END DO

END SUBROUTINE WriteField
