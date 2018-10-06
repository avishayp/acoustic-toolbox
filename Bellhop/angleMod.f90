MODULE anglemod

  USE MathConstantsMod
  USE SubTabulate
  USE SdRdRMod
  USE SortMod

  IMPLICIT NONE
  SAVE

  INTEGER          :: ialpha, ibeta
  INTEGER, PRIVATE :: AllocateStatus
  INTEGER,       PRIVATE, PARAMETER :: ENVFile = 5, PRTFile = 6
  REAL (KIND=8), PRIVATE, PARAMETER :: c0 = 1500.0

  TYPE AnglesStructure
    INTEGER       :: Nalpha = 0, Nbeta = 1, iSingle = 0, iSingle2 = 0
    REAL (KIND=8) :: Dalpha, Dbeta
    REAL (KIND=8), ALLOCATABLE:: alpha( : ), beta( : )
  END TYPE AnglesStructure

  Type( AnglesStructure ) :: Angles

CONTAINS
  SUBROUTINE ReadRayElevationAngles( freq, Depth, TopOpt, RunType )

    ! *** Beam elevation angles ***

    REAL      (KIND=8), INTENT( IN  ) :: freq, Depth
    CHARACTER (LEN= 6), INTENT( IN  ) :: TopOpt, RunType
    REAL      (KIND=8)                :: d_theta_recommended

    IF ( TopOpt( 6 : 6 ) == 'I' ) THEN
       READ( ENVFile, * ) Angles%Nalpha, Angles%iSingle ! option to trace a single beam
    ELSE
       READ( ENVFile, * ) Angles%Nalpha
    END IF

    IF ( Angles%Nalpha == 0 ) THEN   ! automatically estimate Nalpha to use
       IF ( RunType( 1 : 1 ) == 'R' ) THEN
          Angles%Nalpha = 50   ! For a ray trace plot, we don't want too many rays ...
       ELSE
          ! you're letting ME choose? OK: ideas based on an isospeed ocean
          ! limit based on phase of adjacent beams at maximum range
          Angles%Nalpha = MAX( INT( 0.3 * Pos%r( Pos%Nr ) * Freq / c0 ), 300 )

          ! limit based on having beams that are thin with respect to the water depth
          ! assumes also a full 360 degree angular spread of rays
          d_theta_recommended = ATAN( Depth / ( 10.0 * Pos%r( Pos%Nr ) ) )
          Angles%Nalpha = MAX( INT( pi / d_theta_recommended ), Angles%Nalpha )
       END IF
    END IF

    ALLOCATE( Angles%alpha( Angles%Nalpha ), STAT = AllocateStatus )
    IF ( AllocateStatus /= 0 ) CALL ERROUT( PRTFile, 'F', 'READIN', 'Insufficient memory to store beam angles'  )

    IF ( Angles%Nalpha > 2 ) Angles%alpha( 3 ) = -999.9
    READ( ENVFile, * ) Angles%alpha
    CALL SubTab( Angles%alpha, Angles%Nalpha )
    CALL Sort(   Angles%alpha, Angles%Nalpha )

    ! full 360-degree sweep? remove duplicate beam
    IF ( Angles%alpha( Angles%Nalpha ) == Angles%alpha( 1 ) + 360.0 ) Angles%Nalpha = Angles%Nalpha - 1

    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) 'Number of beams in elevation   = ', Angles%Nalpha
    IF ( Angles%iSingle > 0 ) WRITE( PRTFile, * ) 'Trace only beam number ', Angles%iSingle
    WRITE( PRTFile, * ) 'Beam take-off angles (degrees)'

    IF ( Angles%Nalpha >= 1 ) WRITE( PRTFile, "( 5G14.6 )" ) Angles%alpha( 1 : MIN( Angles%Nalpha, Number_to_Echo ) )
    IF ( Angles%Nalpha > Number_to_Echo ) WRITE( PRTFile, * ) ' ... ', Angles%alpha( Angles%Nalpha )
 
    IF ( Angles%Nalpha > 1 .AND. Angles%alpha( Angles%Nalpha ) == Angles%alpha( 1 ) ) &
         CALL ERROUT( PRTFile, 'F', 'BELLHOP: READIN', 'First and last beam take-off angle are identical' )

    IF ( TopOpt( 6 : 6 ) == 'I' ) THEN
       IF ( Angles%iSingle < 1 .OR. Angles%iSingle > Angles%Nalpha ) &
            CALL ERROUT( PRTFile, 'F', 'BELLHOP: READIN', 'Selected beam, iSingl not in [ 1, Angles%Nalpha ]' )
    END IF

  END SUBROUTINE ReadRayElevationAngles

  !**********************************************************************!

  SUBROUTINE ReadRayBearingAngles( freq, ThreeD, TopOpt, RunType )

    ! *** Beam bearing angles ***

    LOGICAL,            INTENT( IN ) :: ThreeD
    REAL      (KIND=8), INTENT( IN ) :: freq
    CHARACTER (LEN= 6), INTENT( IN ) :: TopOpt, RunType

    IF ( ThreeD ) THEN
       IF ( TopOpt( 6 : 6 ) == 'I' ) THEN
          READ( ENVFile, * ) Angles%Nbeta, Angles%iSingle2 ! option to trace a single beam
       ELSE
          READ( ENVFile, * ) Angles%Nbeta
       END IF

       IF ( Angles%Nbeta == 0 ) THEN   ! automatically estimate Nalpha to use
          IF ( RunType( 1 : 1 ) == 'R' ) THEN
             Angles%Nbeta = 50   ! For a ray trace plot, we don't want too many rays ...
          ELSE
             Angles%Nbeta = MAX( INT( 0.1 * Pos%r( Pos%Nr ) * Freq / c0 ), 300 )
          END IF
       END IF

       ALLOCATE( Angles%beta( Angles%Nbeta ), STAT = AllocateStatus )
       IF ( AllocateStatus /= 0 ) CALL ERROUT( PRTFile, 'F', 'READIN', 'Insufficient memory to store beam angles'  )

       IF ( Angles%Nbeta > 2 ) Angles%beta( 3 ) = -999.9
       READ( ENVFile, * ) Angles%beta
       CALL SubTab( Angles%beta, Angles%Nbeta )
       CALL Sort(   Angles%beta, Angles%Nbeta )

       ! full 360-degree sweep? remove duplicate beam
       IF ( Angles%beta( Angles%Nbeta ) == Angles%beta( 1 ) + 360.0D0 ) Angles%Nbeta = Angles%Nbeta - 1

       ! Nx2D CASE: beams must lie on rcvr radials--- replace beta with theta
       IF ( RunType( 6 : 6 ) == '2' .AND. RunType( 1 : 1 ) /= 'R' ) THEN
          WRITE( PRTFile, * )
          WRITE( PRTFile, * ) 'Replacing beam take-off angles, beta, with receiver bearing lines, theta'
          DEALLOCATE( Angles%beta )

          Angles%Nbeta = Pos%Ntheta
          ALLOCATE( Angles%beta( Angles%Nbeta ), STAT = AllocateStatus )
          IF ( AllocateStatus /= 0 ) CALL ERROUT( PRTFile, 'F', 'READIN', 'Insufficient memory to store beam angles'  )
          Angles%beta( 1 : Angles%Nbeta ) = Pos%theta( 1 : Pos%Ntheta )   ! Nbeta should = Ntheta
       END IF

       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Number of beams in bearing   = ', Angles%Nbeta
       IF ( Angles%iSingle2 > 0 ) WRITE( PRTFile, * ) 'Trace only beam number ', Angles%iSingle2
       WRITE( PRTFile, * ) 'Beam take-off angles (degrees)'

       IF ( Angles%Nbeta >= 1 ) WRITE( PRTFile, "( 5G14.6 )" ) Angles%beta( 1 : MIN( Angles%Nbeta, Number_to_Echo ) )
       IF ( Angles%Nbeta > Number_to_Echo ) WRITE( PRTFile, * ) ' ... ', Angles%beta( Angles%Nbeta )

       IF ( Angles%Nbeta > 1 .AND. Angles%beta( Angles%Nbeta ) == Angles%beta( 1 ) ) &
            CALL ERROUT( PRTFile, 'F', 'BELLHOP: READIN', 'First and last beam take-off angle are identical' )

       IF ( TopOpt( 6 : 6 ) == 'I' ) THEN
          IF ( Angles%iSingle2 < 1 .OR. Angles%iSingle2 > Angles%Nbeta ) &
               CALL ERROUT( PRTFile, 'F', 'BELLHOP: READIN', 'Selected beam, iSingl not in [ 1, Angles%Nbeta ]' )
       END IF
       Angles%beta  = DegRad * Angles%beta   ! convert to radians

       Angles%Dbeta = 0.0
       IF ( Angles%Nbeta /= 1 ) Angles%Dbeta = ( Angles%beta( Angles%NBeta ) - Angles%beta( 1 ) ) / ( Angles%Nbeta - 1 )

    END IF

  END SUBROUTINE ReadRayBearingAngles

END MODULE anglemod
