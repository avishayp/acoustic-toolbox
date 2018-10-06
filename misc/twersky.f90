FUNCTION Twersky( Option, omega, HS, kx, rho0, c0 )

  ! Dummy version

  IMPLICIT NONE
  REAL     (KIND=8) :: omega, rho0, c0
  COMPLEX  (KIND=8) :: Twersky, kx
  CHARACTER (LEN=1) :: Option

  ! Halfspace properties
  TYPE HSInfo
     CHARACTER (LEN=1)          :: BC                          ! Boundary condition type
     COMPLEX (KIND=8)           :: cP, cS                      ! P-wave, S-wave speeds
     REAL    (KIND=8)           :: rho, BumpDensity, eta, xi   ! density, boss parameters
  END TYPE HSInfo

  TYPE( HSInfo )         :: HS

  Twersky = 0.0

END FUNCTION Twersky
