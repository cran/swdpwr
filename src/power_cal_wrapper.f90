SUBROUTINE LinearPowertimewrapper(mu, beta, gamma, tau2, II, JJ, KK, a, b, mincomp, maxcomp, GQ, GQX, GQW, X_in,typeone,power)

  integer :: II, JJ, KK
  integer :: X_in(II,JJ)
  ! II : number of clusters
  ! JJ : number of steps
  ! KK : number of subjects per cluster per step
  double precision :: mu, beta, tau2,typeone
  ! true values of mu, beta, and tau2
  double precision :: gamma(JJ)
  double precision :: a, b        ! lower and upper limits of GL integral
  integer :: mincomp(JJ+2), maxcomp(JJ+2)
  integer :: GQ   ! number of GQ points
  double precision :: GQX(GQ), GQW(GQ)
  double precision :: power

double precision, external :: LinearPower_time

power = LinearPower_time(mu, beta, gamma, tau2, II, JJ, KK, a, b, mincomp, maxcomp, GQ, GQX, GQW, X_in,typeone)

END SUBROUTINE

SUBROUTINE LogPowertimewrapper(mu, beta, gamma, tau2, II, JJ, KK, a, b, mincomp, maxcomp, GQ, GQX, GQW, X_in,typeone, power)

  integer :: II, JJ, KK
  integer :: X_in(II,JJ)
  ! II : number of clusters
  ! JJ : number of steps
  ! KK : number of subjects per cluster per step
  double precision :: typeone,mu, beta, tau2
  ! true values of mu, beta, and tau2
  double precision :: gamma(JJ)
  double precision :: a, b        ! lower and upper limits of GL integral
  integer :: mincomp(JJ+2), maxcomp(JJ+2)
  integer :: GQ   ! number of GQ points
  double precision :: GQX(GQ), GQW(GQ)
  double precision :: power

double precision, external :: LogPower_time

power = LogPower_time(mu, beta, gamma, tau2, II, JJ, KK, a, b, mincomp, maxcomp, GQ, GQX, GQW, X_in,typeone)

END SUBROUTINE

SUBROUTINE LogitPowertimewrapper(mu, beta, gamma, tau2, II, JJ, KK, GQ, GQX, GQW, X_in,typeone, power)

  integer :: II, JJ, KK
  integer :: X_in(II,JJ)
  ! II : number of clusters
  ! JJ : number of steps
  ! KK : number of subjects per cluster per step
  double precision :: typeone,mu, beta, tau2
  ! true values of mu, beta, and tau2
  double precision :: gamma(JJ)
  !double precision :: a, b        ! lower and upper limits of GL integral
  !integer :: mincomp(JJ+2), maxcomp(JJ+2)
  integer :: GQ   ! number of GQ points
  double precision :: GQX(GQ), GQW(GQ)
  double precision :: power

double precision, external :: LogitPower_time

power = LogitPower_time(mu, beta, gamma, tau2, II, JJ, KK, GQ, GQX, GQW, X_in,typeone)
END SUBROUTINE

SUBROUTINE LinearPowernotimewrapper(mu, beta, tau2, II, JJ, KK, a, b, GQ, GQX, GQW, X_in,typeone, power)

  integer :: II, JJ, KK
  integer :: X_in(II,JJ)
  ! II : number of clusters
  ! JJ : number of steps
  ! KK : number of subjects per cluster per step
  double precision :: typeone,mu, beta, tau2
  ! true values of mu, beta, and tau2
  !double precision :: gamma(JJ)
  double precision :: a, b        ! lower and upper limits of GL integral
  !integer :: mincomp(JJ+2), maxcomp(JJ+2)
  integer :: GQ   ! number of GQ points
  double precision :: GQX(GQ), GQW(GQ)
  double precision :: power

double precision, external :: LinearPower_notime

power = LinearPower_notime(mu, beta,tau2, II, JJ, KK, a, b, GQ, GQX, GQW, X_in,typeone)

END SUBROUTINE

SUBROUTINE LogPowernotimewrapper(mu, beta, tau2, II, JJ, KK, a, b, GQ, GQX, GQW, X_in,typeone, power)

  integer :: II, JJ, KK
  integer :: X_in(II,JJ)
  ! II : number of clusters
  ! JJ : number of steps
  ! KK : number of subjects per cluster per step
  double precision :: typeone,mu, beta, tau2
  ! true values of mu, beta, and tau2
  !double precision :: gamma(JJ)
  double precision :: a, b        ! lower and upper limits of GL integral
  !integer :: mincomp(JJ+2), maxcomp(JJ+2)
  integer :: GQ   ! number of GQ points
  double precision :: GQX(GQ), GQW(GQ)
  double precision :: power

double precision, external :: LogPower_notime

power = LogPower_notime(mu, beta, tau2, II, JJ, KK, a, b, GQ, GQX, GQW, X_in,typeone)

END SUBROUTINE

SUBROUTINE LogitPowernotimewrapper(mu, beta, tau2, II, JJ, KK, GQ, GQX, GQW, X_in,typeone, power)

  integer :: II, JJ, KK
  integer :: X_in(II,JJ)
  ! II : number of clusters
  ! JJ : number of steps
  ! KK : number of subjects per cluster per step
  double precision :: typeone,mu, beta, tau2
  ! true values of mu, beta, and tau2
  !double precision :: gamma(JJ)
  !double precision :: a, b        ! lower and upper limits of GL integral
!  integer :: mincomp(JJ+2), maxcomp(JJ+2)
  integer :: GQ   ! number of GQ points
  double precision :: GQX(GQ), GQW(GQ)
  double precision :: power

double precision, external :: LogitPower_notime

power = LogitPower_notime(mu, beta,tau2, II, JJ, KK, GQ, GQX, GQW, X_in,typeone)
END SUBROUTINE

SUBROUTINE LinearPowerGEEwrapper(mu, beta, gamma,rho0, II, JJ, KK, &
                     X_in,ICC0,ICC2,typeone,power)

  integer :: II, JJ, KK
  integer :: X_in(II,JJ)
  ! II : number of clusters
  ! JJ : number of steps
  ! KK : number of subjects per cluster per step
  double precision :: typeone,mu, beta, rho0,ICC0,ICC2
  ! true values of mu, beta, and tau2
  double precision :: gamma(JJ)
  !double precision :: a, b        ! lower and upper limits of GL integral
  !integer :: mincomp(JJ+2), maxcomp(JJ+2)
  double precision :: power

double precision, external ::  LinearPower_GEE

power = LinearPower_GEE(mu, beta, gamma,rho0, II, JJ, KK, X_in,ICC0,ICC2,typeone)
END SUBROUTINE

SUBROUTINE LogPowerGEEwrapper(mu, beta, gamma,rho0, II, JJ, KK, &
                     X_in,ICC0,ICC2,typeone,power)

  integer :: II, JJ, KK
  integer :: X_in(II,JJ)
  ! II : number of clusters
  ! JJ : number of steps
  ! KK : number of subjects per cluster per step
  double precision :: typeone,mu, beta, rho0,ICC0,ICC2
  ! true values of mu, beta, and tau2
  double precision :: gamma(JJ)
  !double precision :: a, b        ! lower and upper limits of GL integral
  !integer :: mincomp(JJ+2), maxcomp(JJ+2)
  double precision :: power

double precision, external ::  LogPower_GEE

power = LogPower_GEE(mu, beta, gamma,rho0, II, JJ, KK, X_in,ICC0,ICC2,typeone)
END SUBROUTINE

SUBROUTINE LogitPowerGEEwrapper(mu, beta, gamma,rho0, II, JJ, KK, &
                     X_in,ICC0,ICC2,typeone,power)

  integer :: II, JJ, KK
  integer :: X_in(II,JJ)
  ! II : number of clusters
  ! JJ : number of steps
  ! KK : number of subjects per cluster per step
  double precision :: typeone,mu, beta, rho0,ICC0,ICC2
  ! true values of mu, beta, and tau2
  double precision :: gamma(JJ)
  !double precision :: a, b        ! lower and upper limits of GL integral
  !integer :: mincomp(JJ+2), maxcomp(JJ+2)
  double precision :: power

double precision, external ::  LogitPower_GEE

power = LogitPower_GEE(mu, beta, gamma,rho0, II, JJ, KK, X_in,ICC0,ICC2,typeone)
END SUBROUTINE

SUBROUTINE LinearPowerGEEnotimewrapper(mu, beta, gamma,rho0, II, JJ, KK, &
                     X_in,ICC0,ICC2,typeone,power)

  integer :: II, JJ, KK
  integer :: X_in(II,JJ)
  ! II : number of clusters
  ! JJ : number of steps
  ! KK : number of subjects per cluster per step
  double precision :: typeone,mu, beta, rho0,ICC0,ICC2
  ! true values of mu, beta, and tau2
  double precision :: gamma(JJ)
  !double precision :: a, b        ! lower and upper limits of GL integral
  !integer :: mincomp(JJ+2), maxcomp(JJ+2)
  double precision :: power

double precision, external ::  LinearPower_GEE_notime

power = LinearPower_GEE_notime(mu, beta, gamma,rho0, II, JJ, KK, X_in,ICC0,ICC2,typeone)
END SUBROUTINE

SUBROUTINE LogPowerGEEnotimewrapper(mu, beta, gamma,rho0, II, JJ, KK, &
                     X_in,ICC0,ICC2,typeone,power)

  integer :: II, JJ, KK
  integer :: X_in(II,JJ)
  ! II : number of clusters
  ! JJ : number of steps
  ! KK : number of subjects per cluster per step
  double precision :: typeone,mu, beta, rho0,ICC0,ICC2
  ! true values of mu, beta, and tau2
  double precision :: gamma(JJ)
  !double precision :: a, b        ! lower and upper limits of GL integral
  !integer :: mincomp(JJ+2), maxcomp(JJ+2)
  double precision :: power

double precision, external ::  LogPower_GEE_notime

power = LogPower_GEE_notime(mu, beta, gamma,rho0, II, JJ, KK, X_in,ICC0,ICC2,typeone)
END SUBROUTINE

SUBROUTINE LogitPowerGEEnotimewrapper(mu, beta, gamma,rho0, II, JJ, KK, &
                     X_in,ICC0,ICC2,typeone,power)

  integer :: II, JJ, KK
  integer :: X_in(II,JJ)
  ! II : number of clusters
  ! JJ : number of steps
  ! KK : number of subjects per cluster per step
  double precision :: typeone,mu, beta, rho0,ICC0,ICC2
  ! true values of mu, beta, and tau2
  double precision :: gamma(JJ)
  !double precision :: a, b        ! lower and upper limits of GL integral
!  integer :: mincomp(JJ+2), maxcomp(JJ+2)
  double precision :: power

double precision, external ::  LogitPower_GEE_notime

power = LogitPower_GEE_notime(mu, beta, gamma,rho0, II, JJ, KK, X_in,ICC0,ICC2,typeone)

END SUBROUTINE

SUBROUTINE ContinuousPowerGEEnotimewrapper(mu, beta, gamma,rho0, II, JJ, KK, &
                     X_in,sigma2,ICC0,ICC2,typeone,power)

  integer :: II, JJ, KK
  integer :: X_in(II,JJ)
  ! II : number of clusters
  ! JJ : number of steps
  ! KK : number of subjects per cluster per step
  double precision ::typeone, mu, beta, sigma2,ICC0,ICC2,rho0
  ! true values of mu, beta, and tau2
  double precision :: gamma(JJ)
  !double precision :: a, b        ! lower and upper limits of GL integral
  !integer :: mincomp(JJ+2), maxcomp(JJ+2)
  double precision :: power

double precision, external ::  GEElinear_continous_notime

power = GEElinear_continous_notime(mu, beta, gamma,rho0, II, JJ, KK, &
                     X_in,sigma2,ICC0,ICC2,typeone)

END SUBROUTINE

SUBROUTINE ContinuousPowerGEEtimewrapper(mu, beta, gamma,rho0, II, JJ, KK, &
                     X_in,sigma2,ICC0,ICC2,typeone,power)

  integer :: II, JJ, KK
  integer :: X_in(II,JJ)
  ! II : number of clusters
  ! JJ : number of steps
  ! KK : number of subjects per cluster per step
  double precision :: typeone,mu, beta, sigma2,ICC0,ICC2,rho0
  ! true values of mu, beta, and tau2
  double precision :: gamma(JJ)
  !double precision :: a, b        ! lower and upper limits of GL integral
  !integer :: mincomp(JJ+2), maxcomp(JJ+2)
  double precision :: power

double precision, external :: GEElinear_continous_time

power = GEElinear_continous_time(mu, beta, gamma,rho0, II, JJ, KK, &
                     X_in,sigma2,ICC0,ICC2,typeone)
END SUBROUTINE
