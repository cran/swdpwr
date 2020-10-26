function GEElinear_continous_time(mu, beta, gamma,rho0, II, JJ, KK, &
                     X_in,sigma2,ICC0,ICC2,typeone) result (power)


    implicit none
    ! ---- arg types -----------------------
    integer :: II, JJ, KK,ifault
    integer :: X_in(II,JJ)
    ! II : number of clusters
    ! JJ : number of steps
    ! KK : number of subjects per cluster per step
    double precision :: typeone,ppp
    double precision :: mu, beta,rho0,alpha0,alpha1,alpha2,sigma2,ICC0,ICC2
    ! true values of mu, beta, and tau2
    double precision :: gamma(JJ),B1(JJ*KK,JJ*KK)
    !double precision :: a, b        ! lower and upper limits of GL integral
    !integer :: mincomp(JJ+2), maxcomp(JJ+2)
    !integer :: GQ   ! number of GQ points
    !double precision :: GQX(GQ), GQW(GQ)
    !double precision :: t_power
    double precision :: power,lambda1,lambda2,lambda3,lambda4,vardelta
    ! ----------------------------------------

    double precision, parameter :: PI = 3.14159265358979323846   ! pi
  !  real :: X(JJ*KK,JJ+1)  ! intervention matrix
    !integer :: DD      ! II/(JJ-1). By design, II is a multiple of (JJ-1)
    integer :: i, j, U1,U2,U3,U4,U5
    !logical :: upper

    !double precision, external :: alnorm
    !integer, external :: updatez
    !double precision :: intm4(JJ+1,JJ+1)
    double precision :: z0(JJ*KK),z1(JJ*KK,JJ*KK)   ! z0: # of subjects with outcome 0; z1: # of subjects with outcome 1.
    !double precision :: R(JJ*KK,JJ*KK)
    !double precision :: Ry(JJ+1)
    !double precision :: Wsub(JJ*KK,JJ*KK)
    !double precision :: W(JJ*KK,JJ*KK)
    ! integer :: nZ               ! total number of possible combination of z0 and z1
    !integer :: finish
    !double precision :: derlikelihood(JJ+2)   ! derivative of likelihood
    !double precision :: derlikelihood2(JJ+2, JJ+2)  ! square of derivative of likelihood
    !double precision :: invVar(JJ+2,JJ+2)      ! inverse of variance matrix
    !double precision :: Var(JJ+2,JJ+2)         ! variance matrix
    !double precision :: prob
    !double precision :: sebeta  ! se of beta
    logical :: upper
    double precision, external :: alnorm
    double precision, external :: ppnd16
  !  double precision, external :: diag
    !double precision, external :: matmual
    !double precision, external :: r4_normal_01_cdf_inverse
    upper = .false.
    ifault = 0
    mu=mu
    z0=1
    z1=1

    B1=0
    !allocate(X_in(II,JJ))

    !DD=II/GG
    !DD = II/(JJ-1)   ! II is a multiple of (JJ-1)
    ! assign intervention


    gamma(1)=gamma(1)
    alpha0= ICC0
    alpha1= rho0
    alpha2= ICC2



    !check positive definiteness
    lambda1 = 1- alpha0 + alpha1 - alpha2
    lambda2 = 1- alpha0 - (JJ-1)*(alpha1 - alpha2)
    lambda3 = 1 + (KK-1)*(alpha0 - alpha1) - alpha2
    lambda4 = 1 + (KK-1)*alpha0 + (JJ-1)*(KK-1)*alpha1 + (JJ-1)*alpha2

    if (min(lambda1,lambda2,lambda3,lambda4)<0) then
       !print *, 'not positive-definite, min(lambda1,lambda2,lambda3,lambda4)<0'
      !stop
    end if
    !inverse of R
    !analytical inverse is used here, if the matrix gets large,
    !one could also use the closed-form expression provided in the paper

!calculate variance
   U1=0
   U2=0
   U3=0
   do i=1,II
    ! per cluster
      do j=1,JJ
        U1=U1+X_in(i,j)
      enddo
    enddo

    do j=1,JJ
     ! per cluster
       U4=0
       do i=1,II
        U4=U4+X_in(i,j)
       enddo
     U2=U2+U4*U4
     enddo

     do i=1,II
      ! per cluster
      U5=0
        do j=1,JJ
          U5=U5+X_in(i,j)
        enddo
        U3=U3+U5*U5
      enddo

    ! nZ = 1
    ! do j=1,JJ
    !     nZ = nZ * (KK+1)
    ! end do
    ! run the algorithm
    !invVar = 0.0d0
    !do i=1,II
        ! per cluster
    !    z0 = 0
    !    finish = 0
    !    do while (finish<1)
    !        z1 = KK - z0
            ! prob = prob_z(mu,beta,gamma, tau2, z0, z1, X(i,:), JJ, KK, GQ, GQX, GQW)
    !        call der_likelihood_time(mu,beta,gamma,tau2, z0, z1, X(i,:), JJ, KK, a, b, &
    !                            mincomp, maxcomp, GQ, GQX, GQW, derlikelihood, prob)
    !        call vectorsquare(derlikelihood, JJ+2, derlikelihood2)
    !        invVar = invVar + derlikelihood2 * prob
    !        finish = updatez(z0, JJ, KK)
    !    end do
    !end do
    !call syminverse(invVar,Var,JJ+2)
    !sebeta = sqrt(Var(2,2))
    !power = alnorm(beta/sebeta-1.959964d0,upper) + alnorm(-beta/sebeta-1.959964d0,upper)


    !delta=beta
    vardelta = (sigma2/KK)*II*JJ*lambda3*lambda4/((U1*U1+II*JJ*U1-JJ*U2-II*U3)*lambda4-(U1*U1-II*U3)*lambda3)


     ! verify whether this is close to the assumed delta, which is log(0.25)
    !vardelta = invOmega(1,1)!  # var of the trt effect estimator


    !delta=0
    !vardelta=1
    ! z-test
    !power = alnorm(delta/sqrt(vardelta)-1.959964d0,upper) + alnorm(-delta/sqrt(vardelta)-1.959964d0,upper)
    !if(beta<0) then
    !power =alnorm(-beta/sqrt(vardelta)-1.959964d0,upper)
  !else
    !power =alnorm(beta/sqrt(vardelta)-1.959964d0,upper)
  !end if
!  temerror=r4_normal_01_cdf_inverse(0.95)
  !temerror=ppnd(0.95)
   ppp=ppnd16(typeone/2,ifault)
    power =alnorm(-beta/sqrt(vardelta)+ppp,upper)+alnorm(beta/sqrt(vardelta)+ppp,upper)

  !power =alnorm(-beta/sqrt(vardelta)-1.96,upper)+alnorm(beta/sqrt(vardelta)-1.96,upper)
    !z_power = pnorm(qnorm(0.05/2) + abs(delta)/sqrt(vardelta))
   !print(z_power)
  !power=temerror
    ! t-test df=# clusters minus p
  !  df = II - (JJ+1)
  !  t_power = pt(qt(0.05/2,df=df) + abs(delta)/sqrt(vardelta),df=df)


end function GEElinear_continous_time
