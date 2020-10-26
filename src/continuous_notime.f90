function GEElinear_continous_notime(mu, beta, gamma,rho0, II, JJ, KK, &
                     X_in,sigma2,ICC0,ICC2,typeone) result (power)

    implicit none
    ! ---- arg types -----------------------
    integer :: II, JJ, KK
    integer :: X_in(II,JJ)
    ! II : number of clusters
    ! JJ : number of steps
    ! KK : number of subjects per cluster per step
    double precision :: mu, beta,rho0,alpha0,alpha1,alpha2,sigma2,ICC0,ICC2
    ! true values of mu, beta, and tau2
    double precision :: gamma(JJ),B1(JJ*KK,JJ*KK)
    !double precision :: a, b        ! lower and upper limits of GL integral
    !integer :: mincomp(JJ+2), maxcomp(JJ+2)
    !integer :: GQ   ! number of GQ points
    !double precision :: GQX(GQ), GQW(GQ)
    !double precision :: t_power
    double precision :: power,lambda1,lambda2,lambda3,lambda4,vardelta,typeone,ppp
    ! ----------------------------------------

    double precision, parameter :: PI = 3.14159265358979323846   ! pi
    real :: X(JJ*KK,2)  ! intervention matrix
    !integer :: DD      ! II/(JJ-1). By design, II is a multiple of (JJ-1)
    integer :: i, j, k, ifault
    !logical :: upper
    double precision :: B2(JJ*KK,JJ*KK),B3(JJ*KK,JJ*KK)
    !double precision, external :: alnorm
    !integer, external :: updatez
    double precision :: M1(JJ,JJ),M2(KK,KK), M3(JJ,JJ), M4(KK,KK)
    double precision :: z0(JJ*KK),z1(JJ*KK,JJ*KK)   ! z0: # of subjects with outcome 0; z1: # of subjects with outcome 1.
    double precision :: bm2(JJ*KK,JJ*KK),bm3(JJ*KK,JJ*KK)
    double precision :: bm1(JJ*KK,JJ*KK),bm4(JJ*KK,JJ*KK),R(JJ*KK,JJ*KK)
    double precision :: invR(JJ*KK,JJ*KK),Omega(2,2)
    double precision :: invOmega(2,2)
    double precision :: New1(2,JJ*KK),New2(2,2)
    ! integer :: nZ               ! total number of possible combination of z0 and z1
    !integer :: finish
    !double precision :: derlikelihood(JJ+2)   ! derivative of likelihood
    !double precision :: derlikelihood2(JJ+2, JJ+2)  ! square of derivative of likelihood
    !double precision :: invVar(JJ+2,JJ+2)      ! inverse of variance matrix
    !double precision :: Var(JJ+2,JJ+2)         ! variance matrix
    !double precision :: prob
    !double precision :: sebeta  ! se of beta
    logical :: upper
    integer :: l
    double precision, external :: alnorm
    double precision, external :: ppnd16
  !  double precision, external :: diag
    double precision, external :: matmual
    upper = .false.
    ifault = 0

    gamma=0
    mu=mu


    z0=1
    z1=1

    B1=0



    !allocate(X_in(II,JJ))

    !DD=II/GG
    !DD = II/(JJ-1)   ! II is a multiple of (JJ-1)
    ! assign intervention



    alpha0= ICC0
    alpha1= rho0
    alpha2= ICC2

    !correlation matrix
    do i = 1, JJ*KK
    B1(i,i) = z0(i)
    end do

    M1=1
    M2=0
    do i = 1,KK
      M2(i,i)=1
    end do

    M3=0
    do i = 1,JJ
      M3(i,i)=1
    end do
    M4=1

    B2=0
    B3=0
    bm1 = (1- alpha0 + alpha1 - alpha2) * B1

    Do i=1,JJ  ! tab(i,.)
          Do j=1,JJ  ! tab(i,j)
                Do k=1,KK
                      Do l=1,KK
                             B2((i-1)*KK+k,(j-1)*KK+l)=M1(i,j)*M2(k,l)
                       Enddo
                 Enddo
           Enddo
    enddo






    !call kron_prod(M1,M2,B2)
    bm2 = (alpha2 - alpha1)*B2


        Do i=1,JJ  ! tab(i,.)
              Do j=1,JJ  ! tab(i,j)
                    Do k=1,KK
                          Do l=1,KK
                                 B3((i-1)*KK+k,(j-1)*KK+l)=M3(i,j)*M4(k,l)
                           Enddo
                     Enddo
               Enddo
        enddo



    !call kron_prod(M3,M4,B3)
    bm3 = (alpha0 - alpha1)*B3
    bm4 = alpha1 * z1
    R= bm1+bm2+bm3+bm4



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
  call syminverse(R,invR,JJ*KK)

      X = 0
  Omega = 0
      !calculate variance

 do i=1,II
          ! per cluster
            do j=1,JJ
              do k = 1,KK
              X((j-1)*KK+k,1) =  X_in(i,j)
              enddo
            enddo
           do j=1,JJ
               do k = 1,KK
               X((j-1)*KK+k,2) = 1
            enddo
          enddo

        New1=matmul(TRANSPOSE(X),invR)
        New2=matmul(New1,X)
Omega=Omega+New2

 enddo

   call syminverse(Omega,invOmega,2)
    vardelta = sigma2*invOmega(1,1)
    ppp=ppnd16(typeone/2,ifault)
     power =alnorm(-beta/sqrt(vardelta)+ppp,upper)+alnorm(beta/sqrt(vardelta)+ppp,upper)

!  power =alnorm(-beta/sqrt(vardelta)-1.959964d0,upper)+alnorm(beta/sqrt(vardelta)-1.959964d0,upper)
    !z_power = pnorm(qnorm(0.05/2) + abs(delta)/sqrt(vardelta))
   !print(z_power)

    ! t-test df=# clusters minus p
  !  df = II - (JJ+1)
  !  t_power = pt(qt(0.05/2,df=df) + abs(delta)/sqrt(vardelta),df=df)


end function GEElinear_continous_notime
