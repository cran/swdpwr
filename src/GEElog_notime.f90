function LogPower_GEE_notime(mu, beta, gamma,rho0, II, JJ, KK, &
                     X_in,ICC0,ICC2,typeone) result (power)


    implicit none
    ! ---- arg types -----------------------
    integer :: II, JJ, KK
    integer :: X_in(II,JJ)
    ! II : number of clusters
    ! JJ : number of steps
    ! KK : number of subjects per cluster per step
    double precision :: mu, beta,rho0,alpha0,alpha1,alpha2,ICC0,ICC2,typeone,ppp
    ! true values of mu, beta, and tau2
    double precision :: gamma(JJ),period(2),B1(JJ*KK,JJ*KK)
    !double precision :: a, b        ! lower and upper limits of GL integral
    !integer :: mincomp(JJ+2), maxcomp(JJ+2)
    !integer :: GQ   ! number of GQ points
    !double precision :: GQX(GQ), GQW(GQ)
    !double precision :: t_power
    double precision :: power,lambda1,lambda2,lambda3,lambda4,vardelta
    ! ----------------------------------------

    double precision, parameter :: PI = 3.14159265358979323846   ! pi
    real :: X(JJ,2)  ! intervention matrix
         ! II/(JJ-1). By design, II is a multiple of (JJ-1)
    integer :: i, j, k,m,ifault
    integer :: l
    !logical :: upper

    !double precision, external :: alnorm
    !integer, external :: updatez
    double precision :: M1(JJ,JJ),M2(KK,KK), M3(JJ,JJ), M4(KK,KK)
    double precision :: z0(JJ*KK),z1(JJ*KK,JJ*KK),gmu(JJ),mmu(JJ,JJ)   ! z0: # of subjects with outcome 0; z1: # of subjects with outcome 1.
    double precision :: bm1(JJ*KK,JJ*KK),bm4(JJ*KK,JJ*KK),R(JJ*KK,JJ*KK)
    double precision :: bm2(JJ*KK,JJ*KK),bm3(JJ*KK,JJ*KK)
    double precision :: B2(JJ*KK,JJ*KK),B3(JJ*KK,JJ*KK),Wsub4(2,JJ)
    double precision :: invR(JJ*KK,JJ*KK),Omega(2,2),Ry(2)
    double precision :: invOmega(2,2),Wsub(JJ,2),Wsub3(JJ,JJ),Wsub2(JJ,JJ)
    double precision :: W(JJ,JJ)
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
    double precision, external :: matmual
    upper = .false.
    ifault = 0
    gamma=0
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
    Ry = 0


            period(2)=beta

            period(1)= mu
            mmu=0
             X=1

           Omega=0
           do i=1,II
              do j=1,JJ
                X(j,2)=X_in(i,j)
           end do
           gmu=matmul(X,period)
           do m=1,JJ
             mmu(m,m)=EXP(gmu(m))

             if ((mmu(m,m)*mmu(m,m)+alpha0*mmu(m,m)*(1-mmu(m,m)))>mmu(m,m)) then
                !print *, '(mmu(m,m)*mmu(m,m)+alpha0*mmu(m,m)*(1-mmu(m,m)))>mmu(m,m)'
               !stop
             end if

             if ((mmu(m,m)*mmu(m,m)+alpha0*mmu(m,m)*(1-mmu(m,m)))<max(2*mmu(m,m)-1,0d0)) then
                !print *, '(mmu(m,m)*mmu(m,m)+alpha0*mmu(m,m)*(1-mmu(m,m)))<max(2*mmu(m,m)-1,0d0)'
               !stop
             end if

             mmu(m,m)=(mmu(m,m)*(1-mmu(m,m)))**0.5
             mmu(m,m)=EXP(gmu(m))*1/mmu(m,m)
           enddo

           do m=1,(JJ-1)
             do k=(m+1),JJ
               if (((gmu(m)*gmu(k)+alpha1*((gmu(m)*(1-gmu(m)))**0.5)*((gmu(k)*(1-gmu(k)))**0.5)))>min(gmu(m),gmu(k))) then
                  !print *, '((gmu(m)*gmu(k)+alpha1*((gmu(m)*(1-gmu(m)))**0.5)*((gmu(k)*(1-gmu(k)))**0.5)))>min(gmu(m),gmu(k))'
                 !stop
               end if

               if (((gmu(m)*gmu(k)+alpha1*((gmu(m)*(1-gmu(m)))**0.5)*((gmu(k)*(1-gmu(k)))**0.5)))<max(gmu(m)+gmu(k)-1,0d0)) then
                  !print *, '((gmu(m)*gmu(k)+alpha1*((gmu(m)*(1-gmu(m)))**0.5)*((gmu(k)*(1-gmu(k)))**0.5)))<max(0,gmu(m)+gmu(k)-1)'
                 !stop
               end if
               if (((gmu(m)*gmu(k)+alpha2*((gmu(m)*(1-gmu(m)))**0.5)*((gmu(k)*(1-gmu(k)))**0.5)))>min(gmu(m),gmu(k))) then
                  !print *, '((gmu(m)*gmu(k)+alpha2*((gmu(m)*(1-gmu(m)))**0.5)*((gmu(k)*(1-gmu(k)))**0.5)))>min(gmu(m),gmu(k))'
                 !stop
               end if

               if (((gmu(m)*gmu(k)+alpha2*((gmu(m)*(1-gmu(m)))**0.5)*((gmu(k)*(1-gmu(k)))**0.5)))<max(gmu(m)+gmu(k)-1,0d0)) then
                  !print *, '((gmu(m)*gmu(k)+alpha2*((gmu(m)*(1-gmu(m)))**0.5)*((gmu(k)*(1-gmu(k)))**0.5)))<max(0,gmu(m)+gmu(k)-1)'
                !stop
               end if

             end do
           end do

           Wsub=matmul(mmu,X)
           W=0
           do k=1,JJ
             W(k,k)=KK/lambda3
           end do
           Wsub2=1
           do k=1,JJ
           do m=1,JJ
             Wsub2=(lambda4-lambda3)*KK/(JJ*lambda4*lambda3)
           enddo
           enddo

           Wsub3=W-Wsub2
           Wsub4=matmul(TRANSPOSE(Wsub),Wsub3)
           Omega=Omega+matmul(Wsub4,Wsub)
           end do

           call syminverse(Omega,invOmega,2)
      ! call syminverse(Omega,invOmega,2)


        vardelta = invOmega(2,2)
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
  ppp=ppnd16(typeone/2,ifault)
  power =alnorm(-beta/sqrt(vardelta)+ppp,upper)+alnorm(beta/sqrt(vardelta)+ppp,upper)

  !power =alnorm(-beta/sqrt(vardelta)-1.959964d0,upper)+alnorm(beta/sqrt(vardelta)-1.959964d0,upper)

    !z_power = pnorm(qnorm(0.05/2) + abs(delta)/sqrt(vardelta))
   !print(z_power)

    ! t-test df=# clusters minus p
  !  df = II - (JJ+1)
  !  t_power = pt(qt(0.05/2,df=df) + abs(delta)/sqrt(vardelta),df=df)


end function LogPower_GEE_notime
