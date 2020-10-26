!
!

function LogitPower_time(mu, beta, gamma, tau2, II, JJ, KK, GQ, GQX, GQW, X_in,typeone) result (power)
    implicit none
  ! ---- arg types -----------------------
    integer :: II, JJ, KK,ifault
    integer :: X_in(II,JJ)    ! II : number of clusters
    ! JJ : number of steps
    ! KK : number of subjects per cluster per step
    double precision :: mu, beta, tau2,typeone,ppp
    double precision :: gamma(JJ)
    ! true values of beta0, beta1, and tau2
    integer :: GQ   ! number of GQ points
    double precision :: GQX(GQ), GQW(GQ)
    double precision :: power
    ! ----------------------------------------
    double precision, parameter :: PI = 3.14159265358979323846   ! pi
    integer :: X(II,JJ)   ! intervention matrix
    !integer :: DD      ! II/(JJ-1). By design, II is a multiple of (JJ-1)
    integer :: i, j, k
    logical :: upper
    double precision, external :: alnorm
    double precision, external :: ppnd16
    integer :: z0(JJ), z1(JJ)   ! z0: # of subjects with outcome 0; z1: # of subjects with outcome 1.
    integer :: nZ               ! total number of possible combination of z0 and z1
    double precision :: derlikelihood(JJ+2)   ! derivative of likelihood
    double precision :: derlikelihood2(JJ+2, JJ+2)  ! square of derivative of likelihood
    double precision :: invVar(JJ+2,JJ+2)      ! inverse of variance matrix
    double precision :: Var(JJ+2,JJ+2),invVarPre(JJ+2,JJ+2)         ! variance matrix
    double precision :: prob
    double precision :: sebeta  ! se of beta
    ifault = 0
    upper = .false.

  !  DD = II/(JJ-1)   ! II is a multiple of (JJ-1)
    ! assign intervention
    X = 0
    do i=1,II
      do j = 1,JJ
        X(i,j) =  X_in(i,j)
      enddo
    enddo

    ! compute the possible combinations of (z0,z1)
    nZ = 1
    do j=1,JJ
        nZ = nZ * (KK+1)
    end do
    ! run the algorithm
    invVar = 0.0d0
    do i=1,II
      if(i>1 .and. ALL(X(i,:).EQ. X(i-1,:))) then
      invVar = invVar + invVarPre
      !finish = updatez(z0, JJ, KK)
      else
        ! per cluster
        invVarPre=0
        z0 = 0
        do k=1, nZ
            z1 = KK - z0
            call der_likelihood_timelogit(mu,beta,gamma,tau2, z0, z1, X(i,:), JJ, KK, GQ, GQX, GQW, derlikelihood, prob)
            call vectorsquare(derlikelihood, JJ+2, derlikelihood2)
            invVar = invVar + derlikelihood2 * prob
            invVarPre = invVarPre+derlikelihood2 * prob
            call updatezz(z0, JJ, KK)
        end do
      end if
    end do
    call syminverse(invVar,Var,JJ+2)
    sebeta = sqrt(Var(2,2))
    ppp=ppnd16(typeone/2,ifault)
    power = alnorm(beta/sebeta+ppp,upper) + alnorm(-beta/sebeta+ppp,upper)
  !if(beta<0) then
  !power =alnorm(-beta/sebeta-1.959964d0,upper)
!else
!  power =alnorm(beta/sebeta-1.959964d0,upper)
!end if
end function LogitPower_time

subroutine updatezz(z0, JJ, KK)
    implicit none
    ! ---- arg types -----------------------
    integer :: JJ, KK
    integer :: z0(JJ)
    ! --------------------------------------
    integer :: j
    z0(1) = z0(1) + 1
    do j=1,JJ-1
        if (z0(j)>KK) then
            z0(j) = 0
            z0(j+1) = z0(j+1) + 1
        else
            exit
        end if
    end do
end subroutine updatezz


subroutine der_likelihood_timelogit(mu,beta,gamma,tau2, z0, z1, XX, JJ, KK, GQ, GQX, GQW, derlikelihood, prob)
    implicit none
    ! ---- arg types -----------------------
    integer :: JJ, KK
    double precision :: mu, beta, tau2
    double precision :: gamma(JJ)
    ! true values of mu, beta, gamma, tau2
    integer :: z0(JJ), z1(JJ)   ! z0: # of subjects with outcome 0; z1: # of subjects with outcome 1.
    integer :: XX(JJ)     ! treatment assignment
    integer :: GQ   ! number of GQ points
    double precision :: GQX(GQ), GQW(GQ)
    double precision :: derlikelihood(JJ+2)
    double precision :: prob
    ! --------------------------------------
    double precision, parameter :: PI = 3.14159265358979323846   ! pi
    double precision :: likelihoodf
    double precision :: x
    integer :: i, j, k
    double precision :: ff0, ff1, ff01
    double precision :: ff, ffprob
    double precision :: ff_mu, ff_beta
    double precision :: ff_gamma(JJ-1)
    double precision :: temp

    derlikelihood = 0.0d0
    likelihoodf = 0.0d0
    prob = 0.0d0
    do i=1,GQ
        x = GQX(i)
        ff = 1.0d0
        ffprob = 1.0d0
        ff_mu = 0.0d0
        ff_beta = 0.0d0
        do j=1,JJ
            ff0 = 1.0d0/(1.0d0+exp(mu+beta*XX(j)+gamma(j)+sqrt(2.0d0*tau2)*x))
            ff1 = 1-ff0
            ff = ff*(ff0**z0(j))*(ff1**z1(j))
            temp = z1(j)*ff0 - z0(j)*ff1
            ff_mu = ff_mu + temp
            ff_beta = ff_beta + XX(j)*temp
            k = j-1
            if (k>0) then
                ff_gamma(k) = temp
            end if

            ff01 = ff0*ff1
            ! compute binomial
            ! compute combination number with power of ff0 and ff1
            ! to avoid overflow
            if (z0(j)<z1(j)) then
                ffprob = ffprob*ff1**(z1(j)-z0(j))
                do k=0,(z0(j)-1)
                    ffprob = ffprob*dble(KK-k)/dble(z0(j)-k)*ff01
                end do
            else
                ffprob = ffprob*ff0**(z0(j)-z1(j))
                do k=0,(z1(j)-1)
                    ffprob = ffprob*dble(KK-k)/dble(z1(j)-k)*ff01
                end do
            end if
        end do
        prob = prob + GQW(i)*ffprob
        likelihoodf = likelihoodf + GQW(i)*ff
        derlikelihood(1) = derlikelihood(1) + GQW(i)*ff*ff_mu
        derlikelihood(2) = derlikelihood(2) + GQW(i)*ff*ff_beta
        derlikelihood(3:(JJ+1)) = derlikelihood(3:(JJ+1)) + GQW(i)*ff*ff_gamma
        derlikelihood(JJ+2) = derlikelihood(JJ+2) + GQW(i)*ff*(x*x-0.5d0)/tau2
    enddo
    derlikelihood = derlikelihood / likelihoodf
    prob = prob/sqrt(pi)

end subroutine der_likelihood_timelogit

function LogitPower_notime(mu, beta, tau2, II, JJ, KK, GQ, GQX, GQW,X_in,typeone) result (power)
    implicit none
    ! ---- arg types -----------------------
    double precision :: mu, beta, tau2,typeone,ppp
    integer :: II, JJ, KK,ifault
    integer :: X_in(II,JJ)
    ! true values of beta0, beta1, and tau2

    ! II : number of clusters
    ! JJ : number of steps
    ! KK : number of subjects per cluster per step
    integer :: GQ   ! number of GQ points
    double precision :: GQX(GQ), GQW(GQ)
    double precision :: power
    ! ----------------------------------------
    integer :: NI   ! number of subjects per cluster, NI=JJ*KK
    !integer :: DD   ! II/(JJ-1). By design, II is a multiple of (JJ-1)

    double precision, parameter :: PI = 3.14159265358979323846   ! pi
    integer :: i, j
    integer :: z0(II), z1(II)
    integer :: z00, z01, z10, z11
    double precision :: h11, h12, h13, h22, h23, h33,h11pre,h12pre,h13pre,h22pre,h23pre,h33pre
    double precision :: derlikelihood_mu, derlikelihood_beta, derlikelihood_tau2
    double precision, external :: prob_z
    double precision :: prob
    double precision :: sebeta  ! se of beta1
    logical :: upper
    double precision, external :: ppnd16
    double precision, external :: alnorm

    upper = .false.
    ifault = 0

    NI = JJ * KK
    !DD = II/(JJ-1)


    z0 = 0
    do j=1,II
        z0(j) =(JJ-sum(X_in(j,:)))*KK
    enddo
    z1 = NI - z0

    h11 = 0.0d0
    h12 = 0.0d0
    h13 = 0.0d0
    h22 = 0.0d0
    h23 = 0.0d0
    h33 = 0.0d0
    do i=1,II
      if(i>1 .and. ALL(X_in(i,:).EQ. X_in(i-1,:))) then
        h11=h11+h11pre
        h22=h22+h22pre
        h33=h33+h33pre
        h12=h12+h12pre
        h13=h13+h13pre
        h23=h23+h23pre
      else
        h11pre=0
        h22pre=0
        h33pre=0
        h12pre=0
        h13pre=0
        h23pre=0
        do z00 = 0,z0(i)
            z01 = z0(i) - z00
            do z10 = 0,z1(i)
                z11 = z1(i) - z10
                call der_likelihood_notimelogit(mu, beta, tau2, z00, z01, z10, z11, GQ, GQX, GQW, derlikelihood_mu, &
                derlikelihood_beta, derlikelihood_tau2, prob)
                h11 = h11 + derlikelihood_mu*derlikelihood_mu*prob
                h22 = h22 + derlikelihood_beta*derlikelihood_beta*prob
                h33 = h33 + derlikelihood_tau2*derlikelihood_tau2*prob
                h12 = h12 + derlikelihood_mu*derlikelihood_beta*prob
                h13 = h13 + derlikelihood_mu*derlikelihood_tau2*prob
                h23 = h23 + derlikelihood_beta*derlikelihood_tau2*prob
                h11pre=h11pre+derlikelihood_mu*derlikelihood_mu*prob
                h22pre=h22pre+derlikelihood_beta*derlikelihood_beta*prob
                h33pre = h33pre + derlikelihood_tau2*derlikelihood_tau2*prob
                h12pre = h12pre + derlikelihood_mu*derlikelihood_beta*prob
                h13pre = h13pre + derlikelihood_mu*derlikelihood_tau2*prob
                h23pre = h23pre + derlikelihood_beta*derlikelihood_tau2*prob
            enddo
        enddo
      end if
    enddo
    sebeta = sqrt(abs((h33*h11-h13*h13)/(h11*h22*h33+2.0d0*h12*h23*h13-h13*h13*h22-h12*h12*h33-h23*h23*h11)))

  !  if(beta<0) then
  !  power =alnorm(-beta/sebeta-1.959964d0,upper)
  !else
  !  power =alnorm(beta/sebeta-1.959964d0,upper)
  !end if
  ppp=ppnd16(typeone/2,ifault)
  power = alnorm(beta/sebeta+ppp,upper) + alnorm(-beta/sebeta+ppp,upper)
end function LogitPower_notime


subroutine der_likelihood_notimelogit(mu, beta, tau2, z00, z01, z10, z11, GQ, GQX, GQW, derlikelihood_mu, &
derlikelihood_beta, derlikelihood_tau2, prob)
    implicit none
    ! ---- arg types -----------------------
    double precision :: mu, beta, tau2
    ! true values of beta0, beta1, and tau2
    integer :: z00, z01, z10, z11
    integer :: GQ   ! number of GQ points
    double precision :: GQX(GQ), GQW(GQ)
    double precision :: derlikelihood_mu
    double precision :: derlikelihood_beta
    double precision :: derlikelihood_tau2
    double precision :: prob
    ! ---------------------------------------
    double precision :: likelihoodf

    double precision :: ff00, ff01, ff10, ff11
    double precision :: ff, ff1, ffprob, ff0prob, ff1prob
    integer :: z0, z1
    integer :: i, k
    double precision :: x

    double precision, parameter :: PI = 3.14159265358979323846   ! pi

    z0 = z00 + z01
    z1 = z10 + z11
    derlikelihood_mu = 0.0d0
    derlikelihood_beta = 0.0d0
    derlikelihood_tau2 = 0.0d0
    likelihoodf = 0.0d0
    prob = 0.0d0
    do i=1,GQ
        x = GQX(i)
        ff = 1.0d0
        ffprob = 1.0d0
        ff00 = 1.0d0/(1.0d0+exp(mu+sqrt(2.0d0*tau2)*x))
        ff01 = 1-ff00
        ff10 = 1.0d0/(1.0d0+exp(mu+beta+sqrt(2.0d0*tau2)*x))
        ff11 = 1-ff10

        ff = (ff00**z00)*(ff01**z01)*(ff10**z10)*(ff11**z11)
        likelihoodf = likelihoodf + GQW(i) * ff
        ff1 = ff*(z01*ff00-z00*ff01+z11*ff10-z10*ff11)
        derlikelihood_mu = derlikelihood_mu + GQW(i) * ff1
        ff1 = ff*(z11*ff10-z10*ff11)
        derlikelihood_beta = derlikelihood_beta + GQW(i) * ff1
        ff1 = ff*(x*x-0.5d0)/tau2
        derlikelihood_tau2 = derlikelihood_tau2 + GQW(i) * ff1

        ! prob
        ff0prob = ff00*ff01
        ff1prob = ff10*ff11
        ! compute combination number with power of ff00 and ff01
        ! to avoid overflow
        if (z00<z01) then
            ffprob = ffprob*ff01**(z01-z00)
            do k=0,(z00-1)
                ffprob = ffprob*dble(z0-k)/dble(z00-k)*ff0prob
            end do
        else
            ffprob = ffprob*ff00**(z00-z01)
            do k=0,(z01-1)
                ffprob = ffprob*dble(z0-k)/dble(z01-k)*ff0prob
            end do
        end if
        ! compute combination number with power of ff10 and ff11
        ! to avoid overflow
        if (z10<z11) then
            ffprob = ffprob*ff11**(z11-z10)
            do k=0,(z10-1)
                ffprob = ffprob*dble(z1-k)/dble(z10-k)*ff1prob
            end do
        else
            ffprob = ffprob*ff10**(z10-z11)
            do k=0,(z11-1)
                ffprob = ffprob*dble(z1-k)/dble(z11-k)*ff1prob
            end do
        end if
        prob = prob + GQW(i) * ffprob
    enddo
    derlikelihood_mu = derlikelihood_mu / likelihoodf
    derlikelihood_beta = derlikelihood_beta / likelihoodf
    derlikelihood_tau2 = derlikelihood_tau2 / likelihoodf
    ! likelihoodf = likelihoodf/sqrt(pi)
    prob = prob / sqrt(pi)
end subroutine der_likelihood_notimelogit
