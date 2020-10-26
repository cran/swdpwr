subroutine computeparameterlogit(JJ, mu, beta, gamma, tau2, p0, p11, rho0, GQ, GQX, GQW, convergence)
    implicit none
    ! ---- arg types -----------------------
    integer :: JJ
    double precision :: mu, beta, tau2
    double precision :: gamma(JJ)
    double precision :: p11, rho0
    double precision :: p0(JJ)
    integer :: GQ   ! number of GQ points
    double precision :: GQX(GQ), GQW(GQ)
    integer :: convergence
    ! ---------------------------------------
    ! ---------------------------------------
    double precision :: f2(2)
    double precision :: derf2(2,2)
    double precision :: invderf2(2,2)
    double precision :: f, derf, invderf
    double precision :: gammatemp
    !integer :: iter
    integer :: i, j

    double precision, parameter :: epsilon=1.0d-5
    integer, parameter :: maxiter=100

    convergence = 0
    ! compute mu and tau2
    do i=1,maxiter
        call computef_mutau2(f2, mu, tau2, p0(1), rho0, GQ, GQX, GQW)
        call derivativef_mutau2(derf2, mu, tau2, p0(1), rho0, GQ, GQX, GQW)
        if (maxval(abs(f2))<epsilon) then
            convergence = convergence + 1
            exit
        end if
        call inverse(derf2,invderf2,2)
        mu = mu - dot_product(invderf2(1,:), f2)
        tau2 = tau2 - dot_product(invderf2(2,:), f2)
        if (tau2<epsilon) tau2 = epsilon
    end do
    ! compute beta
    do i=1,maxiter
        call computef(f, mu, beta, tau2, p11, GQ, GQX, GQW)
        call derivativef(derf, mu, beta, tau2, p11, GQ, GQX, GQW)
        if (abs(f)<epsilon) then
            convergence = convergence + 1
            exit
        end if
        invderf = 1.0d0 / derf
        beta = beta - invderf * f
    end do
    ! gamma
    gamma(1) = 0.0d0
    do j=2,JJ
        gammatemp = 0.0d0
        do i=1,maxiter
            call computef(f, mu, gammatemp, tau2, p0(j), GQ, GQX, GQW)
            call derivativef(derf, mu, gammatemp, tau2, p0(j), GQ, GQX, GQW)
            if (abs(f)<epsilon) then
                gamma(j) = gammatemp
                convergence = convergence + 1
                exit
            end if
            invderf = 1/derf
            gammatemp = gammatemp - invderf * f
        end do
    end do

    if (convergence==(JJ+1)) then
        convergence=1
    else
        convergence=0
    end if
end subroutine computeparameterlogit

subroutine computef_mutau2(f, mu, tau2, p01, rho0, GQ, GQX, GQW)
    implicit none
    ! ---- arg types -----------------------
    double precision :: f(2)
    double precision :: mu
    double precision :: tau2
    double precision :: p01, rho0
    integer :: GQ   ! number of GQ points
    double precision :: GQX(GQ), GQW(GQ)
    ! ---------------------------------------
    double precision, parameter :: PI = 3.14159265358979323846   ! pi
    double precision :: intp0
    ! intp0 = int p0
    double precision :: intp02
    ! intp02 = int (p0)^2
    double precision :: intp0p00
    ! intp0p00 = int p0(1-p0)
    integer :: i
    double precision :: x
    double precision :: ff00, ff01
    double precision :: denominator, numerator

    intp0 = 0.0d0
    intp02 = 0.0d0
    intp0p00 = 0.0d0
    do i=1,GQ
        x = GQX(i)
        ff00 = 1.0d0/(1.0d0+exp(mu+sqrt(2.0d0*tau2)*x))
        ff01 = 1-ff00
        intp0 = intp0 + GQW(i) * ff01
        intp02 = intp02 + GQW(i) * ff01 * ff01
        intp0p00 = intp0p00 + GQW(i) * ff00 * ff01
    enddo
    intp0 = intp0 / sqrt(pi)
    intp02 = intp02 / sqrt(pi)
    intp0p00 = intp0p00 / sqrt(pi)

    f(1) = intp0 - p01
    denominator = intp0p00+intp02-intp0*intp0
    numerator = intp02-intp0*intp0
    f(2) = numerator - rho0 * denominator
end subroutine computef_mutau2

subroutine derivativef_mutau2(derf, mu, tau2, p01, rho0, GQ, GQX, GQW)
    implicit none
    ! ---- arg types -----------------------
    double precision :: derf(2,2)
    ! double precision :: f(2)
    double precision :: mu
    double precision :: tau2
    double precision :: p01, rho0
    integer :: GQ   ! number of GQ points
    double precision :: GQX(GQ), GQW(GQ)
    ! ---------------------------------------
    double precision, parameter :: PI = 3.14159265358979323846   ! pi
    double precision :: intp0
    ! intp0 = int p0
    double precision :: intp0tau
    double precision :: intp02tau
    double precision :: intp0p00
    ! intp0p00 = int p0(1-p0
    double precision :: intp0p00tau
    double precision :: intp02p00
    double precision :: intp00temp
    integer :: i
    double precision :: x
    double precision :: ff00, ff01

    intp0 = 0.0d0
    intp0tau = 0.0d0
    intp02tau = 0.0d0
    intp0p00 = 0.0d0
    intp0p00tau = 0.0d0
    intp02p00 = 0.0d0
    intp00temp = 0.0d0
    p01= p01
    do i=1,GQ
        x = GQX(i)
        ff00 = 1.0d0/(1.0d0+exp(mu+sqrt(2.0d0*tau2)*x))
        ff01 = 1-ff00
        intp0 = intp0 + GQW(i) * ff01
        intp0tau = intp0tau + GQW(i) * ff01 * (x*x-0.5d0)/tau2
        intp02tau = intp02tau + GQW(i) * ff01 * ff01 * (x*x-0.5d0)/tau2
        intp0p00 = intp0p00 + GQW(i) * ff00 * ff01
        intp0p00tau = intp0p00tau + GQW(i)*ff00*ff01*(x*x-0.5d0)/tau2
        intp02p00 = intp02p00 + GQW(i) * ff01*ff01*ff00
        intp00temp = intp00temp + GQW(i) * ff01*ff00*(ff00-ff01)
    enddo
    intp0 = intp0 / sqrt(pi)
    intp0tau = intp0tau / sqrt(pi)
    intp02tau = intp02tau / sqrt(pi)
    intp0p00 = intp0p00 / sqrt(pi)
    intp0p00tau = intp0p00tau / sqrt(pi)
    intp02p00 = intp02p00 / sqrt(pi)
    intp00temp = intp00temp / sqrt(pi)

    derf(1,1) = intp0p00
    derf(2,1) = 2*(1-rho0)*intp02p00 -2*(1-rho0)*intp0*intp0p00-rho0*intp00temp
    derf(1,2) = intp0tau
    derf(2,2) = (1-rho0)*intp02tau - 2*(1-rho0)*intp0*intp0tau - rho0*intp0p00tau
end subroutine derivativef_mutau2

subroutine computef(f, mu, beta, tau2, p1, GQ, GQX, GQW)
    implicit none
    ! ---- arg types -----------------------
    double precision :: f
    double precision :: mu
    double precision :: beta
    double precision :: tau2
    double precision :: p1
    integer :: GQ   ! number of GQ points
    double precision :: GQX(GQ), GQW(GQ)
    ! ---------------------------------------
    double precision, parameter :: PI = 3.14159265358979323846   ! pi
    double precision :: intp
    integer :: i
    double precision :: x
    double precision :: ff

    intp = 0.0d0
    do i=1,GQ
        x = GQX(i)
        ff = 1.0d0 - 1.0d0/(1.0d0+exp(mu+beta+sqrt(2.0d0*tau2)*x))
        intp = intp + GQW(i) * ff
    enddo
    intp = intp / sqrt(pi)

    f = intp - p1

end subroutine computef

subroutine derivativef(derf, mu, beta, tau2, p1, GQ, GQX, GQW)
    implicit none
    ! ---- arg types -----------------------
    double precision :: derf
    double precision :: mu
    double precision :: beta
    double precision :: tau2
    double precision :: p1
    integer :: GQ   ! number of GQ points
    double precision :: GQX(GQ), GQW(GQ)
    ! ---------------------------------------
    double precision, parameter :: PI = 3.14159265358979323846   ! pi
    double precision :: intp0p00
    integer :: i
    double precision :: x
    double precision :: ff00, ff01
    p1=0+p1
    intp0p00 = 0.0d0
    do i=1,GQ
        x = GQX(i)
        ff00 = 1.0d0/(1.0d0+exp(mu+beta+sqrt(2.0d0*tau2)*x))
        ff01 = 1-ff00
        intp0p00 = intp0p00 + GQW(i) * ff00 * ff01
    enddo
    derf = intp0p00 / sqrt(pi)

end subroutine derivativef
