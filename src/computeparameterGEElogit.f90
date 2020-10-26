subroutine computeparameterGEElogit(JJ, mu, beta, gamma, tau2, p0, p11, rho0)
    implicit none
    ! ---- arg types -----------------------
    integer :: JJ
    double precision :: mu, beta, tau2
    double precision :: gamma(JJ)
    double precision :: p11, rho0
    double precision :: p0(JJ)
    ! ---------------------------------------
    ! ---------------------------------------
    integer :: j

    mu = log(1/(1-p0(1))-1)
    beta = log(1/(1-p11)-1)-mu
    tau2 = rho0*mu*(1-mu)
    ! gamma
    gamma(1) = 0.0d0
    do j=2,JJ
        gamma(j) =  log(1/(1-p0(j))-1) - mu
    end do
end subroutine computeparameterGEElogit
