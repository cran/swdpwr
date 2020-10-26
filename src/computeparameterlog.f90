subroutine computeparameterlog(JJ, mu, beta, gamma, tau2, p0, p11, rho0)
    implicit none
    ! ---- arg types -----------------------
    integer :: JJ
    double precision :: mu, beta, tau2
    double precision :: gamma(JJ)
    double precision :: p11, rho0
    double precision :: p0(JJ)
    !double precision, parameter:: TwoPi = 8.D0*atan(1.D0)
    ! ---------------------------------------
    ! ---------------------------------------
    integer :: j

    beta=log(p11/p0(1))
    gamma(1)=0
    do j=2,JJ
        gamma(j)=log(p0(j)/p0(1))
    end do

    tau2=log(rho0*(1-p0(1))/p0(1)+1)
    mu=log(p0(1))-tau2/2

end subroutine computeparameterlog
