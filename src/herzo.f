
        SUBROUTINE HERZO(N,X,W)
C
C       ========================================================
C       Purpose : Compute the zeros of Hermite polynomial Ln(x)
C                 in R, and the corresponding
C                 weighting coefficients for Gauss-Hermite
C                 integration by Shanjie Zhang, Jianming Jin.
C       Input :   n    --- Order of the Hermite polynomial
C                 X(n) --- Zeros of the Hermite polynomial
C                 W(n) --- Corresponding weighting coefficients
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION X(N),W(N)
        HN=1.0D0/N
        ZL=-1.1611D0+1.46D0*N**0.5
        DO NR=1,N/2
           IF (NR.EQ.1) Z=ZL
           IF (NR.NE.1) Z=Z-HN*(N/2+1-NR)
           IT=0
10         IT=IT+1
           Z0=Z
           F0=1.0D0
           F1=2.0D0*Z
           DO K=2,N
              HF=2.0D0*Z*F1-2.0D0*(K-1.0D0)*F0
              HD=2.0D0*K*F1
              F0=F1
              F1=HF
           END DO
           P=1.0D0
           DO I=1,NR-1
              P=P*(Z-X(I))
           END DO
           FD=HF/P
           Q=0.0D0
           DO I=1,NR-1
              WP=1.0D0
              DO 25 J=1,NR-1
                 IF (J.EQ.I) GO TO 25
                 WP=WP*(Z-X(J))
25            CONTINUE
              Q=Q+WP
           END DO
           GD=(HD-Q*FD)/P
           Z=Z-FD/GD
           IF (IT.LE.40.AND.DABS((Z-Z0)/Z).GT.1.0D-15) GO TO 10
           X(NR)=Z
           X(N+1-NR)=-Z
           R=1.0D0
           DO K=1,N
              R=2.0D0*R*K
           END DO
           W(NR)=3.544907701811D0*R/(HD*HD)
           W(N+1-NR)=W(NR)
        END DO
        IF (N.NE.2*INT(N/2)) THEN
           R1=1.0D0
           R2=1.0D0
           DO 45 J=1,N
              R1=2.0D0*R1*J
              IF (J.GE.(N+1)/2) R2=R2*J
45         CONTINUE
           W(N/2+1)=0.88622692545276D0*R1/(R2*R2)
           X(N/2+1)=0.0D0
        ENDIF
        RETURN
        END
