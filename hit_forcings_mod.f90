!!!.....................................................................
!!!.......................HIT FORCINGS MOD..........................
!!!.....................................................................

MODULE hit_forcings_mod

USE parameters_mod
USE variables_mod
USE stats_and_probes_mod
IMPLICIT NONE

REAL(KIND=rk)               :: A_linear
REAL(KIND=rk)               :: KE0
CONTAINS


!!!. . . . . . . . . . . LINEAR FORCING . . . . . . . . . . . . . . . . .
SUBROUTINE linear_forcing_init
IMPLICIT NONE
INTEGER(KIND=ik)                  :: xx,yy,zz
REAL(KIND=rk)                     :: u_rand,v_rand,w_rand

KE0=f_amp

CALL random_seed(put=[seed(1),seed(1)])
CALL B_FFT(uu_C,uu)
ZL10: DO zz=1,nzp
YL10: DO yy=1,nyp
XL10: DO xx=1,nxp
        CALL random_number(u_rand)
        CALL random_number(v_rand)
        CALL random_number(w_rand)
        uu(xx,yy,zz,1)=u_rand*0.3
        uu(xx,yy,zz,2)=v_rand*0.3
        uu(xx,yy,zz,3)=w_rand*0.3
ENDDO XL10
ENDDO YL10
ENDDO ZL10
CALL average_energy(.TRUE.)
CALL F_FFT(uu,uu_C)

END SUBROUTINE linear_forcing_init
!!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!!!. . . . . . . . . . . LINEAR FORCING . . . . . . . . . . . . . . . . .
SUBROUTINE linear_forcing
!! Lundgren linear forcing with Carroll and Blanquart modification
!! See "A proposed modification to Lundgren’s physical space
!! velocity forcing method for isotropic turbulence" 2013

IMPLICIT NONE
INTEGER(KIND=ik)                  :: xx,yy,zz

CALL average_energy(.TRUE.)
KE0=f_amp
A_linear=aa*KE0/(pi*KE)


ZL10: DO zz=1,nzp
YL10: DO yy=1,nyp
XL10: DO xx=1,nxp
        hh(xx,yy,zz,1)=hh(xx,yy,zz,1)+uu(xx,yy,zz,1)*A_linear
        hh(xx,yy,zz,2)=hh(xx,yy,zz,2)+uu(xx,yy,zz,2)*A_linear
        hh(xx,yy,zz,3)=hh(xx,yy,zz,3)+uu(xx,yy,zz,3)*A_linear
ENDDO XL10
ENDDO YL10
ENDDO ZL10
END SUBROUTINE linear_forcing
!!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
END MODULE hit_forcings_mod
