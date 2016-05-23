!!!.....................................................................
!!!.......................HIT FORCINGS MOD..........................
!!!.....................................................................

MODULE hit_forcings_mod

USE parameters_mod
USE variables_mod
USE stats_and_probes_mod

IMPLICIT NONE

!! Linear forcing - - - - - - - - - - - - - - - - - - - - -
REAL(KIND=rk)               :: A_linear
REAL(KIND=rk)               :: KE0
!! Alvelious forcing - - - - - - - - - - - - - - - - - - - -
INTEGER(KIND=ik)                                :: ka,kb
REAL(KIND=rk)                                   :: kf,cc
REAL(KIND=rk)                                   :: Power_in
COMPLEX(KIND=rk)                                :: Aran,Bran
REAL(KIND=rk)                                   :: gA,gB
REAL(KIND=rk),DIMENSION(:,:,:),ALLOCATABLE      :: Fk
COMPLEX(KIND=rk),DIMENSION(:,:,:,:),ALLOCATABLE :: e1,e2
COMPLEX(KIND=rk),DIMENSION(:,:,:,:),ALLOCATABLE :: f_alv


CONTAINS
!!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!! Lundgren linear forcing with Carroll and Blanquart modification
!! See "A proposed modification to Lundgren’s physical space
!! velocity forcing method for isotropic turbulence" 2013
!!!. . . . . . . . . . . LINEAR FORCING . . . . . . . . . . . . . . . . .
SUBROUTINE HIT_linear_forcing_init
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
        uu(xx,yy,zz,1)=0.25+u_rand*0.03
        uu(xx,yy,zz,2)=0.25+v_rand*0.03
        uu(xx,yy,zz,3)=0.25+w_rand*0.03
        ! uu(xx,yy,zz,1)=u_rand*0.3
        ! uu(xx,yy,zz,2)=v_rand*0.3
        ! uu(xx,yy,zz,3)=w_rand*0.3
ENDDO XL10
ENDDO YL10
ENDDO ZL10
CALL STATS_average_energy(.TRUE.)
CALL F_FFT(uu,uu_C)

END SUBROUTINE HIT_linear_forcing_init
!!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!!!. . . . . . . . . . . LINEAR FORCING . . . . . . . . . . . . . . . . .
SUBROUTINE HIT_linear_forcing

IMPLICIT NONE
INTEGER(KIND=ik)                  :: xx,yy,zz

CALL STATS_average_energy(.TRUE.)
KE0=f_amp
! A_linear=aa*KE0/(pi*KE)
A_linear=0.0667
print*,'A =',A_linear
ZL10: DO zz=1,nzp
YL10: DO yy=1,nyp
XL10: DO xx=1,nxp
        hh(xx,yy,zz,1)=hh(xx,yy,zz,1)+uu(xx,yy,zz,1)*A_linear
        hh(xx,yy,zz,2)=hh(xx,yy,zz,2)+uu(xx,yy,zz,2)*A_linear
        hh(xx,yy,zz,3)=hh(xx,yy,zz,3)+uu(xx,yy,zz,3)*A_linear
ENDDO XL10
ENDDO YL10
ENDDO ZL10
END SUBROUTINE HIT_linear_forcing
!!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

!!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!! Alvelious isotropic forcing see "Random forcing of three-dimensional
!! homogenous isotropic turbulence" 1999 Physics of fluids

!!! . . . . . . . . . . . ALVELIUS FORCING INIT . . . . . . . . . . . .
SUBROUTINE HIT_alvelius_forcing_init
!! initialize variables and sets the forcing distribution in the wave numbers
!! TODO The exponential distribution doesn't make much sense as it is
!! TODO never properly approximated on the discret wave numbers (to few points).
!! TODO Change it with something simpler
IMPLICIT NONE

INTEGER(KIND=ik)                     :: xx,yy,zz
REAL(KIND=rk)                        :: int_dist
COMPLEX(KIND=rk)                     :: kw


ka=2
kb=6
kf=3
cc=0.05
Power_in=1.


ALLOCATE(Fk(ka:kb,ka:kb,ka:kb))

ALLOCATE(e1(ka:kb,ka:kb,ka:kb,1:3))
ALLOCATE(e2(ka:kb,ka:kb,ka:kb,1:3))
ALLOCATE(f_alv(ka:kb,ka:kb,ka:kb,1:3))




int_dist=0._rk

ZL10: DO zz=ka,kb
YL10: DO yy=ka,kb
XL10: DO xx=ka,kb
        kw=SQRT(kx(xx)**2+ky(yy)**2+kz(zz)**2)
        int_dist=int_dist+exp(-(ABS(kw)-kf)**2/cc)
ENDDO XL10
ENDDO YL10
ENDDO ZL10


ZL20: DO zz=ka,kb
YL20: DO yy=ka,kb
XL20: DO xx=ka,kb
        kw=( kx(xx)**2+ky(yy)**2+kz(zz)**2 )**0.5
        FK(xx,yy,zz)=Power_in*exp( -(ABS(kw)-kf)**2/cc )/(dt*int_dist)

        e1(xx,yy,zz,1) =  ky(yy)/( kx(xx)**2+ky(yy)**2 )**0.5
        e1(xx,yy,zz,2) = -kx(xx)/( kx(xx)**2+ky(yy)**2 )**0.5
        e1(xx,yy,zz,3) =  CMPLX(0._rk,0._rk)

        e2(xx,yy,zz,1) =  kx(xx)*kz(zz)/( kw*(kx(xx)**2+ky(yy)**2)**0.5 )
        e2(xx,yy,zz,2) =  ky(yy)*kz(zz)/( kw*(kx(xx)**2+ky(yy)**2)**0.5 )
        e2(xx,yy,zz,3) =  -( kx(xx)**2+ky(yy)**2 )**0.5/kw
ENDDO XL20
ENDDO YL20
ENDDO ZL20


END SUBROUTINE HIT_alvelius_forcing_init
!!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!!! . . . . . . . .. . . ALVELIUS FORCING UPDATE. . . . . . . . . . . .
SUBROUTINE HIT_alvelius_forcing_update

IMPLICIT NONE

INTEGER(KIND=ik)                            :: xx,yy,zz,ii
REAL(KIND=rk)                               :: Psi,Phi
REAL(KIND=rk)                               :: Theta1,Theta2
REAL(KIND=rk)                               :: num,den
COMPLEX(KIND=rk)                            :: zeta1,zeta2
COMPLEX(KIND=rk)                            :: kw



CALL random_seed(put=[seed(1),seed(1)])

ZL10: DO zz=ka,kb
YL10: DO yy=ka,kb
XL10: DO xx=ka,kb

        DO ii=1,3
        zeta1=uu_C(xx,yy,zz,ii)*e1(xx,yy,zz,ii)
        zeta2=uu_C(xx,yy,zz,ii)*e2(xx,yy,zz,ii)
        CALL random_number(Psi)
        CALL random_number(Phi)
        Psi=Psi*2_rk*pi
        Phi=Phi*2_rk*pi

        gA=SIN(2._rk*Phi)
        gB=COS(2._rk*Phi)

        num =  gA*REAL(DREAL(zeta1),KIND=rk)  + gB*(SIN(Psi)*REAL(DIMAG(zeta2),KIND=rk)&
               +COS(Psi)*REAL(DIMAG(zeta2),KIND=rk))
        den = -gA*REAL(AIMAG(zeta1),KIND=rk)  + gB*(SIN(Psi)*REAL(DIMAG(zeta2),KIND=rk)&
               -COS(Psi)*REAL(DIMAG(zeta2),KIND=rk))

        Theta1 = ATAN(num/den)
        IF (den == 0._rk) THEN
        CALL random_number(Theta1)
        Theta1=Theta1*2_rk*pi
        ENDIF
        Theta2 = Theta1 + Psi

        kw= kx(xx)**2+ky(yy)**2+kz(zz)**2

        Aran=gA*exp(CMPLX(0,Theta1))*(Fk(xx,yy,zz)/(2_rk*pi*kw))**0.5
        Bran=gB*exp(CMPLX(0,Theta2))*(Fk(xx,yy,zz)/(2_rk*pi*kw))**0.5

        f_alv(xx,yy,zz,ii) = Aran*e1(xx,yy,zz,ii) + Bran*e2(xx,yy,zz,ii)

        ENDDO

ENDDO XL10
ENDDO YL10
ENDDO ZL10
        print*,'Forcing max',maxval(abs(f_alv))
END SUBROUTINE HIT_alvelius_forcing_update
!!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!!! . . . . . . . . . . . . . ALVELIUS FORCING  . . . . . . . . . . . .
SUBROUTINE HIT_alvelius_forcing
IMPLICIT NONE
INTEGER(KIND=ik)                            :: xx,yy,zz,jj

ZL21: DO zz=ka,kb
YL21: DO yy=ka,kb
XL21: DO xx=ka,kb

      hh_C(xx,yy,zz,1)=hh_C(xx,yy,zz,1)+f_alv(xx,yy,zz,1)
      hh_C(xx,yy,zz,2)=hh_C(xx,yy,zz,2)+f_alv(xx,yy,zz,2)
      hh_C(xx,yy,zz,3)=hh_C(xx,yy,zz,3)+f_alv(xx,yy,zz,3)

      hh_C(xx,nyp-kb+yy,zz,1)=hh_C(xx,nyp-kb+yy,zz,1)+f_alv(xx,yy,zz,1)
      hh_C(xx,nyp-kb+yy,zz,2)=hh_C(xx,nyp-kb+yy,zz,2)+f_alv(xx,yy,zz,2)
      hh_C(xx,nyp-kb+yy,zz,3)=hh_C(xx,nyp-kb+yy,zz,3)+f_alv(xx,yy,zz,3)

      hh_C(xx,yy,nzp-kb+zz,1)=hh_C(xx,yy,nzp-kb+zz,1)+f_alv(xx,yy,zz,1)
      hh_C(xx,yy,nzp-kb+zz,2)=hh_C(xx,yy,nzp-kb+zz,2)+f_alv(xx,yy,zz,2)
      hh_C(xx,yy,nzp-kb+zz,3)=hh_C(xx,yy,nzp-kb+zz,3)+f_alv(xx,yy,zz,3)

      hh_C(xx,nyp-kb+yy,nzp-kb+zz,1)=hh_C(xx,nyp-kb+yy,nzp-kb+zz,1)+f_alv(xx,yy,zz,1)
      hh_C(xx,nyp-kb+yy,nzp-kb+zz,2)=hh_C(xx,nyp-kb+yy,nzp-kb+zz,2)+f_alv(xx,yy,zz,2)
      hh_C(xx,nyp-kb+yy,nzp-kb+zz,3)=hh_C(xx,nyp-kb+yy,nzp-kb+zz,3)+f_alv(xx,yy,zz,3)
ENDDO XL21
ENDDO YL21
ENDDO ZL21
END SUBROUTINE HIT_alvelius_forcing
!!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
END MODULE hit_forcings_mod
