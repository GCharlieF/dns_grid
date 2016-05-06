!!!.....................................................................
!!!........................FFTW3 MOD....................................
!!!.....................................................................
MODULE time_advancement_mod
USE parameters_mod
USE variables_mod
USE fft_mod
USE grid_forcing_mod
USE stats_and_probes_mod
USE hit_forcings_mod

implicit none

REAL(KIND=rk),DIMENSION(4,3:4)            :: ark,brk
INTEGER(KIND=ik),DIMENSION(0:3)           :: i_rk
INTEGER(KIND=ik)                          :: n_k

contains

!!!. . . . . . . . . . . . . .CROSS. . . . . . . . . . . . . . . . . . .
FUNCTION cross(a,b)
!compute cross product between complex numbers
 implicit none
 COMPLEX(KIND=rk), DIMENSION(3) :: cross
 COMPLEX(KIND=rk), DIMENSION(3), INTENT(IN) :: a, b
 cross(1) = a(2) * b(3) - a(3) * b(2)
 cross(2) = a(3) * b(1) - a(1) * b(3)
 cross(3) = a(1) * b(2) - a(2) * b(1)
END FUNCTION cross

!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
SUBROUTINE rk_initialize
 !! Sets the Runge-Kutta coefficients for the time integration scheme
 !! brk coefficients are permutated by one position in order to match
 !! the correct time level

      !!Runge-Kutta 3rd order
      ark(1,3)=8._rk/15._rk
      ark(2,3)=5._rk/12._rk
      ark(3,3)=3._rk/4._rk
      brk(1,3)=0._rk
      brk(2,3)=-17._rk/60._rk
      brk(3,3)=-5._rk/12._rk
      !!Runge-Kutta 4th order
      ark(1,4)=8._rk/17._rk
      ark(2,4)=17._rk/60._rk
      ark(3,4)=5._rk/12._rk
      ark(4,4)=3._rk/4._rk
      brk(1,4)=0._rk
      brk(2,4)=-15._rk/68._rk
      brk(3,4)=-17._rk/60._rk
      brk(4,4)=-5._rk/12._rk
END SUBROUTINE rk_initialize

!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
SUBROUTINE partial_right_hand_side
!!TODO: trasformare in una funzione che prende variabili di qualsiasi shape
implicit none
COMPLEX(KIND=rk)                             :: k_quad
REAL(KIND=rk)                                :: pnrk,qnrk
INTEGER(KIND=ik)                             :: zz,yy,xx,jj,kk

 pnrk=2_rk*Re/(dt*ark(n_k,rk_steps)+dt*brk(n_k,rk_steps))

 qnrk=dt*brk(n_k,rk_steps)*pnrk

 ZL10: DO zz=1,nz
 YL10: DO yy=1,ny
 XL10: DO xx=1,nx/2
       jj=day(yy)
       kk=daz(zz)
       k_quad=(kx(xx)**2+ky(yy)**2+kz(zz)**2)
       puu_C(xx,yy,zz,:)=(pnrk+k_quad)*uu_C(xx,jj,kk,:)+qnrk*hh_C(xx,jj,kk,:)

 ENDDO XL10
 ENDDO YL10
 ENDDO ZL10

 print *,'prhs',puu_C(7,12,12,1)
 print *,'prhs',uu_C(7,12,12,1)
 xx=7
 yy=12
 zz=12
 jj=day(yy)
 kk=daz(zz)
 k_quad=(kx(xx)**2+ky(yy)**2+kz(zz)**2)
 print *,pnrk,k_quad,ark(n_k,rk_steps),brk(n_k,rk_steps)
 print *,pnrk+k_quad,qnrk*hh_C(xx,jj,kk,1)
END SUBROUTINE partial_right_hand_side

!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
SUBROUTINE nonlinear
 !! Compute non-linear terms in the phisical space and transform them back
 !! in the Fourier space
 !! Enforce zero-divergence on the trasnformed non-linear terms
 implicit none
 INTEGER(KIND=ik)               :: xx,yy,zz
 INTEGER(KIND=ik)               :: jj,kk
 COMPLEX(KIND=rk)               :: k_quad
 COMPLEX(KIND=rk)               :: ac1,ac2,ac3,div
 REAL(KIND=rk)                  :: norm,ff1,ff2,ff3
 REAL(KIND=rk)                  :: a1,a2,a3,b1,b2,b3
 REAL(KIND=rk)                  :: xf_min,xf_max,delta_xf
 REAL(KIND=rk)                  :: coeff,amp_x


 norm=1_rk/REAL(nxp*nyp*nzp,KIND=rk)

 !! Sets x-bounds of the forced region
 delta_xf=0.5_rk*xl*thick
 xf_min=(-delta_xf + xl/2_rk) * REAL(nxp,KIND=rk)/ xl  !left bound of the forced region
 xf_max=(+delta_xf + xl/2_rk) * REAL(nxp,KIND=rk)/ xl  !right bound of the force region
 ! aa=1.5*pi
 ! coeff=1_rk/(sigma*SQRT(2_rk*pi))


 !! Computes vorticity in the fourier space and puts it in hu hv hw
 ZL10 : DO zz=1,nz
 YL10 :    DO yy=1,ny
 XL10 :       DO xx=1,nx/2
             jj=day(yy)
             kk=daz(zz)
             hh_C(xx,jj,kk,:)=cross((/ kx(xx), ky(yy), kz(zz) /),&
             (/ uu_C(xx,jj,kk,1), uu_C(xx,jj,kk,2), uu_C(xx,jj,kk,3) /))
              ENDDO XL10
          ENDDO YL10
       ENDDO ZL10

    CALL B_FFT(hh_C,hh)
    CALL B_FFT(uu_C,uu)

    print *,'u ta   ::',uu(12,12,12,1)
    print *,'hh ta  ::',hh(12,12,12,1)

    CALL average_energy(stats_time)

!! Compute non-linear term in the phisical space
!! If within the forced region adds the forcing term multiplied by
!! a gaussian distribution function of xx
         ff1=0.;ff2=0.;ff3=0.
  ZL20 : DO zz=1,nzp
  YL20 : DO yy=1,nyp
  XL20 : DO xx=1,nxp
         a1=uu(xx,yy,zz,1)
         a2=uu(xx,yy,zz,2)
         a3=uu(xx,yy,zz,3)

         b1=hh(xx,yy,zz,1)
         b2=hh(xx,yy,zz,2)
         b3=hh(xx,yy,zz,3)

         hh(xx,yy,zz,1)=a2*b3-a3*b2
         hh(xx,yy,zz,2)=a3*b1-a1*b3
         hh(xx,yy,zz,3)=a1*b2-a2*b1

        !  hh(xx,yy,zz,:)=cross((/ uu(xx,yy,zz,1), uu(xx,yy,zz,2), uu(xx,yy,zz,3) /),&
        !  (/ hh(xx,yy,zz,1), hh(xx,yy,zz,2), hh(xx,yy,zz,3) /))

        !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !! Forcing input in the physical space for the grid forcing
         IF (REAL(xx,KIND=rk) > xf_min-nxp/8 .AND. REAL(xx,KIND=rk) < xf_max +nxp/8) THEN

      !    amp_x=coeff*EXP(- (REAL(xx-nxp/2,KIND=rk)*xl/REAL(nxp,KIND=rk) )**2/sigma)
         amp_x=0.5*(1 + TANH(aa*(delta_xf-abs(REAL(nxp/2-xx,KIND=rk))*xl/REAL(nxp,KIND=rk) )))

         hh(xx,yy,zz,1)=hh(xx,yy,zz,1) + fu(yy,zz)*amp_x
         hh(xx,yy,zz,2)=hh(xx,yy,zz,2) + fv(yy,zz)*amp_x
         hh(xx,yy,zz,3)=hh(xx,yy,zz,3) + fw(yy,zz)*amp_x
         ff1=ff1+(fu(yy,zz)*amp_x)**2/REAL(nxp*nyp*nzp,KIND=rk)
         ff2=ff2+(fv(yy,zz)*amp_x)**2/REAL(nxp*nyp*nzp,KIND=rk)
         ff3=ff3+(fw(yy,zz)*amp_x)**2/REAL(nxp*nyp*nzp,KIND=rk)
         ENDIF
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ENDDO XL20
  ENDDO YL20
  ENDDO ZL20
         print *,'ff',ff1,ff2,ff3,hh(nxp/2,nyp/2,nzp/2,1)

   CALL compute_CFL

   ! CALL linear_forcing

   CALL F_FFT(hh,hh_C)
   CALL F_FFT(uu,uu_C)

   ! hh_C(4,4,4,1)=hh_C(4,4,4,1)+CMPLX(1_rk,1_rk)
   ! hh_C(4,4,4,2)=hh_C(4,4,4,2)+CMPLX(1_rk,1_rk)
   ! hh_C(4,4,4,3)=hh_C(4,4,4,3)+CMPLX(1_rk,1_rk)

 ZL30 : DO zz=1,nz
 YL30 :    DO yy=1,ny
            xx=1
            jj=day(yy)
            kk=daz(zz)
              IF(jj==1 .AND. kk==1) THEN
                hh_C(xx,jj,kk,1)=CMPLX(0._rk,0._rk)
                hh_C(xx,jj,kk,2)=CMPLX(0._rk,0._rk)
                hh_C(xx,jj,kk,3)=CMPLX(0._rk,0._rk)
              ELSE
                ac1=hh_C(xx,jj,kk,1)
                ac2=hh_C(xx,jj,kk,2)
                ac3=hh_C(xx,jj,kk,3)

                k_quad=kx(xx)**2+ky(yy)**2+kz(zz)**2
                k_quad=1._rk/k_quad
            !     k_quad=1./max(1.0E-10,abs(k_quad))

                div=kx(xx)*hh_C(xx,jj,kk,1)&
                   +ky(yy)*hh_C(xx,jj,kk,2)&
                   +kz(zz)*hh_C(xx,jj,kk,3)


                hh_C(xx,jj,kk,1)=ac1-div*kx(xx)*k_quad
                hh_C(xx,jj,kk,2)=ac2-div*ky(yy)*k_quad
                hh_C(xx,jj,kk,3)=ac3-div*kz(zz)*k_quad
              ENDIF

 XL31 :       DO xx=2,nx/2
 ! XL31 :       DO xx=1,nx/2

              ac1=hh_C(xx,jj,kk,1)
              ac2=hh_C(xx,jj,kk,2)
              ac3=hh_C(xx,jj,kk,3)

               k_quad=kx(xx)**2+ky(yy)**2+kz(zz)**2
               k_quad=1._rk/k_quad
            !    k_quad=1./max(1.0E-10,abs(k_quad))

               div=kx(xx)*hh_C(xx,jj,kk,1)&
                  +ky(yy)*hh_C(xx,jj,kk,2)&
                  +kz(zz)*hh_C(xx,jj,kk,3)

                  hh_C(xx,jj,kk,1)=ac1-div*kx(xx)*k_quad
                  hh_C(xx,jj,kk,2)=ac2-div*ky(yy)*k_quad
                  hh_C(xx,jj,kk,3)=ac3-div*kz(zz)*k_quad


 ENDDO XL31
 ENDDO YL30
 ENDDO ZL30
      ! CALL divfree(hh_C)


END SUBROUTINE nonlinear

!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!! . . . . . . . . . . . . . . . . LINEAR . . . . . . . . . . . . . . . . . . .
SUBROUTINE linear
 implicit none
 INTEGER(KIND=ik)               :: zz,yy,xx,jj,kk
 COMPLEX(KIND=rk)               :: k_quad
 REAL(KIND=rk)                  :: den,pnrk,qnrk

 ! n_k=i_rk(mod(it,rk_steps))

 pnrk=2_rk*Re/(dt*ark(n_k,rk_steps)+dt*brk(n_k,rk_steps))

 qnrk=dt*ark(n_k,rk_steps)*pnrk

 ZL10 : DO zz=1,nz
 YL10 : DO yy=1,ny
 XL10 : DO xx=1,nx/2
               jj=day(yy)
               kk=daz(zz)
               k_quad=kx(xx)**2+ky(yy)**2+kz(zz)**2
               den=1./(pnrk-k_quad)
               uu_C(xx,jj,kk,:)=(puu_C(xx,yy,zz,:)+qnrk*hh_C(xx,jj,kk,:))*den
 ENDDO XL10
 ENDDO YL10
 ENDDO ZL10
END SUBROUTINE linear
!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
SUBROUTINE divfree(vv_C)

IMPLICIT NONE

INTEGER(KIND=IK)                                   ::xx,yy,zz,jj,kk
COMPLEX(KIND=rk)                                   :: k_quad
COMPLEX(KIND=rk)                                   :: ac1,ac2,ac3,div
COMPLEX(KIND=rk),DIMENSION(:,:,:,:), INTENT(INOUT) :: vv_C

ZL30 : DO zz=1,nz
YL30 :    DO yy=1,ny
           xx=1
           jj=day(yy)
           kk=daz(zz)
            !  IF(jj==1 .AND. kk==1) THEN
            !    vv_C(xx,jj,kk,1)=CMPLX(0._rk,0._rk)
            !    vv_C(xx,jj,kk,2)=CMPLX(0._rk,0._rk)
            !    vv_C(xx,jj,kk,3)=CMPLX(0._rk,0._rk)
            !  ELSE
               ac1=vv_C(xx,jj,kk,1)
               ac2=vv_C(xx,jj,kk,2)
               ac3=vv_C(xx,jj,kk,3)

               k_quad=kx(xx)**2+ky(yy)**2+kz(zz)**2
               k_quad=1._rk/k_quad
           !     k_quad=1./max(1.0E-10,abs(k_quad))

               div=kx(xx)*vv_C(xx,jj,kk,1)&
                  +ky(yy)*vv_C(xx,jj,kk,2)&
                  +kz(zz)*vv_C(xx,jj,kk,3)
             IF(jj==1 .AND. kk==1) THEN
                   k_quad=CMPLX(0._rk,0._rk)
             ENDIF

               vv_C(xx,jj,kk,1)=ac1-div*kx(xx)*k_quad
               vv_C(xx,jj,kk,2)=ac2-div*ky(yy)*k_quad
               vv_C(xx,jj,kk,3)=ac3-div*kz(zz)*k_quad
            !  ENDIF

XL31 :       DO xx=2,nx/2
! XL31 :       DO xx=1,nx/2

             ac1=vv_C(xx,jj,kk,1)
             ac2=vv_C(xx,jj,kk,2)
             ac3=vv_C(xx,jj,kk,3)

              k_quad=kx(xx)**2+ky(yy)**2+kz(zz)**2
              k_quad=1._rk/k_quad
           !    k_quad=1./max(1.0E-10,abs(k_quad))

              div=kx(xx)*vv_C(xx,jj,kk,1)&
                 +ky(yy)*vv_C(xx,jj,kk,2)&
                 +kz(zz)*vv_C(xx,jj,kk,3)

                 vv_C(xx,jj,kk,1)=ac1-div*kx(xx)*k_quad
                 vv_C(xx,jj,kk,2)=ac2-div*ky(yy)*k_quad
                 vv_C(xx,jj,kk,3)=ac3-div*kz(zz)*k_quad


ENDDO XL31
ENDDO YL30
ENDDO ZL30
END SUBROUTINE divfree


!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

END MODULE time_advancement_mod

!!!.....................................................................
!!!.....................................................................
!!!.....................................................................
!!TODO routines per statistiche: energia,dissipazione,potenza
