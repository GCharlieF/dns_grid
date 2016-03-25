!!!.....................................................................
!!!........................FFTW3 MOD....................................
!!!.....................................................................
MODULE time_advancement_mod
USE parameters_mod
USE variables_and_IO_mod !, ONLY: i_kutta,Re,dt
USE fft_mod
implicit none
REAL(KIND=rk),DIMENSION(4,3:4)            :: ark,brk
INTEGER(KIND=ik),DIMENSION(0:3)           :: i_rk
INTEGER(KIND=ik)                          :: n_k
contains

!!!. . . . . . . . . . . . . .CROSS. . . . . . . . . . . . . . . . . . .
FUNCTION cross(a,b)
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
      !Runge-Kutta 3rd order
      ark(1,3)=8._rk/15._rk
      ark(2,3)=5._rk/12._rk
      ark(3,3)=3._rk/4._rk
      brk(1,3)=0._rk
      brk(2,3)=-17._rk/60._rk
      brk(3,3)=-5._rk/12._rk
      !Runge-Kutta 4th order
      ark(1,4)=8._rk/17._rk
      ark(2,4)=17._rk/60._rk
      ark(3,4)=5._rk/12._rk
      ark(4,4)=3._rk/4._rk
      brk(1,4)=0._rk
      brk(2,4)=-15._rk/68._rk
      brk(3,4)=-17._rk/60._rk
      brk(4,4)=-5._rk/12._rk
      !Step counter. Update through i_rk(mod(it,i_kutta))
      i_rk(0)=rk_steps
      i_rk(1)=1_ik
      i_rk(2)=2_ik
      i_rk(3)=3_ik
END SUBROUTINE rk_initialize

!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
SUBROUTINE partial_right_hand_side
!!TODO: trasformare in una funzione che prende variabili di qualsiasi shape
implicit none
COMPLEX(KIND=rk)                             :: k_quad
REAL(KIND=rk)                                :: pnrk,qnrk
INTEGER(KIND=ik)                             :: zz,yy,xx,xr,xi

 ! n_k=i_rk(mod(it,i_kutta))

 pnrk=CMPLX(2_rk*Re/(dt*ark(n_k,rk_steps)+dt*brk(n_k,rk_steps)),0._rk)

 qnrk=CMPLX(dt*brk(n_k,rk_steps),0._rk)*pnrk

 ZL10: DO zz=1,nz
 YL10: DO yy=1,ny
 XL10: DO xx=1,nx/2

       k_quad=kx(xx)**2+ky(yy)**2+kz(zz)**2
       puu_C(xx,yy,zz,:)=(pnrk+k_quad)*uu_C(xx,yy,zz,:)+qnrk*hh_C(xx,yy,zz,:)


 ENDDO XL10
 ENDDO YL10
 ENDDO ZL10
END SUBROUTINE partial_right_hand_side

!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
SUBROUTINE nonlinear
 implicit none
 INTEGER(KIND=ik)               :: xx,yy,zz,xr,xi
 COMPLEX(KIND=rk)               :: k_quad,ac1,ac2,ac3,div
 REAL(KIND=rk)                  :: norm
 REAL(KIND=rk)                  :: a1,a2,a3,b1,b2,b3
 REAL(KIND=rk)                  :: xf_min,xf_max,delta_xf   !forcing geom
 REAL(KIND=rk)                  :: coeff,amp_x
 REAL(KIND=rk)                  :: CFL,dx,dy,dz,u_max,v_max,w_max


 norm=1_rk/REAL(nxp*nyp*nzp,KIND=rk)

 ! delta_xf=0.5_rk*xl*thick
 ! xf_min=(-delta_xf + xl/2_rk) * REAL(nxp,KIND=rk)/ xl  !left bound of the forced region
 ! xf_max=(+delta_xf + xl/2_rk) * REAL(nxp,KIND=rk)/ xl  !right bound of the force region
 !
 ! coeff=1_rk/(sigma*SQRT(2_rk*pi))


 !Compute vorticity in the fourier space and put it in hu hv hw
 ZL10 : DO zz=1,nz
 YL10 :    DO yy=1,ny
 XL10 :       DO xx=1,nx/2
             hh_C(xx,yy,zz,:)=cross((/ kx(xx), ky(yy), kz(zz) /),&
             (/ uu_C(xx,yy,zz,1), uu_C(xx,yy,zz,2), uu_C(xx,yy,zz,3) /))
              ENDDO XL10
          ENDDO YL10
       ENDDO ZL10

    CALL B_FFT(hh_C,hh)
    CALL B_FFT(uu_C,uu)

! Compute non-linear term in the phisical space
! If within the forced region adds the forcing term multiplied by
! a gaussian distribution function of xx

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


      !    IF (REAL(xx,KIND=rk) > xf_min .AND. REAL(xx,KIND=rk) < xf_max) THEN
         !
      !    amp_x=coeff*EXP(- (REAL(xx-nxp/2,KIND=rk)*xl/REAL(nxp,KIND=rk) )**2/sigma)
         !
      !    hh(xx,yy,zz,1)=hh(xx,yy,zz,1) + fu(yy,zz)*amp_x
      !    hh(xx,yy,zz,2)=hh(xx,yy,zz,2) + fv(yy,zz)*amp_x
      !    hh(xx,yy,zz,3)=hh(xx,yy,zz,3) + fw(yy,zz)*amp_x
         !
      !    ENDIF

  ENDDO XL20
  ENDDO YL20
  ENDDO ZL20

!!!   CFL




   u_max=maxval(uu(:,:,:,1))
   v_max=maxval(uu(:,:,:,2))
   w_max=maxval(uu(:,:,:,3))
   dx=xl/REAL(nxp,KIND=rk)
   dy=yl/REAL(nyp,KIND=rk)
   dz=zl/REAL(nzp,KIND=rk)
   print *,dx,dy,dz
   CFL = u_max*(dt/dx)+v_max*(dt/dy)+w_max*(dt/dz)
   CFL = CFL*pi
   print *,'CFL =',CFL,u_max,v_max,w_max


   CALL F_FFT(hh,hh_C)


 ZL30 : DO zz=1,nz
 YL30 :    DO yy=1,ny
            xx=1
              IF(yy==1 .AND. zz==1) THEN
                hh_C(xx,yy,zz,1)=0._rk
                hh_C(xx,yy,zz,2)=0._rk
                hh_C(xx,yy,zz,3)=0._rk
              ELSE
                ac1=hh_C(xx,yy,zz,1)
                ac2=hh_C(xx,yy,zz,2)
                ac3=hh_C(xx,yy,zz,3)

                k_quad=kx(xx)**2+ky(yy)**2+kz(zz)**2
                k_quad=1./k_quad

                div=kx(xx)*hh_C(xx,yy,zz,1)&
                   +ky(yy)*hh_C(xx,yy,zz,2)&
                   +kz(zz)*hh_C(xx,yy,zz,3)


                hh_C(xx,yy,zz,1)=ac1-div*kx(xx)*k_quad
                hh_C(xx,yy,zz,2)=ac2-div*ky(yy)*k_quad
                hh_C(xx,yy,zz,3)=ac3-div*kz(zz)*k_quad
              ENDIF

 XL31 :       DO xx=2,nx/2

              ac1=hh_C(xx,yy,zz,1)
              ac2=hh_C(xx,yy,zz,2)
              ac3=hh_C(xx,yy,zz,3)

               k_quad=kx(xx)**2+ky(yy)**2+kz(zz)**2
               k_quad=1./k_quad

               div=kx(xx)*hh_C(xx,yy,zz,1)&
                  +ky(yy)*hh_C(xx,yy,zz,2)&
                  +kz(zz)*hh_C(xx,yy,zz,3)

                  hh_C(xx,yy,zz,1)=ac1-div*kx(xx)*k_quad
                  hh_C(xx,yy,zz,2)=ac2-div*ky(yy)*k_quad
                  hh_C(xx,yy,zz,3)=ac3-div*kz(zz)*k_quad


 ENDDO XL31
 ENDDO YL30
 ENDDO ZL30

END SUBROUTINE nonlinear

!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
SUBROUTINE linear
 implicit none
 INTEGER(KIND=ik)               :: zz,yy,xx,xr,xi
 COMPLEX(KIND=rk)               :: k_quad
 REAL(KIND=rk)                  :: den,pnrk,qnrk

 ! n_k=i_rk(mod(it,rk_steps))

 pnrk=2_rk*Re/(dt*ark(n_k,rk_steps)+dt*brk(n_k,rk_steps))

 qnrk=dt*ark(n_k,rk_steps)*pnrk


 ZL10 : DO zz=1,nz
 YL10 : DO yy=1,ny
 XL10 : DO xx=1,nx/2
               k_quad=kx(xx)**2+ky(yy)**2+kz(zz)**2
               den=1./(pnrk-k_quad)
               uu_C(xx,yy,zz,:)=(puu_C(xx,yy,zz,:)+qnrk*hh_C(xx,yy,zz,:))*den
 ENDDO XL10
 ENDDO YL10
 ENDDO ZL10

END SUBROUTINE linear
!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
END MODULE time_advancement_mod

!!!.....................................................................
!!!.....................................................................
!!!.....................................................................
!!TODO k_quad,k1,k2,k3 vanno definiti come complessi, linear e prhs vanno cambiati
!!  per operazioni con numeri complessi, rename a1,a2,a3,ss complessi
