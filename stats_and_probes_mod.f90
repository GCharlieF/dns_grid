!!!.....................................................................
!!!.......................STATS AND PROBES MOD..........................
!!!.....................................................................

MODULE stats_and_probes_mod
 USE parameters_mod
 USE variables_mod
 USE fft_mod
 USE grid_forcing_mod

IMPLICIT NONE
 REAL(KIND=rk)                           :: KE
 REAL(KIND=rk),DIMENSION(3)              :: KE_ii
 REAL(KIND=rk)                           :: diss
 LOGICAL                                 :: stats_time

CONTAINS
!!! . . . . . . . . . . . . COMPUTE CFL . . . . . . . . . . . . . . . . . . . . . .

SUBROUTINE compute_CFL
 IMPLICIT NONE
 REAL(KIND=rk)                         :: u_max,v_max,w_max
 REAL(KIND=rk)                         :: dx,dy,dz
 REAL(KIND=rk)                         :: CFL

  u_max=maxval(uu(:,:,:,1))
  v_max=maxval(uu(:,:,:,2))
  w_max=maxval(uu(:,:,:,3))

  dx=xl/REAL(nxp,KIND=rk)
  dy=yl/REAL(nyp,KIND=rk)
  dz=zl/REAL(nzp,KIND=rk)
  CFL = u_max*(dt/dx)+v_max*(dt/dy)+w_max*(dt/dz)
  CFL = CFL*pi
  PRINT *,'- - - - - - - -(compute cfl) - - - - - - - - - -'
  PRINT *,'CFL    ::',CFL
  PRINT *,'u max  ::',u_max
  PRINT *,'v max  ::',v_max
  PRINT *,'w max  ::',w_max
  PRINT *,'- - - - - - - - - - - - - - - - - - - - - - - - -'

END SUBROUTINE compute_CFL
!!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

!!! . . . . . . . . . . . . . . ENERGY . . . . . . . . . . . . . . . . .
SUBROUTINE average_energy(n_k)

IMPLICIT NONE
 INTEGER(KIND=ik)                        :: xx,yy,zz,ii,n_k
 IF (stats_time) THEN
 ZL10 :DO zz=1,nzp
 YL10 :DO yy=1,nyp
 XL10 :DO xx=1,nxp
       DO ii=1,3
       KE=KE+0.5_rk*uu(xx,yy,zz,ii)**2
       KE_ii(ii)=KE_ii(ii)+0.5_rk*uu(xx,yy,zz,ii)**2
       ENDDO
 ENDDO XL10
 ENDDO YL10
 ENDDO ZL10
 KE=KE/REAL(nxp*nyp*nzp,KIND=rk)
 KE_ii=KE_ii/REAL(nxp*nyp*nzp,KIND=rk)
 PRINT *,'- - - - - - - -(average energy) - - - - - - - - - -'

 PRINT *,'KE     ::',KE
 PRINT *,'u^2/2  ::',KE_ii(1)
 PRINT *,'v^2/2  ::',KE_ii(2)
 PRINT *,'w^2/2  ::',KE_ii(3)
 PRINT *,'- - - - - - - - - - - - - - - - - - - - - - - - - -'
 ENDIF

!! Stats output
IF (MOD(it-1,it_stat)==0 .AND. n_k==1) THEN
  WRITE(21,2)t,KE,KE_ii(1),KE_ii(2),KE_ii(3)
ENDIF
2 FORMAT(5(e15.6,1x))

END SUBROUTINE average_energy
!!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

!!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
SUBROUTINE average_dissipiation
  IMPLICIT NONE
   INTEGER(KIND=ik)                        :: xx,yy,zz,ii,n_k
!!TODO dichiara dui_dxi in fft_mod
!!TODO controlla che la trasformata funzioni cosi,
!!TODO se necessario definisic nuovo tipo di trasformata
 ZL10: DO zz=1,nz
 YL10: DO zz=1,ny
 XL10: DO zz=1,nx
   dui_dxi_C(xx,yy,zz,1)=uu_C(xx,yy,zz,1)*kx(xx)
   dui_dxi_C(xx,yy,zz,2)=uu_C(xx,yy,zz,1)*ky(yy)
   dui_dxi_C(xx,yy,zz,3)=uu_C(xx,yy,zz,1)*kz(zz)
   dui_dxi_C(xx,yy,zz,4)=uu_C(xx,yy,zz,2)*kx(xx)
   dui_dxi_C(xx,yy,zz,5)=uu_C(xx,yy,zz,2)*ky(yy)
   dui_dxi_C(xx,yy,zz,6)=uu_C(xx,yy,zz,2)*kz(zz)
   dui_dxi_C(xx,yy,zz,7)=uu_C(xx,yy,zz,3)*kx(xx)
   dui_dxi_C(xx,yy,zz,8)=uu_C(xx,yy,zz,3)*ky(yy)
   dui_dxi_C(xx,yy,zz,9)=uu_C(xx,yy,zz,3)*kz(zz)
 ENDDO XL10
 ENDDO YL10
 ENDDO ZL10
 CALL FFT_B(dui_dxi_C(:,:,:,1:3),dui_dxi(:,:,:,1:3))
 CALL FFT_B(dui_dxi_C(:,:,:,4:6),dui_dxi(:,:,:,4:6))
 CALL FFT_B(dui_dxi_C(:,:,:,7:9),dui_dxi(:,:,:,7:9))

 diss=0._rk
 ZL20: DO zz=1,nzp
 YL20: DO yy=1,nyp
 XL20: DO xx=1,nxp
   diss=diss
 ENDDO XL20
 ENDDO YL20
 ENDDO ZL20
END SUBROUTINE average_dissipiation

!!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .


END MODULE stats_and_probes_mod
!!!.....................................................................
!!!.....................................................................
