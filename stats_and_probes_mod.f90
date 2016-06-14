!!!.....................................................................
!!!.......................STATS AND PROBES MOD..........................
!!!.....................................................................

MODULE stats_and_probes_mod
 USE parameters_mod
 USE variables_mod
 ! USE fft_mod
 USE MPI_mod
 USE grid_forcing_mod

IMPLICIT NONE
#INCLUDE 'fpp_macros.h'
 REAL(KIND=rk)                           :: KE
 REAL(KIND=rk),DIMENSION(3)              :: KE_ii
 REAL(KIND=rk)                           :: CFL
 REAL(KIND=rk)                           :: diss
 REAL(KIND=rk)                           :: dt_new
 LOGICAL                                 :: stats_time

CONTAINS
!!! . . . . . . . . . . . . COMPUTE CFL . . . . . . . . . . . . . . . . . . . . . .

SUBROUTINE STATS_compute_CFL
 IMPLICIT NONE
 REAL(KIND=rk)                         :: u_max,v_max,w_max
 REAL(KIND=rk)                         :: u_max_loc,v_max_loc,w_max_loc
 REAL(KIND=rk)                         :: dx,dy,dz
 REAL(KIND=rk)                         :: CFL_max

  CFL_max=0.15
  u_max_loc=maxval(uu(:,:,:,1))
  v_max_loc=maxval(uu(:,:,:,2))
  w_max_loc=maxval(uu(:,:,:,3))

  dx=xl/REAL(nxp,KIND=rk)
  dy=yl/REAL(nyp,KIND=rk)
  dz=zl/REAL(nzp,KIND=rk)
  CALL MPI_ALLREDUCE(u_max_loc,u_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX,mpi_comm_world,ierr)
  CALL MPI_ALLREDUCE(v_max_loc,v_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX,mpi_comm_world,ierr)
  CALL MPI_ALLREDUCE(w_max_loc,w_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX,mpi_comm_world,ierr)

  CFL = u_max*(dt/dx)+v_max*(dt/dy)+w_max*(dt/dz)
  CFL = CFL*pi
!! FIXME check se la u_max locale va bene o se bisogna broadcastarla

  ! CALL MPI_REDUCE(CFL_loc, CFL, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi_comm_world,ierr)
!! FIXME dt_new sees the global maxima
  dt_new=CFL_max/(u_max/dx+v_max/dy+w_max/dz)/pi


  MASTER PRINT *,'- - - - - - - -(compute cfl) - - - - - - - - - -'
  MASTER PRINT *,'CFL    ::',CFL
  MASTER PRINT *,'u max  ::',u_max
  MASTER PRINT *,'v max  ::',v_max
  MASTER PRINT *,'w max  ::',w_max
  MASTER PRINT *,'- - - - - - - - - - - - - - - - - - - - - - - - -'
  ! IF (CFL > 3.0)  STOP
END SUBROUTINE STATS_compute_CFL
!!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

!!! . . . . . . . . . . . . . . ENERGY . . . . . . . . . . . . . . . . .
SUBROUTINE STATS_average_energy(compute_stats)

IMPLICIT NONE
 INTEGER(KIND=ik)                        :: xx,yy,zz,ii
 REAL(KIND=rk)                           :: KE_loc
 REAL(KIND=rk),DIMENSION(3)              :: KE_ii_loc
 LOGICAL                                 :: compute_stats
 IF (compute_stats) THEN
 KE_loc=0.
 KE_ii_loc=0.
 MASTER KE=0.
 MASTER KE_ii=0.
 ZL10 :DO zz=1,Rsize(3)
 YL10 :DO yy=1,Rsize(2)
 XL10 :DO xx=1,Rsize(1)
       DO ii=1,3
       KE_loc=KE_loc+0.5_rk*uu(xx,yy,zz,ii)**2
       KE_ii_loc(ii)=KE_ii_loc(ii)+0.5_rk*uu(xx,yy,zz,ii)**2
       ENDDO
 ENDDO XL10
 ENDDO YL10
 ENDDO ZL10
 KE_loc=KE_loc/REAL(nxp*nyp*nzp,KIND=rk)
 KE_ii_loc=KE_ii_loc/REAL(nxp*nyp*nzp,KIND=rk)

 CALL MPI_ALLREDUCE(KE_loc, KE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world,ierr)
 CALL MPI_ALLREDUCE(KE_ii_loc(1), KE_ii(1), 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world,ierr)
 CALL MPI_ALLREDUCE(KE_ii_loc(2), KE_ii(2), 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world,ierr)
 CALL MPI_ALLREDUCE(KE_ii_loc(3), KE_ii(3), 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world,ierr)


 MASTER PRINT *,'- - - - - - - -(average energy) - - - - - - - - - -'

 MASTER PRINT *,'KE     ::',KE,KE_loc
 MASTER PRINT *,'u^2/2  ::',KE_ii(1)
 MASTER PRINT *,'v^2/2  ::',KE_ii(2)
 MASTER PRINT *,'w^2/2  ::',KE_ii(3)
 MASTER PRINT *,'- - - - - - - - - - - - - - - - - - - - - - - - - -'

ENDIF
!! Stats output
IF (compute_stats) THEN
  MASTER WRITE(21,2)t,KE,KE_ii(1),KE_ii(2),KE_ii(3)
ENDIF
2 FORMAT(5(e15.6,1x))

END SUBROUTINE STATS_average_energy
!!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

!!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! SUBROUTINE average_dissipiation
!   IMPLICIT NONE
!    INTEGER(KIND=ik)                        :: xx,yy,zz,ii,n_k
! !!TODO dichiara dui_dxi in fft_mod
! !!TODO controlla che la trasformata funzioni cosi,
! !!TODO se necessario definisic nuovo tipo di trasformata
!  ZL10: DO zz=1,nz
!  YL10: DO yy=1,ny
!  XL10: DO xx=1,nx
!    dui_dxi_C(xx,yy,zz,1)=uu_C(xx,yy,zz,1)*kx(xx)
!    dui_dxi_C(xx,yy,zz,2)=uu_C(xx,yy,zz,1)*ky(yy)
!    dui_dxi_C(xx,yy,zz,3)=uu_C(xx,yy,zz,1)*kz(zz)
!    dui_dxi_C(xx,yy,zz,4)=uu_C(xx,yy,zz,2)*kx(xx)
!    dui_dxi_C(xx,yy,zz,5)=uu_C(xx,yy,zz,2)*ky(yy)
!    dui_dxi_C(xx,yy,zz,6)=uu_C(xx,yy,zz,2)*kz(zz)
!    dui_dxi_C(xx,yy,zz,7)=uu_C(xx,yy,zz,3)*kx(xx)
!    dui_dxi_C(xx,yy,zz,8)=uu_C(xx,yy,zz,3)*ky(yy)
!    dui_dxi_C(xx,yy,zz,9)=uu_C(xx,yy,zz,3)*kz(zz)
!  ENDDO XL10
!  ENDDO YL10
!  ENDDO ZL10
!  CALL FFT_B(dui_dxi_C(:,:,:,1:3),dui_dxi(:,:,:,1:3))
!  CALL FFT_B(dui_dxi_C(:,:,:,4:6),dui_dxi(:,:,:,4:6))
!  CALL FFT_B(dui_dxi_C(:,:,:,7:9),dui_dxi(:,:,:,7:9))
!
!  diss=0._rk
!  ZL20: DO zz=1,nzp
!  YL20: DO yy=1,nyp
!  XL20: DO xx=1,nxp
!    diss=diss
!  ENDDO XL20
!  ENDDO YL20
!  ENDDO ZL20
! END SUBROUTINE average_dissipiation

!!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .


END MODULE stats_and_probes_mod
!!!.....................................................................
!!!.....................................................................
