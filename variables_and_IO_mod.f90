!!!.....................................................................
!!!..........................VARIABLES AND IO MOD.......................
!!!.....................................................................
MODULE variables_and_IO_mod

  USE parameters_mod
 !~  USE, INTRINSIC :: iso_c_binding
  implicit none

!!Grid forcing variables
  ! REAL(KIND=rk),DIMENSION(:,:),ALLOCATABLE       :: fu_prev,fv_prev,fw_prev      !forcings
  ! REAL(KIND=rk),DIMENSION(:,:),ALLOCATABLE       :: fu_next,fv_next,fw_next
  ! REAL(KIND=rk),DIMENSION(:,:),ALLOCATABLE       :: fu,fv,fw
  REAL(KIND=rk)                                  :: f_amp
  REAL(KIND=rk)                                  :: t_forz
  REAL(KIND=rk)                                  :: thick
  REAL(KIND=RK),DIMENSION(2)                     :: t_weight

  INTEGER(KIND=ik)                               :: nc
  INTEGER(KIND=ik)                               :: nforz
  CHARACTER(LEN=10)                              :: file_gr_rest

!!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
  !! Problem variables
   INTEGER(KIND=ik)                               :: nx,ny,nz   !fourier modes
   INTEGER(KIND=ik)                               :: nxp,nyp,nzp   !real points
   REAL(KIND=rk)                                  :: xl,yl,zl   !domain's dimensions
   COMPLEX(KIND=rk),DIMENSION(:),ALLOCATABLE      :: kx,ky,kz
   COMPLEX(C_DOUBLE_COMPLEX), dimension(:,:,:,:),ALLOCATABLE     ::puu_C
 ! COMPLEX(C_DOUBLE_COMPLEX), dimension(:,:,:,:),ALLOCATABLE      ::pcc_C
   REAL(KIND=rk)                                  :: sigma,al

  !!time variables
   INTEGER(KIND=ik)                               :: it,itmin,itmax
   INTEGER(KIND=ik)                               :: it_out,it_wrt
   INTEGER(KIND=ik)                               :: it_stat
   REAL(KIND=rk)                                  :: t,dt
  !!problem variables
   INTEGER(KIND=ik),DIMENSION(12)                 :: iseed
   INTEGER(KIND=ik)                               :: rk_steps
   INTEGER(KIND=ik)                               :: i_couple
   REAL(KIND=rk)                                  :: Re
   REAL(KIND=rk)                                  :: Pe
   REAL(KIND=rk)                                  :: De
   REAL(KIND=rk)                                  :: eta_p

   !!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
   !addresses and characters

   CHARACTER(LEN=4)                              :: ver
   CHARACTER(LEN=40)                             :: path
   CHARACTER(LEN=10)                             :: file_bin

   !!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .

 CONTAINS
   !.......................................................................
   SUBROUTINE read_input_parameters
   !read case parameters and allocates all the variables
      OPEN(unit=1,file='sim_par.dat',status='unknown')
       READ(1,*)  ver
       READ(1,*)  file_bin
       READ(1,*)  file_gr_rest
       READ(1,*)  path
       READ(1,*)  xl
       READ(1,*)  yl
       READ(1,*)  zl
       READ(1,*)  nx
       READ(1,*)  ny
       READ(1,*)  nz
       READ(1,*)  Re
       READ(1,*)  itmin
       READ(1,*)  itmax
       READ(1,*)  it_out
       READ(1,*)  it_stat
       READ(1,*)  dt
       READ(1,*)  rk_steps
       READ(1,*)  f_amp
       READ(1,*)  iseed(1)
       READ(1,*)  t_forz
       READ(1,*)  nc
       READ(1,*)  thick
       READ(1,*)  i_couple
      CLOSE(1)
   !multiply domain dimensions by pi
       xl=xl*pi; yl=yl*pi; zl=zl*pi
       al=0.
   !assigns n.er of real points in case the dealiasing is on (al=1)
      IF (rk_steps /= 3 .AND. rk_steps /= 4) rk_steps=3
       nxp=nx+al*nx/2 ; nyp=ny+al*ny/2 ; nzp=nz+al*nz/2
      !  itmin=itmin*rk_steps
      !  itmax=itmax*rk_steps
      !  it_out=it_out*rk_steps
      !  it_stat=it_stat*rk_steps
      !  nforz=rk_steps*int(t_forz/dt,KIND=4)

       nforz=int(t_forz/dt,KIND=4)
       it_wrt=itmin+it_out

         write(*,*) 'initial parameters'
         write(*,*) 'xl,yl,zl:',xl,yl,zl
         write(*,*) 'Re:',Re
         write(*,*) 'dt',dt
         write(*,*) 'f_amp',f_amp
         write(*,*) 'sigma',sigma
         write(*,*) 'iseed',iseed(1)
         write(*,*)  'it mx/min', itmax,'/',itmin
         write(*,*)  'it out', it_wrt
         write(*,*) 'n_step', rk_steps
         write(*,*) 'n_forz', nforz
         write(*,*) 'n_celle', nc
         write(*,*) 'thcik', thick
         write(*,*) 'i_couple',i_couple
   !allocation of variables
!! De-aliasing?

END SUBROUTINE read_input_parameters

!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
SUBROUTINE memory_initialization
  ALLOCATE(puu_C(1:nx/2,1:ny,1:nz,1:3))
!  ALLOCATE(pcc_C(1:nx/2,1:ny,1:nz,1:3))
 ALLOCATE(kx(1:nx/2))
 ALLOCATE(ky(1:ny))
 ALLOCATE(kz(1:nz))

END SUBROUTINE memory_initialization
!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .


!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
SUBROUTINE free_memory
 DEALLOCATE(puu_C)
 ! DEALLOCATE(pcc_C)
 DEALLOCATE(kx)
 DEALLOCATE(ky)
 DEALLOCATE(kz)

END SUBROUTINE free_memory
!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

!!!. . . . . . . . . . WAVE NUMBERS. . . . . . . . . . . . . . . . . . .

SUBROUTINE wave_numbers
 implicit none
 INTEGER(KIND=ik)                               ::ii,jj,kk


 DO ii=1,nx/2
	kx(ii)=CMPLX(0._rk,(2._rk*pi/xl)*REAL(ii-1,KIND=rk))
 ENDDO

 ky(1)=(0._rk,0._rk)
 DO jj=2,ny/2+1
	ky(jj)      = CMPLX(0._rk, (2._rk*pi/yl)*REAL(jj-1,KIND=rk))
	ky(ny+2-jj) = CMPLX(0._rk,-(2._rk*pi/yl)*REAL(jj-1,KIND=rk))
 ENDDO

 kz(1)=(0._rk,0._rk)
 DO kk=2,nz/2+1
	kz(kk)      = CMPLX(0._rk, (2._rk*pi/zl)*REAL(kk-1,KIND=rk))
	kz(nz+2-kk) = CMPLX(0._rk,-(2._rk*pi/zl)*REAL(kk-1,KIND=rk))
 ENDDO

END SUBROUTINE wave_numbers
!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .


END MODULE variables_and_IO_mod
!!!.....................................................................
!!!.....................................................................
!!!.....................................................................
