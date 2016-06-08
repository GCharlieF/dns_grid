!!!.....................................................................
!!!..........................VARIABLES MOD..............................
!!!.....................................................................
MODULE variables_mod

USE parameters_mod
IMPLICIT NONE

!!Grid forcing variables
INTEGER(KIND=ik)                               :: nc
INTEGER(KIND=ik)                               :: dt_forc
CHARACTER(LEN=10)                              :: file_gr_rest
REAL(KIND=rk),DIMENSION(:,:),ALLOCATABLE       :: fu_prev,fv_prev,fw_prev      !forcings
REAL(KIND=rk),DIMENSION(:,:),ALLOCATABLE       :: fu_next,fv_next,fw_next
REAL(KIND=rk),DIMENSION(:,:),ALLOCATABLE       :: fu,fv,fw
REAL(KIND=rk)                                  :: f_amp
REAL(KIND=rk)                                  :: t_forz
REAL(KIND=rk)                                  :: thick
REAL(KIND=RK),DIMENSION(2)                     :: t_weight


!!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
!! Problem variables
INTEGER(KIND=ik)                               :: nx,ny,nz   !fourier modes
INTEGER(KIND=ik)                               :: nxp,nyp,nzp   !real points
REAL(KIND=rk)                                  :: xl,yl,zl   !domain's dimensions
COMPLEX(KIND=rk),DIMENSION(:),ALLOCATABLE      :: kx,ky,kz
INTEGER(KIND=ik),DIMENSION(:),ALLOCATABLE      :: day,daz
COMPLEX(C_DOUBLE_COMPLEX), dimension(:,:,:,:),ALLOCATABLE     ::puu_C
! COMPLEX(C_DOUBLE_COMPLEX), dimension(:,:,:,:),ALLOCATABLE      ::pcc_C
REAL(KIND=rk)                                  :: sigma,al

INTEGER(KIND=ik),DIMENSION(12)                 :: seed
INTEGER(KIND=ik)                               :: rk_steps
INTEGER(KIND=ik)                               :: i_couple
REAL(KIND=rk)                                  :: Re
REAL(KIND=rk)                                  :: eta_p
REAL(KIND=rk)                                  :: aa
REAL(KIND=rk)                                  :: Pe
REAL(KIND=rk)                                  :: De

!!time variables
INTEGER(KIND=ik)                               :: it,itmin,itmax
INTEGER(KIND=ik)                               :: it_out,it_wrt
INTEGER(KIND=ik)                               :: it_stat
REAL(KIND=rk)                                  :: t,dt

!!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
!!addresses and characters

CHARACTER(LEN=4)                              :: ver
CHARACTER(LEN=40)                             :: path
CHARACTER(LEN=10)                             :: file_bin

!!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .

 CONTAINS
!.......................................................................
SUBROUTINE VAR_read_input_parameters
!!read case parameters and allocates all the variables
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
       READ(1,*)  seed(1)
       READ(1,*)  t_forz
       READ(1,*)  nc
       READ(1,*)  thick
       READ(1,*)  aa
       READ(1,*)  al
       !  READ(1,*)  i_couple

CLOSE(1)
!!multiply domain dimensions by pi
 xl=xl*pi; yl=yl*pi; zl=zl*pi
 aa=aa*pi
!! set rk_steps to 3 if not properly defined
IF (rk_steps /= 3 .AND. rk_steps /= 4) rk_steps=3

!!assigns n.er of real points in case the dealiasing is on (al=1)
nxp=nx+al*nx/2 ; nyp=ny+al*ny/2 ; nzp=nz+al*nz/2

dt_forc=int(t_forz/dt,KIND=ik)
it_wrt=itmin+it_out

 write(*,*) 'initial parameters'
 write(*,*) 'xl,yl,zl: ',xl,yl,zl
 write(*,*) 'Fourier modes: ',nx,ny,nz
 write(*,*) 'De-aliased Fourier modes:',nxp,nyp,nzp
 write(*,*) 'Re:',Re
 write(*,*) 'dt',dt
 write(*,*) 'f amp',f_amp
 write(*,*) 'seeding',seed(1)
 write(*,*) 'it mx/min', itmax,'/',itmin
 write(*,*) 'it out', it_wrt
 write(*,*) 'n_step', rk_steps
 write(*,*) 'n_forz', dt_forc
 write(*,*) 'n_celle', nc
 write(*,*) 'thcik and aa', thick,aa
        !  write(*,*) 'i_couple',i_couple

END SUBROUTINE VAR_read_input_parameters

!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
SUBROUTINE VAR_memory_initialization
  ALLOCATE(puu_C(1:nx/2,1:ny,1:nz,1:3)) !! FIXME rimuovere +1
!  ALLOCATE(pcc_C(1:nx/2,1:ny,1:nz,1:3))
 ALLOCATE(kx(1:nx/2))  !! FIXME rimuovere +1
 ALLOCATE(ky(1:ny))
 ALLOCATE(kz(1:nz))
 ALLOCATE(day(1:ny))
 ALLOCATE(daz(1:nz))
 !allocation of forcing variables
 ALLOCATE(fu(1:nyp,1:nzp))
 ALLOCATE(fv(1:nyp,1:nzp))
 ALLOCATE(fw(1:nyp,1:nzp))
 ALLOCATE(fu_prev(1:nyp,1:nzp))
 ALLOCATE(fv_prev(1:nyp,1:nzp))
 ALLOCATE(fw_prev(1:nyp,1:nzp))
 ALLOCATE(fu_next(1:nyp,1:nzp))
 ALLOCATE(fv_next(1:nyp,1:nzp))
 ALLOCATE(fw_next(1:nyp,1:nzp))
END SUBROUTINE VAR_memory_initialization
!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .


!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
SUBROUTINE VAR_free_memory
 DEALLOCATE(puu_C)
 ! DEALLOCATE(pcc_C)
 DEALLOCATE(kx)
 DEALLOCATE(ky)
 DEALLOCATE(kz)

END SUBROUTINE VAR_free_memory
!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

!!!. . . . . . . . . . WAVE NUMBERS. . . . . . . . . . . . . . . . . . .

SUBROUTINE VAR_wave_numbers
 implicit none
 INTEGER(KIND=ik)                               ::ii,jj,kk

 kx(1)=CMPLX(0._rk,0._rk)
 DO ii=2,nx/2
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

 ! write(18,*),kx
 ! write(19,*),ky
 ! write(20,*),kz
 ! stop
END SUBROUTINE VAR_wave_numbers
!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
SUBROUTINE VAR_dealiased_indeces
 implicit none
 INTEGER(KIND=ik)                               ::ii
 DO ii=1,ny/2
   day(ii)=ii
  !  print *,ii,day(ii)
 ENDDO
 DO ii=1,nz/2
   daz(ii)=ii
 ENDDO
 DO ii=ny/2+1,ny
   day(ii)=ii+al*ny/2
    ! print *,ii,day(ii)
 ENDDO
 ! stop
 DO ii=nz/2+1,nz
   daz(ii)=ii+al*nz/2
 ENDDO
END SUBROUTINE VAR_dealiased_indeces
!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
END MODULE variables_mod
!!!.....................................................................
!!!.....................................................................
!!!.....................................................................
