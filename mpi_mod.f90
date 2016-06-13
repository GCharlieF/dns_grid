!!!.....................................................................
!!!..............................MPI MOD................................
!!!.....................................................................
MODULE MPI_mod
 USE p3dfft
 USE parameters_mod
 USE variables_mod
 IMPLICIT NONE
#INCLUDE 'fpp_macros.h'   !!#INCLUDE fpp macros
 INCLUDE 'mpif.h'
 INCLUDE 'fftw3.f03'
!!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
!!FFTW Variables
   INTEGER,SAVE                                           ::plan_type=FFTW_ESTIMATE
   TYPE(C_PTR),SAVE                                       ::plan_ub,plan_uf
  !  TYPE(C_PTR),SAVE                                       ::plan_cb,plan_cf
   REAL(C_DOUBLE),pointer,DIMENSION(:,:,:,:)              ::uu
   REAL(C_DOUBLE),pointer,DIMENSION(:,:,:,:)              ::hh
  !  REAL(C_DOUBLE),pointer,DIMENSION(:,:,:,:)              ::cc
  !  REAL(C_DOUBLE),pointer,DIMENSION(:,:,:,:)              ::hcc
  !  REAL(C_DOUBLE),pointer,DIMENSION(:,:,:,:)              ::app_cc

   COMPLEX(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:,:) ::uu_C
   COMPLEX(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:,:) ::hh_C
  !  COMPLEX(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:,:) ::cc_C
  !  COMPLEX(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:,:) ::hcc_C
  !  COMPLEX(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:,:) ::app_cc_C

   TYPE(C_PTR)                                            :: pointer_uu
   TYPE(C_PTR)                                            :: pointer_hh
  !  TYPE(C_PTR)                                            :: pointer_cc
  !  TYPE(C_PTR)                                            :: pointer_hcc
  !  TYPE(C_PTR)                                            :: pointer_app_cc
!!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
!!P3DFFT Variables
INTEGER(KIND=IK),DIMENSION(3)                             :: Rstart,Rend,Rsize
INTEGER(KIND=IK),DIMENSION(3)                             :: Cstart,Cend,Csize
INTEGER(KIND=IK),DIMENSION(2)                             :: dims
!!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
!!MPI Variables
INTEGER(KIND=IK)                                          :: n_proc
INTEGER(KIND=IK)                                          :: proc_id
INTEGER(KIND=IK)                                          :: ierr
LOGICAL,SAVE                                              :: proc0
INTEGER(KIND=ik),DIMENSION(:),ALLOCATABLE      :: xr_loc,yr_loc,zr_loc
INTEGER(KIND=ik),DIMENSION(:),ALLOCATABLE      :: xc_loc,yc_loc,zc_loc
!!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
CONTAINS
!!!. . . . . . . . . . . MPI INITIALIZE . . . . . . . . . . . . . . . . . . . .
SUBROUTINE MPI_initialize
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE (MPI_COMM_WORLD,n_proc,ierr)
      CALL MPI_COMM_RANK (MPI_COMM_WORLD,proc_id,ierr)
      proc0=(proc_id==0)

      CALL MPI_Bcast(nxp,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      CALL MPI_Bcast(nyp,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      CALL MPI_Bcast(nzp,1, MPI_INTEGER,0,mpi_comm_world,ierr)
END SUBROUTINE MPI_initialize
!!!. . . . . . . . . . . MPI INITIALIZE . . . . . . . . . . . . . . . . . . . .
SUBROUTINE MPI_P3DFFT_initialize
!! Sets P3DFFT decomposition to 1D (2D not yet implemented)
      dims(1) = 1
      dims(2) = n_proc
      PPRINT("Number of processes")
      PPRINT(n_proc)
!! TODO p3dfft docs
      ! CALL p3dfft_setup(dims,nxp,nyp,nzp,MPI_COMM_WORLD,nxp,nyp,nzp,.true.)
      CALL p3dfft_setup(dims,nxp,nyp,nzp,MPI_COMM_WORLD)
!!    Gets dimension of local real array
      CALL p3dfft_get_dims(Rstart,Rend,Rsize,1)
!!    Gets dimension of local complex array
      CALL p3dfft_get_dims(Cstart,Cend,Csize,2)

      CALL_BARRIER

!!    Sets C-style pointer for in-place transforms
      pointer_uu = fftw_alloc_complex(int(Rsize(1)*Rsize(2)*Rsize(3)*3,C_SIZE_T))
      pointer_hh = fftw_alloc_complex(int(Rsize(1)*Rsize(2)*Rsize(3)*3,C_SIZE_T))

      CALL c_f_pointer(pointer_uu,uu_C,[Csize(1),Csize(2),Csize(3),3])
      CALL c_f_pointer(pointer_uu,uu,[Rsize(1),Rsize(2),Rsize(3),3])

      CALL c_f_pointer(pointer_hh,hh_C,[Csize(1),Csize(2),Csize(3),3])
      CALL c_f_pointer(pointer_hh,hh,[Rsize(1),Rsize(2),Rsize(3),3])
!! TODO Test this
!! Sets map with local indices
      ALLOCATE(xr_loc(Rsize(1)))
      ALLOCATE(yr_loc(Rsize(2)))
      ALLOCATE(zr_loc(Rsize(3)))

      ALLOCATE(xc_loc(Csize(1)))
      ALLOCATE(yc_loc(Csize(2)))
      ALLOCATE(zc_loc(Csize(3)))

      xr_loc(:)=(/Rstart(1):Rend(1)/)
      yr_loc(:)=(/Rstart(2):Rend(2)/)
      zr_loc(:)=(/Rstart(3):Rend(3)/)

      xc_loc(:)=(/Cstart(1):Cend(1)/)
      yc_loc(:)=(/Cstart(2):Cend(2)/)
      zc_loc(:)=(/Cstart(3):Cend(3)/)
      CALL_BARRIER
      uu=cmplx(0.,0.)
      hh=cmplx(0.,0.)

      CALL_BARRIER
!!    p3dfft seems to require that first transform is r2c, if removed it tries
!!    to re-allocate variables during first c2r transform
      CALL p3dfft_ftran_r2c_many (uu,Rsize(1)*Rsize(2)*Rsize(3),uu_C, &
                           Csize(1)*Csize(2)*Csize(3),3,'fft')
      CALL p3dfft_ftran_r2c_many (hh,Rsize(1)*Rsize(2)*Rsize(3),hh_C, &
                           Csize(1)*Csize(2)*Csize(3),3,'fft')
!! FIXME remove when IO is implemented
                           uu_C=cmplx(0.,0.)
                           hh_C=cmplx(0.,0.)
END SUBROUTINE MPI_P3DFFT_initialize
!!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .


!!! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
END MODULE MPI_mod
!!!.....................................................................
!!!.....................................................................
!!!.....................................................................
