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
!!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
CONTAINS
SUBROUTINE MPI_initialize
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE (MPI_COMM_WORLD,n_proc,ierr)
      CALL MPI_COMM_RANK (MPI_COMM_WORLD,proc_id,ierr)
      proc0=(proc_id==0)

      CALL MPI_Bcast(nxp,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      CALL MPI_Bcast(nyp,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      CALL MPI_Bcast(nzp,1, MPI_INTEGER,0,mpi_comm_world,ierr)

!! Sets P3DFFT decomposition to 1D (2D not yet implemented)
      dims(1) = 1
      dims(2) = n_proc
      PPRINT("'Number of processes ',n_proc")
!! TODO p3dfft docs
      CALL p3dfft_setup(dims,nxp,nyp,nzp,MPI_COMM_WORLD,nxp,nyp,nzp,.true.)
!!    Gets dimension of local real array
      CALL p3dfft_get_dims(Rstart,Rend,Rsize,1)
!!    Gets dimension of local complex array
      CALL p3dfft_get_dims(Cstart,Cend,Csize,2)

      CALL_BARRIER

!!    Sets C-style pointer for in-place transforms
      pointer_uu = fftw_alloc_complex(int(Rsize(1)*Rsize(2)*Rsize(3)*3,C_SIZE_T))
      pointer_hh = fftw_alloc_complex(int(Rsize(1)*Rsize(2)*Rsize(3)*3,C_SIZE_T))

      CALL c_f_pointer(pointer_uu,uu_C,[Cstart(1):Cend(1),Cstart(2):Cend(2),Cstart(3):Cend(3),3])
      CALL c_f_pointer(pointer_uu,uu,[Rstart(1):Rend(1),Rstart(2):Rend(2),Rstart(3):Rend(3),3])

      CALL c_f_pointer(pointer_hh,hh_C,[Cstart(1):Cend(1),Cstart(2):Cend(2),Cstart(3):Cend(3),3])
      CALL c_f_pointer(pointer_hh,hh,[Rstart(1):Rend(1),Rstart(2):Rend(2),Rstart(3):Rend(3),3])

      CALL_BARRIER

      uu_C=cmplx(0._rk,0._rk)
END SUBROUTINE MPI_initialize
!!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .

END MODULE MPI_mod
!!!.....................................................................
!!!.....................................................................
!!!.....................................................................
