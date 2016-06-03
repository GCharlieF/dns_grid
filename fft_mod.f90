!!!.....................................................................
!!!..........................FFT MOD....................................
!!!.....................................................................

MODULE fft_mod
 USE parameters_mod
 USE variables_mod

 implicit none
 INCLUDE 'fftw3.f03'

 !!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
 !!FFTW variables
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

   INTEGER(KIND=ik)                                       :: rank,howmany,idist,odist,istride,ostride
   INTEGER(KIND=ik),DIMENSION(3)                          :: numb,inembed,onembed
   INTEGER(KIND=ik)                                       :: iret
 !!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .

  contains

 SUBROUTINE FFT_initialization

!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


!!~   Allocates C vector for the velocity field
!    pointer_uu = fftw_alloc_complex(int((nxp/2+1)*nyp*nzp*3,C_SIZE_T))
!    pointer_hh = fftw_alloc_complex(int((nxp/2+1)*nyp*nzp*3,C_SIZE_T))
!
! !!~   Generates C pointers uuC and uu pointing at the same memory space
!    CALL c_f_pointer(pointer_uu,uu_C,[nxp/2+1,nyp,nzp,3])
!    CALL c_f_pointer(pointer_uu,uu,[nxp+2,nyp,nzp,3])
!
!    CALL c_f_pointer(pointer_hh,hh_C,[nxp/2+1,nyp,nzp,3])
!    CALL c_f_pointer(pointer_hh,hh,[nxp+2,nyp,nzp,3])

!!    FFTW parameters
   rank=3
   numb(1)=nzp
   numb(2)=nyp
   numb(3)=nxp
   howmany=3
   istride=1
   idist=(nxp/2+1)*nyp*nzp
   ostride=1
   odist=(nxp+2)*nyp*nzp
   inembed(1)=nzp
   inembed(2)=nyp
   inembed(3)=nxp/2+1
   onembed(1)=nzp
   onembed(2)=nyp
   onembed(3)=nxp+2

!     Initialize fft with multithread
!~  	 CALL dfftw_init_threads(iret)
!~ 	     CALL dfftw_plan_with_nthreads(2)
!    Create plan for the FFT

    plan_ub=fftw_plan_many_dft_c2r(rank,numb,howmany,uu_C,&
            inembed,istride,idist,uu,onembed,ostride,odist,plan_type)

    plan_uf=fftw_plan_many_dft_r2c(rank,numb,howmany,uu,&
            onembed,ostride,odist,uu_C,inembed,istride,idist,plan_type)


!!   Allocates C vector for the velocity field
!     pointer_cc     = fftw_alloc_complex(int((nxp/2+1)*nyp*nzp*3,C_SIZE_T))
!     pointer_app_cc = fftw_alloc_complex(int((nxp/2+1)*nyp*nzp*3,C_SIZE_T))
!
! !!~   Generates C ers uuC and uu pointing at the same memory space
!     CALL c_f_pointer(pointer_cc,cc_C,[nxp/2+1,nyp,nzp,3])
!     CALL c_f_pointer(pointer_cc,cc,[nxp+2,nyp,nzp,3])
!
!     CALL c_f_pointer(pointer_hcc,hcc_C,[nxp/2+1,nyp,nzp,3])
!     CALL c_f_pointer(pointer_hcc,hcc,[nxp+2,nyp,nzp,3])
!
!     CALL c_f_pointer(pointer_app_cc,app_cc_C,[nxp/2+1,nyp,nzp,3])
!     CALL c_f_pointer(pointer_app_cc,app_cc,[nxp+2,nyp,nzp,3])
!
! !!  FFTW parameters
!     rank=3
!     numb(1)=nzp
!     numb(2)=nyp
!     numb(3)=nxp
!     howmany=6
!     istride=1
!     idist=(nxp/2+1)*nyp*nzp
!     ostride=1
!     odist=(nxp+2)*nyp*nzp
!     inembed(1)=nzp
!     inembed(2)=nyp
!     inembed(3)=nxp/2+1
!     onembed(1)=nzp
!     onembed(2)=nyp
!     onembed(3)=nxp+2
!
! !     Initialize fft with multithread
! !   	 CALL dfftw_init_threads(iret)
! ! 	   CALL dfftw_plan_with_nthreads(2)
! !!    Create plan for the FFT
!
!      plan_cb=fftw_plan_many_dft_c2r(rank,numb,howmany,cc_C,&
!              inembed,istride,idist,cc,onembed,ostride,odist,plan_type)
!
!      plan_cf=fftw_plan_many_dft_r2c(rank,numb,howmany,cc,&
!              onembed,ostride,odist,cc_C,inembed,istride,idist,plan_type)
 END SUBROUTINE FFT_initialization

 !!!. . . . . . . . . . BACKWARD FFT VELOCITY FIELD. . . . . . . . . . . .
SUBROUTINE B_FFT(vf_C,vf)
 implicit none
 REAL(C_DOUBLE), DIMENSION(:,:,:,:)                  :: vf
 COMPLEX(C_DOUBLE_COMPLEX),DIMENSION(:,:,:,:)        :: vf_C


!!!.   .    .   .    .   .    .   .    .   .    .   .    .   .    .   .
      !!De-aliasing
      IF (al==1) vf_C(nx/2+1:nxp/2,ny/2+1:ny/2+ny/2,nz/2+1:nz/2+nz/2,:)=CMPLX(0._RK,0._RK)

      CALL fftw_execute_dft_c2r(plan_ub,vf_C,vf)

END SUBROUTINE B_FFT
!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!!!. . . . . . . . . . FORWARD FFT VELOCITY FIELD. . . . . . . . . . . . . . .
SUBROUTINE F_FFT(vf,vf_C)
 implicit none
 REAL(C_DOUBLE), DIMENSION(:,:,:,:)                  :: vf
 COMPLEX(C_DOUBLE_COMPLEX),DIMENSION(:,:,:,:)        :: vf_C


!!!.   .    .   .    .   .    .   .    .   .    .   .    .   .    .   .

      CALL fftw_execute_dft_r2c(plan_uf,vf,vf_C)
      vf_C=vf_C/REAL(nxp*nyp*nzp,KIND=rk)

END SUBROUTINE F_FFT
!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

! !!!. . . . . . . . . . BACKWARD FFT CONFORMATION TENSOR . . . . . . . . .
! SUBROUTINE B_FFT_CT(ct_C,ct)
! implicit none
! REAL(C_DOUBLE), DIMENSION(:,:,:,:)          :: ct
! COMPLEX(C_DOUBLE),DIMENSION(:,:,:,:)        :: ct_C
!
!
! !!!.   .    .   .    .   .    .   .    .   .    .   .    .   .    .   .
!
!      CALL fftw_execute_dft_c2r(plan_cb,ct_C,ct)
! END SUBROUTINE B_FFT_CT
! !!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! !!!. . . . . . . . . . FORWARD FFT CONFORMATION TENSOR . . . . . . . . .
! SUBROUTINE F_FFT_CT(ct,ct_C)
! implicit none
! REAL(C_DOUBLE), DIMENSION(:,:,:,:)          :: ct
! COMPLEX(C_DOUBLE),DIMENSION(:,:,:,:)        :: ct_C
!
!
! !!!.   .    .   .    .   .    .   .    .   .    .   .    .   .    .   .
!
!      CALL fftw_execute_dft_r2c(plan_cf,ct,ct_C)
!      ct_C=ct_C/REAL(nxp*nyp*nzp,KIND=rk)
!
! END SUBROUTINE F_FFT_CT
! !!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
SUBROUTINE free_FFT

     CALL fftw_free(pointer_uu)
     CALL fftw_free(pointer_hh)
    !  CALL fftw_free(pointer_cc)
    !  CALL fftw_free(pointer_hcc)

END SUBROUTINE free_FFT
!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

END MODULE fft_mod
!!!.....................................................................
!!!.....................................................................
