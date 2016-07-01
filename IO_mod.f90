!!!.....................................................................
!!!...............................IO MOD................................
!!!.....................................................................
MODULE IO_mod
 USE parameters_mod
 USE variables_mod
 USE MPI_mod
 use hdf5
 USE h5util_mod
 implicit none
 INTEGER(KIND=ik),DIMENSION(:),ALLOCATABLE    	  ::ind_x,ind_y,ind_z
REAL(KIND=rk),DIMENSION(:),ALLOCATABLE 		  ::xc,yc,zc
contains


!!!. . . . . . . . . . RE-INDEXING . . . . . . . . . . . . . . . . . . .

SUBROUTINE IO_re_indexing
!!TODO ELIMINARE RE-INDEXING
 INTEGER(KIND=ik)                       :: ii,jj,kk,xx,yy,zz
!!!.   .    .   .    .   .    .   .    .   .    .   .    .   .    .

 ALLOCATE(ind_x(1:nxp))
 ALLOCATE(ind_y(1:nyp))
 ALLOCATE(ind_z(1:nzp))


 ALLOCATE(xc(1:nxp))
 ALLOCATE(yc(1:nyp))
 ALLOCATE(zc(1:nzp))
! . . . . . . . . . . .
	    DO ii=1,nxp/2
         ind_x(ii)=ii
         xx=ind_x(ii)-1
         xc(ii)=xl*REAL(xx-nxp/2,KIND=rk)/REAL(nxp,KIND=rk)

	 ENDDO
     DO ii=nxp/2+1,nxp
         ind_x(ii)=ii
         xx=ind_x(ii)-1
         xc(ii)=xl*REAL(xx-nxp/2,KIND=rk)/REAL(nxp,KIND=rk)
	 ENDDO
! . . . . . . . . . . .
     DO jj=1,nyp/2
         ind_y(jj)=nyp/2+jj
         yy=ind_y(jj)
         yc(jj)=yl*REAL(yy-nyp-1,KIND=rk)/REAL(nyp,KIND=rk)
	 ENDDO
     DO jj=nyp/2+1,nyp
         ind_y(jj)=jj-nyp/2
         yy=ind_y(jj)
         yc(jj)=yl*REAL(yy-1,KIND=rk)/REAL(nyp,KIND=rk)
	 ENDDO
! . . . . . . . . . . .
     DO kk=1,nzp/2
         ind_z(kk)=nzp/2+kk
         zz=ind_z(kk)
         zc(kk)=zl*REAL(zz-nzp-1,KIND=rk)/REAL(nzp,KIND=rk)
	 ENDDO
     DO kk=nzp/2+1,nzp
         ind_z(kk)=kk-nzp/2
         zz=ind_z(kk)
         zc(kk)=zl*REAL(zz-1,KIND=rk)/REAL(nzp,KIND=rk)
	 ENDDO
 ! . . . . . . . . . . .
END SUBROUTINE IO_re_indexing
!!!. . . . . . . . . . . READ FIELD  . . . . . . . . . . . . . . . . . .

SUBROUTINE IO_read_field
!! Reads the velocity in the Fourier space  from input file


implicit none
 INTEGER(KIND=ik)                  :: xx,yy,zz,ii
 INTEGER(KIND=ik)                  :: xr,xi
 INTEGER(KIND=ik)                  :: it,i_field
 REAL(KIND=rk)                     :: tt
 REAL(KIND=rk),DIMENSION(3)        :: vr,vi
 CHARACTER(LEN=100)                :: address
 CHARACTER(LEN=6)                  :: iteration
!.   .    .   .    .   .    .   .    .   .    .   .    .   .    .   .

   WRITE(iteration,'(I6.6)')itmin
   address=trim(path)//trim(ver)//trim(file_bin)//trim(iteration)//'.bin'
!~    address=trim(path)//trim(in_file)
   print*,'reading : ',address
   OPEN(unit=2,file=address,status='unknown',form='unformatted')

   READ(2) tt
 L00 : DO xx=1,nx/2
   L10 : DO yy=ny/2+1+al*ny/2,ny+al*ny/2
      L20 : DO zz=nz/2+1+al*nz/2,nz+al*nz/2
            READ(2)(vr(ii),vi(ii),ii=1,3)
            DO ii=1,3
            uu_C(xx,yy,zz,ii)=cmplx(vr(ii),vi(ii))
            ENDDO
      ENDDO L20

      L30 : DO zz=1,nz/2
            READ(2)(vr(ii),vi(ii),ii=1,3)
            DO ii=1,3
            uu_C(xx,yy,zz,ii)=cmplx(vr(ii),vi(ii))
            ENDDO
      ENDDO L30
   ENDDO L10

   L11 : DO yy=1,ny/2
      L21 : DO zz=nz/2+1+al*nz/2,nz+al*nz/2
            READ(2)(vr(ii),vi(ii),ii=1,3)
            DO ii=1,3
            uu_C(xx,yy,zz,ii)=cmplx(vr(ii),vi(ii))
            ENDDO
      ENDDO L21

      L31 : DO zz=1,nz/2
            READ(2)(vr(ii),vi(ii),ii=1,3)
            DO ii=1,3
            uu_C(xx,yy,zz,ii)=cmplx(vr(ii),vi(ii))
            ENDDO
      ENDDO L31
   ENDDO L11
 ENDDO L00
   CLOSE(2)

END SUBROUTINE IO_read_field


!!!. . . . . . . . . . WRITE FIELD . . . . . . . . . . . . . . . . . . .

SUBROUTINE IO_write_velocity_field(it)
!! Writes in the output file the velocity field in the Fourier space
implicit none
 INTEGER(KIND=ik)                  :: ii,jj,kk,hh,xx,yy,zz
 INTEGER(KIND=ik)                  :: it
 REAL(KIND=rk),DIMENSION(3)        :: vr,vi
 CHARACTER(LEN=250)                :: header
 CHARACTER(LEN=100)                :: address
 CHARACTER(LEN=6)                  :: iteration

!!!.   .    .   .    .   .    .   .    .   .    .   .    .   .    .   .

IF (it-1 /= itmin) THEN
   WRITE(iteration,'(I6.6)')it-1
   address=trim(path)//trim(ver)//'_'//trim(file_bin)//'_'//trim(iteration)//'.dat'
   print*,'writing : ',address
   OPEN(unit=3,file=address,status='unknown',form='unformatted')

   WRITE(3) t

   DO xx=1,nx/2
      DO yy=ny/2+1,ny
         DO zz=nz/2+1,nz
            vr(:)=dreal(uu_C(xx,yy,zz,:))
            vi(:)=dimag(uu_C(xx,yy,zz,:))
             write(3) (vr(hh),vi(hh),hh=1,3)
         ENDDO
         DO zz=1,nz/2
            vr(:)=dreal(uu_C(xx,yy,zz,:))
            vi(:)=dimag(uu_C(xx,yy,zz,:))
             write(3) (vr(hh),vi(hh),hh=1,3)
         ENDDO
      ENDDO

      DO yy=1,ny/2
         DO zz=nz/2+1,nz
             vr(:)=dreal(uu_C(xx,yy,zz,:))
             vi(:)=dimag(uu_C(xx,yy,zz,:))
             write(3) (vr(hh),vi(hh),hh=1,3)
         ENDDO
         DO zz=1,nz/2
            vr(:)=dreal(uu_C(xx,yy,zz,:))
            vi(:)=dimag(uu_C(xx,yy,zz,:))
            write(3) (vr(hh),vi(hh),hh=1,3)
         ENDDO
       ENDDO

    ENDDO

   CLOSE(3)
ENDIF

END SUBROUTINE IO_write_velocity_field
!!!.....................................................................
!!!.....................................................................
SUBROUTINE IO_write_hdf5_velocity_field(it)

 IMPLICIT NONE
  INTEGER(kind=IK)                         :: it
  integer(hid_t)                           :: file
  integer(hsize_t)                         :: gdims(3), goffset(3)
  integer(hsize_t)                         :: ldims(3), loffset(3)
  integer(hsize_t)                         :: dimm(3), current_dims(3)
  integer                                  :: ndim, nrank, nsize, current_ndim

   CALL p3dfft_btran_c2r_many (uu_C,Csize(1)*Csize(2)*Csize(3),uu, &
               Rsize(1)*Rsize(2)*Rsize(3),3,'tff')


  call h5util_init()
	call h5util_create_file("h5test.h5")
  call h5util_open_file("h5test.h5", file)
  call h5util_put_attribute(file, "it", it)
  call h5util_put_attribute(file, "Dims", (/nxp,nyp,nzp/))
  call h5util_put_attribute(file, "Domain", (/xl,yl,zl/))

  ndim          = 3
  ldims         = shape(uu(:,:,:,2))
  loffset       = 0
  gdims         = ldims
  gdims(ndim)   = ldims(ndim) * n_proc
  goffset       = 0
  goffset(ndim) = ldims(ndim) * proc_id
  call h5util_create_dataset(file, "U", H5T_NATIVE_DOUBLE, ndim, gdims)
  call h5util_create_dataset(file, "V", H5T_NATIVE_DOUBLE, ndim, gdims)
  call h5util_create_dataset(file, "W", H5T_NATIVE_DOUBLE, ndim, gdims)
  call h5util_write_dataset(file, "U", ndim, ldims, ldims, loffset, goffset, &
       & reshape(uu(:,:,:,1), (/size(uu(:,:,:,1))/))  )
  call h5util_write_dataset(file, "V", ndim, ldims, ldims, loffset, goffset, &
       & reshape(uu(:,:,:,2), (/size(uu(:,:,:,2))/))  )
  call h5util_write_dataset(file, "W", ndim, ldims, ldims, loffset, goffset, &
       & reshape(uu(:,:,:,3), (/size(uu(:,:,:,3))/))  )
  call h5util_close_file(file)
  call h5util_finalize()



   CALL p3dfft_ftran_r2c_many (uu,Rsize(1)*Rsize(2)*Rsize(3),uu_C, &
               Csize(1)*Csize(2)*Csize(3),3,'fft')
               uu_C=uu_C/REAL(nxp*nyp*nzp,KIND=rk)
END SUBROUTINE IO_write_hdf5_velocity_field
!!!.....................................................................


END MODULE IO_mod

!!!.....................................................................
!!!.....................................................................
!!!.....................................................................
!!!.....................................................................
