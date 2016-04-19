!!!.....................................................................
!!!...............................IO MOD................................
!!!.....................................................................
MODULE IO_mod
 USE parameters_mod
 USE variables_mod
 USE fft_mod
 implicit none
 INTEGER(KIND=ik),DIMENSION(:),ALLOCATABLE    	  ::ind_x,ind_y,ind_z
REAL(KIND=rk),DIMENSION(:),ALLOCATABLE 		  ::xc,yc,zc
contains


!!!. . . . . . . . . . RE-INDEXING . . . . . . . . . . . . . . . . . . .

SUBROUTINE re_indexing
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
!~      DO ii=1,nxp/2
!~          ind_x(ii)=nxp/2+ii
!~          xx=ind_x(ii)
!~          xc(ii)=xl*REAL(xx-nxp-1,KIND=rk)/REAL(nxp,KIND=rk)
!~
!~ 	 ENDDO
!~      DO ii=nxp/2+1,nxp
!~          ind_x(ii)=ii-nxp/2
!~          xx=ind_x(ii)
!~          xc(ii)=xl*REAL(xx-1,KIND=rk)/REAL(nxp,KIND=rk)
!~ 	 ENDDO
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
END SUBROUTINE re_indexing
!!!. . . . . . . . . . . READ FIELD  . . . . . . . . . . . . . . . . . .

SUBROUTINE read_field
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

END SUBROUTINE read_field


!!!. . . . . . . . . . WRITE FIELD . . . . . . . . . . . . . . . . . . .

SUBROUTINE write_velocity_field(it)
!! Writes in the output file the velocity field in the Fourier space
implicit none
 INTEGER(KIND=ik)                  :: ii,jj,kk,hh,xx,yy,zz
 INTEGER(KIND=ik)                  :: it
 REAL(KIND=rk),DIMENSION(3)        :: vr,vi
 CHARACTER(LEN=250)                :: header
 CHARACTER(LEN=100)                :: address
 CHARACTER(LEN=6)                  :: iteration

!!!.   .    .   .    .   .    .   .    .   .    .   .    .   .    .   .

    WRITE(iteration,'(I6.6)')it
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

END SUBROUTINE write_velocity_field
!!!.....................................................................

END MODULE IO_mod

!!!.....................................................................
!!!.....................................................................
!!!.....................................................................
!!!.....................................................................
