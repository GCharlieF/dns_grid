!!!.....................................................................
!!!...............................IO MOD................................
!!!.....................................................................
MODULE IO_mod
 USE parameters_mod
 USE variables_and_IO_mod
 USE fft_mod
 implicit none
 INTEGER(KIND=ik),DIMENSION(:),ALLOCATABLE    	  ::ind_x,ind_y,ind_z
REAL(KIND=rk),DIMENSION(:),ALLOCATABLE 		  ::xc,yc,zc
contains


!!!. . . . . . . . . . RE-INDEXING . . . . . . . . . . . . . . . . . . .

SUBROUTINE re_indexing

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
!!!reads the input file


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
             if(xx.eq.1) then
              xr=1
             xi=2
             else
             xr=2*xx-1
             xi=2*xx
             endif
   L10 : DO yy=ny/2+1,ny
      L20 : DO zz=nz/2+1,nz
            READ(2)(vr(ii),vi(ii),ii=1,3)
            DO ii=1,3
            uu(xr,yy,zz,ii)=vr(ii)
            uu(xi,yy,zz,ii)=vi(ii)
            ENDDO
      ENDDO L20

      L30 : DO zz=1,nz/2
            READ(2)(vr(ii),vi(ii),ii=1,3)
             DO ii=1,3
            uu(xr,yy,zz,ii)=vr(ii)
            uu(xi,yy,zz,ii)=vi(ii)
            ENDDO
      ENDDO L30
   ENDDO L10

   L11 : DO yy=1,ny/2
      L21 : DO zz=nz/2+1,nz
            READ(2)(vr(ii),vi(ii),ii=1,3)
             DO ii=1,3
            uu(xr,yy,zz,ii)=vr(ii)
            uu(xi,yy,zz,ii)=vi(ii)
            ENDDO
      ENDDO L21

      L31 : DO zz=1,nz/2
            READ(2)(vr(ii),vi(ii),ii=1,3)
             DO ii=1,3
            uu(xr,yy,zz,ii)=vr(ii)
            uu(xi,yy,zz,ii)=vi(ii)
            ENDDO
      ENDDO L31
   ENDDO L11
 ENDDO L00
   CLOSE(2)

    DO zz=1,nz
    DO yy=1,ny
    DO xx=1,nx/2+1
    DO ii=1,3
     if(xx.eq.1) then
     xr=1
     xi=2
     ELSE
     xr=2*xx-1
     xi=2*xx
     ENDIF
    uu_C(xx,yy,zz,ii)=cmplx(uu(xr,yy,zz,ii),uu(xi,yy,zz,ii))
    ENDDO
    ENDDO
    ENDDO
    ENDDO
END SUBROUTINE read_field


!!!. . . . . . . . . . WRITE FIELD . . . . . . . . . . . . . . . . . . .

SUBROUTINE write_field(it)
!!!reads the input file
implicit none
 INTEGER(KIND=ik)                  :: ii,jj,kk,hh,xx,yy,zz
 INTEGER(KIND=ik)                  :: it
 REAL(KIND=rk)                     :: tp
 CHARACTER(LEN=250)                :: header
 CHARACTER(LEN=100)                :: address
 CHARACTER(LEN=6)                  :: iteration

!!!.   .    .   .    .   .    .   .    .   .    .   .    .   .    .   .

    WRITE(iteration,'(I6.6)')it
   address=trim(path)//trim(ver)//'_'//trim(file_bin)//'_'//trim(iteration)//'.dat'
   print*,'writing : ',address
   OPEN(unit=3,file=address,status='unknown',form='formatted')
   header = ' Variables = "x","y","z","u","v","w"'
   WRITE(3,*) header
   header = 'ZONE I=1234 J=1234 K=1234'
   WRITE(header(08:11),1000)  nxp
   WRITE(header(15:18),1000)  nyp
   WRITE(header(22:25),1000)  nzp
1000  format(i4.4)

   WRITE(3,*) header

   do kk=1,nzp
      zz=ind_z(kk)
         do jj=1,nyp
          yy=ind_y(jj)
            do ii=1,nxp
            xx=ind_x(ii)
             write(3,99)xc(ii),yc(jj),zc(kk),&
!~              (rr(xx,yy,zz,hh),hh=1,3)
             (uu(ii,jj,kk,hh),hh=1,3)
            enddo
         enddo
      enddo
   CLOSE(3)
99    format(6(e15.6,1x))
END SUBROUTINE write_field
!!!.....................................................................

END MODULE IO_mod

!!!.....................................................................
!!!.....................................................................
!!!.....................................................................
!!!.....................................................................
