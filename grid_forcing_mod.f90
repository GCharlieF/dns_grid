!!!.....................................................................
!!!..........................GRID FORCING MOD...........................
!!!.....................................................................

MODULE grid_forcing_mod

USE parameters_mod
USE variables_mod
implicit none

CONTAINS

!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
SUBROUTINE grid_forcing_init
 !! Initializes at the first iteration the grid body force and updates it
 !! every time step.
 !! Bicubic interpolation in space and linear in time

 implicit none
 REAL(KIND=rk),DIMENSION(3)              :: ao,an,as,ak
 REAL(KIND=rk)                           :: u_rand,v_rand,w_rand
 REAL(KIND=rk)                           :: dy,dz
 REAL(KIND=rk)                           :: yy,zz !position in the grid
 INTEGER                                 :: y,z,yn,zn !y,z: mesh,yn,zn: grid
  character(len=100)                      :: address !for printing out test forcing distribution
 !grid steps
 dy=REAL(nc,KIND=rk)/REAL(nyp,KIND=rk)
 dz=REAL(nc,KIND=rk)/REAL(nzp,KIND=rk)


            call random_seed(put=[12,12])

 !forcing startup (sets first couple f_prev/f_next)
 IF (it==0) THEN
 ! ZL100: DO zn=1,nzp,nzp/nc
 ! YL100: DO yn=1,nyp,nyp/nc
 !        !fills the grid nodes (nodes with assigned forcing value)
 !        !with a random value between f_amp and - f_amp
 !        CALL random_number(u_rand)
 !        CALL random_number(v_rand)
 !        CALL random_number(w_rand)
 !        fu_prev(yn,zn)=f_amp*(u_rand*2_rk - 1_rk)
 !        fv_prev(yn,zn)=f_amp*(v_rand*2_rk - 1_rk)
 !        fw_prev(yn,zn)=f_amp*(w_rand*2_rk - 1_rk)
 !        ! CALL random_number(u_rand)
 !        ! CALL random_number(v_rand)
 !        ! CALL random_number(w_rand)
 !        fu_next(yn,zn)=f_amp*(u_rand*2_rk - 1_rk)
 !        fv_next(yn,zn)=f_amp*(v_rand*2_rk - 1_rk)
 !        fw_next(yn,zn)=f_amp*(w_rand*2_rk - 1_rk)
 !        ENDDO YL100
 !        ENDDO ZL100
 ! ENDIF

 ZL100: DO zn=1,nzp,nzp/nc
 YL100: DO yn=1,nyp,nyp/nc
        !fills the grid nodes (nodes with assigned forcing value)
        !with a random value between f_amp and - f_amp
        CALL random_number(u_rand)
        CALL random_number(v_rand)
        CALL random_number(w_rand)
        fu_prev(yn,zn)=f_amp*(u_rand*2_rk - 1_rk)
        fv_prev(yn,zn)=f_amp*(v_rand*2_rk - 1_rk)
        fw_prev(yn,zn)=f_amp*(w_rand*2_rk - 1_rk)
        ENDDO YL100
        ENDDO ZL100


  ZL101: DO zn=1,nzp,nzp/nc
  YL101: DO yn=1,nyp,nyp/nc
         !fills the grid nodes (nodes with assigned forcing value)
         !with a random value between f_amp and - f_amp
         CALL random_number(u_rand)
         CALL random_number(v_rand)
         CALL random_number(w_rand)
         fu_next(yn,zn)=f_amp*(u_rand*2_rk - 1_rk)
         fv_next(yn,zn)=f_amp*(v_rand*2_rk - 1_rk)
         fw_next(yn,zn)=f_amp*(w_rand*2_rk - 1_rk)
       ENDDO YL101
     ENDDO ZL101
  ENDIF

        !assigns the value in the remaning points of the plane y-z
        !via bicubic interpolation
 ZN200: DO zn=1,nzp,nzp/nc
 YN200: DO yn=1,nyp,nyp/nc

 ZL201: DO z=zn,zn+nzp/nc-1
 YL201: DO y=yn,yn+nyp/nc-1

        yy=REAL(y-yn,KIND=rk)*dy
        zz=REAL(z-zn,KIND=rk)*dz

        IF ((z .LE. nzp-nzp/nc).AND.(y .LE. nyp-nyp/nc)) THEN
        ao(1)=fu_prev(yn,zn)
        an(1)=fu_prev(yn,zn)-fu_prev(yn,zn+nzp/nc)
        as(1)=fu_prev(yn,zn)-fu_prev(yn+nyp/nc,zn)
        ak(1)=fu_prev(yn,zn)-fu_prev(yn,zn+nzp/nc)-fu_prev(yn+nyp/nc,zn)&
             +fu_prev(yn+nyp/nc,zn+nzp/nc)

        ao(2)=fv_prev(yn,zn)
        an(2)=fv_prev(yn,zn)-fv_prev(yn,zn+nzp/nc)
        as(2)=fv_prev(yn,zn)-fv_prev(yn+nyp/nc,zn)
        ak(2)=fv_prev(yn,zn)-fv_prev(yn,zn+nzp/nc)-fv_prev(yn+nyp/nc,zn)&
             +fv_prev(yn+nyp/nc,zn+nzp/nc)

        ao(3)=fw_prev(yn,zn)
        an(3)=fw_prev(yn,zn)-fw_prev(yn,zn+nzp/nc)
        as(3)=fw_prev(yn,zn)-fw_prev(yn+nyp/nc,zn)
        ak(3)=fw_prev(yn,zn)-fw_prev(yn,zn+nzp/nc)-fw_prev(yn+nyp/nc,zn)&
             +fw_prev(yn+nyp/nc,zn+nzp/nc)

        ENDIF

        IF((z .GT. nzp-nzp/nc) .AND. (y.LE.nyp-nyp/nc)) THEN
        ao(1)=fu_prev(yn,zn)
        an(1)=fu_prev(yn,zn)-fu_prev(yn,1)
        as(1)=fu_prev(yn,zn)-fu_prev(yn+nyp/nc,zn)
        ak(1)=fu_prev(yn,zn)-fu_prev(yn,1)-fu_prev(yn+nyp/nc,zn)&
             +fu_prev(yn+nyp/nc,1)

        ao(2)=fv_prev(yn,zn)
        an(2)=fv_prev(yn,zn)-fv_prev(yn,1)
        as(2)=fv_prev(yn,zn)-fv_prev(yn+nyp/nc,zn)
        ak(2)=fv_prev(yn,zn)-fv_prev(yn,1)-fv_prev(yn+nyp/nc,zn)&
             +fv_prev(yn+nyp/nc,1)

        ao(3)=fw_prev(yn,zn)
        an(3)=fw_prev(yn,zn)-fw_prev(yn,1)
        as(3)=fw_prev(yn,zn)-fw_prev(yn+nyp/nc,zn)
        ak(3)=fw_prev(yn,zn)-fw_prev(yn,1)-fw_prev(yn+nyp/nc,zn)&
             +fw_prev(yn+nyp/nc,1)
        ENDIF

        IF((z .LE. nzp-nzp/nc) .AND. (y.GT.nyp-nyp/nc))then

        ao(1)=fu_prev(yn,zn)
        an(1)=fu_prev(yn,zn)-fu_prev(yn,zn+nzp/nc)
        as(1)=fu_prev(yn,zn)-fu_prev(1,zn)
        ak(1)=fu_prev(yn,zn)-fu_prev(yn,zn+nzp/nc)-fu_prev(1,zn)&
             +fu_prev(1,zn+nzp/nc)

        ao(2)=fv_prev(yn,zn)
        an(2)=fv_prev(yn,zn)-fv_prev(yn,zn+nzp/nc)
        as(2)=fv_prev(yn,zn)-fv_prev(1,zn)
        ak(2)=fv_prev(yn,zn)-fv_prev(yn,zn+nzp/nc)-fv_prev(1,zn)&
             +fv_prev(1,zn+nzp/nc)

        ao(3)=fw_prev(yn,zn)
        an(3)=fw_prev(yn,zn)-fw_prev(yn,zn+nzp/nc)
        as(3)=fw_prev(yn,zn)-fw_prev(1,zn)
        ak(3)=fw_prev(yn,zn)-fw_prev(yn,zn+nzp/nc)-fw_prev(1,zn)&
             +fw_prev(1,zn+nzp/nc)

        ENDIF
        IF((z .GT. nzp-nzp/nc) .AND. (y.GT.nyp-nyp/nc))then

        ao(1)=fu_prev(yn,zn)
        an(1)=fu_prev(yn,zn)-fu_prev(yn,1)
        as(1)=fu_prev(yn,zn)-fu_prev(1,zn)
        ak(1)=fu_prev(yn,zn)-fu_prev(yn,1)-fu_prev(1,zn)&
             +fu_prev(1,1)

        ao(2)=fv_prev(yn,zn)
        an(2)=fv_prev(yn,zn)-fv_prev(yn,1)
        as(2)=fv_prev(yn,zn)-fv_prev(1,zn)
        ak(2)=fv_prev(yn,zn)-fv_prev(yn,1)-fv_prev(1,zn)&
             +fv_prev(1,1)

        ao(3)=fw_prev(yn,zn)
        an(3)=fw_prev(yn,zn)-fw_prev(yn,1)
        as(3)=fw_prev(yn,zn)-fw_prev(1,zn)
        ak(3)=fw_prev(yn,zn)-fw_prev(yn,1)-fw_prev(1,zn)&
             +fw_prev(1,1)

        ENDIF
          fu_prev(y,z)=ao(1)-3_rk*as(1)*yy**2-3_rk*an(1)*zz**2&
                +2_rk*as(1)*yy**3+2_rk*an(1)*zz**3-6_rk*ak(1)*(yy**3)*(zz**2)&
                -6_rk*ak(1)*(yy**2)*(zz**3)+9_rk*ak(1)*(yy**2)*(zz**2)&
                +4_rk*ak(1)*(yy**3)*(zz**3)

          fv_prev(y,z)=ao(2)-3_rk*as(2)*yy**2-3_rk*an(2)*zz**2&
                +2_rk*as(2)*yy**3+2_rk*an(2)*zz**3-6_rk*ak(2)*(yy**3)*(zz**2)&
                -6_rk*ak(2)*(yy**2)*(zz**3)+9_rk*ak(2)*(yy**2)*(zz**2)&
                +4_rk*ak(2)*(yy**3)*(zz**3)

          fw_prev(y,z)=ao(3)-3_rk*as(3)*yy**2-3_rk*an(3)*zz**2&
                +2_rk*as(3)*yy**3+2_rk*an(3)*zz**3-6_rk*ak(3)*(yy**3)*(zz**2)&
                -6_rk*ak(3)*(yy**2)*(zz**3)+9_rk*ak(3)*(yy**2)*(zz**2)&
                +4_rk*ak(3)*(yy**3)*(zz**3)
        ENDDO YL201
        ENDDO ZL201

        ENDDO YN200
        ENDDO ZN200

ZN300: DO zn=1,nzp,nzp/nc
YN300: DO yn=1,nyp,nyp/nc

 ZL301: DO z=zn,zn+nzp/nc-1
 YL301: DO y=yn,yn+nyp/nc-1

        yy=REAL(y-yn,KIND=rk)*dy
        zz=REAL(z-zn,KIND=rk)*dz

        IF ((z .LE. nzp-nzp/nc).AND.(y .LE. nyp-nyp/nc)) THEN
        ao(1)=fu_next(yn,zn)
        an(1)=fu_next(yn,zn)-fu_next(yn,zn+nzp/nc)
        as(1)=fu_next(yn,zn)-fu_next(yn+nyp/nc,zn)
        ak(1)=fu_next(yn,zn)-fu_next(yn,zn+nzp/nc)-fu_next(yn+nyp/nc,zn)&
             +fu_next(yn+nyp/nc,zn+nzp/nc)

        ao(2)=fv_next(yn,zn)
        an(2)=fv_next(yn,zn)-fv_next(yn,zn+nzp/nc)
        as(2)=fv_next(yn,zn)-fv_next(yn+nyp/nc,zn)
        ak(2)=fv_next(yn,zn)-fv_next(yn,zn+nzp/nc)-fv_next(yn+nyp/nc,zn)&
             +fv_next(yn+nyp/nc,zn+nzp/nc)

        ao(3)=fw_next(yn,zn)
        an(3)=fw_next(yn,zn)-fw_next(yn,zn+nzp/nc)
        as(3)=fw_next(yn,zn)-fw_next(yn+nyp/nc,zn)
        ak(3)=fw_next(yn,zn)-fw_next(yn,zn+nzp/nc)-fw_next(yn+nyp/nc,zn)&
             +fw_next(yn+nyp/nc,zn+nzp/nc)

        ENDIF

        IF((z .GT. nzp-nzp/nc) .AND. (y.LE.nyp-nyp/nc)) THEN
        ao(1)=fu_next(yn,zn)
        an(1)=fu_next(yn,zn)-fu_next(yn,1)
        as(1)=fu_next(yn,zn)-fu_next(yn+nyp/nc,zn)
        ak(1)=fu_next(yn,zn)-fu_next(yn,1)-fu_next(yn+nyp/nc,zn)&
             +fu_next(yn+nyp/nc,1)

        ao(2)=fv_next(yn,zn)
        an(2)=fv_next(yn,zn)-fv_next(yn,1)
        as(2)=fv_next(yn,zn)-fv_next(yn+nyp/nc,zn)
        ak(2)=fv_next(yn,zn)-fv_next(yn,1)-fv_next(yn+nyp/nc,zn)&
             +fv_next(yn+nyp/nc,1)

        ao(3)=fw_next(yn,zn)
        an(3)=fw_next(yn,zn)-fw_next(yn,1)
        as(3)=fw_next(yn,zn)-fw_next(yn+nyp/nc,zn)
        ak(3)=fw_next(yn,zn)-fw_next(yn,1)-fw_next(yn+nyp/nc,zn)&
             +fw_next(yn+nyp/nc,1)
        ENDIF

        IF((z .LE. nzp-nzp/nc) .AND. (y.GT.nyp-nyp/nc))then

        ao(1)=fu_next(yn,zn)
        an(1)=fu_next(yn,zn)-fu_next(yn,zn+nzp/nc)
        as(1)=fu_next(yn,zn)-fu_next(1,zn)
        ak(1)=fu_next(yn,zn)-fu_next(yn,zn+nzp/nc)-fu_next(1,zn)&
             +fu_next(1,zn+nzp/nc)

        ao(2)=fv_next(yn,zn)
        an(2)=fv_next(yn,zn)-fv_next(yn,zn+nzp/nc)
        as(2)=fv_next(yn,zn)-fv_next(1,zn)
        ak(2)=fv_next(yn,zn)-fv_next(yn,zn+nzp/nc)-fv_next(1,zn)&
             +fv_next(1,zn+nzp/nc)

        ao(3)=fw_next(yn,zn)
        an(3)=fw_next(yn,zn)-fw_next(yn,zn+nzp/nc)
        as(3)=fw_next(yn,zn)-fw_next(1,zn)
        ak(3)=fw_next(yn,zn)-fw_next(yn,zn+nzp/nc)-fw_next(1,zn)&
             +fw_next(1,zn+nzp/nc)

        ENDIF
        IF((z .GT. nzp-nzp/nc) .AND. (y.GT.nyp-nyp/nc)) THEN

        ao(1)=fu_next(yn,zn)
        an(1)=fu_next(yn,zn)-fu_next(yn,1)
        as(1)=fu_next(yn,zn)-fu_next(1,zn)
        ak(1)=fu_next(yn,zn)-fu_next(yn,1)-fu_next(1,zn)&
             +fu_next(1,1)

        ao(2)=fv_next(yn,zn)
        an(2)=fv_next(yn,zn)-fv_next(yn,1)
        as(2)=fv_next(yn,zn)-fv_next(1,zn)
        ak(2)=fv_next(yn,zn)-fv_next(yn,1)-fv_next(1,zn)&
             +fv_next(1,1)

        ao(3)=fw_next(yn,zn)
        an(3)=fw_next(yn,zn)-fw_next(yn,1)
        as(3)=fw_next(yn,zn)-fw_next(1,zn)
        ak(3)=fw_next(yn,zn)-fw_next(yn,1)-fw_next(1,zn)&
             +fw_next(1,1)

        ENDIF
          fu_next(y,z)=ao(1)-3_rk*as(1)*yy**2-3_rk*an(1)*zz**2&
                +2_rk*as(1)*yy**3+2_rk*an(1)*zz**3-6_rk*ak(1)*(yy**3)*(zz**2)&
                -6_rk*ak(1)*(yy**2)*(zz**3)+9_rk*ak(1)*(yy**2)*(zz**2)&
                +4_rk*ak(1)*(yy**3)*(zz**3)

          fv_next(y,z)=ao(2)-3_rk*as(2)*yy**2-3_rk*an(2)*zz**2&
                +2_rk*as(2)*yy**3+2_rk*an(2)*zz**3-6_rk*ak(2)*(yy**3)*(zz**2)&
                -6_rk*ak(2)*(yy**2)*(zz**3)+9_rk*ak(2)*(yy**2)*(zz**2)&
                +4_rk*ak(2)*(yy**3)*(zz**3)

          fw_next(y,z)=ao(3)-3_rk*as(3)*yy**2-3_rk*an(3)*zz**2&
                +2_rk*as(3)*yy**3+2_rk*an(3)*zz**3-6_rk*ak(3)*(yy**3)*(zz**2)&
                -6_rk*ak(3)*(yy**2)*(zz**3)+9_rk*ak(3)*(yy**2)*(zz**2)&
                +4_rk*ak(3)*(yy**3)*(zz**3)
        ENDDO YL301
        ENDDO ZL301
        ENDDO YN300
        ENDDO ZN300


 !       address = 'Variables = "y","z","f1","f2","f3","f1n","f2n","f3n"'
 !       write(777,*) address
 !       address = 'ZONE I=1234 J=1234'
 !       write(address(08:11),1001) nyp
 !       write(address(15:18),1001) nzp
 ! 1001  format(i4.4)
 !       write(777,*) address
 !       do 4000 z=1,nzp
 !       do 4000 y=1,nyp
 !       write(777,*)float(y),float(z),fu_prev(y,z),fv_prev(y,z),fw_prev(y,z),&
 !            fu_next(y,z),fv_next(y,z),fw_next(y,z)
 ! 4000  continue
 !       stop

END SUBROUTINE grid_forcing_init
!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
SUBROUTINE grid_forcing_update

implicit none

 !! Initialize at the first iteration the grid body force and updates it
 !! every time step.
 !! Bicubic interpolation in space and linear in time

 REAL(KIND=rk),DIMENSION(3)              :: ao,an,as,ak
 REAL(KIND=rk)                           :: u_rand,v_rand,w_rand
 REAL(KIND=rk)                           :: dy,dz
 REAL(KIND=rk)                           :: yy,zz !position in the grid
 INTEGER                                 :: y,z,yn,zn !y,z: mesh,yn,zn: grid
  character(len=100)                      :: address !for printing out test forcing distribution
 !grid steps
 dy=REAL(nc,KIND=rk)/REAL(nyp,KIND=rk)
 dz=REAL(nc,KIND=rk)/REAL(nzp,KIND=rk)

 !After dt_forc iterations updates f_next with new values and put the old ones in f_prev
      print *,'adv'
 IF ((it/=1).AND.(MOD(it-1_ik,dt_forc)==0_ik)) THEN
  fu_prev=fu_next
  fv_prev=fv_next
  fw_prev=fw_next

  !fills the grid nodes (nodes with assigned forcing value)
  !with a random value between f_amp and - f_amp
  ZL100: DO zn=1,nzp,nzp/nc
  YL100: DO yn=1,nyp,nyp/nc
          CALL random_number(u_rand)
          CALL random_number(v_rand)
          CALL random_number(w_rand)
          fu_next(yn,zn)=f_amp*(u_rand*2_rk - 1_rk)
          fv_next(yn,zn)=f_amp*(v_rand*2_rk - 1_rk)
          fw_next(yn,zn)=f_amp*(w_rand*2_rk - 1_rk)
         ENDDO YL100
         ENDDO ZL100
 ENDIF


 ZN200: DO zn=1,nzp,nzp/nc
 YN200: DO yn=1,nyp,nyp/nc

 ZL201: DO z=zn,zn+nzp/nc-1
 YL201: DO y=yn,yn+nyp/nc-1

        yy=REAL(y-yn,KIND=rk)*dy
        zz=REAL(z-zn,KIND=rk)*dz

        IF ((z .LE. nzp-nzp/nc).AND.(y .LE. nyp-nyp/nc)) THEN
        ao(1)=fu_next(yn,zn)
        an(1)=fu_next(yn,zn)-fu_next(yn,zn+nzp/nc)
        as(1)=fu_next(yn,zn)-fu_next(yn+nyp/nc,zn)
        ak(1)=fu_next(yn,zn)-fu_next(yn,zn+nzp/nc)-fu_next(yn+nyp/nc,zn)&
             +fu_next(yn+nyp/nc,zn+nzp/nc)

        ao(2)=fv_next(yn,zn)
        an(2)=fv_next(yn,zn)-fv_next(yn,zn+nzp/nc)
        as(2)=fv_next(yn,zn)-fv_next(yn+nyp/nc,zn)
        ak(2)=fv_next(yn,zn)-fv_next(yn,zn+nzp/nc)-fv_next(yn+nyp/nc,zn)&
             +fv_next(yn+nyp/nc,zn+nzp/nc)

        ao(3)=fw_next(yn,zn)
        an(3)=fw_next(yn,zn)-fw_next(yn,zn+nzp/nc)
        as(3)=fw_next(yn,zn)-fw_next(yn+nyp/nc,zn)
        ak(3)=fw_next(yn,zn)-fw_next(yn,zn+nzp/nc)-fw_next(yn+nyp/nc,zn)&
             +fw_next(yn+nyp/nc,zn+nzp/nc)

        ENDIF

        IF((z .GT. nzp-nzp/nc) .AND. (y.LE.nyp-nyp/nc)) THEN
        ao(1)=fu_next(yn,zn)
        an(1)=fu_next(yn,zn)-fu_next(yn,1)
        as(1)=fu_next(yn,zn)-fu_next(yn+nyp/nc,zn)
        ak(1)=fu_next(yn,zn)-fu_next(yn,1)-fu_next(yn+nyp/nc,zn)&
             +fu_next(yn+nyp/nc,1)

        ao(2)=fv_next(yn,zn)
        an(2)=fv_next(yn,zn)-fv_next(yn,1)
        as(2)=fv_next(yn,zn)-fv_next(yn+nyp/nc,zn)
        ak(2)=fv_next(yn,zn)-fv_next(yn,1)-fv_next(yn+nyp/nc,zn)&
             +fv_next(yn+nyp/nc,1)

        ao(3)=fw_next(yn,zn)
        an(3)=fw_next(yn,zn)-fw_next(yn,1)
        as(3)=fw_next(yn,zn)-fw_next(yn+nyp/nc,zn)
        ak(3)=fw_next(yn,zn)-fw_next(yn,1)-fw_next(yn+nyp/nc,zn)&
             +fw_next(yn+nyp/nc,1)
        ENDIF

        IF((z .LE. nzp-nzp/nc) .AND. (y.GT.nyp-nyp/nc))then

        ao(1)=fu_next(yn,zn)
        an(1)=fu_next(yn,zn)-fu_next(yn,zn+nzp/nc)
        as(1)=fu_next(yn,zn)-fu_next(1,zn)
        ak(1)=fu_next(yn,zn)-fu_next(yn,zn+nzp/nc)-fu_next(1,zn)&
             +fu_next(1,zn+nzp/nc)

        ao(2)=fv_next(yn,zn)
        an(2)=fv_next(yn,zn)-fv_next(yn,zn+nzp/nc)
        as(2)=fv_next(yn,zn)-fv_next(1,zn)
        ak(2)=fv_next(yn,zn)-fv_next(yn,zn+nzp/nc)-fv_next(1,zn)&
             +fv_next(1,zn+nzp/nc)

        ao(3)=fw_next(yn,zn)
        an(3)=fw_next(yn,zn)-fw_next(yn,zn+nzp/nc)
        as(3)=fw_next(yn,zn)-fw_next(1,zn)
        ak(3)=fw_next(yn,zn)-fw_next(yn,zn+nzp/nc)-fw_next(1,zn)&
             +fw_next(1,zn+nzp/nc)
       ENDIF
       IF((z .GT. nzp-nzp/nc) .AND. (y.GT.nyp-nyp/nc)) THEN

             ao(1)=fu_next(yn,zn)
             an(1)=fu_next(yn,zn)-fu_next(yn,1)
             as(1)=fu_next(yn,zn)-fu_next(1,zn)
             ak(1)=fu_next(yn,zn)-fu_next(yn,1)-fu_next(1,zn)&
                  +fu_next(1,1)

             ao(2)=fv_next(yn,zn)
             an(2)=fv_next(yn,zn)-fv_next(yn,1)
             as(2)=fv_next(yn,zn)-fv_next(1,zn)
             ak(2)=fv_next(yn,zn)-fv_next(yn,1)-fv_next(1,zn)&
                  +fv_next(1,1)

             ao(3)=fw_next(yn,zn)
             an(3)=fw_next(yn,zn)-fw_next(yn,1)
             as(3)=fw_next(yn,zn)-fw_next(1,zn)
             ak(3)=fw_next(yn,zn)-fw_next(yn,1)-fw_next(1,zn)&
                  +fw_next(1,1)

        ENDIF


          fu_next(y,z)=ao(1)-3_rk*as(1)*yy**2-3_rk*an(1)*zz**2&
                +2_rk*as(1)*yy**3+2_rk*an(1)*zz**3-6_rk*ak(1)*(yy**3)*(zz**2)&
                -6_rk*ak(1)*(yy**2)*(zz**3)+9_rk*ak(1)*(yy**2)*(zz**2)&
                +4_rk*ak(1)*(yy**3)*(zz**3)

          fv_next(y,z)=ao(2)-3_rk*as(2)*yy**2-3_rk*an(2)*zz**2&
                +2_rk*as(2)*yy**3+2_rk*an(2)*zz**3-6_rk*ak(2)*(yy**3)*(zz**2)&
                -6_rk*ak(2)*(yy**2)*(zz**3)+9_rk*ak(2)*(yy**2)*(zz**2)&
                +4_rk*ak(2)*(yy**3)*(zz**3)

          fw_next(y,z)=ao(3)-3_rk*as(3)*yy**2-3_rk*an(3)*zz**2&
                +2_rk*as(3)*yy**3+2_rk*an(3)*zz**3-6_rk*ak(3)*(yy**3)*(zz**2)&
                -6_rk*ak(3)*(yy**2)*(zz**3)+9_rk*ak(3)*(yy**2)*(zz**2)&
                +4_rk*ak(3)*(yy**3)*(zz**3)
        ENDDO YL201
        ENDDO ZL201
        ENDDO YN200
        ENDDO ZN200

 !Time interpolation between f_prev and f_next
 t_weight(1)=REAL(dt_forc - MOD(it-1_ik,dt_forc) ,KIND=rk)/REAL(dt_forc,KIND=rk)
 t_weight(2)=REAL(MOD(it-1_ik,dt_forc),KIND=rk)/REAL(dt_forc, KIND=rk)
  print *,'weights',t_weight(1),t_weight(2)
  fu=t_weight(1)*fu_prev+t_weight(2)*fu_next
  fv=t_weight(1)*fv_prev+t_weight(2)*fv_next
  fw=t_weight(1)*fw_prev+t_weight(2)*fw_next

 ! if (it==31) then
 !       address = 'Variables = "y","z","f1","f2","f3"'
 !       write(777,*) address
 !       address = 'ZONE I=1234 J=1234'
 !       write(address(08:11),1001) nyp
 !       write(address(15:18),1001) nzp
 ! 1001  format(i4.4)
 !       write(777,*) address
 !       do 4000 z=1,nzp
 !       do 4000 y=1,nyp
 !       write(777,*)float(y),float(z),fu(y,z),fv(y,z),fw(y,z)
 ! 4000  continue
 !       stop
 ! endif


END SUBROUTINE grid_forcing_update
!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

END MODULE grid_forcing_mod
!.......................................................................
!.......................................................................

!TODO test distribuzione nel tempo
!TODO inizializzione forcing e restart
