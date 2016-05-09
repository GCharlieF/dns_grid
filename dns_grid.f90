!=======================================================================
!=======================================================================
PROGRAM dns_grid

USE parameters_mod         !fixed parameters and constants, kinds
USE variables_mod      !program's main variables. Contains input_grid
USE fft_mod
USE grid_forcing_mod
USE stats_and_probes_mod
USE hit_forcings_mod
USE IO_mod
USE time_advancement_mod
implicit none
!!=======================================================================
REAL(KIND=4)          :: t1,t2
REAL(kind=rk)         :: x
INTEGER :: mm

!! Reads case variables and allocates the arrays
CALL read_input_parameters !!var_m
CALL memory_initialization !!var_m
CALL fft_initialization    !!fft_m
CALL wave_numbers          !!var_m
CALL read_field            !!IO_m
 print *,'read field',uu_C(7,12,12,1)
CALL re_indexing
CALL rk_initialize         !!time_m
CALL dealiased_indeces     !!var_m
! CALL grid_forcing_init     !!grid_m
CALL alvelius_forcing_init
! CALL linear_forcing_init
! CALL divfree(uu_C)
t=REAL(itmin,KIND=rk)*dt
it=0

CALL CPU_TIME(t1)
TIMELOOP: DO it=itmin+1,itmax+1
RK_LOOP:  DO n_k=1,rk_steps
          PRINT *,'_________________________________________________________'
          PRINT *,'it , ik =',it,n_k
          PRINT *,'          '
          PRINT *,'dt, t >>>>>>>>>>',dt,t
          PRINT *,'          '

          t=t+dt/REAL(rk_steps,KIND=rk)

          !! Switch off and on the computation of statistics
          stats_time=.false.
          IF (MOD(it-1_ik,it_stat)==0 .and. n_k==1) stats_time=.true.

          CALL partial_right_hand_side  !!time_m

          !! Updates forcing distirbtution every n_k steps
          ! IF (n_k==1_ik) CALL grid_forcing_update !!grid_m

          CALL nonlinear                          !!time_m

          CALL linear                             !!time_m

          !! Write the velocity field every it_out time steps
          IF (MOD(it-1,it_out)==0 .AND. n_k==1) CALL write_velocity_field(it)

ENDDO RK_LOOP
          ! dt=dt_new
ENDDO TIMELOOP
CALL CPU_TIME(t2)
PRINT *,'Computation time:', t2- t1
!=======================================================================
END PROGRAM dns_grid
!=======================================================================
!!TODO aggiustare quando stamapre il campo
