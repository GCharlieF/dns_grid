!=======================================================================
!! TODO Clean,comment and explain
!! TODO Parallel grid_forcing (check random number generation)
!! TODO MPI communication for statistics
!! TODO Parrallel IO
!! TODO De-aliasing
!! TODO remove fft_mod
!! TODO way to avoid first r2c transform in MPI_initialize
!! TODO change x_loc,y_loc,z_loc in something more readable (preprocessor?)
!! TODO PPRINT statement for text and numbers
!! TODO Alvelious?
!=======================================================================
PROGRAM dns_grid
USE parameters_mod                                !! PAR
USE variables_mod                                 !! VAR
USE MPI_mod                                       !! MPI
USE fft_mod                                       !! FFT
USE grid_forcing_mod                              !! GRID
USE stats_and_probes_mod                          !! STATS
! USE hit_forcings_mod                              !! HIT
USE IO_mod                                        !! IO
USE time_advancement_mod                          !! TA
IMPLICIT NONE
#INCLUDE 'fpp_macros.h'
REAL(KIND=4)          :: t1,t2
REAL(kind=rk)         :: x
INTEGER :: mm
!!=======================================================================

!! Reads case variables and allocates the arrays
CALL VAR_read_input_parameters
CALL MPI_initialize
CALL VAR_memory_initialization
CALL VAR_wave_numbers
! CALL IO_read_field
CALL TA_rk_initialize
CALL VAR_dealiased_indeces   !! FIXME To be removed
CALL GRID_forcing_init

t=REAL(itmin,KIND=rk)*dt
it=0
CALL_BARRIER
MASTER CALL CPU_TIME(t1)
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

          CALL TA_partial_right_hand_side

          !! Updates forcing distirbtution every n_k steps
          IF (n_k==1_ik) CALL GRID_forcing_update !!grid_m

          CALL TA_nonlinear

          CALL TA_linear

          !! Write the velocity field every it_out time steps
          ! IF (MOD(it-1,it_out)==0 .AND. n_k==1) CALL IO_write_velocity_field(it)

ENDDO RK_LOOP
          ! dt=dt_new
ENDDO TIMELOOP
MASTER CALL CPU_TIME(t2)
PRINT *,'Computation time:', t2- t1
!=======================================================================
END PROGRAM dns_grid
!=======================================================================
