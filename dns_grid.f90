!=======================================================================
!=======================================================================
PROGRAM dns_grid

USE parameters_mod         !fixed parameters and constants, kinds
USE variables_and_IO_mod      !program's main variables. Contains input_grid
USE fft_mod
USE IO_mod
USE time_advancement_mod
implicit none
!=======================================================================
REAL(KIND=4)          :: t1,t2
REAL(kind=rk)         :: x
LOGICAL               :: stats_time

!reads case variables and allocates the arrays
CALL read_input_parameters
CALL memory_initialization
CALL fft_initialization
CALL wave_numbers
!check newf_switch from sim_par.dat. 0=reads old field 1=generates ampty field
! if (newf_switch==0) CALL read_field
! if (newf_switch==1) CALL new_field
CALL read_field
CALL re_indexing
CALL rk_initialize
t=REAL(itmin,KIND=rk)*dt
CALL CPU_TIME(t1)
TIMELOOP: DO it=itmin+1,itmax
RK_LOOP:  DO n_k=1,rk_steps
          print *,'it ,t , ik =',it,t,n_k,rk_steps
          t=t+dt/REAL(rk_steps,KIND=rk)


          IF (MOD(it-1_ik,it_stat)==0_ik) stats_time=.true.

          CALL partial_right_hand_side



          ! IF (MOD(it-1_ik,rk_steps)==0_ik) CALL grid_forcing_update

          CALL nonlinear

          CALL linear

          IF (it==it_wrt) THEN
          CALL B_FFT(uu_C,uu)
          CALL write_field(it)
                    CALL F_FFT(uu,uu_C)
          ! CALL random_seed_io
          it_wrt=it_wrt+it_out
          ENDIF

ENDDO RK_LOOP

ENDDO TIMELOOP
CALL CPU_TIME(t2)
PRINT *,'Computation time:', t2- t1
!=======================================================================
END PROGRAM dns_grid
!=======================================================================
!!TODO test lettura velocita (fatto)
!!TODO test fftb (fatto)
!!TODO test coefficenti rk (fatto)
!!TODO test prnk qnrk (fatto)
!!TODO test k_quad (fatto) (invertire tutti i segni con k_quad)
!!TODO test p_uu prhs(non tutti i punti coincidono, differenza dell'ordine 10^-53)
