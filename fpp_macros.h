!!Macros for fortran pre-processor

!! Define macro for screen output
#DEFINE PPRINT(input) IF (proc_id==0) PRINT *,input
!! Define macro for MPI_Barrier
#DEFINE CALL_BARRIER CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
!! Macro for assertion
! #ifdef NDEBUG
! #define ASSERT(str,cond)
! #else
! #define ASSERT(str,cond) if (.not. cond) call abortdns(str)
! #endif
