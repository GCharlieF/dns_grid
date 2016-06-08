!!Macros for fortran pre-processor
#DEFINE MASTER IF (proc_id==0)
!! Define macro for screen output
#DEFINE PPRINT(input) IF (proc_id==0) PRINT *,input
!! Define macro for MPI_Barrier
#DEFINE CALL_BARRIER CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
!! Macro for assertion
#ifdef DEBUG
#define ASSERT(str,cond)
#else
#DEFINE ASSERT(str,cond) if (.not. cond) print *,'Failed assert ::',str
#endif
