!!!.....................................................................
!!!..........................MPI MOD....................................
!!!.....................................................................
MODULE mpi_mod
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  logical :: iammaster
  integer(kind=MPI_INTEGER_KIND) :: myid_world, numprocs_world
  integer(kind=MPI_INTEGER_KIND) :: myid, numprocs, master, mpi_err
CONTAINS
!!!. . . . . . . . . . . . . MPI INIT . . . . . . . . . . . . . . . . . .
SUBROUTINE MPI_initialization
  CALL MPI_INIT(mpi_err)
  CALL MPI_Comm_size(MPI_COMM_WORLD,numprocs_world,mpi_err)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,myid_world,mpi_err)
END SUBROUTINE MPI_initialization
!!!. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
END MODULE mpi_mod
!!!.....................................................................
!!!.....................................................................
