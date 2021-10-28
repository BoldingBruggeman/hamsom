c-----------------------------------------------------------------------
c     mpi variables
c-----------------------------------------------------------------------
       include "mpif.h"
       PARAMETER (nPh=3,nPv=6,nprocs=nPh*nPv)
c       PARAMETER (nPh=2,nPv=2,nprocs=nPh*nPv)
       common /mpivar/  myid,ierr,status(100,MPI_STATUS_SIZE)
