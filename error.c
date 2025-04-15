
#include "mpiP.h"

/*
 * Error handling code
 * Just a stub for now to support the MPI interface without actually
 * doing anything
 */

 int MPI_Errhandler_set(MPI_Comm comm, MPI_Errhandler handle)
 {
   return(MPI_SUCCESS);
 }

 int MPI_Error_class(int errorcode, int *errorclass)
 {
   *errorclass = errorcode; // similar to OpenMPI
   return(MPI_SUCCESS);
 }
