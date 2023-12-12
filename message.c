#include "mpiP.h"

int MPI_Mprobe(int source, int tag, MPI_Comm comm, MPI_Message *message, MPI_Status *status)
{
    abort();
    return MPI_SUCCESS;
}

int MPI_Improbe(int source, int tag, MPI_Comm comm, int *flag, MPI_Message *message, MPI_Status *status)
{
    abort();
    return MPI_SUCCESS;
}

int MPI_Mrecv(void *buf, int count, MPI_Datatype type, MPI_Message *message, MPI_Status *status)
{
    abort();
    return MPI_SUCCESS;
}
