#include "mpiP.h"

int MPI_Win_create(void *base, MPI_Aint size, int disp_unit, MPI_Info info, 
                  MPI_Comm comm, MPI_Win *win)
{
    abort();
    return MPI_SUCCESS;
}

int MPI_Win_allocate(MPI_Aint size, int disp_unit, MPI_Info info,
                     MPI_Comm comm, void *baseptr, MPI_Win * win)
{
    abort();
    return MPI_SUCCESS;
}

int MPI_Win_fence(int assert, MPI_Win win)
{
    abort();
    return MPI_SUCCESS;
}

int MPI_Win_lock(int lock_type, int rank, int assert, MPI_Win win)
{
    abort();
    return MPI_SUCCESS;
}

int MPI_Win_lock_all(int assert, MPI_Win win)
{
    abort();
    return MPI_SUCCESS;
}

int MPI_Win_unlock(int rank, MPI_Win win)
{
    abort();
    return MPI_SUCCESS;
}

int MPI_Win_unlock_all(MPI_Win win)
{
    abort();
    return MPI_SUCCESS;
}

int MPI_Win_sync(MPI_Win win)
{
    abort();
    return MPI_SUCCESS;
}

int MPI_Win_flush(int rank, MPI_Win win)
{
    abort();
    return MPI_SUCCESS;
}

int MPI_Win_flush_all(MPI_Win win)
{
    abort();
    return MPI_SUCCESS;
}

int MPI_Win_flush_local(int rank, MPI_Win win)
{
    abort();
    return MPI_SUCCESS;
}

int MPI_Win_flush_local_all(MPI_Win win)
{
    abort();
    return MPI_SUCCESS;
}

int MPI_Win_free(MPI_Win *win)
{
    abort();
    return MPI_SUCCESS;
}

int MPI_Get(void *origin_addr, int origin_count, MPI_Datatype
            origin_datatype, int target_rank, MPI_Aint target_disp,
            int target_count, MPI_Datatype target_datatype, MPI_Win
            win)
{
    abort();
    return MPI_SUCCESS;
}

int MPI_Put(const void *origin_addr, int origin_count, MPI_Datatype
            origin_datatype, int target_rank, MPI_Aint target_disp,
            int target_count, MPI_Datatype target_datatype, MPI_Win win)
{
    abort();
    return MPI_SUCCESS;
}

int MPI_Fetch_and_op(const void *origin_addr, void *result_addr,
        MPI_Datatype datatype, int target_rank, MPI_Aint target_disp,
        MPI_Op op, MPI_Win win)
{
    abort();
    return MPI_SUCCESS;
}
