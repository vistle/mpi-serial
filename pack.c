#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpiP.h"
#include "type.h"

/*
 *
 */


int FC_FUNC( mpi_pack , MPI_PACK )
     ( void *inbuf, int *incount, int *datatype,
       void *outbuf, int *outsize, int *position, int *comm, int *ierror)
{
  *ierror=MPI_Pack(inbuf, *incount,* datatype,
  	           outbuf, *outsize, position, *comm);
  return MPI_SUCCESS;
}



int MPI_Pack(const void *inbuf, int incount, MPI_Datatype datatype,
              void *outbuf, int outsize, int *position, MPI_Comm comm)
{
  int ret;

  Datatype type_ptr = *(Datatype*) mpi_handle_to_datatype(datatype);
  Comm* comm_ptr = mpi_handle_to_ptr(comm);

  ret = Pack(inbuf, incount, type_ptr, outbuf, outsize, position, comm_ptr);

  return ret;
}



int Pack(const void *inbuf, int incount, Datatype type,
              void *outbuf, int outsize, int *position, Comm * comm)
{
  int i, j;
  MPI_Aint extent;
  //check that buffer is large enough
  Type_extent(type, &extent);
  for (i = 0; i < incount; i++)
  {
    for (j = 0; j < type->count; j++)
    {
      if ((*position) + Simpletype_length(type->pairs[j].type) > outsize)
      {
        printf("MPI_Pack: data exceeds buffer size\n");
        exit(1);
      }
      memcpy(((char*) outbuf)+(*position), (char*)inbuf+type->pairs[j].disp + (extent*i),
             Simpletype_length(type->pairs[j].type));
      *position += Simpletype_length(type->pairs[j].type);
    }
  }
  return MPI_SUCCESS;
}

int Pack_size(int incount, Datatype datatype,
              Comm * comm, MPI_Aint * size)
{
    int i;
    *size = 0;
    //sum up all sizes
    for(i = 0; i < datatype->count; i++)
    {
        *size += Simpletype_length(datatype->pairs[i].type);
    }
    *size *= incount;
    printf("Size = %ld\n", (long)*size);
  return MPI_SUCCESS;
}

int FC_FUNC( mpi_pack_size, MPI_PACK_SIZE )(int * incount, int * datatype,
                                          int * comm, int * size, int *ierr)
{
  *ierr = MPI_Pack_size(*incount, *datatype, *comm, size);
  return MPI_SUCCESS;
}

#if 0
int MPI_Pack_size(int incount, MPI_Datatype datatype,
                  MPI_Comm comm, MPI_Aint * size)
{
  int ret;
  Datatype type_ptr = *(Datatype*) mpi_handle_to_datatype(datatype);
  Comm * comm_ptr = mpi_handle_to_ptr(comm);

  ret = Pack_size(incount, type_ptr, comm_ptr, size);

  return ret;
}
#else
int MPI_Pack_size(int incount, MPI_Datatype datatype,
                  MPI_Comm comm, int * size)
{
  int ret;
  Datatype type_ptr = *(Datatype*) mpi_handle_to_datatype(datatype);
  Comm * comm_ptr = mpi_handle_to_ptr(comm);
  MPI_Aint sz;

  ret = Pack_size(incount, type_ptr, comm_ptr, &sz);
  *size = sz;

  return ret;
}
#endif



/*
 *
 */


int FC_FUNC( mpi_unpack , MPI_UNPACK )
     ( void *inbuf, int *insize, int *position,
       void *outbuf, int *outcount, int *datatype,
       int *comm, int *ierror )
{
  *ierror=MPI_Unpack( inbuf, *insize, position,
                      outbuf, *outcount, *datatype, *comm);
  return MPI_SUCCESS;
}


int MPI_Unpack(const void * inbuf, int insize, int * position, void * outbuf,
               int outcount, MPI_Datatype type, MPI_Comm comm)
{
  int ret;
  Datatype type_ptr = *(Datatype*) mpi_handle_to_datatype(type);
  Comm * comm_ptr = mpi_handle_to_ptr(comm);

  ret = Unpack(inbuf, insize, position, outbuf, outcount, type_ptr, comm_ptr);

  return ret;
}

int Unpack(const void * inbuf, int insize, int * position, void *outbuf,
                int outcount, Datatype type, Comm* comm)
{
  int i, j;
  MPI_Aint extent;

  Type_extent(type, &extent);

  for (i = 0; i < outcount; i++)
  {
    for (j = 0; j < type->count; j++)
    {
      if ((*position) + Simpletype_length(type->pairs[j].type) > insize)
      {
        printf("MPI_Unpack: Data exceeds buffer size\n");
	exit(1);
      }
      memcpy((char*)outbuf+type->pairs[j].disp + (extent*i), ((char*) inbuf)+(*position) ,
             Simpletype_length(type->pairs[j].type));
      *position += Simpletype_length(type->pairs[j].type);
    }
  }
  return MPI_SUCCESS;
}


