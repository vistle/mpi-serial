
#include "mpiP.h"



/*
 * COLLECTIVE
 */


int FC_FUNC( mpi_barrier , MPI_BARRIER )(int *comm, int *ierror)
{
  *ierror=MPI_Barrier( *comm );
  return MPI_SUCCESS;
}


int MPI_Barrier(MPI_Comm comm )
{
  return(MPI_SUCCESS);
}

int MPI_Ibarrier(MPI_Comm comm, MPI_Request * request)
{
  Req *req;
  if (comm == MPI_COMM_NULL)
    return MPI_ERR_COMM;

  mpi_alloc_handle(request,(void **) &req);
  req->complete = 1;
  req->source=0;

  return(MPI_SUCCESS);
}


/*********/


int FC_FUNC( mpi_bcast , MPI_BCAST )(void *buffer, int *count, int *datatype,
				   int *root, int *comm, int *ierror )
{
  *ierror=MPI_Bcast(buffer, *count, *datatype, *root, *comm);
  return MPI_SUCCESS;
}



int MPI_Bcast(void* buffer, int count, MPI_Datatype datatype,
	      int root, MPI_Comm comm )
{
  if (root==MPI_ROOT)
    return(MPI_SUCCESS);

  if (root!=0)
    {
      fprintf(stderr,"MPI_Bcast: bad root = %d\n",root);
      abort();
    }


  return(MPI_SUCCESS);
}

int MPI_Ibcast(void *buffer, int count, MPI_Datatype datatype,
    int root, MPI_Comm comm, MPI_Request *request)
{
  Req *req;
  if (comm == MPI_COMM_NULL)
    return MPI_ERR_COMM;

  mpi_alloc_handle(request,(void **) &req);
  req->complete = 1;
  req->source=0;

  if (root==MPI_ROOT)
    return(MPI_SUCCESS);

  if (root!=0)
    {
      fprintf(stderr,"MPI_Bcast: bad root = %d\n",root);
      abort();
    }

  return(MPI_SUCCESS);
}

/*********/


int FC_FUNC( mpi_gather , MPI_GATHER )
                       (void *sendbuf, int *sendcount, int *sendtype,
			void *recvbuf, int *recvcount, int *recvtype,
			int *root, int *comm, int *ierror)
{
  *ierror=MPI_Gather( mpi_c_in_place(sendbuf), *sendcount, *sendtype,
		      recvbuf, *recvcount, *recvtype,
		      *root, *comm);
  return MPI_SUCCESS;
}


int MPI_Gather(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
	       void* recvbuf, int recvcount, MPI_Datatype recvtype,
	       int root, MPI_Comm comm)
{
  if (sendbuf==MPI_IN_PLACE)
    return(MPI_SUCCESS);

  if (root==MPI_ROOT)
    return(MPI_SUCCESS);

  if (root!=0)
    {
      fprintf(stderr,"MPI_Gather: bad root = %d\n",root);
      abort();
    }

  copy_data2(sendbuf, sendcount, sendtype,
             recvbuf, recvcount, recvtype);
//  memcpy(recvbuf,sendbuf,sendcount*sendtype);

  return(MPI_SUCCESS);
}

/*********/



int FC_FUNC( mpi_gatherv , MPI_GATHERV )
                        ( void *sendbuf, int *sendcount, int *sendtype,
			  void *recvbuf, int *recvcounts, int *displs,
			  int *recvtype, int *root, int *comm, int *ierror)
{
  *ierror=MPI_Gatherv( mpi_c_in_place(sendbuf), *sendcount, *sendtype,
		       recvbuf, recvcounts, displs,
		       *recvtype, *root, *comm);
  return MPI_SUCCESS;
}


int MPI_Gatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
		void* recvbuf, const int recvcounts[], const int displs[],
		MPI_Datatype recvtype, int root, MPI_Comm comm)
{
  int offset;
  MPI_Aint rt_extent;

  if (sendbuf==MPI_IN_PLACE)
    return(MPI_SUCCESS);

  if (root==MPI_ROOT)
    return(MPI_SUCCESS);

  if (root!=0)
    {
      fprintf(stderr,"MPI_Gatherv: bad root = %d\n",root);
      abort();
    }

  MPI_Type_extent(recvtype, &rt_extent);
  offset=displs[0]*rt_extent;

  copy_data2(sendbuf, sendcount, sendtype,
             (char*)recvbuf+offset, recvcounts[0], recvtype);

//  memcpy( (char *)recvbuf+offset, sendbuf, recvcounts[0] * recvtype);

  return(MPI_SUCCESS);
}



/*********/


int FC_FUNC( mpi_allgather , MPI_ALLGATHER )
                          ( void *sendbuf, int *sendcount, int *sendtype,
			    void *recvbuf, int *recvcount, int *recvtype,
			    int *comm, int *ierror)
{
  *ierror=MPI_Allgather( mpi_c_in_place(sendbuf), *sendcount, *sendtype,
			 recvbuf, *recvcount, *recvtype,
			 *comm );
  return MPI_SUCCESS;
}


int MPI_Allgather(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
		  void* recvbuf, int recvcount, MPI_Datatype recvtype,
		  MPI_Comm comm)
{
  if (sendbuf==MPI_IN_PLACE)
    return(MPI_SUCCESS);

  copy_data2(sendbuf, sendcount, sendtype,
             recvbuf, recvcount, recvtype);
//  memcpy(recvbuf,sendbuf,sendcount * sendtype);

  return(MPI_SUCCESS);

}


/*********/


int FC_FUNC( mpi_allgatherv , MPI_ALLGATHERV )
                          ( void *sendbuf, int *sendcount, int *sendtype,
			    void *recvbuf, int *recvcounts, int *displs,
                            int *recvtype, int *comm, int *ierror)
{
  *ierror=MPI_Allgatherv( mpi_c_in_place(sendbuf), *sendcount, *sendtype,
			  recvbuf, recvcounts, displs,
                          *recvtype, *comm );
  return MPI_SUCCESS;
}


int MPI_Allgatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
		   void* recvbuf, const int recvcounts[], const int displs[],
                   MPI_Datatype recvtype, MPI_Comm comm)
{
  int offset;
  MPI_Aint rt_extent;

  if (sendbuf==MPI_IN_PLACE)
    return(MPI_SUCCESS);

  MPI_Type_extent(recvtype, &rt_extent);
  offset=displs[0]*rt_extent;

  copy_data2(sendbuf, sendcount, sendtype,
             (char*)recvbuf+offset, recvcounts[0], recvtype);

//  memcpy( (char *)recvbuf+offset, sendbuf, recvcounts[0] * recvtype);

  return(MPI_SUCCESS);
}


/*********/

/* MPI_Scatter
 * Scattering to one proc involves only one copy operation, so copy
 * data from source to dest pointer
 */

int FC_FUNC( mpi_scatter, MPI_SCATTER )
                         ( void *sendbuf, int *sendcount, int *sendtype,
			 void *recvbuf, int *recvcount, int *recvtype,
			 int *root, int *comm, int *ierror)
{
  *ierror = MPI_Scatter(sendbuf, *sendcount, *sendtype,
  			mpi_c_in_place(recvbuf), *recvcount, *recvtype,
			*root, *comm);
  return MPI_SUCCESS;
}

int MPI_Scatter(const void * sendbuf, int sendcount, MPI_Datatype sendtype,
		void * recvbuf, int recvcount, MPI_Datatype recvtype,
		int root, MPI_Comm comm)
{
  if (recvbuf==MPI_IN_PLACE)
    return(MPI_SUCCESS);

  if (root==MPI_ROOT)
    return(MPI_SUCCESS);

  if (root!=0)
    {
      fprintf(stderr,"MPI_Scatter: bad root = %d\n",root);
      abort();
    }

  copy_data2(sendbuf, sendcount, sendtype,
             recvbuf, recvcount, recvtype);

  return(MPI_SUCCESS);
}



/*********/


int FC_FUNC( mpi_scatterv , MPI_SCATTERV )
                         ( void *sendbuf, int *sendcounts, int *displs,
			   int *sendtype, void *recvbuf, int *recvcount,
			   int *recvtype, int *root, int *comm, int *ierror)
{
  *ierror=MPI_Scatterv(sendbuf, sendcounts, displs,
		       *sendtype, mpi_c_in_place(recvbuf), *recvcount,
		       *recvtype, *root, *comm);
  return MPI_SUCCESS;
}



int MPI_Scatterv(const void* sendbuf, const int sendcounts[], const int displs[],
		 MPI_Datatype sendtype, void* recvbuf, int recvcount,
		 MPI_Datatype recvtype, int root, MPI_Comm comm)
{
  int offset;
  MPI_Aint st_extent;

  if (recvbuf==MPI_IN_PLACE)
    return(MPI_SUCCESS);

  if (root==MPI_ROOT)
    return(MPI_SUCCESS);

  if (root!=0)
    {
      fprintf(stderr,"MPI_Scatterv: bad root = %d\n",root);
      abort();
    }
  MPI_Type_extent(sendtype, &st_extent);
  offset=displs[0]*st_extent;

  copy_data2((char*)sendbuf+offset, sendcounts[0], sendtype,
             recvbuf, recvcount, recvtype);
//  memcpy(recvbuf,(char *)sendbuf+offset,sendcounts[0] * sendtype);

  return(MPI_SUCCESS);
}



/*********/


int FC_FUNC( mpi_reduce , MPI_REDUCE )
                       ( void *sendbuf, void *recvbuf, int *count,
			 int *datatype, int *op, int *root, int *comm,
			 int *ierror)
{
  *ierror=MPI_Reduce(sendbuf, recvbuf, *count,
		     *datatype, *op, *root, *comm);
  return MPI_SUCCESS;
}



int MPI_Reduce(const void* sendbuf, void* recvbuf, int count,
	       MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)

{
  if (root!=0)
    {
      fprintf(stderr,"MPI_Reduce: bad root = %d\n",root);
      abort();
    }

  copy_data2(sendbuf, count, datatype, recvbuf, count, datatype);
//  memcpy(recvbuf,sendbuf,count * datatype);

  return(MPI_SUCCESS);
}


/*********/


int FC_FUNC( mpi_allreduce , MPI_ALLREDUCE )
                          ( void *sendbuf, void *recvbuf, int *count,
			    int *datatype, int *op, int *comm, int *ierror)
{
  *ierror=MPI_Allreduce(sendbuf, recvbuf, *count,
			*datatype, *op, *comm);
  return MPI_SUCCESS;

}


int MPI_Allreduce(const void* sendbuf, void* recvbuf, int count,
		  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  if (sendbuf==MPI_IN_PLACE)
    return(MPI_SUCCESS);

  copy_data2(sendbuf, count, datatype, recvbuf, count, datatype);
//  memcpy(recvbuf,sendbuf,count * datatype);

  return(MPI_SUCCESS);

}

extern int MPI_Iallreduce(const void *sendbuf, void *recvbuf, int count,
                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                   MPI_Request *request)
{
  Req *req;
  if (comm == MPI_COMM_NULL)
    return MPI_ERR_COMM;

  mpi_alloc_handle(request,(void **) &req);
  req->complete = 1;
  req->source=0;

  if (sendbuf==MPI_IN_PLACE)
    return(MPI_SUCCESS);

  copy_data2(sendbuf, count, datatype, recvbuf, count, datatype);
//  memcpy(recvbuf,sendbuf,count * datatype);

  return(MPI_SUCCESS);
}

/*********/


/* MPI_Reduce_scatter
 * Performs reduction of n*sum(recvcounts) and distributes to all members
 * in a group. We do this to only one proc, so recvcounts[0] is only used.
 */

int FC_FUNC(mpi_reduce_scatter, MPI_REDUCE_SCATTER)
                (void * sendbuf, void * recvbuf, int *recvcounts,
                 int *datatype, int *op, int *comm, int *ierr)
{
  *ierr = MPI_Reduce_scatter(sendbuf, recvbuf, recvcounts, *datatype, *op, *comm);
  return MPI_SUCCESS;
}


int MPI_Reduce_scatter(const void* sendbuf, void* recvbuf, const int recvcounts[],
                         MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  copy_data2(sendbuf, recvcounts[0], datatype, recvbuf, recvcounts[0], datatype);
  return MPI_SUCCESS;
}


/*********/


int FC_FUNC( mpi_scan , MPI_SCAN)
                       ( void *sendbuf, void *recvbuf, int *count,
                         int *datatype, int *op, int *comm,
                         int *ierror)
{
  *ierror=MPI_Scan( sendbuf, recvbuf, *count,
                    *datatype, *op, *comm);
  return MPI_SUCCESS;
}



int MPI_Scan(const void* sendbuf, void* recvbuf, int count,
             MPI_Datatype datatype, MPI_Op op, MPI_Comm comm )
{
    copy_data2(sendbuf, count, datatype, recvbuf, count, datatype);

    return(MPI_SUCCESS);
}

/*********/


int FC_FUNC( mpi_alltoall , MPI_ALLTOALL )
                        ( void *sendbuf, int *sendcount, int *sendtype,
			  void *recvbuf, int *recvcount, int *recvtype,
                          int *comm, int *ierror )
{
  *ierror=MPI_Alltoall(sendbuf, *sendcount, *sendtype,
		       recvbuf, *recvcount, *recvtype,
		       *comm);
  return MPI_SUCCESS;
}


int MPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
		 void *recvbuf, int recvcount, MPI_Datatype recvtype,
		 MPI_Comm comm)
{
  copy_data2(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype);
//  memcpy(recvbuf,sendbuf,sendcount * sendtype);

  return(MPI_SUCCESS);
}


/*********/


int FC_FUNC( mpi_alltoallv , MPI_ALLTOALLV )
           ( void *sendbuf, int *sendcounts, int *sdispls, int *sendtype,
	     void *recvbuf, int *recvcounts, int *rdispls, int *recvtype,
             int *comm, int *ierror )
{

  *ierror=MPI_Alltoallv(sendbuf, sendcounts, sdispls, *sendtype,
			recvbuf, recvcounts, rdispls, *recvtype,
			*comm);

  return MPI_SUCCESS;
}

int MPI_Alltoallv(const void *sendbuf, const int sendcounts[],
		  const int sdispls[], MPI_Datatype sendtype,
                  void *recvbuf, const int recvcounts[],
		  const int rdispls[], MPI_Datatype recvtype,
                  MPI_Comm comm)

{
  int send_offset;
  int recv_offset;
  MPI_Aint st_extent;
  MPI_Aint rt_extent;

  MPI_Type_extent(sendtype, &st_extent);
  MPI_Type_extent(recvtype, &rt_extent);

  send_offset=sdispls[0]*st_extent;
  recv_offset=rdispls[0]*rt_extent;

  copy_data2((char*)sendbuf+send_offset, sendcounts[0], sendtype,
             (char*)recvbuf+recv_offset, recvcounts[0], recvtype);

//  memcpy( (char *)recvbuf+recv_offset, (char *)sendbuf+send_offset,
//	  sendcounts[0] * sendtype);


  return(MPI_SUCCESS);
}


/*********/


int FC_FUNC( mpi_alltoallw , MPI_ALLTOALLW )
           ( void *sendbuf, int *sendcounts, int *sdispls, int *sendtypes,
	     void *recvbuf, int *recvcounts, int *rdispls, int *recvtypes,
             int *comm, int *ierror )
{

  *ierror=MPI_Alltoallw(sendbuf, sendcounts, sdispls, sendtypes,
			recvbuf, recvcounts, rdispls, recvtypes,
			*comm);

  return MPI_SUCCESS;
}


int MPI_Alltoallw(const void *sendbuf, const int sendcounts[],
		  const int sdispls[], const MPI_Datatype sendtypes[],
                  void *recvbuf, const int recvcounts[],
		  const int rdispls[], const MPI_Datatype recvtypes[],
                  MPI_Comm comm)

{

  copy_data2((char*)sendbuf+sdispls[0], sendcounts[0], sendtypes[0],
             (char*)recvbuf+rdispls[0], recvcounts[0], recvtypes[0]);


  return(MPI_SUCCESS);
}



/*********/

MPI_Op MPI_Op_f2c(MPI_Fint op)
{
  return(op);
}


/*********/


MPI_Fint MPI_Op_c2f(MPI_Op op)
{
  return(op);
}
