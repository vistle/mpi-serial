#ifndef _MPI_H_
#define _MPI_H_

#include "export.h"
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

#define MPI_MAX_LIBRARY_VERSION_STRING (80)

typedef int V_MPIEXPORT MPI_Comm;
typedef int V_MPIEXPORT MPI_Request;
typedef int V_MPIEXPORT MPI_Message;


#define MPI_COMM_WORLD (1)
#define MPI_COMM_NULL (0)      /* handle 0 maps to NULL */


typedef int V_MPIEXPORT MPI_Group;

/* MPI_GROUP_EMPTY and MPI_GROUP_NULL must not conflict with MPI_GROUP_ONE */
#define MPI_GROUP_EMPTY (-1)
#define MPI_GROUP_NULL  (0)


/*
 * Return codes
 *   On error, mpi-serial aborts so the values don't really matter
 *   as long as they are different than MPI_SUCCESS
 *
 */

#define MPI_SUCCESS        (0)
#define MPI_ERR_BUFFER     (-1)
#define MPI_ERR_COUNT      (-1)
#define MPI_ERR_TYPE       (-1)
#define MPI_ERR_TAG        (-1)
#define MPI_ERR_COMM       (-1)
#define MPI_ERR_RANK       (-1)
#define MPI_ERR_REQUEST    (-1)
#define MPI_ERR_ROOT       (-1)
#define MPI_ERR_GROUP      (-1)
#define MPI_ERR_OP         (-1)
#define MPI_ERR_TOPOLOGY   (-1)
#define MPI_ERR_DIMS       (-1)
#define MPI_ERR_ARG        (-1)
#define MPI_ERR_UNKNOWN    (-1)
#define MPI_ERR_TRUNCATE   (-1)
#define MPI_ERR_OTHER      (-1)
#define MPI_ERR_INTERN     (-1)
#define MPI_PENDING        (-1)
#define MPI_ERR_IN_STATUS  (-1)
#define MPI_ERR_LASTCODE   (-1)
#define MPI_ERR_NO_MEM     (-1)
#define MPI_ERR_KEYVAL     (-1)
#define MPI_ERR_PENDING    (-1)
#define MPI_ERR_IO         (-1)

/*
 * MPI_UNDEFINED
 *
 * Uses:
 *   value for "color" in e.g. comm_split
 *   value for rank in Group_translate_ranks
 *
 */


#define MPI_UNDEFINED (-1)

/*
 * Data types etc.
 */

typedef intptr_t V_MPIEXPORT MPI_Aint;
#define MPI_BOTTOM (0)
#define MPI_IN_PLACE (void *)(-1)
typedef int V_MPIEXPORT MPI_Datatype;

/*
 * Topology
 */
#define MPI_CART (1)
#define MPI_GRAPH (2)
#define MPI_DIST_GRAPH (3)


#define MPI_HOST 0
#define MPI_IO 0


/* The type's value is now a handle */

#define MPI_DATATYPE_NULL   (0)

//C types
#define MPI_CHAR            (-1)
#define MPI_SIGNED_CHAR     (-46)
#define MPI_WCHAR           (-47)
#define MPI_SHORT           (-2)
#define MPI_INT             (-3)
#define MPI_LONG            (-4)
#define MPI_UNSIGNED_CHAR   (-5)
#define MPI_UNSIGNED_SHORT  (-6)
#define MPI_UNSIGNED        (-7)
#define MPI_UNSIGNED_LONG   (-8)
#define MPI_FLOAT           (-9)
#define MPI_DOUBLE          (-10)
#define MPI_LONG_DOUBLE     (-11)

#define MPI_CXX_BOOL        MPI_CHAR

//Cross-language
#define MPI_BYTE            (-12)
#define MPI_PACKED          (-13)
#define MPI_LB              (-14)
#define MPI_UB              (-15)

// Fortran types
#define MPI_INTEGER           (-16)     // RML: why not (MPI_INT)
#define MPI_REAL              (-17)     // RML: why not (MPI_FLOAT)
#define MPI_DOUBLE_PRECISION  (-18)     // RML: why not (MPI_DOUBLE)

#define MPI_COMPLEX           (-19)
#define MPI_DOUBLE_COMPLEX    (-20)
#define MPI_LOGICAL           (-21)
#define MPI_CHARACTER         (-22)
#define MPI_2REAL             (-23)
#define MPI_2DOUBLE_PRECISION (-24)
#define MPI_2INTEGER          (-25)

//Reduction function types

#define MPI_FLOAT_INT       (-26)
#define MPI_DOUBLE_INT      (-27)
#define MPI_LONG_INT        (-28)
#define MPI_2INT            (-29)
#define MPI_SHORT_INT       (-30)
#define MPI_LONG_DOUBLE_INT (-31)


/* Fortran size-specific types */

#define MPI_INTEGER1       (-32)
#define MPI_INTEGER2       (-33)
#define MPI_INTEGER4       (-34)
#define MPI_INTEGER8       (-35)
#define MPI_INTEGER16      (-36)

#define MPI_REAL4          (-37)
#define MPI_REAL8          (-38)
#define MPI_REAL16         (-39)

#define MPI_COMPLEX8       (-40)
#define MPI_COMPLEX16      (-41)
#define MPI_COMPLEX32      (-42)

/* Some more types */

#define MPI_LONG_LONG_INT       (-43)
#define MPI_LONG_LONG           MPI_LONG_LONG_INT
#define MPI_UNSIGNED_LONG_LONG  (-44)

#define MPI_OFFSET              (-45)



/*
 * Fortran int size
 *
 */

typedef int V_MPIEXPORT MPI_Fint;



#define MPI_ANY_TAG (-1)
#define MPI_TAG_UB (INT_MAX)

#define MPI_ANY_SOURCE (-1)
#define MPI_PROC_NULL (-2)
#define MPI_ROOT (-3)

#define MPI_REQUEST_NULL (0)

#define MPI_MAX_ERROR_STRING (128)
#define MPI_MAX_PROCESSOR_NAME (128)

#define MPI_THREAD_SINGLE (0)
#define MPI_THREAD_FUNNELED (1)
#define MPI_THREAD_SERIALIZED (2)
#define MPI_THREAD_MULTIPLE (3)

/*
 * MPI_Status
 *
 * definition must be compatible with the mpif.h values for
 * MPI_STATUS_SIZE, MPI_SOURCE, MPI_TAG, and MPI_ERROR.
 *
 * Note: The type used for MPI_Status_int must be chosen to match
 * Fortran INTEGER.
 *
 */

typedef int V_MPIEXPORT MPI_Status_int;

typedef struct                  /* Fortran: INTEGER status(MPI_STATUS_SIZE) */
{
  MPI_Status_int MPI_SOURCE;    /* Fortran: status(MPI_SOURCE) */
  MPI_Status_int MPI_TAG;       /* Fortran: status(MPI_TAG) */
  MPI_Status_int MPI_ERROR;     /* Fortran: status(MPI_ERROR) */
  int            get_count;     /* Number specified for send */

} MPI_Status;


#define MPI_STATUS_IGNORE    ((MPI_Status *)0)
#define MPI_STATUSES_IGNORE  ((MPI_Status *)0)


/*
 * MPI Errhandling stubs (Not functional currently)
 */
typedef int V_MPIEXPORT MPI_Errhandler;

#define MPI_ERRORS_ARE_FATAL ((MPI_Errhandler)0)
#define MPI_ERRORS_RETURN    ((MPI_Errhandler)-1)


/*
 * Collective operations
 */


typedef int V_MPIEXPORT MPI_Op;

typedef void MPI_User_function( void *invec, void *inoutvec, int *len,
                                MPI_Datatype *datatype);

#define MPI_OP_NULL (0)

#define MPI_MAX     (0)
#define MPI_MIN     (0)
#define MPI_SUM     (0)
#define MPI_PROD    (0)
#define MPI_LAND    (0)
#define MPI_BAND    (0)
#define MPI_LOR     (0)
#define MPI_BOR     (0)
#define MPI_LXOR    (0)
#define MPI_BXOR    (0)
#define MPI_MAXLOC  (0)
#define MPI_MINLOC  (0)



#define MPI_STATUS_SIZE       (sizeof(MPI_Status) / sizeof(int))


/* NOTE: the C type MPI_Offset is NOT the same as MPI datatype MPI_OFFSET */
typedef long long int V_MPIEXPORT MPI_Offset;


/* info
 */

typedef int V_MPIEXPORT MPI_Info;         /* handle */

#define MPI_INFO_NULL (0)


/* communicator comparison
 */
#define MPI_IDENT 1
#define MPI_CONGRUENT 2
#define MPI_SIMILAR 3
#define MPI_UNEQUAL 0



/**********************************************************
 *
 * Note: if you need to regenerate the prototypes below,
 * you can use 'protify.awk' and paste the output here.
 *
 */
extern int V_MPIEXPORT MPI_Get_library_version(char *version, int *resultlen);

extern int V_MPIEXPORT MPI_Intercomm_create(MPI_Comm local_comm, int local_leader,
                          MPI_Comm peer_comm, int remote_leader,
                          int tag, MPI_Comm *newintercomm);
extern int V_MPIEXPORT MPI_Intercomm_merge(MPI_Comm intercomm, int high,
			       MPI_Comm *newintercomm);
extern int V_MPIEXPORT MPI_Cart_create(MPI_Comm comm_old, int ndims, int *dims,
                        int *periods, int reorder, MPI_Comm *comm_cart);
extern int V_MPIEXPORT MPI_Cart_get(MPI_Comm comm, int maxdims, int *dims,
                        int *periods, int *coords);
extern int V_MPIEXPORT MPI_Cart_coords(MPI_Comm comm, int rank, int maxdims,
                        int *coords);
extern int V_MPIEXPORT MPI_Dims_create(int nnodes, int ndims, int *dims);

extern int V_MPIEXPORT MPI_Graph_create(MPI_Comm comm_old, int nnodes, const int index[],
    const int edges[], int reorder, MPI_Comm *comm_graph);

extern int V_MPIEXPORT MPI_Barrier(MPI_Comm comm );
extern int V_MPIEXPORT MPI_Ibarrier(MPI_Comm comm, MPI_Request * request);
extern int V_MPIEXPORT MPI_Bcast(void* buffer, int count, MPI_Datatype datatype,
                     int root, MPI_Comm comm );
extern int V_MPIEXPORT MPI_Ibcast(void *buffer, int count, MPI_Datatype datatype,
    int root, MPI_Comm comm, MPI_Request *request);
extern int V_MPIEXPORT MPI_Gather(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
                      void* recvbuf, int recvcount, MPI_Datatype recvtype,
                      int root, MPI_Comm comm);
extern int V_MPIEXPORT MPI_Gatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
                       void* recvbuf, const int *recvcounts, const int *displs,
                       MPI_Datatype recvtype, int root, MPI_Comm comm);
extern int V_MPIEXPORT MPI_Allgather(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
                         void* recvbuf, int recvcount, MPI_Datatype recvtype,
                         MPI_Comm comm);
extern int V_MPIEXPORT MPI_Allgatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
                          void* recvbuf, const int *recvcounts, const int *displs,
                          MPI_Datatype recvtype, MPI_Comm comm);
extern int V_MPIEXPORT MPI_Scatter(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
                        void* recvbuf, int recvcount, MPI_Datatype recvtype,
                        int root, MPI_Comm comm);
extern int V_MPIEXPORT MPI_Scatterv(const void* sendbuf, int *sendcounts, int *displs,
                        MPI_Datatype sendtype, void* recvbuf, int recvcount,
                        MPI_Datatype recvtype, int root, MPI_Comm comm);
extern int V_MPIEXPORT MPI_Reduce(const void* sendbuf, void* recvbuf, int count,
                      MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
extern int V_MPIEXPORT MPI_Reduce_scatter(const void* sendbuf, void* recvbuf, int *recvcounts,
                         MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
extern int V_MPIEXPORT MPI_Allreduce(const void* sendbuf, void* recvbuf, int count,
                         MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
extern int V_MPIEXPORT MPI_Iallreduce(const void *sendbuf, void *recvbuf, int count,
                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                   MPI_Request *request);
extern int V_MPIEXPORT MPI_Scan(const void* sendbuf, void* recvbuf, int count,
                    MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
extern int V_MPIEXPORT MPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                        void *recvbuf, int recvcount, MPI_Datatype recvtype,
                        MPI_Comm comm);
extern int V_MPIEXPORT MPI_Alltoallv(const void *sendbuf, int *sendcounts,
                         int *sdispls, MPI_Datatype sendtype,
                         void *recvbuf, int *recvcounts,
                         int *rdispls, MPI_Datatype recvtype,
                         MPI_Comm comm) ;
extern int V_MPIEXPORT MPI_Alltoallw(const void *sendbuf, int *sendcounts,
                         int *sdispls, MPI_Datatype *sendtypes,
                         void *recvbuf, int *recvcounts,
                         int *rdispls, MPI_Datatype *recvtypes,
                         MPI_Comm comm) ;


extern int V_MPIEXPORT MPI_Op_create(MPI_User_function *function, int commute,
                         MPI_Op *op);
extern MPI_Op V_MPIEXPORT MPI_Op_f2c(MPI_Fint op);
extern MPI_Fint V_MPIEXPORT MPI_Op_c2f(MPI_Op op);
extern MPI_Comm V_MPIEXPORT mpi_comm_new(void);
extern int V_MPIEXPORT MPI_Op_free(MPI_Op *op);
extern int V_MPIEXPORT MPI_Comm_free(MPI_Comm *comm);
extern int V_MPIEXPORT MPI_Comm_size(MPI_Comm comm, int *size);
extern int V_MPIEXPORT MPI_Comm_rank(MPI_Comm comm, int *rank);
extern int V_MPIEXPORT MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm);
extern int V_MPIEXPORT MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm);
extern int V_MPIEXPORT MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm);
extern int V_MPIEXPORT MPI_Comm_group(MPI_Comm comm, MPI_Group *group);
extern int V_MPIEXPORT MPI_Comm_compare(MPI_Comm comm1, MPI_Comm comm2, int *result);
extern int V_MPIEXPORT MPI_Comm_get_attr(MPI_Comm comm, int comm_keyval, void *attribute_val, int *flag);
extern int V_MPIEXPORT MPI_Comm_test_inter(MPI_Comm comm, int *flag);
extern int V_MPIEXPORT MPI_Comm_set_errhandler(MPI_Comm comm, MPI_Errhandler errhandler);
extern int V_MPIEXPORT MPI_Topo_test(MPI_Comm comm, int *status);
extern int V_MPIEXPORT MPI_Attr_get(MPI_Comm comm, int keyval,void *attribute_val, int *flag );

extern MPI_Comm V_MPIEXPORT MPI_Comm_f2c(MPI_Fint comm);
extern MPI_Fint V_MPIEXPORT MPI_Comm_c2f(MPI_Comm comm);

extern int V_MPIEXPORT MPI_Group_incl(MPI_Group group, int n, int *ranks, MPI_Group *newgroup);
extern int V_MPIEXPORT MPI_Group_excl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup);
extern int V_MPIEXPORT MPI_Group_range_incl(MPI_Group group, int n, int ranges[][3],
                         MPI_Group *newgroup);
extern int V_MPIEXPORT MPI_Group_union(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup);
extern int V_MPIEXPORT MPI_Group_intersection(MPI_Group group1, MPI_Group group2,
                                  MPI_Group *newgroup);
extern int V_MPIEXPORT MPI_Group_difference(MPI_Group group1, MPI_Group group2,
                         MPI_Group *newgroup);
extern int V_MPIEXPORT MPI_Group_free(MPI_Group *group);
extern int V_MPIEXPORT MPI_Group_translate_ranks(MPI_Group group1, int n, int *ranks1,
                                     MPI_Group group2, int *ranks2);
extern int V_MPIEXPORT MPI_Group_compare(MPI_Comm group1, MPI_Comm group2, int *result);
extern int V_MPIEXPORT MPI_Group_rank(MPI_Group group, int *rank);
extern int V_MPIEXPORT MPI_Group_size(MPI_Group group, int *size);
extern MPI_Group V_MPIEXPORT MPI_Group_f2c(MPI_Fint group);
extern MPI_Fint V_MPIEXPORT MPI_Group_c2f(MPI_Group group);

extern int V_MPIEXPORT MPI_Init(int *argc, char **argv[]) ;
extern int V_MPIEXPORT MPI_Init_thread( int *argc, char **argv[], int required, int *provided );
extern int V_MPIEXPORT MPI_Query_thread(int *provided);
extern int V_MPIEXPORT MPI_Finalize(void);
extern int V_MPIEXPORT MPI_Finalized(int *flag);
extern int V_MPIEXPORT MPI_Abort(MPI_Comm comm, int errorcode);
extern int V_MPIEXPORT MPI_Error_string(int errorcode, char *string, int *resultlen);
extern int V_MPIEXPORT MPI_Get_processor_name(char *name, int *resultlen);
extern int V_MPIEXPORT MPI_Get_version(int *mpi_vers, int *mpi_subvers);
extern int V_MPIEXPORT MPI_Get_library_version(char *version, int *resultlen);
extern int V_MPIEXPORT MPI_Initialized(int *flag);
extern int V_MPIEXPORT MPI_Is_thread_main(int *flag);

extern int V_MPIEXPORT MPI_Info_create(MPI_Info *info);
extern int V_MPIEXPORT MPI_Info_set(MPI_Info info, char *key, char *value);

extern int V_MPIEXPORT MPI_Pack( void *inbuf, int incount, MPI_Datatype datatype,
                     void *outbuf, int outsize, int *position, MPI_Comm comm);
extern int V_MPIEXPORT MPI_Unpack( void *inbuf, int insize, int *position,
                       void *outbuf, int outcount, MPI_Datatype datatype,
                       MPI_Comm comm );
extern int V_MPIEXPORT MPI_Irecv(void *buf, int count, MPI_Datatype datatype,
                     int source, int tag, MPI_Comm comm, MPI_Request *request);
extern int V_MPIEXPORT MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source,
                    int tag, MPI_Comm comm, MPI_Status *status);

extern int V_MPIEXPORT MPI_Test(MPI_Request *request, int *flag, MPI_Status *status);
extern int V_MPIEXPORT MPI_Wait(MPI_Request *request, MPI_Status *status);
extern int V_MPIEXPORT MPI_Testany(int count,  MPI_Request *array_of_requests,
                       int *index, int *flag, MPI_Status *status);
extern int V_MPIEXPORT MPI_Waitany(int count, MPI_Request *array_of_requests,
                       int *index, MPI_Status *status);
extern int V_MPIEXPORT MPI_Testall(int count, MPI_Request *array_of_requests,
                       int *flag, MPI_Status *array_of_statuses);
extern int V_MPIEXPORT MPI_Waitall(int count, MPI_Request *array_of_requests,
                       MPI_Status *array_of_statuses);
extern MPI_Request V_MPIEXPORT MPI_Request_f2c(MPI_Fint request);
extern MPI_Fint V_MPIEXPORT MPI_Request_c2f(MPI_Request request);
extern int V_MPIEXPORT MPI_Testsome(int incount, MPI_Request *array_of_requests,
                        int *outcount, int *array_of_indices,
                        MPI_Status *array_of_statuses);
extern int V_MPIEXPORT MPI_Waitsome(int incount, MPI_Request *array_of_requests,
                        int *outcount, int *array_of_indices,
                        MPI_Status *array_of_statuses);
extern int V_MPIEXPORT MPI_Request_free(MPI_Request * req);
extern int V_MPIEXPORT MPI_Cancel(MPI_Request * req);
extern int V_MPIEXPORT MPI_Test_cancelled(const MPI_Status *status, int *flag);
extern int V_MPIEXPORT MPI_Isend(const void *buf, int count, MPI_Datatype datatype,
                     int dest, int tag, MPI_Comm comm, MPI_Request *request) ;
extern int V_MPIEXPORT MPI_Send(const void* buf, int count, MPI_Datatype datatype,
                    int dest, int tag, MPI_Comm comm);
extern int V_MPIEXPORT MPI_Ssend(const void* buf, int count, MPI_Datatype datatype,
                     int dest, int tag, MPI_Comm comm);
extern int V_MPIEXPORT MPI_Issend(const void *buf, int count, MPI_Datatype datatype, int dest,
               int tag, MPI_Comm comm, MPI_Request *request);
extern int V_MPIEXPORT MPI_Rsend(const void* buf, int count, MPI_Datatype datatype,
                     int dest, int tag, MPI_Comm comm);
extern int V_MPIEXPORT MPI_Irsend(const void *buf, int count, MPI_Datatype datatype,
                     int dest, int tag, MPI_Comm comm, MPI_Request *request) ;
extern int V_MPIEXPORT MPI_Sendrecv(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
                        int dest, int sendtag,
                        void *recvbuf, int recvcount, MPI_Datatype recvtype,
                        int source, int recvtag,
                        MPI_Comm comm, MPI_Status *status);

extern int V_MPIEXPORT MPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status *status);
extern int V_MPIEXPORT MPI_Iprobe(int source, int tag, MPI_Comm comm,
                      int *flag, MPI_Status *status);

//extern int V_MPIEXPORT MPI_Pack_size(int incount, MPI_Datatype type, MPI_Comm comm, MPI_Aint * size);
extern int V_MPIEXPORT MPI_Pack_size(int incount, MPI_Datatype type, MPI_Comm comm, int * size);

/* Error handling stub, not currently functional */
extern int V_MPIEXPORT MPI_Errhandler_set(MPI_Comm comm, MPI_Errhandler handle);
extern int V_MPIEXPORT MPI_Error_class(int errorcode, int *errorclass);


/* new type functions */
extern int V_MPIEXPORT MPI_Get_count(const MPI_Status *status, MPI_Datatype datatype, int *count);
extern int V_MPIEXPORT MPI_Get_elements(const MPI_Status *status, MPI_Datatype datatype, int *count);
extern int V_MPIEXPORT MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype);

extern int V_MPIEXPORT MPI_Type_vector(int count, int blocklen, int stride, MPI_Datatype oldtype,
                           MPI_Datatype *newtype);

extern int V_MPIEXPORT MPI_Type_hvector(int count, int blocklen, MPI_Aint stride,
                            MPI_Datatype oldtype, MPI_Datatype *newtype);

extern int V_MPIEXPORT MPI_Type_create_hvector(int count, int blocklen, MPI_Aint stride,
                            MPI_Datatype oldtype, MPI_Datatype *newtype);

extern int V_MPIEXPORT MPI_Type_indexed(int count, int *blocklens, int *displacements,
                            MPI_Datatype oldtype, MPI_Datatype *newtype);

extern int V_MPIEXPORT MPI_Type_create_subarray(int ndims, const int array_of_sizes[], const
int array_of_subsizes[], const int array_of_starts[], int order, MPI_Datatype
oldtype, MPI_Datatype *newtype);

extern int V_MPIEXPORT MPI_Type_create_indexed_block(int count, int blocklen, int *displacements,
                                  MPI_Datatype oldtype, MPI_Datatype *newtype);
extern int V_MPIEXPORT MPI_Type_hindexed(int count, int *blocklens, MPI_Aint *displacements,
                             MPI_Datatype oldtype, MPI_Datatype *newtype);
extern int V_MPIEXPORT MPI_Type_size(MPI_Datatype type, int * size);
extern int V_MPIEXPORT MPI_Type_struct(int count, int *blocklens, MPI_Aint *displacements,
                           MPI_Datatype *oldtypes, MPI_Datatype *newtype);
extern int V_MPIEXPORT MPI_Type_dup(MPI_Datatype oldtype, MPI_Datatype *newtype);

extern int V_MPIEXPORT MPI_Type_extent(MPI_Datatype datatype, MPI_Aint * extent);
extern int V_MPIEXPORT MPI_Type_commit(MPI_Datatype * datatype);
extern int V_MPIEXPORT MPI_Type_free(MPI_Datatype * datatype);
extern int V_MPIEXPORT MPI_Type_lb(MPI_Datatype datatype, MPI_Aint * lb);
extern int V_MPIEXPORT MPI_Type_ub(MPI_Datatype datatype, MPI_Aint * ub);
extern int V_MPIEXPORT MPI_Type_create_struct(int count, int array_of_blocklengths[],
     const MPI_Aint array_of_displacements[], const MPI_Datatype array_of_types[],
     MPI_Datatype *newtype);

#define MPI_WTIME_IS_GLOBAL 1
extern double V_MPIEXPORT MPI_Wtime(void);
extern double V_MPIEXPORT MPI_Wtick(void);

extern int V_MPIEXPORT MPI_Alloc_mem(MPI_Aint size, MPI_Info info, void *baseptr);
extern int V_MPIEXPORT MPI_Free_mem(void *baseptr);
extern int V_MPIEXPORT MPI_Get_address(const void *location, MPI_Aint *address);


extern int V_MPIEXPORT MPI_Mprobe(int source, int tag, MPI_Comm comm, MPI_Message *message, MPI_Status *status);
extern int V_MPIEXPORT MPI_Improbe(int source, int tag, MPI_Comm comm, int *flag, MPI_Message *message, MPI_Status *status);
extern int V_MPIEXPORT MPI_Mrecv(void *buf, int count, MPI_Datatype type, MPI_Message *message, MPI_Status *status);


#define MPI_ORDER_C 89
#define MPI_ORDER_FORTRAN 90

typedef int V_MPIEXPORT MPI_Win;
#define MPI_WIN_NULL 0

#define MPI_MODE_NOCHECK 99

#define MPI_NO_OP 222
#define MPI_REPLACE 223

extern int V_MPIEXPORT MPI_Win_create(void *base, MPI_Aint size, int disp_unit, MPI_Info info, 
                  MPI_Comm comm, MPI_Win *win);
extern int V_MPIEXPORT MPI_Win_allocate(MPI_Aint size, int disp_unit, MPI_Info info,
                     MPI_Comm comm, void *baseptr, MPI_Win * win);
extern int V_MPIEXPORT MPI_Win_fence(int assert, MPI_Win win);
extern int V_MPIEXPORT MPI_Win_lock(int lock_type, int rank, int assert, MPI_Win win);
extern int V_MPIEXPORT MPI_Win_lock_all(int assert, MPI_Win win);
extern int V_MPIEXPORT MPI_Win_unlock(int rank, MPI_Win win);
extern int V_MPIEXPORT MPI_Win_unlock_all(MPI_Win win);
extern int V_MPIEXPORT MPI_Win_sync(MPI_Win win) ;
extern int V_MPIEXPORT MPI_Win_flush(int rank, MPI_Win win);
extern int V_MPIEXPORT MPI_Win_flush_all(MPI_Win win);
extern int V_MPIEXPORT MPI_Win_flush_local(int rank, MPI_Win win);
extern int V_MPIEXPORT MPI_Win_flush_local_all(MPI_Win win);
extern int V_MPIEXPORT MPI_Win_free(MPI_Win *win);
extern int V_MPIEXPORT MPI_Get(void *origin_addr, int origin_count, MPI_Datatype
            origin_datatype, int target_rank, MPI_Aint target_disp,
            int target_count, MPI_Datatype target_datatype, MPI_Win
            win);
extern int V_MPIEXPORT MPI_Put(const void *origin_addr, int origin_count, MPI_Datatype
            origin_datatype, int target_rank, MPI_Aint target_disp,
            int target_count, MPI_Datatype target_datatype, MPI_Win win);
extern int V_MPIEXPORT MPI_Fetch_and_op(const void *origin_addr, void *result_addr,
        MPI_Datatype datatype, int target_rank, MPI_Aint target_disp,
        MPI_Op op, MPI_Win win);

typedef int V_MPIEXPORT MPI_File;
#define MPI_FILE_NULL 0

#define MPI_MODE_RDONLY 100
#define MPI_MODE_RDWR 101
#define MPI_MODE_WRONLY 102
#define MPI_MODE_CREATE 103
#define MPI_MODE_EXCL 104
#define MPI_MODE_DELETE_ON_CLOSE 105
#define MPI_MODE_UNIQUE_OPEN 106
#define MPI_MODE_SEQUENTIAL 107
#define MPI_MODE_APPEND 108

extern int V_MPIEXPORT MPI_File_open(MPI_Comm comm, const char *filename,
    int amode, MPI_Info info,
    MPI_File *fh);
extern int V_MPIEXPORT MPI_File_get_size(MPI_File fh, MPI_Offset *size);
extern int V_MPIEXPORT MPI_File_write(MPI_File fh, const void *buf, int count,
                   MPI_Datatype datatype, MPI_Status *status);
extern int V_MPIEXPORT MPI_File_write_at(MPI_File fh, MPI_Offset offset, const void *buf,
                      int count, MPI_Datatype datatype, MPI_Status *status);
extern int V_MPIEXPORT MPI_File_write_all(MPI_File fh, const void *buf, int count,
                       MPI_Datatype datatype, MPI_Status * status);
extern int V_MPIEXPORT MPI_File_write_at_all(MPI_File fh, MPI_Offset offset, const void *buf,
     int count, MPI_Datatype datatype, MPI_Status *status);
extern int V_MPIEXPORT MPI_File_set_size(MPI_File fh, MPI_Offset size);
extern int V_MPIEXPORT MPI_File_read(MPI_File fh, void *buf, int count,
                  MPI_Datatype datatype, MPI_Status *status);
extern int V_MPIEXPORT MPI_File_read_all(MPI_File fh, void *buf, int count,
                      MPI_Datatype datatype, MPI_Status *status);                  
extern int V_MPIEXPORT MPI_File_read_at(MPI_File fh, MPI_Offset offset, void *buf,
                    int count, MPI_Datatype datatype, MPI_Status * status);
                    extern int V_MPIEXPORT MPI_File_write_at_all(MPI_File fh, MPI_Offset offset, const void *buf,
    int count, MPI_Datatype datatype, MPI_Status *status);
extern int V_MPIEXPORT MPI_File_read_at_all(MPI_File fh, MPI_Offset offset, void *buf,
                         int count, MPI_Datatype datatype, MPI_Status * status);    
extern int V_MPIEXPORT MPI_File_close(MPI_File * fh);
extern int V_MPIEXPORT MPI_File_set_view(MPI_File fh, MPI_Offset disp,
     MPI_Datatype etype, MPI_Datatype filetype,
     const char *datarep, MPI_Info info);


#ifdef __cplusplus
}
#endif

#endif
