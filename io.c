#include "mpiP.h"

#include <errno.h>
#include <fcntl.h>
#include <unistd.h>

static int mode(int amode) {
    return amode;
}

static int errno2ret() {
    return MPI_ERR_IO;
}

int MPI_File_open(MPI_Comm comm, const char *filename,
    int amode, MPI_Info info,
    MPI_File *fh)
{
    if (comm == MPI_COMM_NULL)
        return MPI_ERR_COMM;

    int fd = open(filename, mode(amode));
    if (fd < 0) {
        return errno2ret();
    }

    *fh = fd;
    return MPI_SUCCESS;
}

int MPI_File_close(MPI_File * fh)
{
    int fd = *fh;
    int err = close(fd);
    if (err == -1) {
        return errno2ret();
    }
    return MPI_SUCCESS;
}

int MPI_File_get_size(MPI_File fh, MPI_Offset *size);
int MPI_File_set_size(MPI_File fh, MPI_Offset size);

int MPI_File_write(MPI_File fh, const void *buf, int count,
                   MPI_Datatype datatype, MPI_Status *status)
{
    abort();
    return MPI_SUCCESS;
}

int MPI_File_write_at(MPI_File fh, MPI_Offset offset, const void *buf,
                      int count, MPI_Datatype datatype, MPI_Status *status);
int MPI_File_write_all(MPI_File fh, const void *buf, int count,
                       MPI_Datatype datatype, MPI_Status * status)
{
    abort();
    return MPI_SUCCESS;
}

int MPI_File_write_at_all(MPI_File fh, MPI_Offset offset, const void *buf,
     int count, MPI_Datatype datatype, MPI_Status *status);
int MPI_File_read(MPI_File fh, void *buf, int count,
                  MPI_Datatype datatype, MPI_Status *status);
int MPI_File_read_all(MPI_File fh, void *buf, int count,
                      MPI_Datatype datatype, MPI_Status *status);                  
int MPI_File_read_at(MPI_File fh, MPI_Offset offset, void *buf,
                    int count, MPI_Datatype datatype, MPI_Status * status);
                    extern int MPI_File_write_at_all(MPI_File fh, MPI_Offset offset, const void *buf,
    int count, MPI_Datatype datatype, MPI_Status *status);
int MPI_File_read_at_all(MPI_File fh, MPI_Offset offset, void *buf,
                         int count, MPI_Datatype datatype, MPI_Status * status);    
int MPI_File_set_view(MPI_File fh, MPI_Offset disp,
     MPI_Datatype etype, MPI_Datatype filetype,
     const char *datarep, MPI_Info info)
{
    abort();
    return MPI_SUCCESS;
}
