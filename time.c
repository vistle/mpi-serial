
#include <time.h>
#include "mpiP.h"
#include "export.h"
#ifdef _WIN32
#include <windows.h>
#include <sysinfoapi.h>
#include <stdint.h>

struct timeval {
  long tv_sec;  // Seconds
  long tv_usec; // Microseconds
};

struct timezone {
    int tz_minuteswest; // Minutes west of GMT
    int tz_dsttime;     // Daylight saving time flag
};

int gettimeofday(struct timeval *tv, struct timezone *tz) {
    if (tv) {
        FILETIME ft;
        ULARGE_INTEGER li;
        uint64_t epoch_time;

        // Get the current system time
        GetSystemTimeAsFileTime(&ft);

        // Convert FILETIME to a 64-bit integer
        li.LowPart = ft.dwLowDateTime;
        li.HighPart = ft.dwHighDateTime;

        // Convert to Unix epoch time (subtract Windows epoch offset)
        epoch_time = (li.QuadPart - 116444736000000000ULL) / 10;

        tv->tv_sec = (long)(epoch_time / 1000000ULL);
        tv->tv_usec = (long)(epoch_time % 1000000ULL);
    }

    if (tz) {
        // Timezone information is not supported on Windows
        tz->tz_minuteswest = 0;
        tz->tz_dsttime = 0;
    }

    return 0;
}

#else
#include <sys/time.h>
#endif




double MPI_Wtime(void);
double MPI_Wtick(void);



double FC_FUNC( mpi_wtime, MPI_WTIME )(void)
{
  return(MPI_Wtime());
}



double MPI_Wtime(void)
{
  struct timeval tv;

  if (gettimeofday(&tv,0))
    {
      fprintf(stderr,"MPI_Wtime: error calling gettimeofday()\n");
      abort();
    }


  return((double)(tv.tv_sec) + (double)(tv.tv_usec)/1e6) ;
}

double MPI_Wtick(void)
{
  return 1.0f/(1000000.0f);
}
