#include <unistd.h>
//#include <mpi.h>
#include <stdio.h>
#include <pthread.h>
#include <time.h>
#include <errno.h> 

const long t_pause = 150;

int msleep(long msec){
   struct timespec ts;
   int res;

   if (msec < 0){
      errno = EINVAL;
      return -1;
   }

   ts.tv_sec = msec / 1000;
   ts.tv_nsec = (msec % 1000) * 1000000;

   do {
      res = nanosleep(&ts, &ts);
   } while (res && errno == EINTR);
   return res;
}

void * check_mpi(void * arg){
   int flag;
   while(1){
      int flag;
      //MPI_Iprobe(MPI_ANY_SOURCE, 1336, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
      fortran_check_mpi();
      msleep(t_pause);
   }
}

void start_prog_thread(pthread_t *threadId){
   // printf("started progress thread. Pause time: %ld\n", t_pause);
   pthread_create(threadId, NULL, &check_mpi, NULL);
}

void stop_prog_thread(pthread_t *threadId){
   pthread_cancel(*threadId);
   // printf("stopped progress thread\n");
}

