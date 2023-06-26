#include <unistd.h>
//#include <mpi.h>
#include <stdio.h>
#include <pthread.h>
#include <time.h>
#include <errno.h> 

const long t_pause = 150;

int endThread;

void * check_mpi(void * arg){
   int flag;
   while(endThread==0){
      int flag;
      //MPI_Iprobe(MPI_ANY_SOURCE, 1336, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
      fortran_check_mpi();
      
      usleep(1000*t_pause);
   }
}

void start_prog_thread(pthread_t *threadId){
   // printf("started progress thread. Pause time: %ld\n", t_pause);
   endThread = 0;
   pthread_create(threadId, NULL, &check_mpi, NULL);
}

void stop_prog_thread(pthread_t *threadId){
   endThread= 1;
   usleep(3*1000*t_pause);
   pthread_cancel(*threadId);
   // printf("stopped progress thread\n");
}

