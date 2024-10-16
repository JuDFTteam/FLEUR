#include <unistd.h>
//#include <mpi.h>
#include <stdio.h>
#include <pthread.h>
#include <time.h>
#include <errno.h> 

#ifdef CPP_MPI

const long t_pause = 150;
extern void fortran_check_mpi();
int endThread;

void * check_mpi(void * arg){
   while(endThread==0){
      fortran_check_mpi();
      usleep(1000*t_pause);
   }
}

void start_prog_thread(pthread_t *threadId){
   endThread = 0;
   pthread_create(threadId, NULL, &check_mpi, NULL);
}

void stop_prog_thread(pthread_t *threadId){
   endThread= 1;
}

#endif
