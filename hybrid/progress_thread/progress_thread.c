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

