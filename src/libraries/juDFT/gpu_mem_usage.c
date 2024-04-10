
#ifdef CPP_GPU
#include "cuda.h"


double gpu_mem_usage(){ 
 
    size_t free_byte ;
    size_t total_byte ;
    cudaMemGetInfo( &free_byte, &total_byte ) ;

    return ((double)total_byte-(double)free_byte)/1024/1024/1024 ;

}
#endif