/*
!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>

/* Interface for the cusolverDN routines for solving a generalized Eigenvalue problem
Code adopted from the example in the documentation
*/

void cusolver_complex(cuDoubleComplex *H,cuDoubleComplex *S,int n,int ne,double tol,int max_sweeps,double* eig,cuDoubleComplex *z){

  cusolverDnHandle_t cusolverH = NULL;
  cudaStream_t stream = NULL;
  syevjInfo_t syevj_params = NULL;

  cusolverStatus_t status = CUSOLVER_STATUS_SUCCESS;
  cudaError_t cudaStat1 = cudaSuccess;
  cudaError_t cudaStat2 = cudaSuccess;
  cudaError_t cudaStat3 = cudaSuccess;

  double *d_W = NULL; /* eigenvalues on device*/
  int *d_info = NULL; /* error info */
  int  lwork = 0;     /* size of workspace */
  cuDoubleComplex *d_work = NULL; /* device workspace for syevj */
  int info = 0;       /* host copy of error info */

  /* configuration of syevj  */
  const cusolverEigType_t itype = CUSOLVER_EIG_TYPE_1; //Solve H psi=S lambda psi
  const cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR; // compute eigenvectors.
  const cublasFillMode_t  uplo = CUBLAS_FILL_MODE_LOWER;

  /* numerical results of syevj  */
  double residual = 0;
  int executed_sweeps = 0;
  
  /* step 1: create cusolver handle, bind a stream  */
  status = cusolverDnCreate(&cusolverH);
  assert(CUSOLVER_STATUS_SUCCESS == status);

  cudaStat1 = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
  assert(cudaSuccess == cudaStat1);

  status = cusolverDnSetStream(cusolverH, stream);
  assert(CUSOLVER_STATUS_SUCCESS == status);

  /* step 2: configuration of syevj */
  status = cusolverDnCreateSyevjInfo(&syevj_params);
  assert(CUSOLVER_STATUS_SUCCESS == status);
  
  /* default value of tolerance is machine zero */
  status = cusolverDnXsyevjSetTolerance(&syevj_params,tol);
  assert(CUSOLVER_STATUS_SUCCESS == status);
  
  /* default value of max. sweeps is 100 */
  status = cusolverDnXsyevjSetMaxSweeps(&syevj_params,max_sweeps);
  assert(CUSOLVER_STATUS_SUCCESS == status);

  /* step 3: copy A to device */
  cudaStat2 = cudaMalloc ((void**)&d_W, sizeof(cuDoubleComplex) * n);
  cudaStat3 = cudaMalloc ((void**)&d_info, sizeof(int));
  assert(cudaSuccess == cudaStat2);
  assert(cudaSuccess == cudaStat3);
  
  /* step 4: query working space of sygvj */
  status = cusolverDnZhegvj_bufferSize(cusolverH,itype,jobz,uplo,n,H,n,S,n,d_W,&lwork,syevj_params);
  assert(CUSOLVER_STATUS_SUCCESS == status);
  
  cudaStat1 = cudaMalloc((void**)&d_work, sizeof(cuDoubleComplex)*lwork);
  assert(cudaSuccess == cudaStat1);
  
  /* step 5: compute eigen-pair   */
  status = cusolverDnZhegvj(cusolverH,itype,jobz,uplo,n,H,n,S,n,d_W,d_work,lwork,d_info,syevj_params);
  cudaStat1 = cudaDeviceSynchronize();
  assert(CUSOLVER_STATUS_SUCCESS == status);
  assert(cudaSuccess == cudaStat1);

  cudaStat1 = cudaMemcpy(eig, d_W, sizeof(double)*ne, cudaMemcpyDeviceToHost);
  cudaStat2 = cudaMemcpy(z, H, sizeof(cuDoubleComplex)*n*ne, cudaMemcpyDeviceToHost);
  cudaStat3 = cudaMemcpy(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
  assert(cudaSuccess == cudaStat1);
  assert(cudaSuccess == cudaStat2);
  assert(cudaSuccess == cudaStat3);

  if ( 0 == info ){
    printf("sygvj converges \n");
  }else if ( 0 > info ){
    printf("%d-th parameter is wrong \n", -info);
    exit(1);
  }else{
    printf("WARNING: info = %d : sygvj does not converge \n", info );
  }

  status = cusolverDnXsyevjGetSweeps(cusolverH,syevj_params,&executed_sweeps);
  assert(CUSOLVER_STATUS_SUCCESS == status);

  status = cusolverDnXsyevjGetResidual(cusolverH,syevj_params,&residual);
  assert(CUSOLVER_STATUS_SUCCESS == status);

  printf("residual |A - V*W*V**H|_F = %E \n", residual );
  printf("number of executed sweeps = %d \n", executed_sweeps );

  /* free resources */
  if (d_W    ) cudaFree(d_W);
  if (d_info ) cudaFree(d_info);
  if (d_work ) cudaFree(d_work);

  if (cusolverH   ) cusolverDnDestroy(cusolverH);   
  if (stream      ) cudaStreamDestroy(stream);
  if (syevj_params) cusolverDnDestroySyevjInfo(syevj_params);

  //  cudaDeviceReset();

  return ;
}


void cusolver_real(double *H,double *S,int n,int ne,double tol,int max_sweeps,double* eig,double *z){

  cusolverDnHandle_t cusolverH = NULL;
  cudaStream_t stream = NULL;
  syevjInfo_t syevj_params = NULL;

  cusolverStatus_t status = CUSOLVER_STATUS_SUCCESS;
  cudaError_t cudaStat1 = cudaSuccess;
  cudaError_t cudaStat2 = cudaSuccess;
  cudaError_t cudaStat3 = cudaSuccess;

  double *d_W = NULL; /* eigenvalues on device*/
  int *d_info = NULL; /* error info */
  int  lwork = 0;     /* size of workspace */
  double *d_work = NULL; /* device workspace for syevj */
  int info = 0;       /* host copy of error info */

/* configuration of syevj  */
  const cusolverEigType_t itype = CUSOLVER_EIG_TYPE_1; //Solve H psi=S lambda psi
  const cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR; // compute eigenvectors.
  const cublasFillMode_t  uplo = CUBLAS_FILL_MODE_LOWER;

  /* numerical results of syevj  */
  double residual = 0;
  int executed_sweeps = 0;
  
  /* step 1: create cusolver handle, bind a stream  */
  status = cusolverDnCreate(&cusolverH);
  assert(CUSOLVER_STATUS_SUCCESS == status);

  cudaStat1 = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
  assert(cudaSuccess == cudaStat1);

  status = cusolverDnSetStream(cusolverH, stream);
  assert(CUSOLVER_STATUS_SUCCESS == status);

  /* step 2: configuration of syevj */
  status = cusolverDnCreateSyevjInfo(&syevj_params);
  assert(CUSOLVER_STATUS_SUCCESS == status);
  
  /* default value of tolerance is machine zero */
  status = cusolverDnXsyevjSetTolerance(&syevj_params,tol);
  assert(CUSOLVER_STATUS_SUCCESS == status);
  
  /* default value of max. sweeps is 100 */
  status = cusolverDnXsyevjSetMaxSweeps(&syevj_params,max_sweeps);
  assert(CUSOLVER_STATUS_SUCCESS == status);

  /* step 3: copy A to device */
  cudaStat2 = cudaMalloc ((void**)&d_W, sizeof(double) * n);
  cudaStat3 = cudaMalloc ((void**)&d_info, sizeof(int));
  assert(cudaSuccess == cudaStat2);
  assert(cudaSuccess == cudaStat3);
  
  /* step 4: query working space of sygvj */
  status = cusolverDnDsygvj_bufferSize(cusolverH,itype,jobz,uplo,n,H,n,S,n,d_W,&lwork,syevj_params);
  assert(CUSOLVER_STATUS_SUCCESS == status);
  
  cudaStat1 = cudaMalloc((void**)&d_work, sizeof(double)*lwork);
  assert(cudaSuccess == cudaStat1);
  
  /* step 5: compute eigen-pair   */
  status = cusolverDnDsygvj(cusolverH,itype,jobz,uplo,n,H,n,S,n,d_W,d_work,lwork,d_info,syevj_params);
  cudaStat1 = cudaDeviceSynchronize();
  assert(CUSOLVER_STATUS_SUCCESS == status);
  assert(cudaSuccess == cudaStat1);

  cudaStat1 = cudaMemcpy(eig, d_W, sizeof(double)*ne, cudaMemcpyDeviceToHost);
  cudaStat2 = cudaMemcpy(z, H, sizeof(double)*n*ne, cudaMemcpyDeviceToHost);
  cudaStat3 = cudaMemcpy(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
  assert(cudaSuccess == cudaStat1);
  assert(cudaSuccess == cudaStat2);
  assert(cudaSuccess == cudaStat3);

  if ( 0 == info ){
    printf("sygvj converges \n");
  }else if ( 0 > info ){
    printf("%d-th parameter is wrong \n", -info);
    exit(1);
  }else{
    printf("WARNING: info = %d : sygvj does not converge \n", info );
  }

  status = cusolverDnXsyevjGetSweeps(cusolverH,syevj_params,&executed_sweeps);
  assert(CUSOLVER_STATUS_SUCCESS == status);

  status = cusolverDnXsyevjGetResidual(cusolverH,syevj_params,&residual);
  assert(CUSOLVER_STATUS_SUCCESS == status);

  printf("residual |A - V*W*V**H|_F = %E \n", residual );
  printf("number of executed sweeps = %d \n", executed_sweeps );

  /* free resources */
  if (d_W    ) cudaFree(d_W);
  if (d_info ) cudaFree(d_info);
  if (d_work ) cudaFree(d_work);

  if (cusolverH   ) cusolverDnDestroy(cusolverH);   
  if (stream      ) cudaStreamDestroy(stream);
  if (syevj_params) cusolverDnDestroySyevjInfo(syevj_params);

  //  cudaDeviceReset();

  return ;
}
