program test_blas
!  use magma
  use omp_lib
  implicit none
  INTEGER,PARAMETER::Ns(3)=(/5000,10000,15000/)
  INTEGER,PARAMETER::M=1000
  COMPLEX,ALLOCATABLE :: A(:,:),B(:,:),z(:,:)
  REAL,ALLOCATABLE    :: R(:,:,:),eig(:)
  COMPLEX,ALLOCATABLE :: WORK(:)
  REAL,ALLOCATABLE    ::rwork(:)
  INTEGER,ALLOCATABLE ::iwork(:),ifail(:)
  INTEGER :: i,j,lwork,lrwork,liwork,ne,info,nn,n 

  REAL    :: time1,time2,wtemp(1)

  DO nn=1,size(Ns)
        N=Ns(nn)
  time1=omp_get_wtime()
  ALLOCATE (eig(N),ifail(N))
  Allocate(R(N,N,4),A(N,N),b(N,N),z(N,N))
  CALL random_number(r)
  A=cmplx(R(:,:,1),R(:,:,2))
  call zherk("U","N",N,N,1.0,cmplx(R(:,:,3),R(:,:,4)),N,0.0,B,N)
  DO i=1,N
     DO j=1,i
        A(i,j)=conjg(A(j,i))
        B(i,j)=conjg(B(j,i))
     ENDDO
     A(i,i)=real(A(i,i))
  ENDDO
  time2=omp_get_wtime()
  print *, "Init:",time2-time1

  time1=omp_get_wtime()
  allocate(iwork(5*N),rwork(7*N))
  CALL zhegvx (1, "V", "I", "U", N, A, N, B, N, 0.0, 0.0, 1, M, 1E-8, ne, eig, Z, N, Wtemp, -1, RWORK, IWORK, IFAIL, INFO)
  lwork=wtemp(1)
  print *,"lwork:",lwork
   allocate(work(lwork))
  CALL zhegvx (1, "V", "I", "U", N, A, N, B, N, 0.0, 0.0, 1, M, 1E-8, ne, eig, Z, N, WORK, LWORK, RWORK, IWORK, IFAIL, INFO)
  time2=omp_get_wtime()
  print *, N,"MKL:",time2-time1
  deallocate(work,iwork,rwork,A,B,z,r,eig,ifail)
  enddo   
 
 
end
  
