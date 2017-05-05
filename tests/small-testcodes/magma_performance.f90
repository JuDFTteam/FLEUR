program test_blas
  use magma
  use omp_lib
  implicit none
  INTEGER,PARAMETER::Ns(3)=(/5000,10000,15000/)
  INTEGER,PARAMETER::M=1000
  COMPLEX,ALLOCATABLE :: A(:,:),B(:,:),z(:,:)
  COMPLEX,ALLOCATABLE :: AA(:,:),BB(:,:)
  REAL,ALLOCATABLE    :: R(:,:,:)
  COMPLEX,ALLOCATABLE :: WORK(:)
  REAL,ALLOCATABLE    ::rwork(:),eig(:)
  INTEGER,ALLOCATABLE ::iwork(:)
  INTEGER :: i,j,lwork,lrwork,liwork,ne(5),info 
  INTEGER :: nn,n

  REAL    :: time1,time2

  call magmaf_init()
  DO nn=1,3
        N=Ns(nn)
   time1=omp_get_wtime()
  Allocate(eig(N),R(N,N,4),A(N,N),b(N,N),z(N,N))
  CALL random_number(r)
  A=cmplx(R(:,:,1),R(:,:,2))
  call zherk("L","N",N,N,1.0,cmplx(R(:,:,3),R(:,:,4)),N,0.0,B,N)
  DO i=1,N
     DO j=1,i
        A(i,j)=conjg(A(j,i))
 !       B(i,j)=conjg(B(j,i))
     ENDDO
     b(i,i)=b(i,i)+50.0
     A(i,i)=real(A(i,i))
  ENDDO
  time2=omp_get_wtime()
  print *, "Init:",time2-time1

 time1=omp_get_wtime()
 allocate(work(1),rwork(1),iwork(1))
 call magmaf_zhegvdx_m(4,1,'v','i','L',N,A,N,B,N,0.0,0.0,1,M,ne,eig,work,-1,rwork,lrwork,iwork,liwork,info)
 print*,info
 lwork=work(1);lrwork=rwork(1);liwork=iwork(1)
 print *,info,lwork
 deallocate(work,iwork,rwork)
 allocate(work(lwork),rwork(lrwork),iwork(liwork))
 call magmaf_zhegvdx_m(4,1,'v','i','L',N,A,N,B,N,0.0,0.0,1,M,ne,eig,work,lwork,rwork,lrwork,iwork,liwork,info)
 print *,ne,info 
 time2=omp_get_wtime()
 print *, N,"MAGMA:",time2-time1
 deallocate(work,iwork,rwork)
 deallocate(eig,R,A,b,z)

 enddo
 
 
end
  
