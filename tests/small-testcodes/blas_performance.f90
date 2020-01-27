program test_blas
  INTEGER,PARAMETER::M=200,N=20000

  COMPLEX :: A(N,M),B(N,M),C(N,N)
  COMPLEX :: AA(M,M)
  REAL    :: R(N,M,4)
  REAL    :: time1,time2

  CALL cpu_time(time1)
  CALL random_number(r)
  A=cmplx(R(:,:,1),R(:,:,2))
  B=cmplx(R(:,:,3),R(:,:,4))
  C=0.0
  AA=cmplx(R(:M,:,1),R(:M,:,3))
  DO i=1,M
     DO j=1,i
        AA(i,j)=conjg(AA(j,i))
     ENDDO
     AA(i,i)=i*10
  ENDDO

  CALL cpu_time(time2)
  print *, "Init:",time2-time1

  CALL cpu_time(time1)
  CALL zgemm("N","T",N,N,M,cmplx(1.,0),A,N,B,N,cmplx(1.,0),C,N)
  CALL cpu_time(time2)
  print *, "zgemm:",time2-time1

  CALL cpu_time(time1)
  CALL zherk("U","N",N,M,cmplx(1.,0),A,N,cmplx(1.,0),C,N)
  CALL cpu_time(time2)
  print *, "zherk:",time2-time1

  CALL cpu_time(time1)
  CALL zher2k("U","N",N,M,cmplx(1.,0),A,N,B,N,cmplx(1.,0),C,N)
  CALL cpu_time(time2)
  print *, "zher2k:",time2-time1

  CALL cpu_time(time1)
  CALL zpotrf("U",M,AA,M,i)
  CALL cpu_time(time2)
  print *, "zpotrf:",time2-time1
  print *,i
  
end
  
