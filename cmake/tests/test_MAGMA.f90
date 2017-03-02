      program test
      use magma
      IMPLICIT NONE
      integer :: mout,err,iwork(5)
      REAL    :: eigtemp(5),rwork(5)
      COMPLEX :: a(5,5),b(5,5),work(5)

      call magmaf_zhegvdx_2stage_m(1,1,MagmaVec,MagmaRangeI,MagmaLower,5,a,5,b,5,&
           0.0,0.0,1,2,mout,eigTemp,work,5,rwork,5,iwork,5,err)

      end
