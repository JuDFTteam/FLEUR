!--------------------------------------------------------------------------------
! Copyright (c) 2023 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_matmul_dgemm
#ifdef _OPENACC
    use openacc    
    use cublas
#endif
#ifdef CPP_MAGMA
    use magma
#endif
    implicit none
    PRIVATE
    integer,PARAMETER:: blas_select=1,cublas_select=2,magma_select=3,magmablas_select=4

    interface blas_matmul
      module procedure:: blas_matmul_r,blas_matmul_c
    end interface

    public :: blas_matmul

    contains
    subroutine blas_matmul_r(m,n,k,a,b,c,alpha,beta,op_a,op_b)
        INTEGER,INTENT(IN):: n,m,k
        REAL,INTENT(IN)   :: a(:,:),b(:,:)
        REAL,INTENT(INOUT):: c(:,:)
        REAL,INTENT(IN),OPTIONAL:: alpha,beta
        CHARACTER,INTENT(IN),OPTIONAL :: op_a,op_b

        REAL      :: alphaa,betaa
        CHARACTER :: op_aa,op_bb
        INTEGER   :: lda,ldb,ldc

        alphaa=1.0; betaa=0.0
        if (present(alpha)) alphaa=alpha
        if (present(beta))  betaa=beta
      
        call priv_set_defaults(op_aa,op_bb,lda,ldb,ldc,m,n,k,a,b,c,op_a,op_b)
       
        select case (priv_select_multiply_r(a,b,c))   
        case (blas_select)
            call dgemm(op_aa,op_bb,m,n,k, alphaa, a, lda, b,ldb,betaa, c, ldc)
#ifdef _OPENACC
        case (cublas_select)
            if (IS_CONTIGUOUS(a).and.IS_CONTIGUOUS(b).and.IS_CONTIGUOUS(c)) THEN
                !$acc host_data use_device(a,b,c)
                call cublasDgemm(op_aa, op_bb, m, n, k, alphaa,a, lda, b, ldb, betaa, c, ldc)
                !$acc end host_data
            else
                block
                    REAL:: aa(size(a,1),size(a,2)),bb(size(b,1),size(b,2)),cc(size(c,1),size(c,2))
                    !$acc data create(aa,bb,cc)
                    !$acc kernels
                    aa=a
                    bb=b
                    cc=c
                    !$acc end kernels
                    !$acc host_data use_device(aa,bb,cc)
                    call cublasdgemm(op_aa, op_bb, m, n, k, alphaa,aa, lda, bb, ldb, betaa, cc, ldc)
                    !$acc end host_data
                    !$acc kernels
                    c=cc
                    !$acc end kernels
                    !$acc end data
                end block
            ENDIF  
#endif            
#ifdef CPP_MAGMA            
        case (magma_select)
            call magmaf_dgemm(op_aa, op_bb, m, n, k, alphaa,a, lda, b, ldb, betaa, c, ldc)
        case (magmablas_select)    
            call magmablasf_dgemm(op_aa, op_bb, m, n, k, alphaa,a, lda, b, ldb, betaa, c, ldc)
#endif        
        end select 

    end subroutine


    subroutine blas_matmul_c(m,n,k,a,b,c,alpha,beta,op_a,op_b)
        INTEGER,INTENT(IN):: n,m,k
        COMPLEX,INTENT(IN)   :: a(:,:),b(:,:)
        COMPLEX,INTENT(INOUT):: c(:,:)
        COMPLEX,INTENT(IN),OPTIONAL:: alpha,beta
        CHARACTER,INTENT(IN),OPTIONAL :: op_a,op_b

        COMPLEX      :: alphaa,betaa
        CHARACTER :: op_aa,op_bb
        INTEGER   :: lda,ldb,ldc

        alphaa=1.0; betaa=0.0
        if (present(alpha)) alphaa=alpha
        if (present(beta))  betaa=beta
      
        call priv_set_defaults(op_aa,op_bb,lda,ldb,ldc,m,n,k,a,b,c,op_a,op_b)
       
        select case (priv_select_multiply_c(a,b,c))   
        case (blas_select)
            call zgemm(op_aa,op_bb,m,n,k, alphaa, a, lda, b,ldb,betaa, c, ldc)
#ifdef _OPENACC
        case (cublas_select)
            if (IS_CONTIGUOUS(a).and.IS_CONTIGUOUS(b).and.IS_CONTIGUOUS(c)) THEN
                !$acc host_data use_device(a,b,c)
                call cublaszgemm(op_aa, op_bb, m, n, k, alphaa,a, lda, b, ldb, betaa, c, ldc)
                !$acc end host_data
            else
                block
                    complex:: aa(size(a,1),size(a,2)),bb(size(b,1),size(b,2)),cc(size(c,1),size(c,2))
                    !$acc data create(aa,bb,cc)
                    !$acc kernels
                    aa=a
                    bb=b
                    cc=c
                    !$acc end kernels
                    !$acc host_data use_device(aa,bb,cc)
                    call cublaszgemm(op_aa, op_bb, m, n, k, alphaa,aa, lda, bb, ldb, betaa, cc, ldc)
                    !$acc end host_data
                    !$acc kernels
                    c=cc
                    !$acc end kernels
                    !$acc end data
                end block
            ENDIF    
#endif            
#ifdef CPP_MAGMA            
        case (magma_select)
            call magmaf_zgemm(op_aa, op_bb, m, n, k, alphaa,a, lda, b, ldb, betaa, c, ldc)
        case (magmablas_select)    
            call magmablasf_zgemm(op_aa, op_bb, m, n, k, alphaa,a, lda, b, ldb, betaa, c, ldc)
#endif        
        end select 

    end subroutine


    subroutine priv_set_defaults(op_aa,op_bb,lda,ldb,ldc,m,n,k,a,b,c,op_a,op_b)
        INTEGER,INTENT(IN):: n,m,k
        CLASS(*),INTENT(IN)   :: a(:,:),b(:,:),c(:,:)
        CHARACTER,INTENT(IN),OPTIONAL :: op_a,op_b

        CHARACTER,INTENT(OUT) :: op_aa,op_bb
        INTEGER,INTENT(OUT)   :: lda,ldb,ldc

        op_aa='N'; op_bb='N'
        if (present(op_a)) op_aa=op_a
        if (present(op_b)) op_bb=op_b
        lda=size(a,1)
        ldb=size(b,1)
        ldc=size(c,1)

    END subroutine

    integer function priv_select_multiply_r(a,b,c)result(sel)
        REAL,INTENT(IN):: a(:,:),b(:,:),c(:,:)

#ifdef _OPENACC
        
        if (acc_is_present(a).and.acc_is_present(b).and.acc_is_present(c)) THEN
            !All data on GPU

            sel=cublas_select
#ifdef CPP_MAGMA
            sel=magmablas_select
#endif
        ENDIF
#endif        
        sel=blas_select
        return 
    end function                

    integer function priv_select_multiply_c(a,b,c)result(sel)
    COMPLEX,INTENT(IN):: a(:,:),b(:,:),c(:,:)
    sel=blas_select
#ifdef _OPENACC
    if (acc_is_present(a).and.acc_is_present(b).and.acc_is_present(c)) THEN
        !All data on GPU
        sel=cublas_select            
#ifdef CPP_MAGMA
        sel=magmablas_select
#endif
    ENDIF
#endif    
    
end function                

end module