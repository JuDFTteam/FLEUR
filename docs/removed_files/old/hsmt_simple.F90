!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Gr�nberg Institut, Forschungszentrum J�lich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt_simple
  use m_juDFT
  implicit none
CONTAINS
  SUBROUTINE hsmt_simple(jspin,bkpt,input,sym,cell,atoms,lapw,td,noco,usdus,enpara,hmat,smat)
    use m_types
    use m_hsmt_fjgj
    USE m_hsmt_blas
    
    TYPE(t_input),INTENT(IN)      :: input
    TYPE(t_sym),INTENT(IN)        :: sym
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_enpara),INTENT(IN)     :: enpara
    TYPE(t_lapw),INTENT(IN)       :: lapw
    TYPE(t_usdus),INTENT(IN)      :: usdus
    TYPE(t_noco),INTENT(IN)       :: noco
    TYPE(t_tlmplm),INTENT(IN)     :: td
    TYPE(t_lapwmat),INTENT(INOUT) :: smat,hmat
    INTEGER,INTENT(IN)            :: jspin
    REAL,    INTENT (IN)          :: bkpt(3) 

    integer::k,jsp,n,nn,i
    REAL, ALLOCATABLE    :: fj(:,:,:,:),gj(:,:,:,:)
    REAL, ALLOCATABLE    :: gk(:,:,:),vk(:,:,:)
    REAL                 :: v(3),diff_h,diff_s
     REAL, PARAMETER      :: eps = 1.0e-30


    
    ALLOCATE(fj(lapw%dim_nbasfcn(),0:atoms%lmaxd,atoms%ntype,input%jspins))
    ALLOCATE(gj(lapw%dim_nbasfcn(),0:atoms%lmaxd,atoms%ntype,input%jspins))
    DO jsp=jspin,jspin

       !Set up the k+G+qss vectors
       ALLOCATE(vk(lapw%dim_nbasfcn(),3,1),gk(lapw%dim_nbasfcn(),3,1))

       DO k = 1,lapw%nv(jsp)
          v=bkpt+(/lapw%k1(k,jsp),lapw%k2(k,jsp),lapw%k3(k,jsp)/)!-noco%qss/2
          vk(k,:,1) = v
          gk(k,:,1) = MATMUL(TRANSPOSE(cell%bmat),v)/MAX (lapw%rk(k,jsp),eps)
       ENDDO
       

       CALL hsmt_fjgj(input,atoms,jsp,cell,lapw,usdus,fj,gj)
       
       CALL timestart("hsmt_blas")
       CALL  hsmt_blas(sym,atoms,jsp,noco,cell,lapw,td,usdus,gk,vk,fj,gj,smat,hmat)
       CALL timestop("hsmt_blas")
       
    ENDDO

    !CALL hsmt_lo()

  end SUBROUTINE hsmt_simple
end MODULE m_hsmt_simple
