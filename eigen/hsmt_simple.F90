!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt_simple
  use m_juDFT
  implicit none
CONTAINS
  SUBROUTINE hsmt_simple(jspin,jspins,bkpt,dimension,input,sym,cell,atoms,lapw,td,usdus,enpara,hamOvlp)
    use m_types
    use m_hsmt_fjgj
    use m_hsmt_overlap
    use m_hsmt_hamil
    TYPE(t_dimension),INTENT(IN)  :: dimension
    TYPE(t_input),INTENT(IN)      :: input
    TYPE(t_sym),INTENT(IN)        :: sym
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_enpara),INTENT(IN)     :: enpara
    TYPE(t_lapw),INTENT(IN)       :: lapw
    TYPE(t_usdus),INTENT(IN)      :: usdus
    TYPE(t_tlmplm),INTENT(IN)     :: td
    TYPE(t_hamovlp),INTENT(INOUT) :: hamovlp
    INTEGER,INTENT(IN)            :: jspin,jspins
    REAL,    INTENT (IN)          :: bkpt(3) 

    integer::k,jsp,n,nn,i
    REAL, ALLOCATABLE    :: fj(:,:,:,:),gj(:,:,:,:)
    REAL, ALLOCATABLE    :: gk(:,:,:),vk(:,:,:)
    REAL                 :: v(3)
     REAL, PARAMETER      :: eps = 1.0e-30


    hamovlp%l_real=.false.
    allocate(hamovlp%h_c(dimension%nbasfcn,dimension%nbasfcn))
    allocate(hamovlp%s_c(dimension%nbasfcn,dimension%nbasfcn))

    ALLOCATE(fj(dimension%nbasfcn,0:atoms%lmaxd,atoms%ntype,jspins))
    ALLOCATE(gj(dimension%nbasfcn,0:atoms%lmaxd,atoms%ntype,jspins))
    DO jsp=jspin,jspin

       !Set up the k+G+qss vectors
       ALLOCATE(vk(dimension%nbasfcn,3,1),gk(dimension%nbasfcn,3,1))

       DO k = 1,lapw%nv(jsp)
          v=bkpt+(/lapw%k1(k,jsp),lapw%k2(k,jsp),lapw%k3(k,jsp)/)!-noco%qss/2
          vk(k,:,1) = v
          gk(k,:,1) = MATMUL(TRANSPOSE(cell%bmat),v)/MAX (lapw%rk(k,jsp),eps)
       ENDDO
       

       CALL hsmt_fjgj(input,atoms,jsp,cell,lapw,usdus,fj,gj)
       
       CALL timestart("hsmt_overlap")
       CALL  hsmt_overlap(sym,atoms,jsp,cell,lapw,usdus,gk,vk,fj,gj,hamovlp)
       CALL timestop("hsmt_overlap")
       
       CALL timestart("hsmt_overlap")
       CALL  hsmt_hamil(sym,atoms,jsp,cell,lapw,td,gk,vk,fj,gj,hamovlp)
       CALL timestop("hsmt_overlap")
    enddo

    i=0
    DO n=1,lapw%nv(1)
       DO nn=1,n
          i=i+1
          write(61,"(2i5,4f10.7)") n,nn,hamovlp%b_c(i),hamovlp%s_c(nn,n)
          write(62,"(2i5,4f10.7)") n,nn,hamovlp%a_c(i),hamovlp%h_c(nn,n)
       ENDDO
    enddo
    stop "debug"
  end SUBROUTINE hsmt_simple
end MODULE m_hsmt_simple
