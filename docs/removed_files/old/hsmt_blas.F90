!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt_blas
  use m_juDFT
  implicit none
CONTAINS
  SUBROUTINE hsmt_blas(sym,atoms,isp,noco,cell,lapw,td,ud,gk,vk,fj,gj,smat,hmat)
!Calculate overlap matrix
    USE m_hsmt_ab
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_ylm
    IMPLICIT NONE
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_tlmplm),INTENT(IN)   :: td
    TYPE(t_usdus),INTENT(IN)    :: ud
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: isp
    !     ..
    !     .. Array Arguments ..
    REAL,INTENT(IN) :: gk(:,:,:),vk(:,:,:)
    REAL,INTENT(IN) :: fj(:,0:,:,:),gj(:,0:,:,:)
    TYPE(t_lapwmat),INTENT(INOUT)::smat,hmat

    
    INTEGER:: n,nn,na,aboffset,l,ll,m
    COMPLEX,ALLOCATABLE:: ab(:,:),tmpdata(:,:),tmp_s(:,:),tmp_h(:,:),ab1(:,:)


    ALLOCATE(ab(lapw%nv(isp),2*atoms%lmaxd*(atoms%lmaxd+2)+2),ab1(lapw%nv(isp),2*atoms%lmaxd*(atoms%lmaxd+2)+2))
    
    ALLOCATE(tmp_s(smat%matsize1,smat%matsize2),tmp_h(smat%matsize1,smat%matsize2))
    tmp_s=0.0;tmp_h=0.0;ab=0.0;ab1=0.0

    ntyploop: DO n=1,atoms%ntype
       DO nn = 1,atoms%neq(n)
          na = SUM(atoms%neq(:n-1))+nn
          IF ((sym%invsat(na)==0) .OR. (sym%invsat(na)==1)) THEN
             
             !--->   Calculate Overlapp matrix
             CALL timestart("ab-coefficients")
             CALL hsmt_ab(sym,atoms,isp,n,na,cell,lapw,gk,vk,fj,gj,ab,aboffset)
             CALL timestop("ab-coefficients")
             CALL timestart("Overlapp")
             CALL ZHERK("U","N",lapw%nv(isp),aboffset,1.,ab,SIZE(ab,1),1.0,tmp_s,SIZE(tmp_s,1))
             DO l=0,atoms%lmax(n)
                ll=l*(l+1)
                DO m=-l,l
                   ab1(:,1+ll+m)=SQRT(ud%ddn(l,n,isp))*ab(:,aboffset+1+ll+m)
                ENDDO
             ENDDO
             CALL ZHERK("U","N",lapw%nv(isp),aboffset,1.,ab1,SIZE(ab,1),1.0,tmp_s,SIZE(tmp_s,1))
             CALL timestop("Overlapp")
             CALL timestart("Hamiltonian")
         
             !Calculate Hamiltonian
             CALL zgemm("N","N",SIZE(ab,1),2*aboffset,2*aboffset,CMPLX(1.0,0.0),ab,SIZE(ab,1),td%h_loc(:,:,n,isp),SIZE(td%h_loc,1),CMPLX(0.,0.),ab1,SIZE(ab,1))
!             CALL zgemm("N","C",lapw%nv(isp),lapw%nv(isp),2*aboffset,CMPLX(1.0,0.0),ab,SIZE(ab,1),ab1,SIZE(ab,1),CMPLX(1.0,0),tmp_h,SIZE(tmp_h,1))
             CALL ZHERK("U","N",lapw%nv(isp),2*aboffset,1.,ab1,SIZE(ab,1),1.0,tmp_h,SIZE(tmp_h,1))
             CALL timestop("Hamiltonian")
          ENDIF
       END DO
    END DO ntyploop
    !Copy tmp array back
    IF (smat%l_real) THEN
       smat%data_r=smat%data_r+tmp_s
       hmat%data_r=hmat%data_r+tmp_h-td%e_shift*tmp_s
    ELSE
       smat%data_c=smat%data_c+tmp_s
       hmat%data_c=hmat%data_c+tmp_h-td%e_shift*tmp_s
    ENDIF
  END SUBROUTINE hsmt_blas

#if 1==2
    !this version uses zherk for Hamiltonian
    ntyploop: DO n=1,atoms%ntype
       DO nn = 1,atoms%neq(n)
          na = SUM(atoms%neq(:n-1))+nn
          IF ((sym%invsat(na)==0) .OR. (sym%invsat(na)==1)) THEN
             
             !--->   Calculate Overlapp matrix
             CALL timestart("ab-coefficients")
             CALL hsmt_ab(sym,atoms,isp,n,na,cell,lapw,gk,vk,fj,gj,ab,aboffset)
             CALL timestop("ab-coefficients")
             CALL timestart("Overlapp")
             CALL ZHERK("U","N",lapw%nv(isp),aboffset,1.,ab,SIZE(ab,1),1.0,tmp_s,SIZE(tmp_s,1))
             DO l=0,atoms%lmax(n)
                ll=l*(l+1)
                DO m=-l,l
                   ab1(:,1+ll+m)=SQRT(ud%ddn(l,n,isp))*ab(:,aboffset+1+ll+m)
                ENDDO
             ENDDO
             CALL ZHERK("U","N",lapw%nv(isp),aboffset,1.,ab1,SIZE(ab,1),1.0,tmp_s,SIZE(tmp_s,1))
             CALL timestop("Overlapp")
             CALL timestart("Hamiltonian")
         
             !Calculate Hamiltonian
             CALL zgemm("N","N",SIZE(ab,1),2*aboffset,2*aboffset,CMPLX(1.0,0.0),ab,SIZE(ab,1),td%h_loc(:,:,n,isp),SIZE(td%h_loc,1),CMPLX(0.,0.),ab1,SIZE(ab,1))
             CALL ZHERK("U","N",lapw%nv(isp),2*aboffset,1.,ab1,SIZE(ab,1),1.0,tmp_h,SIZE(tmp_h,1))
             CALL timestop("Hamiltonian")
          ENDIF
       END DO
    END DO ntyploop
    !Copy tmp array back
    IF (smat%l_real) THEN
       smat%data_r=smat%data_r+tmp_s
       hmat%data_r=hmat%data_r+tmp_h-td%e_shift*tmp_s
    ELSE
       smat%data_c=smat%data_c+tmp_s
       hmat%data_c=hmat%data_c+tmp_h-td%e_shift*tmp_s
    ENDIF
#endif
  
END MODULE m_hsmt_blas
