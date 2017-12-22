MODULE m_hsmt
  USE m_juDFT
  IMPLICIT NONE
CONTAINS
  SUBROUTINE hsmt(atoms,sphhar,sym,enpara,&
       jsp,input,mpi,noco,cell,lapw,usdus,td,smat,hmat)
    USE m_hsmt_nonsph
    USE m_hsmt_sph
    USE m_hsmt_distspins
    USE m_types
    USE m_hsmt_fjgj
    use m_hsmt_spinor
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)        :: mpi
    TYPE(t_input),INTENT(IN)      :: input
    TYPE(t_noco),INTENT(IN)       :: noco
    TYPE(t_sym),INTENT(IN)        :: sym
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_sphhar),INTENT(IN)     :: sphhar
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_enpara),INTENT(IN)     :: enpara
    TYPE(t_lapw),INTENT(IN)       :: lapw 
    TYPE(t_tlmplm),INTENT(IN)     :: td
    TYPE(t_usdus),INTENT(IN)      :: usdus
    TYPE(t_lapwmat),INTENT(INOUT) :: smat(:,:),hmat(:,:)
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: jsp  
    
    !locals
    REAL, ALLOCATABLE    :: fj(:,:,:),gj(:,:,:)

    INTEGER :: iintsp,jintsp,ispin,n
    LOGICAL :: l_socfirst
    COMPLEX :: chi(2,2),chi0(2,2),chi_one

    TYPE(t_lapwmat)::smat_tmp,hmat_tmp

    !
    l_socfirst= noco%l_soc .AND. noco%l_noco .AND. (.NOT. noco%l_ss)
    IF (noco%l_noco.AND..NOT.noco%l_ss) THEN
       CALL smat_tmp%alloc(smat(1,1)%l_real,smat(1,1)%matsize1,smat(1,1)%matsize2)
       CALL hmat_tmp%alloc(smat(1,1)%l_real,smat(1,1)%matsize1,smat(1,1)%matsize2)
    ENDIF
    
    ALLOCATE(fj(MAXVAL(lapw%nv),0:atoms%lmaxd,MERGE(2,1,noco%l_ss)))
    ALLOCATE(gj(MAXVAL(lapw%nv),0:atoms%lmaxd,MERGE(2,1,noco%l_ss)))

    iintsp=1;jintsp=1;chi_one=1.0 !Defaults in non-noco case
    DO ispin=MERGE(1,jsp,noco%l_noco),MERGE(2,jsp,noco%l_noco) !spin-loop over mt-spin
       DO n=1,atoms%ntype
          CALL hsmt_fjgj(input,atoms,cell,lapw,noco,usdus,n,ispin,fj,gj)
          IF (.NOT.noco%l_noco) THEN
             CALL hsmt_sph(n,atoms,mpi,ispin,input,noco,cell,1,1,chi_one,lapw,&
                  enpara%el0,td%e_shift,usdus,fj,gj,smat(1,1),hmat(1,1))
             CALL hsmt_nonsph(n,mpi,sym,atoms,ispin,1,1,chi_one,noco,cell,lapw,td,fj,gj,hmat(1,1))
             !CALL hsmt_lo()
          ELSEIF(noco%l_noco.AND..NOT.noco%l_ss) THEN
             CALL hsmt_spinor(ispin,n,noco,input,chi0,chi)
             CALL hmat_tmp%clear();CALL smat_tmp%clear()
             CALL hsmt_sph(n,atoms,mpi,ispin,input,noco,cell,1,1,chi_one,lapw,enpara%el0,td%e_shift,usdus,fj,gj,smat_tmp,hmat_tmp)
             CALL hsmt_nonsph(n,mpi,sym,atoms,ispin,1,1,chi_one,noco,cell,lapw,td,fj,gj,hmat_tmp)
             !CALL hsmt_lo()
             CALL hsmt_distspins(chi,smat_tmp,smat)
             CALL hsmt_distspins(chi,hmat_tmp,hmat)
          ELSE !l_ss
             DO iintsp=1,2
                DO jintsp=1,2
                   CALL hsmt_sph(n,atoms,mpi,ispin,input,noco,cell,iintsp,jintsp,chi(iintsp,jintsp),lapw,enpara%el0,td%e_shift,usdus,fj,gj,smat(iintsp,jintsp),hmat(iintsp,jintsp))
                   CALL hsmt_nonsph(n,mpi,sym,atoms,ispin,iintsp,jintsp,chi(iintsp,jintsp),noco,cell,lapw,td,fj,gj,hmat(1,1))
                ENDDO
             ENDDO
          ENDIF
          IF (noco%l_constr.OR.l_socfirst) THEN
             stop "offdiag"
             !CALL hsmt_offdiag()
          ENDIF
       END DO
    END DO

    
    RETURN
  END SUBROUTINE hsmt
END MODULE m_hsmt
