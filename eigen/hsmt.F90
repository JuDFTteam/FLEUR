MODULE m_hsmt
  USE m_juDFT
  IMPLICIT NONE
CONTAINS
  SUBROUTINE hsmt(atoms,sphhar,sym,enpara,&
       ispin,input,mpi,noco,cell,lapw,usdus,td,smat,hmat)
    USE m_hsmt_nonsph
    USE m_hsmt_sph
    use m_hsmt_lo
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
    INTEGER, INTENT (IN) :: ispin  
    
    !locals
    REAL, ALLOCATABLE    :: fj(:,:,:),gj(:,:,:)

    INTEGER :: iintsp,jintsp,n,i,ii
    LOGICAL :: l_socfirst
    COMPLEX :: chi(2,2),chi0(2,2),chi_one

    TYPE(t_lapwmat)::smat_tmp,hmat_tmp

    !
    l_socfirst= noco%l_soc .AND. noco%l_noco .AND. (.NOT. noco%l_ss)
    IF (noco%l_noco.AND..NOT.noco%l_ss) THEN
       CALL smat_tmp%alloc(smat(1,1)%l_real,smat(1,1)%matsize1,smat(1,1)%matsize2)
       CALL hmat_tmp%alloc(smat(1,1)%l_real,smat(1,1)%matsize1,smat(1,1)%matsize2)
    ENDIF
    
    ALLOCATE(fj(MAXVAL(lapw%nv),0:atoms%lmaxd,MERGE(2,1,noco%l_noco)))
    ALLOCATE(gj(MAXVAL(lapw%nv),0:atoms%lmaxd,MERGE(2,1,noco%l_noco)))

    iintsp=1;jintsp=1;chi_one=1.0 !Defaults in non-noco case
       DO n=1,atoms%ntype
          CALL timestart("fjgj coefficients")
          CALL hsmt_fjgj(input,atoms,cell,lapw,noco,usdus,n,ispin,fj,gj)
          CALL timestop("fjgj coefficients")
          IF (.NOT.noco%l_noco) THEN
             CALL timestart("spherical setup")
             CALL hsmt_sph(n,atoms,mpi,ispin,input,noco,cell,1,1,chi_one,lapw,&
                  enpara%el0,td%e_shift,usdus,fj,gj,smat(1,1),hmat(1,1))
             CALL timestop("spherical setup")
             CALL timestart("non-spherical setup")
             CALL hsmt_nonsph(n,mpi,sym,atoms,ispin,1,1,chi_one,noco,cell,lapw,td,fj,gj,hmat(1,1))
             CALL timestop("non-spherical setup")
             CALL hsmt_lo(input,atoms,sym,cell,mpi,noco,lapw,usdus,td,fj,gj,n,chi_one,ispin,iintsp,jintsp,hmat(1,1),smat(1,1))
          ELSEIF(noco%l_noco.AND..NOT.noco%l_ss) THEN
             CALL hsmt_spinor(ispin,n,noco,input,chi0,chi)
             CALL hmat_tmp%clear();CALL smat_tmp%clear()
             CALL hsmt_sph(n,atoms,mpi,ispin,input,noco,cell,1,1,chi_one,lapw,enpara%el0,td%e_shift,usdus,fj,gj,smat_tmp,hmat_tmp)
             CALL hsmt_nonsph(n,mpi,sym,atoms,ispin,1,1,chi_one,noco,cell,lapw,td,fj,gj,hmat_tmp)
             CALL hsmt_lo(input,atoms,sym,cell,mpi,noco,lapw,usdus,td,fj,gj,n,chi_one,ispin,iintsp,jintsp,hmat_tmp,smat_tmp)
             CALL hsmt_distspins(chi,smat_tmp,smat)
             CALL hsmt_distspins(chi,hmat_tmp,hmat)
          ELSE !l_ss
             CALL hsmt_spinor(ispin,n,noco,input,chi0,chi)
             DO iintsp=1,2
                DO jintsp=1,2
                   CALL hsmt_sph(n,atoms,mpi,ispin,input,noco,cell,iintsp,jintsp,chi(iintsp,jintsp),lapw,enpara%el0,td%e_shift,usdus,fj,gj,smat(iintsp,jintsp),hmat(iintsp,jintsp))
                   CALL hsmt_nonsph(n,mpi,sym,atoms,ispin,iintsp,jintsp,chi(iintsp,jintsp),noco,cell,lapw,td,fj,gj,hmat(iintsp,jintsp))
                ENDDO
             ENDDO
          ENDIF
          IF (noco%l_constr.OR.l_socfirst) THEN
             stop "offdiag"
             !CALL hsmt_offdiag()
          ENDIF
    END DO

    
    RETURN
  END SUBROUTINE hsmt
END MODULE m_hsmt
