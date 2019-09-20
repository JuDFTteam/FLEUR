!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_wannier
  USE m_juDFT
  USE m_types_fleurinput_base
  IMPLICIT NONE
  PRIVATE
  !
  ! type for wannier-functions
  !
  TYPE,EXTENDS(t_fleurinput_base):: t_wann
    !New parameters not handled correctly yet...
    LOGICAL :: l_socmatvec
    LOGICAL :: l_socmatvecrs
    LOGICAL :: l_mmn0_unf_to_spn_unf
      LOGICAL :: l_mmn0_to_spn_unf
      LOGICAL :: l_mmn0_to_spn
      LOGICAL :: l_mmn0_to_spn2
      LOGICAL :: l_mmn0_unf_to_spn
      LOGICAL :: l_perpmag_unf_to_tor_unf
      LOGICAL :: l_perpmag_to_tor_unf
      LOGICAL :: l_perpmag_to_tor
      LOGICAL :: l_perpmag_unf_to_tor
      LOGICAL :: l_hsomtxvec_unf_to_lmpzsoc_unf
      LOGICAL :: l_hsomtxvec_to_lmpzsoc_unf
      LOGICAL :: l_hsomtxvec_to_lmpzsoc
      LOGICAL :: l_hsomtxvec_unf_to_lmpzsoc
      LOGICAL :: l_hsomtx_unf_to_hsoc_unf
      LOGICAL :: l_hsomtx_to_hsoc_unf
      LOGICAL :: l_hsomtx_to_hsoc
      LOGICAL :: l_hsomtx_unf_to_hsoc
      INTEGER :: perpmagl
      LOGICAL :: l_perpmagatlres

     INTEGER :: wan90version =3
     INTEGER :: oc_num_orbs =0
     INTEGER, ALLOCATABLE :: oc_orbs(:)
     LOGICAL :: l_unformatted =.FALSE.
     LOGICAL :: l_oc_f=.FALSE.
     LOGICAL :: l_ndegen=.FALSE.
     LOGICAL :: l_orbitalmom=.FALSE.
     LOGICAL :: l_orbcomp=.FALSE.
     LOGICAL :: l_orbcomprs=.FALSE.
     LOGICAL :: l_denmat=.FALSE.
     LOGICAL :: l_perturbrs=.FALSE.
     LOGICAL :: l_perturb=.FALSE.
     LOGICAL :: l_nedrho=.FALSE.
     LOGICAL :: l_anglmomrs=.FALSE.
     LOGICAL :: l_anglmom=.FALSE.
     LOGICAL :: l_spindisp=.FALSE.
     LOGICAL :: l_spindisprs=.FALSE.
     LOGICAL :: l_socspicom=.FALSE.
     LOGICAL :: l_socspicomrs=.FALSE.
     LOGICAL :: l_offdiposoprs=.FALSE.
     LOGICAL :: l_offdiposop=.FALSE.
     LOGICAL :: l_torque=.FALSE.
     LOGICAL :: l_torquers=.FALSE.
     LOGICAL :: l_atomlist=.FALSE.
     INTEGER :: atomlist_num=0
     INTEGER, ALLOCATABLE :: atomlist(:)
     LOGICAL :: l_berry=.FALSE.
     LOGICAL :: l_perpmagrs=.FALSE.
     LOGICAL :: l_perpmag=.FALSE.
     LOGICAL :: l_perpmagat=.FALSE.
     LOGICAL :: l_perpmagatrs=.FALSE.
     LOGICAL :: l_socmatrs=.FALSE.
     LOGICAL :: l_socmat=.FALSE.
     LOGICAL :: l_soctomom=.FALSE.
     LOGICAL :: l_kptsreduc2=.FALSE.
     LOGICAL :: l_nablapaulirs=.FALSE.
     LOGICAL :: l_nablars=.FALSE.
     LOGICAL :: l_surfcurr=.FALSE.
     LOGICAL :: l_updown=.FALSE.
     LOGICAL :: l_ahe=.FALSE.
     LOGICAL :: l_she=.FALSE.
     LOGICAL :: l_rmat=.FALSE.
     LOGICAL :: l_nabla=.FALSE.
     LOGICAL :: l_socodi=.FALSE.
     LOGICAL :: l_pauli=.FALSE.
     LOGICAL :: l_pauliat=.FALSE.
     LOGICAL :: l_potmat=.FALSE.
     LOGICAL :: l_projgen=.FALSE.
     LOGICAL :: l_plot_symm=.FALSE.
     LOGICAL :: l_socmmn0=.FALSE.
     LOGICAL :: l_bzsym=.FALSE.
     LOGICAL :: l_hopping=.FALSE.
     LOGICAL :: l_kptsreduc=.FALSE.
     LOGICAL :: l_prepwan90=.FALSE.
     LOGICAL :: l_plot_umdat=.FALSE.
     LOGICAL :: l_wann_plot=.FALSE.
     LOGICAL :: l_bynumber=.FALSE.
     LOGICAL :: l_stopopt=.FALSE.
     LOGICAL :: l_matrixmmn=.FALSE.
     LOGICAL :: l_matrixamn=.FALSE.
     LOGICAL :: l_projmethod=.FALSE.
     LOGICAL :: l_wannierize=.FALSE.
     LOGICAL :: l_plotw90=.FALSE.
     LOGICAL :: l_byindex=.FALSE.
     LOGICAL :: l_byenergy=.FALSE.
     LOGICAL :: l_proj_plot=.FALSE.
     LOGICAL :: l_bestproj=.FALSE.
     LOGICAL :: l_ikptstart=.FALSE.
     LOGICAL :: l_lapw=.FALSE.
     LOGICAL :: l_plot_lapw=.FALSE.
     LOGICAL :: l_fermi=.FALSE.
     LOGICAL :: l_dipole=.FALSE.
     LOGICAL :: l_dipole2=.FALSE.
     LOGICAL :: l_dipole3=.FALSE.
     LOGICAL :: l_mmn0=.FALSE.
     LOGICAL :: l_mmn0at=.FALSE.
     LOGICAL :: l_manyfiles=.FALSE.
     LOGICAL :: l_collectmanyfiles=.FALSE.
     LOGICAL :: l_ldauwan=.FALSE.
     LOGICAL :: l_lapw_kpts=.FALSE.
     LOGICAL :: l_lapw_gfleur=.FALSE.
     LOGICAL :: l_kpointgen=.FALSE.
     LOGICAL :: l_w90kpointgen=.FALSE.
     LOGICAL :: l_finishnocoplot=.FALSE.
     LOGICAL :: l_finishgwf=.FALSE.
     LOGICAL :: l_skipkov=.FALSE.
     LOGICAL :: l_matrixuHu=.FALSE.
     LOGICAL :: l_matrixuHu_dmi=.FALSE.
     INTEGER :: ikptstart=1
     INTEGER :: band_min(1:2)=-1
     INTEGER :: band_max(1:2)=-1
     INTEGER :: gfthick=0
     INTEGER :: gfcut=0
     INTEGER :: unigrid(6)=0
     INTEGER :: mhp(3)=0
     !---> gwf
     LOGICAL :: l_ms=.FALSE.
     LOGICAL :: l_sgwf=.FALSE.
     LOGICAL :: l_socgwf=.FALSE.
     LOGICAL :: l_gwf=.FALSE.
     LOGICAL :: l_bs_comf=.FALSE.
     LOGICAL :: l_exist=.FALSE.
     LOGICAL :: l_opened=.FALSE.
     LOGICAL :: l_cleverskip=.FALSE.
     LOGICAL :: l_dim(3)=.FALSE.
     REAL    :: scale_param=1.0
     REAL    :: aux_latt_const=8.0
     REAL    :: hdwf_t1=0.0
     REAL    :: hdwf_t2=0.0
     INTEGER :: nparampts=0
     CHARACTER(len=20) :: fn_eig=''
     CHARACTER(len=20) :: param_file='qpts'
     REAL, ALLOCATABLE :: param_vec(:, :)
     REAL, ALLOCATABLE :: param_alpha(:, :)
     CHARACTER(LEN=20), ALLOCATABLE :: jobList(:)
     !---> gwf
   CONTAINS
     PROCEDURE :: read_xml => read_xml_wannier
     PROCEDURE :: mpi_bc => mpi_bc_wannier
  END TYPE t_wann

  PUBLIC t_wann
CONTAINS

  SUBROUTINE mpi_bc_wannier(this,mpi_comm,irank)
    USE m_mpi_bc_tool
    CLASS(t_wann),INTENT(INOUT)::this
    INTEGER,INTENT(IN):: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL::irank
    INTEGER ::rank
    IF (PRESENT(irank)) THEN
       rank=irank
    ELSE
       rank=0
    END IF

    CALL mpi_bc(this%wan90version ,rank,mpi_comm)
    CALL mpi_bc(this%oc_num_orbs ,rank,mpi_comm)
    CALL mpi_bc(this%oc_orbs,rank,mpi_comm)
    CALL mpi_bc(this%l_unformatted ,rank,mpi_comm)
    CALL mpi_bc(this%l_oc_f,rank,mpi_comm)
    CALL mpi_bc(this%l_ndegen,rank,mpi_comm)
    CALL mpi_bc(this%l_orbitalmom,rank,mpi_comm)
    CALL mpi_bc(this%l_orbcomp,rank,mpi_comm)
    CALL mpi_bc(this%l_orbcomprs,rank,mpi_comm)
    CALL mpi_bc(this%l_denmat,rank,mpi_comm)
    CALL mpi_bc(this%l_perturbrs,rank,mpi_comm)
    CALL mpi_bc(this%l_perturb,rank,mpi_comm)
    CALL mpi_bc(this%l_nedrho,rank,mpi_comm)
    CALL mpi_bc(this%l_anglmomrs,rank,mpi_comm)
    CALL mpi_bc(this%l_anglmom,rank,mpi_comm)
    CALL mpi_bc(this%l_spindisp,rank,mpi_comm)
    CALL mpi_bc(this%l_spindisprs,rank,mpi_comm)
    CALL mpi_bc(this%l_socspicom,rank,mpi_comm)
    CALL mpi_bc(this%l_socspicomrs,rank,mpi_comm)
    CALL mpi_bc(this%l_offdiposoprs,rank,mpi_comm)
    CALL mpi_bc(this%l_offdiposop,rank,mpi_comm)
    CALL mpi_bc(this%l_torque,rank,mpi_comm)
    CALL mpi_bc(this%l_torquers,rank,mpi_comm)
    CALL mpi_bc(this%l_atomlist,rank,mpi_comm)
    CALL mpi_bc(this%atomlist_num,rank,mpi_comm)
    CALL mpi_bc(this%atomlist,rank,mpi_comm)
    CALL mpi_bc(this%l_berry,rank,mpi_comm)
    CALL mpi_bc(this%l_perpmagrs,rank,mpi_comm)
    CALL mpi_bc(this%l_perpmag,rank,mpi_comm)
    CALL mpi_bc(this%l_perpmagat,rank,mpi_comm)
    CALL mpi_bc(this%l_perpmagatrs,rank,mpi_comm)
    CALL mpi_bc(this%l_socmatrs,rank,mpi_comm)
    CALL mpi_bc(this%l_socmat,rank,mpi_comm)
    CALL mpi_bc(this%l_soctomom,rank,mpi_comm)
    CALL mpi_bc(this%l_kptsreduc2,rank,mpi_comm)
    CALL mpi_bc(this%l_nablapaulirs,rank,mpi_comm)
    CALL mpi_bc(this%l_nablars,rank,mpi_comm)
    CALL mpi_bc(this%l_surfcurr,rank,mpi_comm)
    CALL mpi_bc(this%l_updown,rank,mpi_comm)
    CALL mpi_bc(this%l_ahe,rank,mpi_comm)
    CALL mpi_bc(this%l_she,rank,mpi_comm)
    CALL mpi_bc(this%l_rmat,rank,mpi_comm)
    CALL mpi_bc(this%l_nabla,rank,mpi_comm)
    CALL mpi_bc(this%l_socodi,rank,mpi_comm)
    CALL mpi_bc(this%l_pauli,rank,mpi_comm)
    CALL mpi_bc(this%l_pauliat,rank,mpi_comm)
    CALL mpi_bc(this%l_potmat,rank,mpi_comm)
    CALL mpi_bc(this%l_projgen,rank,mpi_comm)
    CALL mpi_bc(this%l_plot_symm,rank,mpi_comm)
    CALL mpi_bc(this%l_socmmn0,rank,mpi_comm)
    CALL mpi_bc(this%l_bzsym,rank,mpi_comm)
    CALL mpi_bc(this%l_hopping,rank,mpi_comm)
    CALL mpi_bc(this%l_kptsreduc,rank,mpi_comm)
    CALL mpi_bc(this%l_prepwan90,rank,mpi_comm)
    CALL mpi_bc(this%l_plot_umdat,rank,mpi_comm)
    CALL mpi_bc(this%l_wann_plot,rank,mpi_comm)
    CALL mpi_bc(this%l_bynumber,rank,mpi_comm)
    CALL mpi_bc(this%l_stopopt,rank,mpi_comm)
    CALL mpi_bc(this%l_matrixmmn,rank,mpi_comm)
    CALL mpi_bc(this%l_matrixamn,rank,mpi_comm)
    CALL mpi_bc(this%l_projmethod,rank,mpi_comm)
    CALL mpi_bc(this%l_wannierize,rank,mpi_comm)
    CALL mpi_bc(this%l_plotw90,rank,mpi_comm)
    CALL mpi_bc(this%l_byindex,rank,mpi_comm)
    CALL mpi_bc(this%l_byenergy,rank,mpi_comm)
    CALL mpi_bc(this%l_proj_plot,rank,mpi_comm)
    CALL mpi_bc(this%l_bestproj,rank,mpi_comm)
    CALL mpi_bc(this%l_ikptstart,rank,mpi_comm)
    CALL mpi_bc(this%l_lapw,rank,mpi_comm)
    CALL mpi_bc(this%l_plot_lapw,rank,mpi_comm)
    CALL mpi_bc(this%l_fermi,rank,mpi_comm)
    CALL mpi_bc(this%l_dipole,rank,mpi_comm)
    CALL mpi_bc(this%l_dipole2,rank,mpi_comm)
    CALL mpi_bc(this%l_dipole3,rank,mpi_comm)
    CALL mpi_bc(this%l_mmn0,rank,mpi_comm)
    CALL mpi_bc(this%l_mmn0at,rank,mpi_comm)
    CALL mpi_bc(this%l_manyfiles,rank,mpi_comm)
    CALL mpi_bc(this%l_collectmanyfiles,rank,mpi_comm)
    CALL mpi_bc(this%l_ldauwan,rank,mpi_comm)
    CALL mpi_bc(this%l_lapw_kpts,rank,mpi_comm)
    CALL mpi_bc(this%l_lapw_gfleur,rank,mpi_comm)
    CALL mpi_bc(this%l_kpointgen,rank,mpi_comm)
    CALL mpi_bc(this%l_w90kpointgen,rank,mpi_comm)
    CALL mpi_bc(this%l_finishnocoplot,rank,mpi_comm)
    CALL mpi_bc(this%l_finishgwf,rank,mpi_comm)
    CALL mpi_bc(this%l_skipkov,rank,mpi_comm)
    CALL mpi_bc(this%l_matrixuHu,rank,mpi_comm)
    CALL mpi_bc(this%l_matrixuHu_dmi,rank,mpi_comm)
    CALL mpi_bc(this%ikptstart,rank,mpi_comm)
    CALL mpi_bc(this%band_min(1),rank,mpi_comm)
    CALL mpi_bc(this%band_max(1),rank,mpi_comm)
    CALL mpi_bc(this%band_min(2),rank,mpi_comm)
    CALL mpi_bc(this%band_max(2),rank,mpi_comm)
    CALL mpi_bc(this%gfthick,rank,mpi_comm)
    CALL mpi_bc(this%gfcut,rank,mpi_comm)
    CALL mpi_bc(this%unigrid(1),rank,mpi_comm)
    CALL mpi_bc(this%unigrid(2),rank,mpi_comm)
    CALL mpi_bc(this%unigrid(3),rank,mpi_comm)
    CALL mpi_bc(this%unigrid(4),rank,mpi_comm)
    CALL mpi_bc(this%unigrid(5),rank,mpi_comm)
    CALL mpi_bc(this%unigrid(6),rank,mpi_comm)
    CALL mpi_bc(this%mhp(1),rank,mpi_comm)
    CALL mpi_bc(this%mhp(2),rank,mpi_comm)
    CALL mpi_bc(this%mhp(3),rank,mpi_comm)
    CALL mpi_bc(this%l_ms,rank,mpi_comm)
    CALL mpi_bc(this%l_sgwf,rank,mpi_comm)
    CALL mpi_bc(this%l_socgwf,rank,mpi_comm)
    CALL mpi_bc(this%l_gwf,rank,mpi_comm)
    CALL mpi_bc(this%l_bs_comf,rank,mpi_comm)
    CALL mpi_bc(this%l_exist,rank,mpi_comm)
    CALL mpi_bc(this%l_opened,rank,mpi_comm)
    CALL mpi_bc(this%l_cleverskip,rank,mpi_comm)
    CALL mpi_bc(this%l_dim(1),rank,mpi_comm)
    CALL mpi_bc(this%l_dim(2),rank,mpi_comm)
    CALL mpi_bc(this%l_dim(3),rank,mpi_comm)
    CALL mpi_bc(this%scale_param,rank,mpi_comm)
    CALL mpi_bc(this%aux_latt_const,rank,mpi_comm)
    CALL mpi_bc(this%hdwf_t1,rank,mpi_comm)
    CALL mpi_bc(this%hdwf_t2,rank,mpi_comm)
    CALL mpi_bc(this%nparampts,rank,mpi_comm)
    CALL mpi_bc(this%param_vec,rank,mpi_comm)
    CALL mpi_bc(this%param_alpha,rank,mpi_comm)

    !Not done
    !CHARACTER(len=20) :: fn_eig=''
    !CHARACTER(len=20) :: param_file='qpts'
    !CHARACTER(LEN=20), ALLOCATABLE :: jobList(:)



  END SUBROUTINE mpi_bc_wannier

  SUBROUTINE read_xml_wannier(this,xml)
    USE m_types_xml
    CLASS(t_wann),INTENT(inout):: this
    TYPE(t_xml),INTENT(in)   :: xml
    ! Read in optional Wannier functions parameters

    CHARACTER(len=100):: xPathA
    CHARACTER(len=255):: valueString

    INTEGER           :: numberNodes,i,n,numtokens
    LOGICAL,ALLOCATABLE:: wannAtomList(:)

    xPathA = '/fleurInput/output/wannier'
    numberNodes = xml%getNumberOfNodes(xPathA)

    IF (numberNodes.EQ.1) THEN
       this%l_ms = evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@ms'))
       this%l_sgwf = evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@sgwf'))
       this%l_socgwf = evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@socgwf'))
       this%l_bs_comf = evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@bsComf'))
       this%l_atomlist = evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@atomList'))
    END IF

    xPathA = '/fleurInput/output/wannier/bandSelection'
    numberNodes = xml%getNumberOfNodes(xPathA)

    IF (numberNodes.EQ.1) THEN
       this%l_byindex=.TRUE.
       this%band_min(1) = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@minSpinUp'))
       this%band_max(1) = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@maxSpinUp'))
       xPathA = '/fleurInput/output/wannier/bandSelection/@minSpinDown'
       numberNodes = xml%getNumberOfNodes(xPathA)
       IF (numberNodes.EQ.1) THEN
          this%band_min(2) = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))))
       ELSE
          this%band_min(2) = this%band_min(1)
       END IF
       xPathA = '/fleurInput/output/wannier/bandSelection/@maxSpinDown'
       numberNodes = xml%getNumberOfNodes(xPathA)
       IF (numberNodes.EQ.1) THEN
          this%band_max(2) = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))))
       ELSE
          this%band_max(2) = this%band_max(1)
       END IF
       this%l_byindex = .TRUE.
    END IF

    xPathA = '/fleurInput/output/wannier/jobList'
    numberNodes = xml%getNumberOfNodes(xPathA)

    IF (numberNodes.EQ.1) THEN
       xPathA = '/fleurInput/output/wannier/jobList/text()'

       ! Note: At the moment only 255 characters for the text in this node. Maybe this is not enough.
       valueString = xml%getAttributeValue(TRIM(ADJUSTL(xPathA)))
       numTokens = xml%countStringTokens(valueString)
       ALLOCATE(this%jobList(numTokens))
       DO i = 1, numTokens
          this%jobList(i) = xml%popFirstStringToken(valueString)
       END DO
    END IF

    ALLOCATE(wannAtomList(xml%get_nat()))
    DO i=1,xml%get_nat()
       wannAtomList(i)= evaluateFirstBoolOnly(xml%getAttributeValue(xml%posPath(i)//'/@wannier'))
    ENDDO
    this%atomlist_num = COUNT(wannAtomList)
    n=0
    DO i=1,xml%get_nat()
       IF (wannAtomList(i)) THEN
          n=n+1
          this%atomlist(n) = i
       ENDIF
    ENDDO

    DEALLOCATE(wannAtomList)

  END SUBROUTINE read_xml_wannier
END MODULE m_types_wannier
