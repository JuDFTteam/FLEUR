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
    LOGICAL :: l_socmatvec=.FALSE.
    INTEGER :: socmatvecfmt=1
    LOGICAL :: l_socmatvecrs=.FALSE.
    INTEGER :: socmatvecrsfmt=1
    LOGICAL :: l_mmn0_unf_to_spn_unf=.FALSE.
    LOGICAL :: l_mmn0_to_spn_unf=.FALSE.
    LOGICAL :: l_mmn0_to_spn=.FALSE.
    LOGICAL :: l_mmn0_to_spn2=.FALSE.
    LOGICAL :: l_mmn0_unf_to_spn=.FALSE.
    LOGICAL :: l_perpmag_unf_to_tor_unf=.FALSE.
    LOGICAL :: l_perpmag_to_tor_unf=.FALSE.
    LOGICAL :: l_perpmag_to_tor=.FALSE.
    LOGICAL :: l_perpmag_unf_to_tor=.FALSE.
    LOGICAL :: l_hsomtxvec_unf_to_lmpzsoc_unf=.FALSE.
    LOGICAL :: l_hsomtxvec_to_lmpzsoc_unf=.FALSE.
    LOGICAL :: l_hsomtxvec_to_lmpzsoc=.FALSE.
    LOGICAL :: l_hsomtxvec_unf_to_lmpzsoc=.FALSE.
    LOGICAL :: l_hsomtx_unf_to_hsoc_unf=.FALSE.
    LOGICAL :: l_hsomtx_to_hsoc_unf=.FALSE.
    LOGICAL :: l_hsomtx_to_hsoc=.FALSE.
    LOGICAL :: l_hsomtx_unf_to_hsoc=.FALSE.
    INTEGER :: perpmagl
    LOGICAL :: l_perpmagatlres=.FALSE.

    INTEGER :: wan90version =31
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
     INTEGER :: anglmomrsfmt=1
     LOGICAL :: l_anglmom=.FALSE.
     INTEGER :: anglmomfmt=1
     LOGICAL :: l_spindisp=.FALSE.
     LOGICAL :: l_spindisprs=.FALSE.
     LOGICAL :: l_socspicom=.FALSE.
     LOGICAL :: l_socspicomrs=.FALSE.
     LOGICAL :: l_offdiposoprs=.FALSE.
     LOGICAL :: l_offdiposop=.FALSE.
     LOGICAL :: l_torque=.FALSE.
     INTEGER :: torquefmt=1
     LOGICAL :: l_torquers=.FALSE.
     INTEGER :: torquersfmt=1
     LOGICAL :: l_atomlist=.FALSE.
     INTEGER :: atomlist_num=0
     INTEGER, ALLOCATABLE :: atomlist(:)
     LOGICAL :: l_berry=.FALSE.
     LOGICAL :: l_perpmagrs=.FALSE.
     INTEGER :: perpmagrsfmt=1
     LOGICAL :: l_perpmag=.FALSE.
     INTEGER :: perpmagfmt=1
     LOGICAL :: l_perpmagat=.FALSE.
     INTEGER :: perpmagatfmt=1
     LOGICAL :: l_perpmagatrs=.FALSE.
     INTEGER :: perpmagatrsfmt=1
     LOGICAL :: l_socmatrs=.FALSE.
     INTEGER :: socmatrsfmt=1
     LOGICAL :: l_socmat=.FALSE.
     INTEGER :: socmatfmt=1
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
     INTEGER :: paulifmt=1
     LOGICAL :: l_pauliat=.FALSE.
     INTEGER :: pauliatfmt=1
     LOGICAL :: l_potmat=.FALSE.
     LOGICAL :: l_projgen=.FALSE.
     LOGICAL :: l_plot_symm=.FALSE.
     LOGICAL :: l_socmmn0=.FALSE.
     LOGICAL :: l_bzsym=.FALSE.
     LOGICAL :: l_hopping=.FALSE.
     INTEGER :: hoppingfmt=1
     LOGICAL :: l_kptsreduc=.FALSE.
     LOGICAL :: l_prepwan90=.FALSE.
     LOGICAL :: l_plot_umdat=.FALSE.
     LOGICAL :: l_wann_plot=.FALSE.
     LOGICAL :: l_bynumber=.FALSE.
     LOGICAL :: l_stopopt=.FALSE.
     LOGICAL :: l_stopuhu=.FALSE.
     LOGICAL :: l_stopupdown=.FALSE.
     LOGICAL :: l_matrixmmn=.FALSE.
     INTEGER :: matrixmmnfmt=1
     LOGICAL :: l_matrixamn=.FALSE.
     INTEGER :: matrixamnfmt=1
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
     INTEGER :: mmn0fmt=1
     LOGICAL :: l_mmn0at=.FALSE.
     INTEGER :: mmn0atfmt=1
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
     INTEGER :: matrixuHufmt=1
     LOGICAL :: l_matrixuHu_dmi=.FALSE.
     INTEGER :: matrixuHudmifmt=1
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
    CALL mpi_bc(this%socmatvecfmt,rank,mpi_comm)
    CALL mpi_bc(this%socmatvecrsfmt,rank,mpi_comm)
    CALL mpi_bc(this%l_socmatvec,rank,mpi_comm)
    CALL mpi_bc(this%l_socmatvecrs,rank,mpi_comm)    
    CALL mpi_bc(this%anglmomrsfmt,rank,mpi_comm)
    CALL mpi_bc(this%anglmomfmt,rank,mpi_comm)
    CALL mpi_bc(this%torquefmt,rank,mpi_comm)
    CALL mpi_bc(this%torquersfmt,rank,mpi_comm)
    CALL mpi_bc(this%perpmagrsfmt,rank,mpi_comm)
    CALL mpi_bc(this%perpmagfmt,rank,mpi_comm)
    CALL mpi_bc(this%perpmagatfmt,rank,mpi_comm)
    CALL mpi_bc(this%perpmagatrsfmt,rank,mpi_comm)
    CALL mpi_bc(this%socmatrsfmt,rank,mpi_comm)
    CALL mpi_bc(this%socmatfmt,rank,mpi_comm)
    CALL mpi_bc(this%paulifmt,rank,mpi_comm)
    CALL mpi_bc(this%pauliatfmt,rank,mpi_comm)
    CALL mpi_bc(this%hoppingfmt,rank,mpi_comm)
    CALL mpi_bc(this%matrixmmnfmt,rank,mpi_comm)
    CALL mpi_bc(this%matrixamnfmt,rank,mpi_comm)
    CALL mpi_bc(this%mmn0fmt,rank,mpi_comm)
    CALL mpi_bc(this%mmn0atfmt,rank,mpi_comm)
    CALL mpi_bc(this%matrixuHufmt,rank,mpi_comm)
    CALL mpi_bc(this%matrixuHudmifmt,rank,mpi_comm)
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
    CALL mpi_bc(this%l_perpmag_unf_to_tor_unf,rank,mpi_comm)
    CALL mpi_bc(this%l_perpmag_to_tor_unf,rank,mpi_comm)
    CALL mpi_bc(this%l_perpmag_to_tor,rank,mpi_comm)
    CALL mpi_bc(this%l_perpmag_unf_to_tor,rank,mpi_comm)
    CALL mpi_bc(this%l_hsomtxvec_unf_to_lmpzsoc_unf,rank,mpi_comm)
    CALL mpi_bc(this%l_hsomtxvec_to_lmpzsoc_unf,rank,mpi_comm)
    CALL mpi_bc(this%l_hsomtxvec_to_lmpzsoc,rank,mpi_comm)
    CALL mpi_bc(this%l_hsomtxvec_unf_to_lmpzsoc,rank,mpi_comm)
    CALL mpi_bc(this%l_hsomtx_unf_to_hsoc_unf,rank,mpi_comm)
    CALL mpi_bc(this%l_hsomtx_to_hsoc_unf,rank,mpi_comm)
    CALL mpi_bc(this%l_hsomtx_to_hsoc,rank,mpi_comm)
    CALL mpi_bc(this%l_hsomtx_unf_to_hsoc,rank,mpi_comm)   
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
    CALL mpi_bc(this%l_stopupdown,rank,mpi_comm)
    CALL mpi_bc(this%l_stopuhu,rank,mpi_comm)
    CALL mpi_bc(this%l_matrixmmn,rank,mpi_comm)
    CALL mpi_bc(this%l_matrixamn,rank,mpi_comm)
    CALL mpi_bc(this%l_projmethod,rank,mpi_comm)
    CALL mpi_bc(this%l_wannierize,rank,mpi_comm)
    CALL mpi_bc(this%l_plotw90,rank,mpi_comm)
    CALL mpi_bc(this%l_byindex,rank,mpi_comm)
    CALL mpi_bc(this%l_byenergy,rank,mpi_comm)
    CALL mpi_bc(this%l_mmn0_unf_to_spn_unf,rank,mpi_comm)
    CALL mpi_bc(this%l_mmn0_to_spn_unf,rank,mpi_comm)
    CALL mpi_bc(this%l_mmn0_to_spn,rank,mpi_comm)
    CALL mpi_bc(this%l_mmn0_to_spn2,rank,mpi_comm)
    CALL mpi_bc(this%l_mmn0_unf_to_spn,rank,mpi_comm)
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
    USE m_constants
    CLASS(t_wann),INTENT(inout):: this
    TYPE(t_xml),INTENT(INOUT) ::xml
    ! Read in optional Wannier functions parameters

    CHARACTER(len=100):: xPathA
    CHARACTER(len=255):: valueString
    CHARACTER(len=30):: jobname
    CHARACTER(len=30):: param
    INTEGER           :: parampos
    INTEGER           :: stat
    INTEGER           :: numberNodes,i,n,numtokens
    LOGICAL,ALLOCATABLE:: wannAtomList(:)
    LOGICAL :: l_param
    REAL :: version_real    

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
          IF(this%jobList(i)(1:1).EQ.'!')cycle
          parampos=index(this%jobList(i),'=')
          if(parampos.gt.1)then
             jobname=this%jobList(i)(1:parampos-1)
             param=this%jobList(i)(parampos+1:)
             l_param=.true.
          else   
             jobname=this%jobList(i)
             l_param=.false.
          endif   
          IF(this%jobList(i).EQ.'socmat')THEN
             this%l_socmat=.TRUE.
          ELSEIF(jobname.EQ.'socmatvec')THEN
             this%l_socmatvec=.TRUE.
             if(l_param)then
             read(param,*,iostat=stat) this%socmatvecfmt
             if(stat/=0)then
            CALL juDFT_error("problem with jobparam=",calledby="wann_read_inp")
             endif
            endif
          ELSEIF(jobname.EQ.'socmatvecrs')THEN
             this%l_socmatvecrs=.TRUE.
             if(l_param)then
             read(param,*,iostat=stat) this%socmatvecrsfmt
             if(stat/=0)then
            CALL juDFT_error("problem with jobparam=",calledby="wann_read_inp")
             endif
            endif
          ELSEIF(this%jobList(i).EQ.'unformatted')THEN
             this%l_unformatted=.TRUE.         
          ELSEIF(this%jobList(i).EQ.'ndegen')THEN
             this%l_ndegen=.TRUE.              
          ELSEIF(this%jobList(i).EQ.'socmatrs')THEN
             this%l_socmatrs=.TRUE.
          ELSEIF(this%jobList(i).EQ.'soctomom')THEN
             this%l_soctomom=.TRUE.
          ELSEIF(this%jobList(i).EQ.'surfcurr')THEN
             this%l_surfcurr=.TRUE.
          ELSEIF(this%jobList(i).EQ.'lapw_kpts')THEN
             this%l_lapw_kpts=.TRUE.
          ELSEIF(this%jobList(i).EQ.'updown')THEN
             this%l_updown=.TRUE.
          ELSEIF(this%jobList(i).EQ.'unformatted')THEN
             this%l_unformatted=.TRUE.            
          ELSEIF(this%jobList(i).EQ.'stopopt')THEN
             this%l_stopopt=.TRUE.
          ELSEIF(this%jobList(i).EQ.'stopuhu')THEN
             this%l_stopuhu=.TRUE.
          ELSEIF(this%jobList(i).EQ.'stopupdown')THEN
             this%l_stopupdown=.TRUE.
          ELSEIF(this%jobList(i).EQ.'projgen')THEN
             this%l_projgen=.TRUE.
          ELSEIF(this%jobList(i).EQ.'kpointgen')THEN
             this%l_kpointgen=.TRUE.
          ELSEIF(this%jobList(i).EQ.'potmat')THEN
             this%l_potmat=.TRUE.
          ELSEIF(this%jobList(i).EQ.'w90kpointgen')THEN
             this%l_w90kpointgen=.TRUE.
        !Not done
        !  ELSEIF(this%jobList(i).EQ.'lapw_gfleur')THEN
        !     this%l_lapw_gfleur=.TRUE.
        !     backspace(916)
        !     read(916,*,iostat=ios)task,this%gfthick,this%gfcut
        !     if (ios /= 0) CALL juDFT_error ("error reading gfcut", calledby="wann_read_inp")
        !     if(l_p0)write(oUnit,*)"gfcut=",this%gfthick,this%gfcut
        !Not done
        !  ELSEIF(this%jobList(i).EQ.'lapw')THEN
        !     this%l_lapw=.TRUE.
        !     backspace(916)
        !     read(916,*,iostat=ios)task,this%unigrid(:)
        !     if (ios /= 0) CALL juDFT_error ("error reading unigrid", calledby="wann_read_inp")
        !     if(l_p0)write(oUnit,*)"unigrid=",this%unigrid(:)
          ELSEIF(this%jobList(i).EQ.'plot_lapw')THEN
             this%l_plot_lapw=.TRUE.
          ELSEIF(this%jobList(i).EQ.'bzsym')THEN
             this%l_bzsym=.TRUE.
             !this%l_kpts_fullbz=.false.
          ELSEIF(jobname.EQ.'mmn0')THEN
             this%l_mmn0=.TRUE.
             if(l_param)then
             read(param,*,iostat=stat) this%mmn0fmt
             if(stat/=0)then
            CALL juDFT_error("problem with jobparam=",calledby="wann_read_inp")
             endif
            endif           
          ELSEIF(this%jobList(i).EQ.'mmn0at')THEN
             this%l_mmn0at=.TRUE.
          ELSEIF(this%jobList(i).EQ.'manyfiles')THEN
             this%l_manyfiles=.TRUE.
          ELSEIF(this%jobList(i).EQ.'collectmanyfiles')THEN
             this%l_collectmanyfiles=.TRUE.
          ELSEIF(this%jobList(i).EQ.'bestproj')THEN
             this%l_bestproj=.TRUE.
          ELSEIF(this%jobList(i).EQ.'pauli')THEN
             this%l_pauli=.TRUE.
          ELSEIF(this%jobList(i).EQ.'pauliat')THEN
             this%l_pauliat=.TRUE.
          ELSEIF(this%jobList(i).EQ.'proj_plot')THEN
             this%l_proj_plot=.TRUE.
          ELSEIF(this%jobList(i).EQ.'hopping')THEN
             this%l_hopping=.TRUE.
          ELSEIF(this%jobList(i).EQ.'plot_symm')THEN
             this%l_plot_symm=.TRUE.
          ELSEIF(this%jobList(i).EQ.'kptsreduc')THEN
             this%l_kptsreduc=.TRUE.
          ELSEIF(this%jobList(i).EQ.'fermi')THEN
             this%l_fermi=.TRUE.
          ELSEIF(this%jobList(i).EQ.'prepwan90')THEN
             this%l_prepwan90=.TRUE.
          ELSEIF(this%jobList(i).EQ.'plot_umdat')THEN
             this%l_plot_umdat=.TRUE.
          ELSEIF(this%jobList(i).EQ.'wann_plot')THEN
             this%l_wann_plot=.TRUE.
          ELSEIF(this%jobList(i).EQ.'bynumber')THEN
             this%l_bynumber=.TRUE.
          ELSEIF(jobname.EQ.'matrixmmn')THEN
            this%l_matrixmmn=.TRUE.
            if(l_param)then
             read(param,*,iostat=stat) this%matrixmmnfmt
             if(stat/=0)then
            CALL juDFT_error("problem with jobparam=",calledby="wann_read_inp")
             endif
            endif
          ELSEIF(jobname.EQ.'perpmag')THEN
            this%l_perpmag=.TRUE.
            if(l_param)then
             read(param,*,iostat=stat) this%perpmagfmt
             if(stat/=0)then
            CALL juDFT_error("problem with jobparam=",calledby="wann_read_inp")
             endif
            endif           
          ELSEIF(jobname.EQ.'torquers')THEN
            this%l_torquers=.TRUE.
            if(l_param)then
             read(param,*,iostat=stat) this%torquersfmt
             if(stat/=0)then
            CALL juDFT_error("problem with jobparam=",calledby="wann_read_inp")
             endif
            endif
          ELSEIF(jobname.EQ.'torque')THEN
            this%l_torque=.TRUE.
            if(l_param)then
             read(param,*,iostat=stat) this%torquefmt
             if(stat/=0)then
            CALL juDFT_error("problem with jobparam=",calledby="wann_read_inp")
             endif
          endif                    
            
            
            
          ELSEIF(jobname.EQ.'perpmagrs')THEN
            this%l_perpmagrs=.TRUE.
            if(l_param)then
             read(param,*,iostat=stat) this%perpmagrsfmt
             if(stat/=0)then
            CALL juDFT_error("problem with jobparam=",calledby="wann_read_inp")
             endif
            endif        
          ELSEIF(jobname.EQ.'anglmom')THEN
            this%l_anglmom=.TRUE.
            if(l_param)then
             read(param,*,iostat=stat) this%anglmomfmt
             if(stat/=0)then
            CALL juDFT_error("problem with jobparam=",calledby="wann_read_inp")
             endif
          endif
          ELSEIF(jobname.EQ.'anglmomrs')THEN
            this%l_anglmomrs=.TRUE.
            if(l_param)then
             read(param,*,iostat=stat) this%anglmomrsfmt
             if(stat/=0)then
            CALL juDFT_error("problem with jobparam=",calledby="wann_read_inp")
             endif
          endif                            
            
            
            
            
          ELSEIF(this%jobList(i).EQ.'projmethod')THEN
             this%l_projmethod=.TRUE.
          ELSEIF(this%jobList(i).EQ.'matrixamn')THEN
             this%l_matrixamn=.TRUE.
            if(l_param)then
             read(param,*,iostat=stat) this%matrixamnfmt
             if(stat/=0)then
            CALL juDFT_error("problem with jobparam=",calledby="wann_read_inp")
             endif
            endif
          ELSEIF(this%jobList(i).EQ.'wannierize')THEN
             this%l_wannierize=.TRUE.
          ELSEIF(this%jobList(i).EQ.'plotw90')THEN
             this%l_plotw90=.TRUE.
          ELSEIF(this%jobList(i).EQ.'dipole')THEN
             this%l_dipole=.TRUE.
          ELSEIF(this%jobList(i).EQ.'dipole2')THEN
             this%l_dipole2=.TRUE.      
          ELSEIF(this%jobList(i).EQ.'dipole3')THEN
             this%l_dipole3=.TRUE.
          ELSEIF(this%jobList(i).EQ.'ldauwan')THEN
             this%l_ldauwan=.TRUE.
          ELSEIF(this%jobList(i).EQ.'byenergy')THEN
             this%l_byenergy=.TRUE.
          ELSEIF(this%jobList(i).EQ.'finishnocoplot') THEN
             this%l_finishnocoplot=.TRUE.
          ELSEIF(this%jobList(i).EQ.'mmn0_unf_to_spn_unf') THEN
             this%l_mmn0_unf_to_spn_unf=.TRUE.
          ELSEIF(this%jobList(i).EQ.'mmn0_to_spn_unf') THEN
             this%l_mmn0_to_spn_unf=.TRUE.
          ELSEIF(this%jobList(i).EQ.'mmn0_to_spn') THEN
             this%l_mmn0_to_spn=.TRUE.
          ELSEIF(this%jobList(i).EQ.'mmn0_to_spn2') THEN
             this%l_mmn0_to_spn2=.TRUE.
          ELSEIF(this%jobList(i).EQ.'mmn0_unf_to_spn') THEN
             this%l_mmn0_unf_to_spn=.TRUE.            
          ELSEIF(this%jobList(i).EQ.'permag_unf_to_tor_unf') THEN
             this%l_perpmag_unf_to_tor_unf=.TRUE.
          ELSEIF(this%jobList(i).EQ.'perpmag_to_tor_unf') THEN
             this%l_perpmag_to_tor_unf=.TRUE.
          ELSEIF(this%jobList(i).EQ.'perpmag_to_tor') THEN
             this%l_perpmag_to_tor=.TRUE.
          ELSEIF(this%jobList(i).EQ.'perpmag_unf_to_tor') THEN
             this%l_perpmag_unf_to_tor=.TRUE.                    
          ELSEIF(this%jobList(i).EQ.'hsomtxvec_unf_to_lmpzsoc_unf') THEN
             this%l_hsomtxvec_unf_to_lmpzsoc_unf=.TRUE.
          ELSEIF(this%jobList(i).EQ.'hsomtxvec_to_lmpzsoc_unf') THEN
             this%l_hsomtxvec_to_lmpzsoc_unf=.TRUE.
          ELSEIF(this%jobList(i).EQ.'hsomtxvec_to_lmpzsoc') THEN
             this%l_hsomtxvec_to_lmpzsoc=.TRUE.
          ELSEIF(this%jobList(i).EQ.'hsomtxvec_unf_to_lmpzsoc') THEN
             this%l_hsomtxvec_unf_to_lmpzsoc=.TRUE.  
          ELSEIF(this%jobList(i).EQ.'hsomtx_unf_to_hsoc_unf') THEN
             this%l_hsomtx_unf_to_hsoc_unf=.TRUE.
          ELSEIF(this%jobList(i).EQ.'hsomtx_to_hsoc_unf') THEN
             this%l_hsomtx_to_hsoc_unf=.TRUE.
          ELSEIF(this%jobList(i).EQ.'hsomtx_to_hsoc') THEN
             this%l_hsomtx_to_hsoc=.TRUE.
          ELSEIF(this%jobList(i).EQ.'hsomtx_unf_to_hsoc') THEN
             this%l_hsomtx_unf_to_hsoc=.TRUE.  
          ELSEIF(this%jobList(i).EQ.'finishgwf') THEN
             this%l_finishgwf=.TRUE.
          ELSEIF(this%jobList(i).EQ.'skipkov') THEN
             this%l_skipkov=.TRUE.
          ELSEIF(this%jobList(i).EQ.'matrixuhu') THEN
             this%l_matrixuHu=.TRUE.
          ELSEIF(this%jobList(i).EQ.'matrixuhu-dmi') THEN
             this%l_matrixuHu_dmi=.TRUE.
        
          ELSEIF(jobname.EQ.'wan90version')THEN
             if(l_param)then
                read(param,*,iostat=stat) version_real
                if(stat/=0)then
                   write(*,*)"problem with jobparam=",param
                   CALL juDFT_error ("problem with jobparam", calledby="wann_read_inp")
                endif                                
             else
                  CALL juDFT_error ("parameter needed in wan90version", calledby="wann_read_inp")
             endif   
         
         
             if(abs(version_real-1.1).lt.1.e-9)THEN
                this%wan90version=1
             ELSEIF(abs(version_real-1.2).lt.1.e-9)THEN
                this%wan90version=2
             ELSEIF(abs(version_real-2.0).lt.1.e-9)THEN
                this%wan90version=3
             ELSEIF(abs(version_real-3.0).lt.1.e-9)THEN
                this%wan90version=30
             ELSEIF(abs(version_real-3.1).lt.1.e-9)THEN
                this%wan90version=31                
             ELSE
               CALL judft_error ("chosen w90 version unknown", calledby="wann_read_inp")
             endif
         !Not done
         ! ELSEIF(this%jobList(i).EQ.'ikptstart')THEN
         !    this%l_ikptstart=.TRUE.
         !    backspace(916)
         !    read(916,*,iostat=ios)task,this%ikptstart
         !    if (ios /= 0) CALL juDFT_error ("error reading ikptstart", calledby="wann_read_inp")
         !    if(l_p0)write(oUnit,*)"ikptstart=",this%ikptstart
         ELSEIF(this%jobList(i).EQ.'endjobs')THEN
             Exit
          ELSE
             WRITE(oUnit,*)"unrecognized key: ",this%jobList(i)
             CALL juDFT_error ("unrecognized key in wann_inp", calledby="wann_read_inp")
          END IF
       END DO
    END IF
 
    if(.not.this%l_atomlist)then
        allocate(this%atomlist(xml%get_nat()))
        do n=1,xml%get_nat()
          this%atomlist(n)=n
        enddo
        this%atomlist_num=xml%get_nat()
    else   
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
    endif

    if(this%l_unformatted)then
        this%socmatvecfmt=2
        this%socmatvecrsfmt=2
        this%anglmomrsfmt=2
        this%anglmomfmt=2
        this%torquefmt=2
        this%torquersfmt=2
        this%perpmagrsfmt=2
        this%perpmagfmt=2
        this%perpmagatfmt=2
        this%perpmagatrsfmt=2
        this%socmatrsfmt=2
        this%socmatfmt=2
        this%paulifmt=2
        this%pauliatfmt=2
        this%hoppingfmt=2
        this%matrixmmnfmt=2
        this%matrixamnfmt=2
        this%mmn0fmt=2
        this%mmn0atfmt=2
        this%matrixuHufmt=2
        this%matrixuHudmifmt=2
    endif


  END SUBROUTINE read_xml_wannier
END MODULE m_types_wannier
