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
     INTEGER :: wan90version =3
     INTEGER :: oc_num_orbs =0
     INTEGER, ALLOCATABLE :: oc_orbs(:)
     LOGICAL :: l_unformatted =.false.
     LOGICAL :: l_oc_f=.false.
     LOGICAL :: l_ndegen=.false.
     LOGICAL :: l_orbitalmom=.false.
     LOGICAL :: l_orbcomp=.false.
     LOGICAL :: l_orbcomprs=.false.
     LOGICAL :: l_denmat=.false.
     LOGICAL :: l_perturbrs=.false.
     LOGICAL :: l_perturb=.false.
     LOGICAL :: l_nedrho=.false.
     LOGICAL :: l_anglmomrs=.false.
     LOGICAL :: l_anglmom=.false.
     LOGICAL :: l_spindisp=.false.
     LOGICAL :: l_spindisprs=.false.
     LOGICAL :: l_socspicom=.false.
     LOGICAL :: l_socspicomrs=.false.
     LOGICAL :: l_offdiposoprs=.false.
     LOGICAL :: l_offdiposop=.false.
     LOGICAL :: l_torque=.false.
     LOGICAL :: l_torquers=.false.
     LOGICAL :: l_atomlist=.false.
     INTEGER :: atomlist_num=0
     INTEGER, ALLOCATABLE :: atomlist(:)
     LOGICAL :: l_berry=.false.
     LOGICAL :: l_perpmagrs=.false.
     LOGICAL :: l_perpmag=.false.
     LOGICAL :: l_perpmagat=.false.
     LOGICAL :: l_perpmagatrs=.false.
     LOGICAL :: l_socmatrs=.false.
     LOGICAL :: l_socmat=.false.
     LOGICAL :: l_soctomom=.false.
     LOGICAL :: l_kptsreduc2=.false.
     LOGICAL :: l_nablapaulirs=.false.
     LOGICAL :: l_nablars=.false.
     LOGICAL :: l_surfcurr=.false.
     LOGICAL :: l_updown=.false.
     LOGICAL :: l_ahe=.false.
     LOGICAL :: l_she=.false.
     LOGICAL :: l_rmat=.false.
     LOGICAL :: l_nabla=.false.
     LOGICAL :: l_socodi=.false.
     LOGICAL :: l_pauli=.false.
     LOGICAL :: l_pauliat=.false.
     LOGICAL :: l_potmat=.false.
     LOGICAL :: l_projgen=.false.
     LOGICAL :: l_plot_symm=.false.
     LOGICAL :: l_socmmn0=.false.
     LOGICAL :: l_bzsym=.false.
     LOGICAL :: l_hopping=.false.
     LOGICAL :: l_kptsreduc=.false.
     LOGICAL :: l_prepwan90=.false.
     LOGICAL :: l_plot_umdat=.false.
     LOGICAL :: l_wann_plot=.false.
     LOGICAL :: l_bynumber=.false.
     LOGICAL :: l_stopopt=.false.
     LOGICAL :: l_matrixmmn=.false.
     LOGICAL :: l_matrixamn=.false.
     LOGICAL :: l_projmethod=.false.
     LOGICAL :: l_wannierize=.false.
     LOGICAL :: l_plotw90=.false.
     LOGICAL :: l_byindex=.false.
     LOGICAL :: l_byenergy=.false.
     LOGICAL :: l_proj_plot=.false.
     LOGICAL :: l_bestproj=.false.
     LOGICAL :: l_ikptstart=.false.
     LOGICAL :: l_lapw=.false.
     LOGICAL :: l_plot_lapw=.false.
     LOGICAL :: l_fermi=.false.
     LOGICAL :: l_dipole=.false.
     LOGICAL :: l_dipole2=.false.
     LOGICAL :: l_dipole3=.false.
     LOGICAL :: l_mmn0=.false.
     LOGICAL :: l_mmn0at=.false.
     LOGICAL :: l_manyfiles=.false.
     LOGICAL :: l_collectmanyfiles=.false.
     LOGICAL :: l_ldauwan=.false.
     LOGICAL :: l_lapw_kpts=.false.
     LOGICAL :: l_lapw_gfleur=.false.
     LOGICAL :: l_kpointgen=.false.
     LOGICAL :: l_w90kpointgen=.false.
     LOGICAL :: l_finishnocoplot=.false.
     LOGICAL :: l_finishgwf=.false.
     LOGICAL :: l_skipkov=.false.
     LOGICAL :: l_matrixuHu=.false.
     LOGICAL :: l_matrixuHu_dmi=.false.
     INTEGER :: ikptstart=1
     INTEGER :: band_min(1:2)=-1
     INTEGER :: band_max(1:2)=-1
     INTEGER :: gfthick=0
     INTEGER :: gfcut=0
     INTEGER :: unigrid(6)=0
     INTEGER :: mhp(3)=0
     !---> gwf
     LOGICAL :: l_ms=.false.
     LOGICAL :: l_sgwf=.false.
     LOGICAL :: l_socgwf=.false.
     LOGICAL :: l_gwf=.false.
     LOGICAL :: l_bs_comf=.false.
     LOGICAL :: l_exist=.false.
     LOGICAL :: l_opened=.false.
     LOGICAL :: l_cleverskip=.false.
     LOGICAL :: l_dim(3)=.false.
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
  END TYPE t_wann
  
  PUBLIC t_wann
CONTAINS
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
       wannAtomList(i)= evaluateFirstBoolOnly(xml%getAttributeValue(xml%posPath(i))//'/@wannier')
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
