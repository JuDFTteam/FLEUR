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
     INTEGER :: wan90version
     INTEGER :: oc_num_orbs
     INTEGER, ALLOCATABLE :: oc_orbs(:)
     LOGICAL :: l_unformatted
     LOGICAL :: l_oc_f
     LOGICAL :: l_ndegen
     LOGICAL :: l_orbitalmom
     LOGICAL :: l_orbcomp
     LOGICAL :: l_orbcomprs
     LOGICAL :: l_denmat
     LOGICAL :: l_perturbrs
     LOGICAL :: l_perturb
     LOGICAL :: l_nedrho
     LOGICAL :: l_anglmomrs
     LOGICAL :: l_anglmom
     LOGICAL :: l_spindisp
     LOGICAL :: l_spindisprs
     LOGICAL :: l_socspicom
     LOGICAL :: l_socspicomrs
     LOGICAL :: l_offdiposoprs
     LOGICAL :: l_offdiposop
     LOGICAL :: l_torque
     LOGICAL :: l_torquers
     LOGICAL :: l_atomlist
     INTEGER :: atomlist_num
     INTEGER, ALLOCATABLE :: atomlist(:)
     LOGICAL :: l_berry
     LOGICAL :: l_perpmagrs
     LOGICAL :: l_perpmag
     LOGICAL :: l_perpmagat
     LOGICAL :: l_perpmagatrs
     LOGICAL :: l_socmatrs
     LOGICAL :: l_socmat
     LOGICAL :: l_soctomom
     LOGICAL :: l_kptsreduc2
     LOGICAL :: l_nablapaulirs
     LOGICAL :: l_nablars
     LOGICAL :: l_surfcurr
     LOGICAL :: l_updown
     LOGICAL :: l_ahe
     LOGICAL :: l_she
     LOGICAL :: l_rmat
     LOGICAL :: l_nabla
     LOGICAL :: l_socodi
     LOGICAL :: l_pauli
     LOGICAL :: l_pauliat
     LOGICAL :: l_potmat
     LOGICAL :: l_projgen
     LOGICAL :: l_plot_symm
     LOGICAL :: l_socmmn0
     LOGICAL :: l_bzsym
     LOGICAL :: l_hopping
     LOGICAL :: l_kptsreduc
     LOGICAL :: l_prepwan90
     LOGICAL :: l_plot_umdat
     LOGICAL :: l_wann_plot
     LOGICAL :: l_bynumber
     LOGICAL :: l_stopopt
     LOGICAL :: l_matrixmmn
     LOGICAL :: l_matrixamn
     LOGICAL :: l_projmethod
     LOGICAL :: l_wannierize
     LOGICAL :: l_plotw90
     LOGICAL :: l_byindex
     LOGICAL :: l_byenergy
     LOGICAL :: l_proj_plot
     LOGICAL :: l_bestproj
     LOGICAL :: l_ikptstart
     LOGICAL :: l_lapw
     LOGICAL :: l_plot_lapw
     LOGICAL :: l_fermi
     LOGICAL :: l_dipole
     LOGICAL :: l_dipole2
     LOGICAL :: l_dipole3
     LOGICAL :: l_mmn0
     LOGICAL :: l_mmn0at
     LOGICAL :: l_manyfiles
     LOGICAL :: l_collectmanyfiles
     LOGICAL :: l_ldauwan
     LOGICAL :: l_lapw_kpts
     LOGICAL :: l_lapw_gfleur
     LOGICAL :: l_kpointgen
     LOGICAL :: l_w90kpointgen
     LOGICAL :: l_finishnocoplot
     LOGICAL :: l_finishgwf
     LOGICAL :: l_skipkov
     LOGICAL :: l_matrixuHu
     LOGICAL :: l_matrixuHu_dmi
     INTEGER :: ikptstart
     INTEGER :: band_min(1:2)
     INTEGER :: band_max(1:2)
     INTEGER :: gfthick
     INTEGER :: gfcut
     INTEGER :: unigrid(6)
     INTEGER :: mhp(3)
     !---> gwf
     LOGICAL :: l_ms
     LOGICAL :: l_sgwf
     LOGICAL :: l_socgwf
     LOGICAL :: l_gwf
     LOGICAL :: l_bs_comf
     LOGICAL :: l_exist
     LOGICAL :: l_opened
     LOGICAL :: l_cleverskip
     LOGICAL :: l_dim(3)
     REAL    :: scale_param
     REAL    :: aux_latt_const
     REAL    :: hdwf_t1
     REAL    :: hdwf_t2
     INTEGER :: nparampts
     CHARACTER(len=20) :: fn_eig
     CHARACTER(len=20) :: param_file
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
    CLASS(t_wann),INTENT(out):: this
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
