!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_wannier
  USE m_judft
  USE m_types_fleur_setup
  use m_json_tools
  IMPLICIT NONE

  !This is the old wannier type it needs to be refactored into the new template below..
  TYPE t_wann
     INTEGER :: wan90version
     INTEGER :: oc_num_orbs
     INTEGER,ALLOCATABLE :: oc_orbs(:)
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
     INTEGER,ALLOCATABLE :: atomlist(:)
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
     REAL,ALLOCATABLE :: param_vec(:,:)
     REAL,ALLOCATABLE :: param_alpha(:,:)
     CHARACTER(LEN=20), ALLOCATABLE :: jobList(:)
     !---> gwf

  END TYPE t_wann

  
  TYPE,EXTENDS(t_fleursetup):: t_wannier
     TYPE(t_wannier):: wann
   CONTAINS
     PROCEDURE,PASS :: broadcast=>broadcast_wannier
     PROCEDURE,PASS :: write=>WRITE_wannier
     PROCEDURE,PASS :: read=>READ_wannier
     PROCEDURE,PASS :: read_xml=>read_xml_wannier
  END TYPE t_wannier

CONTAINS
  SUBROUTINE broadcast_wannier(tt,mpi_comm,origin)
    IMPLICIT NONE
    CLASS(t_wannier),INTENT(INOUT):: tt
    INTEGER,INTENT(IN)               :: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL      :: origin

#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER :: pe,ierr
    
    IF (PRESENT(origin)) THEN
       pe=origin
    ELSE
       pe=0
    ENDIF

    call judft_error("Not implemented yet")
    !CALL MPI_BCAST(tt%bbmat,9,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
#endif
      
  END SUBROUTINE broadcast_wannier

  SUBROUTINE write_wannier(tt, unit, iotype, v_list, iostat, iomsg)
    IMPLICIT NONE
    CLASS(t_wannier),INTENT(IN):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    WRITE(unit,*,IOSTAT=iostat) '"wannier":{'

    call judft_error("Not implemented yet")

!    call json_print(unit,"omtil",tt%omtil)
 
    WRITE(unit,*,IOSTAT=iostat) '}'
    
  END SUBROUTINE write_wannier
  SUBROUTINE read_wannier(tt, unit, iotype, v_list, iostat, iomsg)
    IMPLICIT NONE
    CLASS(t_wannier),INTENT(INOUT):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    CHARACTER(len=40)::string
    REAL,allocatable:: rtemp(:)
    CALL json_open_class("wannier",unit,iostat)
    IF (iostat.NE.0)   RETURN

    call judft_error("Not implemented yet")
    
    CALL json_close_class(unit,iostat)
    
  END SUBROUTINE read_wannier

 


  SUBROUTINE read_xml_wannier(tt)
    USE m_xmlIntWrapFort
    USE m_constants
    USE m_calculator
    USE m_inv3
    IMPLICIT NONE
    CLASS(t_wannier),INTENT(OUT):: tt


    LOGICAL::film
    CHARACTER(len=200):: xpath,valueString
    REAL              :: atemp, a1(3),a2(3),a3(3),latticeScale

    call judft_error("Not implemented yet")

    
  END SUBROUTINE read_xml_wannier
  
  
END MODULE m_types_wannier
