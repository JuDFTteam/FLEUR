!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_field
  USE m_judft
  USE m_types_fleur_setup
  USE m_json_tools

  !this module contains actually two  types.
  !t_efield is a subtype of the more general t_fields

  
  IMPLICIT NONE
  TYPE,EXTENDS(t_fleursetup):: t_efield
     REAL    :: zsigma  = 10.0  ! Distance to the charged plates
     REAL,POINTER    :: sigma   ! charge at the plates
     REAL    :: sig_b(2)=  0.0  ! Extra charge for the top/bottom plate
     COMPLEX :: vslope  =  0.0  ! Dirichlet bnd. cond.: Slope
     REAL,    ALLOCATABLE :: sigEF(:,:,:) ! (nx, ny, nvac)
     COMPLEX, ALLOCATABLE :: rhoEF(:,:)   ! (g_||, nvac)
     COMPLEX, ALLOCATABLE :: C1(:), C2(:) ! Coeff. for Dirichlet bnd.cond.
     LOGICAL :: l_segmented = .FALSE.
     LOGICAL :: plot_charge = .FALSE. ! Plot charge as inputted
     LOGICAL :: plot_rho    = .FALSE. ! Plot Fourier-transformed charge
     LOGICAL :: autocomp    = .TRUE.  ! Auto-compensate film charge
     LOGICAL :: dirichlet = .FALSE. ! Dirichlet vs. Neumann boundary cond.
     LOGICAL :: l_dirichlet_coeff = .FALSE. ! For MPI, true if C1/C2 set
   CONTAINS
     PROCEDURE,PASS :: broadcast=>broadcast_efield
     PROCEDURE,PASS :: WRITE=>WRITE_efield
     PROCEDURE,PASS :: READ=>READ_efield
     PROCEDURE,PASS :: read_xml=>read_xml_efield
  END TYPE t_efield

  TYPE,EXTENDS(t_fleursetup):: t_field
     TYPE(t_efield)   :: efield
     LOGICAL          :: l_b_field=.FALSE.
     REAL             :: b_field
     REAL,ALLOCATABLE :: b_field_mt(:)
   CONTAINS
     PROCEDURE,PASS :: broadcast=>broadcast_field
     PROCEDURE,PASS :: WRITE=>WRITE_field
     PROCEDURE,PASS :: READ=>READ_field
     PROCEDURE,PASS :: read_xml=>read_xml_field
  END TYPE t_efield

  
CONTAINS
  SUBROUTINE broadcast_efield(tt,mpi_comm,origin)
#ifdef CPP_MPI
    USE m_mpi_bc_tool
#endif    
    IMPLICIT NONE
    CLASS(t_efield),INTENT(INOUT):: tt
    INTEGER,INTENT(IN)               :: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL      :: origin

#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER :: pe,ierr,ntype,irank

    CALL MPI_COMM_RANK(mpi_comm,irank,ierr)

    IF (PRESENT(origin)) THEN
       pe=origin
    ELSE
       pe=0
    ENDIF
    CALL MPI_BCAST(tt%zsigma,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%sigma,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%sig_b,2,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%vslope,1,MPI_DOUBLE_COMPLEX,pe,mpi_comm,ierr)
 
    CALL MPI_BC(tt%sigef,pe,mpi_comm)
    CALL MPI_BC(tt%rhoef,pe,mpi_comm)
    CALL MPI_BC(tt%c1,pe,mpi_comm)
    CALL MPI_BC(tt%c2,pe,mpi_comm)
    
    CALL MPI_BCAST(tt%l_segmented,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%plot_charge,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%plot_rho,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%autocomp,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%dirichlet,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%l_dirichlet_coeff,1,MPI_LOGICAL,pe,mpi_comm,ierr)

#endif

  END SUBROUTINE broadcast_efield

  SUBROUTINE broadcast_field(tt,mpi_comm,origin)
#ifdef CPP_MPI
    USE m_mpi_bc_tool
#endif    
    IMPLICIT NONE
    CLASS(t_field),INTENT(INOUT):: tt
    INTEGER,INTENT(IN)               :: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL      :: origin

#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER :: pe,ierr,ntype,irank

    CALL MPI_COMM_RANK(mpi_comm,irank,ierr)

    IF (PRESENT(origin)) THEN
       pe=origin
    ELSE
       pe=0
    ENDIF
    CALL tt%efield%broadcast(mpi)
    
    CALL MPI_BCAST(tt%l_b_field,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%b_field,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BC(tt%b_field_mt,pe,mpi_comm)
    
   
#endif

  END SUBROUTINE broadcast_field

  
  SUBROUTINE write_efield(tt, unit, iotype, v_list, iostat, iomsg)
    IMPLICIT NONE
    CLASS(t_efield),INTENT(IN):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    INTEGER:: ntype
    REAL,ALLOCATABLE::d_wgn(:,:)

    WRITE(unit,*,IOSTAT=iostat) '"efield":{'

    CALL json_print(unit,"zsigma",tt%zsigma)
    CALL json_print(unit,"sigma",tt%sigma)
    CALL json_print(unit,"sig_b",tt%sig_b)
    CALL json_print(unit,"vslope",tt%vslope)
 
    CALL json_print(unit,"sigef",tt%sigef)
    CALL json_print(unit,"rhoef",tt%rhoef)
    CALL json_print(unit,"c1",tt%c1)
    CALL json_print(unit,"c2",tt%c2)
    
    CALL json_print(unit,"l_segmented",tt%l_segmented)
    CALL json_print(unit,"plot_charge",tt%plot_charge)
    CALL json_print(unit,"plot_rho",tt%plot_rho)
    CALL json_print(unit,"autocomp",tt%autocomp)
    CALL json_print(unit,"dirichlet",tt%dirichlet)
    CALL json_print(unit,"l_dirichlet_coeff",tt%l_dirichlet_coeff)


    WRITE(unit,*,IOSTAT=iostat) '}'

  END SUBROUTINE write_efield


 SUBROUTINE write_field(tt, unit, iotype, v_list, iostat, iomsg)
    IMPLICIT NONE
    CLASS(t_field),INTENT(IN):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    INTEGER:: ntype
    REAL,ALLOCATABLE::d_wgn(:,:)

    WRITE(unit,*,IOSTAT=iostat) '"field":{'

 
    CALL json_print(unit,"l_b_field",tt%l_b_field)
    CALL json_print(unit,"b_field",tt%b_field)
    CALL json_print(unit,"b_field_mt",tt%b_field_mt)

    call tt%efield%write(unit, iotype, v_list, iostat, iomsg)
    

    WRITE(unit,*,IOSTAT=iostat) '}'

  END SUBROUTINE write_field
 

  SUBROUTINE read_efield(tt, unit, iotype, v_list, iostat, iomsg)
    USE m_json_tools
    IMPLICIT NONE
    CLASS(t_efield),INTENT(INOUT):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    INTEGER :: ntype
    REAL,ALLOCATABLE::rt(:,:)
    CALL json_open_class("efield",unit,iostat)
    IF (iostat.NE.0)   RETURN

    CALL json_read(unit,"zsigma",tt%zsigma)
    CALL json_read(unit,"sigma",tt%sigma)
    CALL json_read(unit,"sig_b",tt%sig_b)
    CALL json_read(unit,"vslope",tt%vslope)
 
    CALL json_read(unit,"sigef",tt%sigef)
    CALL json_read(unit,"rhoef",tt%rhoef)
    CALL json_read(unit,"c1",tt%c1)
    CALL json_read(unit,"c2",tt%c2)
    
    CALL json_read(unit,"l_segmented",tt%l_segmented)
    CALL json_read(unit,"plot_charge",tt%plot_charge)
    CALL json_read(unit,"plot_rho",tt%plot_rho)
    CALL json_read(unit,"autocomp",tt%autocomp)
    CALL json_read(unit,"dirichlet",tt%dirichlet)
    CALL json_read(unit,"l_dirichlet_coeff",tt%l_dirichlet_coeff)


    CALL json_close_class(unit,iostat)

  END SUBROUTINE read_efield


  SUBROUTINE read_field(tt, unit, iotype, v_list, iostat, iomsg)
    USE m_json_tools
    IMPLICIT NONE
    CLASS(t_field),INTENT(INOUT):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    INTEGER :: ntype
    REAL,ALLOCATABLE::rt(:,:)
    CALL json_open_class("field",unit,iostat)
    IF (iostat.NE.0)   RETURN

    CALL json_read(unit,"l_b_field",tt%l_b_field)
    CALL json_read(unit,"b_field",tt%b_field)
    CALL json_read(unit,"b_field_mt",tt%b_field_mt)

    CALL tt%efield%read(unit, iotype, v_list, iostat, iomsg)

    CALL json_close_class(unit,iostat)

  END SUBROUTINE read_field


  SUBROUTINE read_xml_efield(tt)
    USE m_xmlIntWrapFort
    USE m_constants
    USE m_calculator
    USE m_inp_xml
    IMPLICIT NONE
    CLASS(t_efield),INTENT(OUT):: tt

    LOGICAL :: tempBool
    INTEGER:: numberNodes,numtokens
    CHARACTER(len=50):: xpathA,xpathb,valuestring

    xPathA = '/fleurInput/calculationSetup/eField'
    numberNodes = xmlGetNumberOfNodes(xPathA)

    IF (numberNodes.EQ.1) THEN
         tt%zsigma = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@zsigma'))
         tt%sig_b(1) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@sig_b_1'))
         tt%sig_b(2) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@sig_b_2'))
         tt%plot_charge = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@plot_charge'))
         tt%plot_rho = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@plot_rho'))
         tt%autocomp = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@autocomp'))
         tt%dirichlet = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@dirichlet'))
         !l_eV = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@eV'))

         CALL judft_error('Error: Reading input for E-Fields not yet implemented completely!')
         !      ALLOCATE(input%efield%sigEF(3*k1d, 3*k2d, nvac))
         !      input%efield%sigEF = 0.0
         !IF (l_eV) THEN
         !   input%efield%sig_b(:) = input%efield%sig_b/hartree_to_ev_const
         !         input%efield%sigEF(:,:,:) = input%efield%sigEF/hartree_to_ev_const
         !END IF
      END IF


  END SUBROUTINE read_xml_efield

  SUBROUTINE read_xml_field(tt)
    USE m_xmlIntWrapFort
    USE m_constants
    USE m_calculator
    USE m_inp_xml
    IMPLICIT NONE
    CLASS(t_field),INTENT(OUT):: tt

    LOGICAL :: tempBool
    INTEGER:: numberNodes,numtokens
    CHARACTER(len=50):: xpathA,xpathb,valuestring

    CALL judft_error("IO of fields not implemented")
   

  END SUBROUTINE read_xml_field



END MODULE m_types_field
