!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_sphhar
  !Data for the spherical harmonics
  USE m_judft
  USE m_types_fleur_setup
  USE m_json_tools
  IMPLICIT NONE
  TYPE,EXTENDS(t_fleursetup):: t_sphhar
     !No of symmetry types (must
     !equal maxval(atoms%ntypsy)
     INTEGER ::ntypsd
     !Max no of members of sphhar
     INTEGER ::memd
     !max of nlh
     INTEGER ::nlhd
     !No of sphhar (ntypsd)
     INTEGER,ALLOCATABLE ::nlh(:)
     !l's of sphhar (0:nlhd,ntypsd)
     INTEGER,ALLOCATABLE ::llh(:,:)
     !No of members in sphhar (0:nlh
     INTEGER,ALLOCATABLE ::nmem(:,:)
     !lm's of of members (max(nmem),
     INTEGER,ALLOCATABLE ::mlh(:,:,:)
     !phasefactors (max(nmem),0:nlhd
     COMPLEX,ALLOCATABLE ::clnu(:,:,:)
   CONTAINS
     PROCEDURE,PASS :: broadcast=>broadcast_sphhar
     PROCEDURE,PASS :: WRITE=>WRITE_sphhar
     PROCEDURE,PASS :: read=>READ_sphhar
     PROCEDURE,PASS :: read_xml=>read_xml_sphhar
  END TYPE t_sphhar
CONTAINS
  SUBROUTINE broadcast_sphhar(tt,mpi_comm,origin)
#ifdef CPP_MPI
    USE m_bc_tool
#endif    
    IMPLICIT NONE
    CLASS(t_sphhar),INTENT(INOUT):: tt
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

    CALL MPI_BCAST(tt%ntypsd,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%memd,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%nlhd,1,MPI_INTEGER,pe,mpi_comm,ierr)

    CALL MPI_BC(tt%nlh,pe,mpi_comm)
    CALL MPI_BC(tt%llh,pe,mpi_comm)
    CALL MPI_BC(tt%nmem,pe,mpi_comm)
    CALL MPI_BC(tt%mlh,pe,mpi_comm)
    CALL MPI_BC(tt%clnu,pe,mpi_comm)
#endif
  END SUBROUTINE broadcast_sphhar

  SUBROUTINE write_sphhar(tt, unit, iotype, v_list, iostat, iomsg)
    IMPLICIT NONE
    CLASS(t_sphhar),INTENT(IN):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    WRITE(unit,*,IOSTAT=iostat) '"sphhar":{'


    CALL json_print(unit,"ntypsd",tt%ntypsd)
    CALL json_print(unit,"memd",tt%memd)
    CALL json_print(unit,"nlhd",tt%nlhd)

    CALL json_print(unit,"nlh",tt%nlh)
    CALL json_print(unit,"llh",tt%llh)
    CALL json_print(unit,"nmem",tt%nmem)
    CALL json_print(unit,"mlh",tt%mlh)
    CALL json_print(unit,"clnu%re",real(tt%clnu))
    CALL json_print(unit,"clnu%im",aimag(tt%clnu))
    WRITE(unit,*,IOSTAT=iostat) '}'

  END SUBROUTINE write_sphhar
  SUBROUTINE read_sphhar(tt, unit, iotype, v_list, iostat, iomsg)
    use m_json_tools
    IMPLICIT NONE
    CLASS(t_sphhar),INTENT(INOUT):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    INTEGER :: ntype
    REAL,ALLOCATABLE::rt_r(:,:,:),rt_i(:,:,:)
    CALL json_open_class("sphhar",unit,iostat)
    IF (iostat.NE.0)   RETURN


    CALL json_read(unit,"ntypsd",tt%ntypsd)
    CALL json_read(unit,"memd",tt%memd)
    CALL json_read(unit,"nlhd",tt%nlhd)

    CALL json_read(unit,"nlh",tt%nlh)
    CALL json_read(unit,"llh",tt%llh)
    CALL json_read(unit,"nmem",tt%nmem)
    CALL json_read(unit,"mlh",tt%mlh)
    
    CALL json_read(unit,"clnu%re",rt_r)
    CALL json_read(unit,"clnu%im",rt_i)
    tt%clnu=CMPLX(rt_r,rt_i)
    CALL json_close_class(unit,iostat)

  END SUBROUTINE read_sphhar

  SUBROUTINE read_xml_sphhar(tt)
    IMPLICIT NONE
    CLASS(t_sphhar),INTENT(OUT):: tt
  END SUBROUTINE read_xml_sphhar


END MODULE m_types_sphhar
