!--------------------------------------------------------------------------------
! Copyright (c) 2017 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_postprocessInput

CONTAINS

SUBROUTINE postprocessInput(mpi,input,field,sym,stars,atoms,vacuum,kpts,&
     oneD,hybrid,cell,banddos,sliceplot,xcpot,forcetheo,forcetheo_data,&
     noco,DIMENSION,enpara,enparaxml,sphhar,l_kpts)

  USE m_juDFT
  USE m_types
  USE m_constants
  USE m_ylm
  USE m_chkmt
  USE m_dwigner
  USE m_cdn_io
  USE m_prpxcfft
  use m_checks
  use m_lapwdim
  use m_make_stars
  use m_make_sphhar
  use m_make_forcetheo
  use m_make_xcpot
  use m_make_sym
  USE m_convn
  USE m_efield
  USE m_od_kptsgen
  USE m_relaxio


  IMPLICIT NONE

  TYPE(t_mpi)      ,INTENT   (IN) :: mpi
  CLASS(t_forcetheo),ALLOCATABLE,INTENT(OUT):: forcetheo
  TYPE(t_forcetheo_data),INTENT(IN):: forcetheo_data
  TYPE(t_input),    INTENT(INOUT) :: input
  TYPE(t_sym),      INTENT(INOUT) :: sym
  TYPE(t_stars),    INTENT(INOUT) :: stars
  TYPE(t_atoms),    INTENT(INOUT) :: atoms
  TYPE(t_vacuum),   INTENT(INOUT) :: vacuum
  TYPE(t_kpts),     INTENT(INOUT) :: kpts
  TYPE(t_oneD),     INTENT(INOUT) :: oneD
  TYPE(t_hybrid),   INTENT(INOUT) :: hybrid
  TYPE(t_cell),     INTENT(INOUT) :: cell
  TYPE(t_banddos),  INTENT(INOUT) :: banddos
  TYPE(t_sliceplot),INTENT(INOUT) :: sliceplot
  CLASS(t_xcpot),ALLOCATABLE,INTENT(INOUT) :: xcpot
  TYPE(t_noco),     INTENT(INOUT) :: noco
  TYPE(t_dimension),INTENT(INOUT) :: dimension
  TYPE(t_enparaXML)   ,INTENT(IN):: enparaXML
  TYPE(t_enpara)   ,INTENT(OUT)   :: enpara
  TYPE(t_sphhar)   ,INTENT  (OUT) :: sphhar
  TYPE(t_field),    INTENT(INOUT) :: field
  LOGICAL,          INTENT   (IN) :: l_kpts

  INTEGER              :: i, j, n, na, n1, n2, iType, l, ilo, ikpt
  INTEGER              :: minNeigd, nv, nv2, kq1, kq2, kq3, jrc, jsp, ii
  INTEGER              :: ios, ntst, ierr
  REAL                 :: rmtmax, zp, radius, dr
  LOGICAL              :: l_vca, l_test


  INTEGER, ALLOCATABLE :: jri1(:), lmax1(:)
  REAL,    ALLOCATABLE :: rmt1(:), dx1(:)

#ifdef CPP_MPI
  INCLUDE 'mpif.h'
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Start of input postprocessing (calculate missing parameters)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call cell%init(DOT_PRODUCT(atoms%volmts(:),atoms%neq(:)))
  call atoms%init(cell)
  CALL sym%init(cell,input%film)
  call vacuum%init(sym)

  CALL enpara%init_enpara(atoms,input%jspins,input%film,enparaXML)
  CALL make_sym(sym,cell,atoms,noco,oneD,input)
  call make_forcetheo(forcetheo_data,cell,sym,atoms,forcetheo)
  call make_xcpot(xcpot,atoms,input)

  IF (mpi%irank.EQ.0) call check_input_switches(banddos,vacuum,noco,atoms,input)

  ! Generate missing general parameters

  minNeigd = MAX(5,NINT(0.75*input%zelec) + 1)
  IF (noco%l_soc.and.(.not.noco%l_noco)) minNeigd = 2 * minNeigd
  IF (noco%l_soc.and.noco%l_ss) minNeigd=(3*minNeigd)/2
  IF ((dimension%neigd.NE.-1).AND.(dimension%neigd.LT.minNeigd)) THEN
     IF (dimension%neigd>0) THEN
        WRITE(*,*) 'numbands is too small. Setting parameter to default value.'
        WRITE(*,*) 'changed numbands (dimension%neigd) to ',minNeigd
     ENDIF
     dimension%neigd = minNeigd
  END IF

  CALL lapw_dim(kpts,cell,input,noco,oneD,forcetheo,DIMENSION)

  IF(dimension%neigd.EQ.-1) THEN
     dimension%neigd = dimension%nvd + atoms%nlotot
  END IF

  IF (noco%l_noco) dimension%neigd = 2*dimension%neigd



  CALL ylmnorm_init(atoms%lmaxd)


  call oneD%init(atoms)
  ! Initialize missing hybrid functionals arrays
  ALLOCATE (hybrid%nindx(0:atoms%lmaxd,atoms%ntype))
  ! Check muffin tin radii
  l_test = .TRUE. ! only checking, dont use new parameters
  CALL chkmt(atoms,input,vacuum,cell,oneD,l_test)

  !adjust positions by displacements
  CALL apply_displacements(cell,input,vacuum,oneD,sym,noco,atoms)

  call make_sphhar(atoms,sphhar,sym,cell,oneD)
  CALL make_stars(stars,sym,atoms,vacuum,sphhar,input,cell,xcpot,oneD,noco,mpi)

  ! Store structure data
  CALL storeStructureIfNew(input,stars, atoms, cell, vacuum, oneD, sym, mpi,sphhar,noco)

  !Adjust kpoints in case of DOS

  IF ( banddos%dos .AND. banddos%ndir == -3 ) THEN
     WRITE(*,*) 'Recalculating k point grid to cover the full BZ.'
     !CALL gen_bz(kpts,sym)
     kpts%nkpt = kpts%nkptf
     DEALLOCATE(kpts%bk,kpts%wtkpt)
     ALLOCATE(kpts%bk(3,kpts%nkptf),kpts%wtkpt(kpts%nkptf))
     kpts%bk(:,:) = kpts%bkf(:,:)
     kpts%wtkpt = 1.0 / kpts%nkptf
  END IF




  CALL prp_xcfft(stars,input,cell,xcpot)


  IF (.NOT.sliceplot%iplot) THEN
     IF (mpi%irank.EQ.0) THEN
        CALL convn(atoms,stars)
        CALL e_field(atoms,DIMENSION,stars,sym,vacuum,cell,input,field%efield)
     END IF !(mpi%irank.EQ.0)
  END IF

  !At some point this should be enabled for noco as well
#ifdef CPP_MPI
  CALL MPI_BCAST(atoms%nat,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
#endif
  IF (.not.noco%l_noco) &
  CALL transform_by_moving_atoms(mpi,stars,atoms,vacuum, cell, sym, sphhar,input,oned,noco)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End of input postprocessing (calculate missing parameters)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call oneD%init(atoms)

END SUBROUTINE postprocessInput

END MODULE m_postprocessInput
