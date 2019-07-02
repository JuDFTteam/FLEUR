!--------------------------------------------------------------------------------
! Copyright (c) 2017 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_postprocessInput

CONTAINS

SUBROUTINE postprocessInput(mpi,input,field,sym,stars,atoms,vacuum,kpts,&
     oneD,hybrid,cell,banddos,sliceplot,xcpot,forcetheo,forcetheo_data,&
     noco,DIMENSION,enpara,enparaxml,sphhar,l_opti,noel,l_kpts)

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
  use m_make_sym
  USE m_convn
  USE m_efield
  USE m_od_kptsgen
  USE m_types_forcetheo_extended
  USE m_types_xcpot_libxc
  USE m_types_xcpot_inbuild
  USE m_types_xcpot_inbuild_nofunction
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
  LOGICAL,          INTENT  (OUT) :: l_opti
  LOGICAL,          INTENT   (IN) :: l_kpts
  CHARACTER(len=3), ALLOCATABLE, INTENT(IN) :: noel(:)

  INTEGER              :: i, j, n, na, n1, n2, iType, l, ilo, ikpt
  INTEGER              :: minNeigd, nv, nv2, kq1, kq2, kq3, jrc, jsp, ii
  INTEGER              :: ios, ntst, ierr
  REAL                 :: sumWeight, rmtmax, zp, radius, dr
  REAL                 :: kmax1, dtild1, dvac1
  REAL                 :: bk(3)
  LOGICAL              :: l_vca, l_test,l_gga
 
  
  INTEGER, ALLOCATABLE :: jri1(:), lmax1(:)
  REAL,    ALLOCATABLE :: rmt1(:), dx1(:)
  INTEGER              ::func_vxc_id_c,func_vxc_id_x,func_exc_id_c,func_exc_id_x

#ifdef CPP_MPI
  INCLUDE 'mpif.h'
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Start of input postprocessing (calculate missing parameters)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call cell%init()
  cell%volint = cell%vol
  cell%volint = cell%volint - DOT_PRODUCT(atoms%volmts(:),atoms%neq(:))
  CALL sym%init(cell,input%film)
  CALL make_sym(sym,cell,atoms,noco,oneD,input) 
  !Finish setup of xcpot
  IF (xcpot%l_libxc) THEN
     func_vxc_id_c=xcpot%func_vxc_id_c
     func_vxc_id_x=xcpot%func_vxc_id_x
     func_exc_id_c=xcpot%func_exc_id_c
     func_exc_id_x=xcpot%func_exc_id_x
     DEALLOCATE(xcpot)
     ALLOCATE(t_xcpot_libxc::xcpot)
     SELECT TYPE(xcpot)
     CLASS is (t_xcpot_libxc)!just allocated like this
        CALL xcpot%init(func_vxc_id_x,func_vxc_id_c,func_exc_id_x,func_exc_id_c,input%jspins)
     END SELECT
  ELSE
     SELECT TYPE(xcpot)
     CLASS is (t_xcpot_inbuild_nf)
        CALL xcpot%init(atoms%ntype)
     CLASS DEFAULT
        CALL judft_error("Error in setup xcpot")
     END SELECT
  END IF

  !Finish setup of forcetheorem
  SELECT CASE (forcetheo_data%mode)
  CASE(1)
     ALLOCATE(t_forcetheo_mae::forcetheo)
  CASE(2)
     ALLOCATE(t_forcetheo_dmi::forcetheo)
  CASE(3)
     ALLOCATE(t_forcetheo_jij::forcetheo)
  CASE(4)
     ALLOCATE(t_forcetheo_ssdisp::forcetheo)
  CASE default
     ALLOCATE(t_forcetheo::forcetheo)
  END SELECT

  SELECT TYPE(forcetheo)
  TYPE IS(t_forcetheo_mae)
     CALL forcetheo%init(forcetheo_data%theta,forcetheo_data%phi,cell,sym)
  TYPE IS(t_forcetheo_dmi)
     CALL forcetheo%init(forcetheo_data%qvec,forcetheo_data%theta,forcetheo_data%phi)
  TYPE IS(t_forcetheo_jij)
     CALL forcetheo%init(forcetheo_data%qvec,forcetheo_data%theta(1),atoms)
  TYPE IS(t_forcetheo_ssdisp)
     CALL forcetheo%init(forcetheo_data%qvec)
  END SELECT
  

  !Generate enpara datatype
  CALL enpara%init_enpara(atoms,input%jspins,input%film,enparaXML)

  IF (mpi%irank.EQ.0) call check_input_switches(banddos,vacuum,noco,atoms,input)


  ! Check noco stuff and calculate missing noco parameters
  IF (noco%l_noco) THEN
     IF (noco%l_ss) THEN
        !--->    the angle beta is relative to the spiral in a spin-spiral
        !--->    calculation, i.e. if beta = 0 for all atoms in the unit cell
        !--->    that means that the moments are "in line" with the spin-spiral
        !--->    (beta = qss * taual). note: this means that only atoms within
        !--->    a plane perpendicular to qss can be equivalent!
        na = 1
        DO iType = 1,atoms%ntype
           noco%alph(iType) = noco%alphInit(iType) + tpi_const*dot_product(noco%qss,atoms%taual(:,na))
           na = na + atoms%neq(iType)
        END DO
     END IF
  ELSE
     IF (noco%l_ss) THEN
        CALL judft_warn("l_noco=F and l_ss=T is meaningless. Setting l_ss to F.")
        noco%l_ss = .FALSE.
     END IF
  END IF
  
    
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

  ! Generate missing parameters for atoms and calculate volume of the different regions
  CALL ylmnorm_init(atoms%lmaxd)
  dimension%nspd=(atoms%lmaxd+1+mod(atoms%lmaxd+1,2))*(2*atoms%lmaxd+1)
  dimension%msh = 0



  
     
  ! Initialize missing 1D code arrays

  ALLOCATE (oneD%ig1(-oneD%odd%k3:oneD%odd%k3,-oneD%odd%M:oneD%odd%M))
  ALLOCATE (oneD%kv1(2,oneD%odd%n2d),oneD%nstr1(oneD%odd%n2d))
  ALLOCATE (oneD%ngopr1(atoms%nat),oneD%mrot1(3,3,oneD%odd%nop),oneD%tau1(3,oneD%odd%nop))
  ALLOCATE (oneD%invtab1(oneD%odd%nop),oneD%multab1(oneD%odd%nop,oneD%odd%nop))
  ALLOCATE (oneD%igfft1(0:oneD%odd%nn2d-1,2),oneD%pgfft1(0:oneD%odd%nn2d-1))
  
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
  
  !DIMENSION%msh = MAXVAL(DIMENSION%msh,jrc)
   
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
  
  
  ! Other small stuff
  
  
  
  
  input%strho = .FALSE.
  
  INQUIRE(file="cdn1",exist=l_opti)
  if (noco%l_noco) INQUIRE(file="rhomat_inp",exist=l_opti)
  l_opti=.not.l_opti
  IF ((sliceplot%iplot).OR.(input%strho).OR.(input%swsp).OR.&
       (input%lflip).OR.(input%l_bmt)) l_opti = .TRUE.
  
  kpts%wtkpt=kpts%wtkpt/sum(kpts%wtkpt) !Normalize k-point weight
  
  
  CALL prp_xcfft(stars,input,cell,xcpot)
  
  
  IF (.NOT.sliceplot%iplot) THEN   
     IF (mpi%irank.EQ.0) THEN
        CALL convn(atoms,stars)
        CALL e_field(atoms,DIMENSION,stars,sym,vacuum,cell,input,field%efield)
     END IF !(mpi%irank.EQ.0)
  END IF

 
  CALL transform_by_moving_atoms(mpi,stars,atoms,vacuum, cell, sym, sphhar,input,oned,noco)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End of input postprocessing (calculate missing parameters)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE postprocessInput

END MODULE m_postprocessInput
