!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_oneD
  USE m_juDFT
  USE m_types_fleurinput_base
  IMPLICIT NONE
  PRIVATE
  ! types for 1D calculations
  TYPE od_dim
     LOGICAL :: d1=.FALSE.
     INTEGER :: mb=0, M=0, k3=0, m_cyl=0
     INTEGER :: chi=0, rot=0
     LOGICAL :: invs=.FALSE., zrfs=.FALSE.
     INTEGER :: n2d=0, nq2=0, nn2d=0
     INTEGER :: kimax2
     INTEGER :: nop=0, nat=0
  END TYPE od_dim
  
  TYPE od_inp
     LOGICAL :: d1
     INTEGER :: mb, M, k3, m_cyl
     INTEGER :: chi, rot
     LOGICAL :: invs, zrfs
     INTEGER :: n2d, nq2, nn2d
     INTEGER :: kimax2
     INTEGER, POINTER :: ig(:, :)  !(-k3:k3,-M:M)
     INTEGER, POINTER :: kv(:, :)        !(2,n2d)
     INTEGER, POINTER :: nst2(:)        !(n2d)
  END TYPE od_inp
  
  TYPE od_sym
     INTEGER :: nop, nat
     INTEGER, POINTER :: ngopr(:)     !(nat)
     REAL, POINTER :: mrot(:, :, :)  !(3,3,nop)
     REAL, POINTER :: tau(:, :)     !(3,nop)
     INTEGER, POINTER :: invtab(:)    !(nop)
     INTEGER, POINTER :: multab(:, :)  !(nop,nop)
  END TYPE od_sym
  
  TYPE od_lda
     INTEGER :: nn2d
     INTEGER, POINTER :: igf(:, :)  !(0:nn2d-1,2)
     REAL, POINTER :: pgf(:)    !(0:nn2d-1)
  END TYPE od_lda
  
  TYPE od_gga
     INTEGER          :: nn2d
     REAL, POINTER    :: pgfx(:)  ! (0:nn2d-1)
     REAL, POINTER    :: pgfy(:)
     REAL, POINTER    :: pgfxx(:)
     REAL, POINTER    :: pgfyy(:)
     REAL, POINTER    :: pgfxy(:)
  END TYPE od_gga
  TYPE,EXTENDS(t_fleurinput_base):: t_oneD
      TYPE(od_dim) :: odd
      TYPE(od_inp) :: odi
      TYPE(od_sym) :: ods
      TYPE(od_lda) :: odl
      TYPE(od_gga) :: odg
      INTEGER, POINTER :: ig1(:, :)
      INTEGER, POINTER :: kv1(:, :)
      INTEGER, POINTER :: nstr1(:)
      INTEGER, POINTER :: ngopr1(:)
      REAL, POINTER :: mrot1(:, :, :)
      REAL, POINTER :: tau1(:, :)
      INTEGER, POINTER :: invtab1(:)
      INTEGER, POINTER :: multab1(:, :)
      INTEGER, POINTER :: igfft1(:, :)
      REAL, POINTER :: pgfft1(:)
      REAL, POINTER :: pgft1x(:)
      REAL, POINTER :: pgft1y(:)
      REAL, POINTER :: pgft1xx(:)
      REAL, POINTER :: pgft1yy(:)
      REAL, POINTER :: pgft1xy(:)
    contains
      procedure :: read_xml=>read_xml_oneD
      PROCEDURE :: mpi_bc=>mpi_bc_oneD
      procedure :: init=>init_oneD
   END TYPE t_oneD
   PUBLIC::t_oneD,od_dim,od_inp,od_gga,od_lda,od_sym
 CONTAINS
   SUBROUTINE mpi_bc_oneD(this,mpi_comm,irank)
    use m_mpi_bc_tool
    class(t_oneD),INTENT(INOUT)::this
    integer,INTENT(IN):: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL::irank
    INTEGER ::rank
    if (present(irank)) THEN
       rank=irank
    else
       rank=0
    end if
    !Attention only few variables are broadcasted
    CALL mpi_bc(this%odd%d1 ,rank,mpi_comm)
    CALL mpi_bc(this%odd%M ,rank,mpi_comm)
    CALL mpi_bc(this%odd%mb ,rank,mpi_comm)
    CALL mpi_bc(this%odd%m_cyl ,rank,mpi_comm)
    CALL mpi_bc(this%odd%chi ,rank,mpi_comm)
    CALL mpi_bc(this%odd%rot ,rank,mpi_comm)
    CALL mpi_bc(this%odd%invs ,rank,mpi_comm)
    CALL mpi_bc(this%odd%zrfs ,rank,mpi_comm)


  END SUBROUTINE mpi_bc_oneD
   SUBROUTINE read_xml_oneD(this,xml)
     use m_types_xml
     class(t_oned),intent(inout)::this
     type(t_xml),intent(in)   ::xml

     
      ! Read in optional 1D parameters if present
     character(len=100):: xpathA
     integer :: numberNodes
     
      xPathA = '/fleurInput/calculationSetup/oneDParams'
      numberNodes = xml%GetNumberOfNodes(xPathA)

      this%odd%d1 = .FALSE.

      IF (numberNodes.EQ.1) THEN
         this%odd%d1 = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@d1'))
         this%odd%M = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@MM'))
         this%odd%mb = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@vM'))
         this%odd%m_cyl = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@m_cyl'))
         this%odd%chi = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@chi'))
         this%odd%rot = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@rot'))
         this%odd%invs = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@invs1'))
         this%odd%zrfs = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@zrfs1'))
      END IF

    END SUBROUTINE read_xml_oneD

    SUBROUTINE init_oneD(oneD,atoms)
      USE m_types_atoms
      CLASS(t_oned),INTENT(inout)::oneD
      TYPE(t_atoms),INTENT(in)::atoms


      ! Initialize missing 1D code arrays
      
      ALLOCATE (oneD%ig1(-oneD%odd%k3:oneD%odd%k3,-oneD%odd%M:oneD%odd%M))
      ALLOCATE (oneD%kv1(2,oneD%odd%n2d),oneD%nstr1(oneD%odd%n2d))
      ALLOCATE (oneD%ngopr1(atoms%nat),oneD%mrot1(3,3,oneD%odd%nop),oneD%tau1(3,oneD%odd%nop))
      ALLOCATE (oneD%invtab1(oneD%odd%nop),oneD%multab1(oneD%odd%nop,oneD%odd%nop))
      ALLOCATE (oneD%igfft1(0:oneD%odd%nn2d-1,2),oneD%pgfft1(0:oneD%odd%nn2d-1))
    end subroutine init_oneD
 END MODULE m_types_oneD
