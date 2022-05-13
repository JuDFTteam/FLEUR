!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_hybinp
   USE m_judft
   USE m_types_fleurinput_base
   IMPLICIT NONE
   PRIVATE

   TYPE, EXTENDS(t_fleurinput_base):: t_hybinp
      LOGICAL                ::  l_hybrid = .false.
      INTEGER                ::  ewaldlambda = 3
      INTEGER                ::  lexp = 16
      INTEGER                ::  bands1 = -1 !Only read in
      real                   ::  fftcut = 1.0 ! 1.0 is the exact crit (2./3. seem also ok, but less stable)
      INTEGER, ALLOCATABLE   ::  select1(:, :)
      INTEGER, ALLOCATABLE   ::  lcutm1(:)
      INTEGER, ALLOCATABLE   ::  lcutwf(:)
      INTEGER, ALLOCATABLE   ::  map(:, :)
      INTEGER, ALLOCATABLE   ::  tvec(:, :, :)
      !REAL, ALLOCATABLE      ::  radbasfn_mt(:,:,:,:)
      COMPLEX, ALLOCATABLE   ::  d_wgn2(:, :, :, :)

   CONTAINS
      PROCEDURE :: read_xml => read_xml_hybinp
      PROCEDURE :: mpi_bc => mpi_bc_hybinp
      PROCEDURE :: init => init_hybinp
      PROCEDURE :: gen_map => gen_map_hybinp
   END TYPE t_hybinp
   PUBLIC t_hybinp

CONTAINS

   SUBROUTINE mpi_bc_hybinp(this, mpi_comm, irank)
      USE m_mpi_bc_tool
      CLASS(t_hybinp), INTENT(INOUT)::this
      INTEGER, INTENT(IN):: mpi_comm
      INTEGER, INTENT(IN), OPTIONAL::irank
      INTEGER ::rank
      IF (PRESENT(irank)) THEN
         rank = irank
      ELSE
         rank = 0
      END IF
      CALL mpi_bc(this%l_hybrid, rank, mpi_comm)
      CALL mpi_bc(this%ewaldlambda, rank, mpi_comm)
      CALL mpi_bc(this%lexp, rank, mpi_comm)
      CALL mpi_bc(this%bands1, rank, mpi_comm)
      CALL mpi_bc(this%select1, rank, mpi_comm)
      CALL mpi_bc(this%lcutm1, rank, mpi_comm)
      CALL mpi_bc(this%lcutwf, rank, mpi_comm)
      CALL mpi_bc(this%map, rank, mpi_comm)
      CALL mpi_bc(this%tvec, rank, mpi_comm)
      CALL mpi_bc(this%d_wgn2, rank, mpi_comm)
      call mpi_bc(this%fftcut, rank, mpi_comm)
   END SUBROUTINE mpi_bc_hybinp

   SUBROUTINE read_xml_hybinp(this, xml)
      USE m_types_xml
      CLASS(t_hybinp), INTENT(INout):: this
      TYPE(t_xml),INTENT(INOUT) ::xml

      INTEGER::numberNodes, ntype, itype
      CHARACTER(len=100)  :: xPathA
      CHARACTER(len=4), allocatable  :: xc_name

      ntype = xml%GetNumberOfNodes('/fleurInput/atomGroups/atomGroup')
      ALLOCATE (this%lcutm1(ntype), source=0)
      ALLOCATE (this%lcutwf(ntype), source=0)
      ALLOCATE (this%select1(4, ntype), source=0)
      numberNodes = xml%GetNumberOfNodes('/fleurInput/calculationSetup/prodBasis')
      IF (numberNodes == 1) THEN
         this%ewaldlambda = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/prodBasis/@ewaldlambda'))
         this%lexp = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/prodBasis/@lexp'))
         this%bands1 = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/prodBasis/@bands'))
         IF(xml%versionNumber>=34) THEN
            this%fftcut = evaluateFirstOnly (xml%GetAttributeValue('/fleurInput/calculationSetup/prodBasis/@fftcut'))
         ENDIF
      ENDIF

      DO itype = 1, ntype
         xpatha = xml%SpeciesPath(itype)//'/prodBasis'
         numberNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA)))
         IF (numberNodes == 1) THEN
            this%lcutm1(iType) = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@lcutm'))
            this%lcutwf(iType) = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@lcutwf'))
            xPathA = xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@select')
            this%select1(1, iType) = NINT(evaluateFirst(xPathA))
            this%select1(2, iType) = NINT(evaluateFirst(xPathA))
            this%select1(3, iType) = NINT(evaluateFirst(xPathA))
            this%select1(4, iType) = NINT(evaluateFirst(xPathA))
         ELSE
            this%lcutm1(iType) = -1
         ENDIF
      END DO

      xc_name = ''
      IF (xml%GetNumberOfNodes('/fleurInput/calculationSetup/xcFunctional') > 0) THEN
         xc_name = trim(xml%GetAttributeValue('/fleurInput/calculationSetup/xcFunctional/@name'))
      ELSE
         xc_name = trim(xml%GetAttributeValue('/fleurInput/xcFunctional/@name'))
      END IF
      if (trim(xc_name) == "pbe0") then
         this%l_hybrid = .True.
      else
         this%l_hybrid = .False.
      endif
   END SUBROUTINE read_xml_hybinp

   SUBROUTINE init_hybinp(self, atoms, cell, input,   sym, xcpot)
      USE m_dwigner
      use m_types_xcpot
      use m_types_sym
      use m_types_atoms
       
      use m_types_input
      use m_types_cell

      implicit none
      class(t_hybinp), intent(inout) :: self
      type(t_atoms), intent(in)      :: atoms
      type(t_cell), intent(in)       :: cell
      type(t_input), intent(in)      :: input
       
      type(t_sym), intent(in)        :: sym
      class(t_xcpot), intent(in)     :: xcpot

      integer :: isym, iisym, l, m2, m1

      IF (xcpot%is_hybrid() .OR. input%l_rdmft) THEN
         IF (input%film ) THEN
            CALL juDFT_error("2D film and 1D calculations not implemented for HF/EXX/PBE0/HSE", &
                             calledby="fleur", hint="Use a supercell or a different functional")
         END IF

         !             IF( ANY( atoms%l_geo  ) )&
         !                  &     CALL juDFT_error("Forces not implemented for HF/PBE0/HSE ",&
         !                  &                    calledby ="fleur")

         CALL self%gen_map(atoms, sym)

         ! calculate d_wgn
         ALLOCATE (self%d_wgn2(-atoms%lmaxd:atoms%lmaxd, -atoms%lmaxd:atoms%lmaxd, 0:atoms%lmaxd, sym%nsym))
         CALL d_wigner(sym%nop, sym%mrot, cell%bmat, atoms%lmaxd, self%d_wgn2(:, :, 1:, :sym%nop))
         self%d_wgn2(:, :, 0, :) = 1

         DO isym = sym%nop + 1, sym%nsym
            iisym = isym - sym%nop
            DO l = 0, atoms%lmaxd
               DO m2 = -l, l
                  DO m1 = -l, -1
                     self%d_wgn2(m1, m2, l, isym)  = self%d_wgn2(-m1, m2, l, iisym)*(-1)**m1
                     self%d_wgn2(-m1, m2, l, isym) = self%d_wgn2( m1, m2, l, iisym)*(-1)**m1
                  END DO
                  self%d_wgn2(0, m2, l, isym) = self%d_wgn2(0, m2, l, iisym)
               END DO
            END DO
         END DO
      ELSE
         ALLOCATE (self%map(0, 0), self%tvec(0, 0, 0), self%d_wgn2(0, 0, 0, 0))
      ENDIF
   END SUBROUTINE init_hybinp

   SUBROUTINE gen_map_hybinp(hybinp, atoms, sym)
      use m_types_atoms
      use m_types_sym
       
      USE m_juDFT
      use m_map_to_unit
      IMPLICIT NONE
      CLASS(t_hybinp), INTENT(INOUT) :: hybinp
      TYPE(t_atoms), INTENT(IN)      :: atoms
      TYPE(t_sym), INTENT(IN)        :: sym
      ! private scalars
      INTEGER                           :: iatom, first_eq_atom, itype, ieq, isym, iisym, ieq1
      INTEGER                           :: ratom, ok
      ! private arrays
      REAL                              :: rtaual(3), min_val, min_loc, curr_val

      ALLOCATE (hybinp%map(atoms%nat, sym%nsym), stat=ok)
      IF (ok /= 0) call judft_error('gen_map: error during allocation of map')

      ALLOCATE (hybinp%tvec(3, atoms%nat, sym%nsym), stat=ok)
      IF (ok /= 0) call judft_error('gen_map: error during allocation of tvec')

      iatom = 0
      first_eq_atom = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1
            DO isym = 1, sym%nsym

               IF (isym <= sym%nop) THEN
                  iisym = isym
               ELSE
                  iisym = isym - sym%nop
               END IF

               rtaual(:) = matmul(sym%mrot(:, :, iisym), atoms%taual(:, iatom)) + sym%tau(:, iisym)

               ratom = 0
               min_val = 1e33
               min_loc = -1
               DO ieq1 = 1, atoms%neq(itype)
                  curr_val = norm2(map_to_unit(rtaual - atoms%taual(:, first_eq_atom + ieq1)))
                  if(min_val > curr_val ) then 
                     min_loc = ieq1 
                     min_val = curr_val
                  endif 
               END DO
               if(min_val > 1e-6) call judft_error('eigen_hf: ratom not found' // new_line("A") //&
                                                   "iatom = " // int2str(iatom) // new_line("A") //&
                                                   "iisym = " // int2str(iisym) // new_line("A") )
               ratom = first_eq_atom + min_loc
               hybinp%map(iatom, isym) = ratom
               hybinp%tvec(:, iatom, isym) = nint(rtaual - atoms%taual(:, ratom))

            END DO
         END DO
         first_eq_atom = first_eq_atom + atoms%neq(itype)
      END DO

   END SUBROUTINE gen_map_hybinp
END MODULE m_types_hybinp
