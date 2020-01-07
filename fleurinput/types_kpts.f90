!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_kpts
   USE m_judft
   USE m_types_fleurinput_base
   IMPLICIT NONE
   PRIVATE
   TYPE, EXTENDS(t_fleurinput_base):: t_kpts
      character(len=20)     :: name = "default"
      INTEGER               :: nkpt = 0
      INTEGER               :: ntet = 0
      LOGICAL               :: l_gamma = .FALSE.
      !(3,nkpt) k-vectors internal units
      REAL, ALLOCATABLE      :: bk(:, :)
      !(nkpts) weights
      REAL, ALLOCATABLE      :: wtkpt(:)
      INTEGER               :: nkptf = 0   !<k-vectors in full BZ
      REAL, ALLOCATABLE   :: bkf(:, :)
      INTEGER, ALLOCATABLE   :: bkp(:)
      INTEGER, ALLOCATABLE   :: bksym(:)
      INTEGER                       :: numSpecialPoints = 0
      INTEGER, ALLOCATABLE          :: specialPointIndices(:)
      CHARACTER(LEN=50), ALLOCATABLE :: specialPointNames(:)
      REAL, ALLOCATABLE           :: specialPoints(:, :)
      INTEGER, ALLOCATABLE           :: ntetra(:, :)
      REAL, ALLOCATABLE           :: voltet(:)
      REAL, ALLOCATABLE           :: sc_list(:, :) !list for all information about folding of bandstructure (need for unfoldBandKPTS)((k(x,y,z),K(x,y,z),m(g1,g2,g3)),(nkpt),k_original(x,y,z))
   CONTAINS
      PROCEDURE :: add_special_line
      PROCEDURE :: print_xml
      PROCEDURE :: read_xml => read_xml_kpts
      PROCEDURE :: mpi_bc => mpi_bc_kpts
      procedure :: get_nk => kpts_get_nk
      procedure :: to_first_bz => kpts_to_first_bz
      procedure :: is_kpt => kpts_is_kpt
      procedure :: init => init_kpts
   ENDTYPE t_kpts

   PUBLIC :: t_kpts
CONTAINS

   function kpts_get_nk(kpts, kpoint) result(ret_idx)
      ! get the index of a kpoint
      implicit NONE
      class(t_kpts), intent(in)    :: kpts
      real, intent(in)            :: kpoint(3)
      integer                     :: idx, ret_idx

      DO idx = 1, kpts%nkptf
         IF (all(abs(kpoint - kpts%bkf(:, idx)) < 1E-06)) THEN
            ret_idx = idx
            return
         END IF
      END DO
      ret_idx = 0
   end function kpts_get_nk

   function kpts_to_first_bz(kpts, kpoint) result(out_point)
      implicit NONE
      class(t_kpts), intent(in)  :: kpts
      real, intent(in)           :: kpoint(3)
      real                       :: out_point(3)

      out_point = kpoint - floor(kpoint)
   end function kpts_to_first_bz

   function kpts_is_kpt(kpts, kpoint) result(is_kpt)
      implicit none
      class(t_kpts), intent(in)  :: kpts
      real, intent(in)           :: kpoint(3)
      logical                    :: is_kpt

      is_kpt = kpts%get_nk(kpoint) > 0
   end function kpts_is_kpt
   SUBROUTINE mpi_bc_kpts(this, mpi_comm, irank)
      USE m_mpi_bc_tool
      CLASS(t_kpts), INTENT(INOUT)::this
      INTEGER, INTENT(IN):: mpi_comm
      INTEGER, INTENT(IN), OPTIONAL::irank
      INTEGER ::rank
      IF (PRESENT(irank)) THEN
         rank = irank
      ELSE
         rank = 0
      END IF

      CALL mpi_bc(this%nkpt, rank, mpi_comm)
      CALL mpi_bc(this%ntet, rank, mpi_comm)
      CALL mpi_bc(this%l_gamma, rank, mpi_comm)
      CALL mpi_bc(this%bk, rank, mpi_comm)
      CALL mpi_bc(this%wtkpt, rank, mpi_comm)
      CALL mpi_bc(this%nkptf, rank, mpi_comm)
      CALL mpi_bc(this%bkf, rank, mpi_comm)
      CALL mpi_bc(this%bkp, rank, mpi_comm)
      CALL mpi_bc(this%bksym, rank, mpi_comm)
      CALL mpi_bc(this%numSpecialPoints, rank, mpi_comm)
      CALL mpi_bc(this%specialPointIndices, rank, mpi_comm)
      CALL mpi_bc(this%specialPoints, rank, mpi_comm)
      CALL mpi_bc(this%ntetra, rank, mpi_comm)
      CALL mpi_bc(this%voltet, rank, mpi_comm)
      CALL mpi_bc(this%sc_list, rank, mpi_comm)
   END SUBROUTINE mpi_bc_kpts

   SUBROUTINE read_xml_kpts(this, xml)
      USE m_types_xml
      USE m_calculator
      CLASS(t_kpts), INTENT(inout):: this
      TYPE(t_xml),INTENT(INOUT) ::xml

      INTEGER:: number_sets, n
      CHARACTER(len=200)::str, path, path2

      number_sets = xml%GetNumberOfNodes('/fleurInput/calculationSetup/bzIntegration/kPointList')

      DO n = number_sets, 1, -1
         WRITE (path, "(a,i0,a)") '/fleurInput/calculationSetup/bzIntegration/kPointList[', n, ']'
         IF (TRIM(ADJUSTL(this%name)) == xml%GetAttributeValue(TRIM(path)//'/@name')) EXIT
      enddo
      IF (n == 0) CALL judft_error(("No kpoints named:"//TRIM(this%name)//" found"))
      this%nkpt = evaluateFirstOnly(xml%GetAttributeValue(TRIM(path)//'/@count'))

      this%numSpecialPoints = xml%GetNumberOfNodes(TRIM(path)//"/specialPoint")

      IF (this%numSpecialPoints > 0) THEN
         If (this%numSpecialPoints < 2) CALL judft_error(("Giving less than two sepcial points make no sense:"//TRIM(this%name)))
         ALLOCATE (this%specialPoints(3, this%numSpecialPoints))
         ALLOCATE (this%specialPointNames(this%numSpecialPoints))
         DO n = 1, this%numSpecialPoints
            WRITE (path2, "(a,a,i0,a)") trim(path), "/specialPoint[", n, "]"
            this%specialPointNames(n) = xml%getAttributeValue(path2//"/@name")
            str = xml%getAttributeValue(path2)
            this%specialPoints(1, n) = evaluatefirst(str)
            this%specialPoints(2, n) = evaluatefirst(str)
            this%specialPoints(3, n) = evaluatefirst(str)
         ENDDO
      ELSE
         n = xml%GetNumberOfNodes(TRIM(path)//'/kPoint')
         IF (n .NE. this%nkpt) CALL judft_error(("Inconsistent number of k-points in:"//this%name))
         ALLOCATE (this%bk(3, this%nkpt))
         ALLOCATE (this%wtkpt(this%nkpt))
         DO n = 1, this%nkpt
            WRITE (path2, "(a,a,i0,a)") trim(adjustl(path)), "/kPoint[", n, "]"
            this%wtkpt(n) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(path2)//'/@weight'))
            str = xml%getAttributeValue(path2)
            this%bk(1, n) = evaluatefirst(str)
            this%bk(2, n) = evaluatefirst(str)
            this%bk(3, n) = evaluatefirst(str)
         ENDDO
         this%ntet = xml%GetNumberOfNodes(TRIM(path)//'/tetraeder/tet')
         ALLOCATE (this%voltet(this%ntet), this%ntetra(4, this%ntet))
         DO n = 1, this%ntet
            WRITE (path2, "(a,a,i0,a)") trim(path), "/tetraeder/tet[", n, "]"
            this%voltet(n) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(path2)//'/@vol'))
            str = xml%getAttributeValue(path2)
            READ (str, *) this%ntetra(:, n)
         ENDDO
      END IF
      this%wtkpt = this%wtkpt/sum(this%wtkpt) !Normalize k-point weight
   END SUBROUTINE read_xml_kpts

   SUBROUTINE print_xml(kpts, fh, filename)
      CLASS(t_kpts), INTENT(in):: kpts
      INTEGER, INTENT(in)         :: fh
      CHARACTER(len=*), INTENT(in), OPTIONAL::filename

      INTEGER :: n
      LOGICAL :: l_exist

      IF (PRESENT(filename)) THEN
         INQUIRE (file=filename, exist=l_exist)
         IF (l_exist) THEN
            OPEN (fh, file=filename, action="write", position="append")
         ELSE
            OPEN (fh, file=filename, action="write")
         END IF
      ENDIF

205   FORMAT('         <kPointList name="', a, '" count="', i0, '">')
      WRITE (fh, 205) adjustl(trim(kpts%name)), kpts%nkpt
      IF (kpts%numSpecialPoints < 2) THEN
         DO n = 1, kpts%nkpt
206         FORMAT('            <kPoint weight="', f12.6, '">', f12.6, ' ', f12.6, ' ', f12.6, '</kPoint>')
            WRITE (fh, 206) kpts%wtkpt(n), kpts%bk(:, n)
         END DO
         IF (kpts%ntet > 0) THEN
            WRITE (fh, 207) kpts%ntet
207         FORMAT('            <tetraeder ntet="', i0, '">')
            DO n = 1, kpts%ntet
208            FORMAT('          <tet vol="', f12.6, '">', i0, ' ', i0, ' ', i0, ' ', i0, '</tet>')
               WRITE (fh, 208) kpts%voltet(n), kpts%ntetra(:, n)
            END DO
            WRITE (fh, '(a)') '            </tetraeder>'
         ELSE
            DO n = 1, kpts%numSpecialPoints
               WRITE (fh, 209) TRIM(ADJUSTL(kpts%specialPointNames(n))), kpts%specialPoints(:, n)
209            FORMAT('            <specialPoint name="', a, '">', f10.6, ' ', f10.6, ' ', f10.6, '</specialPoint>')
            END DO
         END IF
      END IF
      WRITE (fh, '(a)') ('         </kPointList>')
      IF (PRESENT(filename)) CLOSE (fh)
   END SUBROUTINE print_xml

   SUBROUTINE add_special_line(kpts, point, name)
      CLASS(t_kpts), INTENT(inout):: kpts
      CHARACTER(len=*), INTENT(in):: name
      REAL, INTENT(in)            :: point(3)

      CHARACTER(len=50), ALLOCATABLE:: names(:)
      REAL, ALLOCATABLE             :: points(:, :)

      IF (kpts%numspecialpoints > 0) THEN
         CALL MOVE_ALLOC(kpts%specialPointNames, names)
         CALL MOVE_ALLOC(kpts%specialPoints, points)
         DEALLOCATE (kpts%specialpointindices)
         ALLOCATE (kpts%specialPoints(3, SIZE(names) + 1))
         ALLOCATE (kpts%specialPointIndices(SIZE(names) + 1))
         ALLOCATE (kpts%specialPointNames(SIZE(names) + 1))
         kpts%specialPoints(:, :SIZE(names)) = points
         kpts%specialPointNames(:SIZE(names)) = names
      ELSE
         ALLOCATE (kpts%specialPoints(3, 1))
         ALLOCATE (kpts%specialPointIndices(1))
         ALLOCATE (kpts%specialPointNames(1))
      ENDIF
      kpts%numspecialpoints = kpts%numspecialpoints + 1
      kpts%specialPoints(:, kpts%numspecialpoints) = point
      kpts%specialPointNames(kpts%numspecialpoints) = name
   END SUBROUTINE add_special_line

   SUBROUTINE init_kpts(self, cell, sym, film)
      USE m_types_cell
      USE m_types_sym
      CLASS(t_kpts), INTENT(inout):: self
      TYPE(t_cell), INTENT(IN)    :: cell
      TYPE(t_sym), INTENT(IN)     :: sym
      LOGICAL, INTENT(IN)         :: film

      INTEGER :: n

      DO n = 1, self%nkpt
         self%l_gamma = self%l_gamma .OR. ALL(ABS(self%bk(:, n)) < 1E-9)
      ENDDO
      IF (self%nkptf == 0) CALL gen_bz(self, sym)
   END SUBROUTINE init_kpts

   SUBROUTINE gen_bz(kpts, sym)
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! gen_bz generates the (whole) Brillouin zone from the          !
      ! (irreducible) k-points given                                  !
      !                                     M.Betzinger (09/07)       !
      !                        Refactored in 2017 by G.M.             !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !     bk     ::    irreducible k-points
      !     nkpt   ::    number of irr. k-points
      !     bkf    ::    all k-points
      !     nkptf  ::    number of all k-points
      !     bkp    ::    k-point parent
      !     bksym  ::    symmetry operation, that connects the parent
      !                  k-point with the current one
      USE m_juDFT
      USE m_types_sym
      TYPE(t_kpts), INTENT(INOUT) :: kpts
      TYPE(t_sym), INTENT(IN)     :: sym
      !  - local scalars -
      INTEGER                 ::  ic, iop, ikpt, ikpt1
      LOGICAL                 ::  l_found

      !  - local arrays -
      INTEGER, ALLOCATABLE     ::  iarr(:)
      REAL                    ::  rrot(3, 3, 2*sym%nop), rotkpt(3)
      REAL, ALLOCATABLE        ::  rarr1(:, :)

      INTEGER:: nsym, ID_mat(3, 3)

      nsym = sym%nop
      IF (.NOT. sym%invs) nsym = 2*sym%nop

      ALLOCATE (kpts%bkf(3, nsym*kpts%nkpt))
      ALLOCATE (kpts%bkp(nsym*kpts%nkpt))
      ALLOCATE (kpts%bksym(nsym*kpts%nkpt))

      ! Generate symmetry operations in reciprocal space
      DO iop = 1, nsym
         IF (iop .LE. sym%nop) THEN
            rrot(:, :, iop) = TRANSPOSE(sym%mrot(:, :, sym%invtab(iop)))
         ELSE
            rrot(:, :, iop) = -rrot(:, :, iop - sym%nop)
         END IF
      END DO

      !Add existing vectors to list of full vectors
      id_mat = 0
      ID_mat(1, 1) = 1; ID_mat(2, 2) = 1; ID_mat(3, 3) = 1
      IF (ANY(sym%mrot(:, :, 1) .NE. ID_mat)) CALL judft_error("Identity must be first symmetry operation", calledby="gen_bz")

      ic = 0
      DO iop = 1, nsym
         DO ikpt = 1, kpts%nkpt
            rotkpt = MATMUL(rrot(:, :, iop), kpts%bk(:, ikpt))
            !transform back into 1st-BZ (Do not use nint to deal properly with inaccuracies)
            do while (any(rotkpt < -0.5))
               where (rotkpt < -0.5) rotkpt = rotkpt + 1.0
            enddo
            do while (any(rotkpt > 0.5 + 1E-8))
               where (rotkpt > 0.5 + 1E-8) rotkpt = rotkpt - 1.0
            enddo
            DO ikpt1 = 1, ic
               IF (MAXVAL(ABS(kpts%bkf(:, ikpt1) - rotkpt)) .LE. 1e-07) EXIT
            END DO

            IF (ikpt1 > ic) THEN !new point
               ic = ic + 1
               kpts%bkf(:, ic) = rotkpt
               kpts%bkp(ic) = ikpt
               kpts%bksym(ic) = iop
            END IF
         END DO
      END DO

      kpts%nkptf = ic
      ! Reallocate bkf, bkp, bksym
      ALLOCATE (iarr(kpts%nkptf))
      iarr = kpts%bkp(:kpts%nkptf)
      DEALLOCATE (kpts%bkp)
      ALLOCATE (kpts%bkp(kpts%nkptf))
      kpts%bkp = iarr
      iarr = kpts%bksym(:kpts%nkptf)
      DEALLOCATE (kpts%bksym)
      ALLOCATE (kpts%bksym(kpts%nkptf))
      kpts%bksym = iarr
      DEALLOCATE (iarr)
      ALLOCATE (rarr1(3, kpts%nkptf))
      rarr1 = kpts%bkf(:, :kpts%nkptf)
      DEALLOCATE (kpts%bkf)
      ALLOCATE (kpts%bkf(3, kpts%nkptf))
      kpts%bkf = rarr1
      DEALLOCATE (rarr1)

   END SUBROUTINE gen_bz
END MODULE m_types_kpts
