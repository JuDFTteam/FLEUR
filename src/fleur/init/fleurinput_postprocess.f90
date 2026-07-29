!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_fleurinput_postprocess
  USE m_types_fleurinput
   implicit none
CONTAINS
  SUBROUTINE fleurinput_postprocess(Cell,Sym,Atoms,Input,Noco,Vacuum,&
    Banddos,hybinp ,Xcpot,Kpts,gfinp,wannierlib)
    USE m_juDFT
    USE m_types_fleurinput
    use m_make_sym
    USE m_chkmt
    !use m_make_xcpot
    use m_lapwdim
    use m_checks
    USE m_relaxio
    USE m_types_nococonv
    USE m_constants
    USE m_types_wannierlib
    IMPLICIT NONE

    TYPE(t_cell),INTENT(INOUT)  ::cell
    TYPE(t_sym),INTENT(INOUT)   ::sym
    TYPE(t_atoms),INTENT(INOUT) ::atoms
    TYPE(t_input),INTENT(INOUT) ::input
    TYPE(t_noco),INTENT(INOUT)     ::noco
    TYPE(t_vacuum),INTENT(INOUT)::vacuum
    TYPE(t_banddos),INTENT(IN)  ::banddos
    TYPE(t_hybinp),INTENT(IN)   :: hybinp 
    TYPE(t_wannierlib_wannierize),INTENT(INOUT) ::wannierlib

    CLASS(t_xcpot),ALLOCATABLE,INTENT(INOUT)::xcpot
    TYPE(t_kpts),INTENT(INOUT)     ::kpts
    TYPE(t_gfinp),INTENT(IN)    ::gfinp
    REAL    :: unfold(3,3)  !just for the unfolding
    REAL    :: tempPos(3), tempPosIntern(3), distance ! just for LDA+V
    INTEGER :: i, j, ilo, l
    call cell%init(DOT_PRODUCT(atoms%volmts(:),atoms%neq(:)))
    call atoms%init(cell)
    CALL sym%init(cell,input%film)
    call vacuum%init(sym)

    CALL make_sym(sym,cell,atoms,noco ,input,gfinp)
    !call make_xcpot(xcpot,atoms,input)
    CALL noco%init(atoms,input%ldauSpinoffd)

    IF (noco%l_noco.AND.input%jspins==1) THEN
       CALL judft_error("l_noco=T requires jspins=2",calledby='fleurinput_postprocess')
    END IF

    ! relLO is a FIRST-VARIATION spin-orbit feature only: it is kept in the basis exactly for
    ! noco+SOC (l_soc .AND. l_noco), where SOC is folded directly into the
    ! first-variational Hamiltonian. For every other run mode -- SOC-off, and 2nd-variation SOC
    ! (l_soc,.NOT.l_noco) -- the relLO is removed ENTIRELY by reducing nlo/nlotot. relLO entries
    ! are always at the tail of each type's LO list (see read_xml_atoms), so a simple count
    ! reduction drops them, and every consumer of nlo/nlotot (matrix construction, charge
    ! density, DOS, forces, MPI, radial-function/SOC-integral/energy-parameter generation, ...)
    ! is thereby uniformly blind to the relLO whenever it should be.
    IF (.NOT.(noco%l_soc.AND.noco%l_noco)) THEN
       atoms%nlo     = atoms%nlo     - atoms%nRelLO
       atoms%nlotot = 0
       DO i = 1, atoms%ntype
          DO j = 1, atoms%nlo(i)
             atoms%nlotot = atoms%nlotot + atoms%neq(i) * ( 2*atoms%llo(j,i) + 1 )
          END DO
       END DO
       ! keep lo1l/nlol (first LO index / LO count per l, used e.g. by force_a21_lo.f90) in
       ! sync with the reduced atoms%nlo, mirroring the computation in read_xml_atoms.
       DO i = 1, atoms%ntype
          atoms%lo1l(:,i) = 0
          atoms%nlol(:,i) = 0
          DO ilo = 1, atoms%nlo(i)
             l = atoms%llo(ilo,i)
             IF (ilo==1) THEN
                atoms%lo1l(l,i) = ilo
             ELSE IF (l/=atoms%llo(ilo-1,i)) THEN
                atoms%lo1l(l,i) = ilo
             END IF
             atoms%nlol(l,i) = atoms%nlol(l,i) + 1
          END DO
       END DO
    END IF

    call check_input_switches(banddos,vacuum,noco,atoms,input,sym,kpts,hybinp)
    ! Check muffin tin radii, only checking, dont use new parameters
    CALL chkmt(atoms,input,vacuum,cell ,.TRUE.)
    !adjust positions by displacements
    CALL apply_displacements(cell,input,vacuum ,sym,noco,atoms,gfinp)
!---------------band unfolding ---------------------
    IF (banddos%unfoldband) THEN
      write (*,*) 'input switch unfolding read'
      write (*,*) 'before', kpts%specialPoints(:,2)
      !unfold=banddos%unfoldTransMat
      !unfold(1,1)=banddos%unfoldTransMat(1,1)*banddos%s_cell_x
      !unfold(2,2)=banddos%unfoldTransMat(2,2)*banddos%s_cell_y
      !unfold(3,3)=banddos%unfoldTransMat(3,3)*banddos%s_cell_z
      !Do i= 1,kpts%nkpt
      !  kpts%bk(:,i)=matmul(unfold,kpts%bk(:,i))
      !END DO
      !Do i=1,size(kpts%specialPoints,2)
      !  write (*,*) 'before', kpts%specialPoints(:,i)
      !  kpts%specialPoints(:,i)=matmul(unfold,kpts%specialPoints(:,i))
      !  write (*,*) 'after', kpts%specialPoints(:,i)
      !END DO
      kpts%bk(1,:)=kpts%bk(1,:)*banddos%s_cell_x
      kpts%bk(2,:)=kpts%bk(2,:)*banddos%s_cell_y
      kpts%bk(3,:)=kpts%bk(3,:)*banddos%s_cell_z
      kpts%specialPoints(1,:)=kpts%specialPoints(1,:)*banddos%s_cell_x
      kpts%specialPoints(2,:)=kpts%specialPoints(2,:)*banddos%s_cell_y
      kpts%specialPoints(3,:)=kpts%specialPoints(3,:)*banddos%s_cell_z
      write (*,*) 'after', kpts%specialPoints(:,2)
    END IF
!--------------------------------------------------    

! Temporary output for LDA+V (may be put into an own routine or be deleted)
   IF (atoms%n_v.GT.0) THEN
      WRITE(oUnit,'(a)') 'LDA+V region parameter + distance output:'
      DO i = 1, atoms%n_v
         WRITE(oUnit,'(a,i5,a,i2,a,i2,a,f15.8)') 'refAtom= ', atoms%lda_v(i)%atomIndex, ' refAtomL= ', atoms%lda_v(i)%thisAtomL, ' otherAtomL= ', atoms%lda_v(i)%otherAtomL, ' V= ', atoms%lda_v(i)%V
         DO j = 1, SIZE(atoms%lda_v(i)%otheratomIndices)
            tempPosIntern(:) = atoms%taual(:,atoms%lda_v(i)%otherAtomIndices(j)) + atoms%lda_v(i)%atomShifts(:,j)
            tempPos(:) = MATMUL(cell%amat,tempPosIntern(:))
            tempPos(:) = tempPos(:) - atoms%pos(:,atoms%lda_v(i)%atomIndex)
            distance = norm2(tempPos)
            WRITE(oUnit,'(a,i5,a,3i3,a,f15.8)') 'otherAtom= ', atoms%lda_v(i)%otherAtomIndices(j), ' shift= ', atoms%lda_v(i)%atomShifts(:,j), ' distance= ', distance
         END DO
      END DO
   END IF
!--------------------------------------------------------------------------
  CALL wannierlib%init(atoms, noco)
  END SUBROUTINE fleurinput_postprocess
END MODULE m_fleurinput_postprocess
