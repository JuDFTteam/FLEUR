MODULE m_fleurinput_postprocess
  USE m_types_fleurinput
CONTAINS
  SUBROUTINE fleurinput_postprocess(Cell,Sym,Atoms,Input,Noco,Vacuum,&
    Banddos ,Xcpot,Kpts,gfinp)
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

    TYPE(t_cell),INTENT(INOUT)  ::cell
    TYPE(t_sym),INTENT(INOUT)   ::sym
    TYPE(t_atoms),INTENT(INOUT) ::atoms
    TYPE(t_input),INTENT(INOUT) ::input
    TYPE(t_noco),INTENT(INOUT)     ::noco
    TYPE(t_vacuum),INTENT(INOUT)::vacuum
    TYPE(t_banddos),INTENT(IN)  ::banddos
     
    CLASS(t_xcpot),ALLOCATABLE,INTENT(INOUT)::xcpot
    TYPE(t_kpts),INTENT(INOUT)     ::kpts
    TYPE(t_gfinp),INTENT(IN)    ::gfinp
    REAL    :: unfold(3,3)  !just for the unfolding
    REAL    :: tempPos(3), tempPosIntern(3), distance ! just for LDA+V
    INTEGER :: i, j
    call cell%init(DOT_PRODUCT(atoms%volmts(:),atoms%neq(:)))
    call atoms%init(cell)
    CALL sym%init(cell,input%film)
    call vacuum%init(sym)

    CALL make_sym(sym,cell,atoms,noco ,input,gfinp)
    !call make_xcpot(xcpot,atoms,input)
    CALL noco%init(atoms,input%ldauSpinoffd)

    call check_input_switches(banddos,vacuum,noco,atoms,input,sym,kpts)
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
            distance = SQRT(norm2(tempPos))
            WRITE(oUnit,'(a,i5,a,3i3,a,f15.8)') 'otherAtom= ', atoms%lda_v(i)%otherAtomIndices(j), ' shift= ', atoms%lda_v(i)%atomShifts(:,j), ' distance= ', distance
         END DO
      END DO
   END IF
!--------------------------------------------------------------------------

  END SUBROUTINE fleurinput_postprocess
END MODULE m_fleurinput_postprocess
