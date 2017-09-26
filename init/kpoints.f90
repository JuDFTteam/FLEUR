!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_kpoints
contains
  subroutine kpoints(oneD,jij,sym,cell,input,noco,banddos,kpts,l_kpts)
    USE m_juDFT
    USE m_types
    USE m_julia
    USE m_kptgen_hybrid
    USE m_od_kptsgen
  
    implicit none
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_sym),INTENT(IN)     :: sym
    TYPE(t_oneD),INTENT(IN)     :: oneD
    TYPE(t_Jij),INTENT(IN)      :: Jij
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_banddos),INTENT(IN)  :: banddos
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_kpts),INTENT(INOUT)  :: kpts
    LOGICAL,INTENT(IN)          :: l_kpts

    TYPE(t_sym) :: sym_hlp
    INTEGER     :: i
    REAL        :: recVecLength

    IF (kpts%specificationType.EQ.4) THEN
       DO i = 1, 3
          IF (kpts%kPointDensity(i).LE.0.0) THEN
             CALL juDFT_error('Error: Nonpositive kpointDensity provided', calledby = 'kpoints')
          END IF
          recVecLength = SQRT(cell%bmat(i,1)**2 + cell%bmat(i,2)**2 + cell%bmat(i,3)**2)
          kpts%nkpt3(i) = CEILING(kpts%kPointDensity(i) * recVecLength)
       END DO
       kpts%nkpt = kpts%nkpt3(1) * kpts%nkpt3(2) * kpts%nkpt3(3)
    END IF

    IF (.NOT.l_kpts) THEN
       IF (.NOT.oneD%odd%d1) THEN
          IF (jij%l_J) THEN
             sym_hlp=sym
             sym_hlp%nop=1
             sym_hlp%nop2=1
             CALL julia(sym_hlp,cell,input,noco,banddos,kpts,.FALSE.,.TRUE.)
          ELSE IF(kpts%l_gamma .and. banddos%ndir .eq. 0) THEN
             CALL kptgen_hybrid(kpts,sym%invs,noco%l_soc,sym%nop,sym%mrot,sym%tau)
          ELSE
             CALL julia(sym,cell,input,noco,banddos,kpts,.FALSE.,.TRUE.)
          END IF
       ELSE
          STOP 'Error: No kpoint set generation for 1D systems yet!'
          CALL od_kptsgen (kpts%nkpt)
       END IF
    END IF
    !Rescale weights and kpoints

    kpts%wtkpt(:) = kpts%wtkpt(:) / sum(kpts%wtkpt)
    kpts%bk(:,:) = kpts%bk(:,:) / kpts%posScale
    kpts%posScale = 1.0
    IF (kpts%nkpt3(3).EQ.0) kpts%nkpt3(3) = 1

  
  
end subroutine kpoints
end module m_kpoints
