!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_make_atomic_defaults
  USE m_judft
  IMPLICIT NONE

CONTAINS
  SUBROUTINE make_atomic_defaults(atoms,sym)
    



    CHARACTER(len=1) :: lotype(0:3)=(/'s','p','d','f'/)

     ALLOCATE (atomTypeSpecies(atoms%ntype))
      ALLOCATE (speciesRepAtomType(atoms%nat))
      ALLOCATE (atoms%speciesName(atoms%nat))
      elementNumSpecies = 0
      numSpecies = 0
      speciesRepAtomType = -1
      atomTypeSpecies = -1
      atoms%speciesName = ''
      DO i = 1, atoms%nat
         newSpecies = .TRUE.
         DO j = 1, i-1
            IF(atomid(i).EQ.atomid(j)) THEN
               newSpecies = .FALSE.
               atomTypeSpecies(natype(i)) = atomTypeSpecies(natype(j))
               EXIT
            END IF
         END DO
         IF(newSpecies) THEN
            numSpecies = numSpecies + 1
            speciesRepAtomType(numSpecies) = natype(i)
            atomTypeSpecies(natype(i)) = numSpecies
            element = nint(atoms%zatom(natype(i)))
            elementNumSpecies(element) = elementNumSpecies(element) + 1
            tempNumberString = ''
            WRITE(tempNumberString,'(i0)') elementNumSpecies(element)
            atoms%speciesName(numSpecies) = &
               TRIM(ADJUSTL(namat_const(element))) // '-' // TRIM(ADJUSTL(tempNumberString))
         END IF
      END DO


         DATA xmlCoreRefOccs /2,2,2,4,2,2,4,2,4,6,2,4,2,4,6,2,4,2,6,8,4,&
     &                     6,2,4,2,6,8,4,6/
      xmlElectronStates = noState_const
      xmlPrintCoreStates = .FALSE.
      xmlCoreOccs = 0.0

      atoms%nlod=9  ! This fixed dimensioning might have to be made more dynamical!
      ALLOCATE(atoms%nz(atoms%ntype))
      ALLOCATE(atoms%jri(atoms%ntype))
      ALLOCATE(atoms%dx(atoms%ntype))
      ALLOCATE(atoms%lmax(atoms%ntype))
      ALLOCATE(atoms%nlo(atoms%ntype))
      ALLOCATE(atoms%llo(atoms%nlod,atoms%ntype))
      ALLOCATE(atoms%ncst(atoms%ntype))
      ALLOCATE(atoms%lnonsph(atoms%ntype))
      ALLOCATE(atoms%nflip(atoms%ntype))
      ALLOCATE(atoms%l_geo(atoms%ntype))
      ALLOCATE(atoms%lda_u(atoms%ntype))
      ALLOCATE(atoms%bmu(atoms%ntype))
      ALLOCATE(atoms%relax(3,atoms%ntype))
      ALLOCATE(atoms%ulo_der(atoms%nlod,atoms%ntype))
      
      atoms%nz(:) = NINT(atoms%zatom(:))
      atoms%rmt(:) = 999.9
      atoms%ulo_der = 0
      atoms%l_geo(:) = .TRUE.; atoms%nflip(:) = 1
      atoms%lda_u%l = -1 ; atoms%relax(1:2,:) = 1 ; atoms%relax(:,:) = 1

      !Determine MT-radii
      CALL check_mt_radii(atoms,input,vacuum,cell,oneD,.false.,atoms%rmt)
      !rounding
      atoms%rmt(:)  = real(NINT(atoms%rmt(:)  * 100 ) / 100.)

      !Now set the defaults
      DO n=1,atoms%ntype
         id=NINT(atoms%zatom(n)-atoms%nz(n)*100)
         IF (id>0) THEN
            ap=find_atompar(atoms%nz(n),atoms%rmt(n),id)
            !This specific atom also has a rmt given?
            IF (ap%id==id.AND.ap%rmt>0.0) atoms%rmt(n)=ap%rmt
         ELSE
            ap=find_atompar(atoms%nz(n),atoms%rmt(n))
         ENDIF
         CALL ap%add_defaults()
         atoms%jri(n)=ap%jri
         atoms%dx(n)=ap%dx
         atoms%lmax(n)=ap%lmax
         atoms%lnonsph(n)=ap%lnonsph
         atoms%bmu(n))=ap%bmu
         !local orbitals
         atoms%nlo(n)len_TRIM(ap%lo)/2
         DO i=1,atoms%nlo
            DO l = 0, 3
               IF (ap%lo(2*i:2*) == lotype(l)) atoms%llo(i,n) = l         
            ENDDO
         ENDDO
         atoms%ulo_der(:,n)=.FALSE.
         
         econfig(n)=ap%econfig
         atoms%ncst(n)=econfig_count_core(econfig)
         
         
         atoms%nflip(n)=1
         atoms%lda_u(n)%l=-1
         atoms%l_geo(n)=.FALSE.
         atoms%relax(:,n)=1
         
         ! rounding
         atoms%dx(:)   = REAL(NINT(atoms%dx(:)   * 1000) / 1000.)
      END DO
      



    END SUBROUTINE make_atomic_defaults
  END MODULE m_make_atomic_defaults
