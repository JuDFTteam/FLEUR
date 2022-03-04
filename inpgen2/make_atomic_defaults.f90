!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_make_atomic_defaults
  USE m_judft
  IMPLICIT NONE

CONTAINS
   SUBROUTINE make_atomic_defaults(input,vacuum,profile,cell,oneD,atoms,enpara)
      USE m_check_mt_radii
      USE m_atompar
      USE m_types_atoms
      USE m_types_input
      USE m_types_vacuum
      USE m_types_cell
      USE m_types_oneD
      USE m_constants
      USE m_types_enpara
      USE m_types_profile

      TYPE(t_atoms),INTENT(INOUT)   :: atoms
      TYPE(t_enpara),INTENT(OUT)    :: enpara
      TYPE(t_profile),INTENT(IN)    :: profile
      TYPE(t_input),INTENT(IN)      :: input
      TYPE(t_vacuum),INTENT(IN)     :: vacuum
      TYPE(t_cell),INTENT(IN)       :: cell
      TYPE(t_oneD),INTENT(IN)       :: oneD

      INTEGER :: i,l,id,n,nn,qn,iLO
      INTEGER :: loLCutoff
      INTEGER :: element_species(120)
      INTEGER :: addLOs(atoms%ntype)
      INTEGER :: numLOs(10,0:3)
      INTEGER :: numlLO(0:3)
      INTEGER :: inequivalentNLO

      CHARACTER(len=1) :: lotype(0:3)=(/'s','p','d','f'/)
      TYPE(t_atompar):: ap(atoms%ntype)
      element_species=0

      ALLOCATE(atoms%nz(atoms%ntype))
      ALLOCATE(atoms%jri(atoms%ntype))
      ALLOCATE(atoms%dx(atoms%ntype))
      ALLOCATE(atoms%lmax(atoms%ntype))
      ALLOCATE(atoms%nlo(atoms%ntype))
      ALLOCATE(atoms%lnonsph(atoms%ntype))
      ALLOCATE(atoms%flipSpinPhi(atoms%ntype))
      ALLOCATE(atoms%flipspinTheta(atoms%ntype))
      ALLOCATE(atoms%flipSpinScale(atoms%ntype))
      ALLOCATE(atoms%l_geo(atoms%ntype))
      ALLOCATE(atoms%lda_u(atoms%ntype))
      ALLOCATE(atoms%econf(atoms%ntype))
      ALLOCATE(atoms%relax(3,atoms%ntype))
      ALLOCATE(atoms%rmt(atoms%ntype))
      ALLOCATE(atoms%speciesname(atoms%ntype))
      ALLOCATE(atoms%lapw_l(atoms%ntype))
      ALLOCATE(atoms%llo(99,atoms%ntype));atoms%llo=-1!will be redone later

      addLOs(:) = 0
      atoms%lapw_l=0
      atoms%speciesname=""

      atoms%nz(:) = floor(atoms%zatom(:))
      atoms%rmt(:) = 999.9
      atoms%l_geo(:) = .TRUE.
      atoms%flipSpinPhi=0.0
      atoms%flipSpinTheta=0.0
      atoms%flipSpinScale=.FALSE.

      atoms%lda_u%l = -1 ;  atoms%relax(:,:) = 1

      !Determine MT-radii
      CALL check_mt_radii(atoms,input,vacuum,cell,oneD,profile,.false.,atoms%rmt)

      IF(TRIM(ADJUSTL(profile%profileName)).NE."default") THEN
         atoms%rmt(:) = atoms%rmt(:) * profile%rmtFactor
      END IF

      !rounding
      atoms%rmt(:) = real(NINT(atoms%rmt(:)  * 100 ) / 100.)

      !Now set the defaults
      DO n=1,atoms%ntype
         id=NINT((atoms%zatom(n)-atoms%nz(n))*100)
         IF (id>0) THEN
            ap(n)=find_atompar(atoms%nz(n),atoms%rmt(n),profile,id)
            !This specific atom also has a rmt given?
!            IF (ap(n)%id==id.AND.ap(n)%rmt>0.0) atoms%rmt(n)=ap(n)%rmt
         ELSE
            ap(n)=find_atompar(atoms%nz(n),atoms%rmt(n),profile)   
         ENDIF
         IF (ap(n)%rmt>0.0) atoms%rmt(n)=ap(n)%rmt
         CALL ap(n)%add_defaults()
         atoms%speciesName(n)=ap(n)%desc
         atoms%jri(n)=ap(n)%jri
         atoms%dx(n)=ap(n)%dx
         atoms%lmax(n)=ap(n)%lmax
         atoms%lnonsph(n)=ap(n)%lnonsph
         !atoms%bmu(n))=ap(n)%bmu
         !local orbitals
         atoms%nlo(n)=len_TRIM(ap(n)%lo)/2
         loLCutoff = 4 ! l cutoff for additional LOs. This l quantum number is already excluded.
         IF (atoms%nz(n).EQ.1) loLCutoff = 2
         IF (atoms%nz(n).EQ.2) loLCutoff = 3
         IF ((atoms%nz(n).GE.5).AND.(atoms%nz(n).LE.9)) loLCutoff = 3
         numlLO(:) = 0
         DO i=1,atoms%nlo(n)
            DO l = 0, 3
               !Setting of llo will be redone below
               IF (ap(n)%lo(2*i:2*i) == lotype(l)) THEN
                  atoms%llo(i,n) = l
                  numlLO(l) = numlLO(l) + 1
               END IF
            ENDDO
         ENDDO
         inequivalentNLO = 0
         DO l = 0, 3
            IF(numlLO(l).NE.0) inequivalentNLO = inequivalentNLO + 1
         END DO
         IF ((INDEX(TRIM(ADJUSTL(profile%addLOSetup)),"addHELOs_noSC").NE.0).OR.(INDEX(TRIM(ADJUSTL(profile%addLOSetup)),"addHDLOs_noSC").NE.0)) THEN
            addLOs(n) = loLCutoff - inequivalentNLO
         END IF
         atoms%nlo(n) = atoms%nlo(n) + addLOs(n)
         CALL atoms%econf(n)%init(ap(n)%econfig)
         if (abs(ap(n)%bmu)>1E-8.and.input%jspins.ne.1) call atoms%econf(n)%set_initial_moment(ap(n)%bmu)
         !atoms%ncst(n)=econfig_count_core(econfig)
     ! rounding
         atoms%dx(n)   = REAL(NINT(atoms%dx(n)   * 1000) / 1000.)
         !Generate species-names
         DO nn=1,n
            if (atoms%same_species(n,nn)) atoms%speciesname(n)=atoms%speciesname(nn)
         enddo
         if (len_trim(atoms%speciesname(n))<1) THEN
            element_species(atoms%nz(n))=element_species(atoms%nz(n))+1
            write(atoms%speciesname(n),"(a,a,i0)") namat_const(atoms%nz(n)),"-",element_species(atoms%nz(n))
         endif

         IF(TRIM(ADJUSTL(profile%profileName)).NE."default") THEN
            atoms%lmax(n) = NINT(profile%kmax*profile%lmaxFactor*atoms%rmt(n))
            atoms%lnonsph(n) = MIN( MAX(atoms%lmax(n)-2,3) , 8 )
         END IF

      END DO
      atoms%nlod=MAXVAL(atoms%nlo)
      atoms%lmaxd=MAXVAL(atoms%lmax)
      DEALLOCATE(atoms%llo)
      ALLOCATE(atoms%llo(atoms%nlod,atoms%ntype));atoms%llo=-1
      ALLOCATE(atoms%ulo_der(atoms%nlod,atoms%ntype))
      atoms%ulo_der=0

      CALL enpara%init(atoms%ntype,atoms%nlod,2,.TRUE.,atoms%nz)
      DO n=1,atoms%ntype
         loLCutoff = 4 ! l cutoff for additional LOs. This l quantum number is already excluded.
         IF (atoms%nz(n).EQ.1) loLCutoff = 2
         IF (atoms%nz(n).EQ.2) loLCutoff = 3
         IF ((atoms%nz(n).GE.5).AND.(atoms%nz(n).LE.9)) loLCutoff = 3 
         DO i=1,atoms%nlo(n) - addLOs(n)
            DO l = 0, 3
               IF (ap(n)%lo(2*i:2*i) == lotype(l)) atoms%llo(i,n) = l
            ENDDO
         ENDDO
         CALL enpara%set_quantum_numbers(n,atoms,ap(n)%econfig,ap(n)%lo)

         IF ((INDEX(TRIM(ADJUSTL(profile%addLOSetup)),"addHELOs_noSC").NE.0).OR.(INDEX(TRIM(ADJUSTL(profile%addLOSetup)),"addHDLOs_noSC").NE.0)) THEN
            i = 0
            DO l = 0, loLCutoff - 1
               IF(ANY(atoms%llo(1:atoms%nlo(n)-addLOs(n),n).EQ.l)) CYCLE
               iLO = atoms%nlo(n) - addLOs(n) + 1 + i
               qn = -(enpara%qn_el(l,n,1)+1)
               IF(INDEX(TRIM(ADJUSTL(profile%addLOSetup)),"addHDLOs_noSC").NE.0) THEN
                  qn = enpara%qn_el(l,n,1)
                  atoms%ulo_der(iLO,n) = 2
               END IF
               enpara%qn_ello(iLO,n,:) = qn
               atoms%llo(iLO,n) = l
               i = i + 1
            END DO
         END IF

         numLOs(:,:) = 0
         DO iLO = 1,atoms%nlo(n) - addLOs(n)
            ! If the main quantum number of the LO is larger than that of the
            ! LAPW basis we make it a HELO type LO.
            l = atoms%llo(iLO,n)
            qn = enpara%qn_ello(iLO,n,1)
            IF (qn.GT.enpara%qn_el(l,n,1)) enpara%qn_ello(iLO,n,1) = -qn
            ! Adjust ulo_der for multiple LOs for the same qn, l.
            atoms%ulo_der(iLO,n) = numLOs(ABS(qn),l)
            numLOs(ABS(qn),l) = numLOs(ABS(qn),l) + 1
         END DO

      END DO


    END SUBROUTINE make_atomic_defaults
  END MODULE m_make_atomic_defaults
