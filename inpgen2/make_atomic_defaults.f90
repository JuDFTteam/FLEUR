!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_make_atomic_defaults
  USE m_judft
  IMPLICIT NONE

CONTAINS
  SUBROUTINE make_atomic_defaults(input,vacuum,cell,oneD,atoms,enpara)
    USE m_check_mt_radii
    USE m_atompar
    USE m_types_atoms
    USE m_types_input
    USE m_types_vacuum
    USE m_types_cell
    USE m_types_oneD
    USE m_constants
    USE m_types_enpara
    TYPE(t_atoms),INTENT(INOUT)   :: atoms
    TYPE(t_enpara),INTENT(OUT)    :: enpara

      TYPE(t_input),INTENT(IN)    :: input
      TYPE(t_vacuum),INTENT(IN)   :: vacuum
      TYPE(t_cell),INTENT(IN)     :: cell
      TYPE(t_oneD),INTENT(IN)     :: oneD

      INTEGER :: i,l,id,n,nn
      INTEGER :: element_species(120)

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

      atoms%lapw_l=0
      atoms%speciesname=""

      atoms%nz(:) = NINT(atoms%zatom(:))
      atoms%rmt(:) = 999.9
      atoms%l_geo(:) = .TRUE.;
      atoms%flipSpinPhi=0.0
      atoms%flipSpinTheta=0.0
      atoms%flipSpinScale=.FALSE.

      atoms%lda_u%l = -1 ; atoms%relax(1:2,:) = 1 ; atoms%relax(:,:) = 1

      !Determine MT-radii
      CALL check_mt_radii(atoms,input,vacuum,cell,oneD,.false.,atoms%rmt)
      !rounding
      atoms%rmt(:)  = real(NINT(atoms%rmt(:)  * 100 ) / 100.)

      !Now set the defaults
      DO n=1,atoms%ntype
         id=NINT((atoms%zatom(n)-atoms%nz(n))*100)
         IF (id>0) THEN
            ap(n)=find_atompar(atoms%nz(n),atoms%rmt(n),id)
            !This specific atom also has a rmt given?
            IF (ap(n)%id==id.AND.ap(n)%rmt>0.0) atoms%rmt(n)=ap(n)%rmt
         ELSE
            ap(n)=find_atompar(atoms%nz(n),atoms%rmt(n))
         ENDIF
         CALL ap(n)%add_defaults()
         atoms%jri(n)=ap(n)%jri
         atoms%dx(n)=ap(n)%dx
         atoms%lmax(n)=ap(n)%lmax
         atoms%lnonsph(n)=ap(n)%lnonsph
         !atoms%bmu(n))=ap(n)%bmu
         !local orbitals
         atoms%nlo(n)=len_TRIM(ap(n)%lo)/2
         DO i=1,atoms%nlo(n)
            DO l = 0, 3
               !Setting of llo will be redone below
               IF (ap(n)%lo(2*i:2*i) == lotype(l)) atoms%llo(i,n) = l
            ENDDO
         ENDDO
         CALL atoms%econf(n)%init(ap(n)%econfig)
         if (abs(ap(n)%bmu)>1E-8) call atoms%econf(n)%set_initial_moment(ap(n)%bmu)
         !atoms%ncst(n)=econfig_count_core(econfig)


         atoms%lda_u(n)%l=-1
         atoms%l_geo(n)=.FALSE.
         atoms%relax(:,n)=1

         ! rounding
         atoms%dx(:)   = REAL(NINT(atoms%dx(:)   * 1000) / 1000.)
         !Generate species-names
         DO nn=1,n
            if (atoms%same_species(n,nn)) atoms%speciesname(n)=atoms%speciesname(nn)
         enddo
         if (len_trim(atoms%speciesname(n))<1) THEN
            element_species(atoms%nz(n))=element_species(atoms%nz(n))+1
            write(atoms%speciesname(n),"(a,a,i0)") namat_const(atoms%nz(n)),"-",element_species(atoms%nz(n))
         endif

      END DO
      atoms%nlod=MAXVAL(atoms%nlo)
      atoms%lmaxd=MAXVAL(atoms%lmax)
      DEALLOCATE(atoms%llo)
      ALLOCATE(atoms%llo(atoms%nlod,atoms%ntype));atoms%llo=-1
      ALLOCATE(atoms%ulo_der(atoms%nlod,atoms%ntype))
      atoms%ulo_der=0

      CALL enpara%init(atoms%ntype,atoms%nlod,2,.TRUE.,atoms%nz)
      DO n=1,atoms%ntype
         DO i=1,atoms%nlo(n)
            DO l = 0, 3
               IF (ap(n)%lo(2*i:2*i) == lotype(l)) atoms%llo(i,n) = l
            ENDDO
         ENDDO
         CALL enpara%set_quantum_numbers(n,atoms,ap(n)%econfig,ap(n)%lo)
      END DO


    END SUBROUTINE make_atomic_defaults
  END MODULE m_make_atomic_defaults
