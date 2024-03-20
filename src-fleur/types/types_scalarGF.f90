MODULE m_types_scalarGF

   !This module contains a generic scalarproduct type, which is used to construct
   !the spherical averages for the intersite and l-offdiagonal elements of the GF

   USE m_constants
   USE m_juDFT
   USE m_types_atoms
   USE m_types_input

   IMPLICIT NONE

   PRIVATE

   TYPE t_scalarGF

      LOGICAL :: done

      REAL, ALLOCATABLE :: uun(:,:)
      REAL, ALLOCATABLE :: udn(:,:)
      REAL, ALLOCATABLE :: dun(:,:)
      REAL, ALLOCATABLE :: ddn(:,:)

      REAL, ALLOCATABLE :: uulon(:,:,:)
      REAL, ALLOCATABLE :: uloun(:,:,:)
      REAL, ALLOCATABLE :: dulon(:,:,:)
      REAL, ALLOCATABLE :: ulodn(:,:,:)

      REAL, ALLOCATABLE :: uloulopn(:,:,:,:)

   CONTAINS
      PROCEDURE, PASS :: init             => init_scalarGF
      PROCEDURE       :: addOffdScalarProduct
   END TYPE t_scalarGF

   PUBLIC t_scalarGF

   CONTAINS

   SUBROUTINE init_scalarGF(this,atoms,input)

      CLASS(t_scalarGF),         INTENT(INOUT) :: this
      TYPE(t_atoms),             INTENT(IN)    :: atoms
      TYPE(t_input),             INTENT(IN)    :: input

      this%done =.FALSE.
      IF(ALLOCATED(this%uun)) DEALLOCATE(this%uun)
      IF(ALLOCATED(this%udn)) DEALLOCATE(this%udn)
      IF(ALLOCATED(this%dun)) DEALLOCATE(this%dun)
      IF(ALLOCATED(this%ddn)) DEALLOCATE(this%ddn)

      IF(ALLOCATED(this%uulon)) DEALLOCATE(this%uulon)
      IF(ALLOCATED(this%uloun)) DEALLOCATE(this%uloun)
      IF(ALLOCATED(this%dulon)) DEALLOCATE(this%dulon)
      IF(ALLOCATED(this%ulodn)) DEALLOCATE(this%ulodn)

      IF(ALLOCATED(this%uloulopn)) DEALLOCATE(this%uloulopn)

      ALLOCATE(this%uun(input%jspins,input%jspins),source=0.0)
      ALLOCATE(this%udn(input%jspins,input%jspins),source=0.0)
      ALLOCATE(this%dun(input%jspins,input%jspins),source=0.0)
      ALLOCATE(this%ddn(input%jspins,input%jspins),source=0.0)

      ALLOCATE(this%uulon(atoms%nlod,input%jspins,input%jspins),source=0.0)
      ALLOCATE(this%uloun(atoms%nlod,input%jspins,input%jspins),source=0.0)
      ALLOCATE(this%dulon(atoms%nlod,input%jspins,input%jspins),source=0.0)
      ALLOCATE(this%ulodn(atoms%nlod,input%jspins,input%jspins),source=0.0)

      ALLOCATE(this%uloulopn(atoms%nlod,atoms%nlod,input%jspins,input%jspins),source=0.0)

   END SUBROUTINE init_scalarGF

   SUBROUTINE addOffdScalarProduct(this,l,lp,atomType,atomTypep,l_intersite,l_mperp,atoms,input,f,g,flo)

      USE m_intgr

      CLASS(t_scalarGF),   INTENT(INOUT) :: this
      INTEGER,             INTENT(IN)    :: l,lp
      INTEGER,             INTENT(IN)    :: atomType,atomTypep
      LOGICAL,             INTENT(IN)    :: l_mperp
      LOGICAL,             INTENT(IN)    :: l_intersite !Is there a non-zero interstitial phase
                                                        !(meaning we have to treat r and r' independently)
      TYPE(t_atoms),       INTENT(IN)    :: atoms
      TYPE(t_input),       INTENT(IN)    :: input
      REAL,                INTENT(IN)    :: f(:,:,0:,:,:)
      REAL,                INTENT(IN)    :: g(:,:,0:,:,:)
      REAL,                INTENT(IN)    :: flo(:,:,:,:,:)

      REAL :: uu_tmp(atoms%jmtd),uu_tmp2(atoms%jmtd)
      INTEGER :: j1,j2,j2_start,j2_end,ilo,ilop,jri

      IF(this%done) RETURN !Already calculated

      CALL timestart("Offdiagonal Scalar Product")
      DO j1 = 1, input%jspins
         j2_start = MERGE(1,j1,l_mperp)
         j2_end   = MERGE(input%jspins,j1,l_mperp)
         DO j2 = j2_start, j2_end
            IF(.NOT.l_intersite) THEN
               !Only l/=lp
               uu_tmp(:atoms%jri(atomType)) = f(:atoms%jri(atomType),1,lp,j1,atomTypep)*f(:atoms%jri(atomType),1,l,j2,atomType)&
                                            + f(:atoms%jri(atomType),2,lp,j1,atomTypep)*f(:atoms%jri(atomType),2,l,j2,atomType)
               CALL intgr3(uu_tmp,atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),&
                           this%uun(j1,j2))
               uu_tmp(:atoms%jri(atomType)) = f(:atoms%jri(atomType),1,lp,j1,atomTypep)*g(:atoms%jri(atomType),1,l,j2,atomType)&
                                            + f(:atoms%jri(atomType),2,lp,j1,atomTypep)*g(:atoms%jri(atomType),2,l,j2,atomType)
               CALL intgr3(uu_tmp,atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),&
                           this%udn(j1,j2))
               uu_tmp(:atoms%jri(atomType)) = g(:atoms%jri(atomType),1,lp,j1,atomTypep)*f(:atoms%jri(atomType),1,l,j2,atomType)&
                                            + g(:atoms%jri(atomType),2,lp,j1,atomTypep)*f(:atoms%jri(atomType),2,l,j2,atomType)
               CALL intgr3(uu_tmp,atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),&
                           this%dun(j1,j2))
               uu_tmp(:atoms%jri(atomType)) = g(:atoms%jri(atomType),1,lp,j1,atomTypep)*g(:atoms%jri(atomType),1,l,j2,atomType)&
                                            + g(:atoms%jri(atomType),2,lp,j1,atomTypep)*g(:atoms%jri(atomType),2,l,j2,atomType)
               CALL intgr3(uu_tmp,atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),&
                           this%ddn(j1,j2))

               DO ilo = 1, atoms%nlo(atomType)
                  IF(atoms%llo(ilo,atomType).NE.l) CYCLE
                  uu_tmp(:atoms%jri(atomType)) = f(:atoms%jri(atomType),1,lp,j1,atomTypep)*flo(:atoms%jri(atomType),1,ilo,j2,atomType)&
                                               + f(:atoms%jri(atomType),2,lp,j1,atomTypep)*flo(:atoms%jri(atomType),2,ilo,j2,atomType)
                  CALL intgr3(uu_tmp,atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),&
                              this%uulon(ilo,j1,j2))
                  uu_tmp(:atoms%jri(atomType)) = g(:atoms%jri(atomType),1,lp,j1,atomTypep)*flo(:atoms%jri(atomType),1,ilo,j2,atomType)&
                                               + g(:atoms%jri(atomType),2,lp,j1,atomTypep)*flo(:atoms%jri(atomType),2,ilo,j2,atomType)
                  CALL intgr3(uu_tmp,atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),&
                              this%dulon(ilo,j1,j2))
               ENDDO

               DO ilo = 1, atoms%nlo(atomType)
                  IF(atoms%llo(ilo,atomType).NE.lp) CYCLE
                  uu_tmp(:atoms%jri(atomType)) = flo(:atoms%jri(atomType),1,ilo,j1,atomTypep)*f(:atoms%jri(atomType),1,l,j2,atomType)&
                                               + flo(:atoms%jri(atomType),2,ilo,j1,atomTypep)*f(:atoms%jri(atomType),2,l,j2,atomType)
                  CALL intgr3(uu_tmp,atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),&
                              this%uloun(ilo,j1,j2))
                  uu_tmp(:atoms%jri(atomType)) = flo(:atoms%jri(atomType),1,ilo,j1,atomTypep)*g(:atoms%jri(atomType),1,l,j2,atomType)&
                                               + flo(:atoms%jri(atomType),2,ilo,j1,atomTypep)*g(:atoms%jri(atomType),2,l,j2,atomType)
                  CALL intgr3(uu_tmp,atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),&
                              this%ulodn(ilo,j1,j2))
               ENDDO

               DO ilo = 1, atoms%nlo(atomType)
                  IF(atoms%llo(ilo,atomType).NE.l) CYCLE
                  DO ilop = 1, atoms%nlo(atomType)
                     IF(atoms%llo(ilop,atomType).NE.lp) CYCLE
                     uu_tmp(:atoms%jri(atomType)) = flo(:atoms%jri(atomType),1,ilo,j1,atomTypep)*flo(:atoms%jri(atomType),1,ilop,j2,atomType)&
                                                  + flo(:atoms%jri(atomType),2,ilo,j1,atomTypep)*flo(:atoms%jri(atomType),2,ilop,j2,atomType)
                     CALL intgr3(uu_tmp,atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),&
                                 this%uloulopn(ilo,ilop,j1,j2))
                  ENDDO
               ENDDO
            ELSE
               !Full radial dependence (We need to multiply each term with rmesh(atomtype)*rmesh(atomtypep) to get the right normalization)
               DO jri = 1, atoms%jri(atomTypep)
                  uu_tmp2(:atoms%jri(atomType)) = (f(jri,1,lp,j1,atomTypep)*f(:atoms%jri(atomType),1,l,j2,atomType)&
                                                 + f(jri,2,lp,j1,atomTypep)*f(:atoms%jri(atomType),2,l,j2,atomType)) &
                                                 * atoms%rmsh(jri,atomTypep) * atoms%rmsh(:atoms%jri(atomType),atomType)
                  CALL intgr3(uu_tmp2,atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),uu_tmp(jri))
               ENDDO
               CALL intgr3(uu_tmp,atoms%rmsh(:,atomTypep),atoms%dx(atomTypep),atoms%jri(atomTypep), &
                           this%uun(j1,j2))
               DO jri = 1, atoms%jri(atomTypep)
                  uu_tmp2(:atoms%jri(atomType)) = (f(jri,1,lp,j1,atomTypep)*g(:atoms%jri(atomType),1,l,j2,atomType)&
                                                 + f(jri,2,lp,j1,atomTypep)*g(:atoms%jri(atomType),2,l,j2,atomType)) &
                                                 * atoms%rmsh(jri,atomTypep) * atoms%rmsh(:atoms%jri(atomType),atomType)
                  CALL intgr3(uu_tmp2,atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),uu_tmp(jri))
               ENDDO
               CALL intgr3(uu_tmp,atoms%rmsh(:,atomTypep),atoms%dx(atomTypep),atoms%jri(atomTypep), &
                           this%udn(j1,j2))
               DO jri = 1, atoms%jri(atomTypep)
                  uu_tmp2(:atoms%jri(atomType)) = (g(jri,1,lp,j1,atomTypep)*f(:atoms%jri(atomType),1,l,j2,atomType)&
                                                 + g(jri,2,lp,j1,atomTypep)*f(:atoms%jri(atomType),2,l,j2,atomType)) &
                                                 * atoms%rmsh(jri,atomTypep) * atoms%rmsh(:atoms%jri(atomType),atomType)
                  CALL intgr3(uu_tmp2,atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),uu_tmp(jri))
               ENDDO
               CALL intgr3(uu_tmp,atoms%rmsh(:,atomTypep),atoms%dx(atomTypep),atoms%jri(atomTypep), &
                           this%dun(j1,j2))
               DO jri = 1, atoms%jri(atomTypep)
                  uu_tmp2(:atoms%jri(atomType)) = (g(jri,1,lp,j1,atomTypep)*g(:atoms%jri(atomType),1,l,j2,atomType)&
                                                 + g(jri,2,lp,j1,atomTypep)*g(:atoms%jri(atomType),2,l,j2,atomType)) &
                                                 * atoms%rmsh(jri,atomTypep) * atoms%rmsh(:atoms%jri(atomType),atomType)
                  CALL intgr3(uu_tmp2,atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),uu_tmp(jri))
               ENDDO
               CALL intgr3(uu_tmp,atoms%rmsh(:,atomTypep),atoms%dx(atomTypep),atoms%jri(atomTypep), &
                           this%ddn(j1,j2))
               DO ilo = 1, atoms%nlo(atomType)
                  IF(atoms%llo(ilo,atomType).NE.l) CYCLE

                  DO jri = 1, atoms%jri(atomTypep)
                     uu_tmp2(:atoms%jri(atomType)) = (f(jri,1,lp,j1,atomTypep)*flo(:atoms%jri(atomType),1,ilo,j2,atomType)&
                                                    + f(jri,2,lp,j1,atomTypep)*flo(:atoms%jri(atomType),2,ilo,j2,atomType)) &
                                                    * atoms%rmsh(jri,atomTypep) * atoms%rmsh(:atoms%jri(atomType),atomType)
                     CALL intgr3(uu_tmp2,atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),uu_tmp(jri))
                  ENDDO
                  CALL intgr3(uu_tmp,atoms%rmsh(:,atomTypep),atoms%dx(atomTypep),atoms%jri(atomTypep), &
                              this%uulon(ilo,j1,j2))
                  DO jri = 1, atoms%jri(atomTypep)
                     uu_tmp2(:atoms%jri(atomType)) = (g(jri,1,lp,j1,atomTypep)*flo(:atoms%jri(atomType),1,ilo,j2,atomType)&
                                                    + g(jri,2,lp,j1,atomTypep)*flo(:atoms%jri(atomType),2,ilo,j2,atomType)) &
                                                    * atoms%rmsh(jri,atomTypep) * atoms%rmsh(:atoms%jri(atomType),atomType)
                     CALL intgr3(uu_tmp2,atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),uu_tmp(jri))
                  ENDDO
                  CALL intgr3(uu_tmp,atoms%rmsh(:,atomTypep),atoms%dx(atomTypep),atoms%jri(atomTypep), &
                              this%dulon(ilo,j1,j2))
               ENDDO
               DO ilo = 1, atoms%nlo(atomTypep)
                  IF(atoms%llo(ilo,atomTypep).NE.lp) CYCLE

                  DO jri = 1, atoms%jri(atomTypep)
                     uu_tmp2(:atoms%jri(atomType)) = (flo(jri,1,ilo,j1,atomTypep)*f(:atoms%jri(atomType),1,l,j2,atomType)&
                                                    + flo(jri,2,ilo,j1,atomTypep)*f(:atoms%jri(atomType),2,l,j2,atomType)) &
                                                    * atoms%rmsh(jri,atomTypep) * atoms%rmsh(:atoms%jri(atomType),atomType)
                     CALL intgr3(uu_tmp2,atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),uu_tmp(jri))
                  ENDDO
                  CALL intgr3(uu_tmp,atoms%rmsh(:,atomTypep),atoms%dx(atomTypep),atoms%jri(atomTypep), &
                              this%uloun(ilo,j1,j2))
                  DO jri = 1, atoms%jri(atomTypep)
                     uu_tmp2(:atoms%jri(atomType)) = (flo(jri,1,ilo,j1,atomTypep)*g(:atoms%jri(atomType),1,l,j2,atomType)&
                                                    + flo(jri,2,ilo,j1,atomTypep)*g(:atoms%jri(atomType),2,l,j2,atomType)) &
                                                    * atoms%rmsh(jri,atomTypep) * atoms%rmsh(:atoms%jri(atomType),atomType)
                     CALL intgr3(uu_tmp2,atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),uu_tmp(jri))
                  ENDDO
                  CALL intgr3(uu_tmp,atoms%rmsh(:,atomTypep),atoms%dx(atomTypep),atoms%jri(atomTypep), &
                              this%ulodn(ilo,j1,j2))
               ENDDO

               DO ilo = 1, atoms%nlo(atomType)
                  IF(atoms%llo(ilo,atomType).NE.l) CYCLE
                  DO ilop = 1, atoms%nlo(atomTypep)
                     IF(atoms%llo(ilop,atomTypep).NE.lp) CYCLE
                     DO jri = 1, atoms%jri(atomTypep)
                        uu_tmp2(:atoms%jri(atomType)) = (flo(jri,1,ilop,j1,atomTypep)*flo(:atoms%jri(atomType),1,ilo,j2,atomType)&
                                                       + flo(jri,2,ilop,j1,atomTypep)*flo(:atoms%jri(atomType),2,ilo,j2,atomType)) &
                                                       * atoms%rmsh(jri,atomTypep) * atoms%rmsh(:atoms%jri(atomType),atomType)
                        CALL intgr3(uu_tmp2,atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),uu_tmp(jri))
                     ENDDO
                     CALL intgr3(uu_tmp,atoms%rmsh(:,atomTypep),atoms%dx(atomTypep),atoms%jri(atomTypep), &
                                 this%uloulopn(ilo,ilop,j1,j2))
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO

      this%done = .TRUE.
      CALL timestop("Offdiagonal Scalar Product")

   END SUBROUTINE addOffdScalarProduct

END MODULE m_types_scalarGF