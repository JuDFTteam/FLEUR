!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_flipcdn
!     *******************************************************
!     this subroutine reads the charge density and flips the 
!     magnetic moment within the m.t.sphere for each atom 
!     according to the variable nflip. This variable is read in
!     the main program
!             nflip = -1 : flip spin in sphere
!             nflip = -2 : scale spin by bmu(n)
!             nflip = any: no spin flip
!                            r.pentcheva,kfa,Feb'96
!
!     Extension to multiple U per atom type by G.M. 2017
!     *******************************************************
CONTAINS

SUBROUTINE flipcdn(atoms,input,vacuum,sphhar,stars,sym,noco,oneD,cell)

   USE m_constants
   USE m_cdn_io
   USE m_types

   IMPLICIT NONE

   TYPE(t_stars),INTENT(IN)    :: stars
   TYPE(t_vacuum),INTENT(IN)   :: vacuum
   TYPE(t_atoms),INTENT(IN)    :: atoms
   TYPE(t_sphhar),INTENT(IN)   :: sphhar
   TYPE(t_input),INTENT(INOUT) :: input
   TYPE(t_sym),INTENT(IN)      :: sym
   TYPE(t_noco),INTENT(IN)     :: noco
   TYPE(t_oneD),INTENT(IN)     :: oneD
   TYPE(t_cell),INTENT(IN)     :: cell

   ! Local type instance
   TYPE(t_potden)            :: den

   ! Local Scalars
   REAL                      :: rhodummy,rhodumms,fermiEnergyTemp
   INTEGER                   :: i,nt,j,lh,na,mp,ispin,urec,itype,m,i_u
   INTEGER                   :: archiveType
   LOGICAL                   :: n_exist,l_qfix,l_error

   ! Local Arrays
   CHARACTER(len=80), ALLOCATABLE :: clines(:)

   CALL den%init(stars,atoms,sphhar,vacuum,noco,oneD,input%jspins,.FALSE.,POTDEN_TYPE_DEN)
   IF(noco%l_noco) THEN
      archiveType = CDN_ARCHIVE_TYPE_NOCO_const
   ELSE
      archiveType = CDN_ARCHIVE_TYPE_CDN1_const
   END IF

   ! read the charge density 
   CALL readDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,&
                    CDN_INPUT_DEN_const,0,fermiEnergyTemp,l_qfix,den)

   ! flip cdn for each atom with nflip=-1
   na = 1
   DO itype = 1, atoms%ntype
      IF (atoms%nflip(itype).EQ.-1) THEN
         ! spherical and non-spherical m.t. charge density
         DO lh = 0,sphhar%nlh(atoms%ntypsy(na))
            DO j = 1,atoms%jri(itype)
               rhodummy = den%mt(j,lh,itype,1)
               den%mt(j,lh,itype,1) = den%mt(j,lh,itype,input%jspins)
               den%mt(j,lh,itype,input%jspins) = rhodummy
            END DO
         END DO
      ELSE IF (atoms%nflip(itype).EQ.-2) THEN
         DO lh = 0,sphhar%nlh(atoms%ntypsy(na))
            DO j = 1,atoms%jri(itype)
               rhodummy = den%mt(j,lh,itype,1) + den%mt(j,lh,itype,input%jspins)
               rhodumms = den%mt(j,lh,itype,1) - den%mt(j,lh,itype,input%jspins)
               den%mt(j,lh,itype,1) = 0.5 * (rhodummy + atoms%bmu(itype)*rhodumms)
               den%mt(j,lh,itype,input%jspins) = 0.5 * (rhodummy - atoms%bmu(itype)*rhodumms )
            END DO
         END DO
      END IF
         na = na + atoms%neq(itype)
   END DO

   ! for LDA+U: flip density matrix
   IF (ANY(den%mmpMat(:,:,:,:).NE.0.0).AND.atoms%n_u>0) THEN
      DO i_u = 1, atoms%n_u
         itype = atoms%lda_u(i_u)%atomType
         IF (atoms%nflip(itype).EQ.-1) THEN
            DO m = -3,3
               DO mp = -3,3
                  rhodummy = den%mmpMat(m,mp,i_u,1)
                  den%mmpMat(m,mp,i_u,1) = den%mmpMat(m,mp,i_u,input%jspins)
                  den%mmpMat(m,mp,i_u,input%jspins) = rhodummy
               END DO
            END DO
         ELSE IF (atoms%nflip(itype).EQ.-2) THEN
            DO m = -3,3
               DO mp = -3,3
                  rhodummy = den%mmpMat(m,mp,i_u,1) + den%mmpMat(m,mp,i_u,input%jspins)
                  rhodumms = den%mmpMat(m,mp,i_u,1) - den%mmpMat(m,mp,i_u,input%jspins)
                  den%mmpMat(m,mp,i_u,1) = 0.5 * (rhodummy + atoms%bmu(itype) * rhodumms)
                  den%mmpMat(m,mp,i_u,input%jspins) = 0.5 * (rhodummy - atoms%bmu(itype) * rhodumms)
               END DO
            END DO
         END IF
      END DO
   END IF

   ! write the spin-polarized density
   CALL writeDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,CDN_INPUT_DEN_const,&
                     0,-1.0,0.0,.FALSE.,den)

   ! read enpara and  flip lines
   INQUIRE(file='enpara',exist=n_exist)
   IF (n_exist) THEN
      OPEN(40,file ='enpara',status='old',form='formatted')

      j = 2
      DO itype = 1, atoms%ntype
         j = j + 1
         IF (atoms%nlo(itype)>0) j = j + 2
      END DO
      IF (input%film) j = j + 1
      ALLOCATE (clines(2*j))
      DO i = 1, 2*j
         READ (40,'(a)') clines(i)
      END DO

      REWIND 40
      i = 0 
      DO ispin = 1,input%jspins
         i = i + 2
         WRITE (40,'(a)') TRIM(clines(i-1))
         WRITE (40,'(a)') TRIM(clines(i))
         DO itype = 1, atoms%ntype
            i = i + 1
            m = i
            IF (atoms%nflip(itype)==-1) m = MOD(i+j,2*j)
            IF (m==0) m = 2*j
            WRITE (40,'(a)') TRIM(clines(m))
            IF (atoms%nlo(itype)>0) THEN
               WRITE (40,'(a)') TRIM(clines(m+1))
               WRITE (40,'(a)') TRIM(clines(m+2))
               i = i + 2
            END IF
         END DO
         IF (input%film) THEN
            i = i + 1
            WRITE (40,'(a)') TRIM(clines(i))
         END IF
      END DO
      DEALLOCATE (clines)
      CLOSE(40)
   END IF

END SUBROUTINE flipcdn

END MODULE m_flipcdn
