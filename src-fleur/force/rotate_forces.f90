MODULE m_rotate_forces
   ! This routine writes a file similar to forces.dat, but containing
   ! all atoms in a sequence corresponding to their appearance in the
   ! input file for the input-file generator. forces.dat only
   ! contains forces for a representant of the symmetry equivalent
   ! atoms. This routine also gives files FORCES and POSCAR to be
   ! used with the phon package by Dario Alfe (for a single
   ! displacement),
   ! Reference
   ! Klueppelberg, Jun 2012

   ! Readded to only construct POSCAR and FORCES files in Dec 2020, Neukirchen

   ! Modified to construct a file for use with phonopy instead of PHON.
   ! Neukirchen, Dec 2020 

CONTAINS
   SUBROUTINE rotate_forces(ntypd,ntype,natd,nop,tote,omtil,neq,mrot,amat,bmat,taual,tau,force,label)
      USE m_constants
      USE m_juDFT_string

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: ntypd,ntype,natd,nop
      REAL, INTENT (IN) :: tote,omtil

      INTEGER, INTENT (IN) :: neq(ntypd),mrot(3,3,nop)
      REAL, INTENT (IN) :: amat(3,3),bmat(3,3),taual(3,natd),tau(3,nop)
      REAL, INTENT (IN) :: force(3,ntypd)
      CHARACTER(LEN=20) :: label(natd)

      INTEGER :: iatom,iatom0,itype,ieq,ratom,isym,sym,i,true_atom,atom_map(natd)

      REAL :: pos(3,natd),rtaual(3),forcerot(3,3),forceval(3,natd)
      CHARACTER(len=20) :: string
      LOGICAL :: l_PHON

      l_PHON = .FALSE.
      OPEN (79,file='FORCES')
      IF (l_PHON) THEN
         OPEN (80,file='POSCAR')
      END IF
      OPEN(81,file='FORCES_SORT')

      WRITE (79,'(i1)') 1
      WRITE (79,'(i1,1x,a)') 1,'#'
      WRITE (81,'(i1)') 1
      WRITE (81,'(i1,1x,a)') 1,'#'
      IF (l_PHON) THEN
         WRITE (80,'(a)') "#insert lattice name"
         WRITE (80,'(f20.10)') -omtil*bohr_to_angstrom_const**3
         WRITE (80,FMT=800) amat(1,1:3)*bohr_to_angstrom_const
         WRITE (80,FMT=800) amat(2,1:3)*bohr_to_angstrom_const
         WRITE (80,FMT=800) amat(3,1:3)*bohr_to_angstrom_const
         WRITE (string,'(i3)') ntype
         WRITE (80,'('//string//'(i2,1x))') (neq(itype),itype=1,ntype)
         WRITE (80,'(a)') 'Direct'
      END IF

      iatom  = 0
      iatom0 = 1
      DO itype = 1,ntype
         DO ieq = 1,neq(itype)
            iatom = iatom+1
            true_atom = str2int(TRIM(label(iatom)))
            atom_map(true_atom) = iatom
            pos(:,iatom) = matmul(amat,taual(:,iatom))

            ratom = 0
            sym = 0

            DO isym = 1,nop
               rtaual(:)=matmul(mrot(:,:,isym),taual(:,iatom0))+tau(:,isym)
               IF(all(abs(modulo(rtaual-taual(:,iatom)+0.5,1d0)-0.5).lt.1d-10)) THEN
                  sym = isym
                  EXIT
               END IF
            END DO

            IF (sym.eq.0) THEN
               sym = 1
            END IF

            forcerot = matmul(amat,matmul(mrot(:,:,sym),bmat/tpi_const))
            forceval(:,iatom) = matmul(forcerot,force(:,itype))

            IF (l_PHON) THEN
               WRITE (79,FMT=790) forceval(1:3,iatom)*hartree_to_ev_const/bohr_to_angstrom_const
               WRITE (80,FMT=800) taual(1:3,iatom)
            ELSE
               WRITE (79,*) forceval(1:3,iatom), 'force'
            END IF
         END DO
         iatom0 = iatom0 + neq(itype)
      END DO

      DO iatom = 1, natd
         WRITE (81,*) forceval(1:3,atom_map(iatom)), 'force'
      END DO

790   FORMAT (1x,3(1x,f20.10))
800   FORMAT (1x,3(1x,f20.10))
810   FORMAT (1x,a4,43x,3f14.9)

      CLOSE(79)
      IF (l_PHON) THEN
         CLOSE(80)
      END IF
      CLOSE(81)

   END SUBROUTINE rotate_forces

END MODULE m_rotate_forces
