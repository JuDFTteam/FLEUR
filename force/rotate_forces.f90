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

CONTAINS
   SUBROUTINE rotate_forces(ntypd,ntype,natd,nop,tote,omtil,neq,mrot,amat,bmat,taual,tau,force)
      USE m_constants

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: ntypd,ntype,natd,nop
      REAL, INTENT (IN) :: tote,omtil

      INTEGER, INTENT (IN) :: neq(ntypd),mrot(3,3,nop)
      REAL, INTENT (IN) :: amat(3,3),bmat(3,3),taual(3,natd),tau(3,nop)
      REAL, INTENT (IN) :: force(3,ntypd)

      INTEGER :: iatom,iatom0,itype,ieq,ratom,isym,sym,i
      REAL :: hartreetoev,autoangstrom

      REAL :: pos(3,natd),rtaual(3),forcerot(3,3),forceval(3,natd)
      CHARACTER(len=20) :: string

      hartreetoev = 27.211386245988
      autoangstrom = 0.529177210903

      OPEN (79,file='FORCES')
      OPEN (80,file='POSCAR')

      WRITE (79,'(i1)') 1
      WRITE (79,'(i1,1x,a)') 1,'#insert displacement here (also, edit number of displacement/total displacements'
      WRITE (80,'(a)') "#insert lattice name"
      WRITE (80,'(f20.10)') -omtil*autoangstrom**3
      WRITE (80,FMT=800) amat(1,1:3)*autoangstrom
      WRITE (80,FMT=800) amat(2,1:3)*autoangstrom
      WRITE (80,FMT=800) amat(3,1:3)*autoangstrom
      WRITE (string,'(i3)') ntype
      WRITE (80,'('//string//'(i2,1x))') (neq(itype),itype=1,ntype)
      WRITE (80,'(a)') 'Direct'

      iatom  = 0
      iatom0 = 1
      DO itype = 1,ntype
         DO ieq = 1,neq(itype)
            iatom = iatom+1

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

            WRITE (79,FMT=790) forceval(1:3,iatom)*hartreetoev/autoangstrom
            WRITE (80,FMT=800) taual(1:3,iatom)

         END DO
         iatom0 = iatom0 + neq(itype)
      END DO

790   FORMAT (1x,3(1x,f20.10))
800   FORMAT (1x,3(1x,f20.10))
810   FORMAT (1x,a4,43x,3f14.9)

      CLOSE(79);CLOSE(80)

   END SUBROUTINE rotate_forces

END MODULE m_rotate_forces
