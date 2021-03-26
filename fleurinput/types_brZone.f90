!--------------------------------------------------------------------------------
! Copyright (c) 2020 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_brZone

   IMPLICIT NONE

   INTEGER, PARAMETER :: nop48_const  = 48
   INTEGER, PARAMETER :: mface_const  = 51
   INTEGER, PARAMETER :: nbsz_const   =  3
   INTEGER, PARAMETER :: nv48_const   = (2*3+1)**3+48 !(2*nbsz+1)**3+48

   TYPE :: t_brZone

      INTEGER :: idsyst             ! crystal system identification in MDDFT programs
      INTEGER :: idtype             ! lattice type identification in MDDFT programs
      INTEGER :: nsym               ! number of symmetry operations
      REAL    :: bltv(3,3)          ! cartesian Bravais lattice basis (a.u.)
      REAL    :: rltv(3,3)          ! reciprocal lattice basis (2\pi/a.u.)
      REAL    :: ccr(3,3,nop48_const)     ! rotation matrices in cartesian repr.
      REAL    :: rlsymr(3,3,nop48_const)  ! rotation matrices in reciprocal lattice basis representation
      REAL    :: talfa(3,nop48_const)     ! translation vector associated with (non-symmorphic)
      ! symmetry elements in Bravais lattice representation
      INTEGER :: ncorn,nedge,nface  ! number of corners, faces and edges of the IBZ
      REAL    :: fnorm(3,mface_const)     ! normal vector of the planes bordering the IBZ
      REAL    :: fdist(mface_const)       ! distance vector of the planes bordering t IBZ
      REAL    :: cpoint(3,mface_const)    ! cartesian coordinates of corner points of IBZ
      REAL    :: xvec(3)            ! arbitrary vector lying in the IBZ

      CONTAINS

      PROCEDURE :: initBZone

   END TYPE t_brZone

   PUBLIC t_brZone

   CONTAINS

   SUBROUTINE initBZone(bz, cell, sym, l_soc_or_ss, film, l_onlyIdentitySym)

      USE m_constants
      USE m_types_cell
      USE m_types_sym
      USE m_bravais
      USE m_brzone2

      IMPLICIT NONE

      CLASS(t_brZone), INTENT(INOUT) :: bz
      TYPE(t_cell),    INTENT(IN)    :: cell
      TYPE(t_sym),     INTENT(IN)    :: sym
      LOGICAL,         INTENT(IN)    :: l_soc_or_ss
      LOGICAL,         INTENT(IN)    :: film
      LOGICAL,         INTENT(IN)    :: l_onlyIdentitySym

      INTEGER :: addSym, i
      REAL    :: binv(3,3)

      CALL bravais(cell%amat,bz%idsyst,bz%idtype)

      bz%nsym = MERGE(sym%nop2,sym%nop,film)
      IF(l_onlyIdentitySym) bz%nsym = 1

      ! Lattice information
      bz%bltv=TRANSPOSE(cell%amat)
      binv=TRANSPOSE(cell%bmat)/tpi_const

      bz%talfa(:,:bz%nsym)=MATMUL(bz%bltv,sym%tau(:,:bz%nsym))

      DO i = 1, bz%nsym
         bz%rlsymr(:,:,i)=REAL(sym%mrot(:,:,i))
         bz%ccr(:,:,i) = TRANSPOSE(MATMUL(MATMUL(binv(:,:),TRANSPOSE(bz%rlsymr(:,:,i))),bz%bltv(:,:)))
      END DO

      IF (.NOT.l_OnlyIdentitySym.AND..NOT.l_soc_or_ss.AND.(2*bz%nsym<nop48_const)) THEN
         IF ((film.AND.(.NOT.sym%invs2)).OR.((.NOT.film).AND.(.NOT.sym%invs))) THEN
            addSym = 0
            ! Note: We have to add the negative of each symmetry operation
            !       to exploit time reversal symmetry. However, if the new
            !       symmetry operation is the identity matrix it is excluded.
            !       This is the case iff it is (-Id) + a translation vector.
            DO i = 1, bz%nsym
               ! This test assumes that ccr(:,:,1) is the identity matrix.
               IF(.NOT.ALL(ABS(bz%ccr(:,:,1)+bz%ccr(:,:,i)).LT.10e-10) ) THEN
                  bz%ccr(:,:,bz%nsym+addSym+1 ) = -bz%ccr(:,:,i)
                  bz%rlsymr(:,:,bz%nsym+addSym+1 ) = -bz%rlsymr(:,:,i)
                  addSym = addSym + 1
               END IF
            END DO
            bz%nsym = bz%nsym + addSym
         END IF
      END IF

      ! brzone and brzone2 find the corner-points, the edges, and the
      ! faces of the irreducible wedge of the brillouin zone (IBZ).
      bz%rltv=TRANSPOSE(cell%bmat)
      CALL brzone2(bz%rltv,bz%nsym,bz%ccr,mface_const,nbsz_const,nv48_const,bz%cpoint,bz%xvec,bz%ncorn,bz%nedge,bz%nface,bz%fnorm,bz%fdist)
   END SUBROUTINE initBZone

END MODULE m_types_brZone
