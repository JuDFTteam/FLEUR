!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   This module generates the little group of k and the extended irr. !
!   BZ. Furthermore it calculates the irr. representation             !
!                                                                     !
!   P(R,T)\phi_n,k = \sum_{n'} rep_v(n',n) *\phi_n',k        !
!   where                                                             !
!         P  is an element of the little group of k                   !
!         n' runs over group of degenerat states belonging to n.      !
!                                                                     !
!                                             M.Betzinger (09/07)     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE m_symm_hf
  use m_judft
  USE m_types_hybdat

#define irreps .false.

   USE m_constants
   USE m_types
   USE m_util
   USE m_intgrf
   USE m_io_hybrid

CONTAINS

   SUBROUTINE symm_hf_init(sym, kpts, nk, nsymop, rrot, psym)

      IMPLICIT NONE

      TYPE(t_sym), INTENT(IN)    :: sym
      TYPE(t_kpts), INTENT(IN)    :: kpts
      INTEGER, INTENT(IN)    :: nk
      INTEGER, INTENT(OUT)   :: nsymop
      INTEGER, INTENT(INOUT) :: rrot(:,:,:) ! 3,3,sym%nsym
      INTEGER, INTENT(INOUT) :: psym(:) ! Note: psym is only filled up to index nsymop

      INTEGER :: i
      REAL    :: rotkpt(3)

      ! calculate rotations in reciprocal space
      DO i = 1, sym%nsym
         IF (i <= sym%nop) THEN
            rrot(:, :, i) = transpose(sym%mrot(:, :, sym%invtab(i)))
         ELSE
            rrot(:, :, i) = -rrot(:, :, i - sym%nop)
         END IF
      END DO

      ! determine little group of k., i.e. those symmetry operations
      ! which keep bk(:,nk) invariant
      ! nsymop :: number of such symmetry-operations
      ! psym   :: points to the symmetry-operation

      psym = 0
      nsymop = 0
      DO i = 1, sym%nsym
         rotkpt = matmul(rrot(:, :, i), kpts%bkf(:, nk))

         !transfer rotkpt into BZ
         call judft_error("Missing functionality")
         !rotkpt = modulo1(rotkpt, kpts%nkpt3)

         !check if rotkpt is identical to bk(:,nk)
         IF (maxval(abs(rotkpt - kpts%bkf(:, nk))) <= 1E-07) THEN
            nsymop = nsymop + 1
            psym(nsymop) = i
         END IF
      END DO

      WRITE (6, '(A,i3)') ' nk', nk
      WRITE (6, '(A,3f10.5)') ' kpts%bkf(:,nk):', kpts%bkf(:, nk)
      WRITE (6, '(A,i3)') ' Number of elements in the little group:', nsymop

   END SUBROUTINE symm_hf_init

   SUBROUTINE symm_hf(kpts, nk, sym, hybdat, eig_irr, input,atoms, mpdata, hybrid, cell, &
                      lapw, jsp, rrot, nsymop, psym, nkpt_EIBZ, n_q, parent, &
                      pointer_EIBZ, nsest, indx_sest)

      USE m_olap
      USE m_trafo

      IMPLICIT NONE

      TYPE(t_hybdat), INTENT(IN)   :: hybdat

      TYPE(t_input), INTENT(IN)   :: input
      TYPE(t_mpdata), intent(in) :: mpdata
      TYPE(t_hybrid), INTENT(IN) :: hybrid
      TYPE(t_sym), INTENT(IN)    :: sym
      TYPE(t_cell), INTENT(IN)   :: cell
      TYPE(t_kpts), INTENT(IN)   :: kpts
      TYPE(t_atoms), INTENT(IN)  :: atoms
      TYPE(t_lapw), INTENT(IN)   :: lapw

!     - scalars -
      INTEGER, INTENT(IN)              :: nk
      INTEGER, INTENT(IN)              :: jsp
      INTEGER, INTENT(OUT)             :: nkpt_EIBZ
      INTEGER, INTENT(IN)              :: nsymop

!     - arrays -
      INTEGER, INTENT(IN)              :: rrot(:,:,:)
      INTEGER, INTENT(IN)              :: psym(:)
      INTEGER, INTENT(OUT)             :: parent(kpts%nkptf)
      INTEGER, INTENT(OUT)             :: nsest(hybrid%nbands(nk)), indx_sest(hybrid%nbands(nk), hybrid%nbands(nk))
      INTEGER, ALLOCATABLE, INTENT(OUT) :: pointer_EIBZ(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: n_q(:)

      REAL, INTENT(IN)                 :: eig_irr(:,:)

!     - local scalars -
      INTEGER                         :: ikpt, ikpt1, iop, isym, iisym, m
      INTEGER                         :: itype, ieq, iatom, ratom
      INTEGER                         :: iband, iband1, iband2, iatom0
      INTEGER                         :: i, j, ic, ic1, ic2
      INTEGER                         :: ok
      INTEGER                         :: l, lm
      INTEGER                         :: n1, n2, nn
      INTEGER                         :: ndb, ndb1, ndb2
      INTEGER                         :: nrkpt
      INTEGER                         :: maxndb, nddb

      REAL                            :: tolerance, pi

      COMPLEX                         :: cdum
      COMPLEX, PARAMETER             :: img = (0.0, 1.0)

!     - local arrays -
      INTEGER                         :: neqvkpt(kpts%nkptf)
      INTEGER                         :: list(kpts%nkptf)
      INTEGER                         :: degenerat(hybrid%ne_eig(nk))

      REAL                            :: rotkpt(3), g(3)
      REAL, ALLOCATABLE             :: olapmt(:, :, :, :)

      COMPLEX                         :: cmt(input%neig, hybrid%maxlmindx, atoms%nat)
      COMPLEX                         :: carr1(hybrid%maxlmindx, atoms%nat)
      COMPLEX, ALLOCATABLE             :: carr(:), wavefolap(:, :)
      COMPLEX, ALLOCATABLE             :: cmthlp(:, :, :)
      COMPLEX, ALLOCATABLE             :: cpwhlp(:, :)
      COMPLEX, ALLOCATABLE             :: trace(:, :)

      TYPE(t_mat)                     :: olappw, z
      COMPLEX, ALLOCATABLE             :: rep_d(:, :, :)
      LOGICAL, ALLOCATABLE             :: symequivalent(:, :)

      WRITE (6, '(A)') new_line('n')//new_line('n')//'### subroutine: symm ###'

      ! determine extented irreducible BZ of k ( EIBZ(k) ), i.e.
      ! those k-points, which can generate the whole BZ by
      ! applying the symmetry operations of the little group of k

      neqvkpt = 0

      DO i = 1, kpts%nkptf
         list(i) = i - 1
      END DO

      DO ikpt = 2, kpts%nkptf
         DO iop = 1, nsymop

            rotkpt = matmul(rrot(:, :, psym(iop)), kpts%bkf(:, ikpt))

            !transfer rotkpt into BZ
            call judft_error("Missing functionality")
            !rotkpt = modulo1(rotkpt, kpts%nkpt3)

            !determine number of rotkpt
            nrkpt = 0
            DO ikpt1 = 1, kpts%nkptf
               IF (maxval(abs(rotkpt - kpts%bkf(:, ikpt1))) <= 1E-06) THEN
                  nrkpt = ikpt1
                  EXIT
               END IF
            END DO
            IF (nrkpt == 0) call judft_error('symm: Difference vector not found !')

            IF (list(nrkpt) /= 0) THEN
               list(nrkpt) = 0
               neqvkpt(ikpt) = neqvkpt(ikpt) + 1
               parent(nrkpt) = ikpt
            END IF
            IF (all(list == 0)) EXIT

         END DO
      END DO

      ! for the Gamma-point holds:
      parent(1) = 1
      neqvkpt(1) = 1

      ! determine number of members in the EIBZ(k)
      ic = 0
      DO ikpt = 1, kpts%nkptf
         IF (parent(ikpt) == ikpt) ic = ic + 1
      END DO
      nkpt_EIBZ = ic

      allocate(pointer_EIBZ(nkpt_EIBZ))
      ic = 0
      DO ikpt = 1, kpts%nkptf
         IF (parent(ikpt) == ikpt) THEN
            ic = ic + 1
            pointer_EIBZ(ic) = ikpt
         END IF
      END DO

      WRITE (6, '(A,i5)') ' Number of k-points in the EIBZ', nkpt_EIBZ

      ! determine the factor n_q, that means the number of symmetrie operations of the little group of bk(:,nk)
      ! which keep q (in EIBZ) invariant

      allocate(n_q(nkpt_EIBZ))

      ic = 0
      n_q = 0
      DO ikpt = 1, kpts%nkptf
         IF (parent(ikpt) == ikpt) THEN
            ic = ic + 1
            DO iop = 1, nsymop
               isym = psym(iop)
               rotkpt = matmul(rrot(:, :, isym), kpts%bkf(:, ikpt))

               !transfer rotkpt into BZ
               call judft_error("Missing functionality")
               !rotkpt = modulo1(rotkpt, kpts%nkpt3)

               !check if rotkpt is identical to bk(:,ikpt)
               IF (maxval(abs(rotkpt - kpts%bkf(:, ikpt))) <= 1E-06) THEN
                  n_q(ic) = n_q(ic) + 1
               END IF
            END DO
         END IF
      END DO
      IF (ic /= nkpt_EIBZ) call judft_error('symm: failure EIBZ')

      ! calculate degeneracy:
      ! degenerat(i) = 1 state i  is not degenerat,
      ! degenerat(i) = j state i has j-1 degenerat states at {i+1,...,i+j-1}
      ! degenerat(i) = 0 state i is degenerat

      tolerance = 1E-07 !0.00001

      degenerat = 1

      WRITE (6, '(A,f10.8)') ' Tolerance for determining degenerate states=', tolerance

      DO i = 1, hybrid%nbands(nk)
         DO j = i + 1, hybrid%nbands(nk)
            IF (abs(eig_irr(i, nk) - eig_irr(j, nk)) <= tolerance) THEN
               degenerat(i) = degenerat(i) + 1
            END IF
         END DO
      END DO

      DO i = 1, hybrid%ne_eig(nk)
         IF (degenerat(i) /= 1 .or. degenerat(i) /= 0) THEN
            degenerat(i + 1:i + degenerat(i) - 1) = 0
         END IF
      END DO

      ! maximal number of degenerate bands -> maxndb
      maxndb = maxval(degenerat)

      ! number of different degenerate bands/states
      nddb = count(degenerat >= 1)

      WRITE (6, *) ' Degenerate states:'
      DO iband = 1, hybrid%nbands(nk)/5 + 1
         WRITE (6, '(5i5)') degenerat(iband*5 - 4:min(iband*5, hybrid%nbands(nk)))
      END DO

      IF (irreps) THEN
         ! calculate representation, i.e. the action of an element of
         ! the little group of k on \phi_n,k:
         ! P(R,T)\phi_n,k = \sum_{n'\in degenerat(n)} rep_v(n',n) *\phi_n',k

         ! read in cmt and z at current k-point (nk)
         CALL read_cmt(cmt, nk)
         call read_z(z, kpts%nkptf*(jsp - 1) + nk)

         allocate(rep_d(maxndb, nddb, nsymop), stat=ok)
         IF (ok /= 0) call judft_error('symm: failure allocation rep_v')

         call olappw%alloc(z%l_real, lapw%nv(jsp), lapw%nv(jsp))
         allocate(olapmt(maxval(mpdata%num_radfun_per_l), maxval(mpdata%num_radfun_per_l), 0:atoms%lmaxd, atoms%ntype), stat=ok)
         IF (ok /= 0) call judft_error('symm: failure allocation olapmt')

         olapmt = 0
         CALL wfolap_init(olappw, olapmt, lapw%gvec(:, :, jsp), atoms, mpdata, &
                          cell, hybdat%bas1, hybdat%bas2)

         allocate(cmthlp(hybrid%maxlmindx, atoms%nat, maxndb), cpwhlp(lapw%nv(jsp), maxndb), stat=ok)
         IF (ok /= 0) call judft_error('symm: failure allocation cmthlp/cpwhlp')

         DO isym = 1, nsymop
            iop = psym(isym)

            ic = 0
            DO i = 1, hybrid%nbands(nk)
               ndb = degenerat(i)
               IF (ndb >= 1) THEN
                  ic = ic + 1
                  cmthlp = 0
                  cpwhlp = 0

                  CALL waveftrafo_symm(cmthlp(:, :, :ndb), cpwhlp(:, :ndb), cmt, z%l_real, z%data_r, z%data_c, &
                                       i, ndb, nk, iop, atoms,input, mpdata, hybrid, kpts, sym, jsp, lapw)

                  DO iband = 1, ndb
                     carr1 = cmt(iband + i - 1, :, :)
                     IF (z%l_real) THEN
                        rep_d(iband, ic, isym) = wfolap_inv(carr1, z%data_r(:lapw%nv(jsp), iband + i - 1), cmthlp(:, :, iband), &
                                                            cpwhlp(:, iband), olappw%data_r, olapmt, atoms, mpdata)
                     else
                        rep_d(iband, ic, isym) = wfolap_noinv(carr1, z%data_c(:lapw%nv(jsp), iband + i - 1), cmthlp(:, :, iband), &
                                                              cpwhlp(:, iband), olappw%data_c, olapmt, atoms, mpdata)
                     endif
                  END DO

               END IF
            END DO

         END DO

         deallocate(cmthlp, cpwhlp)

         ! calculate trace of irrecudible representation
         allocate(trace(sym%nsym, nddb), stat=ok)
         IF (ok /= 0) call judft_error('symm: failure allocation trace')

         ic = 0
         trace = 0
         DO iband = 1, hybrid%nbands(nk)
            ndb = degenerat(iband)
            IF (ndb >= 1) THEN
               ic = ic + 1
               !calculate trace
               DO iop = 1, nsymop
                  isym = psym(iop)
                  DO i = 1, ndb
                     trace(isym, ic) = trace(isym, ic) + rep_d(i, ic, iop)
                  END DO
               END DO
            END IF
         END DO

         ! determine symmetry equivalent bands/irreducible representations by comparing the trace

         allocate(symequivalent(nddb, nddb), stat=ok)
         IF (ok /= 0) call judft_error('symm: failure allocation symequivalent')

         ic1 = 0
         symequivalent = .false.
         DO iband1 = 1, hybrid%nbands(nk)
            ndb1 = degenerat(iband1)
            IF (ndb1 >= 1) THEN
               ic1 = ic1 + 1
               ic2 = 0
               DO iband2 = 1, hybrid%nbands(nk)
                  ndb2 = degenerat(iband2)
                  IF (ndb2 >= 1) THEN
                     ic2 = ic2 + 1
                     IF (ndb2 == ndb1) THEN
                        ! note that this criterium is only valid for pure spatial rotations
                        ! if one combines spatial rotations with time reversal symmetry there
                        ! is no unique criteria to identify symequivalent state
                        ! however, also in the latter case the trace of the spatial rotations
                        ! for two symmetry equivalent states must be equivalent
                        IF (all(abs(trace(:sym%nop, ic1) - trace(:sym%nop, ic2)) <= 1E-8))&
                        THEN
                           symequivalent(ic2, ic1) = .true.
                        END IF
                     END IF
                  END IF
               END DO
            END IF
         END DO

      ELSE
         ! read in cmt and z at current k-point (nk)

         CALL read_cmt(cmt, nk)
         !CALL intgrf_init(atoms%ntype,atoms%jmtd,atoms%jri,atoms%dx,atoms%rmsh,hybdat%gridf)

         IF (allocated(olapmt)) deallocate(olapmt)
         allocate(olapmt(maxval(mpdata%num_radfun_per_l), maxval(mpdata%num_radfun_per_l), 0:atoms%lmaxd, atoms%ntype), stat=ok)
         IF (ok /= 0) call judft_error('symm: failure allocation olapmt')
         olapmt = 0

         DO itype = 1, atoms%ntype
            DO l = 0, atoms%lmax(itype)
               nn = mpdata%num_radfun_per_l(l, itype)
               DO n2 = 1, nn
                  DO n1 = 1, nn
                     olapmt(n1, n2, l, itype) = intgrf( &
                                   hybdat%bas1(:, n1, l, itype)*hybdat%bas1(:, n2, l, itype)&
                                  + hybdat%bas2(:, n1, l, itype)*hybdat%bas2(:, n2, l, itype),&
                                   atoms, itype, hybdat%gridf)
                  END DO
               END DO
            END DO
         END DO

         allocate(wavefolap(hybrid%nbands(nk), hybrid%nbands(nk)), carr(maxval(mpdata%num_radfun_per_l)), stat=ok)
         IF (ok /= 0) call judft_error('symm: failure allocation wfolap/maxindx')
         wavefolap = 0

         iatom = 0
         DO itype = 1, atoms%ntype
            DO ieq = 1, atoms%neq(itype)
               iatom = iatom + 1
               lm = 0
               DO l = 0, atoms%lmax(itype)
                  DO M = -l, l
                     nn = mpdata%num_radfun_per_l(l, itype)
                     DO iband1 = 1, hybrid%nbands(nk)
                        carr(:nn) = matmul(olapmt(:nn, :nn, l, itype),&
                                            cmt(iband1, lm + 1:lm + nn, iatom))
                        DO iband2 = 1, iband1
                           wavefolap(iband2, iband1)&
                         = wavefolap(iband2, iband1)&
                         + dot_product(cmt(iband2, lm + 1:lm + nn, iatom), carr(:nn))
                        END DO
                     END DO
                     lm = lm + nn
                  END DO
               END DO
            END DO
         END DO

         DO iband1 = 1, hybrid%nbands(nk)
            DO iband2 = 1, iband1
               wavefolap(iband1, iband2) = conjg(wavefolap(iband2, iband1))
            END DO
         END DO

         allocate(symequivalent(nddb, nddb), stat=ok)
         IF (ok /= 0) call judft_error('symm: failure allocation symequivalent')
         symequivalent = .false.
         ic1 = 0
         DO iband1 = 1, hybrid%nbands(nk)
            ndb1 = degenerat(iband1)
            IF (ndb1 == 0) CYCLE
            ic1 = ic1 + 1
            ic2 = 0
            DO iband2 = 1, hybrid%nbands(nk)
               ndb2 = degenerat(iband2)
               IF (ndb2 == 0) CYCLE
               ic2 = ic2 + 1
               IF (any(abs(wavefolap(iband1:iband1 + ndb1 - 1,&
                                      iband2:iband2 + ndb2 - 1)) > 1E-9)) THEN
!                .and. ndb1 .eq. ndb2 ) THEN
                  symequivalent(ic2, ic1) = .true.
               END IF
            END DO
         END DO
      END IF

      !
      ! generate index field which contain the band combinations (n1,n2),
      ! which are non zero
      !

      ic1 = 0
      indx_sest = 0
      nsest = 0
      DO iband1 = 1, hybrid%nbands(nk)
         ndb1 = degenerat(iband1)
         IF (ndb1 >= 1) ic1 = ic1 + 1
         i = 0
         DO WHILE (degenerat(iband1 - i) == 0)
            i = i + 1
         END DO
         ndb1 = degenerat(iband1 - i)
         ic2 = 0
         DO iband2 = 1, hybrid%nbands(nk)
            ndb2 = degenerat(iband2)
            IF (ndb2 >= 1) ic2 = ic2 + 1
            i = 0
            DO WHILE (degenerat(iband2 - i) == 0)
               i = i + 1
            END DO
            ndb2 = degenerat(iband2 - i)
            ! only upper triangular part
            IF (symequivalent(ic2, ic1) .and. iband2 <= iband1) THEN
!            IF( ndb1 .ne. ndb2 ) call judft_error('symm_hf: failure symequivalent')
               nsest(iband1) = nsest(iband1) + 1
               indx_sest(nsest(iband1), iband1) = iband2
            END IF
         END DO
      END DO

      !
      ! calculate representations for core states
      ! (note that for a core state, these are proportional to the Wigner D matrices)
      !
      ! Definition of the Wigner rotation matrices
      !
      !                     -1                l
      ! P(R) Y  (r)  = Y  (R  r)  =  sum     D    (R)  Y   (r)
      !       lm        lm              m'    m m'      lm'
      !

      pi = pimach()

      IF (hybdat%lmaxcd > atoms%lmaxd) STOP &
       'symm_hf: The very impropable case that hybdat%lmaxcd > atoms%lmaxd occurs'

      iatom = 0
      iatom0 = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1
            DO iop = 1, nsymop
               isym = psym(iop)
               IF (isym <= sym%nop) THEN
                  iisym = isym
               ELSE
                  iisym = isym - sym%nop
               END IF

               ratom = hybrid%map(iatom, isym)
               rotkpt = matmul(rrot(:, :, isym), kpts%bkf(:, nk))
               g = nint(rotkpt - kpts%bkf(:, nk))

               cdum = exp(-2*pi*img*dot_product(rotkpt, sym%tau(:, iisym)))* &
                        exp(2*pi*img*dot_product(g, atoms%taual(:, ratom)))
            END DO
         END DO
         iatom0 = iatom0 + atoms%neq(itype)
      END DO

   END SUBROUTINE symm_hf

   INTEGER FUNCTION symm_hf_nkpt_EIBZ(kpts, nk, sym)

      USE m_util, ONLY: modulo1
      USE m_types
      USE m_juDFT
      IMPLICIT NONE
      TYPE(t_sym), INTENT(IN)   :: sym
      TYPE(t_kpts), INTENT(IN)   :: kpts

!     - scalar input -
      INTEGER, INTENT(IN)   :: nk
!     - array input -

!     - local scalars -
      INTEGER               ::  isym, ic, iop, ikpt, ikpt1
      INTEGER               ::  nsymop, nrkpt
!     - local arrays -
      INTEGER               ::  rrot(3, 3, sym%nsym)
      INTEGER               ::  neqvkpt(kpts%nkptf), list(kpts%nkptf), parent(kpts%nkptf),&
                                symop(kpts%nkptf)
      INTEGER, ALLOCATABLE  ::  psym(:)!,help(:)
      REAL                  ::  rotkpt(3)

      ! calculate rotations in reciprocal space
      DO isym = 1, sym%nsym
         IF (isym <= sym%nop) THEN
            rrot(:, :, isym) = transpose(sym%mrot(:, :, sym%invtab(isym)))
         ELSE
            rrot(:, :, isym) = -rrot(:, :, isym - sym%nop)
         END IF
      END DO

      ! determine little group of k., i.e. those symmetry operations
      ! which keep bk(:,nk,nw) invariant
      ! nsymop :: number of such symmetry-operations
      ! psym   :: points to the symmetry-operation

      ic = 0
      allocate(psym(sym%nsym))

      DO iop = 1, sym%nsym
         rotkpt = matmul(rrot(:, :, iop), kpts%bkf(:, nk))

         !transfer rotkpt into BZ
         call judft_error("Missing functionality")
         !rotkpt = modulo1(rotkpt, kpts%nkpt3)

         !check if rotkpt is identical to bk(:,nk)
         IF (maxval(abs(rotkpt - kpts%bkf(:, nk))) <= 1E-07) THEN
            ic = ic + 1
            psym(ic) = iop
         END IF
      END DO
      nsymop = ic

      ! reallocate psym
!       ALLOCATE(help(ic))
!       help = psym(1:ic)
!       DEALLOCATE(psym)
!       ALLOCATE(psym(ic))
!       psym = help
!       DEALLOCATE(help)

      ! determine extented irreducible BZ of k ( EIBZ(k) ), i.e.
      ! those k-points, which can generate the whole BZ by
      ! applying the symmetry operations of the little group of k

      neqvkpt = 0

!       list = [(ikpt-1, ikpt=1,nkpt) ]
      DO ikpt = 1, kpts%nkptf
         list(ikpt) = ikpt - 1
      END DO

      DO ikpt = 2, kpts%nkptf
         DO iop = 1, nsymop

            rotkpt = matmul(rrot(:, :, psym(iop)), kpts%bkf(:, ikpt))

            !transfer rotkpt into BZ
            call judft_error("Missing functionality")
            !rotkpt = modulo1(rotkpt, kpts%nkpt3)

            !determine number of rotkpt
            nrkpt = 0
            DO ikpt1 = 1, kpts%nkptf
               IF (maxval(abs(rotkpt - kpts%bkf(:, ikpt1))) <= 1E-06) THEN
                  nrkpt = ikpt1
                  EXIT
               END IF
            END DO
            IF (nrkpt == 0) call judft_error('symm: Difference vector not found !')

            IF (list(nrkpt) /= 0) THEN
               list(nrkpt) = 0
               neqvkpt(ikpt) = neqvkpt(ikpt) + 1
               parent(nrkpt) = ikpt
               symop(nrkpt) = psym(iop)
            END IF
            IF (all(list == 0)) EXIT

         END DO
      END DO

      ! for the Gamma-point holds:
      parent(1) = 1
      neqvkpt(1) = 1

      ! determine number of members in the EIBZ(k)
      ic = 0
      DO ikpt = 1, kpts%nkptf
         IF (parent(ikpt) == ikpt) ic = ic + 1
      END DO
      symm_hf_nkpt_EIBZ = ic

   END FUNCTION symm_hf_nkpt_EIBZ

END MODULE m_symm_hf
