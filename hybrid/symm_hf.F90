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
   USE m_types
   USE m_types_hybdat
   USE m_constants
   USE m_util
   USE m_intgrf
   USE m_io_hybrid
#ifdef CPP_MPI 
   use mpi 
#endif
CONTAINS

   SUBROUTINE symm_hf_init(fi, nk, nsymop, rrot, psym)
      use m_juDFT
      IMPLICIT NONE

      type(t_fleurinput), intent(in)    :: fi
      INTEGER, INTENT(IN)    :: nk
      INTEGER, INTENT(INOUT) :: nsymop
      INTEGER, INTENT(INOUT) :: rrot(:, :, :) ! 3,3,fi%sym%nsym
      INTEGER, INTENT(INOUT) :: psym(:) ! Note: psym is only filled up to index nsymop

      INTEGER :: i
      REAL    :: rotkpt(3)

      CALL timestart("symm_hf_init")
      nsymop = 0
      ! calculate rotations in reciprocal space
      DO i = 1, fi%sym%nsym
         IF(i <= fi%sym%nop) THEN
            rrot(:, :, i) = transpose(fi%sym%mrot(:, :, fi%sym%invtab(i)))
         ELSE
            rrot(:, :, i) = -rrot(:, :, i - fi%sym%nop)
         END IF
      END DO

      ! determine little group of k., i.e. those symmetry operations
      ! which keep bk(:,nk) invariant
      ! nsymop :: number of such symmetry-operations
      ! psym   :: points to the symmetry-operation

      psym = 0
      nsymop = 0
      DO i = 1, fi%sym%nsym
         rotkpt = matmul(rrot(:, :, i), fi%kpts%bkf(:, nk))

         !transfer rotkpt into BZ
         rotkpt = fi%kpts%to_first_bz(rotkpt)

         !check if rotkpt is identical to bk(:,nk)
         IF(maxval(abs(rotkpt - fi%kpts%bkf(:, nk))) <= 1E-07) THEN
            nsymop = nsymop + 1
            psym(nsymop) = i
         END IF
      END DO
      CALL timestop("symm_hf_init")
   END SUBROUTINE symm_hf_init

   SUBROUTINE symm_hf(fi, nk, hybdat, submpi, eig_irr, mpdata, cmt, &
                      rrot, nsymop, psym, n_q, parent, nsest, indx_sest, jsp)

      USE m_olap
      USE m_trafo
      use m_calc_cmt
      use m_juDFT

      IMPLICIT NONE

      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_hybdat), INTENT(IN) :: hybdat
      type(t_hybmpi), intent(in) :: submpi
      TYPE(t_mpdata), intent(in) :: mpdata

!     - scalars -
      INTEGER, INTENT(IN)              :: nk, jsp
      INTEGER, INTENT(IN)              :: nsymop

!     - arrays -
      COMPLEX , intent(in)             :: cmt(hybdat%nbands(nk,jsp), hybdat%maxlmindx, fi%atoms%nat)
      INTEGER, INTENT(IN)              :: rrot(:, :, :)
      INTEGER, INTENT(IN)              :: psym(:)
      INTEGER, INTENT(INOUT)           :: parent(fi%kpts%nkptf)
      INTEGER, INTENT(INOUT)           :: nsest(hybdat%nbands(nk,jsp))
      INTEGER, INTENT(INOUT)           :: indx_sest(hybdat%nbands(nk,jsp), hybdat%nbands(nk,jsp))
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: n_q(:)

      REAL, INTENT(IN)                 :: eig_irr(:, :)

!     - local scalars -
      INTEGER                         :: ikpt, iop, isym, iisym, m
      INTEGER                         :: itype, ieq, iatom, ratom, ierr
      INTEGER                         ::  iband1, iband2, iatom0
      INTEGER                         :: i, j, ic, ic1, ic2
      INTEGER                         :: ok, ld_olapmt, ld_cmthlp, ld_tmp, ld_wavefolap
      INTEGER                         :: l, lm
      INTEGER                         :: n1, n2, nn
      INTEGER                         :: ndb1, ndb2
      INTEGER                         :: nrkpt
      INTEGER                         :: maxndb, nddb

      REAL                            :: tolerance, pi

      COMPLEX                         :: cdum
      COMPLEX, PARAMETER             :: img = (0.0, 1.0)

!     - local arrays -
      INTEGER                         :: neqvkpt(fi%kpts%nkptf)
      INTEGER                         :: list(fi%kpts%nkptf)
      INTEGER                         :: degenerat(hybdat%results%neig(nk, jsp))

      REAL                            :: rotkpt(3), g(3)
      complex, ALLOCATABLE            :: olapmt(:, :, :, :)

      COMPLEX, ALLOCATABLE             :: carr(:), wavefolap(:, :), tmp(:,:), carr_tmp(:)
      COMPLEX, ALLOCATABLE             :: cmthlp(:, :)
      LOGICAL, ALLOCATABLE             :: symequivalent(:, :)

      CALL timestart("symm_hf")
      parent = 0; nsest = 0; indx_sest = 0;

      ! determine extented irreducible BZ of k ( EIBZ(k) ), i.e.
      ! those k-points, which can generate the whole BZ by
      ! applying the symmetry operations of the little group of k
      call timestart("calc EIBZ")
      neqvkpt = 0

      DO i = 1, fi%kpts%nkptf
         list(i) = i - 1
      END DO

      DO ikpt = 2, fi%kpts%nkptf
         DO iop = 1, nsymop

            rotkpt = matmul(rrot(:, :, psym(iop)), fi%kpts%bkf(:, ikpt))

            !determine number of rotkpt
            nrkpt = fi%kpts%get_nk(rotkpt)
            IF(nrkpt == 0) call judft_error('symm: Difference vector not found !')

            IF(list(nrkpt) /= 0) THEN
               list(nrkpt) = 0
               neqvkpt(ikpt) = neqvkpt(ikpt) + 1
               parent(nrkpt) = ikpt
            END IF
            IF(all(list == 0)) EXIT

         END DO
      END DO

      ! for the Gamma-point holds:
      parent(1) = 1
      neqvkpt(1) = 1
      call timestop("calc EIBZ")

      ! determine the factor n_q, that means the number of symmetrie operations of the little group of bk(:,nk)
      ! which keep q (in EIBZ) invariant
      call timestart("calc n_q")
      IF(ALLOCATED(n_q)) DEALLOCATE(n_q)
      allocate(n_q(fi%kpts%EIBZ(nk)%nkpt), source=0)

      ic = 0
      n_q = 0
      DO ikpt = 1, fi%kpts%nkptf
         IF(parent(ikpt) == ikpt) THEN
            ic = ic + 1
            DO iop = 1, nsymop
               isym = psym(iop)
               rotkpt = matmul(rrot(:, :, isym), fi%kpts%bkf(:, ikpt))

               !transfer rotkpt into BZ
               rotkpt = fi%kpts%to_first_bz(rotkpt)

               !check if rotkpt is identical to bk(:,ikpt)
               IF(maxval(abs(rotkpt - fi%kpts%bkf(:, ikpt))) <= 1E-06) THEN
                  n_q(ic) = n_q(ic) + 1
               END IF
            END DO
         END IF
      END DO
      IF(ic /= fi%kpts%EIBZ(nk)%nkpt) call judft_error('symm: failure EIBZ')
      call timestop("calc n_q")

      ! calculate degeneracy:
      ! degenerat(i) = 1 state i  is not degenerat,
      ! degenerat(i) = j state i has j-1 degenerat states at {i+1,...,i+j-1}
      ! degenerat(i) = 0 state i is degenerat
      call timestart("calculate degeneracy")
      tolerance = 1E-07 !0.00001

      degenerat = 1

      DO i = 1, hybdat%nbands(nk,jsp)
         DO j = i + 1, hybdat%nbands(nk,jsp)
            IF(abs(eig_irr(i, nk) - eig_irr(j, nk)) <= tolerance) THEN
               degenerat(i) = degenerat(i) + 1
            END IF
         END DO
      END DO

      DO i = 1, hybdat%results%neig(nk, jsp)
         IF(degenerat(i) /= 1 .or. degenerat(i) /= 0) THEN
            degenerat(i + 1:i + degenerat(i) - 1) = 0
         END IF
      END DO

      ! maximal number of degenerate bands -> maxndb
      maxndb = maxval(degenerat)

      ! number of different degenerate bands/states
      nddb = count(degenerat >= 1)
      call timestop("calculate degeneracy")

      call timestart("calc olapmt")
      IF(allocated(olapmt)) deallocate(olapmt)
      allocate(olapmt(maxval(mpdata%num_radfun_per_l), maxval(mpdata%num_radfun_per_l), 0:fi%atoms%lmaxd, fi%atoms%ntype), stat=ok)
      IF(ok /= 0) call judft_error('symm: failure allocation olapmt')
      olapmt = 0

      DO itype = 1, fi%atoms%ntype
         DO l = 0, fi%atoms%lmax(itype)
            nn = mpdata%num_radfun_per_l(l, itype)
            DO n2 = 1, nn
               DO n1 = 1, nn
                  olapmt(n1, n2, l, itype) = intgrf( &
                                             hybdat%bas1(:, n1, l, itype)*hybdat%bas1(:, n2, l, itype) &
                                             + hybdat%bas2(:, n1, l, itype)*hybdat%bas2(:, n2, l, itype), &
                                             fi%atoms, itype, hybdat%gridf)
               END DO
            END DO
         END DO
      END DO
      call timestop("calc olapmt")

      allocate(wavefolap(hybdat%nbands(nk,jsp), hybdat%nbands(nk,jsp)), carr(maxval(mpdata%num_radfun_per_l)), stat=ok)
      IF(ok /= 0) call judft_error('symm: failure allocation wfolap/maxindx')
      wavefolap = 0

      allocate(cmthlp(size(cmt,2), size(cmt,1) ), stat=ierr)
      if(ierr /= 0) call judft_error("can't alloc cmthlp")

      allocate(tmp(maxval(mpdata%num_radfun_per_l), hybdat%nbands(nk,jsp)), stat=ierr)
      if(ierr /= 0 ) call judft_error("cant't alloc tmp")

      call timestart("calc wavefolap")
      ld_olapmt = size(olapmt,1)
      ld_cmthlp = size(cmthlp,1)
      ld_tmp    = size(tmp, 1)
      ld_wavefolap = size(wavefolap,1)

      do iatom = 1+submpi%rank, fi%atoms%nat, submpi%size
         itype = fi%atoms%itype(iatom)
         call timestart("transp cmthlp")
         cmthlp = transpose(cmt(:,:,iatom))
         call timestop("transp cmthlp")
         lm = 0
         DO l = 0, fi%atoms%lmax(itype)
            DO M = -l, l
               nn = mpdata%num_radfun_per_l(l, itype)

               !call zgemm(transa, transb, m,  n,                     k,  alpha,   a,                   lda,      
               call zgemm("N", "N",      nn, hybdat%nbands(nk,jsp), nn, cmplx_1, olapmt(1,1,l,itype), ld_olapmt, &
               !         b,              ldb,       beta,    c,      ldc)
                         cmthlp(lm+1,1), ld_cmthlp, cmplx_0, tmp, ld_tmp)

               !call zgemm(transa, transb, m,              n,                      k,  alpha,   a,                   lda,   
               call zgemm("C", "N", hybdat%nbands(nk,jsp), hybdat%nbands(nk,jsp), nn, cmplx_1, cmthlp(lm+1, 1), ld_cmthlp, &
               !         b,   ldb,    beta,    c,      ldc)
                         tmp, ld_tmp, cmplx_1, wavefolap, ld_wavefolap)
               lm = lm + nn
            END DO
         END DO
      END DO
      
      deallocate(cmthlp)
#ifdef CPP_MPI
      call timestart("allreduce wavefolap")
      call MPI_ALLREDUCE(MPI_IN_PLACE, wavefolap, size(wavefolap), MPI_DOUBLE_COMPLEX, MPI_SUM, submpi%comm, ierr)
      call timestop("allreduce wavefolap")
#endif
      call timestop("calc wavefolap")

      call timestart("calc symmequivalent")

      allocate(symequivalent(nddb, nddb), stat=ok, source=.False.)
      IF(ok /= 0) call judft_error('symm: failure allocation symequivalent')

      !$OMP PARALLEL DO default(none) private(iband1, ndb1, ic1, iband2, ndb2, ic2) &
      !$OMP shared(submpi, hybdat, degenerat, wavefolap, symequivalent, nk, jsp)
      DO iband1 = submpi%rank + 1, hybdat%nbands(nk,jsp), submpi%size
         ndb1 = degenerat(iband1)
         IF(ndb1 /= 0) then
            ic1 = count(degenerat(:iband1) /= 0)
            DO iband2 = 1, hybdat%nbands(nk,jsp)
               ndb2 = degenerat(iband2)
               IF(ndb2 /= 0) then
                  ic2 = count(degenerat(:iband2) /= 0)
                  IF(any(abs(wavefolap(iband1:iband1 + ndb1 - 1, &
                                       iband2:iband2 + ndb2 - 1)) > 1E-9)) THEN
                     symequivalent(ic2, ic1) = .true.
                  END IF
               endif
            END DO
         endif
      END DO
      !$OMP end parallel do
#ifdef CPP_MPI
      call timestart("allreduce symequivalent")
      call MPI_ALLREDUCE(MPI_IN_PLACE, symequivalent, size(symequivalent), MPI_LOGICAL, MPI_LOR, submpi%comm, ierr)
      call timestop("allreduce symequivalent")
#endif
      call timestop("calc symmequivalent")
      !
      ! generate index field which contain the band combinations (n1,n2),
      ! which are non zero
      !
      call timestart("calc bandcombos")
      ic1 = 0
      indx_sest = 0
      nsest = 0
      DO iband1 = 1, hybdat%nbands(nk,jsp)
         ndb1 = degenerat(iband1)
         IF(ndb1 >= 1) ic1 = ic1 + 1
         i = 0
         DO WHILE(degenerat(iband1 - i) == 0)
            i = i + 1
         END DO
         ndb1 = degenerat(iband1 - i)
         ic2 = 0
         DO iband2 = 1, hybdat%nbands(nk,jsp)
            ndb2 = degenerat(iband2)
            IF(ndb2 >= 1) ic2 = ic2 + 1
            i = 0
            DO WHILE(degenerat(iband2 - i) == 0)
               i = i + 1
            END DO
            ndb2 = degenerat(iband2 - i)
            ! only upper triangular part
            IF(symequivalent(ic2, ic1) .and. iband2 <= iband1) THEN
!            IF( ndb1 .ne. ndb2 ) call judft_error('symm_hf: failure symequivalent')
               nsest(iband1) = nsest(iband1) + 1
               indx_sest(nsest(iband1), iband1) = iband2
            END IF
         END DO
      END DO
      call timestop("calc bandcombos")

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

      call timestart("calc core repr")
      IF(hybdat%lmaxcd > fi%atoms%lmaxd) then
         call judft_error('symm_hf: The very impropable case that hybdat%lmaxcd > fi%atoms%lmaxd occurs')
      endif
      iatom = 0
      iatom0 = 0
      DO itype = 1, fi%atoms%ntype
         DO ieq = 1, fi%atoms%neq(itype)
            iatom = iatom + 1
            DO iop = 1, nsymop
               isym = psym(iop)
               IF(isym <= fi%sym%nop) THEN
                  iisym = isym
               ELSE
                  iisym = isym - fi%sym%nop
               END IF

               ratom = fi%hybinp%map(iatom, isym)
               rotkpt = matmul(rrot(:, :, isym), fi%kpts%bkf(:, nk))
               g = nint(rotkpt - fi%kpts%bkf(:, nk))

               cdum = exp(-2*pi*img*dot_product(rotkpt, fi%sym%tau(:, iisym)))* &
                      exp(2*pi*img*dot_product(g, fi%atoms%taual(:, ratom)))
            END DO
         END DO
         iatom0 = iatom0 + fi%atoms%neq(itype)
      END DO
      call timestop("calc core repr")

      CALL timestop("symm_hf")
   END SUBROUTINE symm_hf
END MODULE m_symm_hf
