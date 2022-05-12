MODULE m_checkolap
   use m_ylm
CONTAINS

   SUBROUTINE checkolap(atoms, hybdat, mpdata, hybinp, nkpti, kpts, fmpi, &
                        input, sym, noco, nococonv,   cell, lapw, jsp)
      USE m_util, ONLY: chr, sphbessel, harmonicsr
      use m_intgrf, only: intgrf, intgrf_init
      use m_calc_cmt
      USE m_constants
      USE m_types
      USE m_io_hybrid
      USE m_types_hybdat
      use m_calc_l_m_from_lm
#ifdef CPP_MPI
      use mpi 
#endif

      IMPLICIT NONE

      TYPE(t_hybdat), INTENT(IN)   :: hybdat

      TYPE(t_mpi), INTENT(IN)         :: fmpi
      TYPE(t_mpdata), intent(in)      :: mpdata
      TYPE(t_hybinp), INTENT(IN)      :: hybinp
      TYPE(t_input), INTENT(IN)       :: input
      TYPE(t_noco), INTENT(IN)        :: noco
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_sym), INTENT(IN)         :: sym
      TYPE(t_cell), INTENT(IN)        :: cell
      TYPE(t_kpts), INTENT(IN)        :: kpts
      TYPE(t_atoms), INTENT(IN)       :: atoms
       
      TYPE(t_lapw), INTENT(INOUT)     :: lapw

      ! - scalars -
      INTEGER, INTENT(IN)     :: jsp
      INTEGER, INTENT(IN)     ::  nkpti

      ! - local scalars -
      INTEGER                 ::  i, itype, iatom, ikpt, ineq, igpt, iband
      INTEGER                 ::  j, m
      INTEGER                 ::  l
      INTEGER                 :: lm, lm1
      INTEGER                 ::  n, nbasfcn

      REAL                    ::  rdum, rdum1
      REAL                    ::  qnorm

      COMPLEX                 ::  cexp, cdum, pre_fac
      COMPLEX, PARAMETER     ::  img = (0.0, 1.0)

      ! -local arrays -
      INTEGER                 ::  iarr(2), gpt(3), ierr
      INTEGER, ALLOCATABLE   ::  olapcv_loc(:, :, :, :, :)

      REAL                    ::  sphbes(0:atoms%lmaxd)
      REAL                    ::  q(3)
      REAL                    ::  integrand(atoms%jmtd)
      REAL                    ::  rarr(maxval(hybdat%nbands))
      REAL, ALLOCATABLE   ::  olapcb(:)
      REAL, ALLOCATABLE   :: olapcv_avg(:, :, :, :), olapcv_max(:, :, :, :)
      TYPE(t_mat), ALLOCATABLE :: z(:)

      COMPLEX                 ::  cmt(input%neig, hybdat%maxlmindx, atoms%nat, nkpti)
      COMPLEX                 ::  y((atoms%lmaxd + 1)**2)
      COMPLEX, ALLOCATABLE   ::  olapcv(:, :), c_phase(:)
      COMPLEX, ALLOCATABLE   ::  carr1(:, :), carr2(:, :), carr3(:, :)

      CHARACTER, PARAMETER    ::  lchar(0:38) = &
                                 (/'s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', &
                                   'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', &
                                   'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x'/)
      LOGICAL, parameter      ::  l_mism = .false.

      call timestart("checkolap")
      allocate(z(nkpti))
      DO ikpt = 1, nkpti
         CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, ikpt, cell)
         nbasfcn = MERGE(lapw%nv(1) + lapw%nv(2) + 2*atoms%nlotot, lapw%nv(1) + atoms%nlotot, noco%l_noco)
         call z(ikpt)%alloc(sym%invs, nbasfcn, input%neig)
      ENDDO

      IF(fmpi%irank == 0) WRITE(oUnit, '(//A)') '### checkolap ###'

      cmt = 0

      ! initialize gridf -> was done in eigen_HF_init
      !CALL intgrf_init(atoms%ntype,atoms%jmtd,atoms%jri,atoms%dx,atoms%rmsh,hybdat%gridf)

      ! read in cmt
      DO ikpt = 1, nkpti
         if(allocated(c_phase)) deallocate(c_phase)
         allocate(c_phase(hybdat%nbands(ikpt, jsp)))

         if(ikpt /= kpts%bkp(ikpt)) call juDFT_error("We should be reading the parent z-mat here!")
         call read_z(atoms, cell, hybdat, kpts, sym, noco, nococonv, input, ikpt, &
                     jsp, z(ikpt), c_phase=c_phase)
#ifdef CPP_MPI
         ! call timestart("Post read_z Barrier: checkolap")
         ! call MPI_Barrier(MPI_COMM_WORLD, ierr)
         ! call timestop("Post read_z Barrier: checkolap")
#endif
         call calc_cmt(atoms, cell, input, noco, nococonv, hybinp, hybdat, mpdata, kpts, &
                       sym,   z(kpts%bkp(ikpt)), jsp, ikpt, c_phase, &
                       cmt(:hybdat%nbands(ikpt,jsp), :, :, ikpt))
      END DO

      IF(fmpi%irank == 0) WRITE(oUnit, '(/A)') ' Overlap <core|core>'
      DO itype = 1, atoms%ntype
         IF(atoms%ntype > 1 .AND. fmpi%irank == 0) &
            WRITE(oUnit, '(A,I3)') ' Atom type', itype
         DO l = 0, hybdat%lmaxc(itype)
            DO i = 1, hybdat%nindxc(l, itype)
               IF(fmpi%irank == 0) &
                  WRITE(oUnit, '(1x,I1,A,2X)', advance='no') i + l, lchar(l)
               DO j = 1, i
                  integrand = hybdat%core1(:, i, l, itype)*hybdat%core1(:, j, l, itype) &
                              + hybdat%core2(:, i, l, itype)*hybdat%core2(:, j, l, itype)
                  IF(fmpi%irank == 0) WRITE(oUnit, '(F10.6)', advance='no') &
                     intgrf(integrand, atoms, itype, hybdat%gridf)
               END DO
               IF(fmpi%irank == 0) WRITE(oUnit, *)
            END DO
         END DO
      END DO

      IF(fmpi%irank == 0) WRITE(oUnit, '(/A)') ' Overlap <core|basis>'
      allocate(olapcb(maxval(mpdata%num_radfun_per_l)), olapcv(maxval(hybdat%nbands), nkpti), &
               olapcv_avg(-hybdat%lmaxcd:hybdat%lmaxcd, hybdat%maxindxc, 0:hybdat%lmaxcd, atoms%ntype), &
               olapcv_max(-hybdat%lmaxcd:hybdat%lmaxcd, hybdat%maxindxc, 0:hybdat%lmaxcd, atoms%ntype), &
               olapcv_loc(2, -hybdat%lmaxcd:hybdat%lmaxcd, hybdat%maxindxc, 0:hybdat%lmaxcd, atoms%ntype))

      DO itype = 1, atoms%ntype
         IF(atoms%ntype > 1 .AND. fmpi%irank == 0) &
            WRITE(oUnit, '(A,I3)') ' Atom type', itype
         DO l = 0, hybdat%lmaxc(itype)
            IF(l > atoms%lmax(itype)) EXIT ! very improbable case
            IF(fmpi%irank == 0) &
               WRITE(oUnit, "(9X,'u(',A,')',4X,'udot(',A,')',:,3X,'ulo(',A,"// &
                     "') ...')")(lchar(l), i=1, min(3, mpdata%num_radfun_per_l(l, itype)))
            DO i = 1, hybdat%nindxc(l, itype)
               IF(fmpi%irank == 0) &
                  WRITE(oUnit, '(1x,I1,A,2X)', advance='no') i + l, lchar(l)
               DO j = 1, mpdata%num_radfun_per_l(l, itype)

                  integrand = hybdat%core1(:, i, l, itype)*hybdat%bas1(:, j, l, itype) &
                              + hybdat%core2(:, i, l, itype)*hybdat%bas2(:, j, l, itype)

                  olapcb(j) = &
                     intgrf(integrand, atoms, itype, hybdat%gridf)

                  IF(fmpi%irank == 0) &
                     WRITE(oUnit, '(F10.6)', advance='no') olapcb(j)
               END DO

               lm = sum([(mpdata%num_radfun_per_l(j, itype)*(2*j + 1), j=0, l - 1)])
               iatom = sum(atoms%neq(1:itype - 1)) + 1 ! take first of group of equivalent atoms
               DO m = -l, l
                  olapcv = 0
                  DO j = 1, mpdata%num_radfun_per_l(l, itype)
                     lm = lm + 1
                     olapcv(:, :) = olapcv(:, :) + &
                                    olapcb(j)*cmt(:maxval(hybdat%nbands), lm, iatom, :nkpti)
                  END DO
                  rdum = sum(abs(olapcv(:, :))**2)
                  rdum1 = maxval(abs(olapcv(:, :)))
                  iarr = maxloc(abs(olapcv(:, :)))
                  olapcv_avg(m, i, l, itype) = &
                     sqrt(rdum/nkpti/sum(hybdat%nbands(:nkpti,jsp))*nkpti)
                  olapcv_max(m, i, l, itype) = rdum1
                  olapcv_loc(:, m, i, l, itype) = iarr
               END DO
               IF(fmpi%irank == 0) WRITE(oUnit, *)

            END DO
         END DO
      END DO

      IF(fmpi%irank == 0) THEN
         WRITE(oUnit, '(/A)') ' Average overlap <core|val>'
         DO itype = 1, atoms%ntype
            IF(atoms%ntype > 1) write(oUnit, '(A,I3)') ' Atom type', itype
            DO l = 0, hybdat%lmaxc(itype)
               DO i = 1, hybdat%nindxc(l, itype)
                  WRITE(oUnit, '(1x,I1,A,2X)', advance='no') i + l, lchar(l)
                  WRITE(oUnit, '('//chr(2*l + 1)//'F10.6)') &
                     olapcv_avg(-l:l, i, l, itype)
               END DO
            END DO
         END DO

         WRITE(oUnit, '(/A)') ' Maximum overlap <core|val> at (band/kpoint)'
         DO itype = 1, atoms%ntype
            IF(atoms%ntype > 1) write(oUnit, '(A,I3)') ' Atom type', itype
            DO l = 0, hybdat%lmaxc(itype)
               DO i = 1, hybdat%nindxc(l, itype)
                  WRITE(oUnit, '(1x,I1,A,2X)', advance='no') i + l, lchar(l)
                  WRITE(oUnit, '('//chr(2*l + 1)// &
                        '(F10.6,'' ('',I3.3,''/'',I4.3,'')''))') &
                     (olapcv_max(m, i, l, itype), &
                      olapcv_loc(:, m, i, l, itype), m=-l, l)
               END DO
            END DO
         END DO
      END IF ! fmpi%irank == 0

      deallocate(olapcb, olapcv, olapcv_avg, olapcv_max, olapcv_loc)

      IF(fmpi%irank == 0) WRITE(oUnit, '(/A)') ' Overlap <basis|basis>'

      DO itype = 1, atoms%ntype
         IF(atoms%ntype > 1 .AND. fmpi%irank == 0) &
            WRITE(oUnit, '(A,I3)') ' Atom type', itype
         DO l = 0, atoms%lmax(itype)
            DO i = 1, mpdata%num_radfun_per_l(l, itype)
               IF(fmpi%irank == 0) THEN
                  SELECT CASE(i)
                  CASE(1)
                     WRITE(oUnit, '(1x,''   u('',A,'')'')', advance='no') lchar(l)
                  CASE(2)
                     WRITE(oUnit, '(1x,''udot('',A,'')'')', advance='no') lchar(l)
                  CASE DEFAULT
                     WRITE(oUnit, '(1x,'' ulo('',A,'')'')', advance='no') lchar(l)
                  END SELECT
               END IF
               DO j = 1, i
                  integrand = hybdat%bas1(:, i, l, itype)*hybdat%bas1(:, j, l, itype) &
                              + hybdat%bas2(:, i, l, itype)*hybdat%bas2(:, j, l, itype)

                  IF(fmpi%irank == 0) WRITE(oUnit, '(F10.6)', advance='no') &
                     intgrf(integrand, atoms, itype, hybdat%gridf)
               END DO
               IF(fmpi%irank == 0) WRITE(oUnit, *)
            END DO
         END DO
      END DO

      IF(l_mism)then

         IF(fmpi%irank == 0) WRITE(oUnit, '(/A)') &
            'Mismatch of wave functions at the MT-sphere boundaries'
         allocate(carr1(maxval(hybdat%nbands),(atoms%lmaxd + 1)**2))
         allocate(carr2(maxval(hybdat%nbands),(atoms%lmaxd + 1)**2))
         allocate(carr3(maxval(hybdat%nbands),(atoms%lmaxd + 1)**2))

         ! create lock for race-condition in coulomb
         DO ikpt = 1, nkpti
            iatom = 0
            DO itype = 1, atoms%ntype
               DO ineq = 1, atoms%neq(itype)
                  iatom = iatom + 1            
                  carr1 = 0; carr2 = 0; carr3 = 0

                  ! calculate k1,k2,k3
                  CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, ikpt, cell)
                  call timestart("pw part")
                  ! PW part
                  !$OMP PARALLEL DO default(none) &
                  !$OMP private(igpt, gpt, cexp, q, qnorm, sphbes, y, pre_fac, lm, l, m, iband, cdum) &
                  !$OMP shared(lapw, jsp, atoms, kpts, iatom, ikpt, cell, itype, ineq, z, hybdat) &
                  !$OMP reduction(+:carr2)
                  DO igpt = 1, lapw%nv(jsp)
                     gpt = lapw%gvec(:, igpt, jsp)

                     cexp = exp(img*2*pi_const* &
                              dot_product(kpts%bkf(:, ikpt) + gpt, atoms%taual(:, iatom)))
                     q = matmul(kpts%bkf(:, ikpt) + gpt, cell%bmat)

                     qnorm = norm2(q)
                     call sphbessel(sphbes, atoms%rmt(itype)*qnorm, atoms%lmax(itype))

                     call ylm4(atoms%lmax(itype), q, y)
                     y = conjg(y)
                     
                     pre_fac = fpi_const / sqrt(cell%omtil) * cexp
                     if(z(1)%l_real) THEN
                        do lm = 1, (atoms%lmax(itype)+1)**2
                           call calc_l_m_from_lm(lm, l, m)
                           DO iband = 1, hybdat%nbands(ikpt,jsp)
                              cdum = pre_fac * ImagUnit**l * sphbes(l)
                              carr2(iband, lm) = carr2(iband, lm) + cdum*z(ikpt)%data_r(igpt, iband)*y(lm)
                           enddo
                        enddo
                     else
                        do lm = 1, (atoms%lmax(itype)+1)**2
                           call calc_l_m_from_lm(lm, l, m)
                           DO iband = 1, hybdat%nbands(ikpt,jsp)
                              cdum = pre_fac * ImagUnit**l * sphbes(l)
                              carr2(iband, lm) = carr2(iband, lm) + cdum*z(ikpt)%data_c(igpt, iband)*y(lm)
                           end DO
                        END DO
                     endif
                  enddo
                  !$OMP END PARALLEL DO
                  call timestop("pw part")

                  call timestart("MT-part")
                  ! MT
                  lm1 = 0
                  do lm = 1,(atoms%lmax(itype)+1)**2
                     call calc_l_m_from_lm(lm, l, m)
                     DO n = 1, mpdata%num_radfun_per_l(l, itype)
                        lm1 = lm1 + 1
                        rdum = hybdat%bas1(atoms%jri(itype), n, l, itype)/atoms%rmt(itype)
                        DO iband = 1, hybdat%nbands(ikpt,jsp)
                           carr3(iband, lm) = carr3(iband, lm) + cmt(iband, lm1, iatom, ikpt)*rdum
                        END DO
                     END DO
                  END DO
                  call timestop("MT-part")
                  carr1 = carr2 - carr3


                  rarr = 0
                  lm = 0
                  DO l = 0, atoms%lmax(itype)
                     DO m = -l, l
                        lm = lm + 1
                        rarr = rarr + abs(carr1(:, lm))**2
                     END DO
                  END DO
                  rarr = sqrt(rarr/(4*pi_const))

                  write (oUnit, '(I6,4X,F14.12,''  ('',F14.12,'')'')') ikpt,sum(rarr(:1)**2/hybdat%nbands(ikpt,jsp)),maxval(rarr(:1))
               END DO
            END DO
         END DO
      endif
      call timestop("checkolap")
   END SUBROUTINE checkolap

END MODULE m_checkolap
