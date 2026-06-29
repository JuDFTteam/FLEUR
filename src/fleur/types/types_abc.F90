!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_abc
#ifdef _OPENACC
   use cublas
#define CPP_ACC acc
#define CPP_OMP no_OMP_used
#define zgemm_acc cublaszgemm
#else
#define CPP_ACC No_acc_used
#define CPP_OMP OMP
#define zgemm_acc zgemm
#endif
   use m_judft
   IMPLICIT NONE

   PRIVATE

   
   TYPE t_abc
   !! A derived type for the AB+(C) coeffs
   !! 
   !! This type is used to store and calculate the "large AB" coefficients as 
   !! well as the coefficients for the local orbitals in a unified way.
     
      COMPLEX, ALLOCATABLE :: cof(:, :, :, :) !(nu,lm,iOrd,iAtom)
      !! This array stores coefficients with the following dimensions:
      !! - `nu`: Number of states or bands.
      !! - `lm`: Angular momentum quantum numbers.
      !! - `iOrd`: index of the radial function, 
      !!                   i.e. 1="u", the former "A" coefficient
      !!                        2="\dot{u}", the former "B" coefficient
      !!                        3="LO", the former "C" coefficient     
      !!                        4="LO", the former "C" coefficient      
      !! - `iAtom`: Atom index.
      Integer,allocatable :: n_r(:)
   CONTAINS
      PROCEDURE, PASS :: init => abc_init !!

      PROCEDURE, PASS :: calc_abc 
      !! Calculates the abc coefficients, previously called "abcof".
      
      PROCEDURE, PASS :: calc_force_abc
      !! Calculates the force-related abc coefficients, previously called "to_pulay".      
      
      PROCEDURE, PASS :: rotate
      !! Rotates the coefficients or data within the t_abc type.
      !> \brief Rotates coefficients to the representative atom.
      !! (Currently commented out in the code.)
      !! PROCEDURE, PASS :: rotate_to_rep_atom => rotate_eigveccoeffs_to_rep_atom

      procedure, PASS :: rot_to_unrotated

   END TYPE t_abc


   PUBLIC t_abc

CONTAINS

   SUBROUTINE abc_init(this, input, atoms, n_r_in , noccbd, itype)

      USE m_types_atoms
      USE m_types_input

      IMPLICIT NONE

      CLASS(t_abc), INTENT(INOUT) :: this
      INTEGER, INTENT(IN)   :: n_r_in(0:)
      TYPE(t_atoms), INTENT(IN)    :: atoms
      TYPE(t_input), INTENT(IN)    :: input

      INTEGER, INTENT(IN)    :: itype, noccbd

      IF (ALLOCATED(this%cof)) DEALLOCATE (this%cof,this%n_r)
      allocate(this%n_r(0:size(n_r_in)-1))
      this%n_r(0:)=n_r_in(0:)
      ALLOCATE (this%cof(noccbd, 0:atoms%lmax(itype)*(atoms%lmax(itype) + 2), maxval(this%n_r), atoms%neq(itype)))
      this%cof = CMPLX(0.0, 0.0)

   END SUBROUTINE abc_init

   subroutine calc_abc(this, input, atoms, sym, cell, lapw, ne, usdus, &
                       noco, nococonv, jspin, itype, zMat)
      USE m_juDFT
      USE m_types_atoms
      USE m_types_input
      USE m_types_sym
      USE m_types_cell
      USE m_types_lapw
      USE m_types_usdus
      USE m_types_noco
      USE m_types_nococonv
      USE m_types_enpara
      USE m_constants
      USE m_ylm
      USE m_setabc1lo
      USE m_hsmt_fjgj
      USE m_hsmt_ab
      USE m_types_mat

      IMPLICIT NONE
      CLASS(t_abc), INTENT(INOUT) :: this
      TYPE(t_input), INTENT(IN)             :: input
      TYPE(t_usdus), INTENT(IN)             :: usdus
      TYPE(t_lapw), INTENT(IN)              :: lapw

      TYPE(t_noco), INTENT(IN)              :: noco
      TYPE(t_nococonv), INTENT(IN)          :: nococonv
      TYPE(t_sym), INTENT(IN)               :: sym
      TYPE(t_cell), INTENT(IN)              :: cell
      TYPE(t_atoms), INTENT(IN)             :: atoms
      TYPE(t_mat), INTENT(IN)               :: zMat

! scalar arguments
      INTEGER, INTENT(IN)        :: ne
      INTEGER, INTENT(IN)        :: jspin, itype

! Local objects
      TYPE(t_fjgj) :: fjgj

! Local scalars
      INTEGER :: i, iLAPW, l, lm, nap, jAtom, lmp, m, nkvec, iAtom, acof_size, iAtom_l, jatom_l
      INTEGER :: inv_f, ie, ilo, kspin, iintsp, nintsp, nvmax, lo, inap, abSize, n_l(0:atoms%lmaxd), nbasf
      REAL    :: tmk, qss(3), s2h
      COMPLEX :: phase, c_1, c_2, term1, ctmp
      LOGICAL ::  l_useinversionsym

! Local arrays
      REAL    :: fg(3), fk(3), fkp(3), fkr(3)
      REAL    :: alo1(atoms%nlod, input%jspins), blo1(atoms%nlod, input%jspins)
      REAL    :: clo1(atoms%nlod, input%jspins)
      COMPLEX :: ylm((atoms%lmaxd + 1)**2)
      COMPLEX :: ccchi(2, 2)
      REAL, ALLOCATABLE :: realCoeffs(:, :), imagCoeffs(:, :), workTrans_r(:, :)
      REAL, ALLOCATABLE :: fgpl(:, :)
      COMPLEX, ALLOCATABLE :: s2h_e(:, :)
      COMPLEX, ALLOCATABLE :: work_c(:, :), workTrans_c(:, :), workTrans_cf(:, :),work_lo(:)
      COMPLEX, ALLOCATABLE :: abCoeffs(:, :)
      COMPLEX, ALLOCATABLE :: abTemp(:, :)
      COMPLEX, ALLOCATABLE :: helpMat_force(:, :)

      CALL timestart("abcof")

! Checks
      IF (zmat%l_real) THEN
         IF (noco%l_noco) CALL judft_bug("BUG in abcof, l_noco but real?")
      END IF

! Allocations
      CALL fjgj%alloc(MAXVAL(lapw%nv), atoms%lmaxd, jspin, noco)
      ! abCoeffs is allocated (and filled) inside hsmt_ab.
      ALLOCATE (abTemp(SIZE(this%cof, 1), 0:2*SIZE(this%cof, 2) - 1))
      ALLOCATE (fgpl(3, MAXVAL(lapw%nv)))
      ALLOCATE (work_c(MAXVAL(lapw%nv), ne))
      ALLOCATE (work_lo(ne))

! Initializations
      acof_size = size(this%cof, 1)
!$acc enter data create(abTemp,fjgj,fjgj%fj,fjgj%gj,work_c)
      

!Use inversion symmetry explicitely
      l_useinversionsym = any(sym%invsat == 2)!.and.(.not.noco%l_soc).and.(.not.present(nat_start))

      CALL timestart("fjgj coefficients")
      CALL fjgj%calculate(input, atoms, cell, lapw, noco, usdus, iType, jspin)
!$acc update device (fjgj%fj,fjgj%gj)
      CALL timestop("fjgj coefficients")

      CALL setabc1lo(atoms, iType, usdus, jspin, alo1, blo1, clo1)

! generate the spinors (chi)
      IF (noco%l_noco) ccchi = conjg(nococonv%umat(itype))
      n_l = 2

! loop over atoms
      DO iAtom_l = 1, atoms%neq(itype)
         iAtom = iAtom_l - 1 + atoms%firstAtom(itype)
         if (sym%invsat(iatom) == 2 .and. l_useinversionsym) cycle
         nintsp = 1
         IF (noco%l_ss) nintsp = 2
! loop over the interstitial spin
         DO iintsp = 1, nintsp
            nvmax = lapw%nv(jspin)
            IF (noco%l_ss) nvmax = lapw%nv(iintsp)
            qss = MERGE(-1.0, 1.0, iintsp .EQ. 1)*nococonv%qss/2.0

! Filling of work array (modified zMat)
            CALL timestart("fill work array")
            IF (noco%l_noco) THEN
            IF (noco%l_ss) THEN
               !$acc kernels copyin(zmat,zMat%data_c,ccchi,atoms,lapw,lapw%nv) present(work_c)default(none)
               ! the coefficients of the spin-down basis functions are
               ! stored in the second half of the eigenvector
               kspin = (iintsp - 1)*(lapw%nv(1) + atoms%nlotot)
               work_c(:nvmax, :) = ccchi(iintsp, jspin)*zMat%data_c(kspin + 1:kspin + nvmax, :ne)
               !$acc end kernels
            ELSE
               ! perform sum over the two interstitial spin directions
               ! and take into account the spin boundary conditions
               ! (jspin counts the local spin directions inside each MT)
               !$acc kernels copyin(atoms,zMat,zMat%data_c,ccchi,lapw,lapw%nv) present(work_c) default(none)
               kspin = lapw%nv(1) + atoms%nlotot
               work_c(:nvmax, :) = ccchi(1, jspin)*zMat%data_c(:nvmax, :ne) + &
                                   ccchi(2, jspin)*zMat%data_c(kspin + 1:kspin + nvmax, :ne)
               !$acc end kernels
            END IF
            ELSE
            IF (zmat%l_real) THEN
               !$CPP_OMP PARALLEL DO default(shared) private(i)
               !$acc kernels copyin(zmat,zMat%data_r)present(work_c)default(none)
               DO i = 1, ne
#ifdef _OPENACC
                  work_c(:nvmax, i) = zmat%data_r(:nvmax, i)
#else
                  work_c(:nvmax, i) = 0.0
                  CALL dcopy(nvmax, zMat%data_r(:, i), 1, work_c(:, i), 2)
#endif
               END DO
               !$acc end kernels
               !$CPP_OMP END PARALLEL DO
            ELSE
               !$CPP_OMP PARALLEL DO default(shared) private(i)
               !$acc kernels copyin(zMat,zMat%data_c)present(work_c) default(none)
               DO i = 1, ne
#ifdef _OPENACC
                  work_c(:nvmax, i) = zmat%data_c(:nvmax, i)
#else
                  CALL zcopy(nvmax, zMat%data_c(:, i), 1, work_c(:, i), 1)
#endif
               END DO
               !$acc end kernels
               !$CPP_OMP END PARALLEL DO
            END IF
            END IF

            CALL timestop("fill work array")

! Calculation of a, b coefficients for LAPW basis functions
            CALL timestart("hsmt_ab")
!!$acc data copyin(fjgj,fjgj%fj,fjgj%gj) copyout(abcoeffs)
            CALL hsmt_ab(sym, atoms, noco, nococonv, jspin, iintsp, iType, iAtom, cell, lapw, fjgj, abCoeffs, abSize, .FALSE.)
!!$acc end data
            abSize = abSize/2
            CALL timestop("hsmt_ab")

! Obtaining A, B coefficients for eigenfunctions
            CALL timestart("gemm")

! variant with zgemm

!$acc host_data use_device(work_c,abCoeffs,abTemp)
CALL zgemm_acc("T","T",ne,2*abSize,nvmax,CMPLX(1.0,0.0),work_c,MAXVAL(lapw%nv),abCoeffs,SIZE(abCoeffs,1),CMPLX(0.0,0.0),abTemp,acof_size)
!$acc end host_data
!$acc update self(abTemp)
!stop "DEBUG"
!$OMP PARALLEL DO default(shared) private(i,lm) collapse(2)
            DO lm = 0, absize - 1
            DO i = 1, ne
               this%cof(i, lm, 1, iAtom_l) = this%cof(i, lm, 1, iAtom_l) + abTemp(i, lm)
               this%cof(i, lm, 2, iAtom_l) = this%cof(i, lm, 2, iAtom_l) + abTemp(i, absize + lm)
            END DO
            END DO
!$OMP END PARALLEL DO

            CALL timestop("gemm")
            ! abCoeffs is (re)allocated per call inside hsmt_ab; release the
            ! device copy it created and the host array before the next call.
            !$acc exit data delete(abCoeffs)
            DEALLOCATE(abCoeffs)

            CALL timestart("local orbitals")
! Treatment of local orbitals
!!$acc data copyin(alo1,blo1,clo1,ccchi)create(ylm)
            n_l = 2
            DO lo = 1, atoms%nlo(iType)
               l = atoms%llo(lo, itype)
               n_l(l) = n_l(l) + 1
               DO nkvec = 1, lapw%nkvec(lo, iAtom)
                  iLAPW = lapw%kvec(nkvec, lo, iAtom)
                  fg(:) = MERGE(lapw%gvec(:, iLAPW, iintsp), lapw%gvec(:, iLAPW, jspin), noco%l_ss) + qss + lapw%qPhon
                  fk = lapw%bkpt + fg(:)
                  tmk = tpi_const*DOT_PRODUCT(fk(:), atoms%taual(:, iAtom))
                  phase = CMPLX(COS(tmk), SIN(tmk))

                  nap = sym%ngopr(iAtom)
                  inap = sym%invtab(nap)
                  fkr = MATMUL(TRANSPOSE(sym%mrot(:, :, inap)), fk(:))

                  fkp = MATMUL(fkr, cell%bmat)

                  CALL ylm4(atoms%lmax(iType), fkp, ylm)
    !!$acc update device(ylm)
                  ! Code from previous abclocdn
                  term1 = 2*tpi_const/SQRT(cell%omtil)*((atoms%rmt(itype)**2)/2)*phase
                  !---> the whole program is in hartree units, therefore 1/wronskian is
                  !---> (rmt**2)/2. the factor i**l, which usually appears in the a, b
                  !---> and c coefficients, is included in the t-matrices. thus, it does
                  !---> not show up in the formula above.
                  nbasf = lapw%nv(iintsp) + lapw%index_lo(lo, iatom) + nkvec
                  if (noco%l_noco) Then
                     if (noco%l_ss) THEN
                        work_lo(:ne) = ccchi(iintsp, jspin)*zMat%data_c((iintsp - 1)*(lapw%nv(1) + atoms%nlotot) + nbasf, :ne)
                     else
                        work_lo(:ne) = ccchi(1, jspin)*zMat%data_c(nbasf, :ne) + &
                                         ccchi(2, jspin)*zMat%data_c(lapw%nv(1) + atoms%nlotot + nbasf, :ne)
                     END IF
                  ELSE
                     if (zmat%l_real) Then
                        work_lo(:ne) = zmat%data_r(nbasf, :ne)
                     else
                        work_lo(:ne) = zmat%data_c(nbasf, :ne)
                     end if
                  end if

    !!$acc kernels default(none) present(acof,bcof,ccof,alo1,blo1,clo1,ccchi,ylm)create(ctmp) &
    !!$acc copyin(work,na,term1,l,ne,ll1,noco)
    !!$acc loop seq private(i,m,lm,ctmp,na2,lmp)

       !!$acc loop seq
                  DO m = -l, l
                     lm = l*(l + 1) + m
                     DO i = 1, ne
                        ctmp = term1*conjg(ylm(lm + 1))*work_lo(i)
                        this%cof(i, lm, 1, iatom_l) = this%cof(i, lm, 1, iatom_l) + ctmp*alo1(lo, jspin)
                        this%cof(i, lm, 2, iatom_l) = this%cof(i, lm, 2, iatom_l) + ctmp*blo1(lo, jspin)
                        this%cof(i, lm, n_l(l), iatom_l) = this%cof(i, lm, n_l(l), iatom_l) + ctmp*clo1(lo, jspin)
                     END DO
          !!$acc end loop
                  END DO
    !!$acc end loop
    !!$acc end kernels
               END DO
            END DO ! loop over LOs
!!$acc end data
            CALL timestop("local orbitals")
         END DO ! loop over interstitial spin
      END DO ! loop over atoms
!$acc exit data delete(abTemp,fjgj%fj,fjgj%gj,work_c)
!$acc exit data delete(fjgj)
      DEALLOCATE (work_c)

! Treatment of atoms inversion symmetric to others
      IF (l_useinversionsym) THEN
!Comment on SOC case:
!
!                           -p,n       (l+m)   p,n  *
! Usually, we exploit that A     = (-1)      (A    )  if p and -p are the positions
!                           l,m                l,-m
! of two atoms related by inversion symmetry and the coefficients are considered to
! be in the local frame of the representative atom. This is possible, if z is real.
! After SOC, however, the eigenvectors z are complex and this is no longer possible
! so the z has to enter, not z*. This is done within the k-loop.
!                                    -p,n       m   p,n  *
! When called from hsohelp, we need A     = (-1)  (A    ) because we don't have to
!                                     l,m           l,-m                    l
! rotate, but in the sums in hsoham only products A*  A   enter and the (-1) cancels.
!                                                  lm  lm
         DO iAtom_l = 1, atoms%neq(itype)
            iatom = iatom_l - 1 + atoms%firstAtom(itype)
            IF (sym%invsat(iAtom) .EQ. 1) THEN
               jAtom = sym%invsatnr(iAtom)
               jatom_l = jatom - atoms%firstAtom(itype) + 1
               DO l = 0, atoms%lmax(iType)
                  DO m = -l, l
                     lm = l*(l + 1) + m
                     lmp = l*(l + 1) - m
                     inv_f = (-1)**(m + l)
                     this%cof(:ne, lm, :, jatom_l) = inv_f*CONJG(this%cof(:ne, lmp, :, iatom_l))
                  END DO
               END DO

            END IF
         END DO
      END IF

      CALL timestop("abcof")

   end subroutine calc_abc

   subroutine calc_force_abc(this, input, atoms, sym, cell, lapw, ne, usdus, &
                             noco, nococonv, jspin, itype, zMat,eig,force)
      USE m_juDFT
      USE m_types_atoms
      USE m_types_input
      USE m_types_force
      USE m_types_sym
      USE m_types_cell
      USE m_types_lapw
      USE m_types_usdus
      USE m_types_noco
      USE m_types_nococonv
      USE m_types_enpara
      USE m_constants
      USE m_ylm
      USE m_setabc1lo
      USE m_hsmt_fjgj
      USE m_hsmt_ab
      USE m_types_mat

      IMPLICIT NONE
      CLASS(t_abc), INTENT(INOUT) :: this
      TYPE(t_input), INTENT(IN)             :: input
      TYPE(t_usdus), INTENT(IN)             :: usdus
      TYPE(t_lapw), INTENT(IN)              :: lapw

      TYPE(t_noco), INTENT(IN)              :: noco
      TYPE(t_nococonv), INTENT(IN)          :: nococonv
      TYPE(t_sym), INTENT(IN)               :: sym
      TYPE(t_cell), INTENT(IN)              :: cell
      TYPE(t_atoms), INTENT(IN)             :: atoms
      TYPE(t_mat), INTENT(IN)               :: zMat
      REAL,intent(in)                       :: eig(:)
      TYPE(t_force), INTENT(INOUT)          :: force

! scalar arguments
      INTEGER, INTENT(IN)        :: ne
      INTEGER, INTENT(IN)        :: jspin, itype

! Local objects
      TYPE(t_fjgj) :: fjgj

! Local scalars
      INTEGER :: i, iLAPW, l, lm, nap, jAtom, lmp, m, nkvec, iAtom, acof_size, iAtom_l, jatom_l,j
      INTEGER :: inv_f, ie, ilo, kspin, iintsp, nintsp, nvmax, lo, inap, abSize, n_l(0:atoms%lmaxd), nbasf
      REAL    :: tmk, qss(3), s2h
      COMPLEX :: phase, c_1, c_2, term1, ctmp
      LOGICAL ::  l_useinversionsym

! Local arrays
      REAL    :: fg(3), fk(3), fkp(3), fkr(3), fgr(3), fgp(3)
      REAL    :: alo1(atoms%nlod, input%jspins), blo1(atoms%nlod, input%jspins)
      REAL    :: clo1(atoms%nlod, input%jspins)
      COMPLEX :: ylm((atoms%lmaxd + 1)**2)
      COMPLEX :: ccchi(2, 2)
      REAL, ALLOCATABLE :: realCoeffs(:, :), imagCoeffs(:, :), workTrans_r(:, :)
      REAL, ALLOCATABLE :: fgpl(:, :)
      COMPLEX, ALLOCATABLE :: s2h_e(:, :)
      COMPLEX, ALLOCATABLE :: work_c(:, :), workTrans_c(:, :), workTrans_cf(:, :)
      COMPLEX, ALLOCATABLE :: abCoeffs(:, :)
      COMPLEX, ALLOCATABLE :: abTemp(:, :)
      COMPLEX, ALLOCATABLE :: helpMat_force(:, :)

      CALL timestart("abcof")
      force%acoflo  = CMPLX(0.0,0.0)
      force%bcoflo  = CMPLX(0.0,0.0)
      force%e1cof   = CMPLX(0.0,0.0)
      force%e2cof   = CMPLX(0.0,0.0)
      force%aveccof = CMPLX(0.0,0.0)
      force%bveccof = CMPLX(0.0,0.0)
      force%cveccof = CMPLX(0.0,0.0)
      ALLOCATE(helpMat_force(ne,atoms%lmaxd*(atoms%lmaxd+2)+1))
      ALLOCATE(workTrans_cf(ne,MAXVAL(lapw%nv)))
      ALLOCATE(s2h_e(ne,MAXVAL(lapw%nv)))
! Checks
      IF (zmat%l_real) THEN
         IF (noco%l_noco) CALL judft_bug("BUG in abcof, l_noco but real?")
      END IF

! Allocations
      CALL fjgj%alloc(MAXVAL(lapw%nv), atoms%lmaxd, jspin, noco)
      ! abCoeffs is allocated (and filled) inside hsmt_ab.
      ALLOCATE (abTemp(SIZE(this%cof, 1), 0:2*SIZE(this%cof, 2) - 1))
      ALLOCATE (fgpl(3, MAXVAL(lapw%nv)))
      ALLOCATE (work_c(MAXVAL(lapw%nv), ne))

! Initializations
      acof_size = size(this%cof, 1)
!$acc enter data create(abTemp,fjgj,fjgj%fj,fjgj%gj,work_c)


!Use inversion symmetry explicitely
      l_useinversionsym = any(sym%invsat == 2)!.and.(.not.noco%l_soc).and.(.not.present(nat_start))

      CALL timestart("fjgj coefficients")
      CALL fjgj%calculate(input, atoms, cell, lapw, noco, usdus, iType, jspin)
!$acc update device (fjgj%fj,fjgj%gj)
      CALL timestop("fjgj coefficients")

      CALL setabc1lo(atoms, iType, usdus, jspin, alo1, blo1, clo1)

      ! generate the spinors (chi)
      IF (noco%l_noco) ccchi = conjg(nococonv%umat(itype))
      n_l = 2

! loop over atoms
      DO iAtom_l = 1, atoms%neq(itype)
         iAtom = iAtom_l - 1 + atoms%firstAtom(itype)
         if (sym%invsat(iatom) == 2 .and. l_useinversionsym) cycle
         nintsp = 1
         IF (noco%l_ss) nintsp = 2
! loop over the interstitial spin
         DO iintsp = 1, nintsp
            nvmax = lapw%nv(jspin)
            IF (noco%l_ss) nvmax = lapw%nv(iintsp)
            qss = MERGE(-1.0, 1.0, iintsp .EQ. 1)*nococonv%qss/2.0

            call fill_work_array(zmat,noco,atoms,lapw,ccchi,iintsp,nvmax,jspin,ne,work_c)

! Calculation of a, b coefficients for LAPW basis functions
            CALL timestart("hsmt_ab")
!!$acc data copyin(fjgj,fjgj%fj,fjgj%gj) copyout(abcoeffs)
            CALL hsmt_ab(sym, atoms, noco, nococonv, jspin, iintsp, iType, iAtom, cell, lapw, fjgj, abCoeffs, abSize, .FALSE.)
!!$acc end data
            abSize = abSize/2
            CALL timestop("hsmt_ab")
! Force contributions

            CALL timestart("force contributions")
            DO iLAPW = 1, nvmax
               fg(:) = MERGE(lapw%gvec(:, iLAPW, iintsp), lapw%gvec(:, iLAPW, jspin), noco%l_ss) + qss
               fk = lapw%bkpt + fg(:)
               s2h = 0.5*DOT_PRODUCT(fk, MATMUL(cell%bbmat, fk))
               s2h_e(:ne, iLAPW) = (s2h - eig(:ne))*work_c(iLAPW, :ne)
               nap = sym%ngopr(iAtom)
               inap = sym%invtab(nap)
               fgr = MATMUL(TRANSPOSE(sym%mrot(:, :, inap)), fg(:))
               fgpl(:, iLAPW) = MATMUL(fgr, cell%bmat)
            END DO
      
            workTrans_cf = CMPLX(0.0, 0.0)
            helpMat_force = CMPLX(0.0, 0.0)
      
            CALL zgemm("N", "T", ne, abSize, nvmax, CMPLX(1.0, 0.0), s2h_e, ne, abCoeffs, size(abcoeffs, 1), &
                       CMPLX(1.0, 0.0), force%e1cof(:, :, iAtom), ne)
            CALL zgemm("N", "T", ne, abSize, nvmax, CMPLX(1.0, 0.0), s2h_e, ne, abCoeffs(1 + abSize, 1), &
                       size(abcoeffs, 1), CMPLX(1.0, 0.0), force%e2cof(:, :, iAtom), ne)
            DO i = 1, 3
               DO iLAPW = 1, nvmax
                  workTrans_cf(:, iLAPW) = work_c(iLAPW, :)*fgpl(i, iLAPW)
               END DO
      
               CALL zgemm("N", "T", ne, abSize, nvmax, CMPLX(1.0, 0.0), workTrans_cf, ne, &
                          abCoeffs, size(abCoeffs, 1), CMPLX(0.0, 0.0), helpMat_force, ne)
               force%aveccof(i, :, :, iAtom) = force%aveccof(i, :, :, iAtom) + helpMat_force(:, :)
               CALL zgemm("N", "T", ne, abSize, nvmax, CMPLX(1.0, 0.0), workTrans_cf, ne, &
                          abCoeffs(1 + abSize, 1), size(abcoeffs, 1), CMPLX(0.0, 0.0), helpMat_force, ne)
               force%bveccof(i, :, :, iAtom) = force%bveccof(i, :, :, iAtom) + helpMat_force(:, :)
            END DO
            CALL timestop("force contributions")
            ! abCoeffs is (re)allocated per call inside hsmt_ab; release the
            ! device copy it created and the host array before the next call.
            !$acc exit data delete(abCoeffs)
            DEALLOCATE(abCoeffs)


            CALL timestart("local orbitals")
! Treatment of local orbitals
!!$acc data copyin(alo1,blo1,clo1,ccchi)create(ylm)
            n_l = 2
            DO lo = 1, atoms%nlo(iType)
               l = atoms%llo(lo, itype)
               n_l(l) = n_l(l) + 1
               DO nkvec = 1, lapw%nkvec(lo, iAtom)
                  iLAPW = lapw%kvec(nkvec, lo, iAtom)
                  fg(:) = MERGE(lapw%gvec(:, iLAPW, iintsp), lapw%gvec(:, iLAPW, jspin), noco%l_ss) + qss + lapw%qPhon
                  fk = lapw%bkpt + fg(:)
                  tmk = tpi_const*DOT_PRODUCT(fk(:), atoms%taual(:, iAtom))
                  phase = CMPLX(COS(tmk), SIN(tmk))

                  nap = sym%ngopr(iAtom)
                  inap = sym%invtab(nap)
                  fkr = MATMUL(TRANSPOSE(sym%mrot(:,:,inap)),fk(:))
                  fgr = MATMUL(TRANSPOSE(sym%mrot(:,:,inap)),fg(:))
                  
                  fkp = MATMUL(fkr,cell%bmat)
                  fgp = MATMUL(fgr,cell%bmat)
                  CALL ylm4(atoms%lmax(iType), fkp, ylm)
 !!$acc update device(ylm)
! Code from previous abclocdn
                  term1 = 2*tpi_const/SQRT(cell%omtil)*((atoms%rmt(itype)**2)/2)*phase
!---> the whole program is in hartree units, therefore 1/wronskian is
!---> (rmt**2)/2. the factor i**l, which usually appears in the a, b
!---> and c coefficients, is included in the t-matrices. thus, it does
!---> not show up in the formula above.
                  nbasf = lapw%nv(iintsp) + lapw%index_lo(lo, iatom) + nkvec
                  if (noco%l_noco) Then
                     if (noco%l_ss) THEN
                        work_c(:ne, 1) = ccchi(iintsp, jspin)*zMat%data_c((iintsp - 1)*(lapw%nv(1) + atoms%nlotot) + nbasf, :ne)
                     else
                        work_c(:ne, 1) = ccchi(1, jspin)*zMat%data_c(nbasf, :ne) + &
                                         ccchi(2, jspin)*zMat%data_c(lapw%nv(1) + atoms%nlotot + nbasf, :ne)
                     END IF
                  ELSE
                     if (zmat%l_real) Then
                        work_c(:ne, 1) = zmat%data_r(nbasf, :ne)
                     else
                        work_c(:ne, 1) = zmat%data_c(nbasf, :ne)
                     end if
                  end if

 !!$acc kernels default(none) present(acof,bcof,ccof,alo1,blo1,clo1,ccchi,ylm)create(ctmp) &
 !!$acc copyin(work,na,term1,l,ne,ll1,noco)
 !!$acc loop seq private(i,m,lm,ctmp,na2,lmp)

                     DO i = 1,ne
                       DO m = -l,l
                         lm = l*(l+1) + m
                         ctmp=term1*conjg(ylm(lm+1))*work_c(i,1)
                         force%acoflo(m,i,lo,iatom) = force%acoflo(m,i,lo,iatom) + ctmp*alo1(lo,jspin)
                         force%bcoflo(m,i,lo,iatom) = force%bcoflo(m,i,lo,iatom) + ctmp*blo1(lo,jspin)
                         DO j = 1,3
                          force%aveccof(j,i,lm,iatom)   = force%aveccof(j,i,lm,iatom)   + fgp(j)*ctmp*alo1(lo,jspin)
                          force%bveccof(j,i,lm,iatom)   = force%bveccof(j,i,lm,iatom)   + fgp(j)*ctmp*blo1(lo,jspin)
                          force%cveccof(j,m,i,lo,iatom) = force%cveccof(j,m,i,lo,iatom) + fgp(j)*ctmp*clo1(lo,jspin)
                        END DO
                       END DO
                     END DO

 !!$acc end loop
 !!$acc end kernels
               END DO
            END DO ! loop over LOs
!!$acc end data
            CALL timestop("local orbitals")
         END DO ! loop over interstitial spin
      END DO ! loop over atoms
!$acc exit data delete(abTemp,fjgj%fj,fjgj%gj,work_c)
!$acc exit data delete(fjgj)
      DEALLOCATE (work_c)


! Treatment of atoms inversion symmetric to others
      IF (l_useinversionsym) THEN
!Comment on SOC case:
!
!                           -p,n       (l+m)   p,n  *
! Usually, we exploit that A     = (-1)      (A    )  if p and -p are the positions
!                           l,m                l,-m
! of two atoms related by inversion symmetry and the coefficients are considered to
! be in the local frame of the representative atom. This is possible, if z is real.
! After SOC, however, the eigenvectors z are complex and this is no longer possible
! so the z has to enter, not z*. This is done within the k-loop.
!                                    -p,n       m   p,n  *
! When called from hsohelp, we need A     = (-1)  (A    ) because we don't have to
!                                     l,m           l,-m                    l
! rotate, but in the sums in hsoham only products A*  A   enter and the (-1) cancels.
!                                                  lm  lm
         DO iAtom = atoms%firstatom(itype),atoms%firstAtom(itype)+atoms%neq(itype)-1
            IF (sym%invsat(iAtom) .EQ. 1) THEN
               jAtom = sym%invsatnr(iAtom)
               DO ilo = 1, atoms%nlo(iType)
                  l = atoms%llo(ilo, iType)
                  DO m = -l, l
                     inv_f = (-1)**(m + l)
                        force%acoflo(m, :, ilo, jatom) = inv_f*CONJG(force%acoflo(-m, :, ilo, iatom))
                        force%bcoflo(m, :, ilo, jatom) = inv_f*CONJG(force%bcoflo(-m, :, ilo, iatom))
                        force%cveccof(:, m, :, ilo, jatom) = -inv_f*CONJG(force%cveccof(:, -m, :, ilo, iatom))
                  END DO
               END DO
               DO l = 0, atoms%lmax(iType)
                  DO m = -l, l
                     lm = l*(l + 1) + m
                     lmp = l*(l + 1) - m
                     inv_f = (-1)**(m + l)
                        force%e1cof(:ne, lm, jatom) = inv_f*CONJG(force%e1cof(:ne, lmp, iatom))
                        force%e2cof(:ne, lm, jatom) = inv_f*CONJG(force%e2cof(:ne, lmp, iatom))
                        force%aveccof(:, :ne, lm, jatom) = -inv_f*CONJG(force%aveccof(:, :ne, lmp, iatom))
                        force%bveccof(:, :ne, lm, jatom) = -inv_f*CONJG(force%bveccof(:, :ne, lmp, iatom))
                   
                  END DO
               END DO                   
            END IF
         END DO
      END IF

      CALL timestop("abcof")

   end subroutine calc_force_abc

   function rotate(abc, alpha, beta, gamma, lmax) result(abc_rot)
      USE m_dwigner

      IMPLICIT NONE
      class(t_abc), INTENT(IN)           :: abc
      real, intent(in)                   :: alpha, beta, gamma
      INTEGER, INTENT(IN)                 :: lmax

      TYPE(t_abc) :: abc_rot

      !     ..
      !     .. Local Scalars ..
      INTEGER n, j, l, i
      REAL amx(3, 3, 1), imx(3, 3)
      COMPLEX d_wgn(-lmax:lmax, -lmax:lmax, 1:lmax, 1)

      abc_rot = abc

      CALL euler(alpha, beta, gamma, amx)
      imx(:, :) = 0.; imx(1, 1) = 1.; imx(2, 2) = 1.; imx(3, 3) = 1.
      CALL d_wigner(1, amx, imx, lmax, d_wgn)

      DO n = 1, size(abc%cof, 4)
         DO j = 1, size(abc%cof, 3)
            DO l = 1, lmax
               DO i = 1, size(abc%cof, 1)
                  abc_rot%cof(i, l**2:l*(l + 2), j, n) = MATMUL(CONJG(d_wgn(-l:l, -l:l, l, 1)), &
                                                                  abc%cof(i, l**2:l*(l + 2), j, n))
               end do
            end do
         END DO
      END DO

   END function rotate

   subroutine rot_to_unrotated(abc,hybinp, atoms,  sym,itype)
!     ***************************************************************
!     * This routine transforms a/b/cof which are given wrt rotated *
!     * MT functions (according to invsat/ngopr) into a/b/cof wrt   *
!     * unrotated MT functions. Needed for GW calculations.         *
!     *                                                             *
!     * Christoph Friedrich Mar/2005                                *
!     ***************************************************************
      USE m_types_hybinp
      USE m_types_sym
      USE m_types_atoms
      USE m_juDFT
      IMPLICIT NONE
      CLASS(t_abc), INTENT(INOUT) :: abc
      TYPE(t_hybinp), INTENT(IN) :: hybinp
      TYPE(t_sym), INTENT(IN)    :: sym
      TYPE(t_atoms), INTENT(IN)  :: atoms
      INTEGER, INTENT(IN) :: itype

      INTEGER na, iatom, iop,  i, l, ifac,j
      call timestart("hyb_abcrot")
      IF (.NOT. ALLOCATED(hybinp%d_wgn2)) THEN    !calculate sym%d_wgn only once
         PRINT *, "calculate wigner-matrix"
         call judft_error('WIGNER MATRIX should be available in hybinp part')
      ENDIF

      !$OMP PARALLEL DO default(none) private(na,iatom, iop, ifac, l, i,j) &
      !$OMP shared(atoms, sym, abc, itype, hybinp)
      do na = 1, atoms%neq(itype)
         iatom=atoms%firstAtom(itype)+na-1
         iop = sym%ngopr(iatom)
         !                                    l                        l    l
         ! inversion of spherical harmonics: Y (pi-theta,pi+phi) = (-1)  * Y (theta,phi)
         !                                    m                             m
         ifac = 1
         IF (sym%invsat(iatom) == 2) THEN
            iop = sym%ngopr(sym%invsatnr(iatom))
            ifac = -1
         ENDIF
         DO l = 1, atoms%lmax(itype)
            DO j=1,abc%n_r(l)
               DO i = 1, size(abc%cof,1)
                  abc%cof(i, l**2:l*(l + 2), j,na) = ifac**l*matmul(conjg(hybinp%d_wgn2(-l:l, -l:l, l, iop)), abc%cof(i, l**2:l*(l + 2), j,na))
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !$OMP end parallel do
      call timestop("hyb_abcrot")
   END SUBROUTINE rot_to_unrotated

   subroutine fill_work_array(zmat, noco, atoms, lapw, ccchi, iintsp, nvmax, jspin, ne, work_c)
      use m_types_mat
      use m_types_noco
      use m_types_atoms
      use m_types_lapw
      type(t_mat), intent(in)::zMat
      type(t_noco), intent(in)::noco
      type(t_atoms), intent(in)::atoms
      type(t_lapw),intent(in)::lapw
      complex, intent(in):: ccchi(2, 2)
      integer, intent(in):: iintsp, nvmax, jspin, ne
      complex,INTENT(OUT)            :: work_c(:, :)

      integer:: kspin, i

      CALL timestart("fill work array")
      IF (noco%l_noco) THEN
      IF (noco%l_ss) THEN
!$acc kernels copyin(zmat,zMat%data_c,ccchi,atoms,lapw,lapw%nv) present(work_c)default(none)
! the coefficients of the spin-down basis functions are
! stored in the second half of the eigenvector
         kspin = (iintsp - 1)*(lapw%nv(1) + atoms%nlotot)
         work_c(:nvmax, :) = ccchi(iintsp, jspin)*zMat%data_c(kspin + 1:kspin + nvmax, :ne)
!$acc end kernels
      ELSE
! perform sum over the two interstitial spin directions
! and take into account the spin boundary conditions
! (jspin counts the local spin directions inside each MT)
!$acc kernels copyin(atoms,zMat,zMat%data_c,ccchi,lapw) present(work_c) default(none)
         kspin = lapw%nv(1) + atoms%nlotot
         work_c(:nvmax, :) = ccchi(1, jspin)*zMat%data_c(:nvmax, :ne) + &
                              ccchi(2, jspin)*zMat%data_c(kspin + 1:kspin + nvmax, :ne)
!$acc end kernels
      END IF
      ELSE
      IF (zmat%l_real) THEN
!$CPP_OMP PARALLEL DO default(shared) private(i)
!$acc kernels copyin(zmat,zMat%data_r)present(work_c)default(none)
         DO i = 1, ne
#ifdef _OPENACC
            work_c(:nvmax, i) = zmat%data_r(:nvmax, i)
#else
            work_c(:nvmax, i) = 0.0
            CALL dcopy(nvmax, zMat%data_r(:, i), 1, work_c(:, i), 2)
#endif
         END DO
!$acc end kernels
!$CPP_OMP END PARALLEL DO
      ELSE
!$CPP_OMP PARALLEL DO default(shared) private(i)
!$acc kernels copyin(zMat,zMat%data_c)present(work_c) default(none)
         DO i = 1, ne
#ifdef _OPENACC
            work_c(:nvmax, i) = zmat%data_c(:nvmax, i)
#else
            CALL zcopy(nvmax, zMat%data_c(:, i), 1, work_c(:, i), 1)
#endif
         END DO
!$acc end kernels
!$CPP_OMP END PARALLEL DO
      END IF
      END IF

      CALL timestop("fill work array")
   end subroutine fill_work_array

END MODULE m_types_abc
