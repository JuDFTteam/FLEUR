!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
#ifdef _OPENACC
#define CPP_OMP not_used
#else
#define CPP_OMP $OMP
#endif
MODULE m_hlomat
   use m_matmul_dgemm
  IMPLICIT NONE
  !***********************************************************************
  ! updates the hamiltonian  matrix with the contributions from the local
  ! orbitals.
  ! p.kurz sept. 1996
  !***********************************************************************
CONTAINS
  SUBROUTINE hlomat(input,atoms,fmpi,lapw,tlmplm,sym,cell,noco,nococonv,ilSpinPr,ilSpin,&
       ntyp,na,fjgj,alo1,blo1,clo1, igSpinPr,igSpin,chi,hmat,l_fullj,l_ham,lapwq,fjgjq)

    USE m_hsmt_ab
    USE m_abcoeff_store
    USE m_types
!    USE m_types_mpimat
    USE m_hsmt_fjgj
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)  :: input
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_lapw),INTENT(IN),TARGET   :: lapw
    TYPE(t_mpi),INTENT(IN)    :: fmpi
    TYPE(t_tlmplm),INTENT(IN) :: tlmplm
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_noco),INTENT(IN)   :: noco
    TYPE(t_nococonv),INTENT(IN),TARGET   :: nococonv
    TYPE(t_fjgj),INTENT(IN),TARGET   :: fjgj


    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: na,ntyp
    INTEGER, INTENT (IN) :: ilSpinPr,ilSpin !spins of the local Hamiltonian
    INTEGER, INTENT (IN) :: igSpin,igSpinPr
    COMPLEX, INTENT (IN) :: chi
    !     ..
    !     .. Array Arguments ..
    REAL, INTENT (IN) :: alo1(:,:),blo1(:,:),clo1(:,:)

    CLASS(t_mat),INTENT (INOUT) :: hmat
    LOGICAL, INTENT(IN) :: l_fullj, l_ham

    TYPE(t_lapw), OPTIONAL, INTENT(IN),TARGET :: lapwq
    TYPE(t_fjgj), OPTIONAL, INTENT(IN),TARGET :: fjgjq
    !     ..
    ! Local Scalars
      INTEGER :: invsfct,i,k,locol,nb,ab_size,ab_size_Pr
      LOGICAL :: l_samelapw

      ! Local Arrays
      COMPLEX, ALLOCATABLE :: abCoeffs(:,:), abCoeffsPr(:,:)
      COMPLEX, ALLOCATABLE :: abclo(:,:,:,:), abcloPr(:,:,:,:)
      COMPLEX, ALLOCATABLE :: hmt(:,:), c(:,:), cPr(:,:), x(:,:), y(:,:), q(:,:), w(:,:), z(:,:)
      INTEGER, ALLOCATABLE :: glob(:), globPr(:)


      TYPE(t_lapw) ,POINTER:: lapwPr
      TYPE(t_fjgj) ,POINTER:: fjgjPr

      l_samelapw = .FALSE.
      IF (.NOT.PRESENT(lapwq)) l_samelapw = .TRUE.
      IF (.NOT.l_samelapw) THEN
         lapwPr => lapwq
         fjgjPr => fjgjq
      ELSE
         lapwPr => lapw
         fjgjPr => fjgj
      END IF

      ALLOCATE(abclo(3,-atoms%llod:atoms%llod,2*(2*atoms%llod+1),atoms%nlod))
      ALLOCATE(abcloPr(3,-atoms%llod:atoms%llod,2*(2*atoms%llod+1),atoms%nlod))
      !$acc data create(abclo,abcloPr)
      !$acc data copyin(alo1,blo1,clo1,fjgjPr,fjgjpr%fj,fjgjpr%gj)
      call timestart("hsmt_ab")
      ! Own the abCoeffs mapping in the caller's scope -- see types_abc.F90 for why
      ! hsmt_ab must not do the `enter data` on its own dummy argument.
      IF (.NOT.l_use_abcoeff_store) THEN
         ab_size_Pr = hsmt_ab_size(atoms, ntyp, .TRUE.)
         IF (ALLOCATED(abCoeffsPr)) THEN
            IF (SIZE(abCoeffsPr,1)/=2*ab_size_Pr .OR. SIZE(abCoeffsPr,2)/=lapwPr%nv(igSpinPr)) THEN
               !$acc exit data delete(abCoeffsPr)
               DEALLOCATE(abCoeffsPr)
            END IF
         END IF
         IF (.NOT.ALLOCATED(abCoeffsPr)) THEN
            ALLOCATE(abCoeffsPr(2*ab_size_Pr, lapwPr%nv(igSpinPr)))
            !$acc enter data create(abCoeffsPr)
         END IF
      END IF
      CALL hsmt_ab(sym,atoms,noco,nococonv,ilSpinPr,igSpinPr,ntyp,na,cell,lapwPr,fjgjPr,abCoeffsPr,ab_size_Pr,.TRUE.,abcloPr,alo1(:,ilSpinPr),blo1(:,ilSpinPr),clo1(:,ilSpinPr))
      call timestop("hsmt_ab")
      !we need the "unprimed" abcoeffs
      IF (ilSpin==ilSpinPr.AND.igSpinPr==igSpin.AND.l_samelapw) THEN
         call timestart("abcoeffs copy")
         ! abcoeffs is not produced by hsmt_ab in this branch; allocate it
         ! (matching abcoeffsPr) and put it on the device so it can be copied.
         ALLOCATE(abcoeffs, mold=abcoeffsPr)
         !$acc enter data create(abcoeffs)
         !$acc kernels present(abcoeffs,abcoeffsPr,abclo,abcloPr)
         if (l_fullj) abcoeffs=abcoeffsPr
         abclo=abcloPr
         !$acc end kernels
         call timestop("abcoeffs copy")
      ELSE
         call timestart("hsmt_ab2")
         ! Own the abCoeffs mapping in the caller's scope -- see types_abc.F90 for why
         ! hsmt_ab must not do the `enter data` on its own dummy argument.
         IF (.NOT.l_use_abcoeff_store) THEN
            ab_size = hsmt_ab_size(atoms, ntyp, .TRUE.)
            IF (ALLOCATED(abCoeffs)) THEN
               IF (SIZE(abCoeffs,1)/=2*ab_size .OR. SIZE(abCoeffs,2)/=lapw%nv(igSpin)) THEN
                  !$acc exit data delete(abCoeffs)
                  DEALLOCATE(abCoeffs)
               END IF
            END IF
            IF (.NOT.ALLOCATED(abCoeffs)) THEN
               ALLOCATE(abCoeffs(2*ab_size, lapw%nv(igSpin)))
               !$acc enter data create(abCoeffs)
            END IF
         END IF
         CALL hsmt_ab(sym,atoms,noco,nococonv,ilSpin,igSpin,ntyp,na,cell,lapw,fjgj,abCoeffs,ab_size,.TRUE.,abclo,alo1(:,ilSpin),blo1(:,ilSpin),clo1(:,ilSpin),l_store=.TRUE.)
         call timestop("hsmt_ab2")
      END IF
   

      IF ((sym%invsat(na) == 0) .OR. (sym%invsat(na) == 1)) THEN
         ! If this atom is the first of two atoms related by inversion, the
         ! contributions of both atoms are added simultaneously, using that their
         ! sum is twice the real part of the contribution of each atom. In this
         ! case there are twice as many (2*(2*l+1)) k-vectors (see abccoflo).
         invsfct = MERGE(2,1,sym%invsat(na)==1)
         nb = tlmplm%nbas(ntyp)
         ! LO basis functions of this atom as columns in the unified radial basis
         !$acc update self(abclo,abcloPr)
         CALL lo_coefficients(atoms,tlmplm,ntyp,na,invsfct,abclo,lapw,igSpin,c,glob)
         CALL lo_coefficients(atoms,tlmplm,ntyp,na,invsfct,abcloPr,lapwPr,igSpinPr,cPr,globPr)
         hmt = tlmplm%h(:nb-1,:nb-1,ntyp,ilSpinPr,ilSpin)
         ALLOCATE(x(nb,SIZE(c,2)),y(lapwPr%nv(igSpinPr),SIZE(c,2)),q(SIZE(cPr,2),SIZE(c,2)))
         ALLOCATE(w(SIZE(cPr,2),ab_size_Pr),z(SIZE(cPr,2),lapw%nv(igSpin)))
         !$acc data copyin(hmt,c,cPr,glob,globPr) create(x,y,q,w,z)

         call blas_matmul(nb,SIZE(c,2),nb,hmt,c,x)
         CALL timestart("LAPW-LO")
         call blas_matmul(lapwPr%nv(igSpinPr),SIZE(c,2),ab_size_Pr,abCoeffsPr,x,y,cmplx(1.0,0.0),cmplx(0.0,0.0),'C')
         !$acc kernels present(hmat,hmat%data_c,hmat%data_r,y,glob) copyin(fmpi,fmpi%n_size,fmpi%n_rank,lapwPr,lapwPr%nv)
         DO i = 1,SIZE(c,2)
            IF (MOD(glob(i)-1,fmpi%n_size) /= fmpi%n_rank) CYCLE
            locol = (glob(i)-1)/fmpi%n_size+1
            IF (hmat%l_real) THEN
               hmat%data_r(:lapwPr%nv(igSpinPr),locol) = hmat%data_r(:lapwPr%nv(igSpinPr),locol) + REAL(chi*invsfct*y(:,i))
            ELSE
               hmat%data_c(:lapwPr%nv(igSpinPr),locol) = hmat%data_c(:lapwPr%nv(igSpinPr),locol) + chi*invsfct*y(:,i)
            END IF
         END DO
         !$acc end kernels
         CALL timestop("LAPW-LO")

         IF (l_fullj) THEN
            CALL timestart("LO-LAPW")
            call blas_matmul(SIZE(cPr,2),ab_size_Pr,nb,cPr,hmt,w,cmplx(1.0,0.0),cmplx(0.0,0.0),'C')
            call blas_matmul(SIZE(cPr,2),lapw%nv(igSpin),ab_size_Pr,w,abCoeffs,z)
            !$acc kernels present(hmat,hmat%data_c,hmat%data_r,z,globPr) copyin(fmpi,fmpi%n_size,fmpi%n_rank,lapw,lapw%nv)
            DO k = fmpi%n_rank+1,lapw%nv(igSpin),fmpi%n_size
               locol = (k-1)/fmpi%n_size+1
               IF (hmat%l_real) THEN
                  hmat%data_r(globPr,locol) = hmat%data_r(globPr,locol) + REAL(chi*invsfct*z(:,k))
               ELSE
                  hmat%data_c(globPr,locol) = hmat%data_c(globPr,locol) + chi*invsfct*z(:,k)
               END IF
            END DO
            !$acc end kernels
            CALL timestop("LO-LAPW")
         END IF

         CALL timestart("LO-LO")
         call blas_matmul(SIZE(cPr,2),SIZE(c,2),nb,cPr,x,q,cmplx(1.0,0.0),cmplx(0.0,0.0),'C')
         !$acc kernels present(hmat,hmat%data_c,hmat%data_r,q,glob,globPr) copyin(fmpi,fmpi%n_size,fmpi%n_rank)
         DO i = 1,SIZE(c,2)
            IF (MOD(glob(i)-1,fmpi%n_size) /= fmpi%n_rank) CYCLE
            locol = (glob(i)-1)/fmpi%n_size+1
            DO k = 1,SIZE(cPr,2)
               ! upper triangle only, unless the full matrix is requested
               IF (globPr(k)>glob(i).AND..NOT.l_fullj) CYCLE
               IF (hmat%l_real) THEN
                  hmat%data_r(globPr(k),locol) = hmat%data_r(globPr(k),locol) + REAL(chi*invsfct*q(k,i))
               ELSE
                  hmat%data_c(globPr(k),locol) = hmat%data_c(globPr(k),locol) + chi*invsfct*q(k,i)
               END IF
            END DO
         END DO
         !$acc end kernels
         CALL timestop("LO-LO")
         !$acc end data
      END IF

      ! abCoeffs/abCoeffsPr are allocated and put on the device inside hsmt_ab
      ! (or, for abcoeffs in the copy branch, just above); release them here.
      !$acc exit data delete(abcoeffsPr)
      DEALLOCATE(abcoeffsPr)
      IF (ALLOCATED(abcoeffs)) THEN
         !$acc exit data delete(abcoeffs)
         ! Store only the genuine unprimed abCoeffs produced by hsmt_ab (the ELSE
         ! branch above). The copy branch holds primed data, which must not be
         ! cached under the unprimed (nk,igSpin,ilSpin,na) key.
         IF (.NOT.(ilSpin==ilSpinPr.AND.igSpinPr==igSpin.AND.l_samelapw)) &
              CALL abcoeff_store_save(abcoeffs, lapw%nk, igSpin, ilSpin, na, .TRUE.)
         IF (ALLOCATED(abcoeffs)) DEALLOCATE(abcoeffs)
      END IF

      !$acc end data
      !$acc end data
   END SUBROUTINE hlomat

   SUBROUTINE lo_coefficients(atoms,tlmplm,ntyp,na,invsfct,abclo,lapw,igSpin,c,glob)
      !! LO basis functions of atom na as columns in the unified radial basis of tlmplm%h:
      !! a*u + b*udot + c*u_lo for each m, and their global index in the matrix
      USE m_types
      TYPE(t_atoms),  INTENT(IN) :: atoms
      TYPE(t_tlmplm), INTENT(IN) :: tlmplm
      TYPE(t_lapw),   INTENT(IN) :: lapw
      INTEGER,        INTENT(IN) :: ntyp,na,invsfct,igSpin
      COMPLEX,        INTENT(IN) :: abclo(:,-atoms%llod:,:,:)
      COMPLEX, ALLOCATABLE, INTENT(OUT) :: c(:,:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: glob(:)

      INTEGER :: lo,l,m,lm,nkvec,col

      ALLOCATE(c(tlmplm%nbas(ntyp),SUM(invsfct*(2*atoms%llo(:atoms%nlo(ntyp),ntyp)+1))),source=CMPLX(0.0,0.0))
      ALLOCATE(glob(SIZE(c,2)))
      col = 0
      DO lo = 1,atoms%nlo(ntyp)
         l = atoms%llo(lo,ntyp)
         DO nkvec = 1,invsfct*(2*l+1)
            col = col+1
            glob(col) = lapw%nv(igSpin)+lapw%index_lo(lo,na)+nkvec
            DO m = -l,l
               lm = l*(l+1)+m
               c(tlmplm%ind(1,lm,ntyp)+1,col) = abclo(1,m,nkvec,lo)
               c(tlmplm%ind(2,lm,ntyp)+1,col) = abclo(2,m,nkvec,lo)
               c(tlmplm%ind(atoms%slot_of_lo(lo,ntyp),lm,ntyp)+1,col) = abclo(3,m,nkvec,lo)
            END DO
         END DO
      END DO
   END SUBROUTINE lo_coefficients
END MODULE m_hlomat
