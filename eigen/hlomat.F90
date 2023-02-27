!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
#ifdef _OPENACC
#define CPP_OMP not_used
#else
#define CPP_OMP $OMP
#endif
MODULE m_hlomat
  IMPLICIT NONE
  !***********************************************************************
  ! updates the hamiltonian  matrix with the contributions from the local
  ! orbitals.
  ! p.kurz sept. 1996
  !***********************************************************************
CONTAINS
  SUBROUTINE hlomat(input,atoms,fmpi,lapw,ud,tlmplm,sym,cell,noco,nococonv,ilSpinPr,ilSpin,&
       ntyp,na,fjgj,alo1,blo1,clo1, igSpinPr,igSpin,chi,hmat,l_fullj,l_ham,lapwq,fjgjq)

    USE m_hsmt_ab
    USE m_types
!    USE m_types_mpimat
    USE m_hsmt_fjgj
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)  :: input
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_lapw),INTENT(IN),TARGET   :: lapw
    TYPE(t_mpi),INTENT(IN)    :: fmpi
    TYPE(t_usdus),INTENT(IN)  :: ud
    TYPE(t_tlmplm),INTENT(IN) :: tlmplm
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_noco),INTENT(IN)   :: noco
    TYPE(t_nococonv),INTENT(IN),TARGET   :: nococonv
    TYPE(t_fjgj),INTENT(IN),TARGET   :: fjgj


    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: na,ntyp
    INTEGER, INTENT (IN) :: ilSpinPr,ilSpin !spin for usdus and tlmplm
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
      COMPLEX :: axx,bxx,cxx,dtd,dtu,tdulo,tulod,tulou,tuloulo,utd,utu, tuulo
      INTEGER :: invsfct,l,lm,lmp,lo,lolo,lolop,lop,lp,i
      INTEGER :: mp,nkvec,nkvecp,lmplm,loplo,kp,m,mlo,mlolo,mlolo_new,lolop_new
      INTEGER :: locol,lorow,n,k,ab_size,ab_size_Pr,s
      LOGICAL :: l_samelapw

      ! Local Arrays
      COMPLEX, ALLOCATABLE :: abCoeffs(:,:), ax(:), bx(:), cx(:)
      COMPLEX, ALLOCATABLE :: abclo(:,:,:,:)
      COMPLEX, ALLOCATABLE :: abCoeffsPr(:,:), axPr(:), bxPr(:), cxPr(:)
      COMPLEX, ALLOCATABLE :: abcloPr(:,:,:,:)

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

      ! Synthesize a and b
      
      ALLOCATE(abclo(3,-atoms%llod:atoms%llod,2*(2*atoms%llod+1),atoms%nlod))
      ALLOCATE(ax(MAXVAL(lapw%nv)),bx(MAXVAL(lapw%nv)),cx(MAXVAL(lapw%nv)))      
      ALLOCATE(abCoeffsPr(0:2*atoms%lnonsph(ntyp)*(atoms%lnonsph(ntyp)+2)+1,MAXVAL(lapwPr%nv)))
      ALLOCATE(axPr(MAXVAL(lapwPr%nv)),bxPr(MAXVAL(lapwPr%nv)),cxPr(MAXVAL(lapwPr%nv)))
      ALLOCATE(abcloPr(3,-atoms%llod:atoms%llod,2*(2*atoms%llod+1),atoms%nlod))

      !$acc data create(abcoeffs,abclo,abcoeffsPr,abcloPr)
      !$acc data copyin(alo1,blo1,clo1,fjgjPr,fjgjpr%fj,fjgjpr%gj)
      CALL hsmt_ab(sym,atoms,noco,nococonv,ilSpinPr,igSpinPr,ntyp,na,cell,lapwPr,fjgjPr,abCoeffsPr(:,:),ab_size_Pr,.TRUE.,abcloPr,alo1(:,ilSpinPr),blo1(:,ilSpinPr),clo1(:,ilSpinPr))

      !we need the "unprimed" abcoeffs
      IF (ilSpin==ilSpinPr.AND.igSpinPr==igSpin.AND.l_samelapw) THEN
         !$acc kernels present(abcoeffs,abcoeffsPr,abclo,abcloPr)
         if (l_fullj) abcoeffs=abcoeffsPr !TODO automatic alloc on GPU????
         abclo(:,:,:,:)=abcloPr(:,:,:,:) 
         !$acc end kernels
      ELSE     
         ALLOCATE(abCoeffs(0:2*atoms%lnonsph(ntyp)*(atoms%lnonsph(ntyp)+2)+1,MAXVAL(lapw%nv)))       
         CALL hsmt_ab(sym,atoms,noco,nococonv,ilSpin,igSpin,ntyp,na,cell,lapw,fjgj,abCoeffs(:,:),ab_size,.TRUE.,abclo,alo1(:,ilSpin),blo1(:,ilSpin),clo1(:,ilSpin))
      END IF
   

      mlo=0;mlolo=0;mlolo_new=0
      DO m=1,ntyp-1
         mlo=mlo+atoms%nlo(m)
         mlolo=mlolo+atoms%nlo(m)*(atoms%nlo(m)+1)/2
         mlolo_new=mlolo_new+atoms%nlo(m)**2
      END DO

      IF ((sym%invsat(na) == 0) .OR. (sym%invsat(na) == 1)) THEN
         ! If this atom is the first of two atoms related by inversion, the
         ! contributions to the overlap matrix of both atoms are added simultaneously.
         ! Where it is made use of the fact, that the sum of these contributions
         ! is twice the real part of the contribution of each atom. Note, that
         ! in this case there are twice as many (2*(2*l+1)) k-vectors.
         ! (compare abccoflo and comments there).
         IF (sym%invsat(na) == 0) invsfct = 1
         IF (sym%invsat(na) == 1) invsfct = 2

         !$acc kernels present(hmat,hmat%data_c,hmat%data_r,abcoeffs,abclo,abcoeffsPr,abcloPr) &
         !$acc & copyin(atoms,lapw,lapwPr,tlmplm,tlmplm%tulou,tlmplm%tulod,tlmplm%h_loc(:,:,ntyp,ilSpinPr,ilSpin),lapw%nv(:),lapwPr%nv(:))&
         !$acc & copyin(tlmplm%tdulo(:,:,:,ilSpinPr,ilSpin),tlmplm%tuloulo_newer(:,:,:,:,ntyp,ilSpinPr,ilSpin),atoms%rmt(ntyp))&
         !$acc & copyin(lapw%index_lo(:,na),lapwPr%index_lo(:,na),tlmplm%h_loc2,atoms%llo(:,ntyp),atoms%nlo(ntyp))&
         !$acc & copyin(atoms%lnonsph(ntyp))&
         !$acc & copyin(ud,ud%us(:,ntyp,ilSpin),ud%uds(:,ntyp,ilSpin),ud%dus(:,ntyp,ilSpin),ud%dulos(:,ntyp,ilSpin),ud%duds(:,ntyp,ilSpin))&
         !$acc & copyin(input, input%l_useapw, fmpi, fmpi%n_size, fmpi%n_rank)&
         !$acc & create(ax,bx,cx,axpr,bxpr,cxpr)&
         !$acc & default(none)
         DO lo = 1,atoms%nlo(ntyp)
            l = atoms%llo(lo,ntyp)
            ! Calculate the hamiltonian matrix elements with the regular
            ! Flapw basis-functions
            DO m = -l,l
               lm = l* (l+1) + m
               s = tlmplm%h_loc2_nonsph(ntyp) 
               print*, s,ab_size_pr/2
               axPr = matmul(transpose(conjg(abCoeffsPr(0:2*s-1,:))),tlmplm%h_loc_LO(0:2*s-1,lm,ntyp,ilSpinPr,ilSpin))
               bxPr = matmul(transpose(conjg(abCoeffsPr(0:2*s-1,:))),tlmplm%h_loc_LO(0:2*s-1,s+lm,ntyp,ilSpinPr,ilSpin))
               cxPr = matmul(transpose(conjg(abCoeffsPr(0:2*s-1,:))),tlmplm%h_LO(0:2*s-1,m,lo+mlo,ilSpinPr,ilSpin))
              
               !axpr=0.0
               !bxpr=0.0
               !cxpr=0.0
               !DO kp = 1, lapwPr%nv(igSpinPr)
               !   DO lp = 0, atoms%lnonsph(ntyp)
               !      DO mp = -lp, lp
               !         lmp = lp*(lp+1) + mp
                        !axPr(kp) = axPr(kp) + CONJG(abCoeffsPr(lmp,kp))             *tlmplm%h_loc_LO(lmp,lm,ntyp,ilSpinPr,ilSpin)
                        !axPr(kp) = axPr(kp) + CONJG(abCoeffsPr(ab_size_Pr/2+lmp,kp))*tlmplm%h_loc_LO(s+lmp,lm,ntyp,ilSpinPr,ilSpin)
                        !bxPr(kp) = bxPr(kp) + CONJG(abCoeffsPr(lmp,kp))             *tlmplm%h_loc_LO(lmp,s+lm,ntyp,ilSpinPr,ilSpin)
                        !bxPr(kp) = bxPr(kp) + CONJG(abCoeffsPr(ab_size_Pr/2+lmp,kp))*tlmplm%h_loc_LO(s+lmp,s+lm,ntyp,ilSpinPr,ilSpin)
                        !cxPr(kp) = cxPr(kp) + CONJG(abCoeffsPr(lmp,kp))             *tlmplm%tuulo(lmp,m,lo+mlo,ilSpinPr,ilSpin)
                        !cxPr(kp) = cxPr(kp) + CONJG(abCoeffsPr(ab_size_Pr/2+lmp,kp))*tlmplm%tdulo(lmp,m,lo+mlo,ilSpinPr,ilSpin)
               !      END DO
               !   END DO
               !END DO
               DO nkvec = 1,invsfct*(2*l+1)
                  locol= lapw%nv(igSpin)+lapw%index_lo(lo,na)+nkvec ! This is the column of the matrix
                  IF (MOD(locol-1,fmpi%n_size) == fmpi%n_rank) THEN ! Only this MPI rank calculates this column
                     locol=(locol-1)/fmpi%n_size+1 ! This is the column in local storage
                     IF (hmat%l_real) THEN
                        DO kp = 1,lapwPr%nv(igSpinPr)
                           hmat%data_r(kp,locol) = hmat%data_r(kp,locol) &
                                               & + REAL(chi) * invsfct * (&
                                               & abclo(1,m,nkvec,lo) *  axPr(kp) + &
                                               & abclo(2,m,nkvec,lo) *  bxPr(kp) + &
                                               & abclo(3,m,nkvec,lo) *  cxPr(kp) )
                           IF (input%l_useapw) THEN
                              ! APWlo
                              hmat%data_r(kp,locol) = hmat%data_r(kp,locol) &
                                                  & + 0.25 * atoms%rmt(ntyp)**2 &
                                                  & * REAL(chi) * invsfct * ( &
                                                  & (CONJG(abCoeffsPr(lm,kp))           * ud%us(l,ntyp,ilSpin)   + &
                                                  &  CONJG(abCoeffsPr(ab_size_Pr/2+lm,kp)) * ud%uds(l,ntyp,ilSpin)) * &
                                                  & ( abclo(1,m,nkvec,lo) * ud%dus(l,ntyp,ilSpin)  &
                                                  & + abclo(2,m,nkvec,lo) * ud%duds(l,ntyp,ilSpin) &
                                                  & + abclo(3,m,nkvec,lo) * ud%dulos(lo,ntyp,ilSpin)) )
                           END IF
                       END DO
                     ELSE
                        DO kp = 1,lapwPr%nv(igSpinPr)
                           hmat%data_c(kp,locol) = hmat%data_c(kp,locol) &
                                               & + chi * invsfct * ( &
                                               & abclo(1,m,nkvec,lo) *  axPr(kp) + &
                                               & abclo(2,m,nkvec,lo) *  bxPr(kp) + &
                                               & abclo(3,m,nkvec,lo) *  cxPr(kp) )
                           IF (input%l_useapw) THEN
                              ! APWlo
                              hmat%data_c(kp,locol) = hmat%data_c(kp,locol) &
                                                  & + 0.25 * atoms%rmt(ntyp)**2 &
                                                  & * chi * invsfct * ( &
                                                  & (CONJG(abCoeffsPr(lm,kp))           * ud%us(l,ntyp,ilSpin)   + &
                                                  &  CONJG(abCoeffsPr(ab_size_Pr/2+lm,kp)) * ud%uds(l,ntyp,ilSpin)) * &
                                                  & ( abclo(1,m,nkvec,lo) * ud%dus(l,ntyp,ilSpin)  &
                                                  & + abclo(2,m,nkvec,lo) * ud%duds(l,ntyp,ilSpin) &
                                                  & + abclo(3,m,nkvec,lo) * ud%dulos(lo,ntyp,ilSpin)) )
                           END IF
                           IF (l_ham.AND.l_fullj.AND.ilSpinPr.EQ.ilSpin) THEN
                              ! TODO: Is it 0.25 or 0.5?
                              hmat%data_c(kp,locol) = hmat%data_c(kp,locol) &
                                                  !& + 0.5 * atoms%rmt(ntyp)**2 &
                                                  & + 0.25 * atoms%rmt(ntyp)**2 &
                                                  & * chi * invsfct * &
                                                  & (CONJG(abCoeffsPr(lm,kp))              * ud%us(l,ntyp,ilSpin)   + &
                                                  &  CONJG(abCoeffsPr(ab_size_Pr/2+lm,kp)) * ud%uds(l,ntyp,ilSpin)) * &
                                                  & ( abclo(1,m,nkvec,lo) * ud%dus(l,ntyp,ilSpin)  &
                                                  & + abclo(2,m,nkvec,lo) * ud%duds(l,ntyp,ilSpin) &
                                                  & + abclo(3,m,nkvec,lo) * ud%dulos(lo,ntyp,ilSpin))
                           END IF
                        END DO
                     END IF
                     ! Jump to the last matrix element of the current row
                  END IF
               END DO
               IF (l_fullj) THEN
                  DO k = 1, lapw%nv(igSpin)
                     ax(k) = CMPLX(0.0,0.0)
                     bx(k) = CMPLX(0.0,0.0)
                     cx(k) = CMPLX(0.0,0.0)
                  END DO
                  !CPP_OMP PARALLEL DO DEFAULT(none) &
                  !CPP_OMP& SHARED(ax,bx,cx,ntyp,ilSpin,ilSpinPr,m,lm,lo,mlo) &
                  !CPP_OMP& SHARED(lapw,abCoeffs,ab_size,igSpin) &
                  !CPP_OMP& SHARED(atoms,tlmplm) &
                  !CPP_OMP  PRIVATE(lp,mp,lmp,s)
                  DO k = 1, lapw%nv(igSpin)
                     DO lp = 0, atoms%lnonsph(ntyp)
                        DO mp = -lp, lp
                           lmp = lp*(lp+1) + mp
                           s = tlmplm%h_loc2_nonsph(ntyp)
                           ax(k) = ax(k) + tlmplm%h_loc_LO(lm,lmp,ntyp,ilSpinPr,ilSpin)     * abCoeffs(lmp,k)
                           ax(k) = ax(k) + tlmplm%h_loc_LO(lm,s+lmp,ntyp,ilSpinPr,ilSpin)   * abCoeffs(ab_size/2+lmp,k)
                           bx(k) = bx(k) + tlmplm%h_loc_LO(s+lm,lmp,ntyp,ilSpinPr,ilSpin)   * abCoeffs(lmp,k)
                           bx(k) = bx(k) + tlmplm%h_loc_LO(s+lm,s+lmp,ntyp,ilSpinPr,ilSpin) * abCoeffs(ab_size/2+lmp,k)
                           cx(k) = cx(k) + tlmplm%tulou(lmp,m,lo+mlo,ilSpinPr,ilSpin)    * abCoeffs(lmp,k)
                           cx(k) = cx(k) + tlmplm%tulod(lmp,m,lo+mlo,ilSpinPr,ilSpin)    * abCoeffs(ab_size/2+lmp,k)
                        END DO
                     END DO
                  END DO
                  !CPP_OMP END PARALLEL DO
                  DO nkvec = 1,invsfct*(2*l+1)
                     lorow = lapwPr%nv(igSpinPr)+lapwPr%index_lo(lo,na)+nkvec
                     !IF (MOD(lorow-1,fmpi%n_size) == fmpi%n_rank) THEN
                        !lorow=(lorow-1)/fmpi%n_size+1
                        DO k = 1,lapw%nv(igSpin)
                           hmat%data_c(lorow,k) = hmat%data_c(lorow,k) &
                                              & + chi * invsfct * ( &
                                              & CONJG(abcloPr(1,m,nkvec,lo)) * ax(k) + &
                                              & CONJG(abcloPr(2,m,nkvec,lo)) * bx(k) + &
                                              & CONJG(abcloPr(3,m,nkvec,lo)) * cx(k) )
                           IF (l_ham.AND.ilSpinPr.EQ.ilSpin) THEN
                              ! TODO: Is it 0.25 or 0.5?
                              hmat%data_c(lorow,k) = hmat%data_c(lorow,k) &
                                                 & + 0.25 * atoms%rmt(ntyp)**2 &
                                                 & * chi * invsfct * &
                                                 & ( CONJG(abcloPr(1,m,nkvec,lo)) * ud%dus(l,ntyp,ilSpin)       &
                                                 & + CONJG(abcloPr(2,m,nkvec,lo)) * ud%duds(l,ntyp,ilSpin)      &
                                                 & + CONJG(abcloPr(3,m,nkvec,lo)) * ud%dulos(lo,ntyp,ilSpin)) * &
                                                 & (abCoeffs(lm,k)                * ud%us(l,ntyp,ilSpin)   + &
                                                 &  abCoeffs(ab_size/2+lm,k)      * ud%uds(l,ntyp,ilSpin))

                           END IF
                        END DO
                     !END IF
                  END DO
               END IF
            END DO
            ! Calculate the hamiltonian matrix elements with other local
            ! orbitals at the same atom and with itself
            DO nkvec = 1,invsfct* (2*l+1)
               locol = lapw%nv(igSpin)+lapw%index_lo(lo,na)+nkvec ! This is the column of the matrix
               IF (MOD(locol-1,fmpi%n_size) == fmpi%n_rank) THEN ! Only this MPI rank calculates this column
                  locol=(locol-1)/fmpi%n_size+1 ! This is the column in local storage
                  ! Calculate the Hamiltonian matrix elements with different
                  ! local orbitals at the same atom, if they have the same l
                  DO lop = 1, MERGE(lo-1,atoms%nlo(ntyp),igSpinPr==igSpin.AND..NOT.l_fullj)
                     IF (lop==lo) CYCLE
                     lp = atoms%llo(lop,ntyp)
                     DO nkvecp = 1,invsfct* (2*lp+1)
                        lorow = lapwPr%nv(igSpinPr)+lapwPr%index_lo(lop,na)+nkvecp
                        DO m = -l,l
                           lm = l*(l+1) + m
                           DO mp = -lp,lp
                              lmp = lp* (lp+1) + mp
                              s = tlmplm%h_loc2_nonsph(ntyp)
                              ! Note, that xtx are the t-matrices and NOT their
                              ! respective complex conjugates as in hssphn !TODO: outdated comment?
                              utu = tlmplm%h_loc_LO(lmp,lm,ntyp,ilSpinPr,ilSpin)
                              dtu = tlmplm%h_loc_LO(lmp+s,lm,ntyp,ilSpinPr,ilSpin)
                              utd = tlmplm%h_loc_LO(lmp,lm+s,ntyp,ilSpinPr,ilSpin)
                              dtd = tlmplm%h_loc_LO(lmp+s,lm+s,ntyp,ilSpinPr,ilSpin)

                              tuulo = tlmplm%h_LO(lmp,m,lo+mlo,ilSpinPr,ilSpin)
                              tdulo = tlmplm%h_LO(lmp+s,m,lo+mlo,ilSpinPr,ilSpin)
                              tulou = tlmplm%tulou(lm,mp,lop+mlo,ilSpinPr,ilSpin)
                              tulod = tlmplm%tulod(lm,mp,lop+mlo,ilSpinPr,ilSpin)

                              tuloulo = tlmplm%tuloulo_newer(mp,m,lop,lo,ntyp,ilSpinPr,ilSpin)

                              axx = utu     * abclo(1,m,nkvec,lo) &
                                & + utd     * abclo(2,m,nkvec,lo) &
                                & + tuulo   * abclo(3,m,nkvec,lo)
                              bxx = dtu     * abclo(1,m,nkvec,lo) &
                                & + dtd     * abclo(2,m,nkvec,lo) &
                                & + tdulo   * abclo(3,m,nkvec,lo)
                              cxx = tulou   * abclo(1,m,nkvec,lo) &
                                & + tulod   * abclo(2,m,nkvec,lo) &
                                & + tuloulo * abclo(3,m,nkvec,lo)

                              IF (hmat%l_real) THEN
                                 hmat%data_r(lorow,locol) = hmat%data_r(lorow,locol) &
                                                        & + REAL(chi) * invsfct * ( &
                                                        &  REAL(abcloPr(1,mp,nkvecp,lop))* REAL(axx) + &
                                                        & AIMAG(abcloPr(1,mp,nkvecp,lop))*AIMAG(axx) + &
                                                        &  REAL(abcloPr(2,mp,nkvecp,lop))* REAL(bxx) + &
                                                        & AIMAG(abcloPr(2,mp,nkvecp,lop))*AIMAG(bxx) + &
                                                        &  REAL(abcloPr(3,mp,nkvecp,lop))* REAL(cxx) + &
                                                        & AIMAG(abcloPr(3,mp,nkvecp,lop))*AIMAG(cxx) )
                              ELSE
                                 hmat%data_c(lorow,locol) = hmat%data_c(lorow,locol) &
                                                        & + chi * ( &
                                                        & CONJG(abcloPr(1,mp,nkvecp,lop)) * axx + &
                                                        & CONJG(abcloPr(2,mp,nkvecp,lop)) * bxx + &
                                                        & CONJG(abcloPr(3,mp,nkvecp,lop)) * cxx )
                              END IF
                           END DO
                        END DO
                     END DO
                  END DO
                  ! Calculate the Hamiltonian matrix elements of one local
                  ! orbital with itself
                  lop=lo
                  DO nkvecp = 1,MERGE(nkvec,invsfct* (2*l+1),igSpinPr==igSpin.AND..NOT.l_fullj)
                     lorow=lapwPr%nv(igSpinPr)+lapwPr%index_lo(lop,na)+nkvecp
                     DO m = -l,l
                        lm = l* (l+1) + m
                        DO mp = -l,l
                          lmp = l*(l+1) + mp
                           s = tlmplm%h_loc2_nonsph(ntyp)

                           utu = tlmplm%h_loc_LO(lmp,lm,ntyp,ilSpinPr,ilSpin)
                           dtu = tlmplm%h_loc_LO(lmp+s,lm,ntyp,ilSpinPr,ilSpin)
                           utd = tlmplm%h_loc_LO(lmp,lm+s,ntyp,ilSpinPr,ilSpin)
                           dtd = tlmplm%h_loc_LO(lmp+s,lm+s,ntyp,ilSpinPr,ilSpin)
                           
                           tuulo = tlmplm%h_LO(lmp,m,lo+mlo,ilSpinPr,ilSpin)
                           tdulo = tlmplm%h_LO(s+lmp,m,lo+mlo,ilSpinPr,ilSpin)
                           tulou = tlmplm%tulou(lm,mp,lo+mlo,ilSpinPr,ilSpin)
                           tulod = tlmplm%tulod(lm,mp,lo+mlo,ilSpinPr,ilSpin)

                           tuloulo = tlmplm%tuloulo_newer(mp,m,lo,lo,ntyp,ilSpinPr,ilSpin)

                           axx = utu     * abclo(1,m,nkvec,lo) &
                             & + utd     * abclo(2,m,nkvec,lo) &
                             & + tuulo   * abclo(3,m,nkvec,lo)
                           bxx = dtu     * abclo(1,m,nkvec,lo) &
                             & + dtd     * abclo(2,m,nkvec,lo) &
                             & + tdulo   * abclo(3,m,nkvec,lo)
                           cxx = tulou   * abclo(1,m,nkvec,lo) &
                             & + tulod   * abclo(2,m,nkvec,lo) &
                             & + tuloulo * abclo(3,m,nkvec,lo)

                           IF (hmat%l_real) THEN
                              hmat%data_r(lorow,locol) = hmat%data_r(lorow,locol) &
                                                     & + REAL(chi) * invsfct * ( &
                                                     &  REAL(abcloPr(1,mp,nkvecp,lo))* REAL(axx) + &
                                                     & AIMAG(abcloPr(1,mp,nkvecp,lo))*AIMAG(axx) + &
                                                     &  REAL(abcloPr(2,mp,nkvecp,lo))* REAL(bxx) + &
                                                     & AIMAG(abcloPr(2,mp,nkvecp,lo))*AIMAG(bxx)  + &
                                                     &  REAL(abcloPr(3,mp,nkvecp,lo)) * REAL(cxx) + &
                                                     & AIMAG(abcloPr(3,mp,nkvecp,lo))*AIMAG(cxx) )
                           ELSE
                              hmat%data_c(lorow,locol) = hmat%data_c(lorow,locol) &
                                                     & + chi * invsfct * ( &
                                                     & CONJG(abcloPr(1,mp,nkvecp,lo))*axx + &
                                                     & CONJG(abcloPr(2,mp,nkvecp,lo))*bxx + &
                                                     & CONJG(abcloPr(3,mp,nkvecp,lo))*cxx )
                              !TODO: Should there be an Ekin sufrace term here as well?
                           END IF
                        END DO
                     END DO
                  END DO
               END IF !If this lo to be calculated by fmpi rank
            END DO
         END DO ! end of lo = 1,atoms%nlo loop
         !$acc end kernels
      END IF
      !$acc end data
      !$acc end data
   END SUBROUTINE hlomat
END MODULE m_hlomat
