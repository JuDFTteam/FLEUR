!--------------------------------------------------------------------------------
! Copyright (c) 2020 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_abcof

CONTAINS

  ! The subroutine abcof calculates the A, B, and C coefficients for the
  ! eigenfunctions. Also some force contributions can be calculated.
  SUBROUTINE abcof(input,atoms,sym, cell,lapw,ne,usdus,&
                   noco,nococonv,jspin,oneD, acof,bcof,ccof,zMat,eig,force)

    USE m_juDFT
    USE m_types
    USE m_constants
    USE m_ylm
    USE m_setabc1lo
    USE m_abclocdn
    USE m_hsmt_fjgj
    USE m_hsmt_ab

    IMPLICIT NONE

    TYPE(t_input),INTENT(IN)             :: input
    TYPE(t_usdus),INTENT(IN)             :: usdus
    TYPE(t_lapw),INTENT(IN)              :: lapw
    TYPE(t_oneD),INTENT(IN)              :: oneD
    TYPE(t_noco),INTENT(IN)              :: noco
    TYPE(t_nococonv),INTENT(IN)          :: nococonv
    TYPE(t_sym),INTENT(IN)               :: sym
    TYPE(t_cell),INTENT(IN)              :: cell
    TYPE(t_atoms),INTENT(IN)             :: atoms
    TYPE(t_mat),INTENT(IN)               :: zMat
    TYPE(t_force),OPTIONAL,INTENT(INOUT) :: force

    ! scalar arguments
    INTEGER, INTENT(IN)        :: ne
    INTEGER, INTENT(IN)        :: jspin

    ! array arguments
    COMPLEX, INTENT(OUT)       :: acof(:,0:,:)!(nobd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat)
    COMPLEX, INTENT(OUT)       :: bcof(:,0:,:)!(nobd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat)
    COMPLEX, INTENT(OUT)       :: ccof(-atoms%llod:,:,:,:)!(-llod:llod,nobd,atoms%nlod,atoms%nat)
    REAL, OPTIONAL, INTENT(IN) :: eig(:)!(input%neig)

    ! Local objects
    TYPE(t_fjgj) :: fjgj

    ! Local scalars
    INTEGER :: i,iLAPW,l,ll1,lm,nap,jAtom,lmp,m,nkvec,iAtom,iType
    INTEGER :: inv_f,ie,ilo,kspin,iintsp,nintsp,nvmax,lo,inap,abSize
    REAL    :: tmk, qss(3), s2h, s2h_e(ne)
    COMPLEX :: phase, c_1, c_2
    LOGICAL :: l_force

    ! Local arrays
    REAL    :: fg(3),fgp(3),fgr(3),fk(3),fkp(3),fkr(3)
    REAL    :: alo1(atoms%nlod,input%jspins),blo1(atoms%nlod,input%jspins)
    REAL    :: clo1(atoms%nlod,input%jspins)
    COMPLEX :: ylm((atoms%lmaxd+1)**2)
    COMPLEX :: ccchi(2,2)
    REAL,    ALLOCATABLE :: work_r(:,:), realCoeffs(:,:), imagCoeffs(:,:)
    COMPLEX, ALLOCATABLE :: work_c(:,:)
    COMPLEX, ALLOCATABLE :: abCoeffs(:,:)

    CALL timestart("abcof")

    ! Checks
    IF (zmat%l_real) THEN
       IF (noco%l_soc.AND.sym%invs) CALL judft_error("BUG in abcof, SOC&INVS but real?")
       IF (noco%l_noco) CALL judft_error("BUG in abcof, l_noco but real?")
    ENDIF

    ! Allocations
    CALL fjgj%alloc(MAXVAL(lapw%nv),atoms%lmaxd,jspin,noco)
    ALLOCATE(abCoeffs(MAXVAL(lapw%nv),2*atoms%lmaxd*(atoms%lmaxd+2)+2))

    ! Initializations
    acof(:,:,:)   = CMPLX(0.0,0.0)
    bcof(:,:,:)   = CMPLX(0.0,0.0)
    ccof(:,:,:,:) = CMPLX(0.0,0.0)
    l_force = .FALSE.
    IF(PRESENT(eig).AND.input%l_f) l_force = .TRUE.
    IF(l_force) THEN
       force%acoflo  = CMPLX(0.0,0.0)
       force%bcoflo  = CMPLX(0.0,0.0)
       force%e1cof   = CMPLX(0.0,0.0)
       force%e2cof   = CMPLX(0.0,0.0)
       force%aveccof = CMPLX(0.0,0.0)
       force%bveccof = CMPLX(0.0,0.0)
       force%cveccof = CMPLX(0.0,0.0)
    END IF

    ! loop over atoms
    DO iAtom = 1, atoms%nat
       iType = atoms%itype(iAtom)

       CALL timestart("fjgj coefficients")
       CALL fjgj%calculate(input,atoms,cell,lapw,noco,usdus,iType,jspin)
       CALL timestop("fjgj coefficients")

       CALL setabc1lo(atoms,iType,usdus,jspin,alo1,blo1,clo1)

       IF(noco%l_noco) THEN
          ! generate the spinors (chi)
          ccchi(1,1) =  EXP(ImagUnit*nococonv%alph(iType)/2)*COS(nococonv%beta(iType)/2)
          ccchi(1,2) = -EXP(ImagUnit*nococonv%alph(iType)/2)*SIN(nococonv%beta(iType)/2)
          ccchi(2,1) =  EXP(-ImagUnit*nococonv%alph(iType)/2)*SIN(nococonv%beta(iType)/2)
          ccchi(2,2) =  EXP(-ImagUnit*nococonv%alph(iType)/2)*COS(nococonv%beta(iType)/2)
       END IF

       nintsp = 1
       IF (noco%l_ss) nintsp = 2
       ! loop over the interstitial spin
       DO iintsp = 1,nintsp

          nvmax=lapw%nv(jspin)
          IF (noco%l_ss) nvmax=lapw%nv(iintsp)
          qss = MERGE(-1.0,1.0,iintsp.EQ.1)*nococonv%qss/2.0

          IF ((sym%invsat(iAtom).EQ.0) .OR. (sym%invsat(iAtom).EQ.1)) THEN

             IF (zmat%l_real) THEN
                ALLOCATE ( work_r(ne,nvmax) )
             ELSE
                ALLOCATE ( work_c(ne,nvmax) )
             ENDIF

             ! Filling of work array (modified zMat)
             CALL timestart("fill work array")
             IF (noco%l_noco) THEN
                DO iLAPW = 1,nvmax
                   IF (noco%l_ss) THEN
                      ! the coefficients of the spin-down basis functions are
                      ! stored in the second half of the eigenvector
                      kspin = (iintsp-1)*(lapw%nv(1)+atoms%nlotot)
                      work_c(:ne,iLAPW) = ccchi(iintsp,jspin)*zMat%data_c(kspin+iLAPW,:ne)
                   ELSE
                      ! perform sum over the two interstitial spin directions
                      ! and take into account the spin boundary conditions
                      ! (jspin counts the local spin directions inside each MT)
                      kspin = lapw%nv(1)+atoms%nlotot
                      work_c(:ne,iLAPW) = ccchi(1,jspin)*zMat%data_c(iLAPW,:ne) + ccchi(2,jspin)*zMat%data_c(kspin+iLAPW,:ne)
                   END IF
                END DO
             ELSE
                IF (zmat%l_real) THEN
                   work_r(:ne,:)=TRANSPOSE(zMat%data_r(:nvmax,:ne))
                ELSE
                   work_c(:ne,:)=TRANSPOSE(zMat%data_c(:nvmax,:ne))
                END IF
             END IF
             CALL timestop("fill work array")

             ! Calculation of a, b coefficients for LAPW basis functions
             CALL timestart("hsmt_ab")
             CALL hsmt_ab(sym,atoms,noco,nococonv,jspin,iintsp,iType,iAtom,cell,lapw,fjgj,abCoeffs,abSize,.FALSE.)
             abSize = abSize / 2
             CALL timestop("hsmt_ab")

             ! Obtaining A, B coefficients for eigenfunctions
             CALL timestart("gemm")
             IF (zmat%l_real) THEN
                ! variant with zgemm
!                ALLOCATE ( work_c(ne,nvmax) )
!                work_c(:,:) = work_r(:,:)
!                CALL zgemm("N","N",ne,abSize,nvmax,CMPLX(1.0,0.0),work_c,ne,CONJG(abCoeffs(:nvmax,:abSize)),nvmax,CMPLX(1.0,0.0),acof(:ne,0:abSize-1,iAtom),ne)
!                CALL zgemm("N","N",ne,abSize,nvmax,CMPLX(1.0,0.0),work_c,ne,CONJG(abCoeffs(:nvmax,abSize+1:2*abSize)),nvmax,CMPLX(1.0,0.0),bcof(:ne,0:abSize-1,iAtom),ne)
!                DEALLOCATE(work_c)
                ! variant with dgemm
                ALLOCATE(realCoeffs(ne,0:abSize-1),imagCoeffs(ne,0:abSize-1))
                realCoeffs = 0.0
                imagCoeffs = 0.0
                CALL dgemm("N","N",ne,abSize,nvmax,1.0,work_r,ne,REAL(abCoeffs(:nvmax,:abSize)),nvmax,0.0,realCoeffs,ne)
                CALL dgemm("N","N",ne,abSize,nvmax,-1.0,work_r,ne,AIMAG(abCoeffs(:nvmax,:abSize)),nvmax,0.0,imagCoeffs,ne)
                acof(:ne,0:abSize-1,iAtom) = acof(:ne,0:abSize-1,iAtom) + CMPLX(realCoeffs(:,:),imagCoeffs(:,:))
                realCoeffs = 0.0
                imagCoeffs = 0.0
                CALL dgemm("N","N",ne,abSize,nvmax,1.0,work_r,ne,REAL(abCoeffs(:nvmax,abSize+1:2*abSize)),nvmax,0.0,realCoeffs,ne)
                CALL dgemm("N","N",ne,abSize,nvmax,-1.0,work_r,ne,AIMAG(abCoeffs(:nvmax,abSize+1:2*abSize)),nvmax,0.0,imagCoeffs,ne)
                bcof(:ne,0:abSize-1,iAtom) = bcof(:ne,0:abSize-1,iAtom) + CMPLX(realCoeffs(:,:),imagCoeffs(:,:))
                DEALLOCATE(realCoeffs,imagCoeffs)
             ELSE
                CALL zgemm("N","N",ne,abSize,nvmax,CMPLX(1.0,0.0),work_c,ne,CONJG(abCoeffs(:nvmax,:abSize)),nvmax,CMPLX(1.0,0.0),acof(:ne,0:abSize-1,iAtom),ne)
                CALL zgemm("N","N",ne,abSize,nvmax,CMPLX(1.0,0.0),work_c,ne,CONJG(abCoeffs(:nvmax,abSize+1:2*abSize)),nvmax,CMPLX(1.0,0.0),bcof(:ne,0:abSize-1,iAtom),ne)
             END IF
             CALL timestop("gemm")

             CALL timestart("local orbitals")
             ! Treatment of local orbitals
             DO lo = 1, atoms%nlo(iType)
                DO nkvec = 1, lapw%nkvec(lo,iAtom)
                   iLAPW = lapw%kvec(nkvec,lo,iAtom)
                   fg(:) = MERGE(lapw%gvec(:,iLAPW,iintsp),lapw%gvec(:,iLAPW,jspin),noco%l_ss) + qss
                   fk = lapw%bkpt + fg(:)
                   tmk = tpi_const * DOT_PRODUCT(fk(:),atoms%taual(:,iAtom))
                   phase = CMPLX(COS(tmk),SIN(tmk))

                   IF (oneD%odi%d1) THEN
                      inap = oneD%ods%ngopr(iAtom)
                      fkr = MATMUL(oneD%ods%mrot(:,:,inap),fk(:))
                      fgr = MATMUL(oneD%ods%mrot(:,:,inap),fg(:))
                   ELSE
                      nap = sym%ngopr(iAtom)
                      inap = sym%invtab(nap)
                      fkr = MATMUL(sym%mrot(:,:,inap),fk(:))
                      fgr = MATMUL(sym%mrot(:,:,inap),fg(:))
                   END IF
                   fkp = MATMUL(fkr,cell%bmat)
                   fgp = MATMUL(fgr,cell%bmat)

                   CALL ylm4(atoms%lmax(iType),fkp,ylm)
                   CALL abclocdn(atoms,sym,noco,lapw,cell,ccchi(:,jspin),iintsp,phase,ylm,iType,iAtom,iLAPW,nkvec,&
                                 lo,ne,alo1(:,jspin),blo1(:,jspin),clo1(:,jspin),acof,bcof,ccof,zMat,l_force,fgp,force)
                END DO
             END DO ! loop over LOs
             CALL timestop("local orbitals")

             ! Treatment of inversion symmetric atoms for noco%l_soc.AND.sym%invs
             ! (The complementary case is treated far below)
             IF (noco%l_soc.AND.sym%invs) THEN
                CALL timestart("invsym atoms")
                IF (sym%invsat(iAtom).EQ.1) THEN
                   jatom = sym%invsatnr(iAtom)
                   DO iLAPW = 1,nvmax
                      DO l = 0,atoms%lmax(iType)
                         ll1 = l* (l+1)
                         DO m = -l,l
                            lm = ll1 + m
                            lmp = ll1 - m
                            inv_f = (-1)**(l-m)
                            c_1 = abCoeffs(iLAPW,lm+1) * inv_f
                            c_2 = abCoeffs(iLAPW,lm+1+abSize) * inv_f
                            CALL zaxpy(ne,c_1,work_c(:,iLAPW),1, acof(:,lmp,jatom),1)
                            CALL zaxpy(ne,c_2,work_c(:,iLAPW),1, bcof(:,lmp,jatom),1)
                         END DO
                      END DO
                   END DO ! loop over LAPWs
                END IF ! IF (sym%invsat(iAtom).EQ.1)
                CALL timestop("invsym atoms")
             END IF ! IF (noco%l_soc.AND.sym%invs)

             ! Force contributions
             IF (atoms%l_geo(iType).AND.l_force) THEN
                CALL timestart("force contributions")
                DO iLAPW = 1,nvmax

                   fg(:) = MERGE(lapw%gvec(:,iLAPW,iintsp),lapw%gvec(:,iLAPW,jspin),noco%l_ss) + qss
                   fk = lapw%bkpt + fg(:)
                   s2h = 0.5 * DOT_PRODUCT(fk,MATMUL(cell%bbmat,fk))
                   s2h_e(:ne) = s2h-eig(:ne)
                   IF (oneD%odi%d1) THEN
                      inap = oneD%ods%ngopr(iAtom)
                      fgr = MATMUL(oneD%ods%mrot(:,:,inap),fg(:))
                   ELSE
                      nap = sym%ngopr(iAtom)
                      inap = sym%invtab(nap)
                      fgr = MATMUL(sym%mrot(:,:,inap),fg(:))
                   END IF
                   fgp = MATMUL(fgr,cell%bmat)

                   DO l = 0,atoms%lmax(iType)
                      ll1 = l* (l+1)
                      DO m = -l,l
                         lm = ll1 + m
                         c_1 = CONJG(abCoeffs(iLAPW,lm+1))
                         c_2 = CONJG(abCoeffs(iLAPW,lm+1+abSize))

                         IF (zmat%l_real) THEN
                            force%e1cof(:ne,lm,iAtom) = force%e1cof(:ne,lm,iAtom) + c_1 * work_r(:ne,iLAPW) * s2h_e(:ne)
                            force%e2cof(:ne,lm,iAtom) = force%e2cof(:ne,lm,iAtom) + c_2 * work_r(:ne,iLAPW) * s2h_e(:ne)
                            DO i = 1,3
                               force%aveccof(i,:ne,lm,iAtom) = force%aveccof(i,:ne,lm,iAtom) + c_1 * work_r(:ne,iLAPW) * fgp(i)
                               force%bveccof(i,:ne,lm,iAtom) = force%bveccof(i,:ne,lm,iAtom) + c_2 * work_r(:ne,iLAPW) * fgp(i)
                            END DO
                         ELSE
                            force%e1cof(:ne,lm,iAtom) = force%e1cof(:ne,lm,iAtom) + c_1 * work_c(:ne,iLAPW) * s2h_e(:ne)
                            force%e2cof(:ne,lm,iAtom) = force%e2cof(:ne,lm,iAtom) + c_2 * work_c(:ne,iLAPW) * s2h_e(:ne)
                            DO i = 1,3
                               force%aveccof(i,:ne,lm,iAtom) = force%aveccof(i,:ne,lm,iAtom) + c_1 * work_c(:ne,iLAPW) * fgp(i)
                               force%bveccof(i,:ne,lm,iAtom) = force%bveccof(i,:ne,lm,iAtom) + c_2 * work_c(:ne,iLAPW) * fgp(i)
                            END DO
                         END IF

                         IF (noco%l_soc.AND.sym%invs) THEN
                            IF (sym%invsat(iAtom).EQ.1) THEN
                               jatom = sym%invsatnr(iAtom)
                               lmp = ll1 - m
                               inv_f = (-1)**(l-m)
                               c_1 =  CONJG(c_1) * inv_f
                               c_2 =  CONJG(c_2) * inv_f
                               CALL zaxpy(ne,c_1,work_c(:,iLAPW)*s2h_e(:),1, force%e1cof(1,lmp,jatom),1)
                               CALL zaxpy(ne,c_2,work_c(:,iLAPW)*s2h_e(:),1, force%e2cof(1,lmp,jatom),1)
                               DO i = 1,3
                                  CALL zaxpy(ne,c_1,work_c(:,iLAPW)*fgp(i),1, force%aveccof(i,1,lmp,jatom),3)
                                  CALL zaxpy(ne,c_2,work_c(:,iLAPW)*fgp(i),1, force%bveccof(i,1,lmp,jatom),3)
                               END DO
                            END IF
                         END IF
                      END DO ! loop over m
                   END DO ! loop over l
                END DO ! loop over LAPWs
                CALL timestop("force contributions")
             END IF

             IF (zmat%l_real) THEN
                DEALLOCATE(work_r)
             ELSE
                DEALLOCATE(work_c)
             END IF

          END IF  ! invsatom == ( 0 v 1 )
       END DO ! loop over interstitial spin
    END DO ! loop over atoms

    ! Treatment of atoms inversion symmetric to others
    IF (noco%l_soc.AND.sym%invs) THEN

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
    ELSE
       DO iAtom = 1, atoms%nat
          iType = atoms%itype(iAtom)
          IF (sym%invsat(iAtom).EQ.1) THEN
             CALL timestart("invsym atoms")
             jAtom = sym%invsatnr(iAtom)
             phase = EXP(tpi_const*ImagUnit*DOT_PRODUCT(atoms%taual(:,jAtom) + atoms%taual(:,iAtom),lapw%bkpt))
             DO ilo = 1,atoms%nlo(iType)
                l = atoms%llo(ilo,iType)
                DO m = -l,l
                   inv_f = (-1.0)**(m+l)
                   DO ie = 1,ne
                      ccof(m,ie,ilo,jatom) = inv_f * phase * CONJG( ccof(-m,ie,ilo,iatom))
                      IF(l_force) THEN
                         force%acoflo(m,ie,ilo,jatom) = inv_f * phase * CONJG(force%acoflo(-m,ie,ilo,iatom))
                         force%bcoflo(m,ie,ilo,jatom) = inv_f * phase * CONJG(force%bcoflo(-m,ie,ilo,iatom))
                         force%cveccof(:,m,ie,ilo,jatom) = -inv_f * phase * CONJG(force%cveccof(:,-m,ie,ilo,iatom))
                      END IF
                   END DO
                END DO
             END DO
             DO l = 0,atoms%lmax(iType)
                ll1 = l* (l+1)
                DO m =-l,l
                   lm  = ll1 + m
                   lmp = ll1 - m
                   inv_f = (-1.0)**(m+l)
                   acof(:ne,lm,jAtom) = inv_f * phase * CONJG(acof(:ne,lmp,iAtom))
                   bcof(:ne,lm,jAtom) = inv_f * phase * CONJG(bcof(:ne,lmp,iAtom))
                   IF (atoms%l_geo(iType).AND.l_force) THEN
                      force%e1cof(:ne,lm,jAtom) = inv_f * phase * CONJG(force%e1cof(:ne,lmp,iAtom))
                      force%e2cof(:ne,lm,jAtom) = inv_f * phase * CONJG(force%e2cof(:ne,lmp,iAtom))
                      force%aveccof(:,:ne,lm,jAtom) = -inv_f * phase * CONJG(force%aveccof(:,:ne,lmp,iAtom))
                      force%bveccof(:,:ne,lm,jAtom) = -inv_f * phase * CONJG(force%bveccof(:,:ne,lmp,iAtom))
                   END IF
                END DO
             END DO
             CALL timestop("invsym atoms")
          END IF
       END DO
    END IF

    CALL timestop("abcof")

  END SUBROUTINE abcof
END MODULE m_abcof
