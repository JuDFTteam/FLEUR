!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_symmetrizeh

! symmetrize the Hamiltonian according to the symmetry operations of the little group of k

CONTAINS

SUBROUTINE symmetrizeh(atoms,bk,DIMENSION,jsp,lapw,gpt,sym,kveclo,cell,nsymop,psym,hmat)

   USE m_constants
   USE m_types

   IMPLICIT NONE

   TYPE(t_dimension), INTENT(IN)    :: DIMENSION
   TYPE(t_sym),       INTENT(IN)    :: sym
   TYPE(t_cell),      INTENT(IN)    :: cell
   TYPE(t_atoms),     INTENT(IN)    :: atoms
   TYPE(t_lapw),      INTENT(IN)    :: lapw
   TYPE(T_mat),       INTENT(INOUT) :: hmat

   ! scalars
   INTEGER,           INTENT(IN)    :: nsymop, jsp

   ! arrays
   INTEGER,           INTENT(IN)    :: gpt(:,:)!(3,lapw%nv)
   INTEGER,           INTENT(IN)    :: kveclo(atoms%nlotot)
   INTEGER,           INTENT(IN)    :: psym(nsymop)
   REAL,              INTENT(IN)    :: bk(3)

   ! local scalars
   INTEGER               ::  ilotot,itype,itype1,ilo,ilo1
   INTEGER               ::  iatom,iatom1,iiatom,iiatom1
   INTEGER               ::  i,ieq,ieq1,m
   INTEGER               ::  igpt_lo,igpt_lo1,igpt_lo2,igpt1_lo1
   INTEGER               ::  igpt1_lo2,isym,iop,ic,ic1,ic2
   INTEGER               ::  igpt,igpt1,igpt2,igpt3
   INTEGER               ::  invsfct,invsfct1,idum
   INTEGER               ::  l ,l1,lm,j,ok,ratom,ratom1,nrgpt
   COMPLEX,PARAMETER     ::  img=(0d0,1d0)
   COMPLEX               ::  cdum,cdum2

   ! local arrays
   INTEGER               ::  l_lo(atoms%nlotot)
   INTEGER               ::  itype_lo(atoms%nlotot)
   INTEGER               ::  gpt_lo(3,atoms%nlotot),gpthlp(3),g(3)
   INTEGER               ::  indx(DIMENSION%nbasfcn,DIMENSION%nbasfcn)
   INTEGER               ::  lo_indx(atoms%nlod,atoms%nat)
   INTEGER               ::  rot(3,3,nsymop),rrot(3,3,nsymop)

   INTEGER,ALLOCATABLE   ::  pointer_apw(:,:)
   INTEGER,ALLOCATABLE   ::  ipiv(:)
   INTEGER,ALLOCATABLE   ::  map(:,:)

   REAL                  ::  rtaual(3),kghlp(3)
   REAL                  ::  rotkpthlp(3),rotkpt(3)
   REAL                  ::  trans(3,nsymop)
   COMPLEX,ALLOCATABLE   ::  c_lo(:,:,:,:),c_rot(:,:,:,:,:),y(:)
   COMPLEX,ALLOCATABLE   ::  cfac(:,:),chelp(:,:)

   LOGICAL               ::  ldum(lapw%nv(jsp)+atoms%nlotot,lapw%nv(jsp)+atoms%nlotot)

   ! calculate rotations in reciprocal space
   DO isym = 1,nsymop
      iop = psym(isym)
      IF(iop.LE.sym%nop) THEN
         rrot(:,:,isym) = TRANSPOSE(sym%mrot(:,:,sym%invtab(iop)))
      ELSE
         rrot(:,:,isym) = -TRANSPOSE(sym%mrot(:,:,sym%invtab(iop-sym%nop)))
      END IF
   END DO

   ! calculate rotations in real space (internal coordinates)
   DO isym = 1,nsymop
      iop = psym(isym)
      IF(iop.LE.sym%nop) THEN
         rot(:,:,isym) = sym%mrot(:,:,iop)
         trans(:,isym) = sym%tau(:,iop)
      ELSE
         rot(:,:,isym) = sym%mrot(:,:,iop-sym%nop)
         trans(:,isym) = sym%tau(:,iop-sym%nop)
      END IF
   END DO

   ! caclulate mapping of atoms
   ALLOCATE(map(nsymop,atoms%nat))
   map = 0
   iatom  = 0
   iiatom = 0
   DO itype = 1,atoms%ntype
      DO ieq = 1,atoms%neq(itype)
         iatom = iatom + 1
         DO isym = 1,nsymop
            rtaual = MATMUL(rot(:,:,isym),atoms%taual(:,iatom))+trans(:,isym)
            iatom1 = 0
            DO ieq1 = 1,atoms%neq(itype)
               IF(ALL(ABS(MODULO(rtaual-atoms%taual(:,iiatom + ieq1)+1d-12,1d0)).LT.1d-10)) THEN  !The 1d-12 is a dirty fix.
                  iatom1 = iiatom + ieq1
               END IF
            END DO
            IF(iatom1.EQ.0) STOP 'symmetrizeh_new: error finding rotated atomic position'
            map(isym,iatom) = iatom1
         END DO
      END DO
      iiatom = iiatom + atoms%neq(itype)
   END DO

   ! initialze pointer_apw and the apw part of cfac
   ALLOCATE(pointer_apw(lapw%nv(jsp),nsymop),cfac(lapw%nv(jsp)+atoms%nlotot,nsymop), stat= ok)
   IF(ok.NE.0) STOP 'symmetrizeh_new: failure allocation pointer_apw,cfac'

   pointer_apw = 0
   cfac = 0

   DO isym = 1, nsymop

      ! determine vector g, which map Rk back into BZ
      ! Rk - G = k => G = Rk-k
      rotkpt = MATMUL( rrot(:,:,isym), bk(:) )
      g = NINT(rotkpt-bk)

      DO igpt = 1, lapw%nv(jsp)
         !rotate G vector corresponding to isym
         gpthlp = MATMUL(rrot(:,:,isym),gpt(:,igpt)) + g
         ! determine number of gpthlp
         nrgpt = 0
         DO i = 1, lapw%nv(jsp)
            IF(MAXVAL( ABS( gpthlp - gpt(:,i) ) ) .LE. 1E-06) THEN
               nrgpt = i
               EXIT
            END IF
         END DO
         IF(nrgpt.EQ.0) THEN
            PRINT *,igpt
            PRINT *,gpt(:,igpt)
            PRINT *,gpthlp
            PRINT *,g
            PRINT *,bk
            DO i=1,lapw%nv(jsp)
               WRITE(6,*) i,gpt(:,i)
            ENDDO
            STOP 'symmetrizeh_new: rotated G point not found'
         END IF
         pointer_apw(igpt,isym) = nrgpt
         cfac(igpt,isym) = EXP(-2*pi_const*img* (dot_PRODUCT( bk(:)+gpthlp(:),trans(:,isym) ) ) )
      END DO
   END DO

   ! average apw-part of symmetry-equivalent matrix elements

   ldum = .TRUE.
   DO i = 1, lapw%nv(jsp)
      DO j = 1, i
         cdum = 0
         ic = 0
         DO isym = 1, nsymop
            iop = psym(isym)
            igpt = pointer_apw(i,isym)
            igpt1 = pointer_apw(j,isym)

            IF(iop.LE.sym%nop) THEN
               IF((igpt.NE.0).AND.(igpt1.NE.0)) THEN
                  ic = ic + 1
                  IF (hmat%l_real) THEN
                     cdum = cdum + CONJG(cfac(i,isym))*hmat%data_r(igpt1,igpt)*cfac(j,isym)
                  ELSE
                     cdum = cdum + CONJG(cfac(i,isym))*hmat%data_c(igpt1,igpt)*cfac(j,isym)
                  END IF
               END IF
            ELSE
               IF((igpt.NE.0).AND.(igpt1.NE.0)) THEN
                  ic = ic + 1
                  IF (hmat%l_real) THEN
                     cdum = cdum + CONJG(CONJG(cfac(i,isym))*hmat%data_r(igpt1,igpt)*cfac(j,isym))
                  ELSE
                     cdum = cdum + CONJG(CONJG(cfac(i,isym))*hmat%data_c(igpt1,igpt)*cfac(j,isym))
                  END IF
               END IF
            END IF
         END DO

         cdum = cdum!/ic
         DO isym = 1, nsymop
            iop = psym(isym)
            igpt = pointer_apw(i,isym)
            igpt1 = pointer_apw(j,isym)
            IF ((igpt.EQ.0).OR.(igpt1.EQ.0)) CYCLE
            IF (igpt1.GT.igpt) CYCLE
            IF (ldum(igpt,igpt1)) THEN
               IF (hmat%l_real) THEN
                  IF (iop.LE.sym%nop) THEN
                     hmat%data_r(igpt1,igpt) = cdum/(CONJG(cfac(i,isym))*cfac(j,isym))
                     ldum(igpt,igpt1) = .FALSE.
                  ELSE
                     hmat%data_r(igpt1,igpt) = CONJG(cdum/(CONJG(cfac(i,isym))*cfac(j,isym)))
                     ldum(igpt,igpt1) = .FALSE.
                  END IF
               ELSE
                  IF (iop.LE.sym%nop) THEN
                     hmat%data_c(igpt1,igpt) = cdum/(CONJG(cfac(i,isym))*cfac(j,isym))
                     ldum(igpt,igpt1) = .FALSE.
                  ELSE
                     hmat%data_c(igpt1,igpt) = CONJG(cdum/(CONJG(cfac(i,isym))*cfac(j,isym)))
                     ldum(igpt,igpt1)    = .FALSE.
                  END IF
               END IF
            END IF
         END DO
      END DO
   END DO


   ! average lo-part of matrix elements
   IF (ANY(atoms%nlo.NE.0)) THEN

      ! preparations
      ilotot = 0
      iatom  = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1,atoms%neq(itype)
            iatom = iatom + 1
            IF ((atoms%invsat(iatom).EQ.0).OR.(atoms%invsat(iatom).EQ.1)) THEN
               IF (atoms%invsat(iatom).EQ.0) invsfct = 1
               IF (atoms%invsat(iatom).EQ.1) invsfct = 2
               DO ilo = 1, atoms%nlo(itype)
                  l = atoms%llo(ilo,itype)
                  DO m = 1, invsfct*(2*l+1)
                     ilotot = ilotot + 1
                     l_lo(ilotot) = l
                     itype_lo(ilotot) = itype
                     gpt_lo(:,ilotot) = gpt(:,kveclo(ilotot))
                  END DO
               END DO
            END IF
         END DO
      END DO

      ! calculate expansion coefficients for local orbitals
      IF (hmat%l_real) THEN
         ALLOCATE(c_lo(4*MAXVAL(l_lo)+2,4*MAXVAL(l_lo)+2,atoms%nlod,atoms%nat), stat = ok)
      ELSE
         ALLOCATE(c_lo(2*MAXVAL(l_lo)+1,2*MAXVAL(l_lo)+1,atoms%nlod,atoms%nat), stat = ok)
      END IF

      IF (ok.NE.0) STOP 'symmetrizeh_new: failure allocation c_lo'

      iatom = 0
      ilotot = 0
      lo_indx = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1
            IF ((atoms%invsat(iatom).EQ.0).OR.(atoms%invsat(iatom).EQ.1)) THEN
               IF(atoms%invsat(iatom).EQ.0) invsfct = 1
               IF(atoms%invsat(iatom).EQ.1) invsfct = 2

               DO ilo = 1, atoms%nlo(itype)
                  l = atoms%llo(ilo,itype)
                  ALLOCATE(y((l+1)**2))
                  lo_indx(ilo,iatom) = ilotot + 1

                  DO igpt_lo = 1, invsfct*(2*l+1)
                     ilotot = ilotot + 1
                     kghlp  = bk(:) + gpt_lo(:,ilotot)

                     !generate spherical harmonics
                     CALL harmonicsr(y,MATMUL(kghlp,cell%bmat),l)

                     lm   = l**2
                     cdum = EXP(2*pi_const*img*dot_PRODUCT(kghlp,atoms%taual(:,iatom)))

                     IF (invsfct.EQ.1) THEN
                        DO m = 1,2*l+1
                           lm = lm + 1
                           c_lo(m,igpt_lo,ilo,iatom) = cdum*CONJG(y(lm))
                        END DO
                     ELSE
                        DO m = 1, 2*l+1
                           lm = lm + 1
                           c_lo(m,igpt_lo,ilo,iatom) = cdum *CONJG(y(lm))
                           c_lo(4*l+3-m,igpt_lo,ilo,iatom) = (-1)**(m-1)*CONJG(cdum)*y(lm)
                        END DO
                     END IF
                  END DO
                  DEALLOCATE( y )
               END DO
            END IF
         END DO
      END DO

      IF (ilotot.NE.atoms%nlotot) STOP 'symmetrizeh_new: failure counting local orbitals(ilotot)'

      IF (hmat%l_real) THEN
         ALLOCATE(c_rot(4*MAXVAL(l_lo)+2,4*MAXVAL(l_lo)+2,atoms%nlod,atoms%nat, nsymop))
         ALLOCATE(chelp(4*MAXVAL(l_lo)+2,4*MAXVAL(l_lo)+2))
      ELSE
         ALLOCATE(c_rot(2*MAXVAL(l_lo)+1,2*MAXVAL(l_lo)+1,atoms%nlod,atoms%nat, nsymop))
         ALLOCATE(chelp(2*MAXVAL(l_lo)+1,2*MAXVAL(l_lo)+1))
      END IF

      DO isym = 1,nsymop
         ! determine vector g, which map Rk back into BZ
         ! Rk - G = k => G = Rk-k
         rotkpt = MATMUL( rrot(:,:,isym), bk(:) )
         g = NINT(rotkpt-bk)

         ilotot = 0
         iatom = 0
         DO itype = 1,atoms%ntype
            DO ieq = 1,atoms%neq(itype)
               iatom = iatom + 1
               ratom = map(isym,iatom)

               IF ((atoms%invsat(iatom).EQ.0).OR.(atoms%invsat(iatom).EQ.1)) THEN
                  IF (atoms%invsat(iatom).EQ.0) invsfct = 1
                  IF (atoms%invsat(iatom).EQ.1) THEN
                     IF (atoms%invsat(ratom).EQ.2) THEN
                        ratom = sym%invsatnr(ratom)
                     END IF
                     invsfct = 2
                  END IF

                  DO ilo = 1,atoms%nlo(itype)
                     l = atoms%llo(ilo,itype)
                     ALLOCATE(y((l+1)**2))

                     DO igpt_lo = 1,invsfct*(2*l+1)
                        ilotot = ilotot + 1

                        !rotate G_lo corresponding to iop
                        gpthlp = MATMUL(rrot(:,:,isym),gpt_lo(:,ilotot)) + g
                        kghlp  = bk(:) + MATMUL(rrot(:,:,isym),gpt_lo(:,ilotot)) + g

                        !generate spherical harmonics
                        CALL harmonicsr(y,MATMUL(kghlp,cell%bmat),l)

                        lm = l**2
                        cdum = EXP(2*pi_const*img* dot_PRODUCT(kghlp,atoms%taual(:,ratom)))
                        IF (invsfct.EQ.1) THEN
                           DO m = 1,2*l+1
                              lm = lm + 1
                              c_rot(m,igpt_lo,ilo,ratom,isym) = cdum*CONJG(y(lm))
                           END DO
                        ELSE
                           DO m = 1,2*l+1
                              lm = lm + 1
                              c_rot(m,igpt_lo,ilo,ratom,isym) = cdum *CONJG(y(lm))
                              c_rot(4*l+3-m,igpt_lo,ilo,ratom,isym) = (-1)**(m-1)*(CONJG(cdum))* y(lm)
                           END DO
                        END IF
                        cfac(lapw%nv(jsp)+ilotot,isym) = EXP(-2*pi_const*img* (dot_PRODUCT( kghlp,trans(:,isym) ) ) )
                     END DO

                     idum = invsfct*(2*l+1)

                     ALLOCATE(ipiv(idum))

                     chelp(:idum,:idum) = c_lo(:idum,:idum,ilo,ratom)

                     CALL ZGESV(idum,idum,chelp(:idum,:idum),idum,ipiv, c_rot(:idum,:idum,ilo,ratom,isym),idum,ok )

                     IF(ok.NE.0) STOP 'symmetrizeh_new: failure zgesv'
                     DEALLOCATE(ipiv,y)
                  END DO
               END IF
            END DO
         END DO
      END DO

      ! symmetrize local-orbital-apw-part of the matrix
      i = 0
      iatom = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1
            IF ((atoms%invsat(iatom).EQ.0).OR.(atoms%invsat(iatom).EQ.1)) THEN
               IF (atoms%invsat(iatom).EQ.0) invsfct = 1
               IF (atoms%invsat(iatom).EQ.1) invsfct = 2

               DO ilo = 1, atoms%nlo(itype)
                  l = atoms%llo(ilo,itype)
                  DO igpt = 1, invsfct*(2*l+1)
                     i = i + 1
                     DO j = 1, lapw%nv(jsp)
                        cdum = 0
                        ic = 0
                        DO isym = 1,nsymop
                           ic = ic + 1
                           iop = psym(isym)
                           ratom = map(isym,iatom)
                           IF (invsfct.EQ.2) THEN
                              IF (atoms%invsat(ratom).EQ.2) THEN
                                 ratom = sym%invsatnr(ratom)
                              END IF
                           END IF

                           igpt_lo1 = lo_indx(ilo,ratom)
                           igpt_lo2 = igpt_lo1 + invsfct*2*l
                           IF (invsfct.EQ.2) igpt_lo2 = igpt_lo2 + 1
                           igpt1 = pointer_apw(j,isym)

                           cdum2 = 0
                           ic1 = 0
                           DO igpt2 = igpt_lo1, igpt_lo2
                              ic1 = ic1 + 1
                              IF (hmat%l_real) THEN
                                 cdum2 = cdum2 + CONJG(c_rot(ic1,igpt,ilo,ratom,isym)) * hmat%data_r(igpt1,lapw%nv(jsp)+igpt2)
                              ELSE
                                 cdum2 = cdum2 + CONJG(c_rot(ic1,igpt,ilo,ratom,isym)) * hmat%data_c(igpt1,lapw%nv(jsp)+igpt2)
                              END IF
                           END DO

                           IF (iop.LE.sym%nop) THEN
                              cdum = cdum + cdum2*CONJG(cfac(lapw%nv(jsp)+i,isym)) * cfac(j,isym)
                           ELSE
                              cdum = cdum + CONJG(cdum2*CONJG(cfac(lapw%nv(jsp)+i,isym)) * cfac(j,isym))
                           END IF
                        END DO
                        IF (hmat%l_real) THEN
                           hmat%data_r(j,lapw%nv(jsp)+i) = cdum!/ic
                        ELSE
                           hmat%data_c(j,lapw%nv(jsp)+i) = cdum!/ic
                        END IF
                     END DO
                  END DO
               END DO
            END IF
         END DO
      END DO

      !lo's - lo's
      i = 0
      iatom = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1
            IF ((atoms%invsat(iatom).EQ.0).OR.(atoms%invsat(iatom).EQ.1)) THEN
               IF(atoms%invsat(iatom).EQ.0) invsfct = 1
               IF(atoms%invsat(iatom).EQ.1) invsfct = 2

               DO ilo = 1, atoms%nlo(itype)
                  l = atoms%llo(ilo,itype)
                  DO igpt = 1, invsfct*(2*l+1)
                     i = i + 1
                     j = 0
                     iatom1 = 0
                     DO itype1 = 1, atoms%ntype
                        DO ieq1 = 1, atoms%neq(itype1)
                           iatom1 = iatom1 + 1
                           IF ((atoms%invsat(iatom1).EQ.0).OR.(atoms%invsat(iatom1).EQ.1)) THEN
                              IF(atoms%invsat(iatom1).EQ.0) invsfct1 = 1
                              IF(atoms%invsat(iatom1).EQ.1) invsfct1 = 2

                              DO ilo1 = 1, atoms%nlo(itype1)
                                 l1 = atoms%llo(ilo1,itype1)
                                 DO igpt1 = 1, invsfct1*(2*l1+1)
                                    j = j + 1
                                    IF (j.GT.i) CYCLE
                                    cdum = 0
                                    ic = 0
                                    DO isym = 1, nsymop
                                       iop = psym(isym)
                                       ratom = map(isym,iatom )
                                       ratom1 = map(isym,iatom1)

                                       IF (invsfct.EQ.2) THEN
                                          IF (atoms%invsat(ratom).EQ.2) THEN
                                             ratom = sym%invsatnr(ratom)
                                          END IF
                                       END IF
                                       IF (invsfct1.EQ.2) THEN
                                          IF (atoms%invsat(ratom1).EQ.2) THEN
                                             ratom1 = sym%invsatnr(ratom1)
                                          END IF
                                       END IF

                                       igpt_lo1  = lo_indx(ilo,ratom )
                                       igpt_lo2  = igpt_lo1 + invsfct*2*l
                                       IF (invsfct.EQ.2) igpt_lo2 = igpt_lo2 + 1

                                       igpt1_lo1 = lo_indx(ilo1,ratom1)
                                       igpt1_lo2 = igpt1_lo1 + invsfct1*2*l1
                                       IF (invsfct1.EQ.2) igpt1_lo2 = igpt1_lo2 + 1

                                       cdum2 = 0
                                       ic1 = 0
                                       DO igpt2 = igpt_lo1, igpt_lo2
                                          ic1 = ic1 + 1
                                          ic2 = 0
                                          DO igpt3 = igpt1_lo1, igpt1_lo2
                                             ic2 = ic2 + 1
                                             IF (hmat%l_real) THEN
                                                cdum2 = cdum2 + CONJG(c_rot(ic1,igpt,ilo,ratom,isym)) *&
                                                                hmat%data_r(lapw%nv(jsp)+igpt3,lapw%nv(jsp)+igpt2) *&
                                                                c_rot(ic2,igpt1,ilo1,ratom1,isym)
                                             ELSE
                                                 cdum2 = cdum2 + CONJG(c_rot(ic1,igpt,ilo,ratom,isym)) *&
                                                                 hmat%data_c(lapw%nv(jsp)+igpt3,lapw%nv(jsp)+igpt2) *&
                                                                 c_rot(ic2,igpt1,ilo1,ratom1,isym)
                                             END IF
                                          END DO
                                       END DO
                                       ic = ic + 1
                                       IF (iop.LE.sym%nop) THEN
                                          cdum = cdum + cdum2*CONJG(cfac(lapw%nv(jsp)+i,isym))*cfac(lapw%nv(jsp)+j,isym)
                                       ELSE  
                                          cdum = cdum + CONJG(cdum2*CONJG(cfac(lapw%nv(jsp)+i,isym))*cfac(lapw%nv(jsp)+j,isym))
                                       END IF
                                    END DO
                                    IF (hmat%l_real) THEN
                                       hmat%data_r(lapw%nv(jsp)+j,lapw%nv(jsp)+i) = cdum!/ic
                                    ELSE
                                       hmat%data_c(lapw%nv(jsp)+j,lapw%nv(jsp)+i) = cdum!/ic
                                    END IF
                                 END DO  ! igpt_lo1
                              END DO  ! ilo1
                           END IF
                        END DO  !ieq1
                     END DO  !itype1
                  END DO  ! igpt_lo
               END DO  ! ilo
            END IF
         END DO  !ieq
      END DO  !itype

   END IF ! ANY(atoms%nlo.NE.0)

   CONTAINS

   ! Returns the spherical harmonics Y_lm(^rvec) for l = 0,...,ll in Y(1,...,(ll+1)**2).
   SUBROUTINE harmonicsr(Y,rvec,ll)
      IMPLICIT NONE
      INTEGER,INTENT(IN)    :: ll
      REAL,INTENT(IN)       :: rvec(3)
      COMPLEX,INTENT(OUT)   :: Y((ll+1)**2)
      REAL                  :: stheta,ctheta,sphi,cphi,r,rvec1(3)
      INTEGER               :: l ,lm
      COMPLEX               :: c
      COMPLEX,PARAMETER     :: img=(0d0,1d0)

      Y(1) = 0.282094791773878d0
      IF(ll.EQ.0) RETURN

      stheta = 0
      ctheta = 0
      sphi = 0
      cphi = 0
      r = SQRT(SUM(rvec**2))
      IF (r.GT.1d-16) THEN
         rvec1  = rvec / r
         ctheta = rvec1(3)
         stheta = SQRT(rvec1(1)**2+rvec1(2)**2)
         IF (stheta.GT.1d-16) THEN
            cphi = rvec1(1) / stheta
            sphi = rvec1(2) / stheta
         END IF
      ELSE
         Y(2:) = 0d0
         RETURN
      END IF

      ! define Y,l,-l and Y,l,l
      r = Y(1)
      c = 1
      DO l = 1, ll
         r = r*stheta*SQRT(1d0+1d0/(2*l))
         c = c * (cphi + img*sphi)
         Y(l**2+1) = r*CONJG(c)  ! l,-l
         Y((l+1)**2) = r*c*(-1)**l ! l,l
      END DO

      ! define Y,l,-l+1 and Y,l,l-1
      Y(3) = 0.48860251190292d0*ctheta
      DO l = 2, ll
         r = SQRT(2D0*l+1) * ctheta
         Y(l**2+2) = r*Y((l-1)**2+1) ! l,-l+1
         Y(l*(l+2)) = r*Y(l**2)       ! l,l-1
      END DO

      ! define Y,l,m, |m|<l-1
      DO l = 2, ll
         lm = l**2 + 2
         DO m = -l+2, l-2
            lm = lm + 1
            Y(lm) = SQRT((2d0*l+1)/(l+m)/(l-m)) * (SQRT(2d0*l-1)*ctheta*Y(lm-2*l)- SQRT((l+m-1d0)*(l-m-1)/(2*l-3))*Y(lm-4*l+2) )
         END DO
      END DO
   END SUBROUTINE harmonicsr

END SUBROUTINE symmetrizeh

END MODULE m_symmetrizeh
