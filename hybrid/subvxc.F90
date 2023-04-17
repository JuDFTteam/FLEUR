!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_subvxc

CONTAINS

   SUBROUTINE subvxc(lapw, bk, input, jsp, vr0, atoms, usdus, mpdata, hybdat,&
                     el, ello, sym, cell, sphhar, stars, xcpot, fmpi,   hmat, vx)

      USE m_judft
      USE m_types
      USE m_constants
      USE m_intgr, ONLY: intgr3
      USE m_gaunt, ONLY: gaunt1
      USE m_wrapper
      USE m_loddop
      USE m_radflo
      USE m_radfun
      USE m_abcof3

      IMPLICIT NONE

      CLASS(t_xcpot), INTENT(IN)    :: xcpot
      TYPE(t_mpi), INTENT(IN)    :: fmpi
       
      TYPE(t_mpdata), intent(inout) :: mpdata
      TYPE(t_hybdat), INTENT(IN) :: hybdat
      TYPE(t_input), INTENT(IN)    :: input
      TYPE(t_sym), INTENT(IN)    :: sym
      TYPE(t_stars), INTENT(IN)    :: stars
      TYPE(t_cell), INTENT(IN)    :: cell
      TYPE(t_sphhar), INTENT(IN)    :: sphhar
      TYPE(t_atoms), INTENT(IN)    :: atoms
      TYPE(t_lapw), INTENT(IN)    :: lapw
      TYPE(t_usdus), INTENT(INOUT) :: usdus
      TYPE(t_potden), INTENT(IN)    :: vx
      TYPE(t_mat), INTENT(INOUT) :: hmat

      ! Scalar Arguments
      INTEGER, INTENT(IN) :: jsp

      ! Array Arguments
      REAL, INTENT(IN) :: vr0(:,:,:)               ! just for radial functions
      REAL, INTENT(IN) :: el(0:atoms%lmaxd, atoms%ntype, input%jspins)
      REAL, INTENT(IN) :: ello(:,:,:)
      REAL, INTENT(IN) :: bk(:)

      ! Local Scalars
      INTEGER               ::  ic, indx, m, ig1, ig1_loc, ig2, n, n_loc, nn
      INTEGER               ::  nlharm, nnbas, typsym, lm
      INTEGER               ::  noded, nodeu
      INTEGER               ::  nbasf0
      INTEGER               ::  i, j, j_loc, l, ll, l1, l2, m1, m2, j1, j2
      INTEGER               ::  ok, p1, p2, lh, mh, pp1, pp2
      INTEGER               ::  igrid, itype, ilharm, istar
      INTEGER               ::  ineq, iatom, ilo, ilop, ieq, icentry
      INTEGER               ::  ikvecat, ikvecprevat, invsfct, ikvec, ikvecp
      INTEGER               ::  lp, mp, pp
      REAL                  ::  a_ex
      REAL                  ::  wronk
      COMPLEX               ::  rc, rr

      ! Local Arrays
      INTEGER               ::  gg(3), x, x_loc, y, iband
      INTEGER               ::  pointer_lo(atoms%nlod, atoms%ntype)

      REAL                  ::  integ(0:sphhar%nlhd, maxval(mpdata%num_radfun_per_l), 0:atoms%lmaxd, maxval(mpdata%num_radfun_per_l), 0:atoms%lmaxd)
      REAL                  ::  grid(atoms%jmtd)
      REAL                  ::  vr(atoms%jmtd, 0:sphhar%nlhd)
      REAL                  ::  f(atoms%jmtd, 2, 0:atoms%lmaxd), g(atoms%jmtd, 2, 0:atoms%lmaxd)
      REAL                  ::  flo(atoms%jmtd, 2, atoms%nlod)
      REAL                  ::  uuilon(atoms%nlod, atoms%ntype), duilon(atoms%nlod, atoms%ntype)
      REAL                  ::  ulouilopn(atoms%nlod, atoms%nlod, atoms%ntype)

      REAL                  ::  bas1(atoms%jmtd, maxval(mpdata%num_radfun_per_l), 0:atoms%lmaxd, atoms%ntype)
      REAL                  ::  bas2(atoms%jmtd, maxval(mpdata%num_radfun_per_l), 0:atoms%lmaxd, atoms%ntype)

      COMPLEX               ::  vpw(stars%ng3)
      COMPLEX               ::  carr(hybdat%maxlmindx, lapw%dim_nvd()), carr1(lapw%dim_nvd(), lapw%dim_nvd())
      COMPLEX, ALLOCATABLE  ::  ahlp(:, :, :), bhlp(:, :, :)
      type(t_mat)           ::  vxc, vrmat
      type(t_mat), allocatable ::  bascof(:)
      COMPLEX               ::  bascof_lo(3, -atoms%llod:atoms%llod, 4*atoms%llod + 2, atoms%nlod, atoms%nat)

      CALL timestart("subvxc")
      integ = 0.0
      bas1 = 0.0
      bas2 = 0.0

      call vxc%init(hmat)
      call vxc%reset(cmplx_0)

      ! Calculate radial functions
      call timestart("Calculate radial functions")
      mpdata%num_radfun_per_l = 2
      DO itype = 1, atoms%ntype

         ! Generate the radial basis-functions for each l
         ! WRITE (oUnit, '(a,i3,a)') new_LINE('n')//new_LINE('n')//' wavefunction parameters for atom type', itype, ':'
         ! WRITE (oUnit, '(31x,a,32x,a)') 'radial function', 'energy derivative'
         ! WRITE (oUnit, '(a)') '  l    energy            value        '// &
         !    'derivative    nodes          value        derivative    nodes       norm        wronskian'
         DO l = 0, atoms%lmax(itype)
            CALL radfun(l, itype, jsp, el(l, itype, jsp), vr0(:, itype, jsp), atoms, f(1, 1, l), g(1, 1, l), usdus, nodeu, noded, wronk)
!             WRITE (oUnit, FMT=8010) l, el(l, itype, jsp), usdus%us(l, itype, jsp), &
!                usdus%dus(l, itype, jsp), nodeu, usdus%uds(l, itype, jsp), usdus%duds(l, itype, jsp), noded, &
!                usdus%ddn(l, itype, jsp), wronk
         END DO
! 8010     FORMAT(i3, f10.5, 2(5x, 1p, 2e16.7, i5), 1p, 2e16.7)

         bas1(:, 1, :, itype) = f(:, 1, :)
         bas1(:, 2, :, itype) = g(:, 1, :)
         bas2(:, 1, :, itype) = f(:, 2, :)
         bas2(:, 2, :, itype) = g(:, 2, :)

         ! Generate the extra radial basis-functions for the local orbitals, if there are any.
         IF (atoms%nlo(itype) >= 1) THEN
            CALL radflo(atoms, itype, jsp, ello(:,:, jsp), vr0(:, itype, jsp), f, g, fmpi, &
                        usdus, uuilon, duilon, ulouilopn, flo, .TRUE.)

            DO i = 1, atoms%nlo(itype)
               mpdata%num_radfun_per_l(atoms%llo(i, itype), itype) = mpdata%num_radfun_per_l(atoms%llo(i, itype), itype) + 1
               pointer_lo(i, itype) = mpdata%num_radfun_per_l(atoms%llo(i, itype), itype)
               bas1(:, mpdata%num_radfun_per_l(atoms%llo(i, itype), itype), atoms%llo(i, itype), itype) = flo(:, 1, i)
               bas2(:, mpdata%num_radfun_per_l(atoms%llo(i, itype), itype), atoms%llo(i, itype), itype) = flo(:, 2, i)
            END DO
         END IF
      END DO
      call timestop("Calculate radial functions")

      ! Compute APW coefficients

      ! Calculate bascof
      call timestart("Calculate bascof")
      ALLOCATE(ahlp(lapw%dim_nvd(), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), stat=ok)
      IF (ok /= 0) call judft_error('subvxc: error in allocation of ahlp')
      allocate(bhlp(lapw%dim_nvd(), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), stat=ok)
      IF (ok /= 0) call judft_error('subvxc: error in allocation of bhlp')  
      
      allocate(bascof(atoms%nat))
      do ic = 1, atoms%nat 
         call bascof(ic)%alloc(.false., lapw%dim_nvd(), 2*(atoms%lmaxd*(atoms%lmaxd+2) + 1))
      enddo
      
      CALL abcof3(input, atoms, sym, jsp, cell, bk, lapw, usdus,   &
                  1,lapw%nv(jsp), ahlp, bhlp, bascof_lo)

   
      do iband = 1,lapw%nv(jsp)
         do ic = 1,atoms%nat
            itype = atoms%itype(ic)
            indx = 0
            DO l = 0, atoms%lmax(itype)
               ll = l*(l + 1)
               DO M = -l, l
                  lm = ll + M
                  DO i = 1, 2
                     indx = indx + 1
                     IF (i == 1) THEN
                        bascof(ic)%data_c(iband, indx) = ahlp(iband, lm, ic)
                     ELSE IF (i == 2) THEN
                        bascof(ic)%data_c(iband, indx) = bhlp(iband, lm, ic)
                     END IF
                  END DO
               END DO
            END DO
         END DO
      enddo

      deallocate(ahlp, bhlp)
      call timestop("Calculate bascof")


      ! Loop over atom types
      call vrmat%alloc(.false., hybdat%maxlmindx, hybdat%maxlmindx)
      iatom = 0
      DO itype = 1, atoms%ntype

         typsym = sym%ntypsy(atoms%firstAtom(itype))
         nlharm = sphhar%nlh(typsym)

         ! Calculate vxc = vtot - vcoul
         call timestart("Calculate vxc = vtot - vcoul")
         DO l = 0, nlharm
            DO i = 1, atoms%jri(itype)
               IF (l == 0) THEN
                  ! vr(i,0)= vrtot(i,0,itype)*sfp/rmsh(i,itype) - vrcou(i,0,itype,jsp)
                  vr(i, 0) = vx%mt(i, 0, itype, jsp)! * sfp_const / atoms%rmsh(i,itype)
               ELSE ! vxc = vtot - vcoul
                  ! vr(i,l)= vrtot(i,l,itype)-vrcou(i,l,itype,jsp)
                  vr(i, l) = vx%mt(i, l, itype, jsp)
               END IF
            END DO
         END DO
         call timestop("Calculate vxc = vtot - vcoul")

         ! Calculate MT contribution to vxc matrix elements
         ! Precompute auxiliary radial integrals
         call timestart("calc MT contrib aux radial integral")
         DO ilharm = 0, nlharm
            i = 0
            DO l1 = 0, atoms%lmax(itype)
               DO p1 = 1, 2
                  i = i + 1
                  j = 0
                  DO l2 = 0, atoms%lmax(itype)
                     DO p2 = 1, 2
                        j = j + 1
                        IF (j <= i) THEN
                           DO igrid = 1, atoms%jri(itype)
                              grid(igrid) = vr(igrid, ilharm)*(bas1(igrid, p1, l1, itype)*bas1(igrid, p2, l2, itype) + &
                                                               bas2(igrid, p1, l1, itype)*bas2(igrid, p2, l2, itype))
                           END DO
                           CALL intgr3(grid, atoms%rmsh(:, itype), atoms%dx(itype), atoms%jri(itype), integ(ilharm, p1, l1, p2, l2))
                           integ(ilharm, p2, l2, p1, l1) = integ(ilharm, p1, l1, p2, l2)
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END DO
         call timestop("calc MT contrib aux radial integral")

         ! Calculate muffin tin contribution to vxc matrix
         call vrmat%reset(cmplx_0)
         call timestart("calc MT contrib")

         j1 = 0
         DO l1 = 0, atoms%lmax(itype) ! loop: left basis function
            DO m1 = -l1, l1
               DO p1 = 1, 2
                  j1 = j1 + 1
                  j2 = 0
                  DO l2 = 0, atoms%lmax(itype) ! loop: right basis function
                     DO m2 = -l2, l2
                        DO p2 = 1, 2
                           j2 = j2 + 1
                           rr = 0
                           DO ilharm = 0, nlharm ! loop: lattice harmonics of vxc
                              l = sphhar%llh(ilharm, typsym)
                              DO i = 1, sphhar%nmem(ilharm, typsym)
                                 M = sphhar%mlh(i, ilharm, typsym)
                                 rc = sphhar%clnu(i, ilharm, typsym)*gaunt1(l1, l, l2, m1, M, m2, atoms%lmaxd)
                                 rr = rr + integ(ilharm, p1, l1, p2, l2)*rc
                              END DO
                           END DO
                           rc = CMPLX(0, 1)**(l2 - l1) ! adjusts to a/b/ccof-scaling
                           vrmat%data_c(j1, j2) = rr*rc
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO

         call timestop("calc MT contrib")
         nnbas = j1

         ! Project on bascof
         call timestart("Project on bascof")
         DO ineq = 1, atoms%neq(itype)
            iatom = iatom + 1
            call timestart("zgemm bascof")
            !call zgemm(transa, transb, m,    n,      k,     alpha,   a,                 lda,             
            !           b,                 ldb,           beta,     c,    ldc)
            call zgemm("N", "T", nnbas, lapw%nv(jsp), nnbas, cmplx_1, vrmat%data_c(1,1), hybdat%maxlmindx, &
                        bascof(iatom)%data_c, lapw%dim_nvd(), cmplx_0, carr, hybdat%maxlmindx)
            carr = conjg(carr)

            !call zgemm(transa, transb, m,      n,           k,     alpha,   a,                 lda,            
            !              b,    ldb,              beta,    c,    ldc)
            call zgemm("N", "N", lapw%nv(jsp), lapw%nv(jsp), nnbas, cmplx_1, bascof(iatom)%data_c, lapw%dim_nvd(),&
                        carr, hybdat%maxlmindx, cmplx_0, carr1, lapw%dim_nvd()  )
            call timestop("zgemm bascof")

            call timestart("apply bascof to vxc")
            if(vxc%l_real) then
               do j = fmpi%n_rank+1,lapw%nv(jsp),fmpi%n_size
                  j_loc=(j-1)/fmpi%n_size+1
                  do i = 1,MIN(j,lapw%nv(jsp))  
                     vxc%data_r(i,j_loc) = vxc%data_r(i,j_loc) + real(carr1(i, j))
                  END DO
               END DO
            else
               do j = fmpi%n_rank+1,lapw%nv(jsp),fmpi%n_size
                  j_loc=(j-1)/fmpi%n_size+1
                  do i = 1,MIN(j,lapw%nv(jsp))
                     vxc%data_c(i,j_loc) = vxc%data_c(i,j_loc) + carr1(i, j)
                  END DO
               END DO
            endif
            call timestop("apply bascof to vxc")
         END DO
         call timestop("Project on bascof")
      END DO ! End loop over atom types

      ! Calculate plane wave contribution
      DO i = 1, stars%ng3
         vpw(i) = vx%pw_w(i, jsp)
         ! vpw(i)=vpwtot(i)-vpwcou(i,jsp)
      END DO

      ! Calculate vxc-matrix,  left basis function (ig1)
      !                        right basis function (ig2)
      call timestart("Calculate vxc-matrix")
      ! DO ig1 = 1, lapw%nv(jsp)
      !    DO ig2 = 1, ig1

      do ig1 = fmpi%n_rank+1,lapw%nv(jsp),fmpi%n_size
         ig1_loc = (ig1-1)/fmpi%n_size+1
         do ig2 = 1,min(ig1, lapw%nv(jsp))
            gg = lapw%gvec(:,ig1,jsp) - lapw%gvec(:,ig2,jsp)
            istar = stars%ig(gg(1), gg(2), gg(3))
            IF (istar /= 0) THEN
               if(vxc%l_real) then 
                  vxc%data_r(ig2,ig1_loc) = vxc%data_r(ig2,ig1_loc) + real(stars%rgphs(gg(1), gg(2), gg(3))*vpw(istar))
               else
                  vxc%data_c(ig2,ig1_loc) = vxc%data_c(ig2,ig1_loc) + stars%rgphs(gg(1), gg(2), gg(3))*vpw(istar)
               endif
            ELSE
               IF (fmpi%irank == 0) THEN
                  WRITE (oUnit, '(A,/6I5)') 'Warning: Gi-Gj not in any star:', &
                     lapw%gvec(:,ig1, jsp), lapw%gvec(:,ig2, jsp)
               END IF
            END IF
         END DO
      END DO
      call timestop("Calculate vxc-matrix")

      ! Calculate local orbital contribution
      call timestart("calculate LO contrib")
      IF (ANY(atoms%nlo /= 0)) THEN

         nbasf0 = lapw%nv(jsp)*(lapw%nv(jsp) + 1)/2  ! number of pure APW contributions
         icentry = nbasf0                          ! icentry counts the entry in the matrix vxc
         iatom = 0
         ikvecat = 0
         ikvecprevat = 0

         DO itype = 1, atoms%ntype

            typsym = sym%ntypsy(atoms%firstAtom(itype))
            nlharm = sphhar%nlh(typsym)

            ! Calculate vxc = vtot - vcoul
            call timestart("calc vxc=vtotvcoul")
            DO l = 0, nlharm
               DO i = 1, atoms%jri(itype)
                  IF (l == 0) THEN
                     ! vr(i,0)= vrtot(i,0,itype)*sfp/rmsh(i,itype) -  vrcou(i,0,itype,jsp)
                     vr(i, 0) = vx%mt(i, 0, itype, jsp) !*sfp_const/atoms%rmsh(i, itype)
                  ELSE ! vxc = vtot - vcoul
                     vr(i, l) = vx%mt(i, l, itype, jsp)                    !
                     ! vr(i,l)=  vrtot(i,l,itype)-vrcou(i,l,itype,jsp)
                  END IF
               END DO
            END DO
            call timestop("calc vxc=vtotvcoul")

            ! Precompute auxiliary radial integrals
            call timestart("Precompute aux. rad. integ.")
            DO ilharm = 0, nlharm
               i = 0
               DO l1 = 0, atoms%lmax(itype)
                  DO p1 = 1, mpdata%num_radfun_per_l(l1, itype)
                     i = i + 1
                     j = 0
                     DO l2 = 0, atoms%lmax(itype)
                        DO p2 = 1, mpdata%num_radfun_per_l(l2, itype)
                           j = j + 1
                           IF (j <= i) THEN
                              DO igrid = 1, atoms%jri(itype)
                                 grid(igrid) = vr(igrid, ilharm)*(bas1(igrid, p1, l1, itype)*bas1(igrid, p2, l2, itype) + &
                                                                  bas2(igrid, p1, l1, itype)*bas2(igrid, p2, l2, itype))
                              END DO
                              CALL intgr3(grid, atoms%rmsh(:, itype), atoms%dx(itype), atoms%jri(itype), integ(ilharm, p1, l1, p2, l2))
                              integ(ilharm, p2, l2, p1, l1) = integ(ilharm, p1, l1, p2, l2)
                           END IF
                        END DO
                     END DO
                  END DO
               END DO
            END DO
            call timestop("Precompute aux. rad. integ.")

            DO ieq = 1, atoms%neq(itype)
               iatom = iatom + 1
               IF ((sym%invsat(iatom) == 0) .OR. (sym%invsat(iatom) == 1)) THEN

                  IF (sym%invsat(iatom) == 0) invsfct = 1
                  IF (sym%invsat(iatom) == 1) invsfct = 2

                  DO ilo = 1, atoms%nlo(itype)
#ifdef CPP_OLDINTEL
                     CALL judft_error("no LOs   hybinp with old intel compiler!", calledby="subvxc.F90")
#else
                     l1 = atoms%llo(ilo, itype)
                     DO ikvec = 1, invsfct*(2*l1 + 1)
                        DO m1 = -l1, l1
                           DO p1 = 1, 3
                              IF (p1 == 3) THEN
                                 pp1 = pointer_lo(ilo, itype)
                              ELSE
                                 pp1 = p1
                              END IF

                              IF (mpdata%num_radfun_per_l(l1, itype) <= 2) call judft_error('subvxc: error mpdata%num_radfun_per_l')

                              lm = 0

                              !loop over APW
                              call timestart("loop over APW")
                              DO l2 = 0, atoms%lmax(itype)
                                 DO m2 = -l2, l2
                                    DO p2 = 1, 2
                                       lm = lm + 1

                                       rr = 0
                                       DO ilharm = 0, nlharm
                                          lh = sphhar%llh(ilharm, typsym)
                                          DO i = 1, sphhar%nmem(ilharm, typsym)
                                             mh = sphhar%mlh(i, ilharm, typsym)
                                             rc = sphhar%clnu(i, ilharm, typsym)*gaunt1(l1, lh, l2, m1, mh, m2, atoms%lmaxd)
                                             rr = rr + integ(ilharm, p2, l2, pp1, l1)*rc
                                          END DO
                                       END DO

                                       rc = CMPLX(0.0, 1.0)**(l2 - l1) ! adjusts to a/b/ccof-scaling
                                       ! ic counts the entry in vxc
                                       ic = icentry
                                       DO i = 1, lapw%nv(jsp)
                                          ic = ic + 1
                                          call packed_to_cart(ic, x,y)
                                          if(mod(x-(fmpi%n_rank+1),fmpi%n_size) == 0) then
                                             x_loc = (x-1) / fmpi%n_size+1
                                             IF (hmat%l_real) THEN
                                                vxc%data_r(y,x_loc) = vxc%data_r(y,x_loc) + invsfct*REAL(rr*rc*bascof(iatom)%data_c(i, lm)* &
                                                                                          CONJG(bascof_lo(p1, m1, ikvec, ilo, iatom)))
                                             ELSE
                                                vxc%data_c(y,x_loc) = vxc%data_c(y,x_loc) + rr*rc*bascof(iatom)%data_c(i, lm)* &
                                                                                          CONJG(bascof_lo(p1, m1, ikvec, ilo, iatom))
                                             END IF
                                          endif
                                       END DO
                                    END DO  !p2
                                 END DO  ! m2
                              END DO ! l2 ->  loop over APW
                              call timestop("loop over APW")

                              ! calcualte matrix-elements with local orbitals at the same atom
                              call timestart("calc. matelem with LO at same atm")
                              IF (ic /= icentry + lapw%nv(jsp)) call judft_error('subvxc: error counting ic')

                              ic = ic + ikvecprevat
                              DO ilop = 1, ilo - 1
                                 lp = atoms%llo(ilop, itype)
                                 DO ikvecp = 1, invsfct*(2*lp + 1)
                                    ic = ic + 1
                                    call packed_to_cart(ic, x,y)
                                    if(mod(x-(fmpi%n_rank+1),fmpi%n_size) == 0) then
                                       x_loc = (x-1) / fmpi%n_size+1
                                       DO mp = -lp, lp
                                          DO pp = 1, 3
                                             IF (pp == 3) THEN
                                                pp2 = pointer_lo(ilop, itype)
                                             ELSE
                                                pp2 = pp
                                             END IF

                                             rr = 0
                                             DO ilharm = 0, nlharm
                                                lh = sphhar%llh(ilharm, typsym)
                                                DO i = 1, sphhar%nmem(ilharm, typsym)
                                                   mh = sphhar%mlh(i, ilharm, typsym)
                                                   rc = sphhar%clnu(i, ilharm, typsym)*gaunt1(l1, lh, lp, m1, mh, mp, atoms%lmaxd)
                                                   rr = rr + integ(ilharm, pp2, lp, pp1, l1)*rc
                                                END DO
                                             END DO

                                             rc = CMPLX(0.0, 1.0)**(lp - l1) ! adjusts to a/b/ccof-scaling

                                             IF (hmat%l_real) THEN
                                                vxc%data_r(y,x_loc) = vxc%data_r(y,x_loc) &
                                                   + invsfct*REAL(rr*rc*bascof_lo(pp, mp, ikvecp, ilop, iatom)* &
                                                   CONJG(bascof_lo(p1, m1, ikvec, ilo, iatom)))
                                             ELSE
                                                vxc%data_c(y,x_loc) = vxc%data_c(y,x_loc)&
                                                   + rr*rc*bascof_lo(pp, mp, ikvecp, ilop, iatom)* &
                                                   CONJG(bascof_lo(p1, m1, ikvec, ilo, iatom))
                                             END IF
                                          END DO ! pp
                                       END DO ! mp
                                    endif
                                 END DO !ikvecp
                              END DO ! ilop
                              call timestop("calc. matelem with LO at same atm")

                              ! calculate matrix-elements of one local orbital with itself
                              call timestart("calc. matelem of LO with itself")
                              DO ikvecp = 1, ikvec
                                 ic = ic + 1
                                 call packed_to_cart(ic, x,y)
                                 if(mod(x-(fmpi%n_rank+1),fmpi%n_size) == 0) then
                                    x_loc = (x-1) / fmpi%n_size+1
                                    lp = l1
                                    ilop = ilo
                                    DO mp = -lp, lp
                                       DO pp = 1, 3
                                          IF (pp == 3) THEN
                                             pp2 = pointer_lo(ilop, itype)
                                          ELSE
                                             pp2 = pp
                                          END IF

                                          rr = 0
                                          DO ilharm = 0, nlharm
                                             lh = sphhar%llh(ilharm, typsym)
                                             DO i = 1, sphhar%nmem(ilharm, typsym)
                                                mh = sphhar%mlh(i, ilharm, typsym)
                                                rc = sphhar%clnu(i, ilharm, typsym)*gaunt1(l1, lh, lp, m1, mh, mp, atoms%lmaxd)
                                                rr = rr + integ(ilharm, pp2, lp, pp1, l1)*rc
                                             END DO
                                          END DO

                                          rc = CMPLX(0.0, 1.0)**(lp - l1) ! adjusts to a/b/ccof-scaling

                                          IF (hmat%l_real) THEN
                                             vxc%data_r(y,x_loc) = vxc%data_r(y,x_loc) &
                                                                     + invsfct*REAL(rr*rc*bascof_lo(pp, mp, ikvecp, ilop, iatom)* &
                                                                        CONJG(bascof_lo(p1, m1, ikvec, ilo, iatom)))
                                          ELSE
                                             vxc%data_c(y,x_loc) = vxc%data_c(y,x_loc)&
                                                                  + rr*rc*bascof_lo(pp, mp, ikvecp, ilop, iatom)* &
                                                                  CONJG(bascof_lo(p1, m1, ikvec, ilo, iatom))
                                          END IF
                                       END DO ! pp
                                    END DO ! mp
                                 endif
                              END DO ! ikvecp
                              call timestop("calc. matelem of LO with itself")
                           END DO  ! p1
                        END DO  ! m1
                        icentry = ic
                     END DO !ikvec
                     ikvecat = ikvecat + invsfct*(2*l1 + 1)
#endif
                  END DO  ! ilo
                  ikvecprevat = ikvecprevat + ikvecat
                  ikvecat = 0
               END IF ! sym%invsat(iatom)
            END DO ! ieq
         END DO !itype
      END IF ! if any atoms%llo
      call timestop("calculate LO contrib")

      !initialize weighting factor
      a_ex = xcpot%get_exchange_weight()

      call timestart("apply to hmat")
      do n = fmpi%n_rank+1,hmat%matsize1,fmpi%n_size
         n_loc = (n-1) / fmpi%n_size + 1
         DO nn = 1,MIN(n,hmat%matsize1)  
            IF (hmat%l_real) THEN
               hmat%data_r(nn, n_loc) = hmat%data_r(nn, n_loc) - a_ex * REAL(vxc%data_r(nn,n_loc))
            ELSE
               hmat%data_c(nn, n_loc) = hmat%data_c(nn, n_loc) - a_ex *      vxc%data_c(nn,n_loc)
            ENDIF
         END DO
      END DO

      call timestop("apply to hmat")

      CALL timestop("subvxc")

      deallocate(bascof)
   
   contains
      subroutine packed_to_cart(p_idx, x, y)
         implicit none
         integer, intent(in)    :: p_idx 
         integer, intent(inout) :: x, y 

         x = ceiling(-0.5 + sqrt(0.25 + 2*p_idx)  )
         y = p_idx - ((x - 1)*x)/2
      end subroutine packed_to_cart
   END SUBROUTINE subvxc
END MODULE m_subvxc
