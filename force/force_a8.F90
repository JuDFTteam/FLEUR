MODULE m_forcea8
CONTAINS
   SUBROUTINE force_a8(input,atoms,sym,sphhar,jsp,vr,rho,force,fmpi,results)
      !--------------------------------------------------------------------------
      ! Pulay 1st term force contribution à la Rici et al. 
      ! 
      ! Equation A8, Phys. Rev. B 43, 6411
      !--------------------------------------------------------------------------
      USE m_intgr, ONLY : intgr3
      USE m_constants
      USE m_gaunt, ONLY :gaunt1
      USE m_differentiate,ONLY: difcub
      USE m_types
      USE m_juDFT

      IMPLICIT NONE

      TYPE(t_input),   INTENT(IN)    :: input
      TYPE(t_atoms),   INTENT(IN)    :: atoms
      TYPE(t_sym),     INTENT(IN)    :: sym
      TYPE(t_sphhar),  INTENT(IN)    :: sphhar
      TYPE(t_force),   INTENT(IN)    :: force
      TYPE(t_mpi),     INTENT(IN)    :: fmpi
      TYPE(t_results), INTENT(INOUT) :: results

      INTEGER, INTENT(IN) :: jsp 

      REAL, INTENT(IN)    :: vr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype)
      REAL, INTENT(IN)    :: rho(:,0:,:,:)!(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins)

      ! Local scalars
      COMPLEX, PARAMETER :: czero=CMPLX(0.0,0.0)
      COMPLEX aaa, bbb, ccc, ddd, eee, fff
      REAL    a8_1, a8_2, qval, rr, truerho, truev, xi
      INTEGER i, ir, j, l, l1, l2, lh1, lh2, m1, m2, mem1, mem2, n, nd, na, m

      ! Local arrays
      COMPLEX forc_a8(3), gv(3), f_sum(3)
      REAL    rhoaux(atoms%jmtd), rhodif(atoms%jmtd)

      ! Statement functions
      REAL    alpha, beta, delta, epslon, gamma, phi 
      INTEGER krondel

      ! Kronecker delta for arguments >=0 AND <0
      krondel(i,j) = MIN(ABS(i)+1,ABS(j)+1)/MAX(ABS(i)+1,ABS(j)+1)* (1+SIGN(1,i)*SIGN(1,j))/2
      alpha(l,m) = (l+1)*0.5e0*SQRT(REAL((l-m)* (l-m-1))/ REAL((2*l-1)* (2*l+1)))
      beta(l,m) = l*0.5e0*SQRT(REAL((l+m+2)* (l+m+1))/ REAL((2*l+1)* (2*l+3)))
      GAMMA(l,m) = (l+1)*0.5e0*SQRT(REAL((l+m)* (l+m-1))/ REAL((2*l-1)* (2*l+1)))
      delta(l,m) = l*0.5e0*SQRT(REAL((l-m+2)* (l-m+1))/ REAL((2*l+1)* (2*l+3)))
      epslon(l,m) = (l+1)*SQRT(REAL((l-m)* (l+m))/ REAL((2*l-1)* (2*l+1)))
      phi(l,m) = l*SQRT(REAL((l-m+1)* (l+m+1))/REAL((2*l+1)* (2*l+3)))

      CALL timestart("force_a8")

      WRITE  (oUnit,*)
 
      DO n = 1, atoms%ntype
         na = atoms%firstAtom(n)
         IF (atoms%l_geo(n)) THEN
            nd = sym%ntypsy(na)

            DO i = 1,3
               forc_a8(i) = czero
            END DO

            ! TODO: There is no output for this. Do we want some?

            CALL intgr3(rho(:,0,n,jsp),atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),qval)

            ! Check if the l=0 density is correct.
            ! Note that in general also all l>0 components of the density have
            ! been multiplied by r**2 (see for example checkdop which constructs
            ! the true MT-density at rmt).
            ! The factor sqrt(4*pi) comes from Y_00*\int d\Omega = 1/sqrt(4*pi)*4*pi
8000        FORMAT (' FORCE_A8: valence charge=',1p,e16.8)

            ! PART I of FORCE_A8
            DO lh1 = 0, sphhar%nlh(nd)
               l1 = sphhar%llh(lh1,nd)
               DO lh2 = 0, sphhar%nlh(nd)
                  l2 = sphhar%llh(lh2,nd)

                  DO i = 1, 3
                     gv(i) = czero
                  END DO

                  ! Sum over all m for a particular lattice harmonic.
                  DO mem1 = 1, sphhar%nmem(lh1,nd)
                     m1 = sphhar%mlh(mem1,lh1,nd)
                     DO mem2 = 1, sphhar%nmem(lh2,nd)
                        m2 = sphhar%mlh(mem2,lh2,nd)
                        gv(1) = gv(1) + SQRT(2.e0*pi_const/3.e0) * &
                              sphhar%clnu(mem1,lh1,nd)*sphhar%clnu(mem2,lh2,nd)*&
                              (gaunt1(1,l1,l2,-1,m1,m2,atoms%lmaxd) - &
                               gaunt1(1,l1,l2,1,m1,m2,atoms%lmaxd))
                        gv(2) = gv(2) - ImagUnit*SQRT(2.e0*pi_const/3.e0) * &
                              sphhar%clnu(mem1,lh1,nd)*sphhar%clnu(mem2,lh2,nd)*&
                              (gaunt1(1,l1,l2,-1,m1,m2,atoms%lmaxd) + &
                               gaunt1(1,l1,l2,1,m1,m2,atoms%lmaxd))
                        gv(3) = gv(3) + SQRT(4.e0*pi_const/3.e0) * &
                              sphhar%clnu(mem1,lh1,nd)*sphhar%clnu(mem2,lh2,nd)*&
                              gaunt1(1,l1,l2,0,m1,m2,atoms%lmaxd)
                     END DO
                  END DO

                  ! Note that in general also all l>0 components of the density
                  ! have been multiplied by r**2.
                  ! Here we need the true radial denisity for performing the
                  ! derivative. Rherefore we divide by r**2.
                  DO ir = 1, atoms%jri(n)
                     rhoaux(ir) = rho(ir,lh2,n,jsp)/ (atoms%rmsh(ir,n)**2)
                  END DO

                  ! NOTE: Here we should have: vr  = vtrue
                  ! difcub performs the analytic derivative of Lagrangian of 3rd order
                  xi = atoms%rmsh(1,n)
                  rr = xi*xi
                  rhodif(1) = difcub(atoms%rmsh(1,n),rhoaux(1),xi)*vr(1,lh1,n)*rr
                  DO ir = 2, atoms%jri(n) - 2
                     xi = atoms%rmsh(ir,n)
                     rr = xi*xi
                     rhodif(ir) = difcub(atoms%rmsh(ir-1,n),rhoaux(ir-1),xi)*vr(ir,lh1,n)*rr
                  END DO

                  xi = atoms%rmsh(atoms%jri(n)-1,n)
                  rr = xi*xi
                  rhodif(atoms%jri(n)-1) = difcub(atoms%rmsh(atoms%jri(n)-3,n), &
                                                   rhoaux(atoms%jri(n)-3),xi) * &
                                                   vr(atoms%jri(n)-1,lh1,n)*rr

                  xi = atoms%rmsh(atoms%jri(n),n)
                  rr = xi*xi
                  rhodif(atoms%jri(n)) = difcub(atoms%rmsh(atoms%jri(n)-3,n), &
                                                 rhoaux(atoms%jri(n)-3),xi) * &
                                                 vr(atoms%jri(n),lh1,n)*rr

                  ! NOTE: vr(l=0) is EXPLICITELY multiplied by r/sqrt(4pi) to be
                  ! the TRUE r*V which is needed for the radial Schrödinger eq.
                  ! Here, we need the l=0 component, v_{l=0}(r) which will be
                  ! multiplied by Y_00 in the lm expansion; therefore we MUST
                  ! recorrect vr(l=0) by the inverse factor sqrt(4pi)/r. We do
                  ! the correction for the product array rhodif.
                  IF (lh1.EQ.0) THEN
                     DO ir = 1,atoms%jri(n)
                        rhodif(ir) = rhodif(ir)/atoms%rmsh(ir,n)*sfp_const
                     END DO
                  END IF

                  CALL intgr3(rhodif,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),a8_1)

                  DO i = 1,3
                     forc_a8(i) = forc_a8(i) + a8_1*gv(i)
                  END DO

               END DO !lh2 (0:sphhar%nlh(nd))
            END DO ! lh1 (0:sphhar%nlh(nd))

            ! PART II of FORCE_A8
            DO lh1 = 0, sphhar%nlh(nd)
               l1 = sphhar%llh(lh1,nd)
               DO lh2 = 0, sphhar%nlh(nd)
                  l2 = sphhar%llh(lh2,nd)
                  DO ir = 1, atoms%jri(n)
                     truev = vr(ir,lh1,n)
                     truerho = rho(ir,lh2,n,jsp)/ (atoms%rmsh(ir,n)**2)
                     rhoaux(ir) = truev*truerho*atoms%rmsh(ir,n)
                  END DO

                  IF (lh1.EQ.0) THEN
                     DO ir = 1,atoms%jri(n)
                        rhoaux(ir) = rhoaux(ir)/atoms%rmsh(ir,n)*sfp_const
                     END DO
                  END IF

                  CALL intgr3(rhoaux,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),a8_2)

                  DO i = 1, 3
                     gv(i) = (0.0,0.0)
                  END DO

                  ! Sum over all m for a particular lattice harmonic.
                  DO mem1 = 1, sphhar%nmem(lh1,nd)
                     m1 = sphhar%mlh(mem1,lh1,nd)
                     DO mem2 = 1, sphhar%nmem(lh2,nd)
                        m2 = sphhar%mlh(mem2,lh2,nd)

                        ! NOTE: delta(-m,m')(-1)^m was applied because we have the integrand
                        ! Y nabla Y instead of Y* nabla Y.
                        ! We also know that  Y*(lm) = (-1)^m Y(l,-m)
                        !                and  Y(lm) = (-1)^m Y*(l,-m).
                        ! Therefore  (-1)^m delta(-m,m') appears.

                        aaa = alpha(l2,m2)*sphhar%clnu(mem1,lh1,nd) * &
                                           sphhar%clnu(mem2,lh2,nd) * &
                                           (-1)**m1*krondel(l1,l2-1)* &
                                                    krondel(-m1,m2+1)
                        bbb = beta(l2,m2)*sphhar%clnu(mem1,lh1,nd) *  &
                                          sphhar%clnu(mem2,lh2,nd) *  &
                                          (-1)**m1*krondel(l1,l2+1)*  &
                                                   krondel(-m1,m2+1)
                        ccc = GAMMA(l2,m2)*sphhar%clnu(mem1,lh1,nd) * &
                                           sphhar%clnu(mem2,lh2,nd) * &
                                           (-1)**m1*krondel(l1,l2-1)* &
                                                    krondel(-m1,m2-1)
                        ddd = delta(l2,m2)*sphhar%clnu(mem1,lh1,nd) * &
                                           sphhar%clnu(mem2,lh2,nd) * &
                                           (-1)**m1*krondel(l1,l2+1)* &
                                                    krondel(-m1,m2-1)
                        eee = epslon(l2,m2)*sphhar%clnu(mem1,lh1,nd)* &
                                            sphhar%clnu(mem2,lh2,nd)* &
                                            (-1)**m1*krondel(l1,l2-1)*&
                                                     krondel(-m1,m2)
                        fff = phi(l2,m2)*sphhar%clnu(mem1,lh1,nd) *   &
                                         sphhar%clnu(mem2,lh2,nd) *   &
                                         (-1)**m1*krondel(l1,l2+1)*   &
                                                  krondel(-m1,m2)

                        gv(1) = gv(1) + aaa + bbb - ccc - ddd
                        gv(2) = gv(2) - ImagUnit* (aaa+bbb+ccc+ddd)
                        gv(3) = gv(3) + eee - fff
                     END DO
                  END DO

                  DO i = 1, 3
                     forc_a8(i) = forc_a8(i) + a8_2*gv(i)
                  END DO

               END DO ! lh2 (0:sphhar%nlh(nd))
            END DO ! lh1 (0:sphhar%nlh(nd))

            ! Add onto existing forces.

            DO i = 1, 3
               results%force(i,n,jsp) = results%force(i,n,jsp) + REAL(forc_a8(i))
            END DO

            ! Write out result.
            WRITE (oUnit,FMT=8010) n
            WRITE (oUnit,FMT=8020) (forc_a8(i),i=1,3)
8010        FORMAT (' FORCES: EQUATION A8 FOR ATOM TYPE',i4)
8020        FORMAT (' FX_A8=',2f10.6,' FY_A8=',2f10.6,' FZ_A8=',2f10.6)

         END IF
      END DO

      ! Write out the result of a12, a21, b4 and b8
      ! here as well.

      IF (.NOT.input%l_useapw) THEN

         WRITE  (oUnit,*)
        
         IF (fmpi%isize.EQ.1) THEN
            DO n=1, atoms%ntype
               IF (atoms%l_geo(n)) THEN
                  WRITE  (oUnit,FMT=8030) n
                  WRITE  (oUnit,FMT=8040) (force%f_a12(i,n),i=1,3)
               END IF
8030           FORMAT (' FORCES: EQUATION A12 FOR ATOM TYPE',i4)
8040           FORMAT (' FX_A12=',2f10.6,' FY_A12=',2f10.6,' FZ_A12=',2f10.6)
            END DO
         ELSE
            WRITE (oUnit,*) "If this was a serial calculation, the A12 force component would be written out here. In parallel it holds no meaning."
         END IF
      ELSE

         WRITE  (oUnit,*)

         DO n=1, atoms%ntype
            IF (atoms%l_geo(n)) THEN
               WRITE  (oUnit,FMT=8070) n
               WRITE  (oUnit,FMT=8080) (force%f_b4(i,n),i=1,3)
            END IF
8070        FORMAT (' FORCES: EQUATION B4 FOR ATOM TYPE',i4)
8080        FORMAT (' FX_B4=',2f10.6,' FY_B4=',2f10.6,' FZ_B4=',2f10.6)
         END DO

         WRITE  (oUnit,*)

         DO n=1,atoms%ntype
            IF (atoms%l_geo(n)) THEN
               WRITE  (oUnit,FMT=8090) n
               WRITE  (oUnit,FMT=8100) (force%f_b8(i,n),i=1,3)
            END IF
8090        FORMAT (' FORCES: EQUATION B8 FOR ATOM TYPE',i4)
8100        FORMAT (' FX_B8=',2f10.6,' FY_B8=',2f10.6,' FZ_B8=',2f10.6)
         END DO
      END IF

      WRITE  (oUnit,*)

      IF (fmpi%isize.EQ.1) THEN
         DO n=1,atoms%ntype
            IF (atoms%l_geo(n)) THEN
               WRITE  (oUnit,FMT=8050) n
               WRITE  (oUnit,FMT=8060) (force%f_a21(i,n),i=1,3)
            END IF
8050        FORMAT (' FORCES: EQUATION A21 FOR ATOM TYPE',i4)
8060        FORMAT (' FX_A21=',2f10.6,' FY_A21=',2f10.6,' FZ_A21=',2f10.6)
         END DO
      ELSE
         WRITE (oUnit,*) "If this was a serial calculation, the A21 force component would be written out here. In parallel it holds no meaning."
      END IF


      CALL timestop("force_a8")

   END SUBROUTINE force_a8
END MODULE m_forcea8
