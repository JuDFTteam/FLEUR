!--------------------------------------------------------------------------------
! Copyright (c) 2020 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_forcea21
CONTAINS
   SUBROUTINE force_a21(input,atoms,sym,oneD,cell,we,jsp,epar,ne,eig,usdus,tlmplm,&
                        vtot,eigVecCoeffs,aveccof,bveccof,cveccof,f_a21,f_b4,results)
      !--------------------------------------------------------------------------
      ! Pulay 2nd and 3rd term force contributions à la Rici et al.
      !
      ! Equation A17 and A20 combined, Phys. Rev. B 43, 6411
      !
      ! Note1: We do NOT include the i**l factors in the alm, blm coming from
      ! to_pulay anymore. Therefore, we can use matrix elements from file 28, 38
      ! DIRECTLY.
      !
      ! Note2: The present version only yields forces for the highest energy window
      ! (=valence states). If semicore forces are wanted as well the tmas and tmat
      ! files have to be saved, indexed and properly used here in force_a21.
      !
      ! 22/june/97: Probably found symmetrization error replacing S^-1 by S
      ! (IS instead of isinv)
      !
      ! Force contribution B4 added following
      ! Madsen, Blaha, Schwarz, Sjostedt, Nordstrom
      ! GMadsen FZJ 20/3-01
      !--------------------------------------------------------------------------

      USE m_forcea21lo
      USE m_forcea21U
      USE m_types_setup
      USE m_types_misc
      USE m_types_usdus
      USE m_types_tlmplm
      USE m_types_cdnval
      USE m_types_potden
      USE m_constants
      USE m_juDFT

      IMPLICIT NONE

      TYPE(t_input),        INTENT(IN)    :: input
      TYPE(t_atoms),        INTENT(IN)    :: atoms
      TYPE(t_sym),          INTENT(IN)    :: sym
      TYPE(t_oneD),         INTENT(IN)    :: oneD
      TYPE(t_cell),         INTENT(IN)    :: cell
      TYPE(t_usdus),        INTENT(IN)    :: usdus
      TYPE(t_tlmplm),       INTENT(IN)    :: tlmplm
      TYPE(t_potden),       INTENT(IN)    :: vtot
      TYPE(t_eigVecCoeffs), INTENT(IN)    :: eigVecCoeffs
      TYPE(t_results),      INTENT(INOUT) :: results

      INTEGER, INTENT(IN) :: jsp, ne

      REAL,    INTENT(IN)    :: we(ne), epar(0:atoms%lmaxd,atoms%ntype)
      REAL,    INTENT(IN)    :: eig(input%neig)
      COMPLEX, INTENT(IN)    :: aveccof(3,ne,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat)
      COMPLEX, INTENT(IN)    :: bveccof(3,ne,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat)
      COMPLEX, INTENT(IN)    :: cveccof(3,-atoms%llod:atoms%llod,ne,atoms%nlod,atoms%nat)
      COMPLEX, INTENT(INOUT) :: f_a21(3,atoms%ntype)
      COMPLEX, INTENT(INOUT) :: f_b4(3,atoms%ntype)

      REAL,    PARAMETER :: zero=0.0
      COMPLEX, PARAMETER :: czero=CMPLX(0.,0.)
      COMPLEX dtd, dtu, utd, utu
      INTEGER lo
      INTEGER i, ie, im, l1, l2, ll1, ll2, lm1, lm2, m1, m2, n, natom, m
      INTEGER natrun, is, isinv, j, irinv, it, lmplmd

      REAL, ALLOCATABLE :: a21(:,:), b4(:,:)
      COMPLEX forc_a21(3), forc_b4(3)
      REAL starsum(3), starsum2(3), gvint(3), gvint2(3)
      REAL vec(3), vec2(3), vecsum(3), vecsum2(3)

      CALL timestart("force_a21")

      lmplmd = (atoms%lmaxd*(atoms%lmaxd+2)* (atoms%lmaxd*(atoms%lmaxd+2)+3))/2

      ALLOCATE(a21(3,atoms%nat),b4(3,atoms%nat) )

      natom = 1
      DO  n = 1,atoms%ntype
         IF (atoms%l_geo(n)) THEN
            forc_a21(:) = czero
            forc_b4(:) = czero

            DO natrun = natom,natom + atoms%neq(n) - 1
               a21(:,natrun) = zero
               b4(:,natrun) = zero
            END DO

            DO ie = 1,ne
               DO l1 = 0,atoms%lmax(n)
                  ll1 = l1* (l1+1)
                  DO m1 = -l1,l1
                     lm1 = ll1 + m1
                     DO l2 = 0,atoms%lmax(n)
                        ll2 = l2* (l2+1)
                        DO m2 = -l2,l2
                           lm2 = ll2 + m2
                           DO natrun = natom,natom + atoms%neq(n) - 1
                              utu = CONJG(tlmplm%h_loc(lm2,lm1,n,jsp,jsp))
                              dtd = CONJG(tlmplm%h_loc(lm2+tlmplm%h_loc2(n),lm1+tlmplm%h_loc2(n),n,jsp,jsp))
                              utd = CONJG(tlmplm%h_loc(lm2+tlmplm%h_loc2(n),lm1,n,jsp,jsp))
                              dtu = CONJG(tlmplm%h_loc(lm2,lm1+tlmplm%h_loc2(n),n,jsp,jsp))
                              DO i = 1,3
                                 a21(i,natrun) = a21(i,natrun) + 2.0*&
                                    AIMAG( CONJG(eigVecCoeffs%abcof(ie,lm1,0,natrun,jsp)) *utu*aveccof(i,ie,lm2,natrun)&
                                    +CONJG(eigVecCoeffs%abcof(ie,lm1,0,natrun,jsp)) *utd*bveccof(i,ie,lm2,natrun)&
                                    +CONJG(eigVecCoeffs%abcof(ie,lm1,1,natrun,jsp)) *dtu*aveccof(i,ie,lm2,natrun)&
                                    +CONJG(eigVecCoeffs%abcof(ie,lm1,1,natrun,jsp)) *dtd*bveccof(i,ie,lm2,natrun))*we(ie)/atoms%neq(n)
                              END DO ! i (spatial directions)
                           END DO ! natrun
                        END DO ! m2
                     END DO ! l2

                     ! Correct spherical part
                     utu = -eig(ie)
                     utd = 0.0
                     dtu = 0.0
                     dtd = utu*usdus%ddn(l1,n,jsp)
                     DO i = 1,3
                        DO natrun = natom,natom + atoms%neq(n) - 1
                           a21(i,natrun) = a21(i,natrun) + 2.0*AIMAG(&
                               CONJG(eigVecCoeffs%abcof(ie,lm1,0,natrun,jsp))*utu*aveccof(i,ie,lm1,natrun)&
                              +CONJG(eigVecCoeffs%abcof(ie,lm1,0,natrun,jsp))*utd*bveccof(i,ie,lm1,natrun)&
                              +CONJG(eigVecCoeffs%abcof(ie,lm1,1,natrun,jsp))*dtu*aveccof(i,ie,lm1,natrun)&
                              +CONJG(eigVecCoeffs%abcof(ie,lm1,1,natrun,jsp))*dtd*bveccof(i,ie,lm1,natrun)&
                              )*we(ie) /atoms%neq(n)
                        END DO
                     END DO
                  END DO ! m1
               END DO ! l1
            END DO ! ie

            ! Add the local orbital and U contribution to a21:

            CALL force_a21_lo(atoms,jsp,n,we,eig,ne,eigVecCoeffs,aveccof,bveccof,cveccof,tlmplm,usdus,a21)

            IF (atoms%n_u+atoms%n_hia>0) THEN
               CALL force_a21_U(atoms,n,jsp,we,ne,usdus,vTot%mmpMat(:,:,:,jsp),eigVecCoeffs,aveccof,bveccof,cveccof,a21)
            END IF

            IF (input%l_useapw) THEN
               ! B4 force
               DO ie = 1,ne
                  DO l1 = 0,atoms%lmax(n)
                     ll1 = l1* (l1+1)
                     DO m1 = -l1,l1
                        lm1 = ll1 + m1
                        DO i = 1,3
                           DO natrun = natom,natom + atoms%neq(n) - 1
                              b4(i,natrun) = b4(i,natrun) + 0.5 *&
                                 we(ie)/atoms%neq(n)*atoms%rmt(n)**2*AIMAG(&
                                 CONJG(eigVecCoeffs%abcof(ie,lm1,0,natrun,jsp)*usdus%us(l1,n,jsp)&
                                 +eigVecCoeffs%abcof(ie,lm1,1,natrun,jsp)*usdus%uds(l1,n,jsp))*&
                                 (aveccof(i,ie,lm1,natrun)*usdus%dus(l1,n,jsp)&
                                 +bveccof(i,ie,lm1,natrun)*usdus%duds(l1,n,jsp) )&
                                 -CONJG(aveccof(i,ie,lm1,natrun)*usdus%us(l1,n,jsp)&
                                 +bveccof(i,ie,lm1,natrun)*usdus%uds(l1,n,jsp) )*&
                                 (eigVecCoeffs%abcof(ie,lm1,0,natrun,jsp)*usdus%dus(l1,n,jsp)&
                                 +eigVecCoeffs%abcof(ie,lm1,1,natrun,jsp)*usdus%duds(l1,n,jsp)) )
                           END DO
                        END DO
                     END DO
                  END DO

                  DO lo = 1,atoms%nlo(n)
                     l1 = atoms%llo(lo,n)
                     DO m = -l1,l1
                        lm1 = l1* (l1+1) + m
                        DO i=1,3
                           DO natrun = natom,natom + atoms%neq(n) - 1
                              b4(i,natrun) = b4(i,natrun) + 0.5 *&
                                 we(ie)/atoms%neq(n)*atoms%rmt(n)**2*AIMAG(&
                                 CONJG( eigVecCoeffs%abcof(ie,lm1,0,natrun,jsp)* usdus%us(l1,n,jsp)&
                                 + eigVecCoeffs%abcof(ie,lm1,1,natrun,jsp)* usdus%uds(l1,n,jsp) ) *&
                                 cveccof(i,m,ie,lo,natrun)*usdus%dulos(lo,n,jsp)&
                                 + CONJG(eigVecCoeffs%ccof(m,ie,lo,natrun,jsp)*usdus%ulos(lo,n,jsp)) *&
                                 ( aveccof(i,ie,lm1,natrun)* usdus%dus(l1,n,jsp)&
                                 + bveccof(i,ie,lm1,natrun)* usdus%duds(l1,n,jsp)&
                                 + cveccof(i,m,ie,lo,natrun)*usdus%dulos(lo,n,jsp) )  &
                                 - (CONJG( aveccof(i,ie,lm1,natrun) *usdus%us(l1,n,jsp)&
                                 + bveccof(i,ie,lm1,natrun) *usdus%uds(l1,n,jsp) ) *&
                                 eigVecCoeffs%ccof(m,ie,lo,natrun,jsp)  *usdus%dulos(lo,n,jsp)&
                                 + CONJG(cveccof(i,m,ie,lo,natrun)*usdus%ulos(lo,n,jsp)) *&
                                 ( eigVecCoeffs%abcof(ie,lm1,0,natrun,jsp)*usdus%dus(l1,n,jsp)&
                                 + eigVecCoeffs%abcof(ie,lm1,1,natrun,jsp)*usdus%duds(l1,n,jsp)&
                                 + eigVecCoeffs%ccof(m,ie,lo,natrun,jsp)*usdus%dulos(lo,n,jsp) ) ) )
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END IF

            DO natrun = natom,natom + atoms%neq(n) - 1
               !  to complete summation over stars of k now sum
               !  over all operations which leave (k+G)*R(natrun)*taual(natrun)
               !  invariant. We sum over ALL these operations and not only
               !  the ones needed for the actual star of k. Should be
               !  ok if we divide properly by the number of operations
               !  First, we find operation S where RS=T. T -like R- leaves
               !  the above scalar product invariant (if S=1 then R=T).
               !  R is the operation which generates position of equivalent atom
               !  out of position of representative
               !  S=R^(-1) T
               !  number of ops which leave (k+G)*op*taual invariant: invarind
               !  index of inverse operation of R: irinv
               !  index of operation T: invarop
               !  now, we calculate index of operation S: is
               !
               !  note, that vector in expression A17,A20 + A21 is a
               !  reciprocal lattice vector! other transformation rules
               !
               !  transform recip vector g-g' into internal coordinates

               vec(:) = a21(:,natrun)
               vec2(:) = b4(:,natrun)

               gvint=MATMUL(cell%bmat,vec)/tpi_const
               gvint2=MATMUL(cell%bmat,vec2)/tpi_const
               vecsum(:) = zero
               vecsum2(:) = zero

               !-gb2002
               !            irinv = invtab(ngopr(natrun))
               !            DO it = 1,invarind(natrun)
               !               is = multab(irinv,invarop(natrun,it))
               !c  note, actually we need the inverse of S but -in principle
               !c  because {S} is a group and we sum over all S- S should also
               !c  work; to be lucid we take the inverse:
               !                isinv = invtab(is)
               !!               isinv = is
               ! Rotation is alreadt done in to_pulay, here we work only in the
               ! coordinate system of the representative atom (natom):

               DO it = 1,sym%invarind(natom)
                  is =sym%invarop(natom,it)
                  isinv = sym%invtab(is)
                  IF (oneD%odi%d1) isinv = oneD%ods%ngopr(natom)
                     !-gb 2002
                     !  now we have the wanted index of operation with which we have
                     !  to rotate gv. Note gv is given in cart. coordinates but
                     !  mrot acts on internal ones
                     DO i = 1,3
                        vec(i) = zero
                        vec2(i) = zero
                        DO j = 1,3
                           IF (.NOT.oneD%odi%d1) THEN
                              vec(i) = vec(i) + sym%mrot(i,j,isinv)*gvint(j)
                              vec2(i) = vec2(i) + sym%mrot(i,j,isinv)*gvint2(j)
                           ELSE
                              vec(i) = vec(i) + oneD%ods%mrot(i,j,isinv)*gvint(j)
                              vec2(i) = vec2(i) + oneD%ods%mrot(i,j,isinv)*gvint2(j)
                           END IF
                     END DO
                  END DO
                  DO i = 1,3
                     vecsum(i) = vecsum(i) + vec(i)
                     vecsum2(i) = vecsum2(i) + vec2(i)
                  END DO
               END DO ! end operator loop

               ! Transform from internal to cart. coordinates
               starsum=MATMUL(cell%amat,vecsum)
               starsum2=MATMUL(cell%amat,vecsum2)
               DO i = 1,3
                  forc_a21(i) = forc_a21(i) + starsum(i)/sym%invarind(natrun)
                  forc_b4(i) = forc_b4(i) + starsum2(i)/sym%invarind(natrun)
               END DO
            END DO ! natrun

            ! Add onto existing forces.

            ! NOTE: results%force is real and therefore only the real part of
            ! forc_a21 is added. In general, force must be real after the k-star
            ! summation. Now, we put the proper operations into real space.
            ! Problem: What happens if in real space there is no inversion anymore?
            ! But we have inversion in k-space due to time reversal symmetry:
            ! E(k)=E(-k)
            ! We argue that k-space inversion is automatically taken into account
            ! if force = (1/2)(forc_a21+conjg(forc_a21)), because time reversal
            ! symmetry means that conjg(PSI) is also a solution of Schrödinger eq.
            ! if PSI is one.

            DO i = 1, 3
               results%force(i,n,jsp) = results%force(i,n,jsp) + REAL(forc_a21(i) + forc_b4(i))
               f_a21(i,n)     = f_a21(i,n)     + forc_a21(i)
               f_b4(i,n)      = f_b4(i,n)      + forc_b4(i)
            END DO

         END IF ! l_geo(n)
         natom = natom + atoms%neq(n)
      END DO

      ! The result is written in force_a8.

      CALL timestop("force_a21")

   END SUBROUTINE force_a21
END MODULE m_forcea21
