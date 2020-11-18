MODULE m_forcea12
CONTAINS
   SUBROUTINE force_a12(atoms,nobd,sym,cell,oneD,we,jsp,ne,usdus,eigVecCoeffs, &
                        acoflo,bcoflo,e1cof,e2cof,f_a12,results)
      !--------------------------------------------------------------------------
      ! Pulay 1st term force contribution à la Rici et al. 
      ! 
      ! Equation A12, Phys. Rev. B 43, 6411
      !--------------------------------------------------------------------------
      USE m_types_setup
      USE m_types_misc
      USE m_types_usdus
      USE m_types_cdnval
      USE m_constants
      USE m_juDFT
    
      IMPLICIT NONE

      TYPE(t_atoms),        INTENT(IN)    :: atoms
      TYPE(t_sym),          INTENT(IN)    :: sym
      TYPE(t_cell),         INTENT(IN)    :: cell
      TYPE(t_oneD),         INTENT(IN)    :: oneD
      TYPE(t_usdus),        INTENT(IN)    :: usdus
      TYPE(t_eigVecCoeffs), INTENT(IN)    :: eigVecCoeffs
      TYPE(t_results),      INTENT(INOUT) :: results

      INTEGER, INTENT(IN) :: nobd    
      INTEGER, INTENT(IN) :: jsp, ne 

      REAL,    INTENT(IN)    :: we(nobd)
      COMPLEX, INTENT(IN)    :: acoflo(-atoms%llod:atoms%llod,ne,atoms%nlod,atoms%nat)
      COMPLEX, INTENT(IN)    :: bcoflo(-atoms%llod:atoms%llod,ne,atoms%nlod,atoms%nat)
      COMPLEX, INTENT(IN)    :: e1cof(ne,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat)
      COMPLEX, INTENT(IN)    :: e2cof(ne,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat)
      COMPLEX, INTENT(INOUT) :: f_a12(3,atoms%ntype)

      ! Local scalars
      REAL, PARAMETER :: zero=0.0
      COMPLEX, PARAMETER :: czero=CMPLX(0.0,0.0)
      COMPLEX a12,cil1,cil2
      INTEGER i,ie,irinv,is,isinv,it,j,l,l1,l2,lm1,lm2 ,m1,m2,n,natom,natrun,ilo,m

      ! Local arrays
      COMPLEX forc_a12(3),gv(3)
      COMPLEX acof_flapw(nobd,0:atoms%lmaxd*(atoms%lmaxd+2)),bcof_flapw(nobd,0:atoms%lmaxd*(atoms%lmaxd+2))
      REAL aaa(2),bbb(2),ccc(2),ddd(2),eee(2),fff(2),gvint(3),starsum(3),vec(3),vecsum(3)

      ! Statement functions
      REAL    alpha,beta,delta,epslon,gamma,phi 
      INTEGER krondel

      ! Kronecker delta for arguments >=0 AND <0
      krondel(i,j) = MIN(ABS(i)+1,ABS(j)+1)/MAX(ABS(i)+1,ABS(j)+1)* (1+SIGN(1,i)*SIGN(1,j))/2
      alpha(l,m) = (l+1)*0.5e0*SQRT(REAL((l-m)* (l-m-1))/ REAL((2*l-1)* (2*l+1)))
      beta(l,m) = l*0.5e0*SQRT(REAL((l+m+2)* (l+m+1))/ REAL((2*l+1)* (2*l+3)))
      GAMMA(l,m) = (l+1)*0.5e0*SQRT(REAL((l+m)* (l+m-1))/ REAL((2*l-1)* (2*l+1)))
      delta(l,m) = l*0.5e0*SQRT(REAL((l-m+2)* (l-m+1))/ REAL((2*l+1)* (2*l+3)))
      epslon(l,m) = (l+1)*SQRT(REAL((l-m)* (l+m))/ REAL((2*l-1)* (2*l+1)))
      phi(l,m) = l*SQRT(REAL((l-m+1)* (l+m+1))/REAL((2*l+1)* (2*l+3)))

      CALL timestart("force_a12")

      natom = 1
      DO  n = 1, atoms%ntype
         IF (atoms%l_geo(n)) THEN
            forc_a12(:) = czero

            DO natrun = natom, natom + atoms%neq(n) - 1

               gv(:) = czero

               ! The local orbitals do not contribute to the term a12, because
               ! they vanish at the MT-boundary. Therefore, the LO-contribution
               ! to the a and b coefficients has to be subtracted before 
               ! calculating a12.

               DO l1 = 0, atoms%lmax(n)
                  DO m1 = -l1, l1
                     lm1 = l1*(l1+1) + m1
                     DO ie = 1, ne
                        acof_flapw(ie,lm1) = eigVecCoeffs%acof(ie,lm1,natrun,jsp)
                        bcof_flapw(ie,lm1) = eigVecCoeffs%bcof(ie,lm1,natrun,jsp)
                     END DO
                  END DO
               END DO

               DO ilo = 1, atoms%nlo(n)
                  l1 = atoms%llo(ilo,n)
                  DO m1 = -l1, l1
                     lm1 = l1*(l1+1) + m1
                     DO ie = 1, ne
                        acof_flapw(ie,lm1) = acof_flapw(ie,lm1) - acoflo(m1,ie,ilo,natrun)
                        bcof_flapw(ie,lm1) = bcof_flapw(ie,lm1) - bcoflo(m1,ie,ilo,natrun)
                     END DO
                  END DO
               END DO

               DO l1 = 0, atoms%lmax(n)
                  cil1 = ImagUnit**l1
                  DO m1 = -l1, l1
                     lm1 = l1*(l1+1) + m1
                     DO l2 = 0, atoms%lmax(n)
                        cil2 = ImagUnit**l2
                        DO m2 = -l2, l2
                           lm2 = l2*(l2+1) + m2

                           a12 = czero

                           DO ie = 1, ne

                              a12 = a12 + CONJG(cil1*&
                                 (acof_flapw(ie,lm1)*usdus%us(l1,n,jsp) + bcof_flapw(ie,lm1)*usdus%uds(l1,n,jsp) ))*cil2*&
                                 (e1cof(ie,lm2,natrun)*usdus%us(l2,n,jsp)+ e2cof(ie,lm2,natrun)*usdus%uds(l2,n,jsp))*we(ie)

                           END DO

                           aaa(1) = alpha(l1,m1)*krondel(l2,l1-1)* krondel(m2,m1+1)
                           aaa(2) = alpha(l2,m2)*krondel(l1,l2-1)* krondel(m1,m2+1)
                           bbb(1) = beta(l1,m1)*krondel(l2,l1+1)* krondel(m2,m1+1)
                           bbb(2) = beta(l2,m2)*krondel(l1,l2+1)* krondel(m1,m2+1)
                           ccc(1) = GAMMA(l1,m1)*krondel(l2,l1-1)* krondel(m2,m1-1)
                           ccc(2) = GAMMA(l2,m2)*krondel(l1,l2-1)* krondel(m1,m2-1)
                           ddd(1) = delta(l1,m1)*krondel(l2,l1+1)* krondel(m2,m1-1)
                           ddd(2) = delta(l2,m2)*krondel(l1,l2+1)* krondel(m1,m2-1)
                           eee(1) = epslon(l1,m1)*krondel(l2,l1-1)* krondel(m2,m1)
                           eee(2) = epslon(l2,m2)*krondel(l1,l2-1)* krondel(m1,m2)
                           fff(1) = phi(l1,m1)*krondel(l2,l1+1)* krondel(m2,m1)
                           fff(2) = phi(l2,m2)*krondel(l1,l2+1)* krondel(m1,m2)

                           gv(1) = gv(1) + (aaa(1)+bbb(1)-ccc(1)-ddd(1)+ &
                                            aaa(2)+bbb(2)-ccc(2)-ddd(2)) * &
                                            0.5* atoms%rmt(n)**2*a12

                           gv(2) = gv(2) + ImagUnit* (aaa(1)+bbb(1)+ccc(1)+ddd(1) - &
                                                      aaa(2)-bbb(2)-ccc(2)-ddd(2)) * &
                                                      0.5* atoms%rmt(n)**2*a12

                           gv(3) = gv(3) + (eee(1)+eee(2)-fff(1)-fff(2)) * &
                                            0.5*atoms%rmt(n)**2*a12

                        END DO ! m2 (-l2:l2)
                     END DO ! l2 (0:atoms%lmax(n))
                  END DO ! m1 (-l1:l1)
               END DO ! l1 (0:atoms%lmax(n))

               ! To complete the k summation over the stars now sum over all 
               ! operations which leave (k+G)*R(natrun)*taual(natrun) invariant.
               ! We sum over ALL these operations and not only the ones needed
               ! for the actual star of k. This should be ok if we divide properly
               ! by the number of operations.
               ! First, we find the operation S where RS=T. 
               ! T (like R) leaves the above scalar product invariant (if S=1 then
               ! R=T). R is the operation which generates the position of an equivalent
               ! atom from the position of the representative.

               ! S=R^(-1) T
               ! Number of ops which leave (k+G)*op*taual invariant: invarind
               ! Index of inverse operation of R: irinv
               ! Index of operation T: invarop
               ! Now, we calculate the index of operation S: is

               ! Transform vector gv into internal coordinates:

               vec(:) = REAL(gv(:)) /atoms%neq(n)

               gvint=MATMUL(cell%bmat,vec)/tpi_const
 
               vecsum(:) = zero

               !-gb2002
               !  irinv = invtab(ngopr(natrun))
               !  DO it = 1,invarind(natrun)
               !     is = multab(irinv,invarop(natrun,it))
               !     ! Note, actually we need the inverse of S but -in principle
               !c  because {S} is agroup and we sum over all S- S should also
               !c  work; to be lucid we take the inverse:
               !                isinv = invtab(is)
               !!               isinv = is
               ! Rotation is alreadt done in to_pulay, here we work only in the
               ! coordinate system of the representative atom (natom):
               !!
        
               DO it = 1, sym%invarind(natom)
                  is =sym%invarop(natom,it)
                  isinv = sym%invtab(is)
                  IF (oneD%odi%d1) isinv = oneD%ods%ngopr(natom)
                  !-gb 2002
                  ! Now we have the wanted index of operation with which we have
                  ! to rotate gv. Note gv is given in cart. coordinates but mrot
                  ! acts on internal ones.
                  DO i = 1, 3
                     vec(i) = zero
                     DO j = 1, 3
                        IF (.NOT.oneD%odi%d1) THEN
                           vec(i) = vec(i) + sym%mrot(i,j,isinv)*gvint(j)
                        ELSE
                           vec(i) = vec(i) + oneD%ods%mrot(i,j,isinv)*gvint(j)
                        END IF
                     END DO
                  END DO

                  DO i = 1, 3
                     vecsum(i) = vecsum(i) + vec(i)
                  END DO

               END DO ! End of operator loop

               ! Transform from internal to cart. coordinates:
               starsum=MATMUL(cell%amat,vecsum)
               DO i = 1, 3
                  forc_a12(i) = forc_a12(i) + starsum(i)/sym%invarind(natrun)
               END DO
             
            END DO ! End of natrun loop

            ! Add onto existing forces.

            ! NOTE: results%force is real and therefore only the real part of
            ! forc_a12 is added. In general, force must be real after the k-star
            ! summation. Now, we put the proper operations into real space. 
            ! Problem: What happens if in real space there is no inversion anymore?
            ! But we have inversion in k-space due to time reversal symmetry:
            ! E(k)=E(-k)
            ! We argue that k-space inversion is automatically taken into account
            ! if force = (1/2)(forc_a12+conjg(forc_a12)), because time reversal
            ! symmetry means that conjg(PSI) is also a solution of Schrödinger eq.
            ! if PSI is one.

            DO i = 1, 3
               results%force(i,n,jsp) = results%force(i,n,jsp) + REAL(forc_a12(i))
               f_a12(i,n) = f_a12(i,n) + forc_a12(i)
            END DO

         END IF
         natom = natom + atoms%neq(n)
      END DO

      ! The result is written in force_a8.

      CALL timestop("force_a12")

   END SUBROUTINE force_a12
END MODULE m_forcea12
