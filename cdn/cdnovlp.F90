      MODULE m_cdnovlp
      USE m_juDFT
      CONTAINS
        SUBROUTINE cdnovlp(mpi,&
             &                   sphhar,stars,atoms,sym,&
             &                   DIMENSION,vacuum,cell,&
             &                   input,oneD,&
             &                   jspin,rh,&
             &                   qpw,rhtxy,rho,rht)
          !*****************************************************************
          !     calculates the overlapping core tail density and adds
          !     its contribution to the corresponging components of
          !     valence density.      
          !
          !     OLD VERSION:
          !     The previous version to calculate the overlap of the
          !     core tails was done in real space:
          !     A three dimensional real space grid was set up. On each
          !     of these grid points the charge density of all atoms inside
          !     the unit cell and neigboring unit cells was calculated.
          !     This calculation required a lagrange fit since the core
          !     charge is on a radial mesh. The same was done in the vacuum
          !     region, except on each z-plane we worked on a two dimension
          !     grid. After the charge density was generated on a equidistant
          !     grid the charge density was FFTed into G-space. The set up
          !     of the charge density on the real space grid is rather time
          !     consuming. 3-loops are required for the 3D grid
          !                1-loop over all atoms in the unit cells
          !                Larange interpolation
          !                3D and 2D FFTs
          !     In order to save time the non-spherical contributions inside
          !     the sphere had been ignored. It turns out that the later
          !     approximation is pure in the context of force calculations.
          !     
          !     PRESENT VERSION:
          !     The present version is written from scratch. It is based on the
          !     idea that the FFT of an overlap of spherically symmetric
          !     charges can be expressed by the product of
          !
          !     sum_natype{ F(G,ntype) * sum_atom(atype) {S(\vec{G},atom)}}
          !
          !     of form factor F and structure factor S. The Form factor
          !     depends only G while the structure factor depends on \vec{G}
          !     and can build up in G-space. F of a gaussian chargedensity can
          !     be calculated analytically.

          !     The core-tails to the vacuum region are described by an
          !     exponentially decaying function into the vacuum:

          !     rho(r_||,z,ivac)= sum_n{ rho(n,ivac) * exp(-kappa*(z-z_v))
          !                                          * exp(iG(n)r_||) }

          !     And the plane waves are expanded into lattice harmonics
          !     up to a l_cutoff. Tests of the accuracy inside the sphere
          !     have shown that a reduction of the l_cutoff inside the 
          !     in order to save time leads to sizable errors and should 
          !     be omitted.

          !     rho_L(r) =  4 \pi i^l \sum_{g =|= 0}  \rho_int(g) r_i^{2} \times
          !                              j_{l} (gr_i) \exp{ig\xi_i} Y^*_{lm} (g)

          !     Tests have shown that the present version is about 10 times
          !     faster than the previous one also all nonspherical terms are
          !     included up to l=8 and the previous one included only l=0.
          !     For l=16 it is still faster by factor 2.

          !     coded                  Stefan Bl"ugel, IFF Nov. 1997
          !     tested                 RObert Abt    , IFF Dez. 1997
          !*****************************************************************
          !
          USE m_intgr, ONLY : intgr3,intgz0
          USE m_constants
          USE m_spgrot
          USE m_qpwtonmt
          USE m_cylbes
          USE m_dcylbs
          USE m_rcerf
          USE m_od_cylbes
          USE m_od_chirot
          USE m_diflgr
          USE m_types
#ifdef CPP_MPI
          USE m_mpi_bc_st
#endif
          IMPLICIT NONE
          !
          !     .. Parameters ..
          TYPE(t_mpi),INTENT(IN)     :: mpi
          TYPE(t_sphhar),INTENT(IN)   :: sphhar
          TYPE(t_atoms),INTENT(IN)    :: atoms
          TYPE(t_stars),INTENT(IN)    :: stars
          TYPE(t_cell),INTENT(IN)     :: cell
          TYPE(t_sym),INTENT(IN)      :: sym
          TYPE(t_oneD),INTENT(IN)     :: oneD
          TYPE(t_dimension),INTENT(IN)::DIMENSION
          TYPE(t_vacuum),INTENT(in):: vacuum
          TYPE(t_input),INTENT(in)::input

          INTEGER    method1, method2
          PARAMETER (method1 = 2, method2 = 1)
          !     ..
          !     .. Scalar Arguments ..

          INTEGER,INTENT (IN) :: jspin       
          !     ..
          !     .. Array Arguments ..
          COMPLEX,INTENT (INOUT) :: qpw(stars%n3d,DIMENSION%jspd)
          COMPLEX,INTENT (INOUT) :: rhtxy(vacuum%nmzxyd,oneD%odi%n2d-1,2,DIMENSION%jspd)
          REAL,   INTENT (INOUT) :: rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd,DIMENSION%jspd)
          REAL,   INTENT (INOUT) :: rht(vacuum%nmzd,2,DIMENSION%jspd)
          REAL,   INTENT (INOUT) :: rh(DIMENSION%msh,atoms%ntypd)
          !     ..
          !     .. Local Scalars ..
          COMPLEX czero,sf,carg,VALUE,slope,ci
          REAL    ai,ar,a4,dif,dxx,f11,f12,g,gr,gz,qfin,qfout,dtildh,&
               &        rkappa,sign,signz,tol_14,x,z,zero,zvac,alpha3,&
               &        g2,phi,gamma,qq
          INTEGER ig3,imz,ir,ivac,j,j1,k,kz,k1,k2,l_cutoff,m0,&
               &        n,nz,nrz,nat1,nat2,nzvac,n_out_p,nat,&
               &        irec2,irec3,irec1,m,gzi
          LOGICAL tail
          !     ..
          !     .. Local Arrays ..
          COMPLEX, ALLOCATABLE :: qpwc(:)
          REAL,    ALLOCATABLE :: qf(:)
          REAL    acoff(atoms%ntypd),alpha(atoms%ntypd),rho_out(2)
          COMPLEX phas(sym%nop)
          REAL    rhohelp(DIMENSION%msh),rat(DIMENSION%msh,atoms%ntypd)
          INTEGER kr(3,sym%nop),mshc(atoms%ntypd)
          REAL    kro(3,oneD%ods%nop),fJ(-oneD%odi%M:oneD%odi%M),dfJ(-oneD%odi%M:oneD%odi%M)
          COMPLEX phaso(oneD%ods%nop)
          !     ..
          DATA  czero /(0.0,0.0)/, zero /0.0/, tol_14 /1.0e-10/!-14
          !
          !----> Abbreviation
          !
          !     l_cutoff : cuts off the l-expansion of the nonspherical charge
          !                density due to the core-tails of neigboring atoms
          !     mshc     : maximal radial meshpoint for which the radial coretail
          !                density is larger than tol_14
          !     method1  : two different ways to calculate the derivative of the
          !                charge density at the sphere boundary.
          !                (1) use subroutine diflgr based on lagrange interpol.
          !                (2) use two point formular in real space, 
          !                    see notes of SB.
          !                Tests have shown that (2) is more accurate.
          !     method2  : two different integration routines to calculate form
          !                factor of coretails outside the sphere.
          !                (1) use subroutine intgrz to integrate the tails from
          !                    outside to inside.
          !                (2) use subroutine intgr3 to integrate the tails from
          !                    muffin-tin radius to outside and include correction
          !                    for start up.
          !                Tests have shown that (1) is more accurate.
          !           
          !
          ci = CMPLX(0.0,1.0)

          ALLOCATE (qpwc(stars%n3d),qf(stars%n3d))

          IF (mpi%irank ==0) THEN
             !
             !----> prepare local array to store pw-expansion of pseudo core charge
             !
             DO k = 1 , stars%ng3    
                qpwc(k) = czero
             ENDDO
             ! 
             !----> (1) set up radial mesh beyond muffin-tin radius
             !      (2) cut_off core tails from noise 
             !
             nloop: DO  n = 1 , atoms%ntype
                IF (atoms%ncst(n).GT.0) THEN
                   DO  j = 1 , atoms%jri(n)
                      rat(j,n) = atoms%rmsh(j,n)
                   ENDDO
                   dxx = EXP(atoms%dx(n))
                   DO j = atoms%jri(n) + 1 , DIMENSION%msh
                      rat(j,n) = rat(j-1,n)*dxx
                   ENDDO
                   DO j = atoms%jri(n) - 1 , DIMENSION%msh
                      rh(j,n) = rh(j,n)/ (fpi_const*rat(j,n)*rat(j,n))
                   ENDDO
                   DO j = DIMENSION%msh , atoms%jri(n) , -1
                      IF ( rh(j,n) .GT. tol_14 ) THEN
                         mshc(n) = j
                         CYCLE nloop
                      END IF
                   ENDDO
                   mshc(n) = atoms%jri(n)
                ENDIF
             ENDDO nloop
             !
             !-----> the core density inside the spheres is replaced by a
             !       gaussian-like pseudo density : n(r) = acoff*exp(-alpha*r*r)
             !       acoff and alpha determined to obtain a continous and 
             !       differentiable density at the sphere boundary.
             !       IF mshc = jri  either core tail too small or no core (i.e. H)
             !
             DO  n = 1,atoms%ntype
                IF ((mshc(n).GT.atoms%jri(n)).AND.(atoms%ncst(n).GT.0)) THEN

                   j1 = atoms%jri(n) - 1
                   IF ( method1 .EQ. 1) THEN
                      dif = diflgr(rat(j1,n),rh(j1,n))
                      WRITE (6,FMT=8000) n,rh(atoms%jri(n),n),dif
                      alpha(n) = -0.5 * dif / ( rh(atoms%jri(n),n)*atoms%rmt(n) )
                   ELSEIF ( method1 .EQ. 2) THEN
                      alpha(n) = LOG( rh(j1,n) / rh(atoms%jri(n),n) )
                      alpha(n) = alpha(n)&
                           &                   / ( atoms%rmt(n)*atoms%rmt(n)*( 1.0-EXP( -2.0*atoms%dx(n) ) ) )
                   ELSE
                      WRITE (6,'('' error in choice of method1 in cdnovlp '')')
                      CALL juDFT_error("error in choice of method1 in cdnovlp"&
                           &              ,calledby ="cdnovlp")
                   ENDIF
                   acoff(n) = rh(atoms%jri(n),n) * EXP( alpha(n)*atoms%rmt(n)*atoms%rmt(n) )
                   WRITE (6,FMT=8010) alpha(n),acoff(n)
                   DO j = 1,atoms%jri(n) - 1
                      rh(j,n) = acoff(n) * EXP( -alpha(n)*rat(j,n)**2 )
                   ENDDO

                ELSE
                   alpha(n) = 0.0
                ENDIF
             ENDDO
             !
8000         FORMAT (/,10x,'core density and its first derivative ',&
                  &                 'at sph. bound. for atom type',&
                  &             i2,' is',3x,2e15.7)
8010         FORMAT (/,10x,'alpha=',f10.5,5x,'acoff=',f10.5)
             !
             !=====> calculate the fourier transform of the core-pseudocharge
             !
             !     qpw(\vec{G}) = Sum_{n} [ F(G,n) * Sum_{atm{n}} S(\vec{G},atm) ]
             !                  n = atom_type
             !                  F = Formfactor = F_in_sphere + F_outsphere
             !                  S = Structure factor
             !
             tail = .FALSE.
             !
             !*****> start loop over the atom type
             !
             nat1 = 1
             DO  n = 1,atoms%ntype
                IF ( ( mshc(n) .GT. atoms%jri(n) ).AND.&
                     &        ( alpha(n) .GT. tol_14 ) )    THEN
                   !
                   n_out_p = mshc(n)-atoms%jri(n)+1
                   !
                   ! (1) Form factor for each atom type
                   !
                   f11 = tpi_const * atoms%rmt(n) * rh(atoms%jri(n),n) / alpha(n)
                   f12 = acoff(n) * ( pi_const/alpha(n) ) *SQRT(pi_const/alpha(n))
                   ar  = SQRT( alpha(n) ) * atoms%rmt(n) 
                   !
                   DO  k = 1,stars%ng3
                      g = stars%sk3(k)
                      !    first G=0
                      IF ( k.EQ.1 ) THEN
                         ai = zero
                         !
                         ! ---->     calculate form factor inside the mt-sphere
                         !           (use analytic integration of gaussian)
                         !
                         qfin = - f11 + f12 * rcerf(ar,ai)
                         !
                         ! ---->     calculate form factor outside the mt-sphere
                         !           (do numerical integration of tails)
                         !
                         IF ( method2 .EQ. 1) THEN

                            DO ir = 1 , n_out_p
                               j  = atoms%jri(n)+ir-1
                               rhohelp(mshc(n)+1-j) =  rat(j,n) * rat(j,n) &
                                    &                                     * rat(j,n) *  rh(j,n)
                            END DO
                            CALL intgz0(rhohelp,atoms%dx(n),n_out_p,qfout,tail)

                         ELSE

                            DO ir = 1 , n_out_p
                               j  = atoms%jri(n)+ir-1
                               rhohelp(ir) = rat(j,n) * rat(j,n) * rh(j,n)
                            END DO
                            CALL intgr3(rhohelp,rat(atoms%jri(n),n),atoms%dx(n),&
                                 &                        n_out_p,qfout)
                            !--->             have to remove the small r-correction from intgr3
                            qfout=qfout-atoms%rmt(n)*rhohelp(1)

                         END IF

                         qfout = fpi_const * qfout
                         !
                      ELSE 
                         !    then  G>0
                         ai = 0.5*g/SQRT(alpha(n))
                         gr = g*atoms%rmt(n)
                         a4 = 0.25/alpha(n)
                         !
                         ! ---->     calculate form factor inside the mt-sphere
                         !           (use analytic integration of gaussian)
                         !
                         qfin = - f11 * SIN(gr)/gr &
                              &                + f12 * rcerf(ar,ai) * EXP(-a4*g*g) 
                         !
                         ! ---->     calculate form factor outside the mt-sphere
                         !           (do numerical integration of tails)

                         IF ( method2 .EQ. 1) THEN

                            DO ir = 1 , n_out_p 
                               j  = atoms%jri(n)+ir-1
                               rhohelp(mshc(n)-atoms%jri(n)+2-ir) =  rat(j,n)*rat(j,n) &
                                    &                                     * rh(j,n) * SIN( g*rat(j,n) )
                            END DO
                            !
                            !--->       note we use here the integration routine for vacuum. Because 
                            !           the vacuum integration is made for an inwards integration 
                            !           from outside to inside. Outside the starting value will be 
                            !           nearly zero since the core density is small. if the tail 
                            !           correction (tail=.true.) is included for the integrals, the 
                            !           integrand is from infinity inward. This might not be 
                            !           necessary. Further the integration routine is made for 
                            !           equidistant meshpoints, therefore the term r(i) of
                            !           dr/di = dx*r(i) is included in rhohelp


                            CALL intgz0(rhohelp,atoms%dx(n),n_out_p,qfout,tail)

                         ELSE

                            DO ir = 1 , n_out_p
                               j  = atoms%jri(n)+ir-1
                               rhohelp(ir) = rat(j,n) * rh(j,n) * SIN(g*rat(j,n))
                            END DO
                            CALL intgr3(rhohelp,rat(atoms%jri(n),n),atoms%dx(n),&
                                 &                        n_out_p,qfout)
                            !--->             have to remove the small r-correction from intgr3
                            !roa...correction.from.intgr3.......................
                            IF (rhohelp(1)*rhohelp(2).GT.zero) THEN
                               alpha3 = 1.0 + LOG(rhohelp(2)/rhohelp(1))/atoms%dx(n)
                               IF (alpha3.GT.zero)&
                                    &                 qfout = qfout - rat(atoms%jri(n),n)*rhohelp(1)/alpha3
                            ENDIF
                            !roa...end.correction...............................


                         END IF

                         qfout = fpi_const * qfout / g
                         !
                      END IF
                      !
                      qf(k)    = (qfin + qfout)/cell%omtil
                   ENDDO
                   !
                   ! (2) structure constant for each atom of atom type
                   !
                   !
                   !    first G=0
                   !
                   k=1
                   qpw(k,jspin) = qpw(k,jspin) + atoms%neq(n) * qf(k)
                   qpwc(k)      = qpwc(k)      + atoms%neq(n) * qf(k)
                   !
                   !    then  G>0
                   !
                   DO  k = 2,stars%ng3
                      IF (.NOT.oneD%odi%d1) THEN
                         CALL spgrot(&
                              &                       sym%nop,sym%symor,sym%mrot,sym%tau,sym%invtab,&
                              &                       stars%kv3(:,k),&
                              &                       kr,phas)
                         !
                         ! ----> start loop over equivalent atoms
                         !
                         nat2 = nat1 + atoms%neq(n) - 1
                         DO  nat = nat1,nat2
                            sf = czero
                            DO  j = 1,sym%nop
                               x = -tpi_const* ( kr(1,j) * atoms%taual(1,nat)&
                                    &                           + kr(2,j) * atoms%taual(2,nat)&
                                    &                           + kr(3,j) * atoms%taual(3,nat))
                               sf = sf + CMPLX(COS(x),SIN(x))*phas(j)
                            ENDDO
                            sf = sf / REAL( sym%nop )
                            qpw(k,jspin) = qpw(k,jspin) + sf * qf(k)
                            qpwc(k)      = qpwc(k)      + sf * qf(k)
                         ENDDO
                      ELSE
                         !-odim
                         CALL od_chirot(oneD%odi,oneD%ods,cell%bmat,stars%kv3(:,k),kro,phaso)
                         nat2 = nat1 + atoms%neq(n) - 1
                         DO  nat = nat1,nat2
                            !                  sf = cmplx(1.,0.)
                            sf = czero
                            DO  j = 1,oneD%ods%nop
                               x = -tpi_const* ( kro(1,j)*atoms%taual(1,nat)&
                                    &                          + kro(2,j)*atoms%taual(2,nat)&
                                    &                          + kro(3,j)*atoms%taual(3,nat))
                               sf = sf + CMPLX(COS(x),SIN(x))*phaso(j)
                            ENDDO
                            sf = sf / REAL( oneD%ods%nop )
                            qpw(k,jspin) = qpw(k,jspin) + sf * qf(k)
                            qpwc(k)      = qpwc(k)      + sf * qf(k)
                         ENDDO
                         !+odim
                      END IF
                   ENDDO
                END IF
                nat1 = nat1 + atoms%neq(n)
             ENDDO
             !
             !=====> calculate core-tails to the vacuum region                
             !       Coretails expanded in exponentially decaying functions.
             !       Describe vacuum by: 
             !       rho(r_||,z,ivac)= sum_n{ rho(n,ivac) * exp(-kappa*(z-z_v)) 
             !                                           * exp(iG(n)r_||) }
             IF (input%film .AND. .NOT.oneD%odi%d1) THEN
                !+gu
                dtildh = 0.5 * tpi_const / cell%bmat(3,3)
                IF (vacuum%nvac.EQ.1) THEN
                   rho_out(1) = qpwc(1)*cell%z1
                   DO k = 2,stars%ng3
                      IF ((stars%kv3(1,k).EQ.0).AND.(stars%kv3(2,k).EQ.0)) THEN
                         nz = 1
                         IF (sym%invs.OR.sym%zrfs) nz = 2
                         g = stars%kv3(3,k) * cell%bmat(3,3)
                         rho_out(1) = rho_out(1) + nz*qpwc(k)*SIN(g*cell%z1)/g
                      ENDIF
                   ENDDO
                   rho_out(1) =  qpwc(1) * dtildh - rho_out(1)
                ELSE
                   DO ivac = 1, vacuum%nvac
                      carg = CMPLX(0.0,0.0)
                      DO k = 2,stars%ng3
                         IF ((stars%kv3(1,k).EQ.0).AND.(stars%kv3(2,k).EQ.0)) THEN
                            g = stars%kv3(3,k) * cell%bmat(3,3) * (3. - 2.*ivac)
                            carg = carg -qpwc(k)*(EXP(ci*g*dtildh)-EXP(ci*g*cell%z1))/g
                         ENDIF
                      ENDDO
                      rho_out(ivac) = qpwc(1) * ( dtildh-cell%z1 ) - AIMAG(carg)
                   ENDDO
                ENDIF
                !-gu
                !        nzvac = min(50,nmz)
                m0 = -stars%mx3
                IF (sym%zrfs) m0 = 0
                !
                !---> loop over 2D stars
                !
                DO 280 k = 1,stars%ng2
                   k1 = stars%kv2(1,k)
                   k2 = stars%kv2(2,k)
                   DO 270 ivac = 1,vacuum%nvac
                      VALUE = czero
                      slope = czero
                      sign = 3. - 2.*ivac
                      !
                      ! ---> sum over gz-stars
                      DO 250 kz = m0,stars%mx3
                         ig3 = stars%ig(k1,k2,kz)
                         !        ----> use only stars within the g_max sphere (oct.97 shz)
                         IF (ig3.NE.0) THEN
                            nz = 1
                            IF (sym%zrfs) nz = stars%nstr(ig3)/stars%nstr2(k)
                            gz = kz*cell%bmat(3,3)
                            DO 240 nrz = 1,nz
                               signz = 3. - 2.*nrz
                               carg = ci*sign*signz*gz
                               VALUE = VALUE + qpwc(ig3)* EXP(carg*cell%z1)
                               slope = slope + carg*qpwc(ig3)* EXP(carg*cell%z1)
240                         ENDDO
                         END IF
250                   ENDDO
                      ! roa work-around
                      IF (  ABS(REAL(VALUE)).GT.zero ) THEN
                         ! roa work-around
                         ! gb works also around
                         rkappa = - REAL( slope/VALUE )
                         IF (k.EQ.1) rkappa = VALUE/rho_out(ivac)
                         !               rkappa = - sign * real( slope/value )
                         IF (rkappa.LE.zero) rkappa=MIN(rkappa,-tol_14)
                         IF (rkappa.GT.zero) rkappa=MAX(rkappa,tol_14)
                         ! gb works also around
                         zvac   = - LOG( tol_14/cabs(VALUE) ) / rkappa
                         nzvac  = INT( zvac/vacuum%delz ) + 1
                         !               IF ( rkappa.GT.zero .AND. real(value).GT.zero ) THEN
                         IF ( rkappa.GT.zero ) THEN
                            z = 0. 
                            IF ( k.EQ.1 ) THEN
                               DO imz = 1 , MIN( nzvac,vacuum%nmz )
                                  rht(imz,ivac,jspin) = &
                                       &                  rht(imz,ivac,jspin) + VALUE*EXP(-rkappa*z)
                                  z = z + vacuum%delz
220                            ENDDO
                            ELSE
                               DO imz = 1 , MIN( nzvac,vacuum%nmzxy )
                                  rhtxy(imz,k-1,ivac,jspin) = &
                                       &                  rhtxy(imz,k-1,ivac,jspin) + VALUE*EXP(-rkappa*z)
                                  z = z + vacuum%delz
230                            ENDDO
                            END IF
                         END IF
                         ! roa work-around
                      END IF
                      ! roa work-around
270                ENDDO
280             ENDDO
             ELSEIF (oneD%odi%d1) THEN
                !-odim
                !--->  rho_out is the charge lost between the interstitial and the
                !--->  rectangular unit cell

                rho_out(1) = cell%vol*qpwc(1)

                DO k = 2,stars%ng3
                   IF (stars%kv3(3,k).EQ.0) THEN
                      irec3 = stars%ig(stars%kv3(1,k),stars%kv3(2,k),stars%kv3(3,k))
                      IF (irec3.NE.0) THEN
                         irec2 = stars%ig2(irec3)
                         IF (irec2.NE.1) THEN
                            g2 = stars%sk2(irec2)
                            CALL cylbes(oneD%odi%M,g2*cell%z1,fJ)
                            rho_out(1) = rho_out(1) +&
                                 &                    qpwc(k)*2.*cell%vol*fJ(1)/(cell%z1*g2)
                         END IF
                      END IF
                   END IF
                END DO

                rho_out(1) = cell%bmat(3,3)*(qpwc(1)*cell%omtil - rho_out(1))/(tpi_const*tpi_const)
                !     then we are constructing our radial components of the vacuum
                !     charge density so that the so that they are continuous in the
                !     value and slope on the boundary
                !     value = sum{gpar}[qpw(gz,gpar)exp{-i*m*phi(gpar))J_m(gpar*z1)]
                !     slope = .......the same...................*gpar*dJ_m(gpar*z1)]
                !     for determining the rht we need only continuity, because the
                !     second condition is the charge neutrality of the system
                !                         Y.Mokrousov
                DO irec1 = 1,oneD%odi%nq2
                   VALUE = czero
                   slope = czero
                   gzi = oneD%odi%kv(1,irec1)
                   m = oneD%odi%kv(2,irec1)
                   DO irec2 = 1,stars%ng2
                      irec3 = stars%ig(stars%kv2(1,irec2),stars%kv2(2,irec2),gzi)
                      IF (irec3.NE.0) THEN
                         g2 = stars%sk2(irec2)
                         phi = stars%phi2(irec2)
                         CALL cylbes(oneD%odi%M,g2*cell%z1,fJ)
                         CALL dcylbs(oneD%odi%M,g2*cell%z1,fJ,dfJ)
                         VALUE = VALUE + (ci**m)*qpwc(irec3)*&
                              &                 EXP(CMPLX(0.,-m*phi))*fJ(m)
                         slope = slope + (ci**m)*g2*qpwc(irec3)*&
                              &                 EXP(CMPLX(0.,-m*phi))*dfJ(m)
                      END IF
                   END DO

                   IF (ABS(REAL(VALUE)).GT.zero) THEN
                      IF (irec1.EQ.1) THEN
                         qq = REAL(VALUE/rho_out(1))
                         rkappa = (cell%z1*qq+SQRT(cell%z1*cell%z1*qq*qq + 4*qq))/2.
                         gamma = rkappa
                      ELSE
                         rkappa = gamma
                         !                  rkappa = -real(slope/value)
                      END IF
                      IF (rkappa.LE.zero) rkappa=MIN(rkappa,-tol_14)
                      IF (rkappa.GT.zero) rkappa=MAX(rkappa,tol_14)
                      zvac=-LOG(tol_14/cabs(VALUE))/rkappa
                      nzvac=INT(zvac/vacuum%delz)+1
                      IF (rkappa.GT.zero .AND. REAL(VALUE).GT.zero) THEN
                         z = 0.0
                         IF (irec1.EQ.1) THEN
                            DO imz = 1,MIN(nzvac,vacuum%nmz)
                               rht(imz,vacuum%nvac,jspin)=rht(imz,vacuum%nvac,jspin)+&
                                    &                       VALUE*EXP(-rkappa*z)
                               z = z + vacuum%delz
                            END DO
                         ELSE
                            DO imz = 1,MIN(nzvac,vacuum%nmzxy)
                               rhtxy(imz,irec1-1,vacuum%nvac,jspin)=&
                                    &                          rhtxy(imz,irec1-1,vacuum%nvac,jspin)+&
                                    &                          VALUE*EXP(-rkappa*z)
                               z = z + vacuum%delz
                            END DO
                         END IF
                      END IF
                   END IF
                END DO
                !+odim

             END IF
             !
             !=====> update density inside the spheres                        
             !
             ! ----> (1) subtract on-site contribution, because 
             !           they are contained in the plane wave part 
             !
             DO n = 1,atoms%ntype
                IF  ((mshc(n).GT.atoms%jri(n)).AND.(atoms%ncst(n).GT.0)) THEN
                   DO j = 1,atoms%jri(n)
                      rho(j,0,n,jspin) = rho(j,0,n,jspin)&
                           &                          - sfp_const*rat(j,n)*rat(j,n)*rh(j,n)
                   ENDDO
                ENDIF
             ENDDO
             !
             ! ----> (2) add the plane wave contribution of (core tails + on-site 
             !           contribution) to the m.t. density, include full nonspherical 
             !           components
             !
             l_cutoff = 0
             DO n = 1,atoms%ntype
                l_cutoff = MAX( l_cutoff, atoms%lmax(n) )
             END DO
          ENDIF ! mpi%irank ==0

#ifdef CPP_MPI
          IF ( mpi%isize > 1 ) THEN
             CALL mpi_bc_st(&
                  &               mpi%mpi_comm,mpi%irank,sphhar%memd,sphhar%nlhd,sphhar%ntypsd,atoms%jmtd,atoms%ntypd,stars%n3d,&
                  &               jspin,l_cutoff,stars%ng3,atoms%ntype,sym%nop,atoms%natd,sym%symor,&
                  &               sphhar%clnu,qpwc,atoms%lmax,atoms%ntypsy,atoms%jri,sphhar%nmem,sphhar%nlh,sphhar%mlh,stars%nstr,&
                  &               atoms%neq,stars%kv3,sym%mrot,sym%invtab,sphhar%llh,cell%bmat,sym%tau,atoms%taual,atoms%dx,&
                  &               atoms%rmsh,stars%sk3)
          ENDIF
#endif

          CALL qpw_to_nmt(&
               &                sphhar,atoms,stars,&
               &                sym,cell,oneD,mpi,&
               &                jspin,l_cutoff,qpwc,&
               &                rho)

#ifdef CPP_MPI
          IF ( mpi%isize > 1 ) THEN
             CALL mpi_col_st(                        ! Collect rho()&
             &                 mpi,atoms,sphhar,&
                  &                 rho(1,0,1,jspin))
          ENDIF
#endif

          DEALLOCATE (qpwc,qf)

        END SUBROUTINE cdnovlp
      END MODULE m_cdnovlp
