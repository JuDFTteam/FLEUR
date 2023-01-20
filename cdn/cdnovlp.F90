!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_cdnovlp
   USE m_juDFT
#ifdef CPP_MPI
   USE mpi
#endif
   USE m_force_a4_add
   USE m_sphbes
   USE m_phasy1

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: cdnovlp 
CONTAINS
   SUBROUTINE cdnovlp(fmpi, sphhar, stars, atoms, sym, vacuum,cell, input, &
                        l_st, jspin, rh, qpw, rhtxy, rho, rht, vpw, vr)
      !--------------------------------------------------------------------------
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
      ! 
      !     The core-tails to the vacuum region are described by an
      !     exponentially decaying function into the vacuum:
      ! 
      !     rho(r_||,z,ivac)= sum_n{ rho(n,ivac) * exp(-kappa*(z-z_v))
      !                                          * exp(iG(n)r_||) }
      ! 
      !     And the plane waves are expanded into lattice harmonics
      !     up to a l_cutoff. Tests of the accuracy inside the sphere
      !     have shown that a reduction of the l_cutoff inside the 
      !     in order to save time leads to sizable errors and should 
      !     be omitted.
      ! 
      !     rho_L(r) =  4 \pi i^l \sum_{g =|= 0}  \rho_int(g) r_i^{2} \times
      !                              j_{l} (gr_i) \exp{ig\xi_i} Y^*_{lm} (g)
      ! 
      !     Tests have shown that the present version is about 10 times
      !     faster than the previous one also all nonspherical terms are
      !     included up to l=8 and the previous one included only l=0.
      !     For l=16 it is still faster by factor 2.
      !
      !     coded                  Stefan Bl"ugel, IFF Nov. 1997
      !     tested                 RObert Abt    , IFF Dez. 1997
      !
      !     Added calculation of force contributions from coretails
      !     outside of their native muffin-tin spheres, i.e. in the
      !     interstitial region and other muffin-tins; only for bulk.
      !     Refer to Klüppelberg et al., PRB 91 035105 (2015)
      !     Aaron Klueppelberg, Oct. 2015
      ! 
      !--------------------------------------------------------------------------
         
          USE m_constants
          USE m_qpwtonmt
          USE m_cylbes
          USE m_dcylbs
           
          USE m_diflgr
          USE m_types
          USE m_intgr, ONLY : intgr3
#ifdef CPP_MPI
          USE m_mpi_bc_st
#endif
          !
          !     .. Parameters ..
          TYPE(t_mpi),INTENT(IN)     :: fmpi
          TYPE(t_sphhar),INTENT(IN)   :: sphhar
          TYPE(t_atoms),INTENT(IN)    :: atoms
          TYPE(t_stars),INTENT(IN)    :: stars
          TYPE(t_cell),INTENT(IN)     :: cell
          TYPE(t_sym),INTENT(IN)      :: sym
           
          
          TYPE(t_vacuum),INTENT(in):: vacuum
          TYPE(t_input),INTENT(in)::input

          INTEGER    method1, method2
          PARAMETER (method1 = 2, method2 = 1)
          !     ..
          !     .. Scalar Arguments ..

          INTEGER,INTENT (IN) :: jspin       
          LOGICAL,INTENT (IN) :: l_st
          !     ..
          !     .. Array Arguments ..
          COMPLEX,INTENT(IN),OPTIONAL :: vpw(:,:)
          REAL,INTENT(IN),OPTIONAL :: vr(:,0:,:,:)
          COMPLEX,INTENT (INOUT) :: qpw(stars%ng3,input%jspins)
          COMPLEX,INTENT (INOUT) :: rhtxy(vacuum%nmzxyd,stars%ng2-1,2,input%jspins)
          REAL,   INTENT (INOUT) :: rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins)
          REAL,   INTENT (INOUT) :: rht(vacuum%nmzd,2,input%jspins)
          REAL,   INTENT (INOUT) :: rh(atoms%msh,atoms%ntype)
          !     ..
          !     .. Local Scalars ..
          COMPLEX czero,carg,VALUE,slope,c_ph,sm
          REAL    dif,dxx,g,gz,dtildh,&
               &        rkappa,sign,signz,tol_14,z,zero,zvac,&
               &        g2,phi,gamma,qq,s13,s23,factor,gr
          INTEGER ig3,imz,ivac,j,j1,k,kz,k1,k2,l_cutoff,m0,&
               &        n,nz,nrz,nzvac,&
               &        irec2,irec3,irec1,m,gzi,dir,jm,lh,nat,left,minpar,kp,l,lm,maxl,nd,symint
          !     ..
          !     .. Local Arrays ..
          COMPLEX, ALLOCATABLE :: qpwc(:),ffonat_pT(:,:),pylm2(:,:,:)
          REAL, ALLOCATABLE :: vr2(:,:,:),integrand(:,:),integrandr(:),vrrgrid(:),vrigrid(:),bsl(:,:),force_a4_mt_loc(:,:)
          REAL    acoff(atoms%ntype),alpha(atoms%ntype),rho_out(2)
          REAL    rat(atoms%msh,atoms%ntype)
          COMPLEX gv(3),ycomp1(3,-1:1),ffonat(3,stars%ng3*sym%nop)
          INTEGER mshc(atoms%ntype),ioffset_pT(0:fmpi%isize-1),nkpt_pT(0:fmpi%isize-1)
          LOGICAL l_f2
          INTEGER, ALLOCATABLE :: n1(:),n2(:)
          !     ..
          DATA  czero /(0.0,0.0)/, zero /0.0/, tol_14 /1.0e-10/!-14
#ifdef CPP_MPI
      !EXTERNAL MPI_BCAST
      INTEGER ierr
#endif

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

          ALLOCATE (qpwc(stars%ng3))

          l_f2 = input%l_f.AND.(input%f_level.GE.1).AND.(.NOT.l_st) ! f_level >= 1: coretails completely contained in force calculation
                                                                    ! Klueppelberg (force level 1)
          IF (l_f2) THEN
          ! Allocate the force arrays in the routine force_a4_add.f90
             IF (.NOT.ALLOCATED(force_a4_mt)) THEN
                CALL alloc_fa4_arrays(atoms,input)
             END IF
             ALLOCATE(force_a4_mt_loc,mold=force_a4_mt(:,:,jspin))
             force_a4_mt(:,:,jspin) =  zero
             force_a4_mt_loc(:,:) =  zero
             force_a4_is(:,:,jspin) = czero

             ! Calculate the distribution of Tasks onto processors
             ! Agenda here: every processor in general gets the same number of tasks,
             ! but if there are some left, distribute them on the last processors
             ! reason for that is that within the first few stars, the length of the
             ! stars will change rapidly, while for the last few stars, many will have
             ! the same length. Hence, the first few processors, which will calculate
             ! the spherical Bessel functions quite often, don't get as many stars to
             ! calculate as the last few processors
             ! the first star is \vec{0} and will not contribute to the core forces
             ! thereforce, stars are only considered starting from the second star.
             ! the number of stars is then nq3-1

             ! TODO: Proper parallelization of Klueppelberg force levels.
             
             !ioffset_pT = 0
             !minpar = (stars%ng3-1)/fmpi%isize       ! MINimal number of elements calculated in each PARallel rank
             !nkpt_pT = minpar
             !left = (stars%ng3-1) - minpar * fmpi%isize
             !do j=1,left
             !   nkpt_pT(fmpi%isize-j) = nkpt_pT(fmpi%isize-j)+1
             !end do !j

             !do j=1,fmpi%isize-1
             !   ioffset_pT(j) = sum(nkpt_pT(0:j-1))
             !end do !j

             !ioffset_pT = ioffset_pT+1
             ALLOCATE ( ffonat_pT(3,(stars%ng3-1)*sym%nop) )
             ffonat_pT = czero

             ! lattice/spherical harmonics related variables
             s13 = sqrt(1.0/3.0)
             s23 = sqrt(2.0/3.0)
             ycomp1(1,0) = czero
             ycomp1(2,0) = czero
             ycomp1(3,0) = cmplx(2.0*s13,0.0)
             ycomp1(1,-1) = cmplx(s23,0.0)
             ycomp1(2,-1) = cmplx(0.0,-s23)
             ycomp1(3,-1) = czero
             ycomp1(1,1) = cmplx(-s23,0.0)
             ycomp1(2,1) = cmplx(0.0,-s23)
             ycomp1(3,1) = czero

             ALLOCATE ( vr2(atoms%jmtd,0:sphhar%nlhd,atoms%ntype) )
             vr2=0.0
             ! the l = 0 component of the potential is multiplied by r/sqrt(4 pi), 
             ! for simple use, this is corrected here
             DO n = 1,atoms%ntype
                vr2(:atoms%jri(n),0,n) = sfp_const*vr(:atoms%jri(n),0,n,jspin)/atoms%rmsh(:atoms%jri(n),n)
                vr2(:,1:,n) = vr(:,1:,n,jspin)
             END DO ! n

             ALLOCATE ( integrandr(atoms%jmtd) )
             ALLOCATE ( pylm2( (atoms%lmaxd+1)**2,3,sym%nop ) )
             ALLOCATE ( vrrgrid(atoms%jmtd),vrigrid(atoms%jmtd) )

             ! (f)orce(f)actor(on)(at)oms calculation, parallelization in k
             ffonat = czero
             DO n = 1, atoms%ntype
                nat = atoms%firstAtom(n)
                nd = sym%ntypsy(nat)
                ! find maximal l of the potential for atom (type) n
                ! directly reading max(llh(:,nd)) is only possible if llh is initialized to zero
                ! otherwise, there can be random numbers in it for high lh that are not used by each atom
                maxl = 0
                DO lh = 0,sphhar%nlh(nd)
                   maxl = max(maxl,sphhar%llh(lh,nd))
                END DO ! lh
                ALLOCATE ( bsl(0:maxl,atoms%jmtd),integrand(atoms%jmtd,0:maxl) ) ; bsl=0.0
                g = -0.1 ! g is the norm of a star and can't be negative, this is to initialize a check if the norm between stars has changed

                ! on each processor, calculate a certain consecutive set of k
                kp = 0
                DO k = 2,stars%ng3 ! for k = 1 (G = 0), grad rho_core^alpha is zero
                   IF (abs(g-stars%sk3(k)).gt.tol_14) THEN ! only calculate new spherical Bessel functions if the length of the star vector has changed
                      g = stars%sk3(k)

                      ! generate spherical Bessel functions up to maxl for the radial grid
                      DO j = 1,atoms%jri(n)
                         gr = g * atoms%rmsh(j,n)
                         CALL sphbes(maxl,gr,bsl(:,j))
                         bsl(:,j) = bsl(:,j) * atoms%rmsh(j,n)**2
                      END DO ! j
                   END IF

                   ! as phasy1, but with i\vec{G} in it, i.e. gradient of plane wave, only for atom n and star k
                   CALL phasy2(atoms, stars, sym, cell, k, n, nat, pylm2)

                   ! construct and evaluate radial integral int_0^R_{beta} r^2 j_{l}(Gr) V_{eff,l}^{beta}(r) dr
                   !then, multiply by pylm2 times clnu
                   DO lh = 0,sphhar%nlh(nd)
                      l = sphhar%llh(lh,nd)
                      integrandr(:) = bsl(l,:) * vr2(:,lh,n)
                      CALL intgr3(integrandr,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),factor)
                      DO j = 1,sym%nop
                         symint = kp*sym%nop + j
                         DO dir = 1,3
                            sm = czero
                            DO jm = 1,sphhar%nmem(lh,nd)
                               lm = l*(l+1) + sphhar%mlh(jm,lh,nd) + 1
                               sm = sm + conjg(sphhar%clnu(jm,lh,nd)) * pylm2(lm,dir,j)
                            END DO ! jm
                            ffonat_pT(dir,symint) = ffonat_pT(dir,symint) + factor * sm
                         END DO ! dir
                      END DO ! symint
                   END DO ! lh

                   kp = kp+1
                END DO ! k stars
                DEALLOCATE ( bsl,integrand )
             END DO ! n atom type

             ! collect the entries of ffonat calculated by the different processors
             ! could be all collapsed to irank 0, but if the later part also gets
             ! parallelized over k at some point, now all iranks have ffonat available
!#ifdef CPP_MPI
!             ALLOCATE( n1(0:fmpi%isize-1),n2(0:fmpi%isize-1) )
!             n1(:) = 3*nkpt_pT(:)*sym%nop
!             n2(:) = 3*ioffset_pT(:)*sym%nop
!             CALL MPI_ALLGATHERV(ffonat_pT(1,1),n1(fmpi%irank),MPI_COMPLEX,ffonat(1,1),n1,n2,MPI_COMPLEX,MPI_COMM_WORLD,ierr)
!             DEALLOCATE(ffonat_pT,n1,n2)
!#else
             ffonat(:,(sym%nop+1):) = ffonat_pT
             IF (fmpi%irank==0) THEN
                DEALLOCATE(ffonat_pT)
             END IF
!#endif

          END IF ! l_f2

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
#ifdef CPP_MPI
          CALL MPI_BCAST(rh,atoms%msh*atoms%ntype,MPI_DOUBLE_PRECISION,0,fmpi%mpi_comm,ierr)
#endif
          mshc(:) = 0 ! This initialization is important because there may be atoms without core states.
          nloop: DO  n = 1 , atoms%ntype
              IF ((atoms%econf(n)%num_core_states.GT.0).OR.l_st) THEN
                   DO  j = 1 , atoms%jri(n)
                      rat(j,n) = atoms%rmsh(j,n)
                   ENDDO
                   dxx = EXP(atoms%dx(n))
                   DO j = atoms%jri(n) + 1 , atoms%msh
                      rat(j,n) = rat(j-1,n)*dxx
                   ENDDO
                   DO j = atoms%jri(n) - 1 , atoms%msh
                      rh(j,n) = rh(j,n)/ (fpi_const*rat(j,n)*rat(j,n))
                   ENDDO
                   DO j = atoms%msh , atoms%jri(n) , -1
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
              nat = atoms%firstAtom(n)
              nd = sym%ntypsy(nat)
              IF ((mshc(n).GT.atoms%jri(n)).AND.((atoms%econf(n)%num_core_states.GT.0).OR.l_st)) THEN

                   j1 = atoms%jri(n) - 1
                   IF ( method1 .EQ. 1) THEN
                      dif = diflgr(rat(j1,n),rh(j1,n))
                      WRITE (oUnit,FMT=8000) n,rh(atoms%jri(n),n),dif
                      alpha(n) = -0.5 * dif / ( rh(atoms%jri(n),n)*atoms%rmt(n) )
                   ELSEIF ( method1 .EQ. 2) THEN
                      alpha(n) = LOG( rh(j1,n) / rh(atoms%jri(n),n) )
                      alpha(n) = alpha(n)&
                           &                   / ( atoms%rmt(n)*atoms%rmt(n)*( 1.0-EXP( -2.0*atoms%dx(n) ) ) )
                   ELSE
                      WRITE (oUnit,'('' error in choice of method1 in cdnovlp '')')
                      CALL juDFT_error("error in choice of method1 in cdnovlp"&
                           &              ,calledby ="cdnovlp")
                   ENDIF
                   acoff(n) = rh(atoms%jri(n),n) * EXP( alpha(n)*atoms%rmt(n)*atoms%rmt(n) )
                   !WRITE (oUnit,FMT=8010) alpha(n),acoff(n)
                   DO j = 1,atoms%jri(n) - 1
                      rh(j,n) = acoff(n) * EXP( -alpha(n)*rat(j,n)**2 )
                   ENDDO

                   ! Subtract pseudo density contribution from own mt sphere from mt forces
                   ! Klueppelberg (force level 1)
                   IF (l_f2) THEN
                      ALLOCATE ( integrand(atoms%jmtd,3) )
                      integrandr = zero
                      integrand  = zero
                      DO j = 1,atoms%jri(n)
                         integrandr(j) = -alpha(n) * acoff(n) * atoms%rmsh(j,n)**3 &
                                         *sfp_const * exp(-alpha(n) * atoms%rmsh(j,n)**2) !*2 
                      ! factor of two missing? grad e^{-alpha*r^2} = -2alpha\vec{r}e^{-alpha*r^2}
                      END DO ! j radial mesh

                      DO lh = 0,sphhar%nlh(nd)
                         IF (sphhar%llh(lh,nd).ne.1) CYCLE

                         gv = czero
                         DO jm = 1,sphhar%nmem(lh,nd)
                            m = sphhar%mlh(jm,lh,nd)

                            DO dir = 1,3
                               gv(dir) = gv(dir) + ycomp1(dir,m)* sphhar%clnu(jm,lh,nd) ! why not conjg?
                            END DO ! direction

                         END DO ! jm

                         DO dir = 1,3
                            DO j = 1,atoms%jri(n)
                               integrand(j,dir) = integrand(j,dir) - integrandr(j)* vr2(j,lh,n) * real(gv(dir))
                            END DO ! j radial mesh
                         END DO ! dir ection

                      END DO ! lh lattice harmonics

                      DO dir = 1,3
                         CALL intgr3(integrand(:,dir),atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),force_a4_mt_loc(dir,n))
                      END DO ! direction
                      DEALLOCATE ( integrand )
                   END IF !l_f2
                ELSE
                   alpha(n) = 0.0
              ENDIF
          ENDDO
          !
          IF (fmpi%irank ==0) THEN
8000         FORMAT (/,10x,'core density and its first derivative ',&
                  &                 'at sph. bound. for atom type',&
                  &             i2,' is',3x,2e15.7)
8010         FORMAT (/,10x,'alpha=',f10.5,5x,'acoff=',f10.5)
          END IF          
          !
          !=====> calculate the fourier transform of the core-pseudocharge
          IF (l_f2) THEN
             CALL ft_of_CorePseudocharge(fmpi,input,atoms,mshc,alpha,tol_14,rh, &
                             acoff,stars,method2,rat,cell ,sym,qpwc,jspin,l_f2,vpw,ffonat,force_a4_mt_loc)
          ELSE
             CALL ft_of_CorePseudocharge(fmpi,input,atoms,mshc,alpha,tol_14,rh, &
                             acoff,stars,method2,rat,cell ,sym,qpwc,jspin,l_f2)
          END IF

          DO k = 1 , stars%ng3    
              qpw(k,jspin) = qpw(k,jspin) + qpwc(k) 
          ENDDO

          IF (fmpi%irank ==0) THEN
             !
             !=====> calculate core-tails to the vacuum region                
             !       Coretails expanded in exponentially decaying functions.
             !       Describe vacuum by: 
             !       rho(r_||,z,ivac)= sum_n{ rho(n,ivac) * exp(-kappa*(z-z_v)) 
             !                                           * exp(iG(n)r_||) }
             IF (input%film ) THEN
                !+gu
                dtildh = 0.5 * tpi_const / cell%bmat(3,3)
                IF (vacuum%nvac.EQ.1) THEN
                   rho_out(1) = qpwc(1)*cell%z1
                   DO k = 2,stars%ng3
                      IF ((stars%kv3(1,k).EQ.0).AND.(stars%kv3(2,k).EQ.0)) THEN
                         nz = stars%nstr(k) !1
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
                            carg = carg -qpwc(k)*(EXP(ImagUnit*g*dtildh)-EXP(ImagUnit*g*cell%z1))/g
                         ENDIF
                      ENDDO
                      rho_out(ivac) = qpwc(1) * ( dtildh-cell%z1 ) - AIMAG(carg)
                   ENDDO
                ENDIF
                !-gu
                !        nzvac = min(50,nmz)
               
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
                      DO 250 kz = -stars%mx3,stars%mx3
                         ig3 = stars%ig(k1,k2,kz)
                         c_ph = stars%rgphs(k1,k2,kz) ! phase factor for invs=T & zrfs=F
                         !        ----> use only stars within the g_max sphere (oct.97 shz)
                         IF (ig3.NE.0) THEN
                            gz = kz*cell%bmat(3,3)
                            carg = ImagUnit*sign*gz
                            VALUE = VALUE + c_ph*qpwc(ig3)* EXP(carg*cell%z1)
                            slope = slope + c_ph*carg*qpwc(ig3)* EXP(carg*cell%z1)
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
                         zvac   = MIN (2.*vacuum%nmz,abs(zvac)) ! avoid int-overflow in next line
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
             END IF
             !
             !=====> update density inside the spheres                        
             !
             ! ----> (1) subtract on-site contribution, because 
             !           they are contained in the plane wave part 
             !
             DO n = 1,atoms%ntype
                IF  ((mshc(n).GT.atoms%jri(n)).AND.((atoms%econf(n)%num_core_states.GT.0).OR.l_st)) THEN
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
          ENDIF ! fmpi%irank ==0
          l_cutoff=input%coretail_lmax
#ifdef CPP_MPI
          IF ( fmpi%isize > 1 ) CALL mpi_bc_st(fmpi,stars,qpwc)
#endif

          CALL qpw_to_nmt(&
               &                sphhar,atoms,stars,&
               &                sym,cell ,fmpi,&
               &                jspin,l_cutoff,qpwc,&
               &                rho)

#ifdef CPP_MPI
          IF ( fmpi%isize > 1) CALL mpi_col_st(fmpi,atoms,sphhar,rho(1,0,1,jspin))
#endif

          DEALLOCATE (qpwc)
          IF (l_f2) THEN ! Klueppelberg (force level 1)
             ! Deallocate arrays used specifically during force calculation
             DEALLOCATE ( vr2,integrandr,pylm2 )
             DEALLOCATE ( vrrgrid,vrigrid )
          END IF

        END SUBROUTINE cdnovlp

!***********************************************************************
!     INTERNAL SUBROUTINES
!***********************************************************************

      subroutine ft_of_CorePseudocharge(fmpi,input,atoms,mshc,alpha,&
            tol_14,rh,acoff,stars,method2,rat,cell ,sym,qpwc,jspin,l_f2,vpw,ffonat,force_a4_mt_loc)

      !=====> calculate the fourier transform of the core-pseudocharge
      !
      !     qpw(\vec{G}) = Sum_{n} [ F(G,n) * Sum_{atm{n}} S(\vec{G},atm) ]
      !                  n = atom_type
      !                  F = Formfactor = F_in_sphere + F_outsphere
      !                  S = Structure factor

      USE m_types

      type(t_mpi)      ,intent(in) :: fmpi
      TYPE(t_input),    INTENT(in) ::input
      type(t_atoms)    ,intent(in) :: atoms
      integer          ,intent(in) :: mshc(atoms%ntype),jspin
      real             ,intent(in) :: alpha(atoms%ntype), tol_14
      real             ,intent(in) :: rh(atoms%msh,atoms%ntype)
      real             ,intent(in) :: acoff(atoms%ntype)
      type(t_stars)    ,intent(in) :: stars
      integer          ,intent(in) :: method2
      real             ,intent(in) :: rat(atoms%msh,atoms%ntype)
      type(t_cell)     ,intent(in) :: cell
       
      type(t_sym)      ,intent(in) :: sym
      LOGICAL,         INTENT(IN)  :: l_f2
      COMPLEX,OPTIONAL,INTENT(IN)  :: vpw(:,:),ffonat(:,:)
      REAL,OPTIONAL,INTENT(IN) :: force_a4_mt_loc(:,:)
      complex         ,intent(out) :: qpwc(stars%ng3)

!     ..Local variables
      integer nat1, n, n_out_p, k
      INTEGER :: reducedStarsCutoff ! This is introduced to avoid numerical instabilities.
      complex czero

!     ..Local arrays
      real :: qf(stars%ng3)
      complex qpwc_at(stars%ng3)
#ifdef CPP_MPI
      complex :: qpwc_loc(stars%ng3)
      integer :: ierr
#endif
      czero = (0.0,0.0)
#ifdef CPP_MPI
      DO k = 1, stars%ng3
         qpwc_loc(k) = czero
      ENDDO
#endif
      DO k = 1, stars%ng3
          IF (stars%sk3(k).LE.3.0*input%rkmax) reducedStarsCutoff = k ! The factor 3.0 is arbitrary. One could try going down to 2.0.
          qpwc(k) = czero
      ENDDO

      !
      !*****> start loop over the atom type
      !
      DO  n = 1 + fmpi%irank, atoms%ntype, fmpi%isize
          IF ( ( mshc(n) .GT. atoms%jri(n) ).AND.&
              &        ( alpha(n) .GT. tol_14 ) )    THEN
                   
              n_out_p = mshc(n)-atoms%jri(n)+1
              
              ! (1) Form factor for each atom type
             
              CALL FormFactor_forAtomType(atoms%msh,method2,n_out_p,reducedStarsCutoff,&
                                 atoms%rmt(n),atoms%jri(n),atoms%dx(n),mshc(n),rat(:,n), &
                                 rh(:,n),alpha(n),stars,cell,acoff(n),qf)

              ! (2) structure constant for each atom of atom type
              
              nat1 = atoms%firstAtom(n)
              
              IF (l_f2) THEN     
                 CALL StructureConst_forAtom(nat1,stars ,sym,reducedStarsCutoff,&
                                    atoms%neq(n),atoms%nat,atoms%taual,&
                                    cell,qf,qpwc_at,jspin,l_f2,n,vpw,ffonat)
              ELSE
                 CALL StructureConst_forAtom(nat1,stars ,sym,reducedStarsCutoff,&
                                    atoms%neq(n),atoms%nat,atoms%taual,&
                                    cell,qf,qpwc_at,jspin,l_f2,n)              
              END IF
#ifdef CPP_MPI
              DO k = 1, stars%ng3
                 qpwc_loc(k) = qpwc_loc(k)  + qpwc_at(k)
              END DO
#else
              DO k = 1 , stars%ng3    
                 qpwc(k) = qpwc(k) + qpwc_at(k)
              END DO
#endif

          END IF
       ENDDO
#ifdef CPP_MPI
       CALL mpi_allreduce(qpwc_loc,qpwc,stars%ng3,MPI_DOUBLE_COMPLEX,mpi_sum, &
               fmpi%mpi_comm,ierr)
       IF (l_f2) THEN
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,force_a4_mt,SIZE(force_a4_mt),MPI_DOUBLE,MPI_SUM,fmpi%mpi_comm,ierr)
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,force_a4_is,SIZE(force_a4_is),MPI_DOUBLE_COMPLEX,MPI_SUM,fmpi%mpi_comm,ierr)
       END IF
#endif
       IF (l_f2) THEN
          force_a4_mt(:,:,jspin)=force_a4_mt(:,:,jspin)+force_a4_mt_loc(:,:)
       END IF
      end subroutine ft_of_CorePseudocharge

   SUBROUTINE StructureConst_forAtom(nat1,stars ,sym,reducedStarsCutoff,&
                          neq,natd,taual,cell,qf,qpwc_at,jspin,l_f2,n,vpw,ffonat)
      ! Calculates the structure constant for each atom of atom type

      USE m_types
      USE m_spgrot
      USE m_constants
       

      integer,       intent(in)  :: nat1
      type(t_stars), intent(in)  :: stars
       
      type(t_sym),   intent(in)  :: sym
      INTEGER,       INTENT(IN)  :: reducedStarsCutoff
      integer,       intent(in)  :: neq,natd, jspin, n
      real,          intent(in)  :: taual(3,natd)
      type(t_cell),  intent(in)  :: cell
      real,          intent(in)  :: qf(stars%ng3)
      LOGICAL,       INTENT(IN)  :: l_f2
      COMPLEX,OPTIONAL,INTENT(IN):: vpw(:,:),ffonat(:,:)
      complex,       intent(out) :: qpwc_at(stars%ng3)

      ! ..Local variables
      integer k, nat2, nat, j
      real x
      complex sf, czero

      ! ..Local arrays
      integer kr(3,sym%nop)
      real    force_mt_loc(3)
      complex phas(sym%nop), phase, force_is_loc(3)
      complex  kcmplx(3)

      czero = (0.0,0.0)
      qpwc_at(:) = czero

      !    first G=0
      k=1
      qpwc_at(k)      = qpwc_at(k)      + neq * qf(k)

      !    then  G>0

      force_mt_loc=0.0
      force_is_loc=cmplx(0.0,0.0)
!$OMP PARALLEL DO DEFAULT(none) &
!$OMP SHARED(stars ,sym,reducedStarsCutoff,neq,natd,nat1,taual,cell,qf,qpwc_at,l_f2,ffonat,n,jspin,vpw) &
!$OMP FIRSTPRIVATE(czero) &
!$OMP PRIVATE(k,kr,phas,nat2,nat,sf,j,x,kcmplx,phase) &
!$OMP REDUCTION(+:force_mt_loc,force_is_loc)
      DO  k = 2,reducedStarsCutoff
         
            CALL spgrot(sym%nop, sym%symor, sym%mrot, sym%tau, sym%invtab, &
                        stars%kv3(:,k), kr, phas)
            
            ! ----> start loop over equivalent atoms

            IF (l_f2) THEN ! Klueppelberg (force level 1)
               ! generate phase factors for each G, not only for each star, to incorporate the atomic phase factors
               kcmplx = cmplx(0.0,0.0)
               DO j = 1,sym%nop
                  x = -tpi_const * ( kr(1,j) * taual(1,nat1) &
                                   + kr(2,j) * taual(2,nat1) &
                                   + kr(3,j) * taual(3,nat1) )
                  phase = cmplx(cos(x),sin(x))
                  ! generate muffin-tin part of core force component
                  force_mt_loc(:) = force_mt_loc(:) + qf(k) * &
                                           phase * stars%nstr(k) * ffonat(:,(k-1)*sym%nop+j)
                  kcmplx(:) = kcmplx(:) + kr(:,j) * phase * phas(j) ! should be conjg(phas(j)), but in FLEUR, only real phas(j) are accepted
               END DO !j
               kcmplx = matmul(kcmplx,cell%bmat) * stars%nstr(k) / sym%nop
               ! generate interstitial part of core force component
               force_is_loc(:) = force_is_loc(:) + qf(k) * &
                                        conjg(vpw(k,jspin))*cell%omtil*ImagUnit*kcmplx(:)
            END IF

            nat2 = nat1 + neq - 1
            DO nat = nat1, nat2
               sf = czero
               DO j = 1,sym%nop
                  x = -tpi_const * ( kr(1,j) * taual(1,nat) &
                                   + kr(2,j) * taual(2,nat) &
                                   + kr(3,j) * taual(3,nat) )
                       !gb      sf = sf + CMPLX(COS(x),SIN(x))*phas(j)
                  sf = sf + CMPLX(COS(x),SIN(x))*conjg(phas(j))
               END DO
               sf = sf / REAL( sym%nop )
               qpwc_at(k) = qpwc_at(k) + sf * qf(k)
            END DO
         
      ENDDO
!$OMP END PARALLEL DO

      IF (l_f2) THEN ! Klueppelberg (force level 1)
         force_a4_mt(:,n,jspin) = force_a4_mt(:,n,jspin) + force_mt_loc
         force_a4_is(:,n,jspin) = force_a4_is(:,n,jspin) + force_is_loc
      END IF
   END SUBROUTINE StructureConst_forAtom

   SUBROUTINE FormFactor_forAtomType(msh, method2, n_out_p, reducedStarsCutoff, rmt, jri, dx, &
                                     mshc, rat, rh, alpha, stars, cell, acoff, &
                                     qf)

      USE m_types
      USE m_constants
      USE m_rcerf
      USE m_intgr, ONLY : intgr3, intgz0

      
      integer          ,intent(in) :: msh,method2, n_out_p
      INTEGER,          INTENT(IN) :: reducedStarsCutoff
      real             ,intent(in) :: rmt
      integer          ,intent(in) :: jri
      real             ,intent(in) :: dx
      integer          ,intent(in) :: mshc
      real             ,intent(in) :: rat(msh)
      real             ,intent(in) :: rh(msh)
      real             ,intent(in) :: alpha
      type(t_stars)    ,intent(in) :: stars
      type(t_cell)     ,intent(in) :: cell
      real             ,intent(in) :: acoff
      real            ,intent(out) :: qf(stars%ng3)

      ! ..Local variables
      real f11, f12, ar, g, ai, qfin, qfout, gr, a4, alpha3, zero
      integer k, ir, j
      logical tail

      ! ..Local arrays
      real rhohelp(msh)

      zero = 0.0
      qf(:) = 0.0

      tail = .FALSE.
      f11 = tpi_const * rmt * rh(jri) / alpha
      f12 = acoff * ( pi_const/alpha ) *SQRT(pi_const/alpha)
      ar  = SQRT( alpha ) * rmt 

!$OMP PARALLEL DO DEFAULT(none) & 
!$OMP SHARED(stars,f11,f12,ar,method2,n_out_p,reducedStarsCutoff,jri,rat,rh,dx,tail) &
!$OMP SHARED(alpha,cell,mshc,rmt,qf) &
!$OMP FIRSTPRIVATE(zero) &
!$OMP PRIVATE(k,g,ai,qfin,ir,j,rhohelp,qfout,gr,a4,alpha3)
      DO k = 1, reducedStarsCutoff
         g = stars%sk3(k)
         !    first G=0
         IF ( k.EQ.1 ) THEN
            ai = zero

            ! ---->     calculate form factor inside the mt-sphere
            !           (use analytic integration of gaussian)
 
            qfin = - f11 + f12 * rcerf(ar,ai)

            ! ---->     calculate form factor outside the mt-sphere
            !           (do numerical integration of tails)

            IF ( method2 .EQ. 1) THEN
               DO ir = -6 , n_out_p
                  j = jri+ir-1
                  rhohelp(mshc+1-j) =  rat(j) * rat(j) * rat(j) *  rh(j)
               END DO
               CALL intgz0(rhohelp, dx, n_out_p, qfout, tail)
            ELSE
               DO ir = 1 , n_out_p
                  j = jri+ir-1
                  rhohelp(ir) = rat(j) * rat(j) * rh(j)
               END DO
               CALL intgr3(rhohelp, rat(jri), dx, n_out_p, qfout)
               ! ---->     have to remove the small r-correction from intgr3
               qfout = qfout - rmt*rhohelp(1)
            END IF
 
            qfout = fpi_const * qfout
         ELSE 
            !    then  G>0
            ai = 0.5*g/SQRT(alpha)
            gr = g*rmt
            a4 = 0.25/alpha

            ! ---->     calculate form factor inside the mt-sphere
            !           (use analytic integration of gaussian)

              qfin = - f11 * SIN(gr)/gr + f12 * rcerf(ar,ai) * EXP(-a4*g*g) 

            ! ---->     calculate form factor outside the mt-sphere
            !           (do numerical integration of tails)

            IF ( method2 .EQ. 1) THEN
               DO ir = -6 , n_out_p
                  j  = jri+ir-1
                  rhohelp(mshc-jri+2-ir) =  rat(j)*rat(j) * rh(j) * SIN( g*rat(j) )
               END DO

               ! ---->     note we use here the integration routine for vacuum. Because 
               !           the vacuum integration is made for an inwards integration 
               !           from outside to inside. Outside the starting value will be 
               !           nearly zero since the core density is small. if the tail 
               !           correction (tail=.true.) is included for the integrals, the 
               !           integrand is from infinity inward. This might not be 
               !           necessary. Further the integration routine is made for 
               !           equidistant meshpoints, therefore the term r(i) of
               !           dr/di = dx*r(i) is included in rhohelp

               CALL intgz0(rhohelp,dx,n_out_p,qfout,tail)
            ELSE
               DO ir = 1 , n_out_p
                  j  = jri+ir-1
                  rhohelp(ir) = rat(j) * rh(j) * SIN(g*rat(j))
               END DO
               CALL intgr3(rhohelp, rat(jri), dx, n_out_p, qfout)
               ! ---->     have to remove the small r-correction from intgr3
               !roa...correction.from.intgr3.......................
               IF (rhohelp(1)*rhohelp(2).GT.zero) THEN
                  alpha3 = 1.0 + LOG(rhohelp(2)/rhohelp(1))/dx
                  IF (alpha3.GT.zero) qfout = qfout - rat(jri)*rhohelp(1)/alpha3
               END IF
               !roa...end.correction...............................
            END IF

            qfout = fpi_const * qfout / g
         END IF
         qf(k) = (qfin + qfout)/cell%omtil
      END DO
!$OMP END PARALLEL DO
   END SUBROUTINE FormFactor_forAtomType
END MODULE m_cdnovlp
