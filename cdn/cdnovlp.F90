!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_cdnovlp
      USE m_juDFT
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: cdnovlp 
      CONTAINS
        SUBROUTINE cdnovlp(mpi,&
             &                   sphhar,stars,atoms,sym,&
             &                   DIMENSION,vacuum,cell,&
             &                   input,oneD,l_st,&
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
          USE m_constants
          USE m_qpwtonmt
          USE m_cylbes
          USE m_dcylbs
          USE m_od_cylbes
          USE m_diflgr
          USE m_types
#ifdef CPP_MPI
          USE m_mpi_bc_st
#endif
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
          LOGICAL,INTENT (IN) :: l_st
          !     ..
          !     .. Array Arguments ..
          COMPLEX,INTENT (INOUT) :: qpw(stars%n3d,DIMENSION%jspd)
          COMPLEX,INTENT (INOUT) :: rhtxy(vacuum%nmzxyd,oneD%odi%n2d-1,2,DIMENSION%jspd)
          REAL,   INTENT (INOUT) :: rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,DIMENSION%jspd)
          REAL,   INTENT (INOUT) :: rht(vacuum%nmzd,2,DIMENSION%jspd)
          REAL,   INTENT (INOUT) :: rh(DIMENSION%msh,atoms%ntype)
          !     ..
          !     .. Local Scalars ..
          COMPLEX czero,carg,VALUE,slope,ci
          REAL    dif,dxx,g,gz,dtildh,&
               &        rkappa,sign,signz,tol_14,z,zero,zvac,&
               &        g2,phi,gamma,qq
          INTEGER ig3,imz,ivac,j,j1,k,kz,k1,k2,l_cutoff,m0,&
               &        n,nz,nrz,nzvac,&
               &        irec2,irec3,irec1,m,gzi
          !     ..
          !     .. Local Arrays ..
          COMPLEX, ALLOCATABLE :: qpwc(:)
          REAL    acoff(atoms%ntype),alpha(atoms%ntype),rho_out(2)
          REAL    rat(DIMENSION%msh,atoms%ntype)
          INTEGER mshc(atoms%ntype)
          REAL    fJ(-oneD%odi%M:oneD%odi%M),dfJ(-oneD%odi%M:oneD%odi%M)
          !     ..
          DATA  czero /(0.0,0.0)/, zero /0.0/, tol_14 /1.0e-10/!-14
#ifdef CPP_MPI
      EXTERNAL MPI_BCAST
      INTEGER ierr
#include "cpp_double.h"
      INCLUDE "mpif.h"
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
          ci = CMPLX(0.0,1.0)

          ALLOCATE (qpwc(stars%n3d))
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
          CALL MPI_BCAST(rh,DIMENSION%msh*atoms%ntype,CPP_MPI_REAL,0,mpi%mpi_comm,ierr)
#endif
          nloop: DO  n = 1 , atoms%ntype
              IF ((atoms%ncst(n).GT.0).OR.l_st) THEN
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
              IF ((mshc(n).GT.atoms%jri(n)).AND.((atoms%ncst(n).GT.0).OR.l_st)) THEN

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
          IF (mpi%irank ==0) THEN
8000         FORMAT (/,10x,'core density and its first derivative ',&
                  &                 'at sph. bound. for atom type',&
                  &             i2,' is',3x,2e15.7)
8010         FORMAT (/,10x,'alpha=',f10.5,5x,'acoff=',f10.5)
          END IF          
          !
          !=====> calculate the fourier transform of the core-pseudocharge

          CALL ft_of_CorePseudocharge(mpi,DIMENSION,atoms,mshc,alpha,tol_14,rh, &
                          acoff,stars,method2,rat,cell,oneD,sym,qpwc)

          DO k = 1 , stars%ng3    
              qpw(k,jspin) = qpw(k,jspin) + qpwc(k) 
          ENDDO

          IF (mpi%irank ==0) THEN
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
                IF  ((mshc(n).GT.atoms%jri(n)).AND.((atoms%ncst(n).GT.0).OR.l_st)) THEN
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
          ENDIF ! mpi%irank ==0
          l_cutoff=input%coretail_lmax
#ifdef CPP_MPI
          IF ( mpi%isize > 1 ) CALL mpi_bc_st(mpi,stars,qpwc)
#endif

          CALL qpw_to_nmt(&
               &                sphhar,atoms,stars,&
               &                sym,cell,oneD,mpi,&
               &                jspin,l_cutoff,qpwc,&
               &                rho)

#ifdef CPP_MPI
          IF ( mpi%isize > 1) CALL mpi_col_st(mpi,atoms,sphhar,rho(1,0,1,jspin))
#endif

          DEALLOCATE (qpwc)

        END SUBROUTINE cdnovlp

!***********************************************************************
!     INTERNAL SUBROUTINES
!***********************************************************************

      subroutine ft_of_CorePseudocharge(mpi,DIMENSION,atoms,mshc,alpha,&
            tol_14,rh,acoff,stars,method2,rat,cell,oneD,sym,qpwc)

      !=====> calculate the fourier transform of the core-pseudocharge
      !
      !     qpw(\vec{G}) = Sum_{n} [ F(G,n) * Sum_{atm{n}} S(\vec{G},atm) ]
      !                  n = atom_type
      !                  F = Formfactor = F_in_sphere + F_outsphere
      !                  S = Structure factor

      USE m_types

      type(t_mpi)      ,intent(in) :: mpi
      type(t_dimension),intent(in) :: DIMENSION
      type(t_atoms)    ,intent(in) :: atoms
      integer          ,intent(in) :: mshc(atoms%ntype)
      real             ,intent(in) :: alpha(atoms%ntype), tol_14
      real             ,intent(in) :: rh(DIMENSION%msh,atoms%ntype)
      real             ,intent(in) :: acoff(atoms%ntype)
      type(t_stars)    ,intent(in) :: stars
      integer          ,intent(in) :: method2
      real             ,intent(in) :: rat(DIMENSION%msh,atoms%ntype)
      type(t_cell)     ,intent(in) :: cell
      type(t_oneD)     ,intent(in) :: oneD
      type(t_sym)      ,intent(in) :: sym
      complex         ,intent(out) :: qpwc(stars%n3d)

!     ..Local variables
      integer nat1, n, n_out_p, k
      complex czero

!     ..Local arrays
      real :: qf(stars%n3d)
      complex qpwc_at(stars%n3d)
#ifdef CPP_MPI
      external mpi_bcast
      complex :: qpwc_loc(stars%n3d)
      integer :: ierr
#include "cpp_double.h"
      include "mpif.h"
#endif

      czero = (0.0,0.0)
#ifdef CPP_MPI
      DO k = 1 , stars%n3d
         qpwc_loc(k) = czero
      ENDDO
#endif
      DO k = 1 , stars%n3d    
          qpwc(k) = czero
      ENDDO

      !
      !*****> start loop over the atom type
      !
      DO  n = 1 + mpi%irank, atoms%ntype, mpi%isize
          IF ( ( mshc(n) .GT. atoms%jri(n) ).AND.&
              &        ( alpha(n) .GT. tol_14 ) )    THEN
                   
              n_out_p = mshc(n)-atoms%jri(n)+1
              
              ! (1) Form factor for each atom type
             
              CALL FormFactor_forAtomType(DIMENSION,method2,n_out_p,&
                                 atoms%rmt(n),atoms%jri(n),atoms%dx(n),mshc(n),rat(:,n), &
                                 rh(:,n),alpha(n),stars,cell,acoff(n),qf)

              ! (2) structure constant for each atom of atom type
            
              nat1 = 1
              IF (n>1) THEN
                 DO k = 1, n-1
                    nat1 = nat1 + atoms%neq(k)
                 END DO
              END IF       
              CALL StructureConst_forAtom(nat1,stars,oneD,sym,&
                                 atoms%neq(n),atoms%nat,atoms%taual,&
                                 cell,qf,qpwc_at)
#ifdef CPP_MPI
              DO k = 1, stars%n3d
                 qpwc_loc(k) = qpwc_loc(k)  + qpwc_at(k)
              END DO
#else
              DO k = 1 , stars%n3d    
                 qpwc(k) = qpwc(k) + qpwc_at(k)
              END DO
#endif

          END IF
       ENDDO
#ifdef CPP_MPI
       CALL mpi_allreduce(qpwc_loc,qpwc,stars%n3d,CPP_MPI_COMPLEX,mpi_sum, &
               mpi%mpi_comm,ierr)
#endif

      end subroutine ft_of_CorePseudocharge


!----------------------------------------------------------------------
      subroutine StructureConst_forAtom(nat1,stars,oneD,sym,&
                          neq,natd,taual,cell,qf,qpwc_at)
!       calculates structure constant for each atom of atom type

      USE m_types
      USE m_spgrot
      USE m_constants
      USE m_od_chirot

       integer          ,intent(in) :: nat1
       type(t_stars)    ,intent(in) :: stars
       type(t_oneD)     ,intent(in) :: oneD
       type(t_sym)      ,intent(in) :: sym
       integer          ,intent(in) :: neq,natd
       real             ,intent(in) :: taual(3,natd)
       type(t_cell)     ,intent(in) :: cell
       real             ,intent(in) :: qf(stars%n3d)
       complex         ,intent(out) :: qpwc_at(stars%n3d)

!      ..Local variables
      integer k, nat2, nat, j
      real x
      complex sf, czero

!      ..Local arrays
       integer kr(3,sym%nop)
       real    kro(3,oneD%ods%nop)
       complex phas(sym%nop)
       complex phaso(oneD%ods%nop)

      czero = (0.0,0.0)
      DO k = 1 , stars%n3d    
          qpwc_at(k) = czero
      ENDDO

      !
      !    first G=0
      !
      k=1
      qpwc_at(k)      = qpwc_at(k)      + neq * qf(k)

      !
      !    then  G>0
      !
!$OMP PARALLEL DO DEFAULT(none) &
!$OMP SHARED(stars,oneD,sym,neq,natd,nat1,taual,cell,qf,qpwc_at) &
!$OMP FIRSTPRIVATE(czero) &
!$OMP PRIVATE(k,kr,phas,nat2,nat,sf,j,x,kro,phaso)
      DO  k = 2,stars%ng3
          IF (.NOT.oneD%odi%d1) THEN
              CALL spgrot(&
                   &                       sym%nop,sym%symor,sym%mrot,sym%tau,sym%invtab,&
                   &                       stars%kv3(:,k),&
                   &                       kr,phas)
              !
              ! ----> start loop over equivalent atoms
              !
               nat2 = nat1 + neq - 1
               DO  nat = nat1,nat2
                   sf = czero
                   DO  j = 1,sym%nop
                       x = -tpi_const* ( kr(1,j) * taual(1,nat)&
                           &          + kr(2,j) * taual(2,nat)&
                           &          + kr(3,j) * taual(3,nat))
                       !gb      sf = sf + CMPLX(COS(x),SIN(x))*phas(j)
                       sf = sf + CMPLX(COS(x),SIN(x))*conjg(phas(j))
                   ENDDO
                   sf = sf / REAL( sym%nop )
                   qpwc_at(k)      = qpwc_at(k)      + sf * qf(k)
               ENDDO
          ELSE
               !-odim
               CALL od_chirot(oneD%odi,oneD%ods,cell%bmat,stars%kv3(:,k),kro,phaso)
               nat2 = nat1 + neq - 1
               DO  nat = nat1,nat2
                   !                  sf = cmplx(1.,0.)
                   sf = czero
                   DO  j = 1,oneD%ods%nop
                       x = -tpi_const* ( kro(1,j)*taual(1,nat)&
                           &          + kro(2,j)*taual(2,nat)&
                           &          + kro(3,j)*taual(3,nat))
                       sf = sf + CMPLX(COS(x),SIN(x))*phaso(j)
                   ENDDO
                   sf = sf / REAL( oneD%ods%nop )
                   qpwc_at(k)      = qpwc_at(k)      + sf * qf(k)
               ENDDO
               !+odim
          END IF
      ENDDO
!$OMP END PARALLEL DO
      end subroutine StructureConst_forAtom

!----------------------------------------------------------------------
      subroutine FormFactor_forAtomType(DIMENSION,method2,n_out_p,&
                       rmt,jri,dx,mshc,rat,&
                       rh,alpha,stars,cell,acoff,qf)

      USE m_types
      USE m_constants
      USE m_rcerf
      USE m_intgr, ONLY : intgr3,intgz0

      type(t_dimension),intent(in) :: DIMENSION
      integer          ,intent(in) :: method2, n_out_p
      real             ,intent(in) :: rmt
      integer          ,intent(in) :: jri
      real             ,intent(in) :: dx
      integer          ,intent(in) :: mshc
      real             ,intent(in) :: rat(DIMENSION%msh)
      real             ,intent(in) :: rh(DIMENSION%msh)
      real             ,intent(in) :: alpha
      type(t_stars)    ,intent(in) :: stars
      type(t_cell)     ,intent(in) :: cell
      real             ,intent(in) :: acoff
      real            ,intent(out) :: qf(stars%n3d)


!     ..Local variables
      real f11, f12, ar, g, ai, qfin, qfout, gr, a4, alpha3, zero
      integer k, ir, j
      logical tail

!     ..Local arrays
      real rhohelp(DIMENSION%msh)

      zero = 0.0
      DO k = 1,stars%n3d
        qf(k) = 0.0
      END DO

      tail = .FALSE.
      f11 = tpi_const * rmt * rh(jri) / alpha
      f12 = acoff * ( pi_const/alpha ) *SQRT(pi_const/alpha)
      ar  = SQRT( alpha ) * rmt 

!$OMP PARALLEL DO DEFAULT(none) & 
!$OMP SHARED(stars,f11,f12,ar,method2,n_out_p,jri,rat,rh,dx,tail) &
!$OMP SHARED(alpha,cell,mshc,rmt,qf) &
!$OMP FIRSTPRIVATE(zero) &
!$OMP PRIVATE(k,g,ai,qfin,ir,j,rhohelp,qfout,gr,a4,alpha3)
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
                     j  = jri+ir-1
                     rhohelp(mshc+1-j) =  rat(j) * rat(j) &
                           &                       * rat(j) *  rh(j)
                  END DO
                  CALL intgz0(rhohelp,dx,n_out_p,qfout,tail)
 
              ELSE
 
                  DO ir = 1 , n_out_p
                     j  = jri+ir-1
                     rhohelp(ir) = rat(j) * rat(j) * rh(j)
                  END DO
                  CALL intgr3(rhohelp,rat(jri),dx,&
                           &                        n_out_p,qfout)
                  !--->             have to remove the small r-correction from intgr3
                  qfout=qfout-rmt*rhohelp(1)
 
              END IF
 
              qfout = fpi_const * qfout
              !
          ELSE 
              !    then  G>0
              ai = 0.5*g/SQRT(alpha)
              gr = g*rmt
              a4 = 0.25/alpha
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
                      j  = jri+ir-1
                      rhohelp(mshc-jri+2-ir) =  rat(j)*rat(j) &
                                    &                                     * rh(j) * SIN( g*rat(j) )
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


                  CALL intgz0(rhohelp,dx,n_out_p,qfout,tail)

              ELSE

                  DO ir = 1 , n_out_p
                      j  = jri+ir-1
                      rhohelp(ir) = rat(j) * rh(j) * SIN(g*rat(j))
                  END DO
                  CALL intgr3(rhohelp,rat(jri),dx,&
                              &                        n_out_p,qfout)
                  !--->             have to remove the small r-correction from intgr3
                  !roa...correction.from.intgr3.......................
                  IF (rhohelp(1)*rhohelp(2).GT.zero) THEN
                      alpha3 = 1.0 + LOG(rhohelp(2)/rhohelp(1))/dx
                      IF (alpha3.GT.zero)&
                              &                 qfout = qfout - rat(jri)*rhohelp(1)/alpha3
                      ENDIF
                  !roa...end.correction...............................


                  END IF

                  qfout = fpi_const * qfout / g
                  !
          END IF
          !
          qf(k)    = (qfin + qfout)/cell%omtil
      ENDDO
!$OMP END PARALLEL DO
      end subroutine FormFactor_forAtomType

      END MODULE m_cdnovlp

