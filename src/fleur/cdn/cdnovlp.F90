!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
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
   implicit none

   PRIVATE
   PUBLIC :: cdnovlp, cdnovlp_noco
CONTAINS
   SUBROUTINE cdnovlp(fmpi, sphhar, stars, atoms, sym, vacuum,cell, input, &
                        l_st, jspin, rh, qpw, rho, rhvac, l_CutHighG, vpw, vr)
      !--------------------------------------------------------------------------
      !     Collinear case: the core tail of spin jspin is a single field which is
      !     added to the spin component jspin of the interstitial, muffin-tin and
      !     vacuum densities. All the work is done in add_coretail.
      !--------------------------------------------------------------------------
      USE m_types

      TYPE(t_mpi),INTENT(IN)      :: fmpi
      TYPE(t_sphhar),INTENT(IN)   :: sphhar
      TYPE(t_atoms),INTENT(IN)    :: atoms
      TYPE(t_stars),INTENT(IN)    :: stars
      TYPE(t_cell),INTENT(IN)     :: cell
      TYPE(t_sym),INTENT(IN)      :: sym
      TYPE(t_vacuum),INTENT(IN)   :: vacuum
      TYPE(t_input),INTENT(IN)    :: input

      INTEGER,INTENT (IN) :: jspin
      LOGICAL,INTENT (IN) :: l_st, l_CutHighG
      COMPLEX,INTENT(IN),OPTIONAL :: vpw(:,:)
      REAL,INTENT(IN),OPTIONAL :: vr(:,0:,:,:)
      COMPLEX,INTENT (INOUT) :: qpw(stars%ng3,input%jspins)
      COMPLEX,INTENT (INOUT) :: rhvac(vacuum%nmzd,stars%ng2,2,input%jspins)
      REAL,   INTENT (INOUT) :: rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins)
      REAL,   INTENT (INOUT) :: rh(atoms%msh,atoms%ntype)

      REAL    :: wsrc(1,1,atoms%ntype),wmt(1,1,atoms%ntype)
      COMPLEX :: cone(1,1)
      REAL,ALLOCATABLE :: rh_chan(:,:,:)

      wsrc = 1.0
      wmt  = 1.0
      cone = CMPLX(1.0,0.0)

      !--->    a single spin channel; the copy keeps the (modifying) contract of
      !        cdnovlp for the two dimensional rh
      ALLOCATE(rh_chan(SIZE(rh,1),SIZE(rh,2),1))
      rh_chan(:,:,1) = rh

      CALL add_coretail(fmpi,sphhar,stars,atoms,sym,vacuum,cell,input,l_st,l_CutHighG,jspin,&
                        rh_chan,wsrc,cone,wmt,&
                        qpw(:,jspin:jspin),rho(:,:,:,jspin:jspin),rhvac(:,:,:,jspin:jspin),vpw,vr)

      rh = rh_chan(:,:,1)

   END SUBROUTINE cdnovlp

   SUBROUTINE cdnovlp_noco(fmpi,sphhar,stars,atoms,sym,vacuum,cell,input,noco,nococonv,&
                           l_st,rh,qpw,rho,rhvac,l_CutHighG)
      !--------------------------------------------------------------------------
      !     Non-collinear case.
      !
      !     The core density comes out of the core solver as a spherical,
      !     spin-resolved density that is diagonal in the LOCAL spin frame of the
      !     respective atom type (given by nococonv%alph/beta), while the
      !     interstitial and vacuum density matrices are stored in the GLOBAL
      !     frame. The core tails are therefore carried as the four fields
      !        (n, m_x, m_y, m_z)
      !     with the charge n and the magnetization vector m in global cartesian
      !     coordinates: the core density matrix of atom type n is
      !        rho^n(r) = ( n_n(r) + m_n(r) e_n.sigma ) / 2
      !     with n_n = rh_up + rh_dn, m_n = rh_up - rh_dn and e_n the local
      !     quantization axis of atom type n in global coordinates.
      !
      !     Since nococonv%alph/beta are given per atom TYPE, each of these four
      !     fields is real in real space and has exactly the same form factor,
      !     structure factor and star symmetry as the charge density itself. All
      !     of the machinery of the collinear case can therefore be used
      !     unchanged; only the weights of the spin channels (wsrc) and the
      !     mappings onto the components of the density matrices (cpw, wmt)
      !     differ, and setting these up is all this routine does. In particular
      !     the non-linear steps (tail cutoff, pseudo charge fit) are still
      !     applied to the positive spin-resolved densities and never to a
      !     (possibly small or negative) magnetization.
      !--------------------------------------------------------------------------
      USE m_constants
      USE m_types

      TYPE(t_mpi),INTENT(IN)      :: fmpi
      TYPE(t_sphhar),INTENT(IN)   :: sphhar
      TYPE(t_stars),INTENT(IN)    :: stars
      TYPE(t_atoms),INTENT(IN)    :: atoms
      TYPE(t_sym),INTENT(IN)      :: sym
      TYPE(t_vacuum),INTENT(IN)   :: vacuum
      TYPE(t_cell),INTENT(IN)     :: cell
      TYPE(t_input),INTENT(IN)    :: input
      TYPE(t_noco),INTENT(IN)     :: noco
      TYPE(t_nococonv),INTENT(IN) :: nococonv
      LOGICAL,INTENT(IN)          :: l_st, l_CutHighG
      REAL,INTENT(INOUT)          :: rh(:,:,:)      !(atoms%msh,atoms%ntype,2)
      COMPLEX,INTENT(INOUT)       :: qpw(:,:)       !(stars%ng3,3)
      REAL,INTENT(INOUT)          :: rho(:,0:,:,:)  !(jmtd,0:nlhd,ntype,2 or 4)
      COMPLEX,INTENT(INOUT)       :: rhvac(:,:,:,:) !(nmzd,ng2,2,3)

      INTEGER :: n,i,ndestPw,ndestMt,ndestVac
      REAL    :: mag(0:3),rotm(3,3)
      REAL    :: wsrc(4,2,atoms%ntype)
      REAL    :: wmt(4,SIZE(rho,4),atoms%ntype)
      COMPLEX :: cpw(4,SIZE(qpw,2))

      IF (SIZE(rh,3).LT.2) CALL juDFT_error("cdnovlp_noco needs both spins",calledby="cdnovlp_noco")

      ndestPw  = MIN(3,SIZE(qpw,2))      !the fourth component only exists for DFPT
      ndestMt  = SIZE(rho,4)             !4 components only if l_mperp=T
      ndestVac = MIN(ndestPw,SIZE(rhvac,4))

      !--->    weights of the two spin channels: the charge adds up, the
      !        magnetization enters with the local quantization axis of the atom
      !        type expressed in global cartesian coordinates
      DO n = 1, atoms%ntype
         mag = 0.0
         mag(3) = 1.0
         CALL nococonv%rot_magvec(n,mag,toGlobal=.TRUE.)
         wsrc(1,:,n) = 1.0
         wsrc(2:4,1,n) =  mag(1:3)
         wsrc(2:4,2,n) = -mag(1:3)
      END DO

      ! The transverse (in-plane) part of the core-tail magnetization is not
      ! lattice periodic in the rotating frame of a spin spiral: around each
      ! atom it acquires a factor exp(-i q r), so its form factors would have
      ! to be evaluated at |G+q| for individual G vectors instead of at |G| for
      ! stars (q also breaks the star symmetry). This is not implemented, so
      ! for spin spirals only the charge and the longitudinal (global z) part of
      ! the core-tail magnetization are added. Both are lattice periodic and are
      ! therefore treated exactly. What is dropped is the core magnetization
      ! outside the MT spheres, which is tiny - but note that for a flat spiral
      ! (beta = pi/2) that is the whole core-tail magnetization.
      IF (noco%l_ss) wsrc(2:3,:,:) = 0.0

      !--->    mapping onto the interstitial and vacuum density matrices. Both
      !        are stored in the global frame with the convention (cf.
      !        rotate_int_den_tofrom_local): m_x =  2*REAL(qpw(:,3)),
      !        m_y = -2*AIMAG(qpw(:,3)), m_z = qpw(:,1)-qpw(:,2)
      cpw = CMPLX(0.0,0.0)
      cpw(1,1) = 0.5 ; cpw(4,1) =  0.5                            !rho_11 = (n+m_z)/2
      cpw(1,2) = 0.5 ; cpw(4,2) = -0.5                            !rho_22 = (n-m_z)/2
      IF (ndestPw.GE.3) THEN
         cpw(2,3) = 0.5 ; cpw(3,3) = -0.5*ImagUnit                !rho_21 = (m_x-i*m_y)/2
      END IF
      !--->    mapping onto the muffin-tin density matrix of each atom type. Here
      !        the magnetization has to be rotated into the local frame of the
      !        atom type: rotm(j,i) is the j-th local component of the i-th global
      !        cartesian unit vector. It is determined with the rotation routines
      !        of t_nococonv to be sure that the same conventions as everywhere
      !        else in FLEUR are used.
      wmt = 0.0
      DO n = 1, atoms%ntype
         DO i = 1, 3
            mag = 0.0
            mag(i) = 1.0
            CALL nococonv%rot_magvec(n,mag,toGlobal=.FALSE.)
            rotm(:,i) = mag(1:3)
         END DO
         wmt(1,1,n) = 0.5 ; wmt(2:4,1,n) =  0.5*rotm(3,:)         !rho_11 = (n+m_z^loc)/2
         wmt(1,2,n) = 0.5 ; wmt(2:4,2,n) = -0.5*rotm(3,:)         !rho_22 = (n-m_z^loc)/2
         IF (ndestMt.GE.4) THEN
            wmt(2:4,3,n) = 0.5*rotm(1,:)                          !REAL(rho_21)  = m_x^loc/2
            wmt(2:4,4,n) = 0.5*rotm(2,:)                          !AIMAG(rho_21) = m_y^loc/2
         END IF
      END DO

      CALL add_coretail(fmpi,sphhar,stars,atoms,sym,vacuum,cell,input,l_st,l_CutHighG,1,&
                        rh,wsrc,cpw(:,:ndestPw),wmt,&
                        qpw(:,:ndestPw),rho,rhvac(:,:,:,:ndestVac))

   END SUBROUTINE cdnovlp_noco

   SUBROUTINE add_coretail(fmpi,sphhar,stars,atoms,sym,vacuum,cell,input,l_st,l_CutHighG,jspin,&
                           rh,wsrc,cpw,wmt,qpw,rho,rhvac,vpw,vr)
      !--------------------------------------------------------------------------
      !     calculates the overlapping core tail density and adds
      !     its contribution to the corresponding components of the
      !     valence density.
      !
      !     The subroutine is based on the
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
      !     coded                  Stefan Bl"ugel, IFF Nov. 1997
      !     tested                 RObert Abt    , IFF Dez. 1997
      !
      !     Added calculation of force contributions from coretails
      !     outside of their native muffin-tin spheres, i.e. in the
      !     interstitial region and other muffin-tins; only for bulk.
      !     Refer to Klüppelberg et al., PRB 91 035105 (2015)
      !     Aaron Klueppelberg, Oct. 2015
      !
      !     Generalized to an arbitrary number of source fields and destination
      !     components to cover the non-collinear case as well, see cdnovlp and
      !     cdnovlp_noco for the two entry points.
      !--------------------------------------------------------------------------
      USE m_constants
      USE m_qpwtonmt
      USE m_types

      TYPE(t_mpi),INTENT(IN)      :: fmpi
      TYPE(t_sphhar),INTENT(IN)   :: sphhar
      TYPE(t_atoms),INTENT(IN)    :: atoms
      TYPE(t_stars),INTENT(IN)    :: stars
      TYPE(t_cell),INTENT(IN)     :: cell
      TYPE(t_sym),INTENT(IN)      :: sym
      TYPE(t_vacuum),INTENT(IN)   :: vacuum
      TYPE(t_input),INTENT(IN)    :: input

      LOGICAL,INTENT(IN) :: l_st, l_CutHighG
      INTEGER,INTENT(IN) :: jspin                  !only used to index the force arrays
      REAL,   INTENT(INOUT) :: rh(:,:,:)           !(atoms%msh,atoms%ntype,nchan)
      REAL,   INTENT(IN)    :: wsrc(:,:,:)         !(nsrc,nchan,atoms%ntype)
      COMPLEX,INTENT(IN)    :: cpw(:,:)            !(nsrc,SIZE(qpw,2)); also used for rhvac
      REAL,   INTENT(IN)    :: wmt(:,:,:)          !(nsrc,SIZE(rho,4),atoms%ntype)
      COMPLEX,INTENT(INOUT) :: qpw(:,:)            !(stars%ng3,ndestPw)
      REAL,   INTENT(INOUT) :: rho(:,0:,:,:)       !(jmtd,0:nlhd,ntype,ndestMt)
      COMPLEX,INTENT(INOUT) :: rhvac(:,:,:,:)      !(nmzd,ng2,2,ndestVac)
      COMPLEX,INTENT(IN),OPTIONAL :: vpw(:,:)
      REAL,   INTENT(IN),OPTIONAL :: vr(:,0:,:,:)

      !     method2  : two different integration routines to calculate form
      !                factor of coretails outside the sphere.
      !                (1) use subroutine intgrz to integrate the tails from
      !                    outside to inside.
      !                (2) use subroutine intgr3 to integrate the tails from
      !                    muffin-tin radius to outside and include correction
      !                    for start up.
      !                Tests have shown that (1) is more accurate.
      INTEGER,PARAMETER :: method2 = 1
      !     mshc     : maximal radial meshpoint for which the radial coretail
      !                density is larger than tol_14
      REAL,PARAMETER    :: tol_14 = 1.0e-10 !-14

      INTEGER :: ichan,nchan,nsrc,i,j,k,n
      LOGICAL :: l_f2
      REAL    :: acoff(atoms%ntype),alpha(atoms%ntype)
      REAL    :: rat(atoms%msh,atoms%ntype)
      INTEGER :: mshc(atoms%ntype,SIZE(rh,3))
      COMPLEX,ALLOCATABLE :: qpwc(:,:),ffonat(:,:)
      REAL,ALLOCATABLE    :: force_a4_mt_loc(:,:),rbuf(:,:,:,:)
#ifdef CPP_MPI
      INTEGER :: ierr
#endif

      nchan = SIZE(rh,3)
      nsrc  = SIZE(wsrc,1)

      l_f2 = input%l_f.AND.(input%f_level.GE.1).AND.(.NOT.l_st) & ! f_level >= 1: coretails completely contained in force calculation
             .AND.PRESENT(vpw).AND.PRESENT(vr)                    ! Klueppelberg (force level 1)
      IF (l_f2) THEN
         IF (nchan.GT.1) CALL juDFT_error("Coretail forces are only implemented for the collinear case",&
                                          calledby="add_coretail")
         ! Allocate the force arrays in the routine force_a4_add.f90
         IF (.NOT.ALLOCATED(force_a4_mt)) CALL alloc_fa4_arrays(atoms,input)
         ALLOCATE(force_a4_mt_loc,mold=force_a4_mt(:,:,jspin))
         force_a4_mt(:,:,jspin) = 0.0
         force_a4_mt_loc(:,:)   = 0.0
         force_a4_is(:,:,jspin) = CMPLX(0.0,0.0)
         CALL fa4_formfactors(fmpi,atoms,sphhar,stars,sym,cell,vr(:,:,:,jspin),tol_14,ffonat)
      ELSE
         !dummies: they are passed on but never touched without the forces
         ALLOCATE(ffonat(1,1),force_a4_mt_loc(1,1))
      END IF

      ALLOCATE(qpwc(stars%ng3,nsrc))
      qpwc = CMPLX(0.0,0.0)

#ifdef CPP_MPI
      ! The core density is only calculated on the first rank, but all ranks
      ! need it because the form factors are distributed over the atom types.
      ! rh is copied here to be sure that a contiguous buffer is passed.
      IF (fmpi%isize.GT.1) THEN
         ALLOCATE(rbuf(SIZE(rh,1),SIZE(rh,2),SIZE(rh,3),1))
         rbuf(:,:,:,1) = rh
         CALL MPI_BCAST(rbuf,SIZE(rbuf),MPI_DOUBLE_PRECISION,0,fmpi%mpi_comm,ierr)
         rh = rbuf(:,:,:,1)
         DEALLOCATE(rbuf)
      END IF
#endif

      DO ichan = 1, nchan
         !----> (1) set up radial mesh beyond muffin-tin radius
         !      (2) cut off core tails from noise
         !      (3) replace the core density inside the sphere by a pseudo charge
         CALL prep_coretail(atoms,l_st,rh(:,:,ichan),rat,mshc(:,ichan),alpha,acoff)

         ! Subtract pseudo density contribution from own mt sphere from mt forces
         ! Klueppelberg (force level 1)
         IF (l_f2) CALL fa4_mt_pseudo(atoms,sphhar,sym,l_st,vr(:,:,:,jspin),mshc(:,ichan),&
                                      alpha,acoff,force_a4_mt_loc)

         !=====> calculate the fourier transform of the core-pseudocharge
         CALL ft_of_CorePseudocharge(fmpi,input,atoms,mshc(:,ichan),alpha,tol_14,rh(:,:,ichan),&
                                     acoff,stars,method2,rat,cell,sym,qpwc,jspin,l_f2,&
                                     wgt=wsrc(:,ichan,:),vpw=vpw,ffonat=ffonat,&
                                     force_a4_mt_loc=force_a4_mt_loc)
      END DO
      ! qpwc is identical on all ranks (it comes out of an all-reduce in
      ! ft_of_CorePseudocharge), so no additional broadcast is needed here.

      IF (l_CutHighG) THEN
         DO k = 1, stars%ng3
            IF (stars%sk3(k).GT.(2.0*input%rkmax)) qpwc(k,:) = CMPLX(0.0,0.0)
         END DO
      END IF

      !=====> add the core tails to the interstitial density
      DO i = 1, SIZE(qpw,2)
         qpw(:,i) = qpw(:,i) + MATMUL(qpwc,cpw(:,i))
      END DO

      IF (fmpi%irank.EQ.0) THEN
         !=====> calculate core-tails to the vacuum region
         !the vacuum density matrix is stored with the same convention as the
         !interstitial one, so the same mapping is used
         IF (input%film) CALL coretail_to_vacuum(stars,cell,vacuum,qpwc,&
                                                 cpw(:,:SIZE(rhvac,4)),rhvac)

         !=====> update density inside the spheres
         !
         ! ----> (1) subtract on-site contribution, because they are contained in
         !           the plane wave part. In the non-collinear case the core
         !           density of an atom is diagonal in its own local frame, so
         !           channel ichan is subtracted from component ichan and there
         !           is no off-diagonal term here.
         DO ichan = 1, nchan
            DO n = 1,atoms%ntype
               IF ((mshc(n,ichan).GT.atoms%jri(n)).AND.((atoms%econf(n)%num_core_states.GT.0).OR.l_st)) THEN
                  DO j = 1,atoms%jri(n)
                     rho(j,0,n,ichan) = rho(j,0,n,ichan)&
                          &                          - sfp_const*rat(j,n)*rat(j,n)*rh(j,n,ichan)
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF ! fmpi%irank ==0

      ! ----> (2) add the plane wave contribution of (core tails + on-site
      !           contribution) to the m.t. density, include full nonspherical
      !           components
      CALL qpw_to_nmt(sphhar,atoms,stars,sym,cell,fmpi,input%coretail_lmax,qpwc,rho,wmt)

#ifdef CPP_MPI
      IF (fmpi%isize.GT.1) THEN
         ! qpw_to_nmt zeroed rho on all ranks but the first one, so the sum
         ! contains the valence density once and the core tails of all ranks
         ALLOCATE(rbuf(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3),SIZE(rho,4)))
         rbuf = rho
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,rbuf,SIZE(rbuf),MPI_DOUBLE_PRECISION,MPI_SUM,fmpi%mpi_comm,ierr)
         rho = rbuf
         DEALLOCATE(rbuf)
      END IF
#endif

   END SUBROUTINE add_coretail

!***********************************************************************
!     INTERNAL SUBROUTINES
!***********************************************************************

   SUBROUTINE prep_coretail(atoms,l_st,rh,rat,mshc,alpha,acoff)
      !--------------------------------------------------------------------------
      !     (1) set up radial mesh beyond muffin-tin radius
      !     (2) cut off core tails from noise
      !     (3) the core density inside the spheres is replaced by a
      !         gaussian-like pseudo density : n(r) = acoff*exp(-alpha*r*r)
      !         acoff and alpha determined to obtain a continous and
      !         differentiable density at the sphere boundary.
      !         IF mshc = jri  either core tail too small or no core (i.e. H)
      !
      !     rh enters as 4*pi*r^2*rho_core and leaves as the (pseudized) rho_core
      !--------------------------------------------------------------------------
      USE m_constants
      USE m_types
      USE m_diflgr

      TYPE(t_atoms),INTENT(IN) :: atoms
      LOGICAL,INTENT(IN)       :: l_st
      REAL,INTENT(INOUT)       :: rh(:,:)   !(atoms%msh,atoms%ntype)
      REAL,INTENT(OUT)         :: rat(:,:)  !(atoms%msh,atoms%ntype)
      INTEGER,INTENT(OUT)      :: mshc(:)   !(atoms%ntype)
      REAL,INTENT(OUT)         :: alpha(:),acoff(:)

      !     method1  : two different ways to calculate the derivative of the
      !                charge density at the sphere boundary.
      !                (1) use subroutine diflgr based on lagrange interpol.
      !                (2) use two point formular in real space,
      !                    see notes of SB.
      !                Tests have shown that (2) is more accurate.
      INTEGER,PARAMETER :: method1 = 2
      REAL,PARAMETER    :: tol_14 = 1.0e-10 !-14

      INTEGER :: n,j,j1
      REAL    :: dxx,dif

      rat = 0.0
      alpha = 0.0
      acoff = 0.0
      mshc(:) = 0 ! This initialization is important because there may be atoms without core states.

      nloop: DO n = 1, atoms%ntype
         IF ((atoms%econf(n)%num_core_states.GT.0).OR.l_st) THEN
            DO j = 1, atoms%jri(n)
               rat(j,n) = atoms%rmsh(j,n)
            ENDDO
            dxx = EXP(atoms%dx(n))
            DO j = atoms%jri(n) + 1, atoms%msh
               rat(j,n) = rat(j-1,n)*dxx
            ENDDO
            DO j = atoms%jri(n) - 1, atoms%msh
               rh(j,n) = rh(j,n)/ (fpi_const*rat(j,n)*rat(j,n))
            ENDDO
            DO j = atoms%msh, atoms%jri(n), -1
               IF ( rh(j,n) .GT. tol_14 ) THEN
                  mshc(n) = j
                  CYCLE nloop
               END IF
            ENDDO
            mshc(n) = atoms%jri(n)
         ENDIF
      ENDDO nloop

      DO n = 1, atoms%ntype
         IF ((mshc(n).GT.atoms%jri(n)).AND.((atoms%econf(n)%num_core_states.GT.0).OR.l_st)) THEN
            j1 = atoms%jri(n) - 1
            IF ( method1 .EQ. 1) THEN
               dif = diflgr(rat(j1:j1+2,n),rh(j1:j1+2,n))
               WRITE (oUnit,FMT=8000) n,rh(atoms%jri(n),n),dif
               alpha(n) = -0.5 * dif / ( rh(atoms%jri(n),n)*atoms%rmt(n) )
            ELSEIF ( method1 .EQ. 2) THEN
               alpha(n) = LOG( rh(j1,n) / rh(atoms%jri(n),n) )
               alpha(n) = alpha(n)&
                    &                   / ( atoms%rmt(n)*atoms%rmt(n)*( 1.0-EXP( -2.0*atoms%dx(n) ) ) )
            ELSE
               WRITE (oUnit,'('' error in choice of method1 in cdnovlp '')')
               CALL juDFT_error("error in choice of method1 in cdnovlp"&
                    &              ,calledby ="prep_coretail")
            ENDIF
            acoff(n) = rh(atoms%jri(n),n) * EXP( alpha(n)*atoms%rmt(n)*atoms%rmt(n) )
            DO j = 1, atoms%jri(n) - 1
               rh(j,n) = acoff(n) * EXP( -alpha(n)*rat(j,n)**2 )
            ENDDO
         ENDIF
      ENDDO

8000  FORMAT (/,10x,'core density and its first derivative ',&
           &                 'at sph. bound. for atom type',&
           &             i2,' is',3x,2e15.7)

   END SUBROUTINE prep_coretail

   SUBROUTINE coretail_to_vacuum(stars,cell,vacuum,qpwc,cfac,rhvac)
      !--------------------------------------------------------------------------
      !     Calculate core-tails to the vacuum region.
      !     Coretails expanded in exponentially decaying functions. Describe
      !     vacuum by:
      !     rho(r_||,z,ivac)= sum_n{ rho(n,ivac) * exp(-kappa*(z-z_v))
      !                                          * exp(iG(n)r_||) }
      !
      !     qpwc(:,i) holds the plane wave representation of the i-th field. The
      !     fields are mapped onto the components of the density matrix by cfac:
      !     rhvac(...,j) = sum_i cfac(i,j) * (field i)
      !     and the procedure above is applied to each of these components.
      !
      !     For the diagonal components of the density matrix (j=1,2) that is
      !     straightforward: they are charge-like (positive) densities, so their
      !     logarithmic derivative is a meaningful decay constant. This also
      !     reproduces the collinear case exactly. The off-diagonal component
      !     (j=3) is complex and not sign definite, so its logarithmic derivative
      !     is not useful. Since |rho_21| <= sqrt(rho_11*rho_22), the mean of the
      !     two diagonal decay constants is used for it.
      !--------------------------------------------------------------------------
      USE m_constants
      USE m_types

      TYPE(t_stars),INTENT(IN)  :: stars
      TYPE(t_cell),INTENT(IN)   :: cell
      TYPE(t_vacuum),INTENT(IN) :: vacuum
      COMPLEX,INTENT(IN)        :: qpwc(:,:)      !(stars%ng3,nsrc)
      COMPLEX,INTENT(IN)        :: cfac(:,:)      !(nsrc,ndst)
      COMPLEX,INTENT(INOUT)     :: rhvac(:,:,:,:) !(vacuum%nmzd,stars%ng2,2,ndst)

      REAL,PARAMETER    :: tol_14 = 1.0e-10
      COMPLEX,PARAMETER :: czero = (0.0,0.0)

      INTEGER :: ndst,idst,ndiag,k,k1,k2,kz,ivac,ig3,imz,nz,nzvac
      REAL    :: dtildh,g,gz,sign,z,zvac
      REAL    :: rho_out(2,SIZE(cfac,2)),rkappa(SIZE(cfac,2))
      COMPLEX :: carg,c_ph
      COMPLEX :: vval(SIZE(cfac,2)),slope(SIZE(cfac,2))
      COMPLEX,ALLOCATABLE :: qdst(:,:)

      ndst = SIZE(cfac,2)
      !--->    the diagonal components of the density matrix; only these have a
      !        meaningful logarithmic derivative
      ndiag = MIN(2,ndst)

      !--->    plane wave representation of the components of the density matrix
      ALLOCATE(qdst(stars%ng3,ndst))
      qdst = MATMUL(qpwc,cfac)

      dtildh = 0.5 * tpi_const / cell%bmat(3,3)
      IF (vacuum%nvac.EQ.1) THEN
         rho_out(1,:) = qdst(1,:)*cell%z1
         DO k = 2,stars%ng3
            IF ((stars%kv3(1,k).EQ.0).AND.(stars%kv3(2,k).EQ.0)) THEN
               nz = stars%nstr(k) !1
               g = stars%kv3(3,k) * cell%bmat(3,3)
               rho_out(1,:) = rho_out(1,:) + nz*qdst(k,:)*SIN(g*cell%z1)/g
            ENDIF
         ENDDO
         rho_out(1,:) =  qdst(1,:) * dtildh - rho_out(1,:)
      ELSE
         DO ivac = 1, vacuum%nvac
            DO idst = 1, ndst
               carg = czero
               DO k = 2,stars%ng3
                  IF ((stars%kv3(1,k).EQ.0).AND.(stars%kv3(2,k).EQ.0)) THEN
                     g = stars%kv3(3,k) * cell%bmat(3,3) * (3. - 2.*ivac)
                     carg = carg -qdst(k,idst)*(EXP(ImagUnit*g*dtildh)-EXP(ImagUnit*g*cell%z1))/g
                  ENDIF
               ENDDO
               rho_out(ivac,idst) = qdst(1,idst) * ( dtildh-cell%z1 ) - AIMAG(carg)
            ENDDO
         ENDDO
      ENDIF

      !---> loop over 2D stars
      DO k = 1,stars%ng2
         k1 = stars%kv2(1,k)
         k2 = stars%kv2(2,k)
         DO ivac = 1,vacuum%nvac
            vval = czero
            slope = czero
            sign = 3. - 2.*ivac
            ! ---> sum over gz-stars
            DO kz = -stars%mx3,stars%mx3
               ig3 = stars%ig(k1,k2,kz)
               c_ph = stars%rgphs(k1,k2,kz) ! phase factor for invs=T & zrfs=F
               !        ----> use only stars within the g_max sphere (oct.97 shz)
               IF (ig3.NE.0) THEN
                  gz = kz*cell%bmat(3,3)
                  carg = ImagUnit*sign*gz
                  vval(:)  = vval(:)  + c_ph*qdst(ig3,:)* EXP(carg*cell%z1)
                  slope(:) = slope(:) + c_ph*carg*qdst(ig3,:)* EXP(carg*cell%z1)
               END IF
            ENDDO

            !--->       decay constants of the diagonal components
            rkappa = 0.0
            DO idst = 1, ndiag
               ! roa work-around
               IF (  ABS(REAL(vval(idst))).GT.0.0 ) THEN
                  ! roa work-around
                  ! gb works also around
                  rkappa(idst) = - REAL( slope(idst)/vval(idst) )
                  IF (k.EQ.1) rkappa(idst) = vval(idst)/rho_out(ivac,idst)
                  !               rkappa = - sign * real( slope/value )
                  IF (rkappa(idst).LE.0.0) rkappa(idst)=MIN(rkappa(idst),-tol_14)
                  IF (rkappa(idst).GT.0.0) rkappa(idst)=MAX(rkappa(idst),tol_14)
                  ! gb works also around
               END IF
            ENDDO
            !--->       |rho_21| <= sqrt(rho_11*rho_22), so the off-diagonal
            !           component decays with the mean of the two diagonal
            !           decay constants
            IF (ndst.GT.ndiag) rkappa(ndiag+1:) = SUM(rkappa(:ndiag))/ndiag

            DO idst = 1, ndst
               IF ( ABS(REAL(vval(idst))).LE.0.0 ) CYCLE
               !               IF ( rkappa.GT.zero .AND. real(value).GT.zero ) THEN
               IF ( rkappa(idst).LE.0.0 ) CYCLE
               zvac   = - LOG( tol_14/cabs(vval(idst)) ) / rkappa(idst)
               zvac   = MIN (2.*vacuum%nmz,abs(zvac)) ! avoid int-overflow in next line
               nzvac  = INT( zvac/vacuum%delz ) + 1
               z = 0.
               IF ( k.EQ.1 ) THEN
                  DO imz = 1 , MIN( nzvac,vacuum%nmz )
                     rhvac(imz,1,ivac,idst) = rhvac(imz,1,ivac,idst) + vval(idst)*EXP(-rkappa(idst)*z)
                     z = z + vacuum%delz
                  ENDDO
               ELSE
                  DO imz = 1 , MIN( nzvac,vacuum%nmzxy )
                     rhvac(imz,k,ivac,idst) = rhvac(imz,k,ivac,idst) + vval(idst)*EXP(-rkappa(idst)*z)
                     z = z + vacuum%delz
                  ENDDO
               END IF
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE(qdst)

   END SUBROUTINE coretail_to_vacuum

   SUBROUTINE fa4_vr2(atoms,sphhar,vr,vr2)
      !--------------------------------------------------------------------------
      !     The l = 0 component of the potential is multiplied by r/sqrt(4 pi),
      !     for simple use, this is corrected here.
      !     Klueppelberg (force level 1)
      !--------------------------------------------------------------------------
      USE m_constants
      USE m_types

      TYPE(t_atoms),INTENT(IN)  :: atoms
      TYPE(t_sphhar),INTENT(IN) :: sphhar
      REAL,INTENT(IN)           :: vr(:,0:,:)
      REAL,ALLOCATABLE,INTENT(OUT) :: vr2(:,:,:)

      INTEGER :: n

      ALLOCATE ( vr2(atoms%jmtd,0:sphhar%nlhd,atoms%ntype) )
      vr2=0.0
      DO n = 1,atoms%ntype
         vr2(:atoms%jri(n),0,n) = sfp_const*vr(:atoms%jri(n),0,n)/atoms%rmsh(:atoms%jri(n),n)
         vr2(:,1:,n) = vr(:,1:,n)
      END DO ! n

   END SUBROUTINE fa4_vr2

   PURE FUNCTION fa4_ycomp1() RESULT(ycomp1)
      !--->    lattice/spherical harmonics related constants
      !        Klueppelberg (force level 1)
      COMPLEX :: ycomp1(3,-1:1)
      REAL    :: s13,s23

      s13 = SQRT(1.0/3.0)
      s23 = SQRT(2.0/3.0)
      ycomp1 = CMPLX(0.0,0.0)
      ycomp1(3, 0) = CMPLX(2.0*s13,0.0)
      ycomp1(1,-1) = CMPLX(s23,0.0)
      ycomp1(2,-1) = CMPLX(0.0,-s23)
      ycomp1(1, 1) = CMPLX(-s23,0.0)
      ycomp1(2, 1) = CMPLX(0.0,-s23)

   END FUNCTION fa4_ycomp1

   SUBROUTINE fa4_formfactors(fmpi,atoms,sphhar,stars,sym,cell,vr,tol_14,ffonat)
      !--------------------------------------------------------------------------
      !     (f)orce(f)actor(on)(at)oms calculation
      !
      !     Construct and evaluate the radial integral
      !        int_0^R_{beta} r^2 j_{l}(Gr) V_{eff,l}^{beta}(r) dr
      !     and multiply it by pylm2 times clnu. The entries for the first star
      !     (G = 0) stay zero, because grad rho_core^alpha vanishes there.
      !
      !     TODO: Proper parallelization of Klueppelberg force levels. The
      !           original attempt to distribute the stars over the ranks is
      !           preserved in the repository history.
      !
      !     Klueppelberg (force level 1), Oct. 2015
      !--------------------------------------------------------------------------
      USE m_constants
      USE m_types
      USE m_intgr, ONLY : intgr3

      TYPE(t_mpi),INTENT(IN)    :: fmpi
      TYPE(t_atoms),INTENT(IN)  :: atoms
      TYPE(t_sphhar),INTENT(IN) :: sphhar
      TYPE(t_stars),INTENT(IN)  :: stars
      TYPE(t_sym),INTENT(IN)    :: sym
      TYPE(t_cell),INTENT(IN)   :: cell
      REAL,INTENT(IN)           :: vr(:,0:,:)
      REAL,INTENT(IN)           :: tol_14
      COMPLEX,ALLOCATABLE,INTENT(OUT) :: ffonat(:,:)

      INTEGER :: n,nat,nd,maxl,lh,l,j,k,kp,jm,dir,lm,symint
      REAL    :: g,gr,factor
      COMPLEX :: sm
      REAL,ALLOCATABLE    :: vr2(:,:,:),bsl(:,:),integrandr(:)
      COMPLEX,ALLOCATABLE :: pylm2(:,:,:)

      CALL fa4_vr2(atoms,sphhar,vr,vr2)

      ALLOCATE ( ffonat(3,stars%ng3*sym%nop) )
      ALLOCATE ( integrandr(atoms%jmtd) )
      ALLOCATE ( pylm2( (atoms%lmaxd+1)**2,3,sym%nop ) )
      ffonat = CMPLX(0.0,0.0)

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
         ALLOCATE ( bsl(0:maxl,atoms%jmtd) ) ; bsl=0.0
         g = -0.1 ! g is the norm of a star and can't be negative, this is to initialize a check if the norm between stars has changed

         kp = 1 ! the first star is \vec{0} and does not contribute
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

            DO lh = 0,sphhar%nlh(nd)
               l = sphhar%llh(lh,nd)
               integrandr(:) = bsl(l,:) * vr2(:,lh,n)
               CALL intgr3(integrandr,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),factor)
               DO j = 1,sym%nop
                  symint = kp*sym%nop + j
                  DO dir = 1,3
                     sm = CMPLX(0.0,0.0)
                     DO jm = 1,sphhar%nmem(lh,nd)
                        lm = l*(l+1) + sphhar%mlh(jm,lh,nd) + 1
                        sm = sm + conjg(sphhar%clnu(jm,lh,nd)) * pylm2(lm,dir,j)
                     END DO ! jm
                     ffonat(dir,symint) = ffonat(dir,symint) + factor * sm
                  END DO ! dir
               END DO ! symint
            END DO ! lh

            kp = kp+1
         END DO ! k stars
         DEALLOCATE ( bsl )
      END DO ! n atom type

   END SUBROUTINE fa4_formfactors

   SUBROUTINE fa4_mt_pseudo(atoms,sphhar,sym,l_st,vr,mshc,alpha,acoff,force_a4_mt_loc)
      !--------------------------------------------------------------------------
      !     Subtract pseudo density contribution from own mt sphere from mt forces
      !     Klueppelberg (force level 1)
      !--------------------------------------------------------------------------
      USE m_constants
      USE m_types
      USE m_intgr, ONLY : intgr3

      TYPE(t_atoms),INTENT(IN)  :: atoms
      TYPE(t_sphhar),INTENT(IN) :: sphhar
      TYPE(t_sym),INTENT(IN)    :: sym
      LOGICAL,INTENT(IN)        :: l_st
      REAL,INTENT(IN)           :: vr(:,0:,:)
      INTEGER,INTENT(IN)        :: mshc(:)
      REAL,INTENT(IN)           :: alpha(:),acoff(:)
      REAL,INTENT(INOUT)        :: force_a4_mt_loc(:,:)

      INTEGER :: n,nat,nd,j,lh,jm,m,dir
      COMPLEX :: gv(3),ycomp1(3,-1:1)
      REAL,ALLOCATABLE :: integrand(:,:),integrandr(:),vr2(:,:,:)

      CALL fa4_vr2(atoms,sphhar,vr,vr2)
      ycomp1 = fa4_ycomp1()
      ALLOCATE(integrand(atoms%jmtd,3),integrandr(atoms%jmtd))

      DO n = 1, atoms%ntype
         nat = atoms%firstAtom(n)
         nd = sym%ntypsy(nat)
         IF (.NOT.((mshc(n).GT.atoms%jri(n)).AND.((atoms%econf(n)%num_core_states.GT.0).OR.l_st))) CYCLE

         integrandr = 0.0
         integrand  = 0.0
         DO j = 1,atoms%jri(n)
            integrandr(j) = -alpha(n) * acoff(n) * atoms%rmsh(j,n)**3 &
                            *sfp_const * exp(-alpha(n) * atoms%rmsh(j,n)**2) !*2
            ! factor of two missing? grad e^{-alpha*r^2} = -2alpha\vec{r}e^{-alpha*r^2}
         END DO ! j radial mesh

         DO lh = 0,sphhar%nlh(nd)
            IF (sphhar%llh(lh,nd).ne.1) CYCLE

            gv = CMPLX(0.0,0.0)
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
      END DO ! n

      DEALLOCATE(integrand,integrandr,vr2)

   END SUBROUTINE fa4_mt_pseudo

      subroutine ft_of_CorePseudocharge(fmpi,input,atoms,mshc,alpha,&
            tol_14,rh,acoff,stars,method2,rat,cell ,sym,qpwc,jspin,l_f2,wgt,vpw,ffonat,force_a4_mt_loc)

      !=====> calculate the fourier transform of the core-pseudocharge
      !
      !     qpw(\vec{G}) = Sum_{n} [ F(G,n) * Sum_{atm{n}} S(\vec{G},atm) ]
      !                  n = atom_type
      !                  F = Formfactor = F_in_sphere + F_outsphere
      !                  S = Structure factor
      !
      !     The result is ADDED to qpwc, which may have several components
      !     (SIZE(qpwc,2)): the contribution of atom type n is added to
      !     component i with the weight wgt(i,n). This is used for the
      !     non-collinear case, where the (spherical) core density of each
      !     atom type has to be distributed over the components of the
      !     interstitial density matrix according to the local spin frame
      !     of that atom type. Without wgt all weights are 1.0.

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
      REAL,OPTIONAL,   INTENT(IN)  :: wgt(:,:) !(SIZE(qpwc,2),atoms%ntype)
      COMPLEX,OPTIONAL,INTENT(IN)  :: vpw(:,:),ffonat(:,:)
      REAL,OPTIONAL,INTENT(IN) :: force_a4_mt_loc(:,:)
      complex       ,intent(inout) :: qpwc(:,:) !(stars%ng3,ncomp)

!     ..Local variables
      integer nat1, n, n_out_p, k, icomp, ncomp
      real    w
      complex czero

!     ..Local arrays
      real :: qf(stars%ng3)
      complex qpwc_at(stars%ng3)
      complex, ALLOCATABLE :: qpwc_loc(:,:)
#ifdef CPP_MPI
      integer :: ierr
#endif
      czero = (0.0,0.0)
      ncomp = SIZE(qpwc,2)
      ALLOCATE(qpwc_loc(stars%ng3,ncomp))
      qpwc_loc = czero

      !
      !*****> start loop over the atom type
      !
      DO  n = 1 + fmpi%irank, atoms%ntype, fmpi%isize
          IF ( ( mshc(n) .GT. atoms%jri(n) ).AND.&
              &        ( alpha(n) .GT. tol_14 ) )    THEN
                   
              n_out_p = mshc(n)-atoms%jri(n)+1
              
              ! (1) Form factor for each atom type
             
              CALL FormFactor_forAtomType(atoms%msh,method2,n_out_p,&
                                 atoms%rmt(n),atoms%jri(n),atoms%dx(n),mshc(n),rat(:,n), &
                                 rh(:,n),alpha(n),stars,cell,acoff(n),qf)

              ! (2) structure constant for each atom of atom type
              
              nat1 = atoms%firstAtom(n)
              
              IF (l_f2) THEN
                 CALL StructureConst_forAtom(nat1,stars ,sym,&
                                    atoms%neq(n),atoms%nat,atoms%taual,&
                                    cell,qf,qpwc_at,jspin,l_f2,n,vpw,ffonat)
              ELSE
                 CALL StructureConst_forAtom(nat1,stars ,sym,&
                                    atoms%neq(n),atoms%nat,atoms%taual,&
                                    cell,qf,qpwc_at,jspin,l_f2,n)
              END IF
              DO icomp = 1, ncomp
                 w = 1.0
                 IF (PRESENT(wgt)) w = wgt(icomp,n)
                 IF (w.EQ.0.0) CYCLE
                 DO k = 1 , stars%ng3
                    qpwc_loc(k,icomp) = qpwc_loc(k,icomp) + w*qpwc_at(k)
                 END DO
              END DO

          END IF
       ENDDO
#ifdef CPP_MPI
       CALL mpi_allreduce(MPI_IN_PLACE,qpwc_loc,SIZE(qpwc_loc),MPI_DOUBLE_COMPLEX,mpi_sum, &
               fmpi%mpi_comm,ierr)
       IF (l_f2) THEN
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,force_a4_mt,SIZE(force_a4_mt),MPI_DOUBLE,MPI_SUM,fmpi%mpi_comm,ierr)
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,force_a4_is,SIZE(force_a4_is),MPI_DOUBLE_COMPLEX,MPI_SUM,fmpi%mpi_comm,ierr)
       END IF
#endif
       IF (l_f2) THEN
          force_a4_mt(:,:,jspin)=force_a4_mt(:,:,jspin)+force_a4_mt_loc(:,:)
       END IF
       qpwc = qpwc + qpwc_loc
       DEALLOCATE(qpwc_loc)
      end subroutine ft_of_CorePseudocharge

   SUBROUTINE StructureConst_forAtom(nat1,stars ,sym,&
                          neq,natd,taual,cell,qf,qpwc_at,jspin,l_f2,n,vpw,ffonat)
      ! Calculates the structure constant for each atom of atom type

      USE m_types
      USE m_spgrot
      USE m_constants
       

      integer,       intent(in)  :: nat1
      type(t_stars), intent(in)  :: stars
       
      type(t_sym),   intent(in)  :: sym
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
!$OMP SHARED(stars ,sym,neq,natd,nat1,taual,cell,qf,qpwc_at,l_f2,ffonat,n,jspin,vpw) &
!$OMP FIRSTPRIVATE(czero) &
!$OMP PRIVATE(k,kr,phas,nat2,nat,sf,j,x,kcmplx,phase) &
!$OMP REDUCTION(+:force_mt_loc,force_is_loc)
      DO k = 2, stars%ng3
         
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

   SUBROUTINE FormFactor_forAtomType(msh, method2, n_out_p, rmt, jri, dx, &
                                     mshc, rat, rh, alpha, stars, cell, acoff, &
                                     qf)

      USE m_types
      USE m_constants
      USE m_rcerf
      USE m_intgr, ONLY : intgr3, intgz0

      
      integer          ,intent(in) :: msh,method2, n_out_p
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
!$OMP SHARED(stars,f11,f12,ar,method2,n_out_p,jri,rat,rh,dx,tail) &
!$OMP SHARED(alpha,cell,mshc,rmt,qf) &
!$OMP FIRSTPRIVATE(zero) &
!$OMP PRIVATE(k,g,ai,qfin,ir,j,rhohelp,qfout,gr,a4,alpha3)
      DO k = 1, stars%ng3
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

              qfin = - f11 * SIN(gr)/gr + f12 * rcerfMulExp(ar,ai,-a4*g*g)

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
