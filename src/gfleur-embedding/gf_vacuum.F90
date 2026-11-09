!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_vacuum 
      use m_juDFT
      USE m_constants, ONLY: oUnit
      USE m_gf_types
      USE m_gf_energies
      USE m_gf_iodop, ONLY: gf_iodop_readvacemb
#ifdef old_stuff
      USE m_constants, ONLY: pi_const
      USE m_gf_fft_singleton
      USE m_gf_io2dmat
#endif
      IMPLICIT NONE
      PRIVATE 
#ifdef oldstuff
                                                !matrix for coordinate  
      REAL,SAVE               :: coordconv(2,2) 
                                                !transformation, i.e. bm
                                                !kp-vectors             
      INTEGER,POINTER,SAVE    :: k1p(:),k2p(:) 
                                                !k-vectors              
      REAL,Allocatable,SAVE       :: bk(:,:) 
                                                !size of aux-region     
      REAL,SAVE               :: d_aux 
                                                !d of cell              
      REAL,SAVE               :: d 
                                                !vacuum level           
      REAL,SAVE               :: vac_pot 
                                                !position of image plane
      REAL,SAVE               :: imageplane 
                                                !position of image bound
      REAL,SAVE               :: imageboundary 
                                                !EField applied to surfa
      REAL,SAVE               :: EField 
                                                     !do we have a surfa
      LOGICAL,SAVE            :: surface = .FALSE.,l_noco
                                                                        
                                                 !dt of cell            
      REAL,SAVE               :: dt 
                                                 !FFT box in  z directio
      INTEGER,SAVE            :: mx3 
                                                                        
                                                     !do we have an imag
      LOGICAL,SAVE            :: imagepoti = .FALSE. 
                                                                        
                                                  ! #gridpoints of mesh 
      INTEGER,PARAMETER       :: gridpts=1000000 
                                             !stores V(z)               
      REAL, SAVE              :: V(gridpts) 
                                            ! stores x pos. on mesh     
      REAL, SAVE              :: x(gridpts) 
                                                    !mesh spacing       
      REAL, SAVE              :: xdistance,xdist_sq 
                                        !If V=const set                 
      REAL, SAVE              :: V0=0 
                                                                        
                                                                        
      PUBLIC                  :: gf_vacuum_embpot,gf_vacuum_init 
      PUBLIC                  :: gf_vacuum_Poti_init 
!     PUBLIC                  :: gf_vacuum,gf_vacuum_image_embpot       
      PUBLIC                  :: gf_vacuum 
      PUBLIC                  :: gridpts 
      public                  :: gf_vacuum_check_vacpot
#endif
      public gf_simple_vac
      public gf_vacemb_init, gf_vacemb, gf_vacemb_active

      !vacuum embedding at the outer planes of a film (side 1: below
      !layer 1, side 2: above the last layer) from the planar-averaged
      !potential vz(:,spin) beyond the plane, written by gfprep
      TYPE t_vacside
         LOGICAL :: on = .FALSE.
         REAL    :: delz = 0.0
         REAL, ALLOCATABLE :: vz(:, :)
      END TYPE t_vacside
      TYPE(t_vacside), SAVE   :: vacside(2)
      INTEGER, SAVE           :: vac_nlayers = 0
      LOGICAL, SAVE           :: vac_simple = .FALSE.
      REAL, SAVE              :: vac_energy = 0.0, vac_bmat(2, 2) = 0.0
      REAL, ALLOCATABLE, SAVE :: vac_bk(:, :)
      !RK4 steps per grid interval of vz
      INTEGER, PARAMETER      :: nref = 4
                                                                        
      CONTAINS 

      SUBROUTINE gf_vacemb_init(embinp, cell, kpts, nlayers)
      !Sets up the vacuum embedding of the sides selected by surface (side 2)
      !and vacleft (side 1). With simplevacuum the potential beyond the
      !plane is the constant vacuum_energy, otherwise it is read from
      !gf_pot.hdf.
      TYPE(t_embinp), INTENT(IN) :: embinp
      TYPE(t_cell), INTENT(IN)   :: cell
      TYPE(t_kpts), INTENT(IN)   :: kpts
      INTEGER, INTENT(IN)        :: nlayers

      INTEGER :: side, n
      COMPLEX :: e, s1, s2

      vac_nlayers = nlayers
      vacside(1)%on = embinp%l_vacleft
      vacside(2)%on = embinp%l_surface
      vac_simple = embinp%l_simplevacuum
      vac_energy = embinp%vacuum_energy
      vac_bmat = cell%bmat(1:2, 1:2)
      vac_bk = kpts%bk(1:2, :)
      IF (vac_simple) THEN
         WRITE (oUnit, "(a,f12.6)") "vacuum embedding: constant potential ", vac_energy
         RETURN
      END IF
      DO side = 1, 2
         IF (.NOT. vacside(side)%on) CYCLE
         IF (.NOT. gf_iodop_readvacemb(side, vacside(side)%delz, vacside(side)%vz)) &
            CALL juDFT_error("No vacuum embedding potential in gf_pot.hdf", &
                             hint="generate the layers from a film with gfprep or set simplevacuum=t", &
                             calledby="gf_vacemb_init")
         n = SIZE(vacside(side)%vz, 1)
         WRITE (oUnit, "(a,i1,a,i0,a,f8.4,a,2f12.6)") "vacuum embedding side ", side, ": ", n, &
            " points, delz=", vacside(side)%delz, ", V(plane), V(end):", &
            vacside(side)%vz(1, 1), vacside(side)%vz(n, 1)
         !accuracy of the integration: step doubling at an energy below V(plane)
         e = CMPLX(vacside(side)%vz(1, 1) - 0.2, 0.01)
         s1 = priv_vacsigma(vacside(side)%vz(:, 1), vacside(side)%delz, e, 1.0, nref)
         s2 = priv_vacsigma(vacside(side)%vz(:, 1), vacside(side)%delz, e, 1.0, 2*nref)
         WRITE (oUnit, "(a,i1,a,2f10.6,a,es10.2)") "vacuum embedding side ", side, &
            ": Sigma(E=V(plane)-0.2,|k+G|^2=1)=", s2, ", step-doubling difference:", ABS(s1 - s2)
      END DO
      END SUBROUTINE gf_vacemb_init

      LOGICAL FUNCTION gf_vacemb_active(side, region)
      !is the embedding potential of this side and region the vacuum one?
      INTEGER, INTENT(IN) :: side, region
      gf_vacemb_active = .FALSE.
      IF (side == 1 .AND. region == 1) gf_vacemb_active = vacside(1)%on
      IF (side == 2 .AND. region == vac_nlayers) gf_vacemb_active = vacside(2)%on
      END FUNCTION gf_vacemb_active

      SUBROUTINE gf_vacemb(side, en, nk, jspin, lapw_gf, sigma)
      !Vacuum embedding potential, diagonal in the 2D plane waves:
      !Sigma_GG = 1/2 phi'/phi at the plane, phi the solution decaying (or
      !outgoing) into the vacuum of the planar-averaged potential
      INTEGER, INTENT(IN)          :: side, en, nk, jspin
      TYPE(t_lapw_gf), INTENT(IN)  :: lapw_gf
      COMPLEX, INTENT(OUT)         :: sigma(:, :)

      INTEGER :: n, nn, js
      REAL    :: gk(2), q2
      COMPLEX :: e

      e = gf_Z(en, 0)
      sigma = CMPLX(0.0, 0.0)
      nn = SIZE(sigma, 1)
      IF (lapw_gf%nv2(1) < nn) THEN
         CALL priv_vacemb_noco(side, e, nk, lapw_gf, sigma)
         RETURN
      END IF
      DO n = 1, nn
         gk = MATMUL(vac_bk(:, nk) + [lapw_gf%k1p(n, 1), lapw_gf%k2p(n, 1)], vac_bmat)
         q2 = DOT_PRODUCT(gk, gk)
         IF (vac_simple) THEN
            sigma(n, n) = -0.5*SQRT(priv_retarded(2.0*(vac_energy - e) + q2, e))
         ELSE
            js = MIN(jspin, SIZE(vacside(side)%vz, 2))
            sigma(n, n) = priv_vacsigma(vacside(side)%vz(:, js), vacside(side)%delz, e, q2, nref)
         END IF
      END DO
      END SUBROUTINE gf_vacemb

      SUBROUTINE priv_vacemb_noco(side, e, nk, lapw_gf, sigma)
      !Noco: 2x2 spin blocks of Sigma (order up G's, down G's). With the
      !noco vz (V11,V22,Re V21,Im V21) the matrix equation is integrated,
      !a collinear vz gives diag(Sigma_up,Sigma_down).
      INTEGER, INTENT(IN)          :: side, nk
      COMPLEX, INTENT(IN)          :: e
      TYPE(t_lapw_gf), INTENT(IN)  :: lapw_gf
      COMPLEX, INTENT(OUT)         :: sigma(:, :)

      INTEGER :: n, nn, nc
      REAL    :: gk(2), q2
      COMPLEX :: s2(2, 2)

      sigma = CMPLX(0.0, 0.0)
      nn = SIZE(sigma, 1)/2
      DO n = 1, nn
         gk = MATMUL(vac_bk(:, nk) + [lapw_gf%k1p(n, 1), lapw_gf%k2p(n, 1)], vac_bmat)
         q2 = DOT_PRODUCT(gk, gk)
         s2 = CMPLX(0.0, 0.0)
         IF (vac_simple) THEN
            s2(1, 1) = -0.5*SQRT(priv_retarded(2.0*(vac_energy - e) + q2, e))
            s2(2, 2) = s2(1, 1)
         ELSE
            nc = SIZE(vacside(side)%vz, 2)
            IF (nc >= 4) THEN
               s2 = priv_vacsigma2(vacside(side)%vz, vacside(side)%delz, e, q2, nref)
            ELSE
               s2(1, 1) = priv_vacsigma(vacside(side)%vz(:, 1), vacside(side)%delz, e, q2, nref)
               s2(2, 2) = priv_vacsigma(vacside(side)%vz(:, MIN(2, nc)), vacside(side)%delz, e, q2, nref)
            END IF
         END IF
         sigma(n, n) = s2(1, 1)
         sigma(n + nn, n + nn) = s2(2, 2)
         sigma(n + nn, n) = s2(2, 1)
         sigma(n, n + nn) = s2(1, 2)
      END DO
      END SUBROUTINE priv_vacemb_noco

      FUNCTION priv_vacsigma2(vz, delz, e, q2, nsub) RESULT(sig)
      !2x2 version of priv_vacsigma for the spin matrix potential
      !V = [[V11, conjg(V21)], [V21, V22]]: Sigma = 1/2 Phi' Phi^-1 with
      !Phi'' = (2(V(t)-e)+q2) Phi. Beyond the grid Phi = exp(-K t); the
      !columns of (Phi,Phi') are orthonormalized after each step.
      REAL, INTENT(IN)    :: vz(:, :), delz, q2
      COMPLEX, INTENT(IN) :: e
      INTEGER, INTENT(IN) :: nsub
      COMPLEX             :: sig(2, 2)

      COMPLEX :: u(4, 2), k1(4, 2), k2(4, 2), k3(4, 2), k4(4, 2)
      COMPLEX :: w(2, 2), kap(2), pinv(2, 2), det
      REAL    :: h, t, ev(2)
      INTEGER :: i, n, j

      n = SIZE(vz, 1)
      h = delz/nsub
      !K = W diag(kappa) W^+ from the eigenvectors of V at the end
      CALL priv_heig2(vz(n, 1), vz(n, 2), CMPLX(vz(n, 3), vz(n, 4)), ev, w)
      DO j = 1, 2
         kap(j) = SQRT(priv_retarded(2.0*(ev(j) - e) + q2, e))
      END DO
      u = CMPLX(0.0, 0.0)
      u(1, 1) = 1.0; u(2, 2) = 1.0
      u(3:4, :) = -MATMUL(w, MATMUL(RESHAPE([kap(1), CMPLX(0.0, 0.0), CMPLX(0.0, 0.0), kap(2)], [2, 2]), &
                                    TRANSPOSE(CONJG(w))))
      t = (n - 1)*delz
      DO i = (n - 1)*nsub, 1, -1
         k1 = priv_rhs2(t, u)
         k2 = priv_rhs2(t - 0.5*h, u - 0.5*h*k1)
         k3 = priv_rhs2(t - 0.5*h, u - 0.5*h*k2)
         k4 = priv_rhs2(t - h, u - h*k3)
         u = u - h/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4)
         !Gram-Schmidt of the two columns (a right multiplication)
         u(:, 1) = u(:, 1)/SQRT(SUM(ABS(u(:, 1))**2))
         u(:, 2) = u(:, 2) - DOT_PRODUCT(u(:, 1), u(:, 2))*u(:, 1)
         u(:, 2) = u(:, 2)/SQRT(SUM(ABS(u(:, 2))**2))
         t = t - h
      END DO
      det = u(1, 1)*u(2, 2) - u(1, 2)*u(2, 1)
      pinv = RESHAPE([u(2, 2), -u(2, 1), -u(1, 2), u(1, 1)], [2, 2])/det
      sig = 0.5*MATMUL(u(3:4, :), pinv)
      CONTAINS
         FUNCTION priv_rhs2(tt, uu) RESULT(du)
         REAL, INTENT(IN)    :: tt
         COMPLEX, INTENT(IN) :: uu(4, 2)
         COMPLEX             :: du(4, 2)
         COMPLEX :: a(2, 2)
         REAL    :: x
         x = tt/delz + 1.0
         a(1, 1) = 2.0*(priv_vinterp(vz(:, 1), x) - e) + q2
         a(2, 2) = 2.0*(priv_vinterp(vz(:, 2), x) - e) + q2
         a(2, 1) = 2.0*CMPLX(priv_vinterp(vz(:, 3), x), priv_vinterp(vz(:, 4), x))
         a(1, 2) = CONJG(a(2, 1))
         du(1:2, :) = uu(3:4, :)
         du(3:4, :) = MATMUL(a, uu(1:2, :))
         END FUNCTION priv_rhs2
      END FUNCTION priv_vacsigma2

      SUBROUTINE priv_heig2(a, d, b, ev, w)
      !eigenvalues ev and eigenvectors (columns of w) of [[a,conjg(b)],[b,d]]
      REAL, INTENT(IN)     :: a, d
      COMPLEX, INTENT(IN)  :: b
      REAL, INTENT(OUT)    :: ev(2)
      COMPLEX, INTENT(OUT) :: w(2, 2)
      REAL :: m, r, nrm
      m = 0.5*(a + d)
      r = SQRT(0.25*(a - d)**2 + ABS(b)**2)
      IF (ABS(b) < 1.0E-14) THEN
         ev = [a, d]
         w = RESHAPE([CMPLX(1.0, 0.0), CMPLX(0.0, 0.0), CMPLX(0.0, 0.0), CMPLX(1.0, 0.0)], [2, 2])
         RETURN
      END IF
      ev = [m + r, m - r]
      nrm = SQRT((ev(1) - d)**2 + ABS(b)**2)
      w(:, 1) = [CMPLX(ev(1) - d, 0.0), b]/nrm
      w(:, 2) = [-CONJG(b), CMPLX(ev(1) - d, 0.0)]/nrm
      END SUBROUTINE priv_heig2

      COMPLEX FUNCTION priv_vacsigma(vz, delz, e, q2, nsub)
      !1/2 phi'(0)/phi(0) for phi'' = (2(V(t)-e)+q2) phi, t the distance
      !from the plane into the vacuum. Beyond the last grid point V is
      !constant and phi = exp(-kappa t); the linear system (phi,phi') is
      !integrated inwards with RK4, nsub steps per grid interval (also
      !through nodes of phi).
      REAL, INTENT(IN)    :: vz(:), delz, q2
      COMPLEX, INTENT(IN) :: e
      INTEGER, INTENT(IN) :: nsub

      COMPLEX :: u(2), k1(2), k2(2), k3(2), k4(2)
      REAL    :: h, t
      INTEGER :: i, n

      n = SIZE(vz)
      h = delz/nsub
      u(1) = 1.0
      u(2) = -SQRT(priv_retarded(2.0*(vz(n) - e) + q2, e))
      t = (n - 1)*delz
      DO i = (n - 1)*nsub, 1, -1
         k1 = priv_rhs(t, u)
         k2 = priv_rhs(t - 0.5*h, u - 0.5*h*k1)
         k3 = priv_rhs(t - 0.5*h, u - 0.5*h*k2)
         k4 = priv_rhs(t - h, u - h*k3)
         u = u - h/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4)
         u = u/MAXVAL(ABS(u))
         t = t - h
      END DO
      priv_vacsigma = 0.5*u(2)/u(1)
      CONTAINS
         FUNCTION priv_rhs(tt, uu) RESULT(du)
         REAL, INTENT(IN)    :: tt
         COMPLEX, INTENT(IN) :: uu(2)
         COMPLEX             :: du(2)
         du(1) = uu(2)
         du(2) = (2.0*(priv_vinterp(vz, tt/delz + 1.0) - e) + q2)*uu(1)
         END FUNCTION priv_rhs
      END FUNCTION priv_vacsigma

      COMPLEX FUNCTION priv_retarded(f, e)
      !f = kappa^2; on the real energy axis the branch of the retarded
      !Green function (outgoing waves) is selected
      COMPLEX, INTENT(IN) :: f, e
      priv_retarded = f
      IF (ABS(AIMAG(e)) < 1.0E-12) priv_retarded = CMPLX(REAL(f), -1.0E-12)
      END FUNCTION priv_retarded

      REAL FUNCTION priv_vinterp(f, x)
      !cubic Lagrange interpolation of f at the (real) grid index x
      REAL, INTENT(IN) :: f(:), x
      INTEGER :: i0, i, j
      REAL    :: l
      i0 = MIN(MAX(FLOOR(x) - 1, 1), SIZE(f) - 3)
      priv_vinterp = 0.0
      DO i = i0, i0 + 3
         l = 1.0
         DO j = i0, i0 + 3
            IF (j /= i) l = l*(x - j)/(i - j)
         END DO
         priv_vinterp = priv_vinterp + l*f(i)
      END DO
      END FUNCTION priv_vinterp

      subroutine gf_simple_vac(en,nk,lapw,lapw_gf,kpts,cell,vacuum_energy,sigma,efield)
      integer,intent(in)       :: en,nk
      type(t_lapw),intent(in)  :: lapw
      type(t_lapw_gf),intent(in) :: lapw_gf
      type(t_kpts),intent(in)  :: kpts
      type(t_cell),intent(in)  :: cell
      real,intent(in)          :: vacuum_energy
      complex,intent(out)      :: sigma(:,:)
      real,intent(in),optional :: efield

      integer :: nn,n
      real    :: gk(2),s


      nn=size(sigma,1)
      if (lapw_gf%nv2(1)<nn) nn=nn/2 !noco calculation
      if (present(efield)) THEN
         if (abs(efield)>epsilon(0.0)) THEN
            CALL vacuum_EField(en,nk,lapw,lapw_gf,kpts,cell,vacuum_energy,sigma,efield)
            return
         endif
      endif
        sigma=cmplx(0.0,0.0)
        DO n=1,nn
          gk=MATMUL(kpts%bk(1:2,nk)+(/lapw_gf%k1p(n,1),lapw_gf%k2p(n,1)/),cell%bmat(1:2,1:2))
          s=dot_PRODUCT(gk,gk)
          sigma(n,n)=cmplx(0,0.5)* SQRT(2*(gf_Z(en,0)-vacuum_energy)-s)
        enddo
      if (nn.ne.size(sigma,1)) then !noco
        sigma(nn+1:,nn+1:)=sigma(:nn,:nn)
      endif

      end subroutine


       SUBROUTINE vacuum_EField(en,nk,lapw,lapw_gf,kpts,cell,vacuum_energy,sigma,efield)
!******************************************
! calculates numericaly the embedding potential for
! linear E-field E*z+V_vac applyed to the system
!       01/2007            A. Hanuschkin
!       modified by DW 2012
!******************************************
      IMPLICIT NONE
      integer,intent(in)       :: en,nk
      type(t_lapw),intent(in)  :: lapw
      type(t_lapw_gf),intent(in) :: lapw_gf
      type(t_kpts),intent(in)  :: kpts
      type(t_cell),intent(in)  :: cell
      real,intent(in)          :: vacuum_energy
      complex,intent(out)      :: sigma(:,:)
      real,intent(in),optional :: efield
      !-Locals
      REAL    :: gk(2),s
      COMPLEX :: dphi
      INTEGER :: i,n,nn
      INTEGER,PARAMETER       :: gridpts=100000
      REAL                    :: V(gridpts),x(gridpts)
      COMPLEX                 :: phi(gridpts),ksq(gridpts)
      REAL, SAVE              :: xdistance,xdist_sq !mesh spacing


      !Generate potential in vacuum
      xdistance = REAL(100)/REAL(gridpts)
      xdist_sq = xdistance*xdistance

      DO i=1,gridpts
         x(i)=REAL(i-1)*xdistance;
         V(i) = vacuum_energy+(EField*x(i))
      END DO

      !init Embedding potential
      sigma=0.0

      !CALCULATE PHI and DPHI
      nn=size(sigma,1)
      if (lapw_gf%nv2(1)<nn) nn=nn/2 !noco calculation

      DO n=1,nn
           !initial wavefunction
           !inf. in vacuum region V=vac_pot
           gk=MATMUL(kpts%bk(1:2,nk)+(/lapw_gf%k1p(n,1),lapw_gf%k2p(n,1)/),cell%bmat(:2,:2))
           s=dot_PRODUCT(gk,gk)
           !numerical integration ksq=k**2xh**2
                          !set up k^2*h^2 array
           DO i=1,gridpts
               ksq(i) = ( 2.0*(gf_Z(en,0)-V(i))-s)*xdist_sq
           ENDDO
           phi(gridpts) = exp(I*sqrt(2.0*(gf_Z(en,0)-V(gridpts) )-s)*x(gridpts))
           phi(gridpts-1)= exp(I*sqrt(2.0*(gf_Z(en,0)-V(gridpts-1) )-s)*x(gridpts-1))
           DO i=gridpts-2,1,-1
             ! normalize!
              IF (ABS(phi(i+1))>1000.) THEN
                 phi(i+1)=phi(i+1)/ABS(phi(i+1))
                 phi(i+2)=phi(i+2)/ABS(phi(i+1))
              ENDIF
              !calculate phi
               phi(i)=((2.-(10./12.)*ksq(i+1))*phi(i+1)-(1.+(1./12.)*ksq(i+2))*phi(i+2))/(1.+(1./12.)*ksq(i))
           ENDDO
           dphi = -0.5 * ((phi(1)-phi(1+1))/(xdistance)+((phi(1+1)-phi(1+2)))/(xdistance)  )
           sigma(n,n)= 0.5*dphi/phi(1)
      ENDDO
      if (nn.ne.size(sigma,1)) sigma(nn+1:,nn+1:)=sigma(:nn,:nn) !noco


      END SUBROUTINE

      !>
#ifdef old_stuff
      !<-- S: gf_vacuum_embpot(en,nk,g2)                                

      SUBROUTINE gf_vacuum_embpot(en,nk,embpot) 
!******************************************                             
!     Construct the embedding potential for a surface calculation.      
!     The vacuum potential must be set before by calling gf_vacuum_init 
!     Please note that the potential in the aux-volume might not be zero
!     See Eq. 5.50 from my Thesis!                                      
!                          D. Wortmann                                  
!******************************************                             
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN) :: en,nk 
      COMPLEX,INTENT(INOUT) ::embpot(:,:) 
      !>                                                                
      !<--Locals                                                        
      INTEGER :: n,nn
      REAL    :: gk(2),s 
      COMPLEX :: kl,kr,I 
      I=cmplx(0.0,1.0) 
      !>                                                                
                                                                        
                               !no surface calculation                  
      IF (.NOT.surface) RETURN 
                          !cal with imagepoti                           
      IF (imagepoti) THEN 
          CALL gf_vacuum_image_embpot(en,nk,embpot) 
          RETURN 
      ENDIF 
                              !cal with E-field                         
      IF (EField/=0.0) THEN 
          CALL gf_vacuum_EField(en,nk,embpot) 
          RETURN 
      ENDIF 
                                                                        
                                                                        
                                                                        
                               !init Embedding potential                
      embpot=0.0 

      !calculate embedding potential
      ! EmbPot vacuum V=0
      nn=size(embpot,1)
	  if (l_noco) nn=nn/2
      DO n=1,nn
          gk=MATMUL(coordconv,bk(1:2,nk)+(/k1p(n),k2p(n)/))
          s=dot_PRODUCT(gk,gk)
          embpot(n,n)=I/2* SQRT(2*(gf_Z(en,0)-vac_pot)-s)
      enddo
      if (l_noco) then
      	embpot(nn+1:,nn+1:)=embpot(:nn,:nn)
      endif
!      DO n=1,SIZE(embpot,1)
!         gk=MATMUL(coordconv,bk(1:2,nk)+(/k1p(n),k2p(n)/))
!         s=dot_PRODUCT(gk,gk)
!         kl=SQRT(2*gf_Z(en,0)-s)
!         kr=SQRT(2*(gf_Z(en,0)-vac_pot)-s)
!         embpot(n,n)=I/2*kl*(kr*COS(kl*d_aux)+I*kl*SIN(kl*d_aux))/(kl   &
!     &        *COS(kl*d_aux)+I*kr*SIN(kl*d_aux))
!      WRITE (0,*) en, REAL(gf_Z(en)),nk,"S:",embpot(n,n)               
!     &    ,gk,s,REAL(gf_Z(en)),"k:",REAL(kl)                           
!     & ,REAL(kr)                                                       
                                                                        
!       WRITE (0,*) en, REAL(gf_Z(en)),embpot(n,n)                      
!      ENDDO
! "(5i,f8.4,5i,2f12.5,2f8.4,2f8.4)")                                    
                                                                        
      END SUBROUTINE 

      !>                                                                
                                                                        
                                                                        
                                                                        
      !<-- S: gf_vacuum_image_embpot(en,nk,g2)                          

      SUBROUTINE gf_vacuum_image_embpot(en,nk,embpot) 
!******************************************                             
!       12/2006            A. Hanuschkin                                
!******************************************                             
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN) :: en,nk 
      COMPLEX,INTENT(INOUT) ::embpot(:,:) 
      !>                                                                
      !<--Locals                                                        
      INTEGER :: n 
      REAL    :: gk(2),s 
      COMPLEX :: k_r,I 
                                              ! k**2                    
      COMPLEX, ALLOCATABLE :: ksq(:) 
                                              ! derivation of wf        
      COMPLEX :: dphi 
                                              ! wavefunction on mesh    
      COMPLEX , ALLOCATABLE :: phi(:) 
                                              ! save comp. time         
      REAL    :: one12 = (1/12.0) 
      INTEGER :: ii,index_c,zz 
      CHARACTER (LEN=20) :: filename 
                                                                        
      I=cmplx(0.0,1.0) 
      !>                                                                
                                                                        
!     if (.not.surface) return !no surface calculation                  
      !CHECK THAT IMAGEBOUNDARY >= D/2                                  
      !CHECK THAT IMAGEPLANE  < IMAGEBOUNDARY                           
      !CHECK imageplane < embPlane @ reading gf_inp                     
                                                                        
                                                                        
      !allocate and initialize                                          
      ALLOCATE(ksq(gridpts)) 
      ALLOCATE(phi(gridpts)); 
                               !init Embedding potential                
      embpot=0.0 
                                                                        
#ifdef CPP_DEBUG   
                   !not nessecary because it'll be calculated in necessa
                   !but nice if plotted phi(x<c/2)                      
      phi = 0 
                                                                        
      WRITE (*,*) "xdistance,imageboundary,imageplane"                  &
     &            ,xdistance,imageboundary,imageplane                   
#endif                                                                  
                                                                        
                                                                        
       !phi and dphi at z= c/2 =d-d_aux needed                          
       !up to which point integration is needed?                        
                                                        !check imageplan
       index_c = Int((d-d_aux-imageplane)/(xdistance)) 
!       WRITE (*,*) index_c," index_c ",d-d_aux-imageplane              
!     &         ," d-d_aux-imageplane ",d-d_aux," c"                    
                                                                        
       IF ((index_c>gridpts-2).OR.(index_c<=0)) THEN 
         WRITE (*,*) index_c 
         WRITE (*,*) "Imageboundary closer than Embedding plane?" 
          CALL juDFT_error("PROBLEMS WITH SPACING: INDEX_C",calledby="gf_vacuum.F90")
       ENDIF 
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
      !CALCULATE PHI and DPHI at c/2                                    
      DO n=1,SIZE(embpot,1) 
      !initial wavefunction                                             
      !inf. in vacuum region V=vac_pot                                  
                                                                        
       gk=MATMUL(coordconv,bk(1:2,nk)+(/k1p(n),k2p(n)/)) 
       s=dot_PRODUCT(gk,gk) 
                                                                        
                                                                        
#ifdef CPP_DEBUG                                                        
       WRITE (oUnit,*) "s ",s," vac_pot ",vac_pot 
#endif                                                                  
                                                                        
       phi(gridpts) = exp(I*sqrt(2.0*(gf_Z(en,0)-vac_pot)-s)            &
     &                    *x(gridpts))                                  
       phi(gridpts-1)= exp(I*sqrt(2.0*(gf_Z(en,0)-vac_pot)-s)           &
     &                    *x(gridpts-1))                                
       ksq(gridpts) = (2.0*(gf_Z(en,0)-vac_pot)-s)*xdist_sq 
       ksq(gridpts-1)= ksq(gridpts) 
                                                                        
#ifdef CPP_DEBUG_IMAGEPOTI_V0                                           
! WITHOUT POTENTIAL                                                     
       phi(gridpts) = exp(I*sqrt(2.0*(gf_Z(en,0)-V0)-s)                 &
     &                    *x(gridpts))                                  
       phi(gridpts-1) = exp(I*sqrt(2.0*(gf_Z(en,0)-V0)-s)               &
     &                    *x(gridpts-1))                                
      ksq(gridpts) = (2.0*(gf_Z(en,0))-s)*xdist_sq 
      ksq(gridpts-1)= ksq(gridpts) 
#endif                                                                  
                                                                        
                                                                        
       !numerical integration ksq=k**2xh**2                             
                          !set up k^2*h^2 array                         
       DO ii=1,gridpts-2 
        ksq(ii) = ( 2.0*(gf_Z(en,0)-V(ii)) -S)*xdist_sq 
       ENDDO 
                                                                        
                                   !calculate phi                       
       DO ii=gridpts-2,index_c,-1 
         phi(ii)=( (2-(10.0*one12)*ksq(ii+1))*phi(ii+1)                 &
     &             -(1+(1*one12)*ksq(ii+2))*phi(ii+2)  )                &
     &               /(1+(1*one12)*ksq(ii))                             
!         IF (ABS(REAL(phi(ii))).gt.1000) THEN  ! normalize!            
!           phi(ii+1)=phi(ii+1)/ABS(Real(phi(ii)))                      
!           phi(ii+2)=phi(ii+2)/ABS(Real(phi(ii)))  !index_x+2 ii+3?    
                                                                        
                                         ! normalize!                   
         IF (ABS(phi(ii))>1000) THEN 
           phi(ii+1)=phi(ii+1)/ABS(phi(ii)) 
                                             !index_x+2 ii+3?           
           phi(ii+2)=phi(ii+2)/ABS(phi(ii)) 
                                                                        
                                                                        
#ifdef CPP_DEBUG_IP_NA                                                  
                                             ! normalize all!           
            DO zz=gridpts,ii+3,-1 
!           phi(zz)=phi(zz)/ABS(Real(phi(ii)))                          
             phi(zz)=phi(zz)/ABS(phi(ii)) 
            END DO 
#endif                                                                  
            phi(ii)=( (2-(10.0*one12)*ksq(ii+1))*phi(ii+1)              &
     &             -(1+(1*one12)*ksq(ii+2))*phi(ii+2)  )                &
     &               /(1+(1*one12)*ksq(ii))                             
         ENDIF 
       ENDDO 
                                                                        
       !get phi and dphi at z= c/2 =d-d_aux                             
       !calculate simple derivation of dphi with finite differences     
                                                                        
       !phi(index_c) and                                                
       dphi = -0.5 * ((phi(index_c)-phi(index_c+1))/(xdistance)         &
     &               +((phi(index_c+1)-phi(index_c+2)))/(xdistance)  )  
                                                                        
                                                                        
                                                                        
#ifdef CPP_DEBUG_AP    
!write EmbPot to file                            
!       dphi = 0.5 * ( (phi(1)-phi(2))/(xdistance)                      
!     &               +(phi(2)-phi(3))/(xdistance)  )                   
                                                                        
       OPEN(109,FILE="EmbPot",FORM='formatted',                         &
     &         STATUS='OLD',ACCESS='append')                            
!      WRITE (109,'(A)') "# data Re(E),Im(E):Phi:dPhi:Re(Emb):Im(Emb)"  
       WRITE (109,'(8F10.5)')  gf_Z(en,0),phi(index_c),dphi             &
     &                         ,(dphi/phi(index_c))                     
!       WRITE (109,'(7F11.5)')  gf_Z(en,0),REAL(phi(1)),AIMAG(phi(1))   
!     &                        ,REAL(dphi),AIMAG(dphi)                  
!     &                        ,REAL(dphi/phi(1)),AIMAG(dphi/phi(1))    
       CLOSE(109) 
#endif                                                                  
                                                                        
#ifdef CPP_DEBUG_IP    
!write wavefunction to file                      
      WRITE (filename,'(A,F0.3,A,I0,A)') "wf-",REAL(gf_Z(en,0)),        &
     &          "-",n,".dat"                                            
        IF (n/=1) CLOSE (66) 
        OPEN(66,FILE=filename,                                          &
     &              FORM='formatted',STATUS='REPLACE')                  
                                                                        
       WRITE (66,'(A,I5,A)') "# phi for n=",n,                          &
     &                   "  (z, fft point, Re(phi), Im(phi), V, S)"     
                        !output phi                                     
       DO ii=1,gridpts 
         WRITE (66,'(5F11.5,F7.4)') (imageplane+x(ii)),                 &
     &             (imageplane+x(ii))*(mx3-1)/dt,                       &
     &             phi(ii),V(ii),S                                      
       ENDDO 
      WRITE (filename,'(A,F0.3,A,I0,A)') "wf-",REAL(gf_Z(en,0)),        &
     &          "-",n,"-dphi.dat"                                       
        IF (n/=1) CLOSE (77) 
        OPEN(77,FILE=filename,                                          &
     &              FORM='formatted',STATUS='REPLACE')                  
       WRITE (77,'(A,I5,A)') "# phi for n=",n,                          &
     &                   "  (z, fft point, Re(dphi), Im(dphi), V, S)"   
         WRITE (77,'(5F11.5,F7.4)') (imageplane+x(ii)),                 &
     &             (imageplane+x(ii))*(mx3-1)/dt,                       &
     &             dphi,V(ii),S                                         
!       WRITE (*,*) "d_aux ",d_aux                                      
#endif                                                                  
                                                                        
                                                                        
                            !saves in case d_aux==0 some cpu time       
       IF (d_aux > 0 ) THEN 
        !move trough delta region  & store EmbPoti                      
                ! phi(index_c) , dphi at left side                      
        k_r = sqrt(2.0*gf_Z(en,0)-s) 
        embpot(n,n)= I/2*k_r*                                           &
     &                (                                                 &
     &                         (dphi/(I*k_r)* COS(k_r*d_aux)            &
     &                         +I*phi(index_c)*SIN(k_r*d_aux))          &
     &                  /                                               &
     &                         ( phi(index_c)*COS(k_r*d_aux)            &
     &                        +I*(dphi/(I*k_r))*SIN(k_r*d_aux) )        &
     &                 )                                                
      ELSE 
        embpot(n,n)= 0.5*dphi/phi(index_c) 
      ENDIF 
                                                                        
!      WRITE (0,*) en, REAL(gf_Z(en,0)),embpot(n,n)                     
             ! n                                                        
      ENDDO 
                                                                        
#ifdef CPP_DEBUG_IP    
!close wf. file                                  
      CLOSE (66) 
      CLOSE (77) 
#endif                                                                  
      RETURN 
                                                                        
      END SUBROUTINE 

      !>                                                                
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
      !<-- S: gf_vacuum_EField(en,nk,g2)                                

      SUBROUTINE gf_vacuum_EField(en,nk,embpot) 
!******************************************                             
! calculates numericaly the embedding potential for                     
! linear E-field E*z+V_vac applyed to the system                        
!       01/2007            A. Hanuschkin                                
!******************************************                             
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN) :: en,nk 
      COMPLEX,INTENT(INOUT) ::embpot(:,:) 
      !>                                                                
      !<--Locals                                                        
      INTEGER :: n 
      REAL    :: gk(2),s 
      COMPLEX :: k_r,I 
                                              ! k**2                    
      COMPLEX, ALLOCATABLE :: ksq(:) 
                                              ! derivation of wf        
      COMPLEX :: dphi 
                                              ! wavefunction on mesh    
      COMPLEX , ALLOCATABLE :: phi(:) 
                                              ! save comp. time         
      REAL    :: one12 = (1/12.0) 
      INTEGER :: ii,zz 
      CHARACTER (LEN=20) :: filename 
                                                                        
      I=cmplx(0.0,1.0) 
      !>                                                                
                                                                        
!      if (.not.surface) return !no surface calculation                 
      !CHECK THAT IMAGEBOUNDARY >= D/2                                  
      !CHECK THAT IMAGEPLANE  < IMAGEBOUNDARY                           
      !CHECK imageplane < embPlane @ reading gf_inp                     
                                                                        
                                                                        
      !allocate and initialize                                          
      ALLOCATE(ksq(gridpts)) 
      ALLOCATE(phi(gridpts)); 
                               !init Embedding potential                
      embpot=0.0 
                                                                        
#ifdef CPP_DEBUG   
!not nessecary because it'll be calculated in necessa
                   !but nice if plotted phi(x<c/2)                      
      phi = 0 
      WRITE (*,*) "xdistance,imageboundary,imageplane"                  &
     &            ,xdistance,imageboundary,imageplane                   
#endif                                                                  
                                                                        
       !phi and dphi at z= c/2 =d-d_aux needed                          
       !up to which point integration is needed?                        
!       index_c = Int((d-d_aux)/(xdistance))  !check imageplane < embPla
!       WRITE (*,*) index_c," index_c ",d-d_aux-imageplane              
!     &         ," d-d_aux-imageplane ",d-d_aux," c"                    
                                                                        
!       IF ((index_c.gt.gridpts-2).or.(index_c.le.0)) THEN              
!         WRITE (*,*) index_c                                           
!        WRITE (*,*) "Imageboundary closer than Embedding plane?"       
!          CALL juDFT_error("PROBLEMS WITH SPACING: INDEX_C",calledby="gf_vacuum.F90")
!       ENDIF                                                           
                                                                        
                                                                        
      !CALCULATE PHI and DPHI at c/2                                    
      DO n=1,SIZE(embpot,1) 
      !initial wavefunction                                             
      !inf. in vacuum region V=vac_pot                                  
                                                                        
       gk=MATMUL(coordconv,bk(1:2,nk)+(/k1p(n),k2p(n)/)) 
       s=dot_PRODUCT(gk,gk) 
                                                                        
#ifdef CPP_DEBUG                                                        
       WRITE (oUnit,*) "s ",s," vac_pot ",vac_pot 
#endif                                                                  
                                                                        
       phi(gridpts) = exp(I*sqrt(2.0*(gf_Z(en,0)-V(gridpts) )-s)        &
     &                    *x(gridpts))                                  
       phi(gridpts-1)= exp(I*sqrt(2.0*(gf_Z(en,0)-V(gridpts-1) )-s)     &
     &                    *x(gridpts-1))                                
       ksq(gridpts) = (2.0*(gf_Z(en,0)-V(gridpts))-s)*xdist_sq 
       ksq(gridpts-1)= (2.0*(gf_Z(en,0)-V(gridpts-1))-s)*xdist_sq 
                                                                        
#ifdef CPP_DEBUG_IMAGEPOTI_V0                                           
! WITHOUT POTENTIAL                                                     
       phi(gridpts) = exp(I*sqrt(2.0*(gf_Z(en,0)-V0)-s)                 &
     &                    *x(gridpts))                                  
       phi(gridpts-1) = exp(I*sqrt(2.0*(gf_Z(en,0)-V0)-s)               &
     &                    *x(gridpts-1))                                
      ksq(gridpts) = (2.0*(gf_Z(en,0))-s)*xdist_sq 
      ksq(gridpts-1)= ksq(gridpts) 
#endif                                                                  
                                                                        
                                                                        
       !numerical integration ksq=k**2xh**2                             
                          !set up k^2*h^2 array                         
       DO ii=1,gridpts-2 
        ksq(ii) = ( 2.0*(gf_Z(en,0)-V(ii)) -S)*xdist_sq 
       ENDDO 
                                                                        
                             !calculate phi                             
       DO ii=gridpts-2,1,-1 
         phi(ii)=( (2-(10.0*one12)*ksq(ii+1))*phi(ii+1)                 &
     &             -(1+(1*one12)*ksq(ii+2))*phi(ii+2)  )                &
     &               /(1+(1*one12)*ksq(ii))                             
!         IF (ABS(REAL(phi(ii))).gt.1000) THEN  ! normalize!            
!           phi(ii+1)=phi(ii+1)/ABS(Real(phi(ii)))                      
!           phi(ii+2)=phi(ii+2)/ABS(Real(phi(ii)))  !index_x+2 ii+3?    
                                                                        
                                         ! normalize!                   
         IF (ABS(phi(ii))>1000) THEN 
           phi(ii+1)=phi(ii+1)/ABS(phi(ii)) 
                                             !index_x+2 ii+3?           
           phi(ii+2)=phi(ii+2)/ABS(phi(ii)) 
                                                                        
                                                                        
#ifdef CPP_DEBUG_IP_NA                                                  
                                             ! normalize all!           
            DO zz=gridpts,ii+3,-1 
!           phi(zz)=phi(zz)/ABS(Real(phi(ii)))                          
             phi(zz)=phi(zz)/ABS(phi(ii)) 
            END DO 
#endif                                                                  
            phi(ii)=( (2-(10.0*one12)*ksq(ii+1))*phi(ii+1)              &
     &             -(1+(1*one12)*ksq(ii+2))*phi(ii+2)  )                &
     &               /(1+(1*one12)*ksq(ii))                             
         ENDIF 
       ENDDO 
                                                                        
       !get phi and dphi at z= c/2 =d-d_aux                             
       !calculate simple derivation of dphi with finite differences     
                                                                        
       !phi(index_c) and                                                
       dphi = -0.5 * ((phi(1)-phi(1+1))/(xdistance)                     &
     &               +((phi(1+1)-phi(1+2)))/(xdistance)  )              
                                                                        
                                                                        
                                                                        
#ifdef CPP_DEBUG_AP
    !write EmbPot to file                            
!       dphi = 0.5 * ( (phi(1)-phi(2))/(xdistance)                      
!     &               +(phi(2)-phi(3))/(xdistance)  )                   
                                                                        
       OPEN(109,FILE="EmbPot",FORM='formatted',                         &
     &         STATUS='OLD',ACCESS='append')                            
!      WRITE (109,'(A)') "# data Re(E),Im(E):Phi:dPhi:Re(Emb):Im(Emb)"  
       WRITE (109,'(8F10.5)')  gf_Z(en,0),phi(1),dphi                   &
     &                         ,(dphi/phi(1))                           
!       WRITE (109,'(7F11.5)')  gf_Z(en,0),REAL(phi(1)),AIMAG(phi(1))   
!     &                        ,REAL(dphi),AIMAG(dphi)                  
!     &                        ,REAL(dphi/phi(1)),AIMAG(dphi/phi(1))    
       CLOSE(109) 
#endif                                                                  
                                                                        
#ifdef CPP_DEBUG_IP   
 !write wavefunction to file                      
      WRITE (filename,'(A,F0.3,A,I0,A)') "wf-",REAL(gf_Z(en,0)),        &
     &          "-",n,".dat"                                            
        IF (n/=1) CLOSE (66) 
        OPEN(66,FILE=filename,                                          &
     &              FORM='formatted',STATUS='REPLACE')                  
                                                                        
       WRITE (66,'(A,I5,A)') "# phi for n=",n,                          &
     &                   "  (z, fft point, Re(phi), Im(phi), V, S)"     
                        !output phi                                     
       DO ii=1,gridpts 
         WRITE (66,'(5F11.5,F7.4)')(x(ii)),                            &
     &             (x(ii))*(mx3-1)/dt,                                  &
     &             phi(ii),V(ii),S                                      
       ENDDO 
      WRITE (filename,'(A,F0.3,A,I0,A)') "wf-",REAL(gf_Z(en,0)),        &
     &          "-",n,"-dphi.dat"                                       
        IF (n/=1) CLOSE (77) 
        OPEN(77,FILE=filename,                                          &
     &              FORM='formatted',STATUS='REPLACE')                  
       WRITE (77,'(A,I5,A)') "# phi for n=",n,                          &
     &                   "  (z, fft point, Re(dphi), Im(dphi), V, S)"   
         WRITE (77,'(5F11.5,F7.4)') (x(ii)),                            &
     &             (x(ii))*(mx3-1)/dt,                                  &
     &             dphi,V(ii),S                                         
!       WRITE (*,*) "d_aux ",d_aux                                      
#endif                                                                  
                                                                        
                                                                        
                            !saves in case d_aux==0 some cpu time       
       IF (d_aux > 0 ) THEN 
        !move trough delta region  & store EmbPoti                      
                ! phi(index_c) , dphi at left side                      
        k_r = sqrt(2.0*gf_Z(en,0)-s) 
        embpot(n,n)= I/2*k_r*                                           &
     &                (                                                 &
     &                         (dphi/(I*k_r)* COS(k_r*d_aux)            &
     &                         +I*phi(1)*SIN(k_r*d_aux))                &
     &                  /                                               &
     &                         ( phi(1)*COS(k_r*d_aux)                  &
     &                        +I*(dphi/(I*k_r))*SIN(k_r*d_aux) )        &
     &                 )                                                
      ELSE 
        embpot(n,n)= 0.5*dphi/phi(1) 
      ENDIF 
                                                                        
!      WRITE (0,*) en, REAL(gf_Z(en,0)),embpot(n,n)                     
             ! n                                                        
      ENDDO 
                                                                        
#ifdef CPP_DEBUG_IP    
!close wf. file                                  
      CLOSE (66) 
      CLOSE (77) 
#endif                                                                  
      RETURN 
                                                                        
      END SUBROUTINE 

      !>

                                                                        
      !<-- S: gf_vacuum_init(jspin,cell,kpts,kp,gfinp,vac_pot)          
      SUBROUTINE gf_vacuum_init(jspin,cell,kpts,lapw,lapw_gf,gfinp,     &
     &           vac_pot_in,stars,l_noco_in)
!******************************************                             
!     Initializes the private variables of the module                   
!                          D. Wortmann                                  
!******************************************                             
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)             :: jspin 
      TYPE(t_cell),INTENT(IN)        :: cell 
      TYPE(t_kpts),INTENT(IN)        :: kpts 
      TYPE(t_lapw),INTENT(IN)        :: lapw
      TYPE(t_lapw_gf),INTENT(IN),TARGET :: lapw_gf
      TYPE(t_embinp),INTENT(IN)       :: gfinp 
      TYPE(t_stars),INTENT(IN)       :: stars 
      REAL   ,INTENT(IN)             :: vac_pot_in
      logical,intent(in)             :: l_noco_in
      !>                                                                
      coordconv = cell%bmat(1:2,1:2) 
      allocate(bk(size(kpts%bk,1),size(kpts%bk,2)))
      bk        = kpts%bk 
      k1p       => lapw_gf%k1p(:,jspin) 
      k2p       => lapw_gf%k2p(:,jspin) 
                                 ! d/2                                  
      d         = cell%z1 
                                     ! dt/2                             
      dt        = cell%amat(3,3)/2.0 
      d_aux     = dt-d 
                                 !fft box                               
      mx3       = stars%mx3 
      surface   = gfinp%l_surface 
      imagepoti = gfinp%l_imagepoti 
      imageplane= gfinp%imageplane 
      imageboundary= gfinp%imageboundary 
      Efield    = gfinp%Efield 
      vac_pot   = vac_pot_in 
      l_noco=l_noco_in
      END SUBROUTINE 
      !>                                                                
                                                                        
                                                                        
      !<-- S: gf_vacuum_Poti_init(                                      
      SUBROUTINE gf_vacuum_Poti_init()
      IMPLICIT NONE
!******************************************                             
!     Setup of potential if EmbPot is calculated numericaly             
!         01/24/2007       A. Hanuschkin                                
!******************************************                             

      IF (imagepoti) THEN 
       CALL gf_vacuum_image_init() 
       RETURN 
      ENDIF 
      IF (EField/=0.0) THEN 
       CALL gf_vacuum_EField_init() 
       RETURN 
      ENDIF 
      END SUBROUTINE 
      !>                                                                
                                                                        

      subroutine gf_vacuum_check_vacpot(vpw,stars,c,vac_pot)
      implicit none
      type(t_stars),intent(in)::stars
      complex,intent(in)::vpw(:,:)
      real,intent(in)::c
      real,intent(inout)::vac_pot

      integer::k,n,in
      complex::v(2*stars%mx3+1,size(vpw,2))
      v=0.0
      do k=-stars%mx3,stars%mx3
      	in=stars%ig(0,0,k)
      	if (in==0) cycle
      	if (k<0) then
      	    n=2*stars%mx3+k+2
      	else
      	    n=k+1
      	endif
      	v(n,:)=vpw(in,:)
      enddo
      v(:,1)=fft(v(:,1),inv=.true.)
      if (size(v,2)>1) v(:,2)=fft(v(:,2),inv=.true.)






      write(oUnit,*) "Vacuum potential calculation"

	  v(1,1)=maxval(real(v(:,1)))
	  if (size(v,2)>1) v(1,2)=maxval(real(v(:,2)))

      write(oUnit,*) "Potential:",real(v(1,:))
      k=min(size(v,2),2)
      vac_pot=sum(real(v(1,:k)))/k
      write(oUnit,*) "Vacpot:",vac_pot

      end subroutine
                                                                        
      !<-- S: gf_vacuum_EField_image_init(                              
      SUBROUTINE gf_vacuum_EField_init() 
      IMPLICIT NONE
!******************************************                             
!     Initializes the E*z+V_vac potential                               
!         01/24/2007       A. Hanuschkin                                
!******************************************                             
      INTEGER :: cutoff,ii 
                                                                        
                                                                        
      !calculate gridspacing                                            
      xdistance = REAL(100)/REAL(gridpts) 
      xdist_sq = xdistance*xdistance 
                                                                        
      !set up potential mesh                                            
!      WRITE (oUnit,*) ''                                                   
!      WRITE (oUnit,*) '#EFie potential meshpoint, V, z (a.u.),z (fftbox u)'
      DO ii=1,gridpts 
         x(ii)=REAL(ii)*xdistance; 
         V(ii) = vac_pot+(EField*x(ii)) 
!        WRITE (oUnit,'(I6,3F11.5)') ii,V(ii),(x(ii))                       
!     &         ,(x(ii))*(mx3-1)/dt  !to out file                       
      END DO 
!      WRITE (oUnit,*) '#EField potential END'                              
!      WRITE (oUnit,*) ''                                                   
                                                                        
                                                                        
#ifdef CPP_DEBUG_IMAGEPOTI_V0                                           
!    SET POTENTIAL CONSTANT FOR TESTING                                 
      V0=vac_pot 
      V =vac_pot 
      WRITE (oUnit,*) "SET POTENTIAL CONSTANT=V0=",V0," FOR TESTING" 
#endif                                                                  
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
                                                                        
                                                                        
                                                                        
      !<-- S: gf_vacuum_image_init(                                     
      SUBROUTINE gf_vacuum_image_init() 
      IMPLICIT NONE
!******************************************                             
!     Initializes the 1/z potential                                     
!                          A. Hanuschkin                                
!******************************************                             
      INTEGER :: cutoff,ii 
      REAL :: embplane 
                                                                        
                                                                        
      !calculate gridspacing                                            
      xdistance = REAL(imageboundary-imageplane)/REAL(gridpts) 
      xdist_sq = xdistance*xdistance 
                                                                        
      !set up potential mesh                                            
      IF (EField==0.0) THEN 
       DO ii=1,gridpts 
          x(ii)=REAL(ii)*xdistance; 
          V(ii) = vac_pot-(1/(4.0*x(ii))); 
!         V(1)=vac_pot                                                  
!         V(2)=vac_pot                                                  
       END DO 
      ELSE 
        embplane = d-d_aux 
        DO ii=1,gridpts 
          x(ii)=REAL(ii)*xdistance; 
          V(ii) = vac_pot+EField*(x(ii)+imageplane-embplane)            &
     &            -(1/(4.0*x(ii)))                                      
       END DO 
      ENDIF 
                                                                        
                                                                        
      !set up potential mesh                                            
!      WRITE (oUnit,*) ''                                                   
!      WRITE (oUnit,*) '#IMAGE potential meshpoint, V, z (a.u.),z (fftbox u)
!       DO ii=1,gridpts                                                 
!         WRITE (oUnit,'(I6,3F11.5)') ii,V(ii),(imageplane+x(ii))           
!     &          ,(imageplane+x(ii))*(mx3-1)/dt  !to out file           
!      END DO                                                           
!      WRITE (oUnit,*) '#IMAGE potential END'                               
!      WRITE (oUnit,*) ''                                                   
                                                                        
                                                                        
                                                                        
      !set up mixed potential mesh                                      
!     WRITE (oUnit,*) ''                                                    
!     WRITE (oUnit,*) '#MIXEDIM. potential meshpoint, V, z, fftbox'         
      cutoff = INT(0.16/xdistance) 
      IF (cutoff>gridpts-2) THEN 

          CALL juDFT_error("CUTOFF CHOOSEN TO LARGE !",calledby="gf_vacuum.F90")
      ENDIF 
      DO ii=1,cutoff 
!         V(ii) = (1-IBETA(z))*V_LDA(ii)+IBETA(z)*V(ii);                
          V(ii) = V(cutoff+1) 
!         WRITE (oUnit,'(I6,3F11.5)') ii,V(ii),(imageplane+x(ii))           
!    &         ,(imageplane+x(ii))*(mx3-1)/dt  !to out file             
      END DO 
      DO ii=cutoff,gridpts 
!!        V(ii) = (1-IBETA(z))*V_LDA(ii)+IBETA(z)*V(ii);                
      END DO 
                                                                        
#ifdef CPP_DEBUG_IMAGEPOTI_V0                                           
!    SET POTENTIAL CONSTANT FOR TESTING                                 
      V0=vac_pot 
      V =vac_pot 
      WRITE (oUnit,*) "SET IMAGEPOTENTIAL CONSTANT=V0=",V0," FOR TESTING" 
#endif                                                                  
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
                                                                        
                                                                        
      !<-- S: gf_vacuum(gfinp,lapw,kp(:,jsp),en,nk,jsp)                 
      SUBROUTINE gf_vacuum(gfinp,lapw,lapw_gf,en,nk,jsp)
!******************************************                             
!     Calculate Embedding potential with a potential in Vacuum which    
!     might be different from zero. See Eq. 5.50 from my Thesis!        
!                                                                       
!     (Very ugly hack!!)                                                
!     Two additional parameters are needed:                             
!     The thickness of the auxillary volume ->stored in gfinp%dp1       
!     The value of the vacuum potential ->stored in gfinp%dp2           
!                                                                       
!                          D. Wortmann                                  
!******************************************                             
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_embinp),INTENT(IN) :: gfinp
      TYPE(t_lapw),INTENT(IN)  :: lapw
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      INTEGER,INTENT(IN)       :: en,nk,jsp 
      !>                                                                
      !<--Locals                                                        
      COMPLEX               :: embpot(lapw_gf%nv2d,lapw_gf%nv2d) 
      INTEGER               :: n 
      COMPLEX               :: kl,kr 
      COMPLEX               :: I 
      REAL                  :: d 
      I=CMPLX(0.0,1.0) 
      !>                                                                
      d=gfinp%dp1 
      embpot=0.0 
      DO n=1,lapw_gf%nv2(jsp) 
         kl=SQRT(2*gf_Z(en,0)-lapw_gf%rkp(n,jsp)**2) 
         kr=SQRT(2*(gf_Z(en,0)-gfinp%dp2)-lapw_gf%rkp(n,jsp)**2) 
         embpot(n,n)=I/2*kl*(kr*COS(kl*d)+I*kl*SIN(kl*d))/(kl*COS(kl*d) &
     &        +I*kr*SIN(kl*d))                                          
      ENDDO 
      CALL gf_write2dmat(IO2D_EMB,2,2,en,nk,jsp,lapw_gf,embpot) 
                                                                        
      END SUBROUTINE 
      !>                                                                
#endif
      END                                           
