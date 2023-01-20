!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_force_sf
   USE m_constants
   USE m_types
   !-----------------------------------------------------------------------------
   ! This routine calculates a contribution to the forces stemming
   ! from the discontinuity of density and potential at the muffin-tin
   ! boundary oint n [ rho V (IS) - rho V (MT) ] dS
   ! Klueppelberg May 13
   !
   ! The surface integral is found in the first two lines of equation (14) in 
   ! 
   ! Klüppelberg et al., PRB 91 035105 (2015).
   ! 
   ! Klueppelberg (force level 3) 
   !-----------------------------------------------------------------------------

   IMPLICIT NONE

   COMPLEX, PRIVATE, SAVE, ALLOCATABLE :: force_mt(:,:)
   COMPLEX, PRIVATE, SAVE, ALLOCATABLE :: force_is(:,:)
   LOGICAL, PRIVATE, SAVE :: isdone=.false.,mtdone=.false.

CONTAINS
   SUBROUTINE force_sf_is(atoms_in,stars,sym,jsp,cell,qpw,vpw,excpw,vxcpw)
      !--------------------------------------------------------------------------
      ! This subroutine calculates oint n rho V dS evaluated with quantities 
      ! from the interstital. The Fourier transform of density and potential is 
      ! first calculated as an expansion in spherical harmonics. Then the normal 
      ! vector, which is proportional to Y_1t, connects the l component of rho
      ! with the l+-1 components of V. This is done up to a cutoff lmax.
      ! It is called in a spin loop at the end of vgen.F90-
      !--------------------------------------------------------------------------

      USE m_sphbes
      USE m_phasy1
      USE m_gaunt

      TYPE(t_sym),INTENT(IN)     :: sym
      TYPE(t_stars),INTENT(IN)   :: stars
      TYPE(t_cell),INTENT(IN)    :: cell
      TYPE(t_atoms),INTENT(IN)   :: atoms_in

      INTEGER, INTENT (IN) :: jsp

      COMPLEX, INTENT (IN) :: qpw(:,:) !(stars%ng3,input%jspins)
      COMPLEX, INTENT (IN) :: vpw(:,:) !(n3d,jspd)
      COMPLEX, INTENT (IN) :: excpw(stars%ng3)
      COMPLEX, INTENT (IN) :: vxcpw(:,:) !(stars%ng3,input%jspins)

      TYPE(t_atoms)              :: atoms

      INTEGER :: n,j,itype,s,l ,lm,t,lp,mp,lmp,jp,natom,m
      REAL    :: r,r2
      COMPLEX :: img,rhoprep,Vprep
      LOGICAL :: isthere

      INTEGER :: lmaxb(atoms_in%ntype)
      COMPLEX :: coeff(3,-1:1),qpw2(stars%ng3,size(qpw,2)),qpwcalc(stars%ng3,size(qpw,2))
      REAL   , ALLOCATABLE :: bsl(:,:,:)
      COMPLEX, ALLOCATABLE :: pylm(:,:,:),rho(:),V(:),pylm2(:,:)

      CALL timestart("Force level 3 (IS)")

      atoms=atoms_in
      atoms%lmax = 2*atoms_in%lmax
      atoms%lmaxd = 2*atoms_in%lmaxd
      lmaxb = atoms%lmax
      img = cmplx(0.0,1.0)

      CALL init_sf(sym,cell,atoms)
      isdone = .true.

      ALLOCATE ( bsl(0:atoms%lmaxd,stars%ng3,atoms%ntype) )

      ALLOCATE ( pylm2((atoms%lmaxd+1)**2,atoms%ntype ))
      ALLOCATE ( rho((atoms%lmaxd+1)**2),V((atoms%lmaxd+1)**2) )

      coeff(:, :) =   cmplx(0.0,0.0)
      coeff(1,-1) =     sqrt(tpi_const/3.)
      coeff(1, 1) =    -sqrt(tpi_const/3.)
      coeff(2,-1) = img*sqrt(tpi_const/3.)
      coeff(2, 1) = img*sqrt(tpi_const/3.)
      coeff(3, 0) =  sqrt(2.*tpi_const/3.)

      qpwcalc = qpw

      DO itype = 1,atoms%ntype
         CALL sphbes(atoms%lmax(itype),0.0,bsl(:,1,itype))
         DO s = 2,stars%ng3
            ! Only call sphbes if the length of the star changed.
            IF (abs(stars%sk3(s)-stars%sk3(s-1)).gt.1.0e-14) THEN
               r = stars%sk3(s)*atoms%rmt(itype)
               CALL sphbes(atoms%lmax(itype),r,bsl(:,s,itype))
            ELSE
               bsl(:,s,itype) = bsl(:,s-1,itype)
            END IF
         END DO ! s
      END DO ! itype


      force_is = 0.0

      DO itype = 1,atoms%ntype
         natom = atoms%firstAtom(itype)
         r2  = atoms%rmt(itype)**2
         rho = 0.0
         V   = 0.0

         DO s = 1,stars%ng3
            CALL phasy1(atoms,stars,sym,cell,s,pylm2(:,:))

            DO l = 0,atoms%lmax(itype)-1
               rhoprep = stars%nstr(s) * bsl(l,s,itype) * qpwcalc(s,jsp)
               IF (l.eq.0) THEN
                  Vprep = stars%nstr(s) * bsl(l,s,itype) * &
                                          (vpw(s,jsp)-vxcpw(s,jsp)+excpw(s)) 
                                          ! Switching between Veff and VCoul + exc
               END IF
               ! For l = 0 we calculate rho_00 and V_00:
               DO m = -l,l
                  lm = l*(l+1) + m + 1
                  rho(lm) = rho(lm) + rhoprep * pylm2(lm,itype)!pylm(lm,s,itype)!
                  IF (l.gt.0) CYCLE
                  V(lm) =   V(lm) +   Vprep * pylm2(lm,itype)!pylm(lm,s,itype)!
               END DO ! m
               ! And V_1mp, for l > 0, we calculate rho_lm, V_l+1,mp:
               Vprep = stars%nstr(s) * bsl(l+1,s,itype) * &
                                       (vpw(s,jsp)-vxcpw(s,jsp)+excpw(s))
               DO m = -l-1,l+1
                  lm = (l+1)*(l+2) + m + 1
                  V(lm) =   V(lm) +   Vprep * pylm2(lm,itype)!pylm(lm,s,itype)!
               END DO ! m
            END DO ! l
         END DO ! s

         DO l = 0,atoms%lmax(itype)-1 ! new: altered s and l loop above
            DO m = -l,l
               lm = l*(l+1) + m + 1
               ! Because rho_lm occurs with V_l-1,mp and V_l+1,mp:
               DO lp = abs(l-1),l+1,2
                  DO t = -1,1
                     mp = t-m
                     IF (lp.lt.abs(mp)) CYCLE
                     lmp = lp*(lp+1) + mp + 1
                     force_is(:,itype) = force_is(:,itype) + r2 * rho(lm) * &
                                         V(lmp) * conjg(coeff(:,t)) * &
                                         gaunt2(1,l,lp,t,m,mp,atoms%lmax(itype))
                  END DO ! t
               END DO ! lp
            END DO ! m
         END DO ! l
      END DO ! itype

      DEALLOCATE ( bsl,rho,V )
      DEALLOCATE ( pylm2 )

      CALL timestop("Force level 3 (IS)")

   END SUBROUTINE force_sf_is

   SUBROUTINE force_sf_mt(atoms,sphhar,jspin,ispin,fmpi,vr,excr,vxcr,rho,sym,cell )
      !--------------------------------------------------------------------------
      ! This subroutine calculates the contribution evaluated with quantities
      ! from the Muffin Tins.
      ! 
      ! n rho V = sum(nu,nup) rho(nu)V(nup)
      !           * sum(m,mu,mup) c_1m* Y_1m* c_lnu Y_numu c_lnup Y_nupmup
      !
      ! It is called in a spin loop at the end of cdngen.F90
      !--------------------------------------------------------------------------

      USE m_gaunt
      USE m_ylm

      TYPE(t_mpi),INTENT(IN)   :: fmpi
      TYPE(t_sym),INTENT(IN)   :: sym
      TYPE(t_cell),INTENT(IN)  :: cell
      TYPE(t_sphhar),INTENT(IN):: sphhar
      TYPE(t_atoms),INTENT(IN) :: atoms

      INTEGER, INTENT (IN) :: jspin
      INTEGER, INTENT (IN) :: ispin 

      REAL   , INTENT (IN) :: vr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype) ! 
      REAL   , INTENT (IN) :: rho(:,0:,:,:)!(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins)
      REAL   , INTENT (IN) :: excr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype)
      REAL   , INTENT (IN) :: vxcr(:,0:,:,:)!(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins)

      INTEGER :: natom,itype,nd,lh,l,lhp,lp,mem,m,memp,mp,t,i,lmp
      REAL    :: pot,den
      COMPLEX :: img,factor

      COMPLEX :: coeff(3,-1:1)
      COMPLEX :: d1((atoms%lmaxd+1)**2,atoms%ntype ),d2((atoms%lmaxd+1)**2,atoms%ntype )

      CALL timestart("Force level 3 (MT)")

      CALL init_sf(sym,cell,atoms)
      mtdone = .true.

      img = cmplx(0.0,1.0)
      force_mt = 0.0

      coeff(:, :) =   cmplx(0.0,0.0)
      coeff(1,-1) =     sqrt(tpi_const/3.)
      coeff(1, 1) =    -sqrt(tpi_const/3.)
      coeff(2,-1) = img*sqrt(tpi_const/3.)
      coeff(2, 1) = img*sqrt(tpi_const/3.)
      coeff(3, 0) =  sqrt(2.*tpi_const/3.)

      d1 = 0
      d2 = 0

      ! Calculate forces: For each atom, loop over all lattice harmonics.
      DO itype = 1,atoms%ntype
         natom = atoms%firstAtom(itype)
         nd = sym%ntypsy(natom)

         DO lh = 0,sphhar%nlh(nd)
            l = sphhar%llh(lh,nd)

            ! The l=0 component of the potential array is saved with an additional 
            ! factor r/sfp. For this calculation, we need the pure potential.
            pot = vr(atoms%jri(itype),lh,itype)
            IF (lh.eq.0) THEN
               pot = pot*sfp_const/atoms%rmt(itype)
            END IF

            pot = excr(atoms%jri(itype),lh,itype)-vxcr(atoms%jri(itype),lh,itype,ispin)+pot

            DO mem = 1,sphhar%nmem(lh,nd)
               m = sphhar%mlh(mem,lh,nd)
               lmp = l*(l+1) + m + 1
               d1(lmp,itype) = d1(lmp,itype) + sphhar%clnu(mem,lh,nd) * &
                               rho(atoms%jri(itype),lh,itype,ispin)/atoms%rmt(itype)**2
               d2(lmp,itype) = d2(lmp,itype) + sphhar%clnu(mem,lh,nd) * pot
            END DO ! mem

            DO lhp = 0,sphhar%nlh(nd)
               lp = sphhar%llh(lhp,nd)
               IF (abs(l-lp).ne.1) CYCLE

               den = rho(atoms%jri(itype),lhp,itype,ispin)

               DO mem = 1,sphhar%nmem(lh,nd)
                  m = sphhar%mlh(mem,lh,nd)
                  DO memp = 1,sphhar%nmem(lhp,nd)
                     mp = sphhar%mlh(memp,lhp,nd)
                     IF (abs(m+mp).gt.1) CYCLE

                     ! Due to the normal vector n, the lattice harmonics form a
                     ! Gaunt coefficient with the Y_1m.
                     factor = pot * den * sphhar%clnu(mem,lh,nd) * sphhar%clnu(memp,lhp,nd)&
                                  * gaunt1(1,l,lp,m+mp,m,mp,atoms%lmaxd)
               
                     force_mt(:,itype) = force_mt(:,itype) + factor * conjg(coeff(:,m+mp))

                  END DO ! memp
               END DO ! mem

            END DO ! lhp

         END DO ! lh
      END DO ! itype

      CALL timestop("Force level 3 (MT)")

   END SUBROUTINE force_sf_mt

   SUBROUTINE init_sf(sym,cell,atoms)
      ! Initialize force arrays if neither force_sf_is nor force_sf_mt were
      ! executed yet.
      ! Called at the beginning of cdnval.F90 to once fill the wigner array
      ! Also called in force_sf_is/mt to guarantee that all arrays are allocated.

      TYPE(t_sym),INTENT(IN)   :: sym
      TYPE(t_cell),INTENT(IN)  :: cell
      TYPE(t_atoms),INTENT(IN) :: atoms

      IF (isdone.OR.mtdone) RETURN
      IF (.NOT.ALLOCATED(force_is)) THEN
         ALLOCATE ( force_is(3,atoms%ntype),force_mt(3,atoms%ntype) )
      END IF
      force_is = 0.0
      force_mt = 0.0

   END SUBROUTINE init_sf

   SUBROUTINE exit_sf(isp,atoms,force)
      ! Write out the force contribution from the surface terms and deallocate
      ! arrays if all force_sf_is and force_sf_mt were executed.
      ! Called at the end of totale.f90.

      INTEGER,INTENT(IN)         :: isp
      TYPE(t_atoms),INTENT(IN)   :: atoms
      REAL,INTENT(INOUT)         :: force(:,:,:)

      INTEGER :: itype,dir
      COMPLEX :: force_sf(3,atoms%ntype)

      IF (isdone.AND.mtdone) THEN
         force_sf(:,:) = force_is(:,:) - force_mt(:,:)
         force(:,:,isp) = force(:,:,isp) + real(force_sf(:,:))
         WRITE (oUnit,*)
         DO itype = 1,atoms%ntype
            WRITE (oUnit,FMT=8010) itype
            WRITE (oUnit,FMT=8020) (force_sf(dir,itype),dir=1,3)
         END DO ! itype
         isdone = .false.
         mtdone = .false.
         DEALLOCATE ( force_is,force_mt )
      END IF

8010  FORMAT (' FORCES: SURFACE CORRECTION FOR ATOM TYPE',i4)
8020  FORMAT (' FX_SF=',2f10.6,' FY_SF=',2f10.6,' FZ_SF=',2f10.6)

   END SUBROUTINE exit_sf

END MODULE m_force_sf
