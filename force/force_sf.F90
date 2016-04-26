      MODULE m_force_sf
!     *****************************************************************
!     This routine calculates a contribution to the forces stemming
!     from the discontinuity of density and potential at the muffin-tin
!     boundary oint n [ rho V (IS) - rho V (MT) ] dS
!     Klueppelberg May 13
!     *****************************************************************

!     To enable debug code that compares potential and density on the
!     muffin-tin boundary, uncomment the following line:
!#define debug

      IMPLICIT NONE

      COMPLEX, PRIVATE, SAVE, ALLOCATABLE :: force_mt(:,:)
      COMPLEX, PRIVATE, SAVE, ALLOCATABLE :: force_is(:,:)
      LOGICAL, PRIVATE, SAVE :: isdone=.false.,mtdone=.false.

      CONTAINS

      SUBROUTINE force_sf_is(atoms_in,stars,sym,jsp,cell,qpw,vpw,excpw,vxcpw )
!     *****************************************************************
!     This subroutine calculates the contribution evaluated with
!     quantities from the interstital oint n rho V dS
!     The Fourier transform of density and potential is first calculated
!     as an expansion in spherical harmonics. Then the normal vector,
!     which is proportional to Y_1t connects the l component of rho
!     with the l+-1 components of V. This is done up to a cutoff lmax.
!     It is called in a spin loop at the end of vgen.F
!     *****************************************************************

      USE m_constants, ONLY : tpi_const
      USE m_sphbes
      USE m_phasy1
      USE m_gaunt
#ifdef debug
      USE m_ylm
#endif
    USE m_types
      IMPLICIT NONE
      TYPE(t_sym),INTENT(IN)     :: sym
      TYPE(t_stars),INTENT(IN)   :: stars
      TYPE(t_cell),INTENT(IN)    :: cell
      TYPE(t_atoms),INTENT(IN)   :: atoms_in
      TYPE(t_atoms)              :: atoms !copy with modified data
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jsp

!     .. Array Arguments ..
      COMPLEX, INTENT (IN) :: qpw(:,:) !(stars%n3d,dimension%jspd)
      COMPLEX, INTENT (IN) :: vpw(:,:)!(n3d,jspd)
      COMPLEX, INTENT (IN) :: excpw(stars%n3d)
      COMPLEX, INTENT (IN) :: vxcpw(:,:)!(stars%n3d,dimension%jspd)

!     .. Local Scalars ..
      INTEGER :: n,j,itype,s,l ,lm,t,lp,mp,lmp,jp,natom,m
      REAL    :: r,r2
      COMPLEX :: img,rhoprep,Vprep
      LOGICAL :: isthere

!     .. Local Arrays ..
      INTEGER :: lmaxb(atoms_in%ntypd)
      COMPLEX :: coeff(3,-1:1),qpw2(stars%n3d,size(qpw,2)),qpwcalc(stars%n3d,size(qpw,2))
      REAL   , ALLOCATABLE :: bsl(:,:,:)
      COMPLEX, ALLOCATABLE :: pylm(:,:,:),rho(:),V(:),pylm2(:,:)
!       COMPLEX, ALLOCATABLE :: qpw2(:,:),qpwcalc(:,:)
#ifdef debug
      REAL    :: vec(3)
      COMPLEX :: factorrho,factorv
      COMPLEX, ALLOCATABLE :: ylm(:),testrho(:,:),testV(:,:)
#endif
      atoms=atoms_in
      atoms%lmax = 2*atoms_in%lmaxd!60!
      lmaxb = atoms%lmax
      img = cmplx(0.0,1.0)

      CALL init_sf(sym,cell,atoms)
      isdone = .true.

      ALLOCATE ( bsl(stars%n3d,0:atoms%lmaxd,atoms%ntypd) )

      ALLOCATE ( pylm2((atoms%lmaxd+1)**2,atoms%ntype ))
      ALLOCATE ( rho((atoms%lmaxd+1)**2),V((atoms%lmaxd+1)**2) )

#ifdef debug
      ALLOCATE ( ylm((atoms%lmaxd+1)**2),testrho((atoms%lmaxd+1)**2,atoms%ntypd) )
      ALLOCATE ( testV((atoms%lmaxd+1)**2,atoms%ntypd) )
#endif

      coeff(:, :) =   cmplx(0.0,0.0)
      coeff(1,-1) =     sqrt(tpi_const/3.)
      coeff(1, 1) =    -sqrt(tpi_const/3.)
      coeff(2,-1) = img*sqrt(tpi_const/3.)
      coeff(2, 1) = img*sqrt(tpi_const/3.)
      coeff(3, 0) =  sqrt(2.*tpi_const/3.)
      WRITE (1704,*) 'excpw:',excpw
      WRITE (1704,*) 'vxcpw:',vxcpw
!     load in density without coretails
      INQUIRE (FILE='qpw',EXIST=isthere)
      IF (isthere.and..false.) THEN
        qpw2 = 0.0
        OPEN (15,file='qpw',form='formatted',status='unknown')
        DO jp = 1,size(qpw,2)
        DO s = 1,stars%n3d
          READ (15,'(i2,i10,2f20.14)') n,j,qpw2(s,jp)
        END DO ! s
        END DO ! jp
        CLOSE (15)
        IF (any(abs(qpw2).gt.10**(-6))) THEN
          qpwcalc = qpw2
        ELSE
          qpwcalc = qpw
        END IF
      ELSE
        qpwcalc = qpw
      END IF

!     prepare star quantities
!       DO s = 1,ng3
!         DO itype = 1,ntype
!           r = sk3(s)*rmt(itype)
!           CALL sphbes(lmax,r,bsl(s,:,itype))
!         END DO ! itype
!       END DO ! s
      DO itype = 1,atoms%ntype
        CALL sphbes(atoms%lmax(itype),0.0,bsl(1,:,itype))
        DO s = 2,stars%ng3
!         Only call sphbes if the length of the star changed
          IF (abs(stars%sk3(s)-stars%sk3(s-1)).gt.1.0e-14) THEN
            r = stars%sk3(s)*atoms%rmt(itype)
            CALL sphbes(atoms%lmax(itype),r,bsl(s,:,itype))
          ELSE
            bsl(s,:,itype) = bsl(s-1,:,itype)
          END IF
        END DO ! s
      END DO ! itype


      force_is = 0.0

      natom = 1
      DO itype = 1,atoms%ntype
        r2  = atoms%rmt(itype)**2
        rho = 0.0
        V   = 0.0
!         DO l = 0,lmax-1
!           DO s = 1,ng3
! !           calculate phase factors for the current atom type to prevent
! !           total overhead. it is still calculated lmax times even 
! !           though phasy1 generates pylm for all values of l at once.
! !           but this loop sequence is more convenient.
! !           allocating pylm as pylm(lm,s,itype) and precalculating it
! !           leads to exhaustive use of memory in larger systems
! !             CALL phasy1(
! ! !      >              ntypd,n3d,natd,nop,lmaxd,ntype,neq,lmax,
! !      >              1,n3d,1,nop,lmax,1,neq(itype),lmaxb(itype),0,
! !      >              2.*tpi,taual(1:3,natom),bmat,kv3,!1+sum(neq(1:itype-1))),bmat,kv3,
! !      >              tau,mrot,symor,s,invtab,
! !      <              pylm2(:))
!             pylm2 = 1.
! !             IF ((s.eq.1).and.(l.eq.0)) WRITE (851,*) pylm2(:)
!             rhoprep = nstr(s) * bsl(s,l,itype) * qpwcalc(s,jsp)
!             IF (l.eq.0) THEN
!               Vprep = nstr(s) * bsl(s,l,itype)
!      *              * (1*vpw(s,jsp)-1*vxcpw(s,jsp)+1*excpw(s))
! !      *              * vpw(s,jsp)
!             END IF
! !           for l = 0 we calculate rho_00 and V_00...
!             DO m = -l,l
!               lm = l*(l+1) + m + 1
!               rho(lm) = rho(lm) + rhoprep * pylm2(lm)!pylm(lm,s,itype)!
!               IF (l.gt.0) CYCLE
!                 V(lm) =   V(lm) +   Vprep * pylm2(lm)!pylm(lm,s,itype)!
!             END DO ! m
! !           ... and V_1mp, for l > 0, we calculate rho_lm, V_l+1,mp ...
!               Vprep = nstr(s) * bsl(s,l+1,itype)
!      *              * (1*vpw(s,jsp)-1*vxcpw(s,jsp)+1*excpw(s))
! !      *              * vpw(s,jsp)
!             DO m = -l-1,l+1
!               lm = (l+1)*(l+2) + m + 1
!                 V(lm) =   V(lm) +   Vprep * pylm2(lm)!pylm(lm,s,itype)!
!             END DO ! m
!           END DO ! s

        DO s = 1,stars%ng3 !l = 0,atoms%lmax-1
!         calculate phase factors for the current atom type to prevent
!         total overhead. it is still calculated lmax times even 
!         though phasy1 generates pylm for all values of l at once.
!         but this loop sequence is more convenient.
!         allocating pylm as pylm(lm,s,itype) and precalculating it
!         leads to exhaustive use of memory in larger systems
          CALL phasy1(atoms,stars,sym,cell,s,pylm2(:,:))

          DO l = 0,atoms%lmax(itype)-1 !s = 1,stars%ng3
!           calculate phase factors for the current atom type to prevent
!           total overhead. it is still calculated lmax times even 
!           though phasy1 generates pylm for all values of l at once.
!           but this loop sequence is more convenient.
!           allocating pylm as pylm(lm,s,itype) and precalculating it
!           leads to exhaustive use of memory in larger systems
!             CALL phasy1(
! !      >              ntypd,n3d,natd,nop,lmaxd,ntype,neq,lmax,
!      >              1,n3d,1,nop,lmax,1,neq(itype),lmaxb(itype),0,
!      >              2.*tpi,taual(1:3,natom),bmat,kv3,!1+sum(neq(1:itype-1))),bmat,kv3,
!      >              tau,mrot,symor,s,invtab,
!      <              pylm2(:))
!             IF ((s.eq.1).and.(l.eq.0)) WRITE (851,*) pylm2(:)
            rhoprep = stars%nstr(s) * bsl(s,l,itype) * qpwcalc(s,jsp)
            IF (l.eq.0) THEN
              Vprep = stars%nstr(s) * bsl(s,l,itype) * (1*vpw(s,jsp)-1*vxcpw(s,jsp)+1*excpw(s)) ! Switching between Veff and VCoul + exc
!      *              * vpw(s,jsp)
            END IF
!           for l = 0 we calculate rho_00 and V_00...
            DO m = -l,l
              lm = l*(l+1) + m + 1
              rho(lm) = rho(lm) + rhoprep * pylm2(lm,itype)!pylm(lm,s,itype)!
              IF (l.gt.0) CYCLE
                V(lm) =   V(lm) +   Vprep * pylm2(lm,itype)!pylm(lm,s,itype)!
            END DO ! m
!           ... and V_1mp, for l > 0, we calculate rho_lm, V_l+1,mp ...
              Vprep = stars%nstr(s) * bsl(s,l+1,itype) * (1*vpw(s,jsp)-1*vxcpw(s,jsp)+1*excpw(s))
!      *              * vpw(s,jsp)
            DO m = -l-1,l+1
              lm = (l+1)*(l+2) + m + 1
                V(lm) =   V(lm) +   Vprep * pylm2(lm,itype)!pylm(lm,s,itype)!
            END DO ! m
          END DO ! l
          END DO ! s

!           V = 0.0
!           V(1) = sqrt(2.*tpi)
!           rho = 0.0
!           rho(1) = sqrt(2.*tpi)

        DO l = 0,atoms%lmax(itype)-1 ! new: altered s and l loop above

          DO m = -l,l
            lm = l*(l+1) + m + 1
            WRITE (1705,*) itype,l,m,lm-1,rho(lm),V(lm)
!           ... because rho_lm occurs with V_l-1,mp and V_l+1,mp
            DO lp = abs(l-1),l+1,2
              DO t = -1,1
                mp = t-m
                IF (lp.lt.abs(mp)) CYCLE
                lmp = lp*(lp+1) + mp + 1
                force_is(:,itype) = force_is(:,itype) + r2&
                     * rho(lm) * V(lmp) * conjg(coeff(:,t)) * gaunt1(1,l,lp,t,m,mp,atoms%lmax(itype))
              END DO ! t
            END DO ! lp
          END DO ! m
        END DO ! l
#ifdef debug
        testrho(:,itype) = rho(:)
          testV(:,itype) =   V(:)
#endif
        WRITE (849,'(3(2(f20.14),1x))') (force_is(s,itype),s=1,3)
        natom = natom + atoms%neq(itype)
      END DO ! itype

#ifdef debug
!     test output
!     Evaluate plane wave density in spherical harmonic representation
!     on MT surface along a vector (1,1,1)
      DO itype = 1,atoms%ntype
        WRITE (180,'(2i3)') jsp,itype
        WRITE (180,'(3f20.14)') real(force_is(1,itype)), real(force_is(2,itype)), real(force_is(3,itype))
        vec(1) =-0.730! 1.0
        vec(2) = 0.938! 1.0
        vec(3) =-1.625! 1.0
        CALL ylm4(atoms%lmax,vec,ylm)
        factorrho = 0.0
        factorv   = 0.0
        DO l = 0,atoms%lmax
        DO m = -l,l
          lm = l*(l+1) + m + 1
          factorrho = factorrho + testrho(lm,itype) * ylm(lm)
          factorv   = factorv   +   testV(lm,itype) * ylm(lm)
          WRITE (180,'(3i4,2f20.14)') l,m,lm,testV(lm,itype)
        END DO ! jp
        END DO ! lapw%kp
        WRITE (180,'(a5,2f20.14)') 'den: ',factorrho
        WRITE (180,'(a5,2f20.14)') 'pot: ',factorv
      END DO ! itype
      DEALLOCATE ( ylm,testrho,testV )
#endif

      DEALLOCATE ( bsl,rho,V )
!       DEALLOCATE ( pylm )
      DEALLOCATE ( pylm2 )
!       DEALLOCATE ( qpw2 )
!       DEALLOCATE ( qpwcalc )

      END SUBROUTINE force_sf_is



      SUBROUTINE force_sf_mt(&
                            atoms,sphhar,jspin,&
                            ispin,mpi,&
                            vr,excr,&
                            vxcr,rho,&
                            sym,cell )
!     *****************************************************************
!     This subroutine calculates the contribution evaluated with
!     quantities from the muffin tin
!     n rho V = sum(nu,nup) rho(nu)V(nup)
!             * sum(m,mu,mup) c_1m* Y_1m* c_lnu Y_numu c_lnup Y_nupmup
!     It is called in a spin loop at the end of cdnval.F
!     *****************************************************************

      USE m_constants, ONLY : tpi_const,sfp_const
      USE m_gaunt
      USE m_ylm
      USE m_types
      IMPLICIT NONE
      TYPE(t_mpi),INTENT(IN)   :: mpi
      TYPE(t_sym),INTENT(IN)   :: sym
      TYPE(t_cell),INTENT(IN)  :: cell
      TYPE(t_sphhar),INTENT(IN):: sphhar
      TYPE(t_atoms),INTENT(IN) :: atoms

!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspin
      INTEGER, INTENT (IN) :: ispin 

!     .. Array Arguments ..
      REAL   , INTENT (IN) :: vr(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd) ! 
      REAL   , INTENT (IN) :: rho(:,0:,:,:)!(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd,dimension%jspd)
      REAL   , INTENT (IN) :: excr(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd)
      REAL   , INTENT (IN) :: vxcr(:,0:,:,:)!(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd,dimension%jspd)

!     .. Local Scalars ..
      INTEGER :: natom,itype,nd,lh,l,lhp,lp,mem,m,memp,mp,t,i,lmp
      REAL    :: pot,den
      COMPLEX :: img,factor

!     .. Local Arrays ..
      COMPLEX :: coeff(3,-1:1)
      COMPLEX :: d1((atoms%lmaxd+1)**2,atoms%ntype ),d2((atoms%lmaxd+1)**2,atoms%ntype )
#ifdef debug
      COMPLEX :: testrho((atoms%lmaxd+1)**2 ,atoms%ntype),testv((atoms%lmaxd+1)**2,atoms%ntype )
      COMPLEX :: ylm((atoms%lmaxd+1)**2),factorrho,factorv
      REAL    :: r,Gext(3)
      INTEGER :: lm
      testrho = 0.0
      testv   = 0.0

      WRITE (181,'(a2)') 'MT'
#endif
      WRITE (1704,*) 'vxcr:',vxcr
      WRITE (1704,*) 'excr:',excr

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

!     Calculate forces: For each atom, loop over all lattice harmonics.
      natom = 1
      DO itype = 1,atoms%ntype

        nd = atoms%ntypsy(natom)

        DO lh = 0,sphhar%nlh(nd)
          l = sphhar%llh(lh,nd)

!         The l=0 component of the potential array is saved with a
!         factor r/sfp in front of it. For this calculation, we need
!         the pure potential
          pot = vr(atoms%jri(itype),lh,itype)
          IF (lh.eq.0) THEN
            pot = pot*sfp_const/atoms%rmt(itype)
          END IF
          pot = +1*excr(atoms%jri(itype),lh,itype) -1*vxcr(atoms%jri(itype),lh,itype,ispin) +1*pot
!           pot = 0
!           IF (l.eq.0) pot = sfp

          WRITE (400,'(3(i4,1x),f20.14)') itype,lh,l,vr(atoms%jri(itype),lh,itype)
          WRITE (401,'(3(i4,1x),f20.14)') itype,lh,l,rho(atoms%jri(itype),lh,itype,ispin)
          DO mem = 1,sphhar%nmem(lh,nd)
            m = sphhar%mlh(mem,lh,nd)
            lmp = l*(l+1) + m + 1
            d1(lmp,itype) = d1(lmp,itype) + sphhar%clnu(mem,lh,nd)&
                 *rho(atoms%jri(itype),lh,itype,ispin)/atoms%rmt(itype)**2
            d2(lmp,itype) = d2(lmp,itype) + sphhar%clnu(mem,lh,nd) * pot
          END DO ! mem

          DO lhp = 0,sphhar%nlh(nd)
            lp = sphhar%llh(lhp,nd)
            IF (abs(l-lp).ne.1) CYCLE

            den = rho(atoms%jri(itype),lhp,itype,ispin) ! jspin ?
!             den = 0
!             IF (lp.eq.0) den = sfp*rmt(itype)**2

            DO mem = 1,sphhar%nmem(lh,nd)
              m = sphhar%mlh(mem,lh,nd)
            DO memp = 1,sphhar%nmem(lhp,nd)
              mp = sphhar%mlh(memp,lhp,nd)
              IF (abs(m+mp).gt.1) CYCLE

!             Due to the normal vector n, the lattice harmonics meet
!             with a Y_1m resulting in a gaunt coefficient.
              factor = pot * den * sphhar%clnu(mem,lh,nd) * sphhar%clnu(memp,lhp,nd)&
                    * gaunt1(1,l,lp,m+mp,m,mp,atoms%lmaxd)
              force_mt(:,itype) = force_mt(:,itype) + factor * conjg(coeff(:,m+mp))

            END DO ! memp
            END DO ! mem

          END DO ! lhp


#ifdef debug
!         testrho/v debug code
!         construct spherical harmonic expansion of lattice harmonic
!         density and potential on the muffin-tin boundary
          DO mem = 1,sphhar%nmem(lh,nd)
            m = sphhar%mlh(mem,lh,nd)
          DO lp = 0,atoms%lmaxd
            IF (l.ne.lp) CYCLE
          DO mp = -lp,lp
            IF (m.ne.mp) CYCLE
            lm = lp*(lp+1) + mp + 1
            testrho(lm,itype) = testrho(lm,itype)&
                 + rho(atoms%jri(itype),lh,itype,ispin) &
                 * sphhar%clnu(mem,lh,nd) / atoms%rmt(itype)**2
            testv(lm,itype)   = testv(lm,itype) + pot * sphhar%clnu(mem,lh,nd)
          END DO ! mp
          END DO ! lp
          END DO ! mem
#endif

        END DO ! lh

        DO l = 0,atoms%lmaxd
        DO m = -l,l
          lmp = l*(l+1) + m + 1
          WRITE (1706,*) itype,l,m,lmp-1,d1(lmp,itype),d2(lmp,itype)
        END DO ! m
        END DO ! l

        IF (mpi%irank.eq.0) THEN
          WRITE (850,'(3(2(f20.14),1x))') (force_mt(t,itype),t=1,3)
        END IF
        natom = natom + atoms%neq(itype)
      END DO ! itype

!     debug
      WRITE (*,*) 'sanity check'
      WRITE (*,*) 'look is',force_is
      WRITE (*,*) 'look mt',force_mt

#ifdef debug
!     test output
!     Evaluate lattice harmonic density in spherical harmonic
!     representation in along a reciprocal vector (1,1,1)
      DO itype = 1,atoms%ntype
        WRITE (181,'(2i3)') jspin,itype
        WRITE (181,'(3f20.14)') real(force_mt(1,itype)), real(force_mt(2,itype)), real(force_mt(3,itype))
        Gext(1) =-0.730! 1.0
        Gext(2) = 0.938! 1.0
        Gext(3) =-1.625! 1.0
        CALL ylm4(atoms%lmaxd,Gext,ylm)
        factorrho = 0.0
        factorv   = 0.0
        DO lp = 0,atoms%lmaxd
        DO mp = -lp,lp
          lm = lp*(lp+1) + mp + 1
          factorrho = factorrho + testrho(lm,itype) * ylm(lm)
          factorv   = factorv   + testv(lm,itype)   * ylm(lm)
          WRITE (181,'(3i4,2f20.14)') lp,mp,lm,testv(lm,itype)
        END DO ! mp
        END DO ! lp
        WRITE (181,'(a5,2f20.14)') 'den: ',factorrho
        WRITE (181,'(a5,2f20.14)') 'pot: ',factorv
      END DO ! itype
#endif

      END SUBROUTINE force_sf_mt

      SUBROUTINE init_sf(sym,cell,atoms)
!     Initialize results arrays if neither force_sf_is nor force_sf_mt
!     were executed up until now.
!     Called at the beginning of cdnval.F to once fill the wigner array
!     Also called in force_sf_is/mt to guarantee that
!     all arrays are allocated.

    USE m_types
      IMPLICIT NONE
      TYPE(t_sym),INTENT(IN)   :: sym
      TYPE(t_cell),INTENT(IN)  :: cell
      TYPE(t_atoms),INTENT(IN) :: atoms

!     debug
      IF (ALLOCATED(force_is)) THEN
        WRITE (*,*) 'init is:',force_is
        WRITE (*,*) 'init mt:',force_mt
      ELSE
        WRITE (*,*) 'init: not initialized'
      END IF

      IF (isdone.OR.mtdone) RETURN
      IF (.not.ALLOCATED(force_is)) THEN
        ALLOCATE ( force_is(3,atoms%ntypd),force_mt(3,atoms%ntypd) )
      END IF
      force_is = 0.0
      force_mt = 0.0

        WRITE (*,*) 'init: end'
      END SUBROUTINE init_sf

      SUBROUTINE exit_sf(isp,atoms,force)
!     Write out force contribution from surface and deallocate arrays if
!     all force_sf_is and force_sf_mt were executed
!     Called at the end of cdnval.F and totale.f for writing purposes
!     and to deallocate the arrays.

    USE m_types
      IMPLICIT NONE
      INTEGER,INTENT(IN)         :: isp
      TYPE(t_atoms),INTENT(IN)   :: atoms
      REAL,INTENT(INOUT)         :: force(:,:,:)

      INTEGER :: itype,dir
      COMPLEX :: force_sf(3,atoms%ntypd)

!     debug
      IF (ALLOCATED(force_is)) THEN
        WRITE (*,*) 'exit is:',force_is
        WRITE (*,*) 'exit mt:',force_mt
      ELSE
        WRITE (*,*) 'exit: not initialized'
      END IF

      IF (isdone.AND.mtdone) THEN
        force_sf(:,:) = force_is(:,:) - force_mt(:,:)
        force(:,:,isp) = force(:,:,isp) + real(force_sf(:,:))
        WRITE (6,*)
        WRITE (16,*)
        DO itype = 1,atoms%ntype
          WRITE (6,FMT=8010) itype
          WRITE (16,FMT=8010) itype
          WRITE (6,FMT=8020) (force_sf(dir,itype),dir=1,3)
          WRITE (16,FMT=8020) (force_sf(dir,itype),dir=1,3)
        END DO ! itype
        isdone = .false.
        mtdone = .false.
        DEALLOCATE ( force_is,force_mt )
      END IF

8010  FORMAT (' FORCES: SURFACE CORRECTION FOR ATOM TYPE',i4)
8020  FORMAT (' FX_SF=',2f10.6,' FY_SF=',2f10.6,' FZ_SF=',2f10.6)

      END SUBROUTINE exit_sf

      END MODULE m_force_sf
