MODULE m_atom2
   use m_juDFT
!     *************************************************************
!     fully relativistic atomic program based on the subroutines
!     differ, outint and inwint by d.d.koelling
!     erich wimmer     august 1981
!     modified for use in start-up generator.  m.w. april 1982
!     modified by adding stabilizing well. m. weinert may 1990
!     *************************************************************
CONTAINS
   SUBROUTINE atom2(&
  &                 dimension, atoms, xcpot, input, ntyp, jrc, rnot1,&
  &                 qdel,&
  &                 rhoss, nst, lnum, eig, vbar)

      USE m_intgr, ONLY: intgr1, intgr0
      USE m_constants
      USE m_potl0
      USE m_stpot1
    !  USE m_setcor
      USE m_differ
      USE m_types
      IMPLICIT NONE
!     ..
!     .. Scalar Arguments ..
      TYPE(t_dimension), INTENT(IN)  :: dimension
      TYPE(t_atoms), INTENT(IN)      :: atoms
      CLASS(t_xcpot), INTENT(IN)     :: xcpot
      TYPE(t_input), INTENT(IN)      :: input
      INTEGER, INTENT(IN)  :: jrc, ntyp
      REAL, INTENT(IN)  :: rnot1, qdel
      REAL, INTENT(OUT) :: rhoss(:, :) !(mshd,input%jspins),
      REAL, INTENT(OUT) :: eig(dimension%nstd, input%jspins), vbar(input%jspins)
      INTEGER, INTENT(OUT) :: nst, lnum(dimension%nstd)
!     ..
!     .. Local Scalars ..
      REAL c, d, delrv, dist, distol, e, fisr, fj, fl, fn, h,&
     &     p, p1, pmax, pmin, r, r3, rn, rnot, z, zero, bmu_l, rho
      INTEGER i, inr0, it, itmax, k, l, n, ispin, kk, ierr, msh_l
      LOGICAL conv, lastit, l_start
!     ..
!     .. Local Arrays ..
      REAL a(jrc), b(jrc), dens(jrc), occ(dimension%nstd, input%jspins)
      REAL rad(jrc), rev(dimension%nstd, input%jspins), ahelp(jrc), ain(jrc),&
     &     rh(jrc), vr(jrc), f(0:3),&
     &     vr1(jrc, input%jspins), vr2(jrc, input%jspins), vx(dimension%msh, input%jspins), vxc(dimension%msh, input%jspins)
      INTEGER kappa(dimension%nstd), nprnc(dimension%nstd)
!     ..
!     ..
!     .. Data statements ..
!---->     distol set from 1.0e-6 to 1.0e-3
      DATA zero, distol/0.0e0, 1.0e-3/
!     ..
      c = c_light(1.0)
      vxc(:, :) = 0.0
      vx(:, :) = 0.0
!
      WRITE (6, FMT=8000)
8000  FORMAT(' subroutine atom2 entered')
      z = atoms%zatom(ntyp)
      n = jrc
      rnot = rnot1
      itmax = 100
      pmin = 0.01e0
      pmax = 0.2e0
      h = atoms%dx(ntyp)
      d = exp(h)
      r = rnot
      DO i = 1, n
         rad(i) = r
         r = r*d
      enddo
      rn = rad(n)
      bmu_l = atoms%bmu(ntyp)
      !IF (bmu_l > 0.001 .AND. atoms%numStatesProvided(ntyp) .NE. 0) CALL &
      !   judft_warn("You specified both: inital moment and occupation numbers.", &
      !              hint="The inital moment will be ignored, set magMom=0.0", calledby="atom2.f90")
      !CALL setcor(ntyp, input%jspins, atoms, input, bmu_l, nst, kappa, nprnc, occ)
      CALL atoms%econf(ntyp)%get_core(nst,nprnc,kappa,occ)

!
!--->   for electric field case (sigma.ne.0), add the extra charge
!--->   to the uppermost level; ignore the possible problem that
!--->   the occupations may not be between 0 and 2
      IF (input%jspins == 1) THEN
         occ(nst, 1) = occ(nst, 1) + qdel
      ELSE
         occ(nst, 1) = occ(nst, 1) + qdel/2.
         occ(nst, input%jspins) = occ(nst, input%jspins) + qdel/2.
      ENDIF
!
      CALL stpot1(&
     &           jrc, n, z, rad,&
     &           vr1)
      DO i = 1, n
         vr1(i, input%jspins) = vr1(i, 1)
      ENDDO
!
!     start iterating
!
      lastit = .false.
      conv = .true.
      delrv = 0.100
      inr0 = log(5.0/rnot)/h + 1.5
      DO it = 1, itmax
         DO ispin = 1, input%jspins
!
!---->     load potential
            DO i = 1, n
               vr(i) = vr1(i, ispin)
            ENDDO
!----> adding stabilizing well: m. weinert
            DO i = inr0, n
               vr(i) = vr(i) + rad(i)*delrv*(rad(i) - rad(inr0))
            ENDDO
!---->     note that vr contains r*v(r) in hartree units
            DO i = 1, n
               rhoss(i, ispin) = zero
            ENDDO
            if (lastit) THEN
               inquire (file="startcharges", exist=l_start)
               if (l_start) then
                  OPEN (61, file="startcharges")
                  DO WHILE (.true.)
                     read (61, *, end=888, err=888) i, rho
                     if (i == z) then
                        occ(nst, 1) = occ(nst, 1) + rho
                        goto 888
                     endif
                  enddo
888               continue
                  close (61)
               endif
            endif
            DO 90 k = 1, nst
               fn = nprnc(k)
               fj = iabs(kappa(k)) - 0.5e0
               fl = fj + 0.5e0*isign(1, kappa(k))
               e = -2*(z/(fn + fl))**2
               ierr = -1
               msh_l = jrc
               DO WHILE (ierr .NE. 0)
                  CALL differ(&
       &                      fn, fl, fj, c, z, h, rnot, rn, d, msh_l, vr,&
       &                      e,&
       &                      a, b, ierr)!keep
                  msh_l = msh_l - 1
                  IF (jrc - msh_l > 100) CALL juDFT_error(&
       &               "atom2", calledby="atom2")
               ENDDO
               DO i = msh_l + 1, jrc
                  a(i) = a(msh_l)
                  b(i) = b(msh_l)
               ENDDO

               DO i = 1, n
                  rh(i) = occ(k, ispin)*(a(i)**2 + b(i)**2)
               ENDDO
!+ldau
               IF (lastit) THEN                         ! calculate slater interals
                  l = int(fl)
!                 write(*,*) nprnc(k),l
                  DO kk = 0, 2*l, 2                      ! F0 for s, F0 + F2 for p etc.
                     r = rnot
                     DO i = 1, n
                        ain(i) = a(i)**2*r**(-kk - 1)      ! prepare inner integrand
                        r = r*d
                     ENDDO
                     CALL intgr1(ain, rnot, h, n, &          ! integrate&
                                &ahelp)
                     r = rnot
                     DO i = 1, n - 1
                        ain(i) = a(i)**2*r**kk*(ahelp(n) - ahelp(i))
                        r = r*d
                     ENDDO
                     CALL intgr0(ain, rnot, h, n - 1,      &     ! integrate 2nd r&
                                &f(kk/2))

                  ENDDO
!                 write(*,*) (hartree_to_ev_const*2*f(kk),kk=0,l)
               ENDIF
!-ldau
               eig(k, ispin) = e
!---->       calculate <r>
               DO i = 1, n
                  a(i) = (a(i)**2 + b(i)**2)*rad(i)
               ENDDO
               CALL intgr1(a, rnot, h, n, b)
               rev(k, ispin) = b(n)
               DO i = 1, n
                  rhoss(i, ispin) = rhoss(i, ispin) + rh(i)
               ENDDO
90          ENDDO
         ENDDO
!
!     solve poisson's equation
!
         DO i = 1, n
            dens(i) = rhoss(i, 1)
         ENDDO
         IF (input%jspins == 2) THEN
            DO i = 1, n
               dens(i) = dens(i) + rhoss(i, input%jspins)
            ENDDO
         ENDIF
         CALL intgr1(dens, rnot, h, n, a)
         DO i = 1, n
            rh(i) = dens(i)/rad(i)
         ENDDO
         CALL intgr1(rh, rnot, h, n, b)
         fisr = b(n)
         DO i = 1, n
            vr(i) = (a(i) + rad(i)*(fisr - b(i)) - z)
         ENDDO
!+ta
         DO ispin = 1, input%jspins
            DO i = 1, n
               rhoss(i, ispin) = rhoss(i, ispin)/(fpi_const*rad(i)**2)
            ENDDO
         ENDDO
         IF (xcpot%needs_grad()) THEN
            CALL potl0(xcpot, input%jspins, atoms%dx(ntyp), rad, rhoss, vxc)
         ELSE
            CALL xcpot%get_vxc(input%jspins, rhoss, vxc, vx)
         ENDIF
         DO ispin = 1, input%jspins
            DO i = 1, n
               vr2(i, ispin) = vr(i) + vxc(i, ispin)*rad(i)
            ENDDO
         ENDDO
!-ta
!        determine distance of potentials
!
         r3 = rn**3
         dist = 0.0
         DO ispin = 1, input%jspins
            DO i = 1, n
               a(i) = (vr2(i, ispin) - vr1(i, ispin))**2
            ENDDO
            CALL intgr1(a, rnot, h, n, b)
            dist = dist + sqrt((3.0e0/r3)*b(n))
         ENDDO
         IF (lastit) GO TO 190
         IF (dist < distol) lastit = .true.
!     mix new input potential
         p1 = 1.0e0/dist
         p = min(pmax, p1)
         p = max(p, pmin)
         WRITE (6, FMT=8060) it, dist, p
         p1 = 1.0e0 - p
         DO ispin = 1, input%jspins
            DO i = 1, n
               vr1(i, ispin) = p1*vr1(i, ispin) + p*vr2(i, ispin)
            ENDDO
         ENDDO
      ENDDO
!
! output
!
      WRITE (6, FMT=8030) dist
      conv = .false.
!     list eigenvalues
190   IF (conv) WRITE (6, FMT=8040) it, dist
      DO ispin = 1, input%jspins
         WRITE (6, '(a8,i2)') 'spin No.', ispin
         DO k = 1, nst
            fj = iabs(kappa(k)) - 0.5e0
            l = fj + 0.5e0*isign(1, kappa(k)) + 0.01e0
            lnum(k) = l
            WRITE (6, FMT=8050) nprnc(k), kappa(k), l, fj,&
      &                        occ(k, ispin), eig(k, ispin), rev(k, ispin)
         ENDDO
!
!--->   guess enpara if it doesn't exist, using floating energy parameters
!
         i = atoms%jri(ntyp) - (log(4.0)/atoms%dx(ntyp) + 1.51)
         vbar(ispin) = vr1(i, ispin)/(rnot*exp(atoms%dx(ntyp)*(i - 1)))
         WRITE (6, '(/,'' reference energy = '',2f12.6,/)') vbar(ispin)
      ENDDO

8030  FORMAT(/, /, /, ' $$$ error: not converged, dist=', f10.6,/)
8040  FORMAT(/, /, 3x, 'converged in', i4, ' iterations to a distance of',&
      &       e12.5, ' har', /, /, 3x, 'n  kappa  l    j  ', 5x,&
      &       'occ.   eigenvalue (har)  <r>  ',/)
8050  FORMAT(3x, i1, i5, i5, f6.1, 2(3x, f7.2, 1x, 2f12.6))
8060  FORMAT('it,dist,p=', i4, 2f12.5)

   END SUBROUTINE atom2
END MODULE m_atom2
