module m_structureconstant
   USE m_types
   USE m_juDFT
   USE m_constants
   use m_ylm
   use m_sort
#ifdef CPP_MPI
   use mpi
#endif
contains
   !     -----------------------------------------------------------------------------------------------

   !     Calculates the structure constant
   !                                                        1               *      ^
   !     structconst(lm,ic1,ic2,k) = SUM exp(ikT) -----------------------  Y  ( T + R(ic) )
   !                                  T           | T + R(ic1) - R(ic2) |   lm
   !
   !     with T = lattice vectors
   !
   !     An Ewald summation method devised by O.K. Andersen is used for l<5
   !     (see e.g. H.L. Skriver, "The LMTO method", Springer 1984).
   !     (The real-space function G can be calculated with gfunction.f)
   !

   SUBROUTINE structureconstant(structconst, cell, hybinp, atoms, kpts, fmpi)
      IMPLICIT NONE
      TYPE(t_mpi), INTENT(IN)    :: fmpi
      TYPE(t_hybinp), INTENT(IN) :: hybinp
      TYPE(t_cell), INTENT(IN)   :: cell
      TYPE(t_atoms), INTENT(IN)  :: atoms
      TYPE(t_kpts), INTENT(IN)   :: kpts
      ! - scalars -

      ! - arrays -
      COMPLEX, INTENT(INOUT)   ::  structconst(:, :, :, :)

      ! - local scalars -
      INTEGER                   ::  i, ic1, ic2, lm, ikpt, l, ishell, nshell
      INTEGER                   ::  m
      INTEGER                   ::  nptsh, maxl

      REAL                      ::  rad, rrad, rdum
      REAL                      ::  a, a1, aa
      REAL                      ::  pref, rexp
      REAL                      ::  scale

      COMPLEX                   ::  cexp

      LOGICAL, SAVE             ::  first = .TRUE.
      logical                   ::  run_loop
      ! - local arrays -
      INTEGER                   ::  conv(0:2*hybinp%lexp), ierr, buf_sz, root
      INTEGER, ALLOCATABLE     ::  ptsh(:, :)

      REAL                      ::  k(3), ki(3), ka(3)
      REAL                      ::  convpar(0:2*hybinp%lexp), g(0:2*hybinp%lexp)
      REAL, ALLOCATABLE     ::  radsh(:)

      COMPLEX                   ::  y((2*hybinp%lexp + 1)**2)
      REAL, PARAMETER           :: CONVPARAM = 1e-18
      ! Do some additional shells ( real-space and Fourier-space sum )
      INTEGER, PARAMETER        :: ADDSHELL2 = 0
      
      call timestart("calc struc_const.")

      IF (fmpi%irank /= 0) first = .FALSE.

      rdum = cell%vol**(1.0/3) ! define "average lattice parameter"

      ! ewaldlambda = ewaldscale
      scale = hybinp%ewaldlambda/rdum

      !       lambda = ewaldlambda / rdum

      pref = fpi_const/(scale**3*cell%vol)

      DO l = 0, 2*hybinp%lexp
         convpar(l) = CONVPARAM/scale**(l + 1)
      END DO

      IF (first) THEN
         WRITE (oUnit, '(//A)') '### subroutine: structureconstant ###'
         WRITE (oUnit, '(/A)') 'Real-space sum:'
      END IF

      !
      !     Determine cutoff radii for real-space and Fourier-space summation
      ! (1) real space
      call timestart("determine cutoff radii")

      a = 0
      run_loop = .True.
      do while(run_loop)
         a = a + 1
         rexp = EXP(-a)
         g(0) = rexp/a*(1 + a*11/16*(1 + a*3/11*(1 + a/9)))
         g(1) = rexp/a**2*(1 + a*(1 + a/2*(1 + a*7/24*(1 + a/7))))
         g(2) = rexp/a**3*(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a*3/16 &
                                                            *(1 + a/9))))))
         g(3) = rexp/a**4*(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a/5*(1 + a/6 &
                                                                     *(1 + a/8)))))))
         g(4) = rexp/a**5*(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a/5*(1 + a/6 &
                                                                     *(1 + a/7*(1 + a/8*(1 + a/10)))))))))
         g(5) = rexp/a**6*(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a/5*(1 + a/6 &
                                                                     *(1 + a/7*(1 + a/8*(1 + a/9*(1 + a/10))))))))))
         g(6) = rexp/a**7*(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a/5*(1 + a/6 &
                                                                     *(1 + a/7*(1 + a/8*(1 + a/9*(1 + a/10*(1 + a/11*(1 + a/12))))))))))))
         g(7) = rexp/a**8*(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a/5*(1 + a/6 &
                                                                     *(1 + a/7*(1 + a/8*(1 + a/9*(1 + a/10*(1 + a/11*(1 + a/12*(1 + a/13)))))))))))))
         DO l = 8, 2*hybinp%lexp
            g(l) = a**(-l - 1)
         END DO
         run_loop = ANY(g > convpar/10)
      enddo
      rad = a/scale
      call timestop("determine cutoff radii")

      ! (2) Fourier space
      call timestart("fourier space")
      a = 0 
      run_loop = .True.
      do while(run_loop)
         a = a + 1
         aa = (1 + a**2)**(-1)
         g(0) = pref*aa**4/a**2
         g(1) = pref*aa**4/a
         g(2) = pref*aa**5/3
         g(3) = pref*aa**5*a/15
         g(4) = pref*aa**6*a**2/105
         g(5) = pref*aa**6*a**3/945
         g(6) = pref*aa**7*a**4/10395
         g(7) = pref*aa**7*a**5/135135
         run_loop = ANY(g > convpar)
      enddo
      ! IF (ANY(g > convpar)) THEN
      !    a = a + 1
      !    GOTO 2
      ! END IF
      rrad = a*scale
      call timestop("fourier space")

      IF (first) THEN
         WRITE (oUnit, '(/A,2F10.5)') 'Cutoff radii: ', rad, rrad
         WRITE (oUnit, '(/A)') 'Real-space sum'
      END IF

      call realspace_sum(atoms, cell, hybinp, fmpi, kpts, first, scale, convpar, g, a, a1, rad, structconst)
      
      IF (first) WRITE (oUnit, '(/A)') 'Fourier-space sum'

      !
      !     Determine reciprocal shells
      !
      call timestart("determince reciproc. shell")
      CALL getshells(ptsh, nptsh, radsh, nshell, rrad, cell%bmat, first)
      call timestop("determince reciproc. shell")
      ! minimum nonzero reciprocal-shell radius (needed in routines concerning the non-local hartree-fock exchange)
      !hybinp%radshmin = radsh(2)
      !
      !     Fourier-space sum
      !
      call timestart("fourierspace sum")
      DO ikpt = 1, kpts%nkpt
         k = kpts%bk(:, ikpt)
         maxl = MIN(7, hybinp%lexp*2)
         ishell = 1
         conv = HUGE(i)
         DO i = 1, nptsh
            IF (i > 1) THEN
               IF (ABS(radsh(i) - radsh(i - 1)) > 1e-10) ishell = ishell + 1
            ENDIF
            ki = ptsh(:, i) + k - NINT(k) ! -nint(...) transforms to Wigner-Seitz cell ( i.e. -0.5 <= x,y,z < 0.5 )
            ka = MATMUL(ki, cell%bmat)
            a = norm2(ka)/scale
            aa = (1 + a**2)**(-1)
            IF (ABS(a - a1) > 1e-10) THEN
               a1 = a
               IF (abs(a) < 1e-12) THEN
                  g(0) = pref*(-4)
                  g(1) = 0
               ELSE
                  IF (ishell <= conv(0)) g(0) = pref*aa**4/a**2
                  IF (ishell <= conv(1)) g(1) = pref*aa**4/a
               END IF
               IF (ishell <= conv(2)) g(2) = pref*aa**5/3
               IF (ishell <= conv(3)) g(3) = pref*aa**5*a/15
               IF (ishell <= conv(4)) g(4) = pref*aa**6*a**2/105
               IF (ishell <= conv(5)) g(5) = pref*aa**6*a**3/945
               IF (ishell <= conv(6)) g(6) = pref*aa**7*a**4/10395
               IF (ishell <= conv(7)) g(7) = pref*aa**7*a**5/135135
               IF (ishell > 1) THEN
                  DO l = 0, 7
                     IF (conv(l) == HUGE(i) .AND. g(l) < convpar(l)) conv(l) = ishell + ADDSHELL2
                  END DO
               END IF
            END IF

            IF (ishell > conv(maxl) .AND. maxl /= 0) maxl = maxl - 1
            call ylm4(maxl, ka, y)
            IF (norm2(ka(:)) .LT. 1.0e-16) y(2:(maxl + 1)**2) = CMPLX(0.0, 0.0)
            lm = 0
            !$OMP PARALLEL default(none) &
            !$OMP private(l, M, lm, ic1, ic2, cexp) &
            !$OMP shared(ishell, conv, g, y, maxl, structconst, atoms, ikpt, ki, fmpi)

            !$OMP DO schedule(dynamic)
            DO l = 0, maxl
               lm = l**2
               IF (ishell <= conv(l)) THEN
                  DO M = -l, l
                     lm = lm + 1
                     y(lm) = CONJG(y(lm))*g(l)*ImagUnit**l
                  END DO
               ELSE
                  y(lm + 1:lm + 2*l + 1) = 0
               END IF
            END DO
            !$OMP END DO

            !$OMP DO schedule(dynamic) collapse(2)
            DO ic2 = 1+fmpi%irank, atoms%nat, fmpi%isize
               DO ic1 = 1, atoms%nat
                  IF (ic2 /= 1 .AND. ic1 == ic2) CYCLE
                  cexp = EXP(ImagUnit*tpi_const*dot_PRODUCT(ki, atoms%taual(:, ic1) - atoms%taual(:, ic2)))
                  DO lm = 1, (maxl + 1)**2
                     structconst(lm, ic1, ic2, ikpt) = structconst(lm, ic1, ic2, ikpt) + cexp*y(lm)
                  END DO
               END DO
            END DO
            !$OMP END DO
            !$OMP END PARALLEL

         END DO
#ifdef CPP_MPI
         call timestart("bcast fouriersum")
         buf_sz = size(structconst,1) * size(structconst,2)
         DO ic2 = 1, atoms%nat
            root = mod(ic2-1, fmpi%isize)
            call MPI_Bcast(structconst(1,1,ic2,ikpt), buf_sz, MPI_DOUBLE_COMPLEX, root, fmpi%mpi_comm, ierr)
         enddo
         call timestop("bcast fouriersum")
#endif
      END DO
      call timestop("fourierspace sum")
      !
      !     Add contribution for l=0 to diagonal elements and rescale structure constants
      !
      structconst(1, 1, 1, :) = structconst(1, 1, 1, :) - 5.0/16/SQRT(fpi_const)
      DO i = 2, atoms%nat
         structconst(:, i, i, :) = structconst(:, 1, 1, :)
      END DO
      DO l = 0, 2*hybinp%lexp
         structconst(l**2 + 1:(l + 1)**2, :, :, :) = structconst(l**2 + 1:(l + 1)**2, :, :, :)*scale**(l + 1)
      END DO

      rad = (cell%vol*3/4/pi_const)**(1.0/3) ! Wigner-Seitz radius (rad is recycled)

      !     Calculate accuracy of Gamma-decomposition
      IF (ALL(abs(kpts%bk) > 1e-12)) THEN
         a = 1e30 ! ikpt = index of shortest non-zero k-point
         DO i = 2, kpts%nkpt
            rdum = SUM(MATMUL(kpts%bk(:, i), cell%bmat)**2)
            IF (rdum < a) THEN
               ikpt = i
               a = rdum
            END IF
         END DO
         rdum = norm2(MATMUL(kpts%bk(:, ikpt), cell%bmat))
         a = 0
         DO ic2 = 1, atoms%nat
            DO ic1 = 1, MAX(1, ic2 - 1)
               a = a + ABS(structconst(1, ic1, ic2, ikpt) - &
                           (structconst(1, ic1, ic2, 1) + SQRT(fpi_const)/cell%vol/rdum**2* &
                            EXP(-CMPLX(0.0, 1.0)*tpi_const*dot_PRODUCT( &
                                kpts%bk(:, ikpt), atoms%taual(:, ic2) - atoms%taual(:, ic1)))))**2
            END DO
         END DO
         a = SQRT(a/atoms%nat**2)
         aa = SQRT(SUM(ABS(structconst(1, :, :, ikpt))**2)/atoms%nat**2)
         IF (first) WRITE (oUnit, '(/A,F8.5,A,F8.5,A)') 'Accuracy of Gamma-decomposition (structureconstant):', a, ' (abs)', a/aa, ' (rel)'
      ENDIF

      deallocate (ptsh, radsh)

      first = .FALSE.

      call timestop("calc struc_const.")
   END SUBROUTINE structureconstant

      !     -----------------

   !     Determines all shells of the crystal defined by lat and vol with radii smaller than rad.
   !     The lattice points (number = nptsh) are stored in ptsh, their corresponding lengths (shell radii) in radsh.

   SUBROUTINE getshells(ptsh, nptsh, radsh, nshell, rad, lat, lwrite)

      USE m_juDFT
      USE m_constants

      IMPLICIT NONE

      LOGICAL, INTENT(IN)    :: lwrite
      INTEGER, INTENT(INOUT)   :: nptsh, nshell
      INTEGER, ALLOCATABLE   :: ptsh(:, :)
      REAL, ALLOCATABLE   :: radsh(:)
      REAL, INTENT(IN)    :: rad, lat(:, :)
      REAL                    :: r(3), rdum
      INTEGER, ALLOCATABLE   :: pnt(:)
      INTEGER                 :: n, i, ix, iy, iz, ok
      LOGICAL                 :: found
      INTEGER, ALLOCATABLE   :: ihelp(:, :)
      REAL, ALLOCATABLE   :: rhelp(:)

      allocate (ptsh(3, 100000), radsh(100000), stat=ok)
      IF (ok /= 0) call judft_error('getshells: failure allocation ptsh/radsh')

      ptsh = 0
      radsh = 0

      n = 0
      i = 0
      DO
         found = .FALSE.
         DO ix = -n, n
            DO iy = -(n - ABS(ix)), n - ABS(ix)
               iz = n - ABS(ix) - ABS(iy)
1              r = ix*lat(:, 1) + iy*lat(:, 2) + iz*lat(:, 3)
               rdum = SUM(r**2)
               IF (rdum < rad**2) THEN
                  found = .TRUE.
                  i = i + 1
                  IF (i > SIZE(radsh)) THEN
                     allocate (rhelp(SIZE(radsh)), ihelp(3, SIZE(ptsh, 2)), stat=ok)
                     IF (ok /= 0) call judft_error('getshells: failure allocation rhelp/ihelp')
                     rhelp = radsh
                     ihelp = ptsh
                     deallocate (radsh, ptsh)
                     allocate (radsh(SIZE(rhelp) + 100000), ptsh(3, SIZE(ihelp, 2) + 100000), stat=ok)
                     IF (ok /= 0) call judft_error('getshells: failure re-allocation ptsh/radsh')
                     radsh(1:SIZE(rhelp)) = rhelp
                     ptsh(:, 1:SIZE(ihelp, 2)) = ihelp
                     deallocate (rhelp, ihelp)
                  END IF
                  ptsh(:, i) = [ix, iy, iz]
                  radsh(i) = SQRT(rdum)
               END IF
               IF (iz > 0) THEN
                  iz = -iz
                  GOTO 1
               END IF
            END DO
         END DO
         IF (.NOT. found) EXIT
         n = n + 1
      END DO
      nptsh = i

      allocate (pnt(nptsh))

      !reallocate radsh ptsh
      allocate (rhelp(nptsh), ihelp(3, nptsh))
      rhelp = radsh(1:nptsh)
      ihelp = ptsh(:, 1:nptsh)
      deallocate (radsh, ptsh)
      allocate (radsh(nptsh), ptsh(3, nptsh))
      radsh = rhelp
      ptsh = ihelp
      deallocate (rhelp, ihelp)

      call sort(pnt, radsh, [(1.0*i,i=1,nptsh)])

      radsh = radsh(pnt)
      ptsh = ptsh(:, pnt)
      nshell = 1
      DO i = 2, nptsh
         IF (radsh(i) - radsh(i - 1) > 1e-10) nshell = nshell + 1
      END DO

      IF (lwrite) &
         WRITE (oUnit, '(A,F10.5,A,I7,A,I5,A)') &
         '  Sphere of radius', rad, ' contains', &
         nptsh, ' lattice points and', nshell, ' shells.'

   END SUBROUTINE getshells

   subroutine realspace_sum(atoms, cell, hybinp, fmpi, kpts, first, scale, convpar, g, a, a1, rad, structconst)
      use ieee_arithmetic
      implicit none 
      type(t_atoms), intent(in) :: atoms 
      type(t_cell), intent(in)  :: cell 
      type(t_hybinp), intent(in):: hybinp
      TYPE(t_mpi), INTENT(IN)    :: fmpi
      type(t_kpts), intent(in)  :: kpts
      logical, intent(in)       :: first
      real, intent(in)          :: rad, scale, convpar(0:2*hybinp%lexp)
      real, intent(inout)       :: g(0:2*hybinp%lexp), a, a1
      complex, intent(inout)    :: structconst(:,:,:,:)
      
      integer :: ic2, ic1, i, ishell, l, m, maxl, lm, ikpt, nptsh, nshell, ierr, s
      integer ::  conv(0:2*hybinp%lexp)
      integer, allocatable ::  pnt(:), ptsh(:,:)
      INTEGER, PARAMETER        :: ADDSHELL1 = 40

      real :: rdum, rexp, ra(3),  rc(3), tmp_vec(3)
      real, allocatable    ::  radsh(:)

      complex  ::  shlp((2*hybinp%lexp + 1)**2, kpts%nkpt)
      COMPLEX  ::  cdum, cexp, y((2*hybinp%lexp + 1)**2)

      rdum = cell%vol**(1.0/3) 

      !
      !     Determine atomic shells
      call timestart("determine atomic shell")
      CALL getshells(ptsh, nptsh, radsh, nshell, rad, cell%amat, first)
      call timestop("determine atomic shell")

      allocate (pnt(nptsh))
      structconst = 0
      a1 = 0

      !
      !     Real-space sum
      !
      call timestart("realspace sum")
      DO ic2 = 1+fmpi%irank, atoms%nat, fmpi%isize
!         !$OMP PARALLEL DO default(none) &
!         !$OMP shared(ic2, atoms, cell, nptsh, structconst, hybinp, kpts, scale, convpar) &
!         !$OMP private(ic1, tmp_vec, i, ra, rc, a, pnt, maxl, l, conv, shlp, ishell, rexp, g, y) &
!         !$OMP private(rdum, cexp, lm, cdum)&
!         !$OMP firstprivate(ptsh, radsh) &
!         !$OMP reduction(max:a1)
         DO ic1 = 1, atoms%nat
            IF (ic2 /= 1 .AND. ic1 == ic2) CYCLE
            !MATMUL(cell%amat, (atoms%taual(:, ic2) - atoms%taual(:, ic1)))
            tmp_vec = atoms%taual(:, ic2) - atoms%taual(:, ic1)
            call dgemv("N", 3, 3, 1.0, cell%amat, 3, tmp_vec, 1, 0.0, rc, 1)
            DO i = 1, nptsh
               !ra = MATMUL(cell%amat, ptsh(:, i)) + rc
               tmp_vec = real(ptsh(:, i))
               call dgemv("N", 3, 3, 1.0, cell%amat, 3, tmp_vec, 1, 0.0, ra, 1)
               ra = ra + rc
               radsh(i) = norm2(ra)
            END DO
            call sort(pnt, radsh, [(1.0*s, s=1,nptsh)])
            ptsh = ptsh(:, pnt)
            radsh = radsh(pnt)
            maxl = 2*hybinp%lexp
            a1 = 1e30  ! stupid initial value
            ishell = 1
            conv = HUGE(i)
            shlp = 0
            DO i = 1, nptsh
               IF (ALL(conv /= HUGE(i))) EXIT
               IF (i /= 1) THEN
                  IF (ABS(radsh(i) - radsh(i - 1)) > 1e-10) ishell = ishell + 1
               ENDIF
               !ra = MATMUL(cell%amat, ptsh(:, i)) + rc
               tmp_vec = real(ptsh(:, i))
               call dgemv("N", 3, 3, 1.0, cell%amat, 3, tmp_vec, 1, 0.0, ra, 1)
               ra = ra + rc

               a = scale*norm2(ra)
               IF (abs(a) < 1e-12) THEN
                  CYCLE
               ELSE IF (ABS(a - a1) > 1e-10) THEN
                  a1 = a
                  rexp = EXP(-a)
                  IF (ishell <= conv(0)) g(0) = rexp/a &
                                                *(1 + a*11/16*(1 + a*3/11*(1 + a/9)))
                  IF (ishell <= conv(1)) g(1) = rexp/a**2 &
                                                *(1 + a*(1 + a/2*(1 + a*7/24*(1 + a/7))))
                  IF (ishell <= conv(2)) g(2) = rexp/a**3 &
                                                *(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a*3/16*(1 + a/9))))))
                  IF (ishell <= conv(3)) g(3) = rexp/a**4 &
                                                *(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a/5*(1 + a/6*(1 + a/8)))))))
                  IF (ishell <= conv(4)) g(4) = rexp/a**5 &
                                                *(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a/5*(1 + a/6*(1 + a/7*(1 + a/8 &
                                                                                                               *(1 + a/10)))))))))
                  IF (ishell <= conv(5)) g(5) = rexp/a**6 &
                                                *(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a/5*(1 + a/6*(1 + a/7*(1 + a/8*(1 + a/9 &
                                                                                                                        *(1 + a/10))))))))))
                  IF (ishell <= conv(6)) g(6) = rexp/a**7 &
                                                *(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a/5*(1 + a/6*(1 + a/7*(1 + a/8*(1 + a/9 &
                                                                                                                        *(1 + a/10*(1 + a/11*(1 + a/12))))))))))))
                  IF (ishell <= conv(7)) g(7) = rexp/a**8 &
                                                *(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a/5*(1 + a/6*(1 + a/7*(1 + a/8*(1 + a/9 &
                                                                                                                        *(1 + a/10*(1 + a/11*(1 + a/12*(1 + a/13)))))))))))))
                  DO l = 8, maxl
                     IF (ishell <= conv(l)) g(l) = a**(-l - 1)
                  END DO
                  DO l = 0, maxl
                     IF (conv(l) == HUGE(i) .AND. g(l) < convpar(l)/10) conv(l) = ishell + ADDSHELL1
                  END DO
               END IF
               IF (ishell > conv(maxl) .AND. maxl /= 0) maxl = maxl - 1
               call ylm4(maxl, ra, y)
               y = CONJG(y)
               
               DO ikpt = 1, kpts%nkpt
                  DO l = 0, maxl
                     rdum = dot_product(kpts%bk(:, ikpt), ptsh(:, i))
                     cexp = EXP(ImagUnit*tpi_const*rdum)
                     lm = l**2
                     IF (ishell <= conv(l)) THEN
                        cdum = cexp*g(l)
                        DO M = -l, l
                           lm = lm + 1
                           shlp(lm, ikpt) = shlp(lm, ikpt) + cdum*y(lm)
                        END DO
                     END IF
                  END DO
               END DO
            END DO
            structconst(:, ic1, ic2, :) = shlp
         END DO
!         !$OMP END PARALLEL DO
      END DO
#ifdef CPP_MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE, a1, 1, MPI_DOUBLE_PRECISION, MPI_MAX, fmpi%mpi_comm, ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, structconst, size(structconst), MPI_DOUBLE_COMPLEX,MPI_SUM,fmpi%mpi_comm,ierr)
#endif
      call timestop("realspace sum")
      deallocate (ptsh, radsh)
   end subroutine realspace_sum
end module m_structureconstant
