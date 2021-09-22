module m_slater

   use m_types
   use m_constants
   use m_juDFT
   use m_differ
   use m_intgr
   use m_xmlOutput

   implicit none

contains
   subroutine slater(input, jspin, atoms, vr, l_write, slater_parameters)

      ! Calculate the slater integrals for the orbitals specified
      ! to be corrected with LDA+OP
      !
      ! Adapted from core/cored.F90
      ! Author: Henning Janssen 2021

      type(t_input),     intent(in) :: input
      type(t_atoms),     intent(in) :: atoms
      integer,           intent(in) :: jspin
      real,              intent(in) :: vr(:, :)
      logical, optional, intent(in) :: l_write
      real, allocatable, optional, intent(out) :: slater_parameters(:,:)

      real    :: eig, fj, fl, fn, t2
      real    :: d, dxx, rn, rnot, z, t1, rr, r, c
      integer :: i, n, ncmsh, l, ierr, kk
      integer :: i_opc, atomType, ipm

      real :: vrd(atoms%msh), f(0:3, 2, atoms%n_opc)
      real :: a(atoms%msh), b(atoms%msh), ain(atoms%msh), ahelp(atoms%msh)
      real :: lambda
      character(LEN=20) :: attributes(6)

      logical :: l_write_arg

      l_write_arg = .false.
      if (present(l_write)) l_write_arg = l_write

      if (l_write_arg) call openXMLElementNoAttributes('slaterIntegrals')

      c = c_light(1.0)

      do i_opc = 1, atoms%n_opc

         l = atoms%lda_opc(i_opc)%l
         n = atoms%lda_opc(i_opc)%n
         atomType = atoms%lda_opc(i_opc)%atomType
         lambda = atoms%lda_opc(i_opc)%lambda

         z = atoms%zatom(atomType)
         dxx = atoms%dx(atomType)
         rnot = atoms%rmsh(1, atomType)
         d = exp(atoms%dx(atomType))
         ncmsh = nint(log((atoms%rmt(atomType) + 10.0)/rnot)/dxx + 1)
         ncmsh = min(ncmsh, atoms%msh)
         rn = rnot*(d**(ncmsh - 1))
         if (l_write_arg) write (oUnit, fmt=8000) z, rnot, dxx, atoms%jri(atomType)
         vrd(:atoms%jri(atomType)) = vr(:atoms%jri(atomType), atomType)

         !
         if (input%l_core_confpot) then
            !--->    linear extension of the potential with slope t1 / a.u.
            t1 = 0.125
            t1 = max((vrd(atoms%jri(atomType)) - vrd(atoms%jri(atomType) - 1)*d)* &
                     d/(atoms%rmt(atomType)**2*(d - 1)), t1)
            t2 = vrd(atoms%jri(atomType))/atoms%rmt(atomType) - atoms%rmt(atomType)*t1
            rr = atoms%rmt(atomType)
         else
            t2 = vrd(atoms%jri(atomType))/(atoms%jri(atomType) - ncmsh)
         end if
         if (atoms%jri(atomType) < ncmsh) then
            do i = atoms%jri(atomType) + 1, ncmsh
               if (input%l_core_confpot) then
                  rr = d*rr
                  vrd(i) = rr*(t2 + rr*t1)
                  !               vrd(i) = 2*vrd(jri(jatom)) - rr*( t2 + rr*t1 )
               else
                  vrd(i) = vrd(atoms%jri(atomType)) + t2*(i - atoms%jri(atomType))
               end if
               !
            end do
         end if

         fl = real(l)
         fn = real(n)
         do ipm = 1, 2
            fj = fl + (ipm - 1.5)

            eig = -2*(z/(fn + fl))**2
            call differ(fn, fl, fj, c, z, dxx, rnot, rn, d, ncmsh, vrd, eig, a, b, ierr)
            if (l_write_arg) write (oUnit, fmt=8010) fn, fl, fj, eig

            if (ierr /= 0) call juDFT_error("error in slater routine", calledby="slater")

            do kk = 0, 2*l, 2  ! F0 for s, F0 + F2 for p etc.
               !lambda = 3.5449*sqrt((l+1.0)/100)     ! screening (TF) sqrt(4pi N(ef))
               !IF (kk.GT.0) lambda = 2*lambda
               r = rnot
               do i = 1, ncmsh
                  ain(i) = a(i)**2*r**(-kk - 1)      ! prepare inner integrand
                  if (kk == 0) then
                     ain(i) = ain(i)*exp(-r*lambda)
                  end if
                  r = r*d
               end do
               call intgr1(ain, rnot, dxx, ncmsh, ahelp)

               r = rnot
               do i = 1, ncmsh - 1
                  ain(i) = a(i)**2*r**kk*(ahelp(ncmsh) - ahelp(i))
                  if (kk == 0) then
                     ain(i) = ain(i)*exp(r*lambda)
                  end if
                  r = r*d
               end do
               call intgr0(ain, rnot, dxx, ncmsh - 1, f(kk/2, ipm, i_opc))
            end do
         end do
         f(0:l, :, i_opc) = f(0:l, :, i_opc)*2
         if (l_write_arg) then
            write (oUnit, fmt=8020) n, l, jspin, lambda
            write (oUnit, '(12x,i3,f7.1,4f20.10)') l, fj, (f(kk, 1, i_opc)*hartree_to_ev_const, kk=0, l)
            write (oUnit, '(12x,i3,f7.1,4f20.10)') l, fj, (f(kk, 2, i_opc)*hartree_to_ev_const, kk=0, l)
         end if

         if (l_write_arg) then
            attributes = ''
            write (attributes(1), '(i0)') atomType
            write (attributes(2), '(i0)') nint(z)
            write (attributes(3), '(i0)') jspin
            write (attributes(4), '(i0)') l
            write (attributes(5), '(i0)') n
            write (attributes(6), '(f10.7)') lambda
            call openXMLElementForm('slaterIntegral', (/'atomType     ', 'atomicNumber ', 'spin         ', 'l            ', &
                                                        'n            ', 'screening    '/), &
                                    attributes, reshape((/8, 12, 4, 1, 1, 9, 6, 3, 1, 1, 1, 18/), (/6, 2/)))
            do ipm = 1, 2
               fj = fl + (ipm - 1.5)
               attributes = ''
               write (attributes(1), '(i0)') n
               write (attributes(2), '(i0)') nint(fl)
               write (attributes(3), '(f4.1)') fj
               write (attributes(4), '(a)') 'eV'
               call writeXMLElementPoly('state', (/'n    ', 'l    ', 'j    ','units'/), &
                                        attributes(1:4), contentList=f(0:l, ipm, i_opc)*hartree_to_ev_const)
            end do
            call closeXMLElement('slaterIntegral')
         end if
      end do

      if (l_write_arg) call closeXMLElement('slaterIntegrals')

      if (present(slater_parameters)) then
         allocate(slater_parameters(0:lmaxU_const,atoms%n_opc), source=0.0)
         do i_opc = 1, atoms%n_opc
            l = atoms%lda_opc(i_opc)%l
         !Is simply averaging the right approach??
            slater_parameters(0:l, i_opc) =  (f(0:l,1,i_opc) + f(0:l,2,i_opc))/2.0
         enddo
      endif

8000  format(/, /, 10x, 'z=', f4.0, 5x, 'r(1)=', e14.6, 5x, 'dx=', f9.6, 5x,&
          &       'm.t.index=', i4, /, 15x, 'n', 4x, 'l', 5x, 'j', 4x, 'energy')
8010  format(12x, 2f5.0, f6.1, f10.4, f12.4)
8020  format(/, /, 10x, 'Slater integrals: n=', i2, 5x, 'l=', i2, 5x, 'spin=', i2, 5x, 'screening=', f9.6, 5x,&
          &       /, 14x, 'l', 4x, 'j', 9x, 'F(0,2,4,6)')
   end subroutine slater
end module m_slater
