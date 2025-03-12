! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_l_like
   !!calculate the "l-like charge"
contains

   subroutine print_l_like_charge(input, atoms,usdus, hub1inp, denCoeffs,hub1data)
            !!Use the density matrix in denCoeffs to calculate the l-like charge
            !!Output to out and out.xml
      use m_constants
      use m_types
      use m_xmlOutput
      implicit none
      type(t_input), intent(IN)    :: input
      type(t_atoms), intent(IN)    :: atoms
      type(t_usdus), intent(in)    :: usdus
      type(t_denCoeffs), intent(IN)    :: denCoeffs
      TYPE(t_hub1inp),           INTENT(IN)    :: hub1inp
      TYPE(t_hub1data), OPTIONAL, INTENT(INOUT) :: hub1data

      integer :: itype, ispin, l, llp, lp, lo, lop, i_hia, i_exc
      real    :: qmtl(0:atoms%lmaxd)
      complex :: ctemp
      character(LEN=20) :: attributes(6)

      write (oUnit, FMT=8000)
8000  format(/, 5x, 'l-like charge', /, t6, 'atom', t15, 's', t24, 'p', &
              t33, 'd', t42, 'f', t51, 'total')

      do itype = 1, atoms%ntype
         do ispin = 1, input%jspins
            !l-decomposed density for each atom type
            do l = 0, atoms%lmax(itype)
               llp = l*(atoms%lmax(itype) + 1) + l
               qmtl(l) = real(denCoeffs%nmt_coeff(llp, 0, itype, 0, 0, ispin, ispin) +&
                              denCoeffs%nmt_coeff(llp, 0, itype, 1, 1, ispin, ispin) &
                              *usdus%ddn(l, itype, ispin))/atoms%neq(itype)
            end do
            !LO contribution
            do lo = 1, atoms%nlo(itype)
               l = atoms%llo(lo, itype)
               ctemp = (denCoeffs%mt_ulo_coeff(lo, itype, 0, iSpin, ispin) + denCoeffs%mt_lou_coeff(lo, itype, 0, iSpin, ispin)) &
                         & *usdus%uulon(lo, itype, ispin) &
                   & + (denCoeffs%mt_ulo_coeff(lo, itype, 1, ispin, ispin) + denCoeffs%mt_lou_coeff(lo, itype, 1, ispin, ispin)) &
                         & *usdus%dulon(lo, itype, ispin)
               qmtl(l) = qmtl(l) + real(ctemp)/atoms%neq(itype)
               do lop = 1, atoms%nlo(itype)
                  if (atoms%llo(lop, itype) .eq. l) then
                     ctemp = denCoeffs%mt_lolo_coeff(lop, lo, itype, ispin, ispin)*usdus%uloulopn(lop, lo, itype, ispin)
                     qmtl(l) = qmtl(l) + real(ctemp)/atoms%neq(itype)
                  end if
               end do
            end do


            !Get the magnetic moment for the shells where we defined additional exchange splittings for DFT+Hubbard 1
            IF(PRESENT(hub1data)) THEN
               DO i_hia = 1, atoms%n_hia
                  IF(atoms%lda_u(atoms%n_u+i_hia)%atomType/=itype) CYCLE
                  DO i_exc = 1, hub1inp%n_exc(i_hia)
                     hub1data%mag_mom(i_hia,i_exc) = hub1data%mag_mom(i_hia,i_exc) + (-1)**(ispin-1) &
                                                   *  qmtl(hub1inp%exc_l(i_hia,i_exc))
                  END DO
               END DO
            END IF

            write (oUnit, FMT=8100) itype, (qmtl(l), l=0, 3), sum(qmtl)
8100        format(' -->', i3, 2x, 4f9.5, 2x, f9.5)
            attributes = ''
            write (attributes(1), '(i0)') itype
            write (attributes(2), '(f12.7)') sum(qmtl)
            write (attributes(3), '(f12.7)') qmtl(0)
            write (attributes(4), '(f12.7)') qmtl(1)
            write (attributes(5), '(f12.7)') qmtl(2)
            write (attributes(6), '(f12.7)') qmtl(3)
  call writeXMLElementForm('mtCharge', (/'atomType', 'total   ', 's       ', 'p       ', 'd       ', 'f       '/), attributes(:6), &
                                     reshape((/8, 5, 1, 1, 1, 1, 6, 12, 12, 12, 12, 12/), (/6, 2/)))
         end do
      end do

   end subroutine
end module m_l_like
