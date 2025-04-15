! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_l_like
   !!calculate the "l-like charge"
contains

   subroutine print_l_like_charge(spin0,atoms,radfun,denmatrix,itype)
            !!Use the density matrix in denCoeffs to calculate the l-like charge
            !!Output to out and out.xml
      use m_types_denmatrix
      use m_types_radfun
      use m_types
      use m_xmlOutput
      implicit none
      integer,intent(in)           :: spin0
      type(t_radfun), intent(IN)   :: radfun
      type(t_atoms), intent(IN)    :: atoms
      type(t_denmatrix), intent(in)    :: denmatrix(spin0:,spin0:,:)
      integer,intent(in):: itype
      real    :: qmtl(0:atoms%lmax(itype))
      character(LEN=20) :: attributes(6)
      integer:: l,ispin

      if (itype==1) write (oUnit, FMT=8000) !! write header if this is the first atom
8000  format(/, 5x, 'l-like charge', /, t6, 'atom', t15, 's', t24, 'p', &
              t33, 'd', t42, 'f', t51, 'total')

      do ispin=lbound(denmatrix,1),ubound(denmatrix,1)        
         qmtl=denmatrix(ispin,ispin,itype)%l_like_charge(radfun,ispin,size(qmtl)-1)/atoms%neq(itype)        
      
         write (oUnit, FMT=8100) itype, (qmtl(l), l=0, 3), sum(qmtl)
8100        format(' -->', i3, 2x, 4f9.5, 2x, f9.5)
         attributes = ''
         write (attributes(1), '(i0)') itype
         write (attributes(2), '(f12.7)') sum(qmtl)
         write (attributes(3), '(f12.7)') qmtl(0)
         write (attributes(4), '(f12.7)') qmtl(1)
         write (attributes(5), '(f12.7)') qmtl(2)
         write (attributes(6), '(f12.7)') qmtl(3)
         call writeXMLElementForm('mtCharge', (/'atomType', 'total   ', 's       ', 'p       ', 'd       ', 'f       '/), &
         attributes(:6),reshape((/8, 5, 1, 1, 1, 1, 6, 12, 12, 12, 12, 12/), (/6, 2/)))
      enddo      
   end subroutine
end module m_l_like
