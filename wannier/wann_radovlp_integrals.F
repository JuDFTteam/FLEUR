!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_radovlp_integrals
      contains
      subroutine wann_radovlp_integrals(
     >               rmsh,jri,dx,radial1,radial2,
     <               integral)
c*************************************************
c     Calculate the radial integral needed for the
c     overlap.
c              Frank Freimuth
c*************************************************
      USE m_intgr, ONLY : intgr3      
      implicit none
      real, intent(in)    :: rmsh(:)
      integer,intent(in)  :: jri
      real, intent(in)    :: dx
      real, intent(in)    :: radial1(:,:)
      real, intent(in)    :: radial2(:,:)
      real, intent(out)   :: integral

      integer :: ir
      real    :: rr
      real    :: integrand(jri)
      integer :: slc

      REAL     :: difcub
      EXTERNAL :: difcub

      integrand=0.0
      do slc=1,2
         DO ir = 1,jri
            integrand(ir) = integrand(ir)+
     &         radial1(ir,slc)*radial2(ir,slc)
         END DO 
      enddo !slc

      CALL intgr3(
     &  integrand,rmsh(:),dx,jri,integral)      

      end subroutine wann_radovlp_integrals
      end module m_wann_radovlp_integrals
