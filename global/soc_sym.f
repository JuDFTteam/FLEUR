!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_socsym
      use m_juDFT
!-----------------------------------------------------------------------!
! tests the compatibility of the symmetry elements with the SQA defined !
! by theta and phi in case of spin-orbit coupling.                gb`02 !
!-----------------------------------------------------------------------!
      CONTAINS
      SUBROUTINE soc_sym(
     >                   nop,mrot,theta,phi,amat,
     <                   error)
      
      USE m_constants, ONLY : pimach
      USE m_inv3
      IMPLICIT NONE

      INTEGER, INTENT (IN)  :: nop, mrot(3,3,nop)
      REAL,    INTENT (IN)  :: theta, phi, amat(3,3)
      LOGICAL, INTENT (OUT) :: error(nop)

      INTEGER iop
      REAL st,ct,sp,cp,stt,ctt,pih,test,nn1n,nn2n
      REAL sqa(3),n1(3),n2(3),nn1(3),nn2(3),ssqa(3)
      REAL ainv(3,3),rrot(3,3),rrot1(3,3)

      pih= 0.5*pimach()
      CALL inv3(amat,ainv,test)
      st = sin(theta) ; sp = sin(phi) ; stt = sin(theta + pih)
      ct = cos(theta) ; cp = cos(phi) ; ctt = cos(theta + pih)
!
! --> find two vectors n1,n2 normal to the spin-quantizatiopn axis (sqa)
!
      sqa(1) = st*cp ; n1(1) = stt*cp
      sqa(2) = st*sp ; n1(2) = stt*sp
      sqa(3) = ct    ; n1(3) = ctt

      CALL cross(sqa,n1,
     <           n2)
!dbg      write(*,*) n1
!dbg      write(*,*) n2
!
! --> loop over symmetry elements
!
      error(:) = .false.
      DO iop = 1, nop
       
!        CALL matmul3(mrot(1,1,iop),ainv,rrot1)
        rrot1=matmul(mrot(:,:,iop),ainv)
!        CALL matmul2(amat,rrot1,rrot)
        rrot=matmul(amat,rrot1)
!
! ----> rotate n1 and n2 by symmetry element and form the cross-product
!       of the rotated vectors (nn1,nn2) --> ssqa
!
         nn1=matmul(rrot,n1)
         nn2=matmul(rrot,n2)

        CALL cross(nn1,nn2,
     <             ssqa)
!
! ----> if ssqa is identical with sqa accept this symmetry element
!
        test = sqa(1)*ssqa(1) + sqa(2)*ssqa(2) + sqa(3)*ssqa(3) 
        IF (abs(test-1.0).GT.0.00001) THEN
          error(iop) = .true.
          WRITE (6,100) iop,test
        ENDIF
      ENDDO
 100  FORMAT ('Symmetry element no.',i3,' incompatible with SQA',f20.10)
      
      IF ( ANY(error(:)) ) THEN
!         CALL juDFT_warn
!     +        ("symmetry incompatible with SOC - Spin Quantization Axis"
!     +        ,calledby ="soc_sym")
         WRITE (6,*) 'symmetry incompatible with SOC - Spin Quant Axis'
      ENDIF
      END SUBROUTINE soc_sym
!---------------------------------------------------------------
      SUBROUTINE cross(
     >                 a,b,
     <                 c)

      IMPLICIT NONE
      REAL, INTENT  (IN) :: a(3),b(3)
      REAL, INTENT (OUT) :: c(3)

      c(1) = a(2) * b(3) - a(3) * b(2)
      c(2) = a(3) * b(1) - a(1) * b(3)
      c(3) = a(1) * b(2) - a(2) * b(1)

      END SUBROUTINE cross

      END MODULE m_socsym

