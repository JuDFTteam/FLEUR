!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_outcdn
   implicit none

CONTAINS

   SUBROUTINE outcdn(p, n, na, iv, iflag, jsp, l_potential, stars, vacuum, &
                     sphhar, atoms, sym, cell,   potDen, xdnout,potDenIm,lattvec_index)
      USE m_types
      USE m_constants
      USE m_angle
      USE m_starf, ONLY : starf2,starf3
      USE m_ylm
      USE m_checkdopall


      !--------------------------------------------------------------------------
      ! Calculates the charge density at a given point p(i=1,3).
      !--------------------------------------------------------------------------

      IMPLICIT NONE

      TYPE(t_stars),INTENT(IN)     :: stars
      TYPE(t_vacuum),INTENT(IN)    :: vacuum
      TYPE(t_sphhar),INTENT(IN)    :: sphhar
      TYPE(t_atoms),INTENT(IN)     :: atoms
      TYPE(t_sym),INTENT(IN)       :: sym
      TYPE(t_cell),INTENT(IN)      :: cell
       
      TYPE(t_potden),INTENT(IN)    :: potDen


      ! Scalar Arguments
      INTEGER, INTENT (IN) :: iflag,jsp,n,na,iv
      REAL,    INTENT (OUT) :: xdnout
      ! -odim
      ! +odim

      ! Array Arguments
      REAL,    INTENT (INOUT) :: p(3)

      ! Logical argument
      LOGICAL, INTENT (IN) :: l_potential

      !for responses:
      TYPE(t_potden), OPTIONAL, INTENT(IN)    :: potDenIm
      !TYPE(t_input), OPTIONAL ,   INTENT(IN) :: input

      INTEGER, OPTIONAL, INTENT(IN)       :: lattvec_index(3)

      LOGICAL                         :: l_dfpt

      ! Local scalars
      REAL delta,sx,xd1_r,xd2_r,xx1,xx2,rrr,phi !removed s and made it complex 
      INTEGER i,j,jp3,jr,k,lh,mem,nd,nopa,ivac,ll1,lm ,gzi,m
      real check_mt

      ! Local arrays
      COMPLEX, ALLOCATABLE :: sf2(:),sf3(:),ylm(:)
      REAL rcc(3),x(3)

      REAL :: lattvec(3)
      COMPLEX :: xd1, xd2,s
      l_dfpt = .FALSE. 
      IF (PRESENT(potDenIm)) THEN
         l_dfpt = .TRUE.
      END IF
      ALLOCATE( sf2(stars%ng2),sf3(stars%ng3),ylm((atoms%lmaxd+1)**2))
      ivac=iv

      IF (iflag.NE.1) THEN
         IF (iflag.NE.0) THEN
            ! Interstitial part:
            rcc=matmul(cell%bmat,p)/tpi_const
            IF (l_dfpt) THEN
               !print*,"Im doing starsf"
               CALL starf3(sym%nop, stars%ng3, sym%symor, stars%kv3, sym%mrot, &
                        sym%tau, rcc, sym%invtab, sf3,stars%center)
            ELSE
               CALL starf3(sym%nop, stars%ng3, sym%symor, stars%kv3, sym%mrot, &
               sym%tau, rcc, sym%invtab, sf3)
            END IF
            !WIP make correct for DFPT, same as for MT part
            xdnout=dot_product(real(potDen%pw(:,jsp)*sf3(:)),stars%nstr)
            RETURN

         ENDIF

         ! Vacuum part:
         xdnout = 0.

         
         
            IF (p(3).LT.0.0) THEN
               ivac = vacuum%nvac
               IF (sym%invs) THEN
                  p(1:2) = -p(1:2)
               END IF
               p(3) = abs(p(3))
            END IF
            rcc=matmul(cell%bmat,p)/tpi_const
            IF (l_dfpt) THEN
               CALL starf2(sym%nop2, stars%ng2, stars%kv2, sym%mrot, sym%symor, &
                        sym%tau,rcc,sym%invtab,sf2,stars%center)
            ELSE
               CALL starf2(sym%nop2, stars%ng2, stars%kv2, sym%mrot, sym%symor, &
               sym%tau,rcc,sym%invtab,sf2)
            END IF
            jp3 = (p(3)-cell%z1)/vacuum%delz
            delta = (p(3)-cell%z1)/vacuum%delz - jp3
            ! We count 0 as point 1.
            jp3 = jp3 + 1
            IF (jp3.LT.vacuum%nmz) THEN
               xdnout = REAL(potDen%vac(jp3,1,ivac,jsp)) + &
                        delta*(REAL(potDen%vac(jp3+1,1,ivac,jsp)) - &
                        REAL(potDen%vac(jp3,1,ivac,jsp)))
               IF (jp3.LT.vacuum%nmzxy) THEN
                  xx1 = 0.
                  xx2 = 0.
                  DO  k = 2,stars%ng2
                     xx1 = xx1 + REAL(potDen%vac(jp3,k,ivac,jsp)*sf2(k))* &
                           stars%nstr2(k)
                     xx2 = xx2 + REAL(potDen%vac(jp3+1,k,ivac,jsp)*sf2(k))* &
                           stars%nstr2(k)
                  ENDDO
                  xdnout = xdnout + xx1 + delta* (xx2-xx1)
               END IF
            ELSE
               xdnout = 0.0
            END IF
         ! Vacuum part finished.
         
         RETURN
      ENDIF
      ! MT part:

      nd = sym%ntypsy(na)
      nopa = sym%ngopr(na)
      x = p - atoms%pos(:, na)
      sx = sqrt(sum(x * x))
      IF (nopa.NE.1) THEN
         ! Switch to internal units.
         rcc=matmul(cell%bmat,x)/tpi_const
         ! Rotate into representative.
         p = matmul(sym%mrot(:,:,nopa), rcc)
         ! Switch back to cartesian units.
         x=matmul(cell%amat,p)
      END IF
      DO j = atoms%jri(n),2,-1
         IF (sx.GE.atoms%rmsh(j,n)) EXIT
      ENDDO
      jr = j
      CALL ylm4(atoms%lmax(n),x,ylm)
      xd1 = 0.0
      xd2 = 0.0
      DO  lh = 0, sphhar%nlh(nd)
         ll1 = sphhar%llh(lh,nd) * ( sphhar%llh(lh,nd) + 1 ) + 1
         s = cmplx(0.0)
         DO mem = 1,sphhar%nmem(lh,nd)
            lm = ll1 + sphhar%mlh(mem,lh,nd)
            s = s + sphhar%clnu(mem,lh,nd)*ylm(lm) 
         ENDDO
         IF (l_potential) THEN
            xd1 = xd1 + potDen%mt(jr,lh,n,jsp)*s
            IF (l_dfpt) THEN
               xd1 = xd1 + cmplx(0.0,1.0)*potDenIm%mt(jr,lh,n,jsp)*s
            END IF
         ELSE
            xd1 = xd1 + potDen%mt(jr,lh,n,jsp)*s/(atoms%rmsh(jr,n)**2)
            IF (l_dfpt) THEN
               xd1 = xd1 + cmplx(0.0,1.0)*potDenIm%mt(jr,lh,n,jsp)*s/(atoms%rmsh(jr,n)**2)
            END IF 
         END IF
         IF (jr.EQ.atoms%jri(n)) CYCLE
         IF (l_potential) THEN
            xd2 = xd2 + potDen%mt(jr+1,lh,n,jsp)*s
            IF (l_dfpt) THEN
               xd2 = xd2 + cmplx(0.0,1.0)*potDenIm%mt(jr+1,lh,n,jsp)*s
            END IF
         ELSE
            xd2 = xd2 + potDen%mt(jr+1,lh,n,jsp)*s/(atoms%rmsh(jr+1,n)**2)
            IF (l_dfpt) THEN
               xd2 = xd2 + cmplx(0.0,1.0)*potDenIm%mt(jr+1,lh,n,jsp)*s/(atoms%rmsh(jr+1,n)**2)
            END IF
         END IF
      ENDDO
      If (l_dfpt) THEN
         lattvec = MATMUL(cell%amat,(/-lattvec_index(1),-lattvec_index(2),-lattvec_index(3)/))
         xd1 = xd1*exp(CMPLX(0.0,1.0)*tpi_const*dot_product( stars%center,(/-lattvec_index(1),-lattvec_index(2),-lattvec_index(3)/)))
         xd2 = xd2*exp(CMPLX(0.0,1.0)*tpi_const*dot_product( stars%center,(/-lattvec_index(1),-lattvec_index(2),-lattvec_index(3)/)))
      END IF
      !WIP: introduce some kind of switch in case dfpt real and imag part
      xd1_r = real(xd1)
      xd2_r = real(xd2)
      IF (jr.EQ.atoms%jri(n)) THEN
         xdnout = xd1_r
      ELSE
         xdnout = xd1_r + (xd2_r-xd1_r)* (sx-atoms%rmsh(jr,n))/ &
                  (atoms%rmsh(jr+1,n)-atoms%rmsh(jr,n))
      END IF
   END SUBROUTINE outcdn
END MODULE m_outcdn
