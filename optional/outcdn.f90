      MODULE m_outcdn
      use m_constants
      USE m_types
!     ********************************************************
!     calculates the charge density at given point p(i=1,3)
!     ********************************************************
      CONTAINS
      !SUBROUTINE outcdn(&
     !&                  p,n,na,iv,iflag,jsp,l_potential,stars,&
     !&                  vacuum,sphhar,atoms,sym,cell,oneD,&
     !&                  qpw,rhtxy,rho,rht,&
     !&                  xdnout)
 SUBROUTINE outcdn(&
     &                  p,n,na,iv,iflag,jsp,l_potential,stars,&
     &                  vacuum,sphhar,atoms,sym,cell,oneD,&
     &                  potDen,&
     &                  xdnout)
!
      use m_constants
      USE m_angle
      USE m_starf, ONLY : starf2,starf3
      USE m_ylm
      IMPLICIT NONE
!
!      TYPE(t_sliceplot),INTENT(IN) :: sliceplot TODO:Remove
      TYPE(t_stars),INTENT(IN)     :: stars
      TYPE(t_vacuum),INTENT(IN)    :: vacuum
      TYPE(t_sphhar),INTENT(IN)    :: sphhar
      TYPE(t_atoms),INTENT(IN)     :: atoms
      TYPE(t_sym),INTENT(IN)       :: sym
      TYPE(t_cell),INTENT(IN)      :: cell
      TYPE(t_oneD),INTENT(IN)      :: oneD
      TYPE(t_potden),INTENT(IN)    :: potDen


!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: iflag,jsp,n,na,iv
      REAL,    INTENT (OUT) :: xdnout
!-odim
!+odim
!     ..
!     .. Array Arguments ..
      COMPLEX qpw(:,:) !(stars%ng3,input%jspins)
      COMPLEX  rhtxy(:,:,:,:) !(vacuum%nmzxyd,oneD%odi%n2d-1,2,input%jspins)
      REAL  rho(:,0:,:,:) !(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins)
      REAL  rht(:,:,:) !(vacuum%nmzd,2,input%jspins)
      REAL,    INTENT (INOUT) :: p(3)
!     ..
!     .. Local Scalars ..
      REAL delta,s,sx,xd1,xd2,xx1,xx2,rrr,phi
      INTEGER i,j,jp3,jr,k,lh,mem,nd,nopa,ivac,ll1,lm ,gzi,m

!     .. Logical Argumens  ..
      LOGICAL, INTENT (IN) :: l_potential 
      
!     ..
!     .. Local Arrays ..
      COMPLEX sf2(stars%ng2),sf3(stars%ng3),ylm((atoms%lmaxd+1)**2)
      REAL rcc(3),x(3)
      
      !Temp. Variable Assignment TODO: Can be removed as soon as the potDen datatype has been implemented in the whole file.
      qpw=potDen%pw
      rhtxz=potDen%vacxy
      rho=potDen%mt
      rht=potDen%vacz
       

      ivac=iv
     
      if (iflag.ne.1) THEN
      if (iflag.ne.0) THEN
!     ---> interstitial part
      !CALL cotra1(p(1),rcc,cell%bmat)
      rcc=matmul(cell%bmat,p)/tpi_const
      CALL starf3(&
     &            sym%nop,stars%ng3,sym%symor,stars%kv3,sym%mrot,sym%tau,rcc,sym%invtab,&
     &            sf3)
!
      xdnout=dot_product(real(qpw(:,jsp)*sf3(:)),stars%nstr)
      RETURN
!     ---> vacuum part
      ENDIF
      xdnout = 0.
!-odim
      IF (oneD%odi%d1) THEN
         rrr = sqrt( p(1)**2 + p(2)**2 )
         phi = angle(p(1),p(2))
         jp3 = (rrr-cell%z1)/vacuum%delz
         delta = (rrr-cell%z1)/vacuum%delz - jp3
!*     we count 0 as point 1
         jp3 = jp3 + 1
         IF (jp3.LT.vacuum%nmz) THEN
            xdnout = rht(jp3,ivac,jsp) + delta*&
     &           (rht(jp3+1,ivac,jsp)-rht(jp3,ivac,jsp))
            IF (jp3.LT.vacuum%nmzxy) THEN
               xx1 = 0.
               xx2 = 0.
               DO  k = 2,oneD%odi%nq2
                  m = oneD%odi%kv(2,k)
                  gzi = oneD%odi%kv(1,k)
                  xx1 = xx1 + real(rhtxy(jp3,k-1,ivac,jsp)*&
     &                 exp(ImagUnit*m*phi)*exp(ImagUnit*gzi*cell%bmat(3,3)*p(3)))*&
     &                 oneD%odi%nst2(k)
                  xx2 = xx2 + real(rhtxy(jp3+1,k-1,ivac,jsp)*&
     &                 exp(ImagUnit*m*phi)*exp(ImagUnit*gzi*cell%bmat(3,3)*p(3)))*&
     &                 oneD%odi%nst2(k)
            ENDDO
               xdnout = xdnout + xx1 + delta* (xx2-xx1)
            END IF
         ELSE
            xdnout = 0.0
         END IF

      ELSE
!+odim      
         IF (p(3).LT.0.0) THEN
            ivac = vacuum%nvac
            IF (sym%invs) THEN
               p(1:2) = -p(1:2)
            END IF
            p(3) = abs(p(3))
         END IF
         !CALL cotra1(p,rcc,cell%bmat)
         rcc=matmul(cell%bmat,p)/tpi_const
         CALL starf2(&
     &            sym%nop2,stars%ng2,stars%kv3,sym%mrot,sym%symor,sym%tau,rcc,sym%invtab,&
     &            sf2)
!
         jp3 = (p(3)-cell%z1)/vacuum%delz
         delta = (p(3)-cell%z1)/vacuum%delz - jp3
!*     we count 0 as point 1
         jp3 = jp3 + 1
         IF (jp3.LT.vacuum%nmz) THEN
             xdnout = rht(jp3,ivac,jsp) + delta*&
     &               (rht(jp3+1,ivac,jsp)-rht(jp3,ivac,jsp))
            IF (jp3.LT.vacuum%nmzxy) THEN
               xx1 = 0.
               xx2 = 0.
              DO  k = 2,stars%ng2
               xx1 = xx1 + real(rhtxy(jp3,k-1,ivac,jsp)*sf2(k))*stars%nstr2(k)
               xx2 = xx2 + real(rhtxy(jp3+1,k-1,ivac,jsp)*sf2(k))*&
     &               stars%nstr2(k)
   enddo
              xdnout = xdnout + xx1 + delta* (xx2-xx1)
            END IF
         ELSE
            xdnout = 0.0
         END IF
!----> vacuum finishes
      ENDIF

      RETURN
      ENDIF
!     ----> m.t. part
      
      nd = atoms%ntypsy(na)
      nopa = atoms%ngopr(na)
      IF (oneD%odi%d1) nopa = oneD%ods%ngopr(na)
      sx = 0.0
      DO  i = 1,3
         x(i) = p(i) - atoms%pos(i,na)
         sx = sx + x(i)*x(i)
   enddo
      sx = sqrt(sx)
      IF (nopa.NE.1) THEN
!... switch to internal units
         !CALL cotra1(x,rcc,cell%bmat)
         rcc=matmul(cell%bmat,x)/tpi_const
!... rotate into representative
         DO  i = 1,3
            p(i) = 0.
            DO  j = 1,3
              IF (.NOT.oneD%odi%d1) THEN
               p(i) = p(i) + sym%mrot(i,j,nopa)*rcc(j)
              ELSE
               p(i) = p(i) + oneD%ods%mrot(i,j,nopa)*rcc(j)
              END IF
   enddo
   enddo
!... switch back to cartesian units
         !CALL cotra0(p,x,cell%amat)
         x=matmul(cell%amat,p)
      END IF
      DO j = atoms%jri(n),2,-1
         IF (sx.GE.atoms%rmsh(j,n)) EXIT
      ENDDO
      jr = j
      CALL ylm4(&
     &          atoms%lmax(n),x,&
     &          ylm)
      xd1 = 0.0
      xd2 = 0.0
      DO  lh = 0, sphhar%nlh(nd)
         ll1 = sphhar%llh(lh,nd) * ( sphhar%llh(lh,nd) + 1 ) + 1
         s = 0.0
         DO mem = 1,sphhar%nmem(lh,nd)
           lm = ll1 + sphhar%mlh(mem,lh,nd)
           s = s + real( sphhar%clnu(mem,lh,nd)*ylm(lm) )
         ENDDO
         IF (l_potential) THEN
            xd1 = xd1 + rho(jr,lh,n,jsp)*s
         ELSE
            xd1 = xd1 + rho(jr,lh,n,jsp)*s/ (atoms%rmsh(jr,n)*atoms%rmsh(jr,n))
         END IF
         IF (jr.EQ.atoms%jri(n)) CYCLE
         IF (l_potential) THEN
            xd2 = xd2 + rho(jr+1,lh,n,jsp)*s
         ELSE
            xd2 = xd2 + rho(jr+1,lh,n,jsp)*s/&
     &            (atoms%rmsh(jr+1,n)*atoms%rmsh(jr+1,n))
         END IF
   ENDDO
      IF (jr.EQ.atoms%jri(n)) THEN
         xdnout = xd1
      ELSE
         xdnout = xd1 + (xd2-xd1) *&
     &                  (sx-atoms%rmsh(jr,n)) / (atoms%rmsh(jr+1,n)-atoms%rmsh(jr,n))
      END IF
 8000 FORMAT (2f10.6)
!
      RETURN
      END SUBROUTINE outcdn
      END MODULE m_outcdn
