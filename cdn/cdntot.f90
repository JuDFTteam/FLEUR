      MODULE m_cdntot
!     ********************************************************
!     calculate the total charge density in the interstial.,
!     vacuum, and mt regions      c.l.fu
!     ********************************************************
      CONTAINS
      SUBROUTINE cdntot(&
     &                  stars,atoms,sym,&
     &                  vacuum,input,cell,oneD,&
     &                  qpw,rho,rht,&
     &                  qtot,qistot)

      USE m_intgr, ONLY : intgr3
      USE m_constants
      USE m_qsf
      USE m_pwint
      USE m_types
      use m_juDFT
      IMPLICIT NONE
!     ..
!     .. Scalar Arguments ..
      TYPE(t_stars),INTENT(IN) :: stars
      TYPE(t_atoms),INTENT(IN) :: atoms
      TYPE(t_sym),INTENT(IN)   :: sym
      TYPE(t_vacuum),INTENT(IN):: vacuum
      TYPE(t_input),INTENT(IN) :: input
      TYPE(t_oneD),INTENT(IN)  :: oneD
      TYPE(t_cell),INTENT(IN)  :: cell
      REAL,    INTENT (OUT):: qtot,qistot
!     ..
!     .. Array Arguments ..
      COMPLEX, INTENT (IN) :: qpw(stars%n3d,input%jspins)
      REAL,    INTENT (IN) :: rho(:,0:,:,:) !(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd,input%jspins)
      REAL,    INTENT (IN) :: rht(vacuum%nmzd,2,input%jspins)
!-odim
!+odim
!     ..
!     .. Local Scalars ..
    ! COMPLEX x
      COMPLEX x(stars%ng3)
      REAL q,qis,w
      INTEGER i,ivac,j,jspin,n,nz
!     ..
!     .. Local Arrays ..
      REAL qmt(atoms%ntypd),qvac(2),q2(vacuum%nmz),rht1(vacuum%nmzd,2,input%jspins)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC real
!     ..
!
      CALL timestart("cdntot")
      qtot = 0.e0
      qistot = 0.e0
      DO 40 jspin = 1,input%jspins
         q = 0.e0
!     -----mt charge
         CALL timestart("MT")
         DO 10 n = 1,atoms%ntype
            CALL intgr3(rho(:,0,n,jspin),atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),w)
            qmt(n) = w*sfp_const
            q = q + atoms%neq(n)*qmt(n)
   10    CONTINUE
         CALL timestop("MT")
!     -----vacuum region
         IF (input%film) THEN
            DO 20 ivac = 1,vacuum%nvac
               DO nz = 1,vacuum%nmz
                  IF (oneD%odi%d1) THEN
                     rht1(nz,ivac,jspin) = (cell%z1+(nz-1)*vacuum%delz)*&
     &                    rht(nz,ivac,jspin)
                  ELSE
                     rht1(nz,ivac,jspin) =  rht(nz,ivac,jspin)
                  END IF
               END DO
               CALL qsf(vacuum%delz,rht1(1,ivac,jspin),q2,vacuum%nmz,0)
               qvac(ivac) = q2(1)*cell%area
               IF (.NOT.oneD%odi%d1) THEN
                  q = q + qvac(ivac)*2./real(vacuum%nvac)
               ELSE
                  q = q + cell%area*q2(1)
               END IF
   20       CONTINUE
         END IF
!     -----is region
         qis = 0.
!         DO 30 j = 1,nq3
!            CALL pwint(
!     >                 k1d,k2d,k3d,n3d,ntypd,natd,nop,invtab,odi,
!     >                 ntype,neq,volmts,taual,z1,vol,volint,
!     >                 symor,tau,mrot,rmt,sk3,bmat,ig2,ig,
!     >                 kv3(1,j),
!     <                 x)
!            qis = qis + qpw(j,jspin)*x*nstr(j)
!   30    CONTINUE
         CALL pwint_all(&
     &                 stars,atoms,sym,oneD,&
     &                 cell,&
     &                 x)
         DO j = 1,stars%ng3
             qis = qis + qpw(j,jspin)*x(j)*stars%nstr(j)
         ENDDO
         qistot = qistot + qis
         q = q + qis
         WRITE (6,FMT=8000) jspin,q,qis, (qmt(n),n=1,atoms%ntype)
         IF (input%film) WRITE (6,FMT=8010) (i,qvac(i),i=1,vacuum%nvac)
         WRITE (16,FMT=8000) jspin,q,qis, (qmt(n),n=1,atoms%ntype)
         IF (input%film) WRITE (16,FMT=8010) (i,qvac(i),i=1,vacuum%nvac)
         qtot = qtot + q
   40 CONTINUE
      WRITE (6,FMT=8020) qtot
      WRITE (16,FMT=8020) qtot
 8000 FORMAT (/,10x,'total charge for spin',i3,'=',f12.6,/,10x,&
     &       'interst. charge =   ',f12.6,/,&
     &       (10x,'mt charge=          ',4f12.6,/))
 8010 FORMAT (10x,'vacuum ',i2,'  charge=  ',f12.6)
 8020 FORMAT (/,10x,'total charge  =',f12.6)

      CALL timestop("cdntot")
      END SUBROUTINE cdntot
      END MODULE m_cdntot
