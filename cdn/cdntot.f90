      MODULE m_cdntot
!     ********************************************************
!     calculate the total charge density in the interstial.,
!     vacuum, and mt regions      c.l.fu
!     ********************************************************
      CONTAINS
      SUBROUTINE cdntot(&
     &                  stars,atoms,sym,&
     &                  vacuum,input,cell,oneD,&
     &                  qpw,rho,rht,l_printData,&
     &                  qtot,qistot)

      USE m_intgr, ONLY : intgr3
      USE m_constants
      USE m_qsf
      USE m_pwint
      USE m_types
      USE m_juDFT
      USE m_convol
      USE m_xmlOutput
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
      LOGICAL,INTENT(IN)       :: l_printData
      REAL,    INTENT (OUT):: qtot,qistot
!     ..
!     .. Array Arguments ..
      COMPLEX, INTENT (IN) :: qpw(stars%ng3,input%jspins)
      REAL,    INTENT (IN) :: rho(:,0:,:,:) !(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins)
      REAL,    INTENT (IN) :: rht(vacuum%nmzd,2,input%jspins)
!-odim
!+odim
!     ..
!     .. Local Scalars ..
    ! COMPLEX x
      COMPLEX x(stars%ng3)
      REAL q,qis,w,mtCharge
      INTEGER i,ivac,j,jspin,n,nz
!     ..
!     .. Local Arrays ..
      REAL qmt(atoms%ntype),qvac(2),q2(vacuum%nmz),rht1(vacuum%nmzd,2,input%jspins)
      INTEGER, ALLOCATABLE :: lengths(:,:)
      CHARACTER(LEN=20) :: attributes(6), names(6)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC real
!     ..
!
      IF (input%film) THEN
         ALLOCATE(lengths(4+vacuum%nvac,2))
      ELSE
         ALLOCATE(lengths(4,2))
      END IF
      CALL timestart("cdntot")
      qtot = 0.e0
      qistot = 0.e0
      DO jspin = 1,input%jspins
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
         IF (.not.judft_was_Argument("-oldfix")) THEN
            CALL convol(stars,x,qpw(:,jspin),stars%ufft)
            qis = x(1)*cell%omtil
         ELSE
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
         endif
         qistot = qistot + qis
         q = q + qis
         WRITE (6,FMT=8000) jspin,q,qis, (qmt(n),n=1,atoms%ntype)
         IF (input%film) WRITE (6,FMT=8010) (i,qvac(i),i=1,vacuum%nvac)
         WRITE (16,FMT=8000) jspin,q,qis, (qmt(n),n=1,atoms%ntype)
         IF (input%film) WRITE (16,FMT=8010) (i,qvac(i),i=1,vacuum%nvac)
         mtCharge = SUM(qmt(1:atoms%ntype) * atoms%neq(1:atoms%ntype))
         names(1) = 'spin'         ; WRITE(attributes(1),'(i0)') jspin       ; lengths(1,1)=4  ; lengths(1,2)=1
         names(2) = 'total'        ; WRITE(attributes(2),'(f14.7)') q        ; lengths(2,1)=5  ; lengths(2,2)=14
         names(3) = 'interstitial' ; WRITE(attributes(3),'(f14.7)') qis      ; lengths(3,1)=12 ; lengths(3,2)=14
         names(4) = 'mtSpheres'    ; WRITE(attributes(4),'(f14.7)') mtCharge ; lengths(4,1)=9  ; lengths(4,2)=14
         IF(l_printData) THEN
            IF(input%film) THEN
               DO i = 1, vacuum%nvac
                  WRITE(names(4+i),'(a6,i0)') 'vacuum', i
                  WRITE(attributes(4+i),'(f14.7)') qvac(i)
                  lengths(4+i,1)=7
                  lengths(4+i,2)=14
               END DO
               CALL writeXMLElementFormPoly('spinDependentCharge',names(1:4+vacuum%nvac),&
                                            attributes(1:4+vacuum%nvac),lengths)
            ELSE
               CALL writeXMLElementFormPoly('spinDependentCharge',names(1:4),attributes(1:4),lengths)
            END IF
         END IF
         qtot = qtot + q
      END DO ! loop over spins
      DEALLOCATE (lengths)
      WRITE (6,FMT=8020) qtot
      WRITE (16,FMT=8020) qtot
      IF(l_printData) THEN
         CALL writeXMLElementFormPoly('totalCharge',(/'value'/),(/qtot/),reshape((/5,20/),(/1,2/)))
      END IF
 8000 FORMAT (/,10x,'total charge for spin',i3,'=',f12.6,/,10x,&
     &       'interst. charge =   ',f12.6,/,&
     &       (10x,'mt charge=          ',4f12.6,/))
 8010 FORMAT (10x,'vacuum ',i2,'  charge=  ',f12.6)
 8020 FORMAT (/,10x,'total charge  =',f12.6)

      CALL timestop("cdntot")
      END SUBROUTINE cdntot
      END MODULE m_cdntot
