MODULE m_cdntot
!     ********************************************************
!     calculate the total charge density in the interstial.,
!     vacuum, and mt regions      c.l.fu
!     ********************************************************
CONTAINS
   SUBROUTINE cdntot(mpi,stars,atoms,sym,vacuum,input,cell,oneD,&
                     den,l_printData,qtot,qistot)

      USE m_intgr, ONLY : intgr3
      USE m_constants
      USE m_qsf
      USE m_pwint
      USE m_types
      USE m_juDFT
      USE m_convol
      USE m_xmlOutput
      IMPLICIT NONE

!     .. Scalar Arguments ..
      TYPE(t_mpi),INTENT(IN)    :: mpi
      TYPE(t_stars),INTENT(IN)  :: stars
      TYPE(t_atoms),INTENT(IN)  :: atoms
      TYPE(t_sym),INTENT(IN)    :: sym
      TYPE(t_vacuum),INTENT(IN) :: vacuum
      TYPE(t_input),INTENT(IN)  :: input
      TYPE(t_oneD),INTENT(IN)   :: oneD
      TYPE(t_cell),INTENT(IN)   :: cell
      TYPE(t_potden),INTENT(IN) :: den
      LOGICAL,INTENT(IN)        :: l_printData
      REAL,INTENT(OUT)          :: qtot,qistot

!     .. Local Scalars ..
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
         DO n = 1,atoms%ntype
            CALL intgr3(den%mt(:,0,n,jspin),atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),w)
            qmt(n) = w*sfp_const
            q = q + atoms%neq(n)*qmt(n)
         ENDDO
         CALL timestop("MT")
!     -----vacuum region
         IF (input%film) THEN
            DO ivac = 1,vacuum%nvac
               DO nz = 1,vacuum%nmz
                  IF (oneD%odi%d1) THEN
                     rht1(nz,ivac,jspin) = (cell%z1+(nz-1)*vacuum%delz)*&
                                           den%vacz(nz,ivac,jspin)
                  ELSE
                     rht1(nz,ivac,jspin) =  den%vacz(nz,ivac,jspin)
                  END IF
               END DO
               CALL qsf(vacuum%delz,rht1(1,ivac,jspin),q2,vacuum%nmz,0)
               qvac(ivac) = q2(1)*cell%area
               IF (.NOT.oneD%odi%d1) THEN
                  q = q + qvac(ivac)*2./real(vacuum%nvac)
               ELSE
                  q = q + cell%area*q2(1)
               END IF
            ENDDO
         END IF
!     -----is region
         qis = 0.

         CALL pwint_all(stars,atoms,sym,oneD,cell,1,stars%ng3,x)
         DO j = 1,stars%ng3
            qis = qis + den%pw(j,jspin)*x(j)*stars%nstr(j)
         ENDDO

         qistot = qistot + qis
         q = q + qis
         WRITE (6,FMT=8000) jspin,q,qis, (qmt(n),n=1,atoms%ntype)
         IF (input%film) WRITE (6,FMT=8010) (i,qvac(i),i=1,vacuum%nvac)
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
      IF(l_printData) THEN
         CALL writeXMLElementFormPoly('totalCharge',(/'value'/),(/qtot/),reshape((/5,20/),(/1,2/)))
      END IF
8000  FORMAT (/,10x,'total charge for spin',i3,'=',f12.6,/,10x,&
               'interst. charge =   ',f12.6,/,&
               (10x,'mt charge=          ',4f12.6,/))
8010  FORMAT (10x,'vacuum ',i2,'  charge=  ',f12.6)
8020  FORMAT (/,10x,'total charge  =',f12.6)

      CALL timestop("cdntot")
   END SUBROUTINE cdntot
END MODULE m_cdntot
