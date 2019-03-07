MODULE m_cdntot
!     ********************************************************
!     calculate the total charge density in the interstial.,
!     vacuum, and mt regions      c.l.fu
!     ********************************************************
CONTAINS
   SUBROUTINE integrate_grid(xcpot, stars, atoms, sym, vacuum, input, cell, oneD, sphhar,noco,&
                             is_inte, mt_inte, &
                             q, qis, qmt, qvac, qtot, qistot)
      USE m_pw_tofrom_grid
      USE m_mt_tofrom_grid
      USE m_types
      USE m_constants
      !USE m_types_xcpot
      IMPLICIT NONE
      CLASS(t_xcpot),INTENT(IN) :: xcpot
      TYPE(t_stars),INTENT(IN)  :: stars
      TYPE(t_atoms),INTENT(IN)  :: atoms
      TYPE(t_sym),INTENT(IN)    :: sym
      TYPE(t_vacuum),INTENT(IN) :: vacuum
      TYPE(t_input),INTENT(IN)  :: input
      TYPE(t_cell),INTENT(IN)   :: cell
      TYPE(t_oneD),INTENT(IN)   :: oneD
      TYPE(t_sphhar), INTENT(IN):: sphhar
      TYPE(t_noco), INTENT(INOUT)   :: noco
      TYPE(t_grid), INTENT(in)  :: is_inte, mt_inte(:)
      REAL, INTENT(out)  :: q(input%jspins), qis(input%jspins), qmt(atoms%ntype,input%jspins),&
                            qvac(2,input%jspins), qtot, qistot

      TYPE(t_potden)     :: integrand
      
      TYPE(t_grid)       :: is_inte_mut
      INTEGER            :: n,i 
      
      call init_pw_grid(xcpot, stars, sym, cell)
      call init_mt_grid(input%jspins, atoms, sphhar, xcpot, sym)

      is_inte_mut = is_inte

      !allocate potden type
      call integrand%init(stars, atoms, sphhar, vacuum, noco, input%jspins, POTDEN_TYPE_DEN)
      allocate(integrand%pw_w, mold=integrand%pw)

      !put is in potden-basis
      call pw_from_grid(xcpot, stars,.True., is_inte_mut%grid, integrand%pw, integrand%pw_w)

      !put mt in potden-basis
      do n = 1,atoms%ntype
         call mt_from_grid(atoms,sphhar,n,input%jspins,mt_inte(n)%grid,integrand%mt(:,0:,n,:))
         do i =1,atoms%jri(n)
            integrand%mt(i,0:,n,:) = integrand%mt(i,0:,n,:) * atoms%rmsh(i,n)**2
         enddo
      enddo

      ! integrate my integrand
      call integrate_cdn(stars, atoms, sym, vacuum, input, cell, oneD, integrand,&
                        q, qis, qmt, qvac, qtot, qistot)
      call finish_pw_grid()
      call finish_mt_grid()
   END SUBROUTINE integrate_grid

   SUBROUTINE integrate_cdn(stars,atoms,sym,vacuum,input,cell,oneD, integrand, &
                                   q, qis, qmt, qvac, qtot, qistot)
      USE m_intgr, ONLY : intgr3
      USE m_constants
      USE m_qsf
      USE m_pwint
      USE m_types
      USE m_juDFT
      USE m_convol
      IMPLICIT NONE
      TYPE(t_stars),INTENT(IN)  :: stars
      TYPE(t_atoms),INTENT(IN)  :: atoms
      TYPE(t_sym),INTENT(IN)    :: sym
      TYPE(t_vacuum),INTENT(IN) :: vacuum
      TYPE(t_input),INTENT(IN)  :: input
      TYPE(t_cell),INTENT(IN)   :: cell
      TYPE(t_oneD),INTENT(IN)   :: oneD
      TYPE(t_potden),INTENT(IN) :: integrand
      REAL, INTENT(out)         :: q(input%jspins), qis(input%jspins), qmt(atoms%ntype,input%jspins),&
                                   qvac(2,input%jspins), qtot, qistot
      INTEGER                   :: jsp, j, ivac, nz, n
      REAL                      :: q2(vacuum%nmz), w, rht1(vacuum%nmzd,2,input%jspins)
      COMPLEX                   :: x(stars%ng3)
      
      qtot = 0.0
      qistot = 0.0
      DO jsp = 1,input%jspins
         q(jsp) = 0.0
!     -----mt charge
         CALL timestart("MT")
         DO n = 1,atoms%ntype
            CALL intgr3(integrand%mt(:,0,n,jsp),atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),w)
            qmt(n, jsp) = w*sfp_const
            q(jsp) = q(jsp) + atoms%neq(n)*qmt(n,jsp)
         ENDDO
         CALL timestop("MT")
!     -----vacuum region
         IF (input%film) THEN
            DO ivac = 1,vacuum%nvac
               DO nz = 1,vacuum%nmz
                  IF (oneD%odi%d1) THEN
                     rht1(nz,ivac,jsp) = (cell%z1+(nz-1)*vacuum%delz)*&
                                           integrand%vacz(nz,ivac,jsp)
                  ELSE
                     rht1(nz,ivac,jsp) =  integrand%vacz(nz,ivac,jsp)
                  END IF
               END DO
               CALL qsf(vacuum%delz,rht1(1,ivac,jsp),q2,vacuum%nmz,0)
               qvac(ivac,jsp) = q2(1)*cell%area
               IF (.NOT.oneD%odi%d1) THEN
                  q(jsp) = q(jsp) + qvac(ivac,jsp)*2./real(vacuum%nvac)
               ELSE
                  q(jsp) = q(jsp) + cell%area*q2(1)
               END IF
            ENDDO
         END IF
!     -----is region
         qis(jsp) = 0.

         CALL pwint_all(stars,atoms,sym,oneD,cell,1,stars%ng3,x)
         DO j = 1,stars%ng3
            qis(jsp) = qis(jsp) + integrand%pw(j,jsp)*x(j)*stars%nstr(j)
         ENDDO

         qistot = qistot + qis(jsp)
         q(jsp) = q(jsp) + qis(jsp)
         qtot = qtot + q(jsp)
      END DO ! loop over spins
   END SUBROUTINE integrate_cdn

   SUBROUTINE cdntot(stars,atoms,sym,vacuum,input,cell,oneD,&
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
      REAL q(input%jspins),qis(input%jspins),w,mtCharge
      INTEGER i,ivac,j,jsp,n,nz
!     ..
!     .. Local Arrays ..
      REAL qmt(atoms%ntype,input%jspins),qvac(2,input%jspins)
      INTEGER, ALLOCATABLE :: lengths(:,:)
      CHARACTER(LEN=20) :: attributes(6), names(6)
      
      CALL timestart("cdntot")
      call integrate_cdn(stars,atoms,sym,vacuum,input,cell,oneD, den, &
                                   q, qis, qmt, qvac, qtot, qistot)
 
      IF (input%film) THEN
         ALLOCATE(lengths(4+vacuum%nvac,2))
      ELSE
         ALLOCATE(lengths(4,2))
      END IF

      DO jsp = 1,input%jspins
         WRITE (6,FMT=8000) jsp,q(jsp),qis(jsp), (qmt(n,jsp),n=1,atoms%ntype)
         IF (input%film) WRITE (6,FMT=8010) (i,qvac(i,jsp),i=1,vacuum%nvac)
         mtCharge = SUM(qmt(1:atoms%ntype,jsp) * atoms%neq(1:atoms%ntype))
         names(1) = 'spin'         ; WRITE(attributes(1),'(i0)')    jsp      ; lengths(1,1)=4  ; lengths(1,2)=1
         names(2) = 'total'        ; WRITE(attributes(2),'(f14.7)') q(jsp)   ; lengths(2,1)=5  ; lengths(2,2)=14
         names(3) = 'interstitial' ; WRITE(attributes(3),'(f14.7)') qis(jsp) ; lengths(3,1)=12 ; lengths(3,2)=14
         names(4) = 'mtSpheres'    ; WRITE(attributes(4),'(f14.7)') mtCharge ; lengths(4,1)=9  ; lengths(4,2)=14
         IF(l_printData) THEN
            IF(input%film) THEN
               DO i = 1, vacuum%nvac
                  WRITE(names(4+i),'(a6,i0)') 'vacuum', i
                  WRITE(attributes(4+i),'(f14.7)') qvac(i,jsp)
                  lengths(4+i,1)=7
                  lengths(4+i,2)=14
               END DO
               CALL writeXMLElementFormPoly('spinDependentCharge',names(1:4+vacuum%nvac),&
                                            attributes(1:4+vacuum%nvac),lengths)
            ELSE
               CALL writeXMLElementFormPoly('spinDependentCharge',names(1:4),attributes(1:4),lengths)
            END IF
         END IF
      END DO ! loop over spins
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
