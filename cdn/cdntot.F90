MODULE m_cdntot
#ifdef CPP_MPI
   use mpi
#endif
!     ********************************************************
!     calculate the total charge density in the interstial.,
!     vacuum, and mt regions      c.l.fu
!     ********************************************************
CONTAINS
   SUBROUTINE integrate_cdn(stars,nococonv,atoms,sym,vacuum,input,cell , integrand, &
                                   q, qis, qmt, qvac, qtot, qistot, fmpi)
      ! if called with fmpi variable, distribute the calculation of the pwint 
      ! over fmpi processes in fmpi%mpi_comm
      USE m_intgr, ONLY : intgr3
      USE m_constants
      USE m_qsf
      USE m_pwint
      USE m_types
      USE m_juDFT
      IMPLICIT NONE
      TYPE(t_stars),INTENT(IN)  :: stars
      TYPE(t_nococonv),INTENT(IN):: nococonv
      TYPE(t_atoms),INTENT(IN)  :: atoms
      TYPE(t_sym),INTENT(IN)    :: sym
      TYPE(t_vacuum),INTENT(IN) :: vacuum
      TYPE(t_input),INTENT(IN)  :: input
      TYPE(t_cell),INTENT(IN)   :: cell
       
      TYPE(t_potden),INTENT(IN) :: integrand
      REAL, INTENT(OUT)         :: q(input%jspins), qis(input%jspins), qmt(atoms%ntype,input%jspins),&
                                   qvac(2,input%jspins), qtot, qistot
      TYPE(t_mpi),INTENT(IN),OPTIONAL :: fmpi
      INTEGER                   :: jsp, j, ivac, nz, n, irank, nsize, intstart, intstop, chunk_size, leftover
      REAL                      :: q2(vacuum%nmz), w(4), rht1(vacuum%nmzd,2,input%jspins)
      REAL                      :: sum_over_ng3
      COMPLEX,ALLOCATABLE       :: x(:) !(1:stars%ng3), may be distributed over fmpi ranks
      COMPLEX                   :: w_off
#ifdef CPP_MPI
      INTEGER ierr
#endif
      IF (PRESENT(fmpi)) THEN
         irank = fmpi%irank
         nsize = fmpi%isize
      ELSE
         irank = 0
         nsize = 1
      ENDIF

      qtot = 0.0
      qistot = 0.0
      qvac=0.0
      q=0.0
      qis=0.0
      qmt=0.0
      IF (irank.EQ.0) THEN   
!     -----mt charge
            DO n = 1,atoms%ntype
               DO jsp=1,size(integrand%mt,4)
                   CALL intgr3(integrand%mt(:,0,n,jsp),atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),w(jsp))
               enddo
               if (size(integrand%pw,2)>2) THEN 
                  !this is a noco-calculation
                  if (size(integrand%mt,4)==4) THEN
                     w_off=cmplx(w(3),w(4))
                  else  
                     w_off=0.0
                  endif
                  !rotate into global frame
                  call nococonv%rotdenmat(n,w(1),w(2),w_off,toGlobal=.true.)
               endif
               DO jsp=1,input%jspins       
                  qmt(n, jsp) = w(jsp)*sfp_const
                  q(jsp) = q(jsp) + atoms%neq(n)*qmt(n,jsp)
               ENDDO   
            ENDDO
         IF (input%film) THEN
!     -----vacuum region
            DO jsp = 1,input%jspins        
               DO ivac = 1,vacuum%nvac
                  DO nz = 1,vacuum%nmz
                     rht1(nz,ivac,jsp) =  REAL(integrand%vac(nz,1,ivac,jsp))
                    END DO
                  CALL qsf(vacuum%delz,rht1(1,ivac,jsp),q2,vacuum%nmz,0)
                  qvac(ivac,jsp) = q2(1)*cell%area
                  
                     q(jsp) = q(jsp) + qvac(ivac,jsp)*2./real(vacuum%nvac)
               ENDDO
            enddo
         END IF
      END IF ! irank = 0

      DO jsp = 1,input%jspins
!     -----is region
         chunk_size = stars%ng3/nsize
         leftover = stars%ng3 - chunk_size*nsize
         IF ( leftover > irank ) THEN
            chunk_size = chunk_size + 1
            intstart = irank * chunk_size + 1
         ELSE
            intstart = leftover * (chunk_size+1) + (irank - leftover) * chunk_size + 1
         ENDIF 
         intstop = intstart + chunk_size -1
         ALLOCATE(x(intstart:intstop))
         CALL pwint_all(stars,atoms,sym ,cell,intstart,intstop,x)
         sum_over_ng3 = 0.0
         DO j = intstart,intstop
            sum_over_ng3 = sum_over_ng3 + integrand%pw(j,jsp)*x(j)*stars%nstr(j)
         ENDDO
         DEALLOCATE(x)
#ifdef CPP_MPI
         IF (PRESENT(fmpi)) THEN
            CALL MPI_reduce(sum_over_ng3,qis(jsp),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,fmpi%mpi_comm,ierr)
         ELSE
            qis(jsp) = sum_over_ng3
         ENDIF
#else
         qis(jsp) = sum_over_ng3
#endif
         qistot = qistot + qis(jsp)
         q(jsp) = q(jsp) + qis(jsp)
         qtot = qtot + q(jsp)
      END DO ! loop over spins
   END SUBROUTINE integrate_cdn

   SUBROUTINE integrate_realspace(xcpot, atoms, sym, sphhar, input, &
                                  stars, cell,   vacuum, noco, mt, is, hint)
      use m_types
      use m_mt_tofrom_grid
      use m_pw_tofrom_grid
      use m_constants
      implicit none
      CLASS(t_xcpot), INTENT(inout)   :: xcpot
      TYPE(t_atoms),INTENT(IN)      :: atoms
      TYPE(t_sym), INTENT(in)       :: sym
      TYPE(t_sphhar), INTENT(IN)    :: sphhar
      TYPE(t_input), INTENT(IN)     :: input
      TYPE(t_stars), INTENT(IN)     :: stars
      TYPE(t_cell), INTENT(IN)      :: cell
       
      TYPE(t_vacuum), INTENT(in)    :: vacuum
      TYPE(t_noco), INTENT(in)      :: noco
      real, intent(inout)           :: mt(:,:,:), is(:,:)
      character(len=*), intent(in), optional :: hint
      integer                       :: n_atm, i

      TYPE(t_potden)                :: tmp_potden
      REAL                          :: q(input%jspins), qis(input%jspins), &
                                       qmt(atoms%ntype,input%jspins), qvac(2,input%jspins),&
                                       qtot, qistot

      call tmp_potden%init(stars, atoms, sphhar, vacuum, noco, input%jspins, POTDEN_TYPE_DEN)
      call init_mt_grid(input%jspins, atoms, sphhar, xcpot%needs_grad(), sym)
      do n_atm =1,atoms%ntype
         call mt_from_grid(atoms, sym, sphhar, n_atm, input%jspins, mt(:,:,n_atm), &
                           tmp_potden%mt(:,0:,n_atm,:))

         do i=1,atoms%jri(n_atm)
            tmp_potden%mt(i,:,n_atm,:) = tmp_potden%mt(i,:,n_atm,:) * atoms%rmsh(i,n_atm)**2
         enddo
      enddo
      call finish_mt_grid()

      call init_pw_grid(stars, sym, cell,xcpot)
      call pw_from_grid( stars, is, tmp_potden%pw)  !THIS CODE SEEMS TO BE BROKEN!!
      call finish_pw_grid()

      call judft_error("Bug, integrate_realspace in cdntot")
      !call integrate_cdn(stars,atoms,sym,vacuum,input,cell , tmp_potden, &
      !                             q, qis, qmt, qvac, qtot, qistot)

      call print_cdn_inte(q, qis, qmt, qvac, qtot, qistot, hint)
   END SUBROUTINE integrate_realspace

   SUBROUTINE cdntot(stars,nococonv,atoms,sym,vacuum,input,cell ,&
                     den,l_printData,qtot,qistot,fmpi,l_par)

      USE m_types
      USE m_juDFT
      IMPLICIT NONE

!     .. Scalar Arguments ..
      TYPE(t_stars),INTENT(IN)  :: stars
      TYPE(t_nococonv),INTENT(IN):: nococonv
      TYPE(t_atoms),INTENT(IN)  :: atoms
      TYPE(t_sym),INTENT(IN)    :: sym
      TYPE(t_vacuum),INTENT(IN) :: vacuum
      TYPE(t_input),INTENT(IN)  :: input
       
      TYPE(t_cell),INTENT(IN)   :: cell
      TYPE(t_potden),INTENT(IN) :: den
      LOGICAL,INTENT(IN)        :: l_printData,l_par
      REAL,INTENT(OUT)          :: qtot,qistot
      TYPE(t_mpi),INTENT(IN)    :: fmpi

!     .. Local Scalars ..
      REAL q(input%jspins),qis(input%jspins),w,mtCharge
!     ..
!     .. Local Arrays ..
      REAL qmt(atoms%ntype,input%jspins),qvac(2,input%jspins)

      CALL timestart("cdntot")
      IF (l_par) THEN
         CALL integrate_cdn(stars,nococonv,atoms,sym,vacuum,input,cell , den, &
                                   q, qis, qmt, qvac, qtot, qistot, fmpi)
      ELSE
         CALL integrate_cdn(stars,nococonv,atoms,sym,vacuum,input,cell , den, &
                                   q, qis, qmt, qvac, qtot, qistot)
      ENDIF

      IF (fmpi%irank.EQ.0) CALL cdntot_writings(atoms,vacuum,input,l_printData,q,qis,qmt,qvac,qtot)
      CALL timestop("cdntot")
   END SUBROUTINE cdntot

   SUBROUTINE cdntot_writings(atoms,vacuum,input,l_printData,q,qis,qmt,qvac,qtot)

      USE m_constants
      USE m_types
      USE m_juDFT
      USE m_xmlOutput
      IMPLICIT NONE

!     .. Scalar Arguments ..
      TYPE(t_atoms),INTENT(IN)  :: atoms
      TYPE(t_vacuum),INTENT(IN) :: vacuum
      TYPE(t_input),INTENT(IN)  :: input
      LOGICAL,INTENT(IN)        :: l_printData
      REAL,INTENT(IN)           :: q(input%jspins),qis(input%jspins)
      REAL,INTENT(IN)           :: qmt(atoms%ntype,input%jspins),qvac(2,input%jspins)
      REAL,INTENT(IN)           :: qtot

!     .. Local Scalars ..
      REAL mtCharge
      INTEGER i,jsp,n
!     ..
!     .. Local Arrays ..
      INTEGER, ALLOCATABLE :: lengths(:,:)
      CHARACTER(LEN=20) :: attributes(6), names(6)


      IF (input%film) THEN
         ALLOCATE(lengths(4+vacuum%nvac,2))
      ELSE
         ALLOCATE(lengths(4,2))
      END IF

      DO jsp = 1,input%jspins
         WRITE (oUnit,FMT=8000) jsp,q(jsp),qis(jsp), (qmt(n,jsp),n=1,atoms%ntype)
            IF (input%film) WRITE (oUnit,FMT=8010) (i,qvac(i,jsp),i=1,vacuum%nvac)
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
         WRITE (oUnit,FMT=8020) qtot
         IF(l_printData) THEN
            CALL writeXMLElementFormPoly('totalCharge',(/'value'/),(/qtot/),reshape((/5,20/),(/1,2/)))
         END IF
8000     FORMAT (/,10x,'total charge for spin',i3,'=',f12.6,/,10x,&
               'interst. charge =   ',f12.6,/,&
               (10x,'mt charge=          ',4f12.6,/))
8010     FORMAT (10x,'vacuum ',i2,'  charge=  ',f12.6)
8020     FORMAT (/,10x,'total charge  =',f12.6)

   END SUBROUTINE cdntot_writings

   SUBROUTINE print_cdn_inte(q, qis, qmt, qvac, qtot, qistot, hint)
      use  ieee_arithmetic
      implicit none
      REAL, INTENT(in)                       :: q(:), qis(:), qmt(:,:), qvac(:,:), qtot, qistot
      character(len=*), intent(in), optional :: hint
      integer                                :: n_mt


      if(present(hint)) write (*,*) "DEN of ", hint
      write (*,*) "q   = ", q
      write (*,*) "qis = ", qis

      write (*,*) "qmt"
      do n_mt = 1,size(qmt, dim=1)
         write (*,*) "mt = ", n_mt, qmt(n_mt,:)
      enddo

      if(.not. any(ieee_is_nan(qvac))) then
         write (*, *) "qvac",    qvac
      endif
      write (*, *) "qtot",    qtot
      write (*, *) "qis_tot", qistot
      write (*, *) "-------------------------"
   END SUBROUTINE print_cdn_inte
END MODULE m_cdntot
