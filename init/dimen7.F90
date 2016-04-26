      MODULE m_dimen7
      use m_juDFT
      CONTAINS
      SUBROUTINE dimen7(&
     &                  input,sym,stars,&
     &                  atoms,sphhar,dimension,vacuum,&
     &                  obsolete,kpts,oneD,hybrid,Jij,cell)

!
! This program reads the input files of the flapw-programm (inp & kpts)
! and creates a file 'fl7para' that contains dimensions 
! for the main flapw-programm.
!

      USE m_localsym
      USE m_socsym
      USE m_sssym
      USE m_spg2set
      USE m_constants
      USE m_rwinp
      USE m_inpnoco
      USE m_julia
      USE m_od_kptsgen
      USE m_types
      USE m_firstglance
      USE m_inv3
      USE m_rwsymfile
      USE m_strgndim
      USE m_convndim
      USE m_inpeigdim
      USE m_kptgen_hybrid
      USE m_ylm
      IMPLICIT NONE
!
! dimension-parameters for flapw:
!
      TYPE(t_input),INTENT(INOUT)   :: input
      TYPE(t_sym),INTENT(INOUT)     :: sym
      TYPE(t_stars),INTENT(INOUT)   :: stars 
      TYPE(t_atoms),INTENT(INOUT)   :: atoms
      TYPE(t_sphhar),INTENT(INOUT)  :: sphhar
      TYPE(t_dimension),INTENT(INOUT) :: dimension
      TYPE(t_vacuum),INTENT(INOUT)   :: vacuum
      TYPE(t_obsolete),INTENT(INOUT) :: obsolete
      TYPE(t_kpts),INTENT(INOUT)     :: kpts
      TYPE(t_oneD),INTENT(INOUT)     :: oneD
      TYPE(t_hybrid),INTENT(INOUT)   :: hybrid
      TYPE(t_Jij),INTENT(INOUT)      :: Jij
      TYPE(t_cell),INTENT(INOUT)     :: cell
 
      TYPE(t_noco)      :: noco
      TYPE(t_sliceplot) :: sliceplot
      TYPE(t_banddos)   :: banddos
      TYPE(t_xcpot)     :: xcpot

!
!
!-------------------------------------------------------------------
! ..  Local Scalars ..
      REAL   :: thetad,xa,epsdisp,epsforce ,rmtmax,arltv1,arltv2,arltv3   
      REAL   :: s,r,d ,idsprs,scale
      INTEGER :: ok,ilo,n,nstate,i,j,na,n1,n2,jrc,nopd,symfh
      INTEGER :: nmop(3) ,nmopq(3)
      
      CHARACTER(len=1) :: rw
      CHARACTER(len=4) :: namex 
      CHARACTER(len=7) :: symfn
      CHARACTER(len=12):: relcor
      LOGICAL  ::l_kpts,l_qpts,l_inpexist,l_tmp(2)
! ..
      REAL    :: a1(3),a2(3),a3(3)  
      REAL    :: q(3)

      CHARACTER(len=3), ALLOCATABLE :: noel(:)
      LOGICAL, ALLOCATABLE :: error(:) 
     
      INTEGER ntp1,ii
      INTEGER, ALLOCATABLE :: lmx1(:), nq1(:), nlhtp1(:)

!     added for HF and hybrid functionals
      LOGICAL               ::  l_gamma=.false.

      EXTERNAL prp_xcfft_box,parawrite
!     ..
      
    
!---> First, check whether an inp-file exists
!
      INQUIRE (file='inp',exist=l_inpexist)
      IF (.not.l_inpexist) THEN
         CALL juDFT_error("no inp- or input-file found!",calledby&
     &        ="dimen7")
      ENDIF
!
!---> determine ntype,nop,natd,nwdd,nlod and layerd
!
      CALL first_glance(atoms%ntype,sym%nop,atoms%natd,obsolete%nwdd,atoms%nlod,vacuum%layerd,&
                        input%itmax,l_kpts,l_qpts,l_gamma,kpts%nkpt,nmop,jij%nqpt,nmopq)
      atoms%ntypd=atoms%ntype
      atoms%nlod = max(atoms%nlod,1)

      ALLOCATE (&
     & atoms%lmax(atoms%ntype),atoms%ntypsy(atoms%natd),atoms%neq(atoms%ntype),atoms%nlhtyp(atoms%ntype),&
     & atoms%rmt(atoms%ntype),atoms%zatom(atoms%ntype),atoms%jri(atoms%ntype),atoms%dx(atoms%ntype), &
     & atoms%nlo(atoms%ntype),atoms%llo(atoms%nlod,atoms%ntype),atoms%nflip(atoms%ntype),atoms%bmu(atoms%ntype),&
     & noel(atoms%ntype),vacuum%izlay(vacuum%layerd,2),atoms%ncst(atoms%ntype),atoms%lnonsph(atoms%ntype),&
     & atoms%taual(3,atoms%natd),atoms%pos(3,atoms%natd),&
     & atoms%nz(atoms%ntype),atoms%relax(3,atoms%ntype),&
     & atoms%l_geo(atoms%ntype),noco%soc_opt(atoms%ntype+2),noco%alph(atoms%ntype),noco%beta(atoms%ntype),&
     & atoms%lda_u(atoms%ntype),noco%l_relax(atoms%ntype),jij%l_magn(atoms%ntype),jij%M(atoms%ntype),&
     & jij%magtype(atoms%ntype),jij%nmagtype(atoms%ntype),noco%b_con(2,atoms%ntype),&
     & sphhar%clnu(1,1,1),sphhar%nlh(1),sphhar%llh(1,1),sphhar%nmem(1,1),sphhar%mlh(1,1,1),&
     & hybrid%select1(4,atoms%ntype),hybrid%lcutm1(atoms%ntype),hybrid%select2(4,atoms%ntype),hybrid%lcutm2(atoms%ntype),&
     & hybrid%lcutwf(atoms%ntype), STAT=ok)
!
!---> read complete input and calculate nvacd,llod,lmaxd,jmtd,neigd and 
!
      CALL rw_inp('r',&
     &            atoms,obsolete,vacuum,input,stars,sliceplot,banddos,&
     &                  cell,sym,xcpot,noco,jij,oneD,hybrid,kpts,&
     &                  noel,namex,relcor,a1,a2,a3,scale)

!---> pk non-collinear
!---> read the angle and spin-spiral information from nocoinp
      noco%qss = 0.0
      noco%l_ss = .false.
      IF (noco%l_noco) THEN 
        IF (.not.jij%l_J) THEN
         CALL inpnoco(atoms,input,vacuum,jij,noco)
        ELSE
         noco%l_ss= .true.
        ENDIF
      ENDIF

      vacuum%nvacd = 2
      IF (sym%zrfs .OR. sym%invs .OR. oneD%odd%d1) vacuum%nvacd = 1
      atoms%llod  = 0
      atoms%lmaxd = 0
      atoms%jmtd  = 0
      rmtmax      = 0.0
      dimension%neigd = 0
      dimension%nstd  = maxval(atoms%ncst)
      atoms%lmaxd = maxval(atoms%lmax)
      atoms%jmtd  = maxval(atoms%jri)
      rmtmax      = maxval(atoms%rmt)
      DO n = 1,atoms%ntype
        DO ilo = 1,atoms%nlo(n)
!+apw
          IF (atoms%llo(ilo,n).LT.0) THEN
             atoms%llo(ilo,n) = -atoms%llo(ilo,n) - 1
#ifndef CPP_APW
             atoms%llo(ilo,n) = mod(atoms%llo(ilo,n),10)
#endif
          ELSE
             dimension%neigd = dimension%neigd + atoms%neq(n)*(2*abs(atoms%llo(ilo,n)) +1)
          ENDIF
!-apw
          atoms%llod = max(abs(atoms%llo(ilo,n)),atoms%llod)
        ENDDO
        nstate = 4
        IF ((atoms%nz(n).GE.21.AND.atoms%nz(n).LE.29) .OR. &
     &      (atoms%nz(n).GE.39.AND.atoms%nz(n).LE.47) .OR.&
     &      (atoms%nz(n).GE.57.AND.atoms%nz(n).LE.79)) nstate = 9
        IF ((atoms%nz(n).GE.58.AND.atoms%nz(n).LE.71) .OR.&
     &      (atoms%nz(n).GE.90.AND.atoms%nz(n).LE.103)) nstate = 16
        dimension%neigd = dimension%neigd + nstate*atoms%neq(n)
!
      ENDDO
      CALL ylmnorm_init(atoms%lmaxd)
!      IF (mod(lmaxd,2).NE.0) lmaxd = lmaxd + 1
      IF (2*dimension%neigd.LT.input%zelec) THEN
        WRITE(6,*) dimension%neigd,' states estimated in dimen7 ...'
        dimension%neigd = NINT(0.75*input%zelec)
        WRITE(6,*) 'changed dimension%neigd to ',dimension%neigd
      ENDIF
      IF (noco%l_soc .and. (.not. noco%l_noco)) dimension%neigd=2*dimension%neigd 
      IF (noco%l_soc .and. noco%l_ss) dimension%neigd=(3*dimension%neigd)/2  
       ! not as accurate, but saves much time

      rmtmax = rmtmax*stars%gmax
      CALL convn_dim(rmtmax,dimension%ncvd)
!
! determine core mesh
!
      dimension%msh = 0
      DO n = 1,atoms%ntype
         r = atoms%rmt(n)
         d = exp(atoms%dx(n))
         jrc = atoms%jri(n)
         DO WHILE (r < atoms%rmt(n) + 20.0)
            jrc = jrc + 1
            r = r*d
         ENDDO
         dimension%msh = max( dimension%msh, jrc ) 
      ENDDO
!
! ---> now, set the lattice harmonics, determine nlhd
!
      cell%amat(:,1) = a1(:)*scale
      cell%amat(:,2) = a2(:)*scale
      cell%amat(:,3) = a3(:)*scale
      CALL inv3(cell%amat,cell%bmat,cell%omtil)
      IF (input%film) cell%omtil = cell%omtil/cell%amat(3,3)*vacuum%dvac
!-odim
      IF (oneD%odd%d1) cell%omtil = cell%amat(3,3)*pimach()*(vacuum%dvac**2)/4.
!+odim
      cell%bmat=tpi_const*cell%bmat
    
      na = 0
      DO n = 1,atoms%ntype
        DO n1 = 1,atoms%neq(n)
            na = na + 1
            IF (input%film) atoms%taual(3,na) = atoms%taual(3,na)/a3(3)
            atoms%pos(:,na) = matmul(cell%amat,atoms%taual(:,na))
        ENDDO
        atoms%zatom(n) = real( atoms%nz(n) )
      ENDDO
      ALLOCATE (sym%mrot(3,3,sym%nop),sym%tau(3,sym%nop))
      IF (sym%namgrp.EQ.'any ') THEN
         nopd = sym%nop ; rw = 'R'
         symfh = 94 ; symfn = 'sym.out'
         CALL rw_symfile(rw,symfh,symfn,nopd,cell%bmat,sym%mrot,sym%tau,sym%nop,sym%nop2,sym%symor)
      ELSE
         CALL spg2set(sym%nop,sym%zrfs,sym%invs,sym%namgrp,cell%latnam,sym%mrot,sym%tau,sym%nop2,sym%symor)
      ENDIF
      sphhar%ntypsd = 0
      IF (.NOT.oneD%odd%d1) THEN
        CALL local_sym(atoms%lmaxd,atoms%lmax,sym%nop,sym%mrot,sym%tau,&
                       atoms%natd,atoms%ntype,atoms%neq,cell%amat,cell%bmat,&
                       atoms%taual,sphhar%nlhd,sphhar%memd,sphhar%ntypsd,.true.,&
                       atoms%nlhtyp,atoms%ntypsy,sphhar%nlh,sphhar%llh,&
                       sphhar%nmem,sphhar%mlh,sphhar%clnu)
!-odim
      ELSEIF (oneD%odd%d1) THEN
        ntp1 = atoms%natd
        ALLOCATE (nq1(ntp1),lmx1(ntp1),nlhtp1(ntp1))
        ii = 1
        nq1=1
        DO i = 1,atoms%ntype
          DO j = 1,atoms%neq(i)
            lmx1(ii) = atoms%lmax(i)
            ii = ii + 1
          END DO
        END DO
        CALL local_sym(atoms%lmaxd,lmx1,sym%nop,sym%mrot,sym%tau,&
              atoms%natd,ntp1,nq1,cell%amat,cell%bmat,atoms%taual,&
              sphhar%nlhd,sphhar%memd,sphhar%ntypsd,.true.,nlhtp1,&
              atoms%ntypsy,sphhar%nlh,sphhar%llh,sphhar%nmem,&
              sphhar%mlh,sphhar%clnu)        
        ii = 1
        DO i = 1,atoms%ntype
          atoms%nlhtyp(i) = nlhtp1(ii)
          ii = ii + atoms%neq(i)
        END DO
        DEALLOCATE (nq1,lmx1,nlhtp1)
      END IF
!+odim
!
! Check if symmetry is compatible with SOC or SSDW
!
      IF (noco%l_soc .and. (.not.noco%l_noco)) THEN  
        ! test symmetry for spin-orbit coupling
        ALLOCATE ( error(sym%nop) )
        CALL soc_sym(sym%nop,sym%mrot,noco%theta,noco%phi,cell%amat,error)
        IF ( ANY(error(:)) ) THEN
          WRITE(*,fmt='(1x)')
          WRITE(*,fmt='(A)')&
     &     'Symmetry incompatible with SOC spin-quantization axis ,'  
          WRITE(*,fmt='(A)')&
     &     'do not perform self-consistent calculations !'    
          WRITE(*,fmt='(1x)')
          IF ( input%eonly .or. (noco%l_soc.and.noco%l_ss) .or. input%gw.ne.0 ) THEN  ! .or. .
            CONTINUE 
          ELSE 
            IF (input%itmax>1) THEN
               CALL juDFT_error("symmetry & SOC",calledby&
     &              ="dimen7")
            ENDIF 
          ENDIF 
        ENDIF           
        DEALLOCATE ( error )
      ENDIF
!--- J<
      IF(.not.jij%l_J) THEN
!--- J>
      IF (noco%l_ss) THEN  ! test symmetry for spin-spiral
        ALLOCATE ( error(sym%nop) )
        CALL ss_sym(sym%nop,sym%mrot,noco%qss,error)
        IF ( ANY(error(:)) )  CALL juDFT_error("symmetry & SSDW",&
     &       calledby="dimen7")
        DEALLOCATE ( error )
      ENDIF
!--- J<
      ENDIF
!--- J>
!
! Dimensioning of the stars
!
      IF (input%film.OR.(sym%namgrp.ne.'any ')) THEN
         CALL strgn1_dim(stars%gmax,cell%bmat,sym%invs,sym%zrfs,sym%mrot,&
                    sym%tau,sym%nop,sym%nop2,stars%k1d,stars%k2d,stars%k3d,&
                    stars%n3d,stars%n2d,oneD%odd)

      ELSE
         CALL strgn2_dim(stars%gmax,cell%bmat,sym%invs,sym%zrfs,sym%mrot,&
                    sym%tau,sym%nop,stars%k1d,stars%k2d,stars%k3d,&
                    stars%n3d,stars%n2d)
         oneD%odd%n2d = stars%n2d
         oneD%odd%nq2 = stars%n2d
         oneD%odd%nop = sym%nop
      ENDIF

      IF ( xcpot%gmaxxc .le. 10.0**(-6) ) THEN
         WRITE (6,'(" xcpot%gmaxxc=0 : xcpot%gmaxxc=stars%gmax choosen as default",&
     &              " value")')
         WRITE (6,'(" concerning memory, you may want to choose",&
     &              " a smaller value for stars%gmax")')
         xcpot%gmaxxc=stars%gmax
      END IF

      CALL prp_xcfft_box(xcpot%gmaxxc,cell%bmat,stars%kxc1d,stars%kxc2d,stars%kxc3d)
!
! k-point generator provides kpts-file, if it's missing:
!
      IF (.not.l_kpts) THEN
       IF (.NOT.oneD%odd%d1) THEN
         IF (jij%l_J) THEN
         n1=sym%nop
         n2=sym%nop2
         sym%nop=1
         sym%nop2=1
         CALL julia(&
     &              sym,cell,input,noco,banddos,&
     &              kpts,.false.)
         sym%nop=n1
         sym%nop2=n2
         ELSE IF(l_gamma .and. banddos%ndir .eq. 0) THEN
         CALL kptgen_hybrid(kpts%nmop(1),kpts%nmop(2),kpts%nmop(3),&
                            kpts%nkpt,sym%invs,noco%l_soc,sym%nop,&
                            sym%mrot,sym%tau)
         ELSE
         CALL julia(&
     &              sym,cell,input,noco,banddos,&
     &              kpts,.false.)
         ENDIF
       ELSE
        CALL od_kptsgen (kpts%nkpt)
       ENDIF
      ELSE
        IF(input%gw.eq.2) THEN
          INQUIRE(file='QGpsi',exist=l_kpts) ! Use QGpsi if it exists ot
          IF(l_kpts) THEN
            WRITE(6,*)&
     &        'QGpsi exists and will be used to generate kpts-file'
            OPEN (15,file='QGpsi',form='unformatted',status='old',&
     &        action='read')
            OPEN (41,file='kpts',form='formatted',status='unknown')
            REWIND(41)
            READ (15) kpts%nkpt
            WRITE (41,'(i5,f20.10)') kpts%nkpt,1.0
            DO n = 1, kpts%nkpt
              READ (15) q
              WRITE (41,'(4f10.5)') MATMUL(TRANSPOSE(cell%amat),q)/scale,1.0
              READ (15)
            ENDDO
            CLOSE (15)
            CLOSE (41)
          ENDIF
        ENDIF
      ENDIF
      
      dimension%neigd = max(dimension%neigd,input%gw_neigd)

!
! Using the k-point generator also for creation of q-points for the
! J-constants calculation:
      IF(.not.l_qpts)THEN
        kpts%nmop=nmopq
        l_tmp=(/noco%l_ss,noco%l_soc/)
        noco%l_ss=.false.
        noco%l_soc=.false.
        CALL julia(&
     &             sym,cell,input,noco,banddos,&
     &             kpts,.true.)
        noco%l_ss=l_tmp(1); noco%l_soc=l_tmp(2)
      ENDIF

!
! now proceed as usual
!
      CALL inpeig_dim(input,obsolete,cell,noco,oneD,jij,&
     &                kpts,dimension,stars)
      vacuum%layerd = max(vacuum%layerd,1)
      dimension%nstd = max(dimension%nstd,30)
      atoms%ntypd = atoms%ntype
      IF (noco%l_noco) dimension%neigd = 2*dimension%neigd

      atoms%nlod = max(atoms%nlod,2) ! for chkmt
      dimension%jspd=input%jspins
      CALL parawrite(&
     &               sym,stars,atoms,sphhar,dimension,vacuum,obsolete,&
     &               kpts,oneD)

!
      DEALLOCATE( &
     & atoms%lmax,atoms%ntypsy,atoms%neq,atoms%nlhtyp,atoms%rmt,atoms%zatom,atoms%jri,atoms%dx,atoms%nlo,atoms%llo,atoms%nflip,atoms%bmu,noel,&
     & vacuum%izlay,atoms%ncst,atoms%lnonsph,atoms%taual,atoms%pos,atoms%nz,atoms%relax,&
     & atoms%l_geo,noco%soc_opt,noco%alph,noco%beta,atoms%lda_u,noco%l_relax,jij%l_magn,jij%M,noco%b_con,sphhar%clnu,sphhar%nlh,&
     & sphhar%llh,sphhar%nmem,sphhar%mlh,jij%magtype,jij%nmagtype,hybrid%select1,hybrid%lcutm1,hybrid%select2,hybrid%lcutm2,&
     & hybrid%lcutwf)
!
      RETURN
      END SUBROUTINE dimen7
      END MODULE m_dimen7
