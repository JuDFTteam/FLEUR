!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_setinp
      use m_juDFT
!---------------------------------------------------------------------
!  Check muffin tin radii and determine a reasonable choice for MTRs.
!  Derive also other parameters for the input file, to provide some
!  help in the out-file.                                        gb`02
!---------------------------------------------------------------------
      CONTAINS
      SUBROUTINE set_inp(&
     &                   infh,nline,xl_buffer,bfh,buffer,l_hyb,&
     &                   atoms,sym,cell,title,idlist,&
     &                   input,vacuum,noco,&
     &                   atomTypeSpecies,speciesRepAtomType,numSpecies,&
     &                   a1,a2,a3)

      USE m_types
      USE iso_c_binding
      USE m_chkmt
      USE m_constants
      USE m_atominput
      USE m_lapwinput
      USE m_rwinp
      USE m_winpXML
      USE m_juDFT_init
      USE m_kpoints
      USE m_inv3
      USE m_types_xcpot_inbuild

      IMPLICIT NONE
      TYPE(t_input),INTENT(INOUT)    :: input
      TYPE(t_vacuum),INTENT(INOUT)   :: vacuum
      TYPE(t_noco),INTENT(INOUT)     :: noco
      TYPE(t_sym),INTENT(INOUT)      :: sym
      TYPE(t_cell),INTENT(INOUT)     :: cell
      TYPE(t_atoms),INTENT(INOUT)    :: atoms

      INTEGER, INTENT (IN) :: infh,xl_buffer,bfh,numSpecies
      INTEGER, INTENT (INOUT) :: nline
      INTEGER, INTENT (IN) :: atomTypeSpecies(atoms%ntype)
      INTEGER, INTENT (IN) :: speciesRepAtomType(atoms%nat)
      CHARACTER(len=xl_buffer) :: buffer
      LOGICAL, INTENT (IN) :: l_hyb  
      REAL,    INTENT (IN) :: idlist(:)
      REAL,    INTENT (INOUT) :: a1(3),a2(3),a3(3)
      CHARACTER(len=80), INTENT (IN) :: title
 
      INTEGER nel,i,j, nkptOld
      REAL    kmax,dtild,dvac1,n1,n2,gam,kmax0,dtild0,dvac0,sumWeight
      REAL    recVecLength, kPointDen(3)
      LOGICAL l_test,l_gga,l_exists, l_explicit, l_kpts
      REAL     dx0(atoms%ntype), rmtTemp(atoms%ntype)
      REAL     a1Temp(3),a2Temp(3),a3Temp(3) 
      INTEGER  div(3)
      INTEGER jri0(atoms%ntype),lmax0(atoms%ntype),nlo0(atoms%ntype),llo0(atoms%nlod,atoms%ntype)
      CHARACTER(len=1)  :: ch_rw
      CHARACTER(len=4)  :: namex
      CHARACTER(len=3)  :: noel(atoms%ntype)
      CHARACTER(len=12) :: relcor
      CHARACTER(len=3)  :: latnamTemp
      CHARACTER(LEN=20) :: filename
      INTEGER  nu,iofile
      INTEGER  iggachk
      INTEGER  n ,iostat, errorStatus
      REAL     scpos ,zc

      TYPE(t_banddos)::banddos
      TYPE(t_obsolete)::obsolete
      TYPE(t_sliceplot)::sliceplot
      TYPE(t_oneD)::oneD
      TYPE(t_stars)::stars
      TYPE(t_hybrid)::hybrid
      TYPE(t_xcpot_inbuild)::xcpot
      TYPE(t_kpts)::kpts
      TYPE(t_enpara)::enpara
      TYPE(t_forcetheo)::forcetheo

    !-odim
!+odim
!      REAL, PARAMETER :: eps=0.00000001
!     ..
!HF   added for HF and hybrid functionals
      REAL     ::  taual_hyb(3,atoms%nat)
      INTEGER  ::  bands
      LOGICAL  ::  l_gamma
      INTEGER  :: nkpt3(3)
!HF

      INTEGER :: xmlElectronStates(29,atoms%ntype)
      LOGICAL :: xmlPrintCoreStates(29,atoms%ntype)
      REAL    :: xmlCoreOccs(2,29,atoms%ntype)
      REAL    :: xmlCoreRefOccs(29)

      interface
         function dropInputSchema() bind(C, name="dropInputSchema")
            use iso_c_binding
            INTEGER(c_int) dropInputSchema
         end function dropInputSchema
      end interface

      DATA xmlCoreRefOccs /2,2,2,4,2,2,4,2,4,6,2,4,2,4,6,2,4,2,6,8,4,&
     &                     6,2,4,2,6,8,4,6/
      xmlElectronStates = noState_const
      xmlPrintCoreStates = .FALSE.
      xmlCoreOccs = 0.0

      l_test = .false.
      l_gga  = .true.
      atoms%nlod=9
      ALLOCATE(atoms%nz(atoms%ntype))
      ALLOCATE(atoms%jri(atoms%ntype))
      ALLOCATE(atoms%dx(atoms%ntype))
      ALLOCATE(atoms%lmax(atoms%ntype))
      ALLOCATE(atoms%nlo(atoms%ntype))
      ALLOCATE(atoms%llo(atoms%nlod,atoms%ntype))
      ALLOCATE(atoms%ncst(atoms%ntype))
      ALLOCATE(atoms%lnonsph(atoms%ntype))
      ALLOCATE(atoms%nflip(atoms%ntype))
      ALLOCATE(atoms%l_geo(atoms%ntype))
      ALLOCATE(atoms%lda_u(atoms%ntype))
      ALLOCATE(atoms%bmu(atoms%ntype))
      ALLOCATE(atoms%relax(3,atoms%ntype))
      ALLOCATE(atoms%ulo_der(atoms%nlod,atoms%ntype))
      
      atoms%nz(:) = NINT(atoms%zatom(:))
      DO i = 1, atoms%ntype
       noel(i) = namat_const(atoms%nz(i))
      ENDDO
      atoms%rmt(:) = 999.9
      atoms%pos(:,:) = matmul( cell%amat , atoms%taual(:,:) )
      atoms%ulo_der = 0
      ch_rw = 'w'
      sym%namgrp= 'any ' 
      banddos%dos   = .false. ; banddos%l_mcd = .false. ; banddos%unfoldband = .FALSE. ; input%secvar = .false.
      input%vchk = .false. ; input%cdinf = .false. 
      input%l_bmt= .false. ; input%eonly  = .false.
      input%gauss= .false. ; input%tria  = .false. 
      sliceplot%slice= .false. ;  input%swsp  = .false.
      input%lflip= .false. ; banddos%vacdos= .false. ; input%integ = .false.
      sliceplot%iplot= .false. ; input%score = .false. ; sliceplot%plpot = .false.
      input%pallst = .false. ; obsolete%lwb = .false. ; vacuum%starcoeff = .false.
      input%strho  = .false.  ; input%l_f = .false. ; atoms%l_geo(:) = .true.
      noco%l_noco = noco%l_ss ;   input%jspins = 1
      input%itmax = 9 ; input%maxiter = 99 ; input%imix = 7 ; input%alpha = 0.05
      input%preconditioning_param = 0.0 ; input%minDistance = 1.0e-5
      input%spinf = 2.0 ; obsolete%lepr = 0 ; input%coretail_lmax = 0
      sliceplot%kk = 0 ; sliceplot%nnne = 0  ; vacuum%nstars = 0 ; vacuum%nstm = 0 
      nu = 5 ; vacuum%layerd = 1 ; iofile = 6
      ALLOCATE(vacuum%izlay(vacuum%layerd,2))
      banddos%ndir = 0 ; vacuum%layers = 0 ; atoms%nflip(:) = 1 ; vacuum%izlay(:,:) = 0
      banddos%e_mcd_lo = -10.0 ; banddos%e_mcd_up = 0.0
      atoms%lda_u%l = -1 ; atoms%relax(1:2,:) = 1 ; atoms%relax(:,:) = 1
      input%epsdisp = 0.00001 ; input%epsforce = 0.00001 ; input%xa = 2.0 ; input%thetad = 330.0
      sliceplot%e1s = 0.0 ; sliceplot%e2s = 0.0 ; banddos%e1_dos = 0.5 ; banddos%e2_dos = -0.5 ; input%tkb = 0.001
      banddos%sig_dos = 0.015 ; vacuum%tworkf = 0.0 ; input%scaleCell = 1.0 ; scpos = 1.0
      input%scaleA1 = 1.0 ; input%scaleA2 = 1.0 ; input%scaleC = 1.0
      zc = 0.0 ; vacuum%locx(:) = 0.0 ;  vacuum%locy(:) = 0.0
      kpts%numSpecialPoints = 0
      input%ldauLinMix = .FALSE. ; input%ldauMixParam = 0.05 ; input%ldauSpinf = 1.0
      input%l_wann = .FALSE.

!+odim
      oneD%odd%mb = 0 ; oneD%odd%M = 0 ; oneD%odd%m_cyl = 0 ; oneD%odd%chi = 0 ; oneD%odd%rot = 0
      oneD%odd%k3 = 0 ; oneD%odd%n2d= 0 ; oneD%odd%nq2 = 0 ; oneD%odd%nn2d = 0 
      oneD%odd%nop = 0 ; oneD%odd%kimax2 = 0 ; oneD%odd%nat = 0
      oneD%odd%invs = .false. ; oneD%odd%zrfs = .false. ; oneD%odd%d1 = .false.
!-odim
! check for magnetism
      atoms%bmu(:) = 0.0
      DO n = 1, atoms%ntype
        IF (atoms%nz(n).EQ.24) atoms%bmu(n) = 1.0  ! Cr - Ni
        IF (atoms%nz(n).EQ.25) atoms%bmu(n) = 3.5
        IF (atoms%nz(n).EQ.26) atoms%bmu(n) = 2.2
        IF (atoms%nz(n).EQ.27) atoms%bmu(n) = 1.6
        IF (atoms%nz(n).EQ.28) atoms%bmu(n) = 1.1
        IF (atoms%nz(n).EQ.59) atoms%bmu(n) = 2.1  ! Pr - Tm
        IF (atoms%nz(n).EQ.60) atoms%bmu(n) = 3.1
        IF (atoms%nz(n).EQ.61) atoms%bmu(n) = 4.1
        IF (atoms%nz(n).EQ.62) atoms%bmu(n) = 5.1
        IF (atoms%nz(n).EQ.63) atoms%bmu(n) = 7.1
        IF (atoms%nz(n).EQ.64) atoms%bmu(n) = 7.1 
        IF (atoms%nz(n).EQ.65) atoms%bmu(n) = 6.1
        IF (atoms%nz(n).EQ.66) atoms%bmu(n) = 5.1
        IF (atoms%nz(n).EQ.67) atoms%bmu(n) = 4.1
        IF (atoms%nz(n).EQ.68) atoms%bmu(n) = 3.1
        IF (atoms%nz(n).EQ.69) atoms%bmu(n) = 2.1
      ENDDO

      input%delgau = input%tkb ; atoms%ntype = atoms%ntype ; atoms%nat = atoms%nat
      DO i = 1, 10
        j = (i-1) * 8 + 1
        input%comment(i) = title(j:j+7)
      ENDDO 
      IF (noco%l_noco) input%jspins = 2
       
      a1(:) = cell%amat(:,1) ; a2(:) = cell%amat(:,2) ; a3(:) = cell%amat(:,3) 

      CALL chkmt(&
     &           atoms,input,vacuum,cell,oneD,&
     &           l_gga,noel,l_test,&
     &           kmax,dtild,vacuum%dvac,atoms%lmax,atoms%jri,atoms%rmt,atoms%dx)

! --> read in (possibly) atomic info

      stars%gmax = 3.0 * kmax ; xcpot%gmaxxc = 2.5 * kmax ; input%rkmax = kmax
      atoms%lnonsph(:) = min( max( (atoms%lmax(:)-2),3 ), 8 )

      CALL atom_input(&
     &                infh,xl_buffer,bfh,buffer,&
     &                input,idlist,xmlCoreRefOccs,&
     &                nline,&
     &                xmlElectronStates,xmlPrintCoreStates,xmlCoreOccs,&
     &                atomTypeSpecies,numSpecies,&
     &                nel,atoms,enpara)

      DO n = 1, atoms%ntype
         IF (atoms%lnonsph(n).GT.atoms%lmax(n)) THEN
            WRITE(*,'(a20,i5,a25,i3,a4,i3,a1)')& 
               'NOTE: For atom type ', n,' lnonsph is reduced from ',& 
               atoms%lnonsph(n),' to ', atoms%lmax(n), '.'
            WRITE(6,'(a20,i5,a25,i3,a4,i3,a1)')&
               'NOTE: For atom type ', n, ' lnonsph is reduced from ',& 
               atoms%lnonsph(n),' to ', atoms%lmax(n), '.'
            atoms%lnonsph(n) = atoms%lmax(n)
         END IF
      END DO

      input%zelec = nel

! --> check once more
      rmtTemp = 999.0
      l_test = .true.
      CALL chkmt(&
     &           atoms,input,vacuum,cell,oneD,&
     &           l_gga,noel,l_test,&
     &           kmax0,dtild0,dvac0,lmax0,jri0,rmtTemp,dx0)

      IF ( ANY(atoms%nlo(:).NE.0) ) THEN
        input%ellow = -1.8
      ELSE
        input%ellow = -0.8  
      ENDIF
      IF (input%film) THEN
         input%elup = 0.5
      ELSE
         input%elup = 1.0
      ENDIF 

      IF (.not.input%film) THEN
         vacuum%dvac = a3(3) ; dtild = vacuum%dvac
      ENDIF
      IF ( (abs(a1(3)).GT.eps).OR.(abs(a2(3)).GT.eps).OR.&
     &     (abs(a3(1)).GT.eps).OR.(abs(a3(2)).GT.eps) ) THEN          
        cell%latnam = 'any'
      ELSE
        IF ( (abs(a1(2)).LT.eps).AND.(abs(a2(1)).LT.eps) ) THEN
          IF (abs(a1(1)-a2(2)).LT.eps) THEN
            cell%latnam = 'squ'
          ELSE
            cell%latnam = 'p-r'
          ENDIF
        ELSE
          n1 = sqrt(a1(1)**2 + a1(2)**2); n2 = sqrt(a2(1)**2 + a2(2)**2)
          IF (abs(n1-n2).LT.eps) THEN
            gam = ( a1(1)*a2(1) + a1(2)*a2(2) ) / (n1 * n2)
            gam = 57.295779512*acos(gam)
            IF (abs(gam-60.).LT.eps) THEN
               cell%latnam = 'hex'
               a1(2) = n1 * 0.5
               a1(1) = a1(2) * sqrt(3.0)
            ELSEIF (abs(gam-120.).LT.eps) THEN
               cell%latnam = 'hx3'
               a1(1) = n1 * 0.5
               a1(2) = a1(1) * sqrt(3.0)
            ELSE
               cell%latnam = 'c-r'
               gam = 0.5 * gam / 57.295779512
               a1(1) =  n1 * cos(gam)
               a1(2) = -n1 * sin(gam)
            ENDIF
            a2(1) =   a1(1)
            a2(2) = - a1(2)
          ELSE
            cell%latnam = 'obl'
          ENDIF
        ENDIF
      ENDIF

!HF   added for HF and hybrid functionals
      hybrid%gcutm1       = input%rkmax - 0.5
      hybrid%tolerance1   = 1e-4
      taual_hyb   = atoms%taual
      ALLOCATE(hybrid%lcutwf(atoms%ntype))
      ALLOCATE(hybrid%lcutm1(atoms%ntype))
      ALLOCATE(hybrid%select1(4,atoms%ntype))
      hybrid%lcutwf      = atoms%lmax - atoms%lmax / 10
      hybrid%ewaldlambda = 3
      hybrid%lexp        = 16
      hybrid%lcutm1 = 4
      hybrid%select1(1,:) = 4
      hybrid%select1(2,:) = 0
      hybrid%select1(3,:) = 4
      hybrid%select1(4,:) = 2
      bands       = max( nint(input%zelec)*10, 60 )
      l_gamma     = .false.
      hybrid%l_hybrid = l_hyb
      IF (l_hyb) THEN
         input%ellow = input%ellow -  2.0
         input%elup  = input%elup  + 10.0
         input%gw_neigd = bands
         l_gamma = .true.
         input%minDistance = 1.0e-5
      ELSE
        input%gw_neigd = 0
      END IF
!HF

! rounding
      atoms%rmt(:)  = real(NINT(atoms%rmt(:)  * 100 ) / 100.)
      atoms%dx(:)   = real(NINT(atoms%dx(:)   * 1000) / 1000.)
      stars%gmax    = real(NINT(stars%gmax    * 10  ) / 10.)
      input%rkmax   = real(NINT(input%rkmax   * 10  ) / 10.)
      xcpot%gmaxxc  = real(NINT(xcpot%gmaxxc  * 10  ) / 10.)
      hybrid%gcutm1 = real(NINT(hybrid%gcutm1 * 10  ) / 10.)
      IF (input%film) THEN
       vacuum%dvac = real(NINT(vacuum%dvac*100)/100.)
       dtild = real(NINT(dtild*100)/100.)
      ENDIF
!
! read some lapw input
!
      CALL lapw_input(&
     &                infh,nline,xl_buffer,bfh,buffer,&
     &                input%jspins,input%kcrel,obsolete%ndvgrd,kpts%nkpt,div,kpts%kPointDensity,&
     &                input%frcor,input%ctail,obsolete%chng,input%tria,input%rkmax,stars%gmax,xcpot%gmaxxc,&
     &                vacuum%dvac,dtild,input%tkb,namex,relcor)

      stars%gmaxInit = stars%gmax
!
      IF (input%film) atoms%taual(3,:) = atoms%taual(3,:) * a3(3) / dtild

      INQUIRE(file="inp",exist=l_exists)
      IF (l_exists) THEN
         CALL juDFT_error("inp-file exists. Cannot write another input file in this directory.",calledby="set_inp")
      ENDIF
      INQUIRE(file="inp.xml",exist=l_exists)
      IF (l_exists) THEN
         CALL juDFT_error("inp.xml-file exists. Cannot write another input file in this directory.",calledby="set_inp")
      ENDIF

      nu = 8 
      input%gw = 0

      IF (kpts%nkpt == 0) THEN     ! set some defaults for the k-points
        IF (input%film) THEN
          cell%area = cell%omtil / vacuum%dvac
          kpts%nkpt = MAX(nint((3600/cell%area)/sym%nop2),1)
        ELSE
          kpts%nkpt = MAX(nint((216000/cell%omtil)/sym%nop),1)
        ENDIF
      ENDIF

      kpts%specificationType = 0
      IF((ANY(div(:).NE.0)).AND.(ANY(kpts%kPointDensity(:).NE.0.0))) THEN
         CALL juDFT_error('Double specification of k point set', calledby = 'set_inp')
      END IF
      IF (ANY(div(:).NE.0)) THEN
         kpts%specificationType = 2
      ELSE IF (ANY(kpts%kPointDensity(:).NE.0.0)) THEN
         kpts%specificationType = 4
      ELSE
         kpts%specificationType = 1
      END IF
      l_kpts = .FALSE.

      IF(TRIM(ADJUSTL(sym%namgrp)).EQ.'any') THEN
         sym%symSpecType = 1
      ELSE
         sym%symSpecType = 2
      END IF

      ! set vacuum%nvac
      vacuum%nvac = 2
      IF (sym%zrfs.OR.sym%invs) vacuum%nvac = 1
      IF (oneD%odd%d1) vacuum%nvac = 1
      
      ! Set defaults for noco  types
      ALLOCATE(noco%l_relax(atoms%ntype),noco%b_con(2,atoms%ntype))
      ALLOCATE(noco%alphInit(atoms%ntype),noco%alph(atoms%ntype),noco%beta(atoms%ntype))
   
      IF (noco%l_ss) input%ctail = .FALSE.
      noco%l_mperp = .FALSE.
      noco%l_constr = .FALSE.
      noco%mix_b = 0.0
      noco%qss = 0.0

      noco%l_relax(:) = .FALSE.
      noco%alphInit(:) = 0.0
      noco%alph(:) = 0.0
      noco%beta(:) = 0.0
      noco%b_con(:,:) = 0.0

     
      CALL inv3(cell%amat,cell%bmat,cell%omtil)
      cell%bmat=tpi_const*cell%bmat
      kpts%nkpt3(:) = div(:)

      IF (kpts%specificationType.EQ.4) THEN
         DO i = 1, 3
            IF (kpts%kPointDensity(i).LE.0.0) THEN
               CALL juDFT_error('Error: Nonpositive kpointDensity provided', calledby = 'set_inp')
            END IF
            recVecLength = SQRT(cell%bmat(i,1)**2 + cell%bmat(i,2)**2 + cell%bmat(i,3)**2)
            kpts%nkpt3(i) = CEILING(kpts%kPointDensity(i) * recVecLength)
         END DO
         kpts%nkpt = kpts%nkpt3(1) * kpts%nkpt3(2) * kpts%nkpt3(3)
      END IF

      IF (l_hyb) THEN
         ! Changes for hybrid functionals
         namex = 'pbe0'
         input%ctail = .false. ; atoms%l_geo = .false.! ; input%frcor = .true.
         input%itmax = 15 ; input%maxiter = 25!; input%imix  = 17
         IF (ANY(kpts%nkpt3(:).EQ.0)) kpts%nkpt3(:) = 4
         div(:) = kpts%nkpt3(:)
         kpts%specificationType = 2
      END IF

         nkptOld = kpts%nkpt
         latnamTemp = cell%latnam

         l_explicit = juDFT_was_argument("-explicit")

         a1Temp(:) = a1(:)
         a2Temp(:) = a2(:)
         a3Temp(:) = a3(:)

         IF(l_explicit) THEN
            ! kpts generation
            kpts%l_gamma = l_gamma

            CALL kpoints(oneD,sym,cell,input,noco,banddos,kpts,l_kpts)

            kpts%specificationType = 3
            IF (l_hyb) kpts%specificationType = 2
         END IF

         IF(l_explicit) THEN
            sym%symSpecType = 3
            !set latnam to any
            cell%latnam = 'any'

            a1Temp(:) = cell%amat(:,1)
            a2Temp(:) = cell%amat(:,2)
            a3Temp(:) = cell%amat(:,3)
         END IF

         errorStatus = 0
         errorStatus = dropInputSchema()
         IF(errorStatus.NE.0) THEN
            STOP 'Error: Cannot print out FleurInputSchema.xsd'
         END IF
         filename = 'inp.xml'

         CALL w_inpXML(&
     &                 atoms,obsolete,vacuum,input,stars,sliceplot,forcetheo,banddos,&
     &                 cell,sym,xcpot,noco,oneD,hybrid,kpts,div,l_gamma,&
     &                 noel,namex,relcor,a1Temp,a2Temp,a3Temp,dtild,input%comment,&
     &                 xmlElectronStates,xmlPrintCoreStates,xmlCoreOccs,&
     &                 atomTypeSpecies,speciesRepAtomType,.FALSE.,filename,&
     &                 l_explicit,numSpecies,enpara)

         IF(juDFT_was_argument("-explicit")) THEN
            sumWeight = 0.0
            WRITE(6,*) ''
            WRITE(6,'(a,(a3,i10))') 'k-point count:','', kpts%nkpt
            DO i = 1, kpts%nkpt
               sumWeight = sumWeight + kpts%wtkpt(i)
            END DO
            DO i = 1, kpts%nkpt
               kpts%wtkpt(i) = kpts%wtkpt(i) / sumWeight
               kpts%wtkpt(i) = kpts%wtkpt(i)
            END DO
         END IF

         kpts%nkpt = nkptOld
         cell%latnam = latnamTemp
 
      DEALLOCATE (noco%l_relax,noco%b_con,noco%alphInit,noco%alph,noco%beta)
      DEALLOCATE (atoms%ulo_der)

      IF (ANY(kpts%nkpt3(:).NE.0)) THEN
         DO i = 1, 3
            recVecLength = SQRT(cell%bmat(i,1)**2 + cell%bmat(i,2)**2 + cell%bmat(i,3)**2)
            kPointDen(i) = kpts%nkpt3(i) / recVecLength
         END DO
         WRITE(6,*) ''
         WRITE(6,'(a,3(a4,i10))')   'k-point mesh:'   , '', kpts%nkpt3(1),'', kpts%nkpt3(2),'', kpts%nkpt3(3)
         WRITE(6,'(a,3(a4,f10.6))') 'k-point density:', '', kPointDen(1),'', kPointDen(2),'', kPointDen(3)
         WRITE(6,*) ''
      END IF

      CLOSE (6)


      END SUBROUTINE set_inp
      END MODULE m_setinp
