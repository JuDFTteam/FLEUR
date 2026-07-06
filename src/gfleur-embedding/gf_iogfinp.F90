!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_iogfinp
      IMPLICIT NONE
      CONTAINS
      !<--S: GF_RGFINP

      SUBROUTINE GF_RGFINP(embinp,gfmix,layer_info)
!***********************************************************************
!     This subroutine:
!     - reads the gf_inp file
!     - sets some variables to defaults
!
!     In the port to current FLEUR the gf_inp file only carries the
!     GF-specific control parameters. Everything that describes the
!     electronic structure of a layer (basis cutoffs, LOs, lmax,
!     LDA+U, noco, SOC, xc-potential) now comes from the per-layer
!     inp.xml files; the corresponding gf_inp sections were removed
!     and produce an error with a hint. The layer decomposition itself
!     is given in the new &LAYERS namelist.
!
!     Daniel Wortmann, Juelich 2003
!***********************************************************************
      USE m_gf_types
      use m_juDFT
      IMPLICIT NONE
      !<-- Arguments
      TYPE(t_embinp), INTENT(INOUT) :: embinp
      TYPE(t_gfmix), INTENT(INOUT)  :: gfmix
      TYPE(t_layers), INTENT(INOUT) :: layer_info
      !>
      !<-- locals!

      !layer decomposition
      INTEGER :: nlayers
      CHARACTER(LEN=20) :: prefix
      REAL :: c1(1000),c2(1000),d(1000),dt(1000),c(1000)
      NAMELIST /LAYERS/ nlayers,prefix,c1,c2,d,dt,c

      !mode switches
      LOGICAL :: eigen,tmat,gmat,charge,dos,totalmix,writeT,potmix      &
     &     ,inp2plot,savemem,writehs,IEC,gproj,band,spectral,fullgreen  &
     &     ,nohelpregion,overlap,intdos,hdfio
      INTEGER :: curr,nz_dos
      REAL    :: z_min,z_max
      CHARACTER(LEN=4)::kpts
      LOGICAL :: doslayer(1000)
      NAMELIST /MODE/ intdos,totalmix,eigen,tmat,gmat,curr,charge,dos   &
     &     ,writeT,potmix,inp2plot,kpts,savemem,writehs,IEC,gproj,band  &
     &     ,spectral,fullgreen,nz_dos,z_min,z_max,nohelpregion,overlap,hdfio,doslayer
      !CBS
      REAL   ::eps_current,eps_non_bloch,CBS_bz,d1,d2,kappa_max
      LOGICAL:: cbs
      NAMELIST /CBSMODE/ CBS,eps_current,eps_non_bloch,CBS_bz,d1,d2     &
     &     ,kappa_max
      !BASIS (only the GF-specific region-II description remains;
      !rkmax is taken from the per-layer inp.xml)
      LOGICAL::solwil
      INTEGER::npw
      INTEGER :: napw(1000)
      NAMELIST /BASIS/solwil,npw,napw
      !EMB
      LOGICAL::addemb,advanced,surface,embmt,embspectral,simplevacuum
      REAL  ::charge_limit,imag_broad,EField,vacuum_energy
      NAMELIST /EMB/addemb,advanced,surface,charge_limit,embmt,         &
     &       imag_broad,Efield,embspectral,simplevacuum,vacuum_energy
      !SC
      REAL::trans,c_sc
      INTEGER::imix,maxiter,iter,precond
      REAL::alpha,spinf
      LOGICAL::gauss,tria
      REAL::tkb,delgau
      REAL::gmax_pot,gmax_decouple,g0max,k_kerker,g0scale
      NAMELIST /SC/iter,trans,c_sc,imix,maxiter,alpha,spinf,gauss,tria  &
     &     ,tkb,delgau,gmax_pot,gmax_decouple,k_kerker,g0max,precond,g0scale

      !for free input
      CHARACTER(LEN = 10) :: teststring
      LOGICAL             :: free_input
      LOGICAL             :: required(2)

      !>
      !<--set defaults
      CALL timestart("Reading the gf_inp")
      !LAYERS
      nlayers = 0
      prefix = "layer"
      c1 = 0.0;c2 = 0.0;d = 0.0;dt = 0.0;c = 0.0
      !CBS
      embinp%l_CBS = .FALSE.;embinp%eps_current = 1E-4
      embinp%eps_non_bloch = 1E-4;embinp%CBS_bz =-1
      embinp%dp1     = 0.0;embinp%dp2 = 0.0
      embinp%kappa_max=1E20
      !EMB
      embinp%l_embmt = .FALSE.
      embinp%charge_limit = 0.0
      embinp%l_addemb = .FALSE.
      embinp%l_adv = .FALSE.
      embinp%l_surface  = .FALSE.
      embinp%l_embspectral = .FALSE.
      embinp%imag_broad = 0.0
      embinp%l_simplevacuum=.FALSE.
      embinp%Efield   = 0.0
      !SC
      embinp%trans = 0.0;embinp%c_sc = 0.0
      gfmix%imix = 0;gfmix%maxiter = 99;gfmix%iter =0;
      gfmix%alpha = 0.01;gfmix%spinf = 1.0
      gfmix%precond=0;gfmix%k_kerker=0.0;gfmix%g0max=0.0;gfmix%g0scale=1.0
      embinp%l_gauss = .FALSE.;embinp%l_tria = .FALSE.
      embinp%tkb = 0.005;embinp%delgau = 0.0
      embinp%gmax_pot=50.0
      embinp%gmax_decouple=10.0
      !BASIS
      embinp%l_solwil = .FALSE.;embinp%npw = 0
      napw = 15
      required = .FALSE.

      !>

      !<--READ the gf_inp-file

      OPEN(92,FILE='gf_inp',FORM='formatted',STATUS='OLD')
      free_input = .TRUE.
      DO WHILE(free_input)
         !<--read a teststring
         READ(92,'(a10)',END=999) teststring
                                          !remove heading blanks
         teststring = ADJUSTL(teststring)
         !>
                                            !skip comment lines
         IF (teststring(1:1) =='!') CYCLE
         IF (LEN_TRIM(teststring)==0) CYCLE
                       !The current line has to be read again
         BACKSPACE(92)
         !<--Now test for the different possible input-lines

         IF (teststring(1:4)=='&LAY') THEN
            READ(92,NML = LAYERS)
            IF (nlayers<1) CALL juDFT_error("nlayers<1 in &LAYERS",     &
     &           calledby="gf_iogfinp")
            layer_info%num_layers = nlayers
            layer_info%prefix = prefix
            IF (ALLOCATED(layer_info%c1)) DEALLOCATE(layer_info%c1,layer_info%c2,   &
     &           layer_info%d,layer_info%dt,layer_info%c)
            ALLOCATE(layer_info%c1(nlayers),layer_info%c2(nlayers),             &
     &           layer_info%d(nlayers),layer_info%dt(nlayers),layer_info%c(nlayers))
            layer_info%c1 = c1(:nlayers);layer_info%c2 = c2(:nlayers)
            layer_info%d = d(:nlayers);layer_info%dt = dt(:nlayers)
            layer_info%c = c(:nlayers)
            required(1) = .TRUE.
         ELSEIF (teststring(1:4)=='gmax') THEN
            CALL juDFT_error("The 'gmax' line was removed from gf_inp"  &
     &           ,hint="gmax and gmaxxc are taken from the per-layer "//&
     &           "inp.xml files now",calledby="gf_iogfinp")
         ELSEIF (teststring(1:4)=='&MOD') THEN
            !READ the mode-switches
            spectral = .TRUE.;fullgreen=.FALSE.
            inp2plot = .FALSE.;band = .FALSE.
            eigen = .FALSE.;tmat=.TRUE.;curr=0;
            charge = .FALSE.;dos=.FALSE.;IEC=.FALSE.;intdos=.FALSE.
            writeT = .FALSE.;gmat=.TRUE.;potmix=.FALSE.;kpts="none"
            writeHS = .FALSE.;saveMEM=.FALSE.;gproj=.FALSE.             &
     &           ;nohelpregion = .TRUE.;overlap=.FALSE.
            nz_dos = 1;z_min=0.0;z_max=1.0;totalmix= .FALSE.;hdfio=.false.
            doslayer=.true.
            READ(92,NML = MODE)
            embinp%l_hdfio=hdfio
            embinp%l_intdos=intdos
            embinp%l_totalmix=totalmix
            embinp%l_spectral=spectral;embinp%l_fullgreen=fullgreen
            embinp%l_eigen = eigen;embinp%l_tmat=tmat;embinp%curr=curr
            embinp%l_charge = charge;embinp%l_dos=dos;embinp%l_band=band
            embinp%l_writeT =writeT;embinp%l_gmat=gmat;
            embinp%l_potmix = potmix;embinp%l_inp2plot=inp2plot
            embinp%l_addselfen=.FALSE.
            required(2)    = .TRUE.;embinp%kpts =                       &
     &           kpts;embinp%l_nohelpregion = nohelpregion
            embinp%l_writeHS = writeHS;embinp%l_savemem=savemem
            embinp%l_IEC=IEC;embinp%l_gproj=gproj
            embinp%nz_dos=nz_dos;embinp%z_min=z_min;embinp%z_max=z_max
            IF (.NOT.ALLOCATED(embinp%l_doslayer))                      &
     &           ALLOCATE(embinp%l_doslayer(1000))
            embinp%l_doslayer=doslayer
         ELSEIF (teststring(1:4) =='&CBS') THEN
            !READ the CBS
            kappa_max = 1E20;
            CBS = .FALSE.;eps_current = 1E-4;eps_non_bloch=1E-4
            CBS_bz =-1.0;d1 = 0.0;d2 = 0.0
            READ(92,NML = CBSMODE)
            embinp%kappa_max=kappa_max
            embinp%l_CBS = CBS;embinp%eps_current=eps_current
            embinp%eps_non_bloch = eps_non_bloch;embinp%CBS_bz=CBS_bz
            embinp%dp1 = d1;embinp%dp2=d2
         ELSEIF (teststring(1:4) =='&BAS') THEN
            !READ the BASIS
            solwil=.FALSE.;npw=0;napw=15
            READ(92,NML = BASIS)
            embinp%l_solwil = solwil;embinp%npw = npw
            IF (.NOT.ALLOCATED(embinp%napw))                            &
     &           ALLOCATE(embinp%napw(1000))
            embinp%napw = napw
         ELSEIF (teststring(1:4) =='&EMB') THEN
            !READ EMB
            embmt=.FALSE.
            addemb = .FALSE.
            simplevacuum=.false.
            embspectral=.FALSE.
            advanced = .FALSE.
            surface  = .FALSE.
            imag_broad = 0.0
            charge_limit=0.0
            Efield = 0.0
            vacuum_energy=1.0
            READ(92,NML = EMB)
            embinp%vacuum_energy=vacuum_energy
            embinp%l_simplevacuum=simplevacuum
            embinp%l_embspectral=embspectral
            embinp%charge_limit = charge_limit
            embinp%l_addemb = addemb
            embinp%l_adv    = advanced
            embinp%l_embmt  = embmt
            embinp%l_surface = surface
            embinp%imag_broad = imag_broad
            embinp%Efield   = Efield
         ELSEIF (teststring(1:3) =='&SC') THEN
            !READ SC
            gmax_pot = 50.0
            gmax_decouple = 10.0
            trans = 0.0;c_sc=0.0;imix=1;maxiter=99;alpha=0.05;spinf=0.2
            gauss = .FALSE.;tria=.FALSE.;precond=0;k_kerker=0
            tkb = 0.005;delgau=0.0;g0scale=1.0;g0max=0.0
            READ(92,NML = SC)
            embinp%trans = trans;embinp%c_sc=c_sc
            gfmix%imix    = imix;gfmix%maxiter = maxiter;gfmix%iter=iter
            gfmix%alpha   = alpha;gfmix%spinf = spinf
            gfmix%precond=precond;gfmix%k_kerker=k_kerker
            gfmix%g0max=g0max;gfmix%g0scale=g0scale
            embinp%l_gauss = gauss;embinp%l_tria=tria
            embinp%tkb   = tkb;embinp%delgau = delgau
            embinp%gmax_pot=gmax_pot
            embinp%gmax_decouple=gmax_decouple
         ELSEIF (teststring(1:4) =='&SPI') THEN
            CALL juDFT_error("&SPINORBIT was removed from gf_inp",      &
     &           hint="SOC is specified in the per-layer inp.xml now",  &
     &           calledby="gf_iogfinp")
         ELSEIF (teststring(1:4) =='&NON'.OR.teststring(1:4)=='&NOC')THEN
            CALL juDFT_error("&NONCOLLINEAR/&NOCOL were removed from "//&
     &           "gf_inp",hint="noco is specified in the per-layer "//  &
     &           "inp.xml now",calledby="gf_iogfinp")
         ELSEIF (teststring(1:7) =='*****XC') THEN
            CALL juDFT_error("The XC section was removed from gf_inp",  &
     &           hint="The xc-potential and the atom setup (lmax, LOs, "&
     &           //"LDA+U) are taken from the per-layer inp.xml files " &
     &           //"now - delete everything from *****XC on",           &
     &           calledby="gf_iogfinp")
         ELSE
            WRITE(*,*) "Unknown entry in gf_inp-file"
            WRITE(*,*) "Line starting with:",teststring
            CALL juDFT_error("ERROR in gf_inp")
         ENDIF

         !>
      ENDDO
  999 CONTINUE
      CLOSE(92)

      IF (.NOT.required(1)) CALL juDFT_error                            &
     &     ("&LAYERS section missing in gf_inp",calledby="gf_iogfinp")
      IF (.NOT.required(2)) CALL juDFT_error                            &
     &     ("&MODE section missing in gf_inp",calledby="gf_iogfinp")

      !resize the per-layer arrays now that num_layers is known
      CALL priv_resize_log(embinp%l_doslayer,layer_info%num_layers,.TRUE.)
      CALL priv_resize_int(embinp%napw,layer_info%num_layers,15)

      IF ( embinp%curr/=0 .AND. embinp%l_cbs )                          &
     &     CALL juDFT_error('embtmat: l_CBS&curr')

      CALL timestop("Reading the gf_inp")
      !>
      END SUBROUTINE

      SUBROUTINE priv_resize_log(arr,n,default)
      LOGICAL,ALLOCATABLE,INTENT(INOUT) :: arr(:)
      INTEGER,INTENT(IN) :: n
      LOGICAL,INTENT(IN) :: default
      LOGICAL,ALLOCATABLE :: tmp(:)
      IF (.NOT.ALLOCATED(arr)) THEN
         ALLOCATE(arr(n))
         arr = default
      ELSE
         ALLOCATE(tmp(n))
         tmp = default
         tmp(:MIN(n,SIZE(arr))) = arr(:MIN(n,SIZE(arr)))
         CALL MOVE_ALLOC(tmp,arr)
      ENDIF
      END SUBROUTINE

      SUBROUTINE priv_resize_int(arr,n,default)
      INTEGER,ALLOCATABLE,INTENT(INOUT) :: arr(:)
      INTEGER,INTENT(IN) :: n
      INTEGER,INTENT(IN) :: default
      INTEGER,ALLOCATABLE :: tmp(:)
      IF (.NOT.ALLOCATED(arr)) THEN
         ALLOCATE(arr(n))
         arr = default
      ELSE
         ALLOCATE(tmp(n))
         tmp = default
         tmp(:MIN(n,SIZE(arr))) = arr(:MIN(n,SIZE(arr)))
         CALL MOVE_ALLOC(tmp,arr)
      ENDIF
      END SUBROUTINE

      END
