!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_io2dmat 
#ifdef CPP_MPI
      USE mpi
#endif
      USE,INTRINSIC::iso_c_binding
!-----------------------------------------------                        
! New version of the subroutines needed to read/write 2D matrices       
!                 Daniel Wortmann, (05-11-18)                           
!-----------------------------------------------                        
      USE hdf5 
      USE m_hdf_tools
      USE m_juDFT
      USE m_constants, ONLY: oUnit
      USE m_gf_types
      USE m_gf_energies, ONLY: gf_NoEn, gf_Z, gf_allz
      USE m_hdf_accessprp
      USE m_constants
      USE m_gf_math
      USE m_gf_out

      IMPLICIT NONE
      PRIVATE 
      !<--These are public variables for the different types            
                                                                        
      INTEGER,PARAMETER :: IO2D_EMB      = 1 
      INTEGER,PARAMETER :: IO2D_GMAT     = 2 
      INTEGER,PARAMETER :: IO2D_TMAT     = 3 
      INTEGER,PARAMETER :: IO2D_CBS      = 4 
      INTEGER,PARAMETER :: IO2D_EMBADV   = 5 
      INTEGER,PARAMETER :: IO2D_G12      = 6 
      PUBLIC            :: IO2D_EMB,IO2D_GMAT,IO2D_TMAT,IO2D_CBS        &
     &     ,IO2D_EMBADV,IO2D_G12                                        
      INTEGER,PARAMETER :: IO2D_READ       = 2 
      INTEGER,PARAMETER :: IO2D_WRITE      = 1 
      INTEGER,PARAMETER :: IO2D_READWRITE  = 0 
      PUBLIC            :: IO2D_READ,IO2D_WRITE,IO2D_READWRITE 
                                                                        
      !>                                                                
      TYPE handle 
        INTEGER :: fileno 
        INTEGER :: mode 

        INTEGER :: region 
        INTEGER :: side 
        INTEGER :: rw 
        REAL    :: theta,phi 
        REAL    :: shift(2) 
        LOGICAL :: spinrot 
      END TYPE 
      TYPE t_file 
        CHARACTER(LEN = 10) :: filename 
        LOGICAL             :: writeable,no_hdf 
        INTEGER             :: matsize,spins,matscale 
        REAL                :: energyshift 
        INTEGER(HID_T)      :: fid,varid,mapid
        INTEGER,POINTER     :: emap(:) 
        INTEGER,POINTER     :: kmap(:) 
        INTEGER,POINTER     :: gmap(:) 
        INTEGER             :: supercell 
        REAL                :: eps_kpts 
        COMPLEX             :: eps_energy
      END TYPE

#ifdef CPP_MPI
      TYPE t_mem
        LOGICAL             :: in_mem=.false.
        LOGICAL             :: l_noco
        INTEGER             :: matsize,spins
        INTEGER             :: buffersize_sigma
        INTEGER             :: buffersize_tmat
        INTEGER,ALLOCATABLE :: mem_map(:,:)
        INTEGER             :: mem_handle_sigma
        INTEGER             :: mem_handle_tmat
        COMPLEX,ALLOCATABLE :: buffer_tmat(:),buffer_sigma(:)
        !COMPLEX,POINTER     :: buffer_tmat(:),buffer_sigma(:)
        INTEGER             :: stride_tmat,stride_sigma
      END TYPE
#endif
                                                                        
      TYPE(handle),ALLOCATABLE,SAVE   :: iohandle(:)
      TYPE(t_file),ALLOCATABLE,SAVE   :: iofile(:) 
#ifdef CPP_MPI
      TYPE(t_mem),save :: mem
#endif

      !cell and k-set metadata written as attributes of the 'kpts' dataset
      TYPE t_io2dmeta
        LOGICAL :: l_set = .FALSE.
        !cartesian in-plane lattice vectors (columns, a.u.)
        REAL    :: amat2d(2,2) = 0.0
        !cutoff of the 2D basis |k+G|<=rkmax
        REAL    :: rkmax = 0.0
        LOGICAL :: l_noco = .FALSE., l_soc = .FALSE.
        !time reversal used in the reduction of the k-set
        LOGICAL :: l_timerev = .FALSE.
        !distance of the plane to the outermost atom plane behind it
        !(region 0 side 2: layer 1, region nl+1 side 1: layer nl), <0 unknown
        REAL    :: zplane(2) = -1.0
        INTEGER :: nlayers = 0
        !energy reference: planar-averaged potential at the end of the vacuum
        !of the run (vacuum embedding), in the energy scale of the file
        LOGICAL :: l_vref = .FALSE.
        REAL    :: vref = 0.0
      END TYPE
      TYPE(t_io2dmeta),SAVE :: meta
      !number of handles from gf_io (the automatic handles follow)
      INTEGER,SAVE :: nuser_handles = 0
      PUBLIC :: t_io2dmeta

      PUBLIC init,finalize,gf_write2dmat,gf_read2dmat,gf_io2dmatFID     &
     &     ,gf_io2dmatshift,gf_io2dstatus,gf_io2dmat_memorysync         &
     &     ,gf_io2dmat_create,gf_io2dmat_kmesh,gf_io2dmat_userwrite
      CONTAINS 
                                                                        
      !<-- S: init                                                      
      SUBROUTINE init(kpts,gfinp,lapw_gf,sym,layers,jspins,l_noco,gmpi    &
     &     ,l_tmatmemory,cell,l_soc,zplanes,vref)
!-----------------------------------------------
!
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments

      INTEGER,INTENT(IN)        :: jspins
      TYPE(t_kpts),INTENT(IN)   :: kpts
      TYPE(t_embinp),INTENT(IN)  :: gfinp
      TYPE(t_lapw_gf),INTENT(IN)   :: lapw_gf
      TYPE(t_sym),INTENT(IN)    :: sym
      TYPE(t_layers),INTENT(IN) :: layers
      LOGICAL,INTENT(IN)        :: l_noco
      TYPE(t_gfmpi),intent(in)    :: gmpi
      LOGICAL,INTENT(OUT)       :: l_tmatmemory
      TYPE(t_cell),INTENT(IN),OPTIONAL :: cell
      LOGICAL,INTENT(IN),OPTIONAL      :: l_soc
      REAL,INTENT(IN),OPTIONAL         :: zplanes(2)
      REAL,INTENT(IN),OPTIONAL         :: vref
                                                                        
      !>                                                                
      !<-- Locals                                                       
                                                                        
      LOGICAL             :: lexist,bias 
      INTEGER             :: i,ii,n 
      INTEGER             :: iohandles 
      CHARACTER           :: iochar1,iochar2 
      REAL                :: kpts_eps,energy_shift,shift(2) 
      COMPLEX             :: en_eps 
      CHARACTER(LEN = 30) :: name 
      CHARACTER(LEN = 1)  :: mode,rw 
      CHARACTER(LEN = 7)  :: default
      LOGICAL             :: lead(2)
      INTEGER             :: region,side,id 
      LOGICAL             :: spinrot 
      REAL                :: theta,phi 
      INTEGER             :: foundfileno 
      !new for automatic io                                             
      INTEGER             :: layer,new_handles,fileno 
      TYPE(handle),ALLOCATABLE :: iohandle_tmp(:) 
      TYPE(t_file),ALLOCATABLE   :: iofile_tmp(:) 

      LOGICAL :: NO_HDF_TMP = .FALSE.
                                                                        
      NAMELIST /FILE/ name,energy_shift,bias,en_eps,kpts_eps 
      NAMELIST /IO/ mode,region,side,rw,id,spinrot,theta,phi,shift 
      if (l_noco) NO_HDF_TMP=.false.
      !>
      l_tmatmemory = .FALSE.
      meta%l_set = PRESENT(cell)
      IF (PRESENT(cell)) meta%amat2d = cell%amat(1:2,1:2)
      meta%rkmax = lapw_gf%rkmax
      meta%l_noco = l_noco
      meta%l_soc = .FALSE.
      IF (PRESENT(l_soc)) meta%l_soc = l_soc
      !gfprep reduces the k-set with time reversal unless SOC is on
      meta%l_timerev = .NOT. meta%l_soc
      IF (PRESENT(zplanes)) meta%zplane = zplanes
      meta%l_vref = PRESENT(vref)
      IF (PRESENT(vref)) meta%vref = vref
      meta%nlayers = layers%num_layers
      INQUIRE(FILE ='gf_io',EXIST = lexist)


      IF (.NOT.(lexist)) THEN
        !No IO wanted
        if (gfinp%l_addemb.and..not.(gfinp%l_vacleft.and.gfinp%l_surface)) &
           CALL juDFT_error("addemb set without gf_io-file")
        if (gfinp%l_writeup) CALL juDFT_error("writeup needs a gf_io file", &
           hint="add a handle &IO mode='E', region=0, side=2, rw='w' /")
        IF (layers%num_layers==1) RETURN
        !several layers: the temporary T-matrix/embedding storage is still needed
        iohandles = 0
        ALLOCATE(iohandle(0),iofile(0))
        GOTO 998
      ENDIF
      !<-- Now read the gf_io file                                      
      OPEN(120,FILE ='gf_io',FORM ='formatted',STATUS ='old') 

      !Check if defaults should be used
      read(120,"(a)") default
      rewind 120
      if (default.eq."default") then
         close(120)
         !the defaults are written to a scratch unit, the user's gf_io stays untouched
         open(120,FORM ='formatted',STATUS ='scratch',DELIM="APOSTROPHE")
         !a side with vacuum embedding has no lead file
         lead = (/.NOT.gfinp%l_vacleft,.NOT.gfinp%l_surface/)
         write(120,"(i10)") COUNT(lead)
         mode ="E";rw="r";id=0;spinrot=.FALSE.;theta=0.0;phi = 0.0;shift = 0.0
         if (gfinp%l_cbs) rw="w"
         if (lead(1)) then
            region=1;side=1;id=id+1
            write(120,IO)
         endif
         if (lead(2)) then
            region=layers%num_layers;side=2;id=id+1
            write(120,IO)
         endif
         energy_shift = 0.0 
         bias         = .FALSE. 
         kpts_eps     = 1E-5 
         en_eps       = CMPLX(1E-5,1E-5) 
         name         ='embpot11' 
         if (lead(1)) write(120,FILE)
         name         ='embpotN2' 
         if (lead(2)) write(120,FILE)
         rewind(120)
      endif


      !<-- first the IO-handles                                         
                                                                        
      READ(120,*) iohandles
      IF (iohandles == 0) THEN
         CLOSE(120)
         IF (layers%num_layers==1) RETURN
         ALLOCATE(iohandle(0),iofile(0))
         GOTO 998
      ENDIF
      ALLOCATE(iohandle(iohandles)) 
                                                                        
      DO i = 1,iohandles 
         !mode=E,region=1,side                                          
         !=1,io=r,file=12,spinrot=T,theta=0.000,phi=0.000               
         mode ="E";region= 1;side=1;rw="r";id=1;spinrot=.FALSE.;theta=0.0;phi = 0.0;shift = 0.0
         READ(120,IO) 
          iochar1 = mode;iohandle(i)%region = region;iohandle(i)%side   &
     &        = side;iochar2 = rw;iohandle(i)%fileno =id;iohandle(i     &
     &        )%spinrot = spinrot;iohandle(i)%theta=theta;iohandle(i    &
     &        )%phi = phi;iohandle(i)%shift=shift                       
         IF (iohandle(i)%fileno<1) CALL juDFT_error('illegal fileNO in gf_io')
         IF (iochar1 =='E'.OR.iochar1 =='e') THEN 
            iohandle(i)%mode = IO2D_EMB 
         ELSE IF (iochar1    =='A'.OR.iochar1 =='a') THEN 
            iohandle(i)%mode = IO2D_EMBADV 
         ELSE IF (iochar1    =='T'.OR.iochar1 =='t') THEN 
            iohandle(i)%mode = IO2D_TMAT 
         ELSE IF (iochar1    =='G'.OR.iochar1 =='g') THEN 
            iohandle(i)%mode = IO2D_GMAT 
         ELSE IF (iochar1    =='C'.OR.iochar1 =='c') THEN 
            iohandle(i)%mode = IO2D_CBS 
         ELSE IF (iochar1    =='N'.OR.iochar1 =='n') THEN
            iohandle(i)%mode = IO2D_G12 
         ELSE 
            WRITE (*,*) 'Unkown IO-mode line:',i,'Mode:',iochar1 
            CALL juDFT_error('gf_io2dmat')
         ENDIF 
         IF (iochar2 =='r'.OR.iochar2 =='R') THEN 
            iohandle(i)%rw = IO2D_READ 
         ELSE IF (iochar2    =='w'.OR.iochar2 =='W') THEN 
            iohandle(i)%rw = IO2D_WRITE 
         ELSE IF (iochar2 =='B'.OR.iochar2 =='b') THEN 
            iohandle(i)%rw = IO2D_READWRITE 
         ELSE 
            WRITE (*,*) 'Unkown IO-mode line:',i,'rw:',iochar2 
            CALL juDFT_error('gf_io2dmat')
         ENDIF 
      ENDDO 
                                                                        
      !>
      nuser_handles = iohandles
      !a plane between two layers is also read by the stack: a write handle
      !there becomes read/write
      DO i = 1,iohandles
         IF (iohandle(i)%mode/=IO2D_EMB.OR.iohandle(i)%rw/=IO2D_WRITE) CYCLE
         IF ((iohandle(i)%side==1.AND.iohandle(i)%region>=2.AND.iohandle(i)%region<=layers%num_layers).OR. &
     &       (iohandle(i)%side==2.AND.iohandle(i)%region>=1.AND.iohandle(i)%region<layers%num_layers)) THEN
            iohandle(i)%rw = IO2D_READWRITE
            IF (gmpi%pe0) WRITE(oUnit,"(a,i0,a,i0,a)") "gf_io: handle of region ",iohandle(i)%region, &
     &           ", side ",iohandle(i)%side," is between two layers and is used for reading as well"
         ENDIF
      ENDDO
      !check for required io-handles
      if (gfinp%l_writeup) THEN
         if (count(iohandle%mode==IO2D_EMB.and.iohandle%rw/=IO2D_READ.and.iohandle%region==0.and.iohandle%side==2)/=1) &
            CALL juDFT_error("writeup needs one gf_io write handle for region=0, side=2", &
                             hint="&IO mode='E', region=0, side=2, rw='w', id=<n> /")
      endif
      if (gfinp%l_addemb) THEN
         !is there an left embedding file?
         if (gfinp%l_vacleft) then
            if (count(iohandle%mode==IO2D_EMB.and.iohandle%rw==IO2D_READ.and.iohandle%region==1.and.iohandle%side==1)>0) &
                CALL juDFT_error("Left embedding potential supplied with addemb=T & vacleft=T")
         else
            if (count(iohandle%mode==IO2D_EMB.and.iohandle%rw==IO2D_READ.and.iohandle%region==1.and.iohandle%side==1)<1) &
               CALL juDFT_error("No left embedding potential supplied with addemb=T")
            if (count(iohandle%mode==IO2D_EMB.and.iohandle%rw==IO2D_READ.and.iohandle%region==1.and.iohandle%side==1)>1) &
               CALL juDFT_error("Multiple left embedding potentials supplied with addemb=T")
         endif
         if (gfinp%l_surface) then
            if (count(iohandle%mode==IO2D_EMB.and.iohandle%rw==IO2D_READ.and.iohandle%region==layers%num_layers.and.iohandle%side==2)>0) &
                CALL juDFT_error("Right embedding potential supplied with addemb=T & surface=T")
         else
            if (count(iohandle%mode==IO2D_EMB.and.iohandle%rw==IO2D_READ.and.iohandle%region==layers%num_layers.and.iohandle%side==2)>1) &
                CALL juDFT_error("Multiple right embedding potentials supplied with addemb=T")
            if (count(iohandle%mode==IO2D_EMB.and.iohandle%rw==IO2D_READ.and.iohandle%region==layers%num_layers.and.iohandle%side==2)<1) &
                CALL juDFT_error("No right embedding potentials supplied with addemb=T")
         endif
      endif
      !<-- Now read the filenames, the matrix-sizes and the energyshift!
      ALLOCATE(iofile(MAXVAL(iohandle%fileno))) 
      iofile%matsize=0 
      iofile%no_hdf = .FALSE. 
      DO i = 1,MAXVAL(iohandle%fileno) 
         !<-- defaults                                                  
                                                                        
         name         ='IO2D_2DMAT' 
         energy_shift = 0.0 
         bias         = .FALSE. 
         kpts_eps     = 1E-5 
         en_eps       = CMPLX(1E-5,1E-5) 
                                                                        
         !>                                                             
         READ(120,FILE,ERR = 999,END = 999) 
         iofile(i)%filename    = name 
         iofile(i)%energyshift = energy_shift 
         iofile(i)%eps_kpts=kpts_eps 
         iofile(i)%eps_energy=en_eps 
         IF (bias) THEN 
            CALL juDFT_error("Bias in gf_io2dmat not used anymore")
            IF (energy_shift /= 0.0)                                    &
     &           CALL juDFT_error("Both bias+energy_shift in gf_io")
!            iofile(i)%energyshift = -1.*gfinp%bias ! bias has to be app
                                                   ! right side         
         ENDIF 
         !<-- Check if file must be writeable ...                       
         ! ... and perform consistency check.                           
         iofile(i)%writeable = .FALSE. 
         foundfileno=0 
         DO ii = 1,SIZE(iohandle) 
            IF (iohandle(ii)%fileno == i) THEN 
               foundfileno=foundfileno+1 
               IF(iohandle(ii)%rw == IO2D_WRITE.OR.iohandle(ii)%rw == IO2D_READWRITE)     &
     &             iofile(i)%writeable = .TRUE.                         
                  !fileno=i                                             
            ENDIF 
               !loop over iohandles                                     
         ENDDO 
         IF(foundfileno==0) CALL juDFT_error("foundnofileno")
!         if(foundfileno.ne.1)call juDFT_error("foundmanyfileno")         
                                                                        
                                                                        
                                                                        
         !>                                                             
         !<-- determine the size-factor of the matrix                   
         iofile(i)%matsize=0 
         DO ii = 1,SIZE(iohandle) 
          IF(iohandle(ii)%fileno == i)THEN 
            IF (iohandle(ii)%mode == IO2D_TMAT.OR.iohandle(ii)%mode== IO2D_CBS) THEN
               IF (iofile(i)%matsize /= 0.AND.iofile(i)%matsize /= 2)   &
     &              CALL juDFT_error("Conflicting file sizes in gf_io")
               iofile(i)%matsize = 2 
                  !transfer matrix: large matrices                      
            ENDIF 
            IF (iohandle(ii)%mode == IO2D_EMB .OR.                      &
     &           iohandle(ii)%mode == IO2D_GMAT .OR.                    &
     &           iohandle(ii)%mode == IO2D_EMBADV .OR.                  &
     &           iohandle(ii)%mode == IO2D_G12) THEN                    
               IF (iofile(i)%matsize /= 0.AND.iofile(i)%matsize /= 1)   &
     &              CALL juDFT_error("Conflicting file sizes in gf_io")
               iofile(i)%matsize = 1 
                  !embedding matrix, surface-green, ....: small matrices
            ENDIF 
               !fileno=i                                                
          ENDIF 
               !loop over iohandles                                     
         ENDDO 
         !>                                                             
      ENDDO 
  999 CLOSE(120) 
  998 CONTINUE
                                                                        
      !>                                                                
                                                                        
      !<-- Some checks...                                               
      IF (ANY(iofile%matsize==0)) THEN 
         WRITE(*,*) 'ERROR in reading gf_io' 
         WRITE(*,*) 'You tried to assign an IO to an unkown filenumber' 
         CALL juDFT_error('gf_io')
      ENDIF 
      !>                                                                
      !>                                                                
                                                                        
      !<-- add IO for tempoary files in case of multiple layers         
                                                                        
      fileno = SIZE(iofile)
      !writeup and the output of the top plane of a surface propagate also
      !through a single layer: keep its T-matrix
      IF (layers%num_layers==1) THEN
         IF (gfinp%l_writeup.OR.(gfinp%l_surface.AND.ANY(iohandle%mode==IO2D_EMB.AND.iohandle%region==2 &
     &        .AND.iohandle%side==1.AND.iohandle%rw/=IO2D_READ))) &
     &      l_tmatmemory = (COUNT(iohandle%mode==IO2D_TMAT)==0)
      ENDIF
      IF (layers%num_layers>1) THEN
         !<-- new handles for t-matrix

         new_handles = layers%num_layers-COUNT(iohandle%mode ==         &
     &        IO2D_TMAT)
         !check if t-matrix should be stored in memory
         l_tmatmemory = (new_handles == layers%num_layers)
         IF (new_handles>0.AND.(new_handles /= layers%num_layers)) THEN 
            ALLOCATE(iohandle_tmp(iohandles)) 
            iohandle_tmp = iohandle 
            DEALLOCATE(iohandle) 
            ALLOCATE(iohandle(iohandles+new_handles)) 
            iohandle(:iohandles) = iohandle_tmp 
            DEALLOCATE(iohandle_tmp) 
            DO layer = 1,layers%num_layers 
               IF (ANY(iohandle(:iohandles)%mode == IO2D_TMAT.AND.       &
     &              iohandle(:iohandles)%region == layer)) CYCLE
               iohandles = iohandles+1 
               iohandle(iohandles)%mode    = IO2D_TMAT 
               iohandle(iohandles)%rw      = IO2D_READWRITE 
               iohandle(iohandles)%region  = layer 
               iohandle(iohandles)%side    = 0 
               fileno                      = fileno+1 
               iohandle(iohandles)%fileno  = fileno 
               iohandle(iohandles)%spinrot = .FALSE. 
               iohandle(iohandles)%theta   = 0.0 
               iohandle(iohandles)%phi     = 0.0 
               iohandle(iohandles)%shift   = 0.0 
            ENDDO 
            !<-- new filehandles for t-matrices                         
                                                                        
            n = SIZE(iofile) 
            ALLOCATE(iofile_tmp(n)) 
            iofile_tmp = iofile 
            DEALLOCATE(iofile) 
            ALLOCATE(iofile(new_handles+n)) 
            iofile(:n)=iofile_tmp 
            DEALLOCATE(iofile_tmp) 
            DO i = n+1,n+new_handles 
               WRITE(iofile(i)%filename,"(a3,i7.7)") "_t_",i 
               iofile(i)%energyshift = 0.0 
               iofile(i)%eps_kpts    = 1.0E-5 
               iofile(i)%eps_energy  = CMPLX(1.0E-5,1.0E-5) 
               iofile(i)%writeable   = .TRUE. 
               iofile(i)%matsize     = 2 
               iofile(i)%no_hdf      = .FALSE. 
!               iofile(i)%no_hdf = NO_HDF_TMP                           
            ENDDO 
                                                                        
            !>                                                          
         ENDIF 
                                                                        
         !>                                                             
                                                                        
         !<-- new handles for embpot

         !count the missing handles: user handles of other planes (region 0,
         !nlayers+1) must not reduce the number
         new_handles = 0
         DO layer = 1,layers%num_layers
            DO side = 1,2
               IF (gfinp%l_surface.AND.layer == layers%num_layers.AND.side == 2) CYCLE
               IF (ANY(iohandle%mode == IO2D_EMB.AND.iohandle%region == layer &
     &              .AND.iohandle%side == side)) CYCLE
               new_handles = new_handles+1
            ENDDO
         ENDDO
         IF (gfinp%l_surface) THEN
            IF (.NOT.ANY(iohandle%mode == IO2D_EMB.AND.iohandle%region == layers%num_layers+1 &
     &           .AND.iohandle%side == 1)) new_handles = new_handles+1
         ENDIF

#ifdef CPP_MPI
         IF(gfinp%l_surface) THEN
           mem%in_mem=(COUNT(iohandle%mode ==IO2D_EMB)==1)
         ELSE
           mem%in_mem=(count(iohandle%mode==IO2D_EMB)==2)
         ENDIF
         if (mem%in_mem)THEN
           IF(priv_create_mpi_memory(gmpi,layers,kpts,lapw_gf,jspins,l_noco)) THEN
             new_handles=0
             l_tmatmemory=.false.
           ENDIF
         ENDIF
#endif


         IF (new_handles>0) THEN
               ALLOCATE(iohandle_tmp(iohandles)) 
               iohandle_tmp = iohandle 
               DEALLOCATE(iohandle) 
               ALLOCATE(iohandle(iohandles+new_handles)) 
               iohandle(:iohandles) = iohandle_tmp 
               DO layer = 1,layers%num_layers 
                  DO side  = 1,2 
                     IF (gfinp%l_surface.AND.layer == layers%num_layers &
     &                    .AND.side == 2) CYCLE                         
                     IF (ANY(iohandle(:iohandles)%mode == IO2D_EMB      &
     &                    .AND.iohandle(:iohandles)%region == layer     &
     &                    .AND.iohandle(:iohandles)%side      == side)) &
     &                    CYCLE                                         
                     iohandles = iohandles+1 
                     iohandle(iohandles)%mode = IO2D_EMB 
                     iohandle(iohandles)%rw   = IO2D_READWRITE 
                     iohandle(iohandles)%region = layer 
                     iohandle(iohandles)%side   = side 
                     fileno                     = fileno+1 
                     iohandle(iohandles)%fileno = fileno 
                     iohandle(iohandles)%spinrot = .FALSE. 
                     iohandle(iohandles)%theta   = 0.0 
                     iohandle(iohandles)%phi     = 0.0 
                     iohandle(iohandles)%shift   = 0.0 
                  ENDDO 
               ENDDO 

               if (gfinp%l_surface) THEN
               		!Surface embedding potential must be saved as well
                    if (.not.(any(iohandle(:iohandles)%mode == IO2D_EMB      &
     &                    .AND.iohandle(:iohandles)%region == layers%num_layers+1     &
     &                    .AND.iohandle(:iohandles)%side      == 1))) THEN
                       iohandles=iohandles+1
                       iohandle(iohandles)%mode = IO2D_EMB
                       iohandle(iohandles)%rw   = IO2D_READWRITE
                       iohandle(iohandles)%region = layers%num_layers+1
                       iohandle(iohandles)%side   = 1
                       fileno                     = fileno+1
                       iohandle(iohandles)%fileno = fileno
                       iohandle(iohandles)%spinrot = .FALSE.
                       iohandle(iohandles)%theta   = 0.0
                       iohandle(iohandles)%phi     = 0.0
                       iohandle(iohandles)%shift   = 0.0
                    endif
                endif
               !<-- new filehandles for emb-matrices                    
               n = SIZE(iofile) 
               ALLOCATE(iofile_tmp(n)) 
               iofile_tmp = iofile 
               DEALLOCATE(iofile) 
               ALLOCATE(iofile(new_handles+n)) 
               iofile(:n) = iofile_tmp 
               DEALLOCATE(iofile_tmp) 
               DO i = n+1,n+new_handles 
                  WRITE(iofile(i)%filename,"(a3,i7.7)") "_e_",i 
                  iofile(i)%energyshift = 0.0 
                  iofile(i)%eps_kpts    = 1.0E-5 
                  iofile(i)%eps_energy  = CMPLX(1.0E-5,1.0E-5) 
                  iofile(i)%writeable   = .TRUE. 
                  iofile(i)%matsize     = 1 
                  iofile(i)%no_hdf     = NO_HDF_TMP 
               ENDDO 
               !>                                                       
         ENDIF 
      ENDIF 
                                                                        
      !>                                                                
      !<-- Write out some information                                   
      IF (gmpi%pe0) THEN
         WRITE(oUnit,*) "io2Dmat-Setup" 
         WRITE(oUnit,*) "Handles:",SIZE(iohandle) 
          WRITE(oUnit,*) "No File Region Side RW Spinrot     Theta   Phi" 
         DO i = 1,SIZE(iohandle) 
            WRITE(oUnit,'(i2,2x,i2,4x,i2,5x,i2,3x,i2,3x,l1,3x,2(f10.8,1x))')&
     &           i,iohandle(i)%fileno,iohandle(i)%region,iohandle(i     &
     &           )%side,iohandle(i)%rw,iohandle(i)%spinrot,iohandle(i   &
     &           )%theta,iohandle(i)%phi                                
                                                                        
         ENDDO 
         WRITE(oUnit,*) 
         WRITE(oUnit,*) "Files:",SIZE(iofile) 
         WRITE(oUnit,*) "No  Filename  Write           Size" 
         DO i = 1,SIZE(iofile) 
            WRITE(oUnit,'(i2,3x,a10,2x,l1,6x,i8)')i,iofile(i)%filename      &
     &           ,iofile(i)%writeable,iofile(i)%matsize                 
         ENDDO 
      ENDIF 
      !>                                                                
      !<-- Now open all              
      DO i = 1,SIZE(iofile) 
         CALL priv_open_file(iofile(i),lapw_gf,kpts,sym,jspins,l_noco) 
      ENDDO 
      !>

      END SUBROUTINE 
      !>                                                                





                                                                  
      !<-- S: priv_open_file(io_file,lapw_gf,kpts,sym,jspins,l_noco)       
                                                                        
      SUBROUTINE priv_open_file(io_file,lapw_gf,kpts,sym,jspins,l_noco) 
!-----------------------------------------------                        
!                                                                       
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_file),INTENT(INOUT) :: io_file 
      TYPE(t_lapw_gf),INTENT(IN)  :: lapw_gf 
      TYPE(t_kpts),INTENT(IN)  :: kpts 
      TYPE(t_sym),INTENT(IN)   :: sym 
      INTEGER,INTENT(IN)       :: jspins 
      LOGICAL,INTENT(IN)       :: l_noco 
                                                                        
      !>                                                                
      !<-- Locals                                                       
      LOGICAL              :: l_exist 
      INTEGER              :: n,k,nn,kk,nk 
      COMPLEX,ALLOCATABLE  :: en_in(:) 
      REAL,ALLOCATABLE     :: kpts_in(:,:) 
      INTEGER,ALLOCATABLE  :: gvecs_in(:,:) 
                                                                        
                                                                        
      INTEGER(hid_t)   :: varid,fspace
      INTEGER          :: access_mode
      INTEGER(hsize_t) :: dims(7),maxdims(7)
      INTEGER          :: hdferr,irank
      TYPE(t_io2dmeta) :: fmeta
      !>
      !<--Access properties for MPI

#ifdef CPP_MPI                                                          
      integer :: ierr
      CALL  MPI_COMM_RANK (MPI_COMM_WORLD,irank,hdferr)
#else                                                                   
      irank=0 
#endif
                                                                        
      CALL timestart("io2dmat open:"//io_file%filename)
      INQUIRE(FILE ="supercell_io",EXIST= l_exist)
      IF (l_exist) CALL juDFT_error("supercell_io is no longer supported", &
     &     hint="use incomm_inp (&INCOMM) to map an embedding potential to a supercell")
                                                                        
                                                                        
      INQUIRE(FILE = TRIM(io_file%filename)//'.hdf',EXIST=l_exist) 
                                                                        
      !<--delete file if exist and filename starts with an underscore   
      IF (l_exist.AND.io_file%filename(1:1)=='_') THEN 
         IF(irank==0) THEN 
            OPEN(99,FILE=TRIM(io_file%filename)//'.hdf') 
            CLOSE(99,STATUS="delete") 
            WRITE(*,*) "Deleted "//TRIM(io_file%filename)//".hdf" 
         ENDIF 
         l_exist=.FALSE. 
      ENDIF 

      !>                                                                
      IF (.NOT.l_exist.AND..NOT.io_file%writeable) CALL juDFT_error("gf_io: file for reading does not exist")
      IF (l_exist) THEN 
         !<--Check if file can be opened read-only and open file        
         IF (io_file%writeable) THEN 
            access_mode = H5F_ACC_RDWR_F 
         ELSE 
            access_mode = H5F_ACC_RDONLY_F 
         ENDIF 
         CALL io_hdfopen (TRIM(io_file%filename)//'.hdf',access_Mode,    &
     &        io_file%fid, hdferr,hdf_access_prp(trim(io_file%filename)))
         ! Open dataset
         CALL timestart("open dataset")
         CALL io_dopen(io_file%fid, '2Dmat', io_file%varid, hdferr)
         CALL timestop("open dataset")
         !>
         !Open dataset for gmap
         CALL timestart("open gmap")
         !CALL io_dopen(io_file%fid, 'mapping', io_file%mapid, hdferr)
         CALL timestop("open gmap")

         !>
         !<--Read energies                                              
         CALL timestart("energy+k mapping")
                                                                        
         CALL io_dopen(io_file%fid, 'energies', varid, hdferr) 
         CALL h5dget_space_f(varid,fspace,hdferr) 
         CALL h5sget_simple_extent_dims_f(fspace,dims,maxdims,hdferr) 
         ALLOCATE(en_in(dims(2))) 
         CALL h5sclose_f(fspace,hdferr) 
         CALL io_READ(varid,(/-1,1/),(/1,SIZE(en_in)/),"en_in",en_in)
         CALL io_dclose(varid,hdferr) 
         !generate mapping index                                        
         ALLOCATE(io_file%emap(gf_NoEN())) 
         DO n = 1,gf_NoEN() 
            !find correct energy                                        
            io_file%emap(n) =-1 
            DO nn = 1,SIZE(en_in) 
               IF (ABS(real(gf_Z(n,0))+io_file%energyshift-real(en_in(nn&
     &              )))<REAL(io_file%eps_energy).AND.ABS(aimag(gf_z(n,0)&
     &              )-AIMAG(en_in(nn)))<AIMAG(io_file%eps_energy)       &
     &              )io_file%emap(n) = nn                               
            ENDDO 
            IF (io_file%emap(n) ==-1) THEN 
               WRITE(*,*) io_file%filename 
               WRITE(*,"('Energy:',2f10.5,'->',2f10.5)") gf_Z(n,0)      &
     &              ,gf_Z(n,0)+io_file%energyshift                      
               WRITE(*,*) 'In hdf-file:' 
               WRITE(*,"(999(2f10.5,/))") (en_in(nn),nn = 1,SIZE(en_in))
               write(*,*) "Energy-epsilon:", io_file%eps_energy
               CALL juDFT_error('No matching energy found')
            ENDIF 
         ENDDO 
         DEALLOCATE(en_in) 
                                                                        
         !>                                                             
         !<-- Read kpts                                                 
                                                                        
         CALL io_dopen(io_file%fid, 'kpts', varid, hdferr) 
         CALL h5dget_space_f(varid,fspace,hdferr) 
         CALL h5sget_simple_extent_dims_f(fspace,dims,maxdims,hdferr) 
         ALLOCATE(kpts_in(2,dims(2))) 
         CALL h5sclose_f(fspace,hdferr) 
         CALL io_READ(varid,(/1,1/),(/2,SIZE(kpts_in,2)/),"kpts_in",kpts_in)
         CALL io_dclose(varid,hdferr) 
         !generate mapping index                                        
         ALLOCATE(io_file%kmap(kpts%nkpt)) 
         DO k = 1,kpts%nkpt 
            !find kpt                                                   
            io_file%kmap(k) =-1 
            DO kk = 1,SIZE(kpts_in,2) 
               IF(ALL(ABS(kpts%bk(1:2,k)-kpts_in(:,kk))<io_file%eps_kpts&
     &              ))io_file%kmap(k) = kk                              
            ENDDO 
            IF (io_file%kmap(k) ==-1) THEN 
               WRITE(*,*) io_file%filename 
               WRITE(*,*) 'kpt:',kpts%bk(:,k) 
               WRITE(*,*) 'In hdf-file:' 
               DO kk = 1,SIZE(kpts_in,2) 
                  WRITE(*,*) kpts_in(:,kk) 
               ENDDO 
               CALL juDFT_error('No matching k-point found')
            ENDIF 
         ENDDO 
         DEALLOCATE(kpts_in) 
                                                                        
         !>                                                             
         !<--Now read the gmap and create the mapping index             
                                                                        
         CALL io_dopen(io_file%fid, 'gvecs', varid, hdferr) 
         CALL h5dget_space_f(varid,fspace,hdferr) 
         CALL h5sget_simple_extent_dims_f(fspace,dims,maxdims,hdferr) 
         ALLOCATE(gvecs_in(2,dims(2))) 
         CALL h5sclose_f(fspace,hdferr) 
         CALL io_READ(varid,(/1,1/),(/2,SIZE(gvecs_in,2)/),"gvecs_in",gvecs_in)
         CALL io_dclose(varid,hdferr) 
         !create mapping index                                          
         !gmap maps the index in the current global 2D list onto the index in
         !the file (0 if the G-vector is not in the file)
         ALLOCATE(io_file%gmap(SIZE(lapw_gf%global2Dlist,2))) 
         io_file%gmap = 0 
         DO n = 1,SIZE(gvecs_in,2) 
            IF (ALL(ABS(gvecs_in(1:2,n)) <= lapw_gf%g_max(1:2))) THEN
               nn = lapw_gf%global2Dmap(gvecs_in(1,n),gvecs_in(2,n))
               IF (nn>0) io_file%gmap(nn) = n
            ENDIF
         ENDDO 
         !up to now io_file%matsize was the scaling factor              
         io_file%matscale=io_file%matsize 
         io_file%matsize=io_file%matsize*SIZE(gvecs_in,2) 
         DEALLOCATE(gvecs_in) 
                                                                        
         !>                                                             
         !<-- now check the available amount of data                    
         CALL h5dget_space_f(io_file%varid,fspace,hdferr) 
         CALL h5sget_simple_extent_dims_f(fspace,dims,maxdims,hdferr) 
         CALL h5sclose_f(fspace,hdferr) 
         IF (io_file%matsize/=dims(1)) THEN 
            WRITE (*,*) 'In file:',io_file%filename 
            WRITE (*,*) 'Matrix dimension missmatch' 
            CALL juDFT_error('gf_io2dmat')
         ENDIF 
         !check how many spin's there are...                            
         io_file%spins=dims(4)
         CALL timestop("energy+k mapping")

         !>                                                             
      ELSE !file does not exist
         !<-- create the dimensions                                     
         io_file%matscale=io_file%matsize 
         io_file%matsize=io_file%matsize*SIZE(lapw_gf%global2Dlist,2) 
         io_file%spins=jspins 
         IF (l_noco) io_file%spins=4 
         !>                                                             
         !<-- Create the file and dataset                               
         IF (io_file%no_hdf) THEN 
            !find a new unit                                            
            l_exist = .TRUE. 
            io_file%fid = 777 
            DO WHILE(l_exist) 
               io_file%fid = io_file%fid+1 
               INQUIRE(UNIT = io_file%fid,OPENED=l_exist) 
            ENDDO 
            OPEN(io_file%fid,ACCESS ='direct',RECL = io_file%matsize**2 &
     &           *16,STATUS ='scratch')                                 
         ELSE
            nk = kpts%nkpt
#ifndef CPP_MPI
            if (io_file%filename(1:1)=='_') nk=1
#endif
            !the plane this file describes fixes zplane
            fmeta = meta
            fmeta%zplane = -1.0
            DO n = 1,SIZE(iohandle)
               IF (iohandle(n)%fileno<1.OR.iohandle(n)%fileno>SIZE(iofile)) CYCLE
               IF (iofile(iohandle(n)%fileno)%filename/=io_file%filename) CYCLE
               IF (iohandle(n)%mode/=IO2D_EMB) CYCLE
               IF (iohandle(n)%region==0.AND.iohandle(n)%side==2) fmeta%zplane(1) = meta%zplane(1)
               IF (iohandle(n)%region==meta%nlayers+1.AND.iohandle(n)%side==1) fmeta%zplane(1) = meta%zplane(2)
            ENDDO
            CALL gf_io2dmat_create(TRIM(io_file%filename),lapw_gf,kpts,sym,  &
     &           io_file%spins,nk,io_file%matsize,fmeta,io_file%fid,io_file%varid)
            ALLOCATE(io_file%emap(gf_NoEN()))
            io_file%emap = (/(n,n = 1,gf_NoEn())/)
            ALLOCATE(io_file%kmap(kpts%nkpt))
            io_file%kmap = (/ (n,n = 1,kpts%nkpt ) /)
            if (nk/=kpts%nkpt) io_file%kmap=1
            ALLOCATE(io_file%gmap(SIZE(lapw_gf%global2Dlist,2)))
            io_file%gmap = (/(n,n = 1,SIZE(io_file%gmap))/)
         ENDIF
      ENDIF 
      CALL timestop("io2dmat open:"//io_file%filename)
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S:finalize()                                                 
                                                                        
      SUBROUTINE finalize() 
!********************************************************************** 
!     All hdf-files are closed, the arrays are not deallocated so this  
!     will not allow to OPEN the files again                            
!********************************************************************** 
      IMPLICIT NONE
      INTEGER :: i,hdferr
      !io2dmat_init returns early (and leaves iofile unallocated) when no
      !gf_io file / no io-handles are present - nothing to close then
      IF (.NOT.ALLOCATED(iofile)) RETURN
      DO i = 1,size(iofile)
         IF (iofile(i)%no_hdf) THEN 
            CLOSE(iofile(i)%fid) 
         ELSE 
            !CALL io_dclose(iofile(i)%mapid,hdferr)
            CALL io_dclose(iofile(i)%varid,hdferr) 
            CALL io_hdfclose(iofile(i)%fid,hdferr) 
         ENDIF 
      ENDDO 
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S:gf_write2dmat(mode,region,side,en,nk,jspin,lapw_gf,mat)       
      SUBROUTINE gf_write2dmat(mode,region,side,en,nk,jspin   &
     &     ,lapw_gf,mat)
!********************************************************************** 
!     * This SUBROUTINE saves the 2DMatrix for the                      
!     * energy,kpoint,spin                                              
!     *                           Daniel Wortmann, Tokyo, 2002          
!********************************************************************** 
      IMPLICIT NONE 
      INTEGER,INTENT(IN)          :: mode,region,side 
      INTEGER,INTENT(IN)          :: jspin,nk,en 
      TYPE(t_lapw_gf),INTENT(IN)     :: lapw_gf 
      COMPLEX,INTENT(IN)          :: mat(:,:) 
                                                                        
      INTEGER :: fileid
      REAL   ,ALLOCATABLE :: A(:,:,:,:)

#ifdef CPP_MPI
      if (mem%in_mem.and.(mode==IO2D_tmat.or.mode==IO2D_EMB)) CALL  mpi_memory_write(mode,region,side,en,nk,jspin,mat)
#endif

                                                                        
      !First get the fileid                                             
      fileid = priv_getFILEHANDLE(IO2D_WRITE,mode,region,side) 
      IF (fileid<1) RETURN
      CALL timestart("gf_write2dmat")

      !Check if this is a noco calculation
      if (jspin==1.and.size(mat,1)>iofile(fileid)%matscale*lapw_gf%nv2(1)) then
          allocate(A(iofile(fileid)%matsize,iofile(fileid)%matsize,2,4))
      else
          allocate(A(iofile(fileid)%matsize,iofile(fileid)%matsize,2,1))
      endif
      A=0
      !Now map data to output buffer

      call priv_map_out(jspin,lapw_gf,fileid,mat,A)

      !Write Data to file
       if (size(a,4)==1) then
            CALL io_WRITE(iofile(fileid)%varid,(/1,1,1,jspin,iofile(fileid)%emap(en),iofile(fileid)%kmap(nk)/),&
     &     (/SIZE(a,1),SIZE(a,2),2,1,1,1/),"A",&
     &     A(:,:,:,1))

       else !in noco-case write all spins at once
            CALL io_WRITE(iofile(fileid)%varid,(/1,1,1,1,iofile(fileid)%emap(en),iofile(fileid)%kmap(nk)/),&
     &     (/SIZE(a,1),SIZE(a,2),2,4,1,1/),"A",&
     &     A(:,:,:,:))

       endif
       deallocate (A)
       CALL timestop("gf_write2dmat")

      end subroutine

      RECURSIVE SUBROUTINE priv_map_out(jspin,lapw_gf,fileid,mat,A)
      !This subroutine maps the 2D-matrix to te output buffer
      IMPLICIT NONE
      INTEGER,INTENT(IN)          :: jspin,fileid
      TYPE(t_lapw_gf),INTENT(IN)     :: lapw_gf
      COMPLEX,INTENT(IN)          :: mat(:,:)
      real,intent(inout)          :: a(:,:,:,:)
                                                                        
      INTEGER :: nv2,na,jsp,i,ii,n1,n2

      !--check if matrix cosists of several blocks, e.g. T-Matrix
        IF (iofile(fileid)%matscale == 2.AND.SIZE(mat,1)>lapw_gf%nv2_tot) THEN
         nv2 = lapw_gf%nv2_tot
         na=size(a,1)/2
         call priv_map_out(jspin,lapw_gf,fileid,mat(1:nv2,1:nv2),a(1:na,1:na,:,:))
         call priv_map_out(jspin,lapw_gf,fileid,mat(nv2+1:,1:nv2),a(1:na,na+1:,:,:))
         call priv_map_out(jspin,lapw_gf,fileid,mat(1:nv2,nv2+1:),a(na+1:,1:na,:,:))
         call priv_map_out(jspin,lapw_gf,fileid,mat(nv2+1:,nv2+1:),a(na+1:,na+1:,:,:))
         RETURN 
      ENDIF 
                                                                        
                                                                        
      !<--Matrix might be a noco-matrix containing four components check
      !and call writing subroutine recursive in this case               
      IF (jspin == 1) THEN 
         IF (SIZE(mat,1)>lapw_gf%nv2(1)) THEN 
            !This is a noco calculation                                 
            nv2 = lapw_gf%nv2_tot/2
            call priv_map_out(1,lapw_gf,fileid,mat(1:nv2,1:nv2),a)
            call priv_map_out(2,lapw_gf,fileid,mat(nv2+1:,nv2+1:),a)
            call priv_map_out(3,lapw_gf,fileid,mat(1:nv2,nv2+1:),a)
            call priv_map_out(4,lapw_gf,fileid,mat(nv2+1:,1:nv2),a)
            RETURN 
         ENDIF 
      ENDIF 
      !>                                                                
                                                                        

      IF (iofile(fileid)%no_hdf) &
         CALL juDFT_error("Non-HDF-file no longer supported in gf_io2dmat")

                                                                        
                                                                        
      IF (iofile(fileid)%spins<jspin) THEN 
         WRITE(*,*) "In file jspin=",iofile(fileid)%spins 
         WRITE(*,*) "In this setup jspin=",jspin 
         CALL juDFT_error('gf_io2dmat: Wrong no of spins in writing')
      ENDIF 
      !<-- Check size                                                   
      IF (SIZE(mat,2) /= SIZE(mat,1))                                   &
     &     CALL juDFT_error('gf_io2dmat: No square matrix')
                                                                        
      !>                                                                
                                                                        
      !<-- now map data to temp array                                   


      jsp=1
      if (size(a,4)>2) jsp=jspin

      DO i = 1,lapw_gf%nv2(min(jspin,2))
         n1=lapw_gf%g2map(i)
         if (n1==0) CALL juDFT_error("Uups")
         IF (iofile(fileid)%gmap(n1) == 0.OR.iofile(fileid)%gmap(n1)>SIZE(A,1)) CYCLE
         DO ii = 1,lapw_gf%nv2(min(jspin,2))
            n2=lapw_gf%g2map(ii)
            if (n2==0) cycle
            IF (iofile(fileid)%gmap(n2) == 0.OR.iofile(fileid)%gmap(n2  &
     &           )>SIZE(a,2)) CYCLE
             A(iofile(fileid)%gmap(n1),iofile(fileid)%gmap(n2),1,jsp)=real(mat(i,ii))
             A(iofile(fileid)%gmap(n1),iofile(fileid)%gmap(n2),2,jsp)=aimag(mat(i,ii))
         ENDDO 
      ENDDO 
      !>                                                                
                                                                        
      END SUBROUTINE 
                                                                        
      !>
      SUBROUTINE priv_read_files_to_mem(layers,jspins,nkpts,lapw_gf,l_noco)
      IMPLICIT NONE
      Type(t_layers),intent(in) :: layers
      TYPE(t_lapw_gf),intent(in)   :: lapw_gf
      integer,intent(in)        :: jspins,nkpts
      logical,intent(in)        :: l_noco

      INTEGER:: jsp,l,s,jspin,en,nk
      COMPLEX,ALLOCATABLE:: mat(:,:)

      ALLOCATE(mat(lapw_gf%nv2d,lapw_gf%nv2d))
      jsp=jspins
      if (l_noco) jsp=1
#ifdef CPP_MPI
      DO l=1,layers%num_layers
         DO s=1,2
            IF (priv_getFILEHANDLE(IO2D_READ,IO2D_EMB,l,s)<1) cycle
            !OK we will read this embedding potential
            DO jspin=1,jsp
                DO en=1,gf_noEN()
                    DO nk=1,nkpts
                          mem%in_mem=.false.
                          mem%in_mem=gf_read2dmat(IO2D_EMB,l,s,en,nk,jspin,lapw_gf,mat)
                          CALL gf_write2dmat(IO2D_EMB,l,s,en,nk,jspin,lapw_gf,mat)
                    ENDDO
                ENDDO
            ENDDO
         ENDDO
      ENDDO
      mem%in_mem=.true.
#endif
      END subroutine

      !<-- F:gf_read2dmat(mode,region,side,en,nk,jspin,lapw_gf,mat)        
      RECURSIVE FUNCTION gf_read2dmat(mode,region,side,en,nk,jspin,lapw_gf &
     &     ,mat,pos1,pos2) RESULT(read2dmat)                            
!********************************************************************** 
!     * This SUBROUTINE reads the 2DMatrix for the                      
!     * energy,kpoint,spin                                              
!     * It also performs various transformations if needed:             
!     * Constructs a magnetic matrix from a non-magnetic                
!     * Rotates the spin quantisation axis                              
!     *                           Daniel Wortmann, Tokyo, 2002          
!********************************************************************** 
      IMPLICIT NONE 
      INTEGER,INTENT(IN)         :: mode,region,side 
      INTEGER,INTENT(IN)         :: jspin,nk,en 
      TYPE(t_lapw_gf),INTENT(IN)    :: lapw_gf 
      COMPLEX,INTENT(OUT)        :: mat(:,:) 
      INTEGER,INTENT(IN),OPTIONAL :: pos1,pos2 
      LOGICAL                     :: read2dmat 
                                                                        
      INTEGER :: fileid,nv2,spin,i,ii,start1,start2,n1,n2
      REAL    :: s(2) 
      REAL   ,ALLOCATABLE :: A(:,:,:) 
      COMPLEX,ALLOCATABLE :: phase(:,:)
      LOGICAL             :: l 
      read2dmat=.FALSE. 

#ifdef CPP_MPI
      if (mem%in_mem.and.(mode==IO2D_EMB.or.mode==IO2D_tmat)) THEN
           read2dmat=.true.
           CALL  mpi_memory_read(mode,region,side,en,nk,jspin,mat)
      ENDIF
#endif
      !First get the fileid                                             
      fileid = priv_getFILEHANDLE(IO2D_READ,mode,region,side) 

      IF (fileid<1) RETURN 
      CALL timestart("gf_read2dmat")
      !<-- check if matrix cosists of several blocks                    
      IF (.NOT.PRESENT(pos1).AND.iofile(fileid)%matscale == 2           &
     &     .AND.SIZE(mat,1)>lapw_gf%nv2_tot) THEN                          
         nv2 = lapw_gf%nv2_tot 
         l = gf_read2dmat(mode,region,side,en,nk,jspin                  &
     &           ,lapw_gf,mat(1:nv2,1:nv2),1,1)                            
         l = gf_read2dmat(mode,region,side,en,nk,jspin                  &
     &           ,lapw_gf,mat(nv2+1:2*nv2,nv2+1:2*nv2),2,2)                
         l = gf_read2dmat(mode,region,side,en,nk,jspin                  &
     &           ,lapw_gf,mat(1:nv2,nv2+1:2*nv2),2,1)                      
         l = gf_read2dmat(mode,region,side,en,nk,jspin,                 &
     &        lapw_gf,mat(nv2+1:2*nv2,1:nv2),1,2)                          
         read2dmat = l 
         CALL timestop("gf_read2dmat")
         RETURN 
      ENDIF 
      IF (PRESENT(pos1)) THEN 
         start1 = (pos1-1)*iofile(fileid)%matsize/iofile(fileid         &
     &        )%matscale+1                                              
         start2 = (pos2-1)*iofile(fileid)%matsize/iofile(fileid         &
     &        )%matscale+1                                              
      ELSE 
         start1 = 1 
         start2 = 1 
      ENDIF 
      !>                                                                
                                                                        
                                                                        
      !<-- Matrix might be a noco-matrix containing four components chec
      !and call reading subroutine recursive in this case               
      IF (jspin == 1) THEN 
         IF (SIZE(mat,1)>lapw_gf%nv2(1)) THEN 
            !This is a noco calculation                                 
            nv2 = lapw_gf%nv2_tot/2 
            l = gf_read2dmat(mode,region,side,en,nk,1                   &
     &           ,lapw_gf,mat(1:nv2,1:nv2),pos1,pos2)                      
            l = gf_read2dmat(mode,region,side,en,nk,2                   &
     &           ,lapw_gf,mat(nv2+1:2*nv2,nv2+1:2*nv2),pos1,pos2)          
            l = gf_read2dmat(mode,region,side,en,nk,3                   &
     &           ,lapw_gf,mat(1:nv2,nv2+1:2*nv2),pos1,pos2)                
            l = gf_read2dmat(mode,region,side,en,nk,4,lapw_gf              &
     &           ,mat(nv2+1:2*nv2,1:nv2),pos1,pos2)                     
            read2dmat = l 
            CALL priv_rotspin(mode,region,side,nv2,mat)
            CALL timestop("gf_read2dmat")
            RETURN 
         ENDIF 
      ENDIF 
      !>                                                                
      !<-- read from simple direct access file                          
      IF (iofile(fileid)%no_hdf) THEN 
                IF (jspin>2) CALL juDFT_error("HDF-file needed for noco in gf_io2dmat")
         IF (PRESENT(pos1)) THEN 
            start1 = 4 
            start2 = pos1+2*(pos2-1) 
         ELSE 
            start1 = 1 
            start2 = 1 
         ENDIF 
         READ(iofile(fileid)%fid,REC = ((en-1)*start1)+start2) mat 
!         WRITE(*,*) "r",iofile(fileid)%fid,en,MINVAL(ABS(mat))         
!     +        ,MAXVAL(ABS(mat))                                        
         read2dmat = .TRUE.
         CALL timestop("gf_read2dmat")
         RETURN 
      ENDIF 
      !>                                                                
                                                                        
                                                                        
      mat = CMPLX(0.0,0.0) 
                                                                        
      !<-- Now determine the spin to use                                
      spin = jspin 
                                      !Not enough spin data             
      IF (iofile(fileid)%spins<jspin)                                   &
     &     THEN                                                         
      IF (jspin == 2) THEN 
                  !This is a magnetic calculation but only one spin in  
         spin = 1 
                  !file, simply use this data for both spins now        
      ELSE 
                              !set off-diagonal elements to zero if     
         mat = CMPLX(0.0,0.0) 
                              !not present
         read2dmat=.true.
         CALL timestop("gf_read2dmat")
         RETURN 
      ENDIF 
      ENDIF 
      !>                                                                
                                                                        
      !<-- Check size                                                   
      IF (SIZE(mat,2) /= SIZE(mat,1))                                   &
     &     CALL juDFT_error('gf_io2dmat: No square matrix')
      !>                                                                
                                                                        
                                                                        
                                                                        
      !<-- read from hdf-file                                           
      ALLOCATE(A(iofile(fileid)%matsize/iofile(fileid)%matscale         &
     &     ,iofile(fileid)%matsize/iofile(fileid)%matscale,2))          
      CALL io_READ(iofile(fileid)%varid,(/start1,start2,1,spin,iofile(fileid)%emap(en),iofile(fileid      )%kmap(nk)/),&
     &     (/SIZE(a,1),SIZE(a,2),1,1,1,1/),"A",&
     &     A(:,:,1))
      CALL io_READ(iofile(fileid)%varid,(/start1,start2,2,spin,iofile(fileid)%emap(en),iofile(fileid      )%kmap(nk)/),&
     &     (/SIZE(a,1),SIZE(a,2),1,1,1,1/),"A",&
     &     A(:,:,2))
      !>                                                                
                                                                        
                                                                        
      !<-- create the phase matrix for the shift                        
      s = iohandle(priv_getiohandle(IO2D_READ,mode,region,side))%shift 
      allocate(phase(iofile(fileid)%matsize/iofile(fileid)%matscale         &
     &     ,iofile(fileid)%matsize/iofile(fileid)%matscale))
      IF (ANY(s /= 0.0)) THEN 
         DO i = 1,SIZE(phase,1)
            DO ii = 1,SIZE(phase,2)
               phase(i,ii) = EXP(CMPLX(0.0,2.*pi_const)*dot_PRODUCT(s     &
     &              ,lapw_gf%global2Dlist(:,i)-lapw_gf%global2Dlist(:,ii)))   
            ENDDO 
         ENDDO 
      ELSE 
         phase=1.0
      ENDIF 
      !>                                                                
      !<-- now map data to temp array
      mat=0.0
      DO i = 1,min(SIZE(A,1),size(lapw_gf%g2map))
         n1=lapw_gf%g2map(i)
         if (n1==0) cycle
         IF (iofile(fileid)%gmap(n1) == 0.OR.iofile(fileid)%gmap(n1       &
     &        )>SIZE(A,1)) cycle
         DO ii = 1,min(SIZE(A,2),size(lapw_gf%g2map))
            n2=lapw_gf%g2map(ii)
            if (n2==0) cycle
            IF (iofile(fileid)%gmap(n2) == 0.OR.iofile(fileid)%gmap(n2  &
     &           )>SIZE(A,2)) cycle
            mat(i,ii)=phase(iofile(fileid)%gmap(n1),iofile(fileid)%gmap(n2))* &
     &       cmplx(A(iofile(fileid)%gmap(n1),iofile(fileid)%gmap(n2),1), &
     &        A(iofile(fileid)%gmap(n1),iofile(fileid)%gmap(n2),2))
!            mat(iofile(fileid)%gmap(n1),iofile(fileid)%gmap(n2)) =       &
!     &           mat(iofile(fileid)%gmap(n1),iofile(fileid)%gmap(n2)     &
!     &           )*CMPLX(A(i,ii,1),A(i,ii,2))
         ENDDO 
      ENDDO 
                                                                        
      DEALLOCATE(A,phase)
      !>                                                                
                                                                        
      read2dmat = .TRUE. 
      CALL timestop("gf_read2dmat")
      END FUNCTION 
                                                                        
      !>                                                                
                                                                        
      !<-- S:priv_rotspin(mode,region,side,nv2,mat)                     
      SUBROUTINE priv_rotspin(mode,region,side,nv2,mat) 
!-----------------------------------------------                        
!      Rotates the spin quantisation axis                               
!           (last modified: 05-12-09) D. Wortmann                       
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER,INTENT(IN)     :: region,side,nv2,mode 
      COMPLEX,INTENT(INOUT)  :: mat(:,:) 
      !>                                                                
      !<-- Locals                                                       
      COMPLEX    :: ci,u(2,2),ut(2,2),tmp(2,2) 
      REAL       :: phi,theta 
      INTEGER    :: n,nn,i(2),ii(2),ioh 
      ci = CMPLX(0.0,1.0) 
      !>                                                                
      ioh = priv_getiohandle(IO2D_READ,mode,region,side) 
      IF (ioh<1) RETURN 
      IF (.NOT.iohandle(ioh)%spinrot) RETURN 
      theta=iohandle(ioh)%theta 
      phi=iohandle(ioh)%phi 
      !<-- Now rotate the spin-quantisation axis if needed              
      u(1,1) =  EXP(-ci*phi/2)*COS(theta/2) 
      u(1,2) = -EXP(-ci*phi/2)*SIN(theta/2) 
      u(2,1) =  EXP(ci*phi/2)*SIN(theta/2) 
      u(2,2) =  EXP(ci*phi/2)*COS(theta/2) 
      ut     = TRANSPOSE(CONJG(u)) 
      DO n = 1,nv2 
         i = (/n,n+nv2/) 
         DO nn = 1,nv2 
            ii = (/nn,nn+nv2/) 
            tmp = mat(i,ii) 
            tmp = MATMUL(MATMUL(u,tmp),ut) 
            mat(i,ii) = tmp 
         ENDDO 
      ENDDO 
      !>                                                                
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
                                                                        
      !<-- F: priv_getiohandle(rw,mode,region,side)                     
      FUNCTION priv_getiohandle(rw,mode,region,side) 
!-----------------------------------------------                        
!                                                                       
!             (last modified: 2004-00-00) D. Wortmann                   
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)     :: rw,mode,region,side 
      INTEGER                :: priv_getiohandle 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: n 
      !>                                                                
      priv_getiohandle = -1 
      IF (.NOT.ALLOCATED(iohandle)) RETURN 
      DO n = 1,SIZE(iohandle) 
         !region must match exactly: region 0 is the plane below layer 1
         IF (iohandle(n)%mode == mode.AND.                              &
     &        iohandle(n)%region==region                                &
     &       .AND.                                                      &
     &        (iohandle(n)%side==side.OR.iohandle(n)%side==0).AND.      &
     &        (iohandle(n)%rw==rw.OR.iohandle(n)%rw==IO2D_READWRITE))   &
     &        THEN                                                      
            priv_getiohandle = n 
            RETURN 
         ENDIF 
      ENDDO 
                                                                        
      END FUNCTION 
      !>                                                                
                                                                        
      !<-- F: priv_getFilehandle(rw,mode,region,side)                   
      FUNCTION priv_getFilehandle(rw,mode,region,side) 
!-----------------------------------------------                        
!                                                                       
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)     :: rw,mode,region,side 
      INTEGER                :: priv_getfilehandle 
      !>                                                                
      priv_getFilehandle = priv_getiohandle(rw,mode,region,side) 
      IF (priv_getFilehandle<1) THEN 
           !write(*,*) "No handle:",rw,mode,region,side                 
           RETURN 
      ENDIF 
      IF (.NOT.ALLOCATED(iofile)) CALL juDFT_error("internal error in gf_io2dmat")
      priv_getFilehandle = iohandle(priv_getFilehandle)%fileno 
                                                                        
      END FUNCTION 
      !>                                                                
                                                                        
      !<-- F:gf_io2dmatFID(mode,region,side,readwrite)                  
      FUNCTION gf_io2dmatFID(mode,region,side,readwrite) 
!******************************************                             
!     FUNCTION to deliver the file-handle for USE in outside modules    
!     returns a -1 IF no filehandle is available                        
!     D. Wortmann                                                       
!******************************************                             
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN) :: mode,region,side,readwrite 
      INTEGER(HID_T)     :: gf_io2dmatFID 
      !>                                                                
      INTEGER :: IOH 
      IOH = priv_getfileHANDLE(readwrite,mode,region,side) 
      IF (IOH<0) THEN 
         gf_io2dmatFID =-1 
      ELSE 
         gf_io2dmatFID = iofile(ioh)%fid 
      ENDIF 
                                                                        
      END FUNCTION 
                                                                        
      !>                                                                
                                                                        
                                                                        
      !<-- F:gf_io2dmatshift(mode,region,side,readwrite)                
      FUNCTION gf_io2dmatshift(mode,region,side,readwrite) 
!******************************************                             
!     FUNCTION to deliver the shift for USE in outside modules          
!     returns a (0,0) IF no iohandle is available                       
!     D. Wortmann                                                       
!******************************************                             
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN) :: mode,region,side,readwrite 
      REAL               :: gf_io2dmatshift(2) 
      !>                                                                
      INTEGER :: IOH 
      IOH = priv_getioHANDLE(readwrite,mode,region,side) 
      IF (IOH<0) THEN 
         gf_io2dmatshift = (/0.0,0.0/) 
      ELSE 
         gf_io2dmatshift = iohandle(ioh)%shift 
      ENDIF 
                                                                        
      END FUNCTION 
                                                                        
      !>                                                                
                                                                        
      !<-- S: gf_io2dmat_create
      SUBROUTINE gf_io2dmat_create(filename,lapw_gf,kpts,sym,spins,nk,matsize,fmeta,fid,varid,l_serial)
!-----------------------------------------------
!  Creates an HDF5 file of 2D matrices: dataset 2Dmat(matsize,matsize,2,spins,nE,nk),
!  energies, kpts (with symmetry and cell attributes) and gvecs (global 2D list)
!-----------------------------------------------
      IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(IN)  :: filename
      TYPE(t_lapw_gf),INTENT(IN)   :: lapw_gf
      TYPE(t_kpts),INTENT(IN)      :: kpts
      TYPE(t_sym),INTENT(IN)       :: sym
      INTEGER,INTENT(IN)           :: spins,nk,matsize
      TYPE(t_io2dmeta),INTENT(IN)  :: fmeta
      INTEGER(HID_T),INTENT(OUT)   :: fid,varid
      !file used by this PE only: no MPI-IO access properties
      LOGICAL,INTENT(IN),OPTIONAL  :: l_serial

      INTEGER(HID_T) :: did
      INTEGER        :: hdferr
      LOGICAL        :: l_ser

      l_ser = .FALSE.
      IF (PRESENT(l_serial)) l_ser = l_serial
      CALL timestart("open file")
      IF (l_ser) THEN
         CALL h5fcreate_f(filename//'.hdf',H5F_ACC_TRUNC_F,fid,hdferr)
      ELSE
         CALL h5fcreate_f(filename//'.hdf',H5F_ACC_TRUNC_F,fid,hdferr,H5P_DEFAULT_f, &
                          hdf_access_prp(filename))
      ENDIF
      IF (hdferr/=0) CALL juDFT_error("Could not create "//filename//".hdf",calledby="gf_io2dmat_create")
      CALL timestop("open file")
      CALL io_createvar(fid,"2Dmat",H5T_NATIVE_DOUBLE,                  &
     &     (/matsize,matsize,2,spins,gf_NoEN(),nk/),varid,chunk=4,fill=.false.)
      CALL io_createvar(fid,"energies",H5T_NATIVE_DOUBLE,(/2,gf_NoEn()/),did)
      CALL io_WRITE(did,(/-1,1/),(/1,gf_NoEn()/),"gf_AllZ",gf_AllZ(0))
      CALL io_dclose(did,hdferr)
      CALL io_createvar(fid,"kpts",H5T_NATIVE_DOUBLE,(/2,kpts%nkpt/),did)
      CALL io_WRITE(did,(/1,1/),(/2,kpts%nkpt/),"kpts_bk",kpts%bk(1:2,:))
      CALL priv_write_cellmeta(did,sym,kpts,fmeta)
      CALL io_dclose(did,hdferr)
      CALL io_createvar(fid,"gvecs",H5T_NATIVE_INTEGER,(/2,SIZE(lapw_gf%global2Dlist,2)/),did)
      CALL io_WRITE(did,(/1,1/),(/2,SIZE(lapw_gf%global2Dlist,2)/),"lapw_global2Dlist",lapw_gf%global2Dlist(:,:))
      CALL io_dclose(did,hdferr)
      END SUBROUTINE gf_io2dmat_create

      !<-- S: priv_write_cellmeta
      SUBROUTINE priv_write_cellmeta(did,sym,kpts,fmeta)
!-----------------------------------------------
!  Symmetry, cell and k-set attributes of the kpts dataset. symop are the
!  in-plane parts of the operations x -> symop x + symtau (fractional).
!-----------------------------------------------
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)   :: did
      TYPE(t_sym),INTENT(IN)      :: sym
      TYPE(t_kpts),INTENT(IN)     :: kpts
      TYPE(t_io2dmeta),INTENT(IN) :: fmeta

      INTEGER :: nmesh(2)
      LOGICAL :: l_full

      CALL io_WRITE_att(did,"nsym",sym%nop2)
      CALL io_WRITE_att(did,"symop",sym%mrot(1:2,1:2,1:sym%nop2))
      CALL io_WRITE_att(did,"symtau",sym%tau(1:2,1:sym%nop2))
      CALL gf_io2dmat_kmesh(kpts%bk(1:2,:),nmesh,l_full)
      CALL io_WRITE_att(did,"nkmesh",nmesh)
      CALL io_WRITE_att(did,"kfull",l_full)
      CALL io_WRITE_att(did,"l_noco",fmeta%l_noco)
      CALL io_WRITE_att(did,"l_soc",fmeta%l_soc)
      !without SOC time reversal holds per spin; with SOC only if the k-set was
      !reduced with it (gfprep soc_timerev, non-magnetic systems)
      CALL io_WRITE_att(did,"l_timerev",fmeta%l_timerev.OR.priv_trneeded(kpts%bk(1:2,:),sym))
      IF (.NOT.fmeta%l_set) RETURN
      CALL io_WRITE_att(did,"amat2d",fmeta%amat2d)
      CALL io_WRITE_att(did,"rkmax",fmeta%rkmax)
      IF (fmeta%zplane(1)>=0.0) CALL io_WRITE_att(did,"zplane",fmeta%zplane(1))
      IF (fmeta%l_vref) CALL io_WRITE_att(did,"vref",fmeta%vref)
      END SUBROUTINE priv_write_cellmeta

      !<-- F: gf_io2dmat_userwrite
      LOGICAL FUNCTION gf_io2dmat_userwrite(mode,region,side)
!-----------------------------------------------
!  is there a write handle of gf_io for this matrix (not one of the
!  automatic handles of the temporary files)
!-----------------------------------------------
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: mode,region,side
      gf_io2dmat_userwrite = .FALSE.
      IF (.NOT.ALLOCATED(iohandle)) RETURN
      gf_io2dmat_userwrite = ANY(iohandle(:nuser_handles)%mode==mode.AND.iohandle(:nuser_handles)%region==region &
     &     .AND.iohandle(:nuser_handles)%side==side.AND.iohandle(:nuser_handles)%rw/=IO2D_READ)
      END FUNCTION gf_io2dmat_userwrite

      !<-- F: priv_trneeded
      LOGICAL FUNCTION priv_trneeded(bk,sym)
!-----------------------------------------------
!  T if the k-points cover their mesh only with the in-plane operations
!  together with time reversal (k -> -k)
!-----------------------------------------------
      IMPLICIT NONE
      REAL,INTENT(IN)        :: bk(:,:)
      TYPE(t_sym),INTENT(IN) :: sym
      INTEGER, ALLOCATABLE :: cover(:,:)
      INTEGER :: nmesh(2),ik,iop,it,m(2)
      LOGICAL :: l_full,l_ops
      priv_trneeded = .FALSE.
      CALL gf_io2dmat_kmesh(bk,nmesh,l_full)
      IF (l_full.OR.ANY(nmesh<1)) RETURN
      ALLOCATE(cover(0:nmesh(1)-1,0:nmesh(2)-1))
      DO it = 1,2
         cover = 0
         DO ik = 1,SIZE(bk,2)
            DO iop = 1,MAX(sym%nop2,1)
               m = MODULO(NINT(MATMUL(bk(:,ik),REAL(sym%mrot(1:2,1:2,iop)))*nmesh),nmesh)
               cover(m(1),m(2)) = 1
               m = MODULO(-m,nmesh)
               IF (it==2) cover(m(1),m(2)) = 1
            ENDDO
         ENDDO
         IF (it==1) l_ops = ALL(cover==1)
      ENDDO
      priv_trneeded = (.NOT.l_ops).AND.ALL(cover==1)
      END FUNCTION priv_trneeded

      !<-- S: gf_io2dmat_kmesh
      SUBROUTINE gf_io2dmat_kmesh(bk,nmesh,l_full)
!-----------------------------------------------
!  Smallest Gamma-centered mesh N1xN2 containing all k-points (nmesh=0 if
!  none up to nmax); l_full: the list is the complete mesh
!-----------------------------------------------
      IMPLICIT NONE
      REAL,INTENT(IN)     :: bk(:,:)
      INTEGER,INTENT(OUT) :: nmesh(2)
      LOGICAL,INTENT(OUT) :: l_full
      INTEGER,PARAMETER   :: nmax = 1024
      REAL,PARAMETER      :: eps = 1.0E-6
      INTEGER :: i,n
      nmesh = 0
      DO i = 1,2
         DO n = 1,nmax
            IF (ALL(ABS(bk(i,:)*n-NINT(bk(i,:)*n))<eps*n)) THEN
               nmesh(i) = n
               EXIT
            ENDIF
         ENDDO
      ENDDO
      l_full = ALL(nmesh>0)
      IF (l_full) l_full = (SIZE(bk,2)==nmesh(1)*nmesh(2))
      END SUBROUTINE gf_io2dmat_kmesh
      !>
      !<-- S: gf_io2dstatus(rw,mode,region,side)                        
      SUBROUTINE gf_io2dstatus(rw,mode,region,side) 
!-----------------------------------------------                        
!                                                                       
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)     :: rw,mode,region,side 
      !>                                                                
      !<--Locals                                                        
      INTEGER             :: n 
                                                                        
      !>                                                                
      WRITE(*,"(a,4(i3))") "Status for requested handle",rw,mode,region &
     &     ,side                                                        
      IF (.NOT.ALLOCATED(iohandle)) THEN 
         WRITE(*,*) "No iohandles defined at all" 
         RETURN 
      ENDIF 
      WRITE(*,*) "List of handles" 
      DO n = 1,SIZE(iohandle) 
         WRITE(*,"(a,i3,a,4(i3))") "No:",n," -- ",iohandle(n)%rw        &
     &        ,iohandle(n)%mode,iohandle(n)%region,iohandle(n)%side     
         IF (iohandle(n)%mode == mode.AND.                              &
     &        iohandle(n)%region == region                              &
     &        .AND.                                                     &
     &        (iohandle(n)%side == side.OR.iohandle(n)%side == 0).AND.  &
     &        (iohandle(n)%rw == rw.OR.iohandle(n)%rw==IO2D_READWRITE)) &
     &        THEN                                                      
            WRITE(*,*) "Match" 
         ENDIF 
      ENDDO 
                                                                        
      END SUBROUTINE 
      !>                                                                
#ifdef CPP_MPI
      LOGICAL FUNCTION priv_create_MPI_memory(mpi_d,layer,kpts,lapw_gf,jspins,l_noco) result(ok)
      !use gmpi
      IMPLICIT NONE
      TYPE(t_gfmpi),INTENT(IN)     :: mpi_d
      TYPE(t_layers),INTENT(IN)  :: layer
      TYPE(t_lapw_gf),INTENT(IN)    :: lapw_gf
      TYPE(t_kpts),INTENT(IN)    :: kpts
      INTEGER,INTENT(IN)         :: jspins
      LOGICAL,INTENT(IN)         :: l_noco
      INTEGER :: n,l,irank,nk,nl,rank_index,e,cmplx_size,old_errhandler
      LOGICAL :: win_ok(2),win_ok_all(2),size_ok
      INTEGER(KIND=MPI_ADDRESS_KIND) :: size_tmat,size_sigma
      INTEGER :: k_list(kpts%nkpt),l_list(layer%num_layers)
      INTEGER :: stat(MPI_STATUS_SIZE)
      TYPE(c_ptr)::ptr,ptr1
      ok=.false.
      mem%matsize=SIZE(lapw_gf%global2Dlist,2)
      mem%spins=jspins
      IF (l_noco)  mem%spins=4
      mem%l_noco=l_noco
      !Determine mapping of k-points and layers
      ALLOCATE(mem%mem_map(kpts%nkpt,layer%num_layers))
      mem%mem_map=-1
      IF (mpi_d%fmpi%irank>0) THEN
          CALL MPI_SEND(mpi_d%k_kperpe,1,MPI_INTEGER,0,1,mpi_d%fmpi%mpi_comm,e)
          CALL MPI_SEND(mpi_d%kl_layerperpe,1,MPI_INTEGER,0,2,mpi_d%fmpi%mpi_comm,e)
          CALL MPI_SEND(mpi_d%k_kpts,mpi_d%k_kperpe,MPI_INTEGER,0,3,mpi_d%fmpi%mpi_comm,e)
          CALL MPI_SEND(mpi_d%kl_layers,mpi_d%kl_layerperpe,MPI_INTEGER,0,4,mpi_d%fmpi%mpi_comm,e)
       ELSE
          nk=mpi_d%k_kperpe
          nl=mpi_d%kl_layerperpe
          DO n=1,nk
               DO l=1,nl
                 mem%mem_map(mpi_d%k_kpts(n),mpi_d%kl_layers(l))=0
               ENDDO
          ENDDO
          DO irank=0,mpi_d%fmpi%isize-1
             If (irank==0) THEN
                k_list(:nk)=mpi_d%k_kpts
                l_list(:nl)=mpi_d%kl_layers
             ELSE
                CALL MPI_RECV(nk,1,MPI_INTEGER,irank,1,mpi_d%fmpi%mpi_comm,stat,e)
                CALL MPI_RECV(nl,1,MPI_INTEGER,irank,2,mpi_d%fmpi%mpi_comm,stat,e)
                CALL MPI_RECV(k_list,nk,MPI_INTEGER,irank,3,mpi_d%fmpi%mpi_comm,stat,e)
                CALL MPI_RECV(l_list,nl,MPI_INTEGER,irank,4,mpi_d%fmpi%mpi_comm,stat,e)
             ENDIF
             rank_index=irank
             DO n=1,nk
               DO l=1,nl
                 mem%mem_map(k_list(n),l_list(l))=rank_index
                 rank_index=rank_index+mpi_d%fmpi%isize
               ENDDO
             ENDDO
           ENDDO
       ENDIF

       IF (mpi_d%fmpi%irank==0) THEN
          write(*,*) "Memory map"
          write(*,*) "nk  layer  map"
          DO n=1,kpts%nkpt
            do l=1,layer%num_layers
               write(*,*) n,l,mem%mem_map(n,l)
            enddo
          enddo
       endif

       CALL MPI_BCAST(mem%mem_map,size(mem%mem_map),MPI_INTEGER,0,mpi_d%fmpi%mpi_comm,e)

       !Determine the amount of local datasets
       n=count(mod(mem%mem_map,mpi_d%fmpi%isize)==mpi_d%fmpi%irank)
       call MPI_Type_size(MPI_DOUBLE_COMPLEX, cmplx_size, e)
       !Sizes are computed in 64 bit: buffer sizes and displacement units (in bytes) are
       !default integers and must not overflow
       size_tmat=2_MPI_ADDRESS_KIND*mem%matsize*mem%matsize*mem%spins*gf_noEN()*2*n
       size_sigma=INT(mem%matsize,MPI_ADDRESS_KIND)*mem%matsize*mem%spins*gf_noEN()*2*n
       size_ok=MAX(size_tmat,size_sigma,4_MPI_ADDRESS_KIND*mem%matsize*mem%matsize*cmplx_size)<=HUGE(0)
       !all PEs must agree, the windows are created collectively
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,size_ok,1,MPI_LOGICAL,MPI_LAND,mpi_d%fmpi%mpi_comm,e)
       IF (.NOT.size_ok) THEN
          IF (mpi_d%fmpi%irank==0) WRITE(0,*) "MPI windows too large, embedding potential and t-matrix are stored in files"
          mem%in_mem=.false.
          RETURN
       END IF
       mem%buffersize_tmat=INT(size_tmat)
       mem%buffersize_sigma=INT(size_sigma)
       mem%stride_sigma=mem%matsize**2
       mem%stride_tmat=4*mem%matsize**2
       !write Info
       write(0,"(a,i5,a,f10.3,a,i8,a,i8)") "Rank:",mpi_d%fmpi%irank," created a window of ",REAL(size_sigma)*cmplx_size/1048576.,"Mbyte for sigma ->",&
                           mem%buffersize_sigma/mem%stride_sigma," Slots of ",mem%stride_sigma*cmplx_size
       write(0,"(a,i5,a,f10.3,a,i8,a,i8)") "Rank:",mpi_d%fmpi%irank," created a window of ",REAL(size_tmat)*cmplx_size/1048576.,"Mbyte for tmat ->",&
                           mem%buffersize_tmat/mem%stride_tmat," Slots of ",mem%stride_tmat*cmplx_size


       !CALL MPI_ALLOC_MEM(int(mem%buffersize_sigma*cmplx_size,MPI_ADDRESS_KIND),MPI_INFO_NULL,ptr,e)
       !CALL c_f_pointer(ptr,mem%buffer_sigma,(/mem%buffersize_sigma/))
       ALLOCATE(mem%buffer_sigma(mem%buffersize_sigma),stat=e)
       IF (e.ne.0) CALL juDFT_error("Could not allocated MPI-Data (sigma) in gf_io2dmat")
       !CALL MPI_ALLOC_MEM(int(mem%buffersize_tmat*cmplx_size,MPI_ADDRESS_KIND),MPI_INFO_NULL,ptr1,e)
       !CALL c_f_pointer(ptr1,mem%buffer_tmat,(/mem%buffersize_tmat/))
       ALLOCATE(mem%buffer_tmat(mem%buffersize_tmat),stat=e)
       IF (e.ne.0) CALL juDFT_error("Could not allocated MPI-Data (tmat) in gf_io2dmat")
       !create windows
       !Not every MPI provides one-sided communication for every setup, e.g. OpenMPI
       !with a single process has no osc component for MPI_WIN_CREATE. Errors are
       !therefore returned instead of aborting, and the caller falls back to files.
       CALL MPI_COMM_GET_ERRHANDLER(mpi_d%fmpi%mpi_comm,old_errhandler,e)
       CALL MPI_COMM_SET_ERRHANDLER(mpi_d%fmpi%mpi_comm,MPI_ERRORS_RETURN,e)
       CALL MPI_WIN_CREATE(mem%buffer_sigma, size_sigma*cmplx_size,mem%stride_sigma*cmplx_size,  &
               Mpi_INFO_NULL, mpi_d%fmpi%mpi_comm, mem%mem_handle_sigma, e)
       win_ok(1)=(e==MPI_SUCCESS)
       CALL MPI_WIN_CREATE(mem%buffer_tmat, size_tmat*cmplx_size,mem%stride_tmat*cmplx_size,  &
               MPI_INFO_NULL, mpi_d%fmpi%mpi_comm, mem%mem_handle_tmat, e)
       win_ok(2)=(e==MPI_SUCCESS)
       CALL MPI_COMM_SET_ERRHANDLER(mpi_d%fmpi%mpi_comm,old_errhandler,e)
       !all PEs must agree, the windows are used collectively
       CALL MPI_ALLREDUCE(win_ok,win_ok_all,2,MPI_LOGICAL,MPI_LAND,mpi_d%fmpi%mpi_comm,e)
       IF (.NOT.ALL(win_ok_all)) THEN
          !Free windows that exist on all PEs (collective). A window that exists on some
          !PEs only cannot be freed, its buffer is kept on those PEs.
          IF (win_ok_all(1)) CALL MPI_WIN_FREE(mem%mem_handle_sigma,e)
          IF (win_ok_all(2)) CALL MPI_WIN_FREE(mem%mem_handle_tmat,e)
          IF (win_ok_all(1).OR..NOT.win_ok(1)) DEALLOCATE(mem%buffer_sigma)
          IF (win_ok_all(2).OR..NOT.win_ok(2)) DEALLOCATE(mem%buffer_tmat)
          IF (mpi_d%fmpi%irank==0) WRITE(0,*) "MPI_WIN_CREATE failed, embedding potential and t-matrix are stored in files"
          mem%in_mem=.false.
          RETURN
       END IF


       CALL  gf_io2dmat_memorysync()
       ok=.true.
      END FUNCTION
#endif
      SUBROUTINE gf_io2dmat_memorysync()
#ifdef CPP_MPI
      !USE MPI
      IMPLICIT NONE
      INTEGER::e,irank
      IF (.not.mem%in_mem) RETURN
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,e)
      if (any(mem%buffer_tmat(::mem%stride_tmat)==cmplx(0.,0.))) write(*,*) "hmm :",irank
      !CALL MPI_WIN_FENCE(0,mem%mem_handle_sigma,e)
      !CALL MPI_WIN_FENCE(0,mem%mem_handle_tmat,e)
#endif
      END SUBROUTINE
#ifdef CPP_MPI
      RECURSIVE SUBROUTINE mpi_memory_read(mode,layer,side,en,nk,jspin,data,rec)
      !USE MPI
      IMPLICIT NONE
      INTEGER,INTENT(IN)      ::mode,layer,side,en,nk,jspin
      COMPLEX,INTENT(OUT)     ::data(:,:)
      LOGICAL,INTENT(IN),OPTIONAL::rec

      INTEGER:: pe,s,isize,e,handle,slots
      INTEGER:: mype
      COMPLEX,POINTER::localdata(:)
      INTEGER        ::localstride
      COMPLEX,ALLOCATABLE:: tmp(:)
      if (layer>size(mem%mem_map,2)) return !Do not read in vacuum
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,mype,e)

      !in noco case call subroutine recursively for all spins
      IF (mem%l_noco.and..not.present(rec)) THEN
         e=size(data,1)/2
         CALL mpi_memory_read(mode,layer,side,en,nk,1,data(:e,:e),.false.)
         CALL mpi_memory_read(mode,layer,side,en,nk,2,data(e+1:,e+1:),.false.)
         CALL mpi_memory_read(mode,layer,side,en,nk,3,data(:e,e+1:),.false.)
         CALL mpi_memory_read(mode,layer,side,en,nk,4,data(e+1:,:e),.false.)
         RETURN
      ENDIF
      if ((mode==IO2d_emb.and.size(data)>mem%matsize**2).or.(mode==IO2d_tmat.and.size(data)>4*mem%matsize**2)) &
            CALL juDFT_error("BUG: size error in mpi_memory_read")
      !Determine PE that holds the data
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,isize,e)
      pe=mod(mem%mem_map(nk,layer),isize)

      !Slot in which we find the data
      s=mem%mem_map(nk,layer)/isize
      IF (mode==IO2D_tmat) THEN
           s=s*mem%spins*gf_noen()  &
             +(en-1)*mem%spins &
             +(jspin-1)
            handle=mem%mem_handle_tmat
            !localdata=>mem%buffer_tmat
            localstride=mem%stride_tmat
            slots=mem%buffersize_tmat/localstride
      ELSE
           s=s*2*mem%spins*gf_noen()  &
             +(side-1)*mem%spins*gf_noen()  &
             +(en-1)*mem%spins &
             +(jspin-1)
           handle=mem%mem_handle_sigma
           !localdata=>mem%buffer_sigma
           localstride=mem%stride_sigma
           slots=mem%buffersize_sigma/localstride
      ENDIF
      !Memory access
      if (pe<0.or.pe>isize-1.or.s<0.or.s>slots) then
           write(0,*) mype,"->",pe
           write(0,*) s, mem%mem_map(nk,layer)
           call juDFT_error("Uups, parallelization broken ")
      endif
      if (size(data)>localstride) call juDFT_error("Uups, too much data")

!      IF (mype==pe) THEN
!           data=reshape(localdata(s*localstride+1:s*localstride+size(data)),shape(data))
!           CALL MPI_WIN_UNLOCK(pe,handle,e)
!      ELSE
        ALLOCATE(tmp(size(data)))
        tmp=0.0001
        CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,pe,0,handle,e)
        CALL MPI_GET(tmp,size(data),MPI_DOUBLE_COMPLEX,pe,&
            int(s,MPI_ADDRESS_KIND),size(data),MPI_DOUBLE_COMPLEX,handle, e)
        CALL MPI_WIN_UNLOCK(pe,handle,e)
        data=reshape(tmp,shape(data))
        DEALLOCATE(tmp)
!      ENDIF

      ! write(0,"('MPI_GET:',i4,'PE',i1,'m,',i2,'l,',i1,'s,',i3,'en,',i3,'nk,',i1,'jspin->',i4,'PE,',i4,'S,',i0,'#,',f15.4,i0)") mype,mode,layer,side,en,nk,jspin,pe,s,size(data),MAXVAL(ABS(data)),handle

      !CALL gf_io2dmat_memorysync()
      END SUBROUTINE

      RECURSIVE SUBROUTINE mpi_memory_write(mode,layer,side,en,nk,jspin,data,rec)
!      USE MPI
      IMPLICIT NONE
      INTEGER,INTENT(IN)      ::mode,layer,side,en,nk,jspin
      COMPLEX,INTENT(IN)      ::data(:,:)
      logical,intent(in),optional::rec

      INTEGER:: pe,s,isize,e,handle,mype,slots
      COMPLEX,POINTER::localdata(:)
      INTEGER        ::localstride
      COMPLEX,ALLOCATABLE:: tmp(:)
      if (layer>size(mem%mem_map,2)) return !Do not write vacuum side
      !in noco case call subroutine recursively for all spins
      IF (mem%l_noco.and..not.present(rec))THEN
         e=size(data,1)/2
         CALL mpi_memory_write(mode,layer,side,en,nk,1,data(:e,:e),.false.)
         CALL mpi_memory_write(mode,layer,side,en,nk,2,data(e+1:,e+1:),.false.)
         CALL mpi_memory_write(mode,layer,side,en,nk,3,data(:e,e+1:),.false.)
         CALL mpi_memory_write(mode,layer,side,en,nk,4,data(e+1:,:e),.false.)
         RETURN
      ENDIF
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,mype,e)
      if ((mode==IO2d_emb.and.size(data)>mem%matsize**2).or.(mode==IO2d_tmat.and.size(data)>4*mem%matsize**2)) THEN
            write(0,*) size(data),mem%matsize
            CALL juDFT_error("BUG: size error in mpi_memory_write")
      endif
      !Determine PE that holds the data
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,isize,e)
      pe=mod(mem%mem_map(nk,layer),isize)

      !Slot in which we find the data
      s=mem%mem_map(nk,layer)/isize
      IF (mode==IO2D_tmat) THEN
           s=s*mem%spins*gf_noen()  &
             +(en-1)*mem%spins &
             +(jspin-1)
            handle=mem%mem_handle_tmat
            !localdata=>mem%buffer_tmat
            localstride=mem%stride_tmat
            slots=mem%buffersize_tmat/localstride
      ELSE
           s=s*2*mem%spins*gf_noen()  &
             +(side-1)*mem%spins*gf_noen()  &
             +(en-1)*mem%spins &
             +(jspin-1)
           handle=mem%mem_handle_sigma
           !localdata=>mem%buffer_sigma
           localstride=mem%stride_sigma
           slots=mem%buffersize_sigma/localstride
      ENDIF
      !Memory access

      if (pe<0.or.pe>isize-1.or.int(s,MPI_ADDRESS_KIND)<0.or.s>slots) then
           write(0,*) "BUG:",mype,"->",pe
           write(0,*) "BUG:",s, mem%mem_map(nk,layer)
           write(0,*) "BUG:",nk,layer
           write(0,*) "BUG:",shape(mem%mem_map)
           write(0,*) "BUG:",mem%mem_map
           call juDFT_error("Uups, parallelization broken ")
      endif
      if (size(data)>localstride) call juDFT_error("Uups, too much data")
      !write(0,"('MPI_PUT:',i4,'PE',i1,'m,',i2,'l,',i1,'s,',i3,'en,',i3,'nk,',i1,'jspin->',i4,'PE,',i4,'S,',i0,'#,',f15.4,i0)") mype,mode,layer,side,en,nk,jspin,pe,s,size(data),MAXVAL(ABS(data)),handle
      CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,pe,0,handle,e)
      ALLOCATE(tmp(size(data)))
      tmp=reshape(data,(/size(data)/))
      CALL MPI_PUT(tmp,size(data),MPI_DOUBLE_COMPLEX,pe,&
            int(s,MPI_ADDRESS_KIND),size(data),MPI_DOUBLE_COMPLEX,handle, e)
      !CALL MPI_ACCUMULATE(tmp,size(data),MPI_DOUBLE_COMPLEX,pe,&
      !      int(s,MPI_ADDRESS_KIND),size(data),MPI_DOUBLE_COMPLEX,MPI_REPLACE,handle, e)
      CALL MPI_WIN_UNLOCK(pe,handle,e)!CALL gf_io2dmat_memorysync()
      DEALLOCATE(tmp)
      !write(*,*) "OK:",jspin,mem%spins
      END SUBROUTINE
#endif
      END                                           
