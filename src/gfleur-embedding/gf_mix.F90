!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_mix
#ifdef CPP_MPI
      USE mpi
#endif
      USE m_constants, ONLY: oUnit
!-------------------------------------------------------------
!Module provides an interface to the mixing part of FLEUR
!as used in the GF-code
!------------------------------------------------------------
      USE m_gf_iodop
      IMPLICIT NONE
      integer,parameter:: nmz=250,nmzxy=100
      PRIVATE
      PUBLIC gf_mix
      CONTAINS

      SUBROUTINE gf_mix(gmpi,embinp,gfmix,layers,ld,jspins)
!*****************************************************************
! Mix the potentials/charge densities of a self-consistent embedding
! calculation. Works on the per-layer t_gflayer containers; the plain
! array views needed by the leaf routines (gf_broyden, gf_apply_metric)
! are gathered from ld(:).
!                          Daniel Wortmann, Mon May 27 11:39:02 2002
!*****************************************************************
      USE m_gf_types
      use m_juDFT
      USE m_gf_broyden
      USE m_gf_potdis
      USE m_gf_metric
      IMPLICIT NONE
      TYPE(t_gfmpi),INTENT(IN)     :: gmpi
      TYPE(t_embinp),INTENT(IN)    :: embinp
      TYPE(t_gfmix),INTENT(IN)     :: gfmix
      TYPE(t_layers),INTENT(IN)    :: layers
      TYPE(t_gflayer),INTENT(IN)   :: ld(:)
      INTEGER,INTENT(IN)           :: jspins

      !<-- locals
      INTEGER               :: fileid
      INTEGER               :: len,layer,mmap,layer_loop,ierr,ds,dl,i,nl
      REAL                  :: volume,distance,c_layer
      INTEGER,DIMENSION(layers%num_layers) :: datastart,dataend
      REAL, ALLOCATABLE     :: sm(:),fsm(:),metric(:),sm_old(:)
      !gathered plain-array views for the leaf routines
      TYPE(t_atoms), ALLOCATABLE  :: atoms(:)
      TYPE(t_stars), ALLOCATABLE  :: stars(:)
      TYPE(t_sphhar),ALLOCATABLE  :: sphhar(:)
      TYPE(t_cell),  ALLOCATABLE  :: cell(:)
      TYPE(t_sym),   ALLOCATABLE  :: sym(:)

      nl = layers%num_layers
#ifdef CPP_MPI
      CALL MPI_BARRIER(gmpi%fmpi%mpi_comm,ierr)
#endif

      IF (embinp%l_potmix) THEN
         fileID = GF_POTFILE
      ELSE
         fileID = GF_CDNFILE
      ENDIF

      IF (.NOT.(gmpi%k_kpts(1) == 1)) RETURN
      IF (.NOT.(gmpi%pe0)) RETURN

      !gather plain arrays
      ALLOCATE(atoms(nl),stars(nl),sphhar(nl),cell(nl),sym(nl))
      DO i=1,nl
         atoms(i)=ld(i)%fi%atoms; stars(i)=ld(i)%stars
         sphhar(i)=ld(i)%sphhar;  cell(i)=ld(i)%fi%cell
         sym(i)=ld(i)%fi%sym
      ENDDO

      !<-- distances
      distance = 0.0;volume = 0.0
      DO layer = 1,nl
         c_layer = MERGE(layers%c(layer),cell(layer)%amat(3,3),layers%c(layer)>0.0)
         CALL gf_potdis(jspins,atoms(layer),stars(layer),sphhar(layer),      &
     &        gmpi%fmpi,sym(layer),cell(layer),c_layer,embinp%l_potmix,       &
     &        layer,distance,volume,ld(layer)%enpara)
      ENDDO
      WRITE(oUnit,*) "Total distance of all layers:",1000*SQRT(distance/volume)

      WRITE(oUnit,*) "Mixing:",gfmix%imix
      WRITE(oUnit,*) "Alpha:",gfmix%alpha
      IF (.NOT.ANY((/0,3,5,7/) == gfmix%imix))                              &
     &     CALL juDFT_error('mix: imix =/= 0,3,5,7 ')

      !<-- length of the mixing vector
      len = 0
      DO layer = 1,nl
         len = len+2*stars(layer)%ng3 +                                     &
     &        atoms(layer)%ntype*(MAXVAL(sphhar(layer)%nlh)+1)              &
     &        *MAXVAL(atoms(layer)%jri)
      ENDDO
      IF (embinp%l_surface) len=len+nmz+2*nmz*stars(1)%ng2-1
      len = len*jspins

      ALLOCATE (sm_old(len),sm(len),fsm(len),metric(len))

      !<-- construct the large mixing vector
      CALL priv_makevector(fileid,layers,gmpi,gfmix,ld,sphhar,atoms,stars,   &
     &     cell,sym,jspins,dataend,datastart,sm,fsm,embinp%l_surface)
      sm_old=sm

      !<-- mixing
      IF (gfmix%imix == 0) THEN
         !simple linear (straight) mixing:  s <- s + alpha * F[s]
         !(the old stmix has been replaced by this inline mixing; the
         ! magnetization spin factor is not applied separately)
         sm = sm + gfmix%alpha*fsm
      ELSE
         CALL gf_broyden("gf_broyd",embinp%l_potmix,embinp%l_surface,gmpi,   &
     &        gfmix%imix,gfmix%maxiter,gfmix%alpha,fsm,stars,atoms,sphhar,   &
     &        cell,sym,jspins,nl,sm)
      END IF

      CALL priv_decomposevector(fileid,layers,gmpi,sphhar,atoms,stars,sym,   &
     &     embinp,jspins,sm,embinp%l_surface)

      !<-- distances via the metric
      metric = gf_apply_metric(embinp%l_surface,embinp%l_potmix,gmpi,stars,  &
     &     atoms,cell,sphhar,sym,jspins,nl,fsm)
      WRITE(oUnit,*) "Distance for all layers:"
      DO layer=1,nl
         WRITE(oUnit,"('Layer:',i5,' distance=',f0.8)") layer,              &
     &        1000.*SQRT(ABS(dot_PRODUCT(metric(datastart(layer):dataend(layer)) &
     &        ,fsm(datastart(layer):dataend(layer))))/cell(layer)%vol)
      ENDDO
      WRITE(oUnit,*) "Mix: total distance =",                               &
     &     1000.*SQRT(ABS(dot_PRODUCT(metric,fsm))/SUM(cell%vol))

      !<-- stepsize
      sm=sm-sm_old
      metric = gf_apply_metric(embinp%l_surface,embinp%l_potmix,gmpi,stars,  &
     &     atoms,cell,sphhar,sym,jspins,nl,sm)
      WRITE(oUnit,*) "Stepsize for all layers:"
      DO layer=1,nl
         WRITE(oUnit,"('Layer:',i5,' step=',f0.8)") layer,                  &
     &        1000.*SQRT(ABS(dot_PRODUCT(metric(datastart(layer):dataend(layer)) &
     &        ,sm(datastart(layer):dataend(layer)))))/cell(layer)%vol
      ENDDO
      WRITE(oUnit,*) "Mix: total stepsize =",                               &
     &     1000.*SQRT(ABS(dot_PRODUCT(metric,sm))/SUM(cell%vol))

      DEALLOCATE (sm,fsm,metric,sm_old)
      END SUBROUTINE


      SUBROUTINE priv_makevector(fileid,layers,gmpi,mix,ld,sphhar,atoms,     &
     &     stars,cell,sym,jspins,dataend,datastart,sout,fsout,l_surface)
!-----------------------------------------------
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      USE m_gf_types
      USE m_gf_fft
      USE m_gf_iodop
      use m_gf_precond
      IMPLICIT NONE
      INTEGER,INTENT(IN)        :: fileid
      TYPE(t_layers),INTENT(IN) :: layers
      TYPE(t_gfmpi),INTENT(IN)  :: gmpi
      TYPE(t_gfmix),INTENT(IN)  :: mix
      TYPE(t_gflayer),INTENT(IN):: ld(:)
      TYPE(t_sphhar),INTENT(IN) :: sphhar(:)
      TYPE(t_atoms),INTENT(IN)  :: atoms(:)
      TYPE(t_stars),INTENT(IN)  :: stars(:)
      TYPE(t_cell),INTENT(IN)   :: cell(:)
      TYPE(t_sym),INTENT(IN)    :: sym(:)
      INTEGER,INTENT(IN)        :: jspins
      INTEGER,INTENT(OUT)       :: dataend(:),datastart(:)
      REAL   ,INTENT(OUT)       :: sout(:),fsout(:)
      logical,intent(in)        :: l_surface

      INTEGER             :: offset,layer,js,na,n,l,i
      REAL   ,ALLOCATABLE :: rho(:,:,:,:,:),vz(:,:,:)
      COMPLEX,ALLOCATABLE :: vpw(:,:,:),vxy(:,:,:,:)
      LOGICAL             :: fixed(layers%num_layers+1)

      !<-- keep some layers fixed?
      fixed = .FALSE.
      OPEN(99,FILE ="fixed_layers")
      DO
         READ(99,*,END = 100,ERR=100) layer
         fixed(layer) = .TRUE.
         WRITE(*,*) "Fixed layer:",layer
      ENDDO
  100 CLOSE(99)

      datastart = 0;dataend = 0;sout = 0.0
      offset = 0
      DO layer=1,layers%num_layers
         datastart(layer) = offset+1
         ALLOCATE(vpw(stars(layer)%ng3,jspins,2))
         ALLOCATE(rho(MAXVAL(atoms(layer)%jri),0:MAXVAL(sphhar(layer)%nlh),  &
     &        atoms(layer)%ntype,jspins,2))
         CALL gf_loddop(fileid,layer,jspins,atoms(layer),stars(layer),       &
     &        sphhar(layer),rho(:,:,:,:,1),vpw(:,:,1),old=.true.)
         CALL gf_loddop(fileid,layer,jspins,atoms(layer),stars(layer),       &
     &        sphhar(layer),rho(:,:,:,:,2),vpw(:,:,2),old=.false.)

         if (fixed(layer)) then
            rho(:,:,:,:,2)=0.0
            vpw(:,:,2)=0.0
         else
            rho(:,:,:,:,2) = rho(:,:,:,:,2)-rho(:,:,:,:,1)
            vpw(:,:,2)   = vpw(:,:,2)-vpw(:,:,1)
         endif

         if (mix%precond>0) CALL gf_precond(layer,ld(layer),gmpi,mix,        &
     &        jspins,vpw(:,:,2),rho(:,:,:,:,2))

         DO js = 1,jspins
            sout(offset+1:offset+stars(layer)%ng3) = REAL(vpw(:stars(layer)%ng3,js,1))
            fsout(offset+1:offset+stars(layer)%ng3) = REAL(vpw(:stars(layer)%ng3,js,2))
            offset = offset+stars(layer)%ng3
            sout(offset+1:offset+stars(layer)%ng3) = AIMAG(vpw(:stars(layer)%ng3,js,1))
            fsout(offset+1:offset+stars(layer)%ng3) = AIMAG(vpw(:stars(layer)%ng3,js,2))
            offset = offset+stars(layer)%ng3
            na = 1
            DO n = 1,atoms(layer)%ntype
               DO l = 0,sphhar(layer)%nlh(sym(layer)%ntypsy(na))
                  DO i = 1,atoms(layer)%jri(n)
                     offset       = offset +1
                     sout(offset) = rho(i,l,n,js,1)
                     fsout(offset) = rho(i,l,n,js,2)
                  END DO
               END DO
               na = na + atoms(layer)%neq(n)
            END DO
         ENDDO
         dataend(layer) = offset
         DEALLOCATE(vpw,rho)
      ENDDO

      if (l_surface) then
         allocate(vz(nmz,jspins,2),vxy(nmzxy,stars(1)%ng2-1,jspins,2))
         call gf_iodop_readvacuum(GF_CDNFILE,vxy(:,:,:,1),vz(:,:,1),old=.true.)
         call gf_iodop_readvacuum(GF_CDNFILE,vxy(:,:,:,2),vz(:,:,2))
         if (fixed(layers%num_layers+1)) then
            vz(:,:,2)=0.0
            vxy(:,:,:,2)=0.0
         else
            vxy(:,:,:,2) = vxy(:,:,:,2)-vxy(:,:,:,1)
            vz(:,:,2)   = vz(:,:,2)-vz(:,:,1)
         endif
         DO js = 1,jspins
            sout(offset+1:offset+nmz) = vz(:,js,1)
            fsout(offset+1:offset+nmz) = vz(:,js,2)
            offset=offset+nmz
            DO n=1,stars(1)%ng2-1
                sout(offset+1:offset+nmzxy) = real(vxy(:,n,js,1))
                sout(offset+nmzxy+1:offset+2*nmzxy)= aimag(vxy(:,n,js,1))
                fsout(offset+1:offset+nmzxy) = real(vxy(:,n,js,2))
                fsout(offset+nmzxy+1:offset+2*nmzxy)= aimag(vxy(:,n,js,2))
                offset=offset+2*nmzxy
            enddo
        enddo
      endif
      END SUBROUTINE


      SUBROUTINE priv_decomposevector(fileid,layers,gmpi,sphhar,atoms,       &
     &     stars,sym,gfinp,jspins,smix,l_surface)
!-----------------------------------------------
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      USE m_gf_types
      USE m_gf_fft
      USE m_gf_iodop
      use m_juDFT
      IMPLICIT NONE
      INTEGER,INTENT(IN)        :: fileid
      TYPE(t_layers),INTENT(IN) :: layers
      TYPE(t_gfmpi),INTENT(IN)  :: gmpi
      TYPE(t_sphhar),INTENT(IN) :: sphhar(:)
      TYPE(t_atoms),INTENT(IN)  :: atoms(:)
      TYPE(t_stars),INTENT(IN)  :: stars(:)
      TYPE(t_sym),INTENT(IN)    :: sym(:)
      TYPE(t_embinp),INTENT(IN) :: gfinp
      INTEGER,INTENT(IN)        :: jspins
      REAL   ,INTENT(IN)        :: smix(:)
      logical,intent(in)        :: l_surface

      INTEGER             :: offset,layer,js,na,n,l,i
      REAL   ,ALLOCATABLE :: rho(:,:,:,:),vz(:,:)
      COMPLEX,ALLOCATABLE :: vpw(:,:),vxy(:,:,:)

      offset = 0
      DO layer = 1,layers%num_layers
         ALLOCATE(vpw(stars(layer)%ng3,jspins))
         ALLOCATE(rho(MAXVAL(atoms(layer)%jri),0:MAXVAL(sphhar(layer)%nlh),  &
     &        atoms(layer)%ntype,jspins))
         DO js = 1,jspins
            vpw(:stars(layer)%ng3,js) = CMPLX(smix(offset+1:offset          &
     &           +stars(layer)%ng3),smix(offset+1+stars(layer)%ng3:offset   &
     &           +2*stars(layer)%ng3))
            offset = offset+2*stars(layer)%ng3
            na = 1
            DO n = 1,atoms(layer)%ntype
               DO l = 0,sphhar(layer)%nlh(sym(layer)%ntypsy(na))
                  DO i = 1,atoms(layer)%jri(n)
                     offset = offset +1
                     rho(i,l,n,js) = smix(offset)
                  END DO
               END DO
               na = na + atoms(layer)%neq(n)
            END DO
         ENDDO
         CALL gf_renamepot(fileID,gmpi%iodop_subcom,layer)
         CALL gf_wrtdop(fileid,layer,jspins,gfinp,atoms(layer),stars(layer), &
     &        sphhar(layer),rho,vpw,.FALSE.,gmpi%self_subcom)
         DEALLOCATE(vpw,rho)
      ENDDO
      if (l_surface) then
         allocate(vz(nmz,jspins),vxy(nmzxy,stars(1)%ng2-1,jspins))
         DO js = 1,jspins
            vz(:,js)=smix(offset+1:offset+nmz)
            offset=offset+nmz
            DO n=1,stars(1)%ng2-1
                vxy(:,n,js)=cmplx(smix(offset+1:offset+nmzxy),smix(offset+nmzxy+1:offset+2*nmzxy))
                offset=offset+2*nmzxy
            enddo
        enddo
        call gf_iodop_writevacuum(GF_CDNFILE,vxy(:,:,:),vz(:,:),gmpi%self_subcom)
      endif
      END SUBROUTINE
      END
