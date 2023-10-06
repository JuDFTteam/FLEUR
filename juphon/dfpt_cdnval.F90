!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_dfpt_cdnval

USE m_juDFT
#ifdef CPP_MPI
use mpi
#endif

CONTAINS

SUBROUTINE dfpt_cdnval(eig_id, dfpt_eig_id, fmpi,kpts,jspin,noco,nococonv,input,banddosdummy,cell,atoms,enpara,stars,&
                  vacuum,sphhar,sym,vTot ,cdnvalJob,den,dosdummy,vacdosdummy,&
                  hub1inp, cdnvalJob1, resultsdummy, resultsdummy1, bqpt, iDtype, iDir, denIm, l_real,&
                  qm_eid_id,dfpt_eigm_id,starsmq,resultsdummy1m,cdnvalJob1m)

   USE m_types
   USE m_constants
   USE m_eig66_io
   USE m_abcof
   USE m_pwden
   USE m_vacden
   USE m_cdnmt       ! calculate the density and orbital moments etc.
   USE m_types_dos
   USE m_types_vacdos
#ifdef CPP_MPI
   USE m_mpi_col_den ! collect density data from parallel nodes
#endif
   USE m_dfpt_rhomt
   USE m_dfpt_rhonmt
   USE m_genMTBasis
   USE m_npy

   IMPLICIT NONE

   TYPE(t_mpi),           INTENT(IN)    :: fmpi

   TYPE(t_enpara),        INTENT(IN)    :: enpara
   TYPE(t_banddos),       INTENT(IN)    :: banddosdummy
   TYPE(t_dos),           INTENT(INOUT) :: dosdummy
   TYPE(t_vacdos),        INTENT(INOUT) :: vacdosdummy
   TYPE(t_input),         INTENT(IN)    :: input
   TYPE(t_vacuum),        INTENT(IN)    :: vacuum
   TYPE(t_noco),          INTENT(IN)    :: noco
   TYPE(t_nococonv),      INTENT(IN)    :: nococonv
   TYPE(t_sym),           INTENT(IN)    :: sym
   TYPE(t_stars),         INTENT(IN)    :: stars
   TYPE(t_cell),          INTENT(IN)    :: cell
   TYPE(t_kpts),          INTENT(IN)    :: kpts
   TYPE(t_sphhar),        INTENT(IN)    :: sphhar
   TYPE(t_atoms),         INTENT(IN)    :: atoms
   TYPE(t_hub1inp),       INTENT(IN)    :: hub1inp
   TYPE(t_potden),        INTENT(IN)    :: vTot
   TYPE(t_cdnvalJob),     INTENT(IN)    :: cdnvalJob, cdnvalJob1
   TYPE(t_results),       INTENT(INOUT) :: resultsdummy, resultsdummy1
   TYPE(t_potden),        INTENT(INOUT) :: den, denIm

   ! Scalar Arguments
   INTEGER,               INTENT(IN)    :: eig_id, dfpt_eig_id, jspin, iDtype, iDir
   LOGICAL,               INTENT(IN)    :: l_real

   REAL, INTENT(IN) :: bqpt(3)

   INTEGER, OPTIONAL, INTENT(IN)    :: qm_eid_id, dfpt_eigm_id
   TYPE(t_stars), OPTIONAL, INTENT(IN)         :: starsmq
   TYPE(t_results), OPTIONAL, INTENT(INOUT)    :: resultsdummy1m
   TYPE(t_cdnvalJob), OPTIONAL, INTENT(IN)    :: cdnvalJob1m

   ! Local Scalars
   INTEGER :: ikpt,ikpt_i,jsp_start,jsp_end,ispin,jsp,iType,ikG,iqdir
   INTEGER :: iErr,nbands,noccbd,nbands1,iLo,l,imLo,ikLo,ikGLo,nbands1m
   INTEGER :: skip_t,skip_tt,nbasfcn,nbasfcnq,nbasfcnmq
   REAL    :: gExt(3), q_loop(3), bkpt(3)

   ! Local Arrays
   COMPLEX ::  f_b8_dummy(3, atoms%ntype)
   REAL,ALLOCATABLE :: we(:),eig(:),we1(:),eig1(:),we1m(:),eig1m(:)
   INTEGER,ALLOCATABLE :: ev_list(:)
   REAL,    ALLOCATABLE :: f(:,:,:,:),g(:,:,:,:),flo(:,:,:,:) ! radial functions

   TYPE (t_lapw)              :: lapw, lapwq, lapwmq
   TYPE (t_orb)               :: orbdummy
   TYPE (t_denCoeffs)         :: denCoeffs
   TYPE (t_denCoeffsOffdiag)  :: denCoeffsOffdiag
   TYPE (t_eigVecCoeffs)      :: eigVecCoeffs, eigVecCoeffs1, eigVecCoeffsPref, eigVecCoeffs1m
   TYPE (t_usdus)             :: usdus
   TYPE (t_mat)               :: zMat, zMat1, zMatPref, zMat1m
   TYPE(t_kpts)               :: kpts_mod

   LOGICAL :: l_minusq

   CALL timestart("dfpt_cdnval")

   call timestart("init")

   l_minusq = PRESENT(qm_eid_id)
   IF (noco%l_mperp) THEN
      ! when the off-diag. part of the density matrix, i.e. m_x and
      ! m_y, is calculated inside the muffin-tins (l_mperp = T), cdnval
      ! is called only once. therefore, several spin loops have been
      ! added. if l_mperp = F, these loops run only from jspin - jspin.
      jsp_start = 1
      jsp_end   = 2
   ELSE
      jsp_start = jspin
      jsp_end   = jspin
   END IF

   ALLOCATE (f(atoms%jmtd,2,0:atoms%lmaxd,input%jspins)) ! Deallocation before mpi_col_den
   ALLOCATE (g(atoms%jmtd,2,0:atoms%lmaxd,input%jspins))
   ALLOCATE (flo(atoms%jmtd,2,atoms%nlod,input%jspins))

   ! Initializations
   CALL usdus%init(atoms,input%jspins)
   CALL denCoeffs%init(atoms,sphhar,jsp_start,jsp_end)
   ! The last entry in denCoeffsOffdiag%init is l_fmpl. It is meant as a switch to a plot of the full magnet.
   ! density without the atomic sphere approximation for the magnet. density.
   CALL denCoeffsOffdiag%init(atoms,noco,sphhar,banddosdummy%l_jDOS,any(noco%l_unrestrictMT).OR.noco%l_mperp)
   !CALL orbdummy%init(atoms,noco,jsp_start,jsp_end)

   IF (denCoeffsOffdiag%l_fmpl.AND.(.NOT.noco%l_mperp)) CALL juDFT_error("for fmpl set noco%l_mperp = T!" ,calledby ="cdnval")

   DO iType = 1, atoms%ntype
      DO ispin = 1, input%jspins
         CALL genMTBasis(atoms,enpara,vTot,fmpi,iType,ispin,usdus,f(:,:,0:,ispin),g(:,:,0:,ispin),flo(:,:,:,ispin))
      END DO
   END DO
   DEALLOCATE (f,g,flo)

   skip_tt = dot_product(enpara%skiplo(:atoms%ntype,jspin),atoms%neq(:atoms%ntype))
   IF (noco%l_soc.OR.noco%l_noco) skip_tt = 2 * skip_tt

   jsp = MERGE(1,jspin,noco%l_noco)

   ! TODO: There was the idea that some problems stemmed from k+q>0.5, so we tried implementing
   !       a backfolding option. This turned out to be unnecessary, but I leave it here for possible
   !       future application
   kpts_mod = kpts
   DO ikpt_i = 1, size(cdnvalJob%k_list)
      ikpt=fmpi%k_list(ikpt_i)
      bkpt = kpts%bk(:, ikpt)
      DO iqdir = 1, 3
         !IF (bkpt(iqdir)+bqpt(iqdir)>=0.5) bkpt(iqdir) = bkpt(iqdir) - 1.0
         !IF (bkpt(iqdir)+bqpt(iqdir)<-0.5) bkpt(iqdir) = bkpt(iqdir) + 1.0
         !IF (bkpt(iqdir)+bqpt(iqdir)>=0.5.AND.ABS(bqpt(iqdir))>1e-8) bkpt(iqdir) = bkpt(iqdir) - 1.0
         !IF (bkpt(iqdir)+bqpt(iqdir)<-0.5.AND.ABS(bqpt(iqdir))>1e-8) bkpt(iqdir) = bkpt(iqdir) + 1.0
      END DO
      kpts_mod%bk(:, ikpt) = bkpt
   END DO

   call timestop("init")

   DO ikpt_i = 1,size(cdnvalJob%k_list)
      ikpt=cdnvalJob%k_list(ikpt_i)

      CALL lapw%init(input,noco,nococonv, kpts,atoms,sym,ikpt,cell, fmpi)
      CALL lapwq%init(input,noco,nococonv, kpts_mod,atoms,sym,ikpt,cell, fmpi, bqpt)

      IF (l_minusq) CALL lapwmq%init(input,noco,nococonv, kpts_mod,atoms,sym,ikpt,cell, fmpi, -bqpt)

      skip_t = skip_tt
      ev_list=cdnvaljob%compact_ev_list(ikpt_i,.FALSE.)
      noccbd = SIZE(ev_list)

      we  = cdnvalJob%weights(ev_list,ikpt)
      we1  = cdnvalJob1%weights(ev_list,ikpt)
      eig = resultsdummy%eig(ev_list,ikpt,jsp)
      eig1 = resultsdummy1%eig(ev_list,ikpt,jsp)

      IF (l_minusq) THEN
         we1m = cdnvalJob1m%weights(ev_list,ikpt)
         eig1m = resultsdummy1m%eig(ev_list,ikpt,jsp)
      END IF

      IF (cdnvalJob%l_evp) THEN
         IF (minval(ev_list) > skip_tt) skip_t = 0
         IF (maxval(ev_list) <= skip_tt) skip_t = noccbd
         IF ((minval(ev_list) <= skip_tt).AND.(maxval(ev_list) > skip_tt)) skip_t = mod(skip_tt,noccbd)
      END IF

      nbasfcn = MERGE(lapw%nv(1)+lapw%nv(2)+2*atoms%nlotot,lapw%nv(1)+atoms%nlotot,noco%l_noco)
      nbasfcnq = MERGE(lapwq%nv(1)+lapwq%nv(2)+2*atoms%nlotot,lapwq%nv(1)+atoms%nlotot,noco%l_noco)

      IF (l_minusq) nbasfcnmq = MERGE(lapwmq%nv(1)+lapwmq%nv(2)+2*atoms%nlotot,lapwmq%nv(1)+atoms%nlotot,noco%l_noco)

      CALL zMat%init(l_real,nbasfcn,noccbd)
      CALL zMat1%init(.FALSE.,nbasfcnq,noccbd)
      CALL zMatPref%init(.FALSE.,nbasfcn,noccbd)

      IF (l_minusq) THEN
         CALL zMat1m%init(.FALSE.,nbasfcnmq,noccbd)
      END IF

      CALL read_eig(eig_id,ikpt,jsp,list=ev_list,neig=nbands,zmat=zMat)
      CALL read_eig(dfpt_eig_id,ikpt,jsp,list=ev_list,neig=nbands1,zmat=zMat1)

      IF (l_minusq) CALL read_eig(dfpt_eigm_id,ikpt,jsp,list=ev_list,neig=nbands1m,zmat=zMat1m)

      ! TODO: Implement correct spin logic here! Only collinear operational for now!
      DO ikG = 1, lapw%nv(jsp)
         gExt = MATMUL(lapw%vk(:, ikG, jsp),cell%bmat)
         IF (zMat%l_real) THEN
            zMatPref%data_c(ikG,:) = ImagUnit * gExt(idir) * zMat%data_r(ikG, :)
         ELSE
            zMatPref%data_c(ikG,:) = ImagUnit * gExt(idir) * zMat%data_c(ikG, :)
         END IF
      END DO

      ! TODO: LOs matching coefficients are unperturbed for now, because they derailed
      !       the calculation. Find out why; forces can use the perturbation!
      !DO ikG = lapw%nv(jsp) + 1, lapw%nv(jsp) + atoms%nlo(iDtype)
      !   iLo = ikG-lapw%nv(jsp)
      !   l = atoms%llo(iLo, iDtype)
      !   DO imLo = 1, 2*l+1
      !      ikLo = lapw%kvec(imLo,iLo,iDtype)
      !      ikGLo = lapw%nv(jsp) + lapw%index_lo(iLo,iDtype) + imLo
      !      !gExt = MATMUL(cell%bmat,lapw%vk(:,ikLo, jsp))
      !      gExt = MATMUL(cell%bmat,lapw%bkpt)
      !      IF (zMat%l_real) THEN
      !         zMatPref%data_c(ikGLo,:) = ImagUnit * gExt(idir) * zMat%data_r(ikGLo, :)
      !      ELSE
      !         zMatPref%data_c(ikGLo,:) = ImagUnit * gExt(idir) * zMat%data_c(ikGLo, :)
      !      END IF
      !   END DO
      !END DO

      !IF (.NOT.(nbands==nbands1)) TODO: Can this ever be a problem?
#ifdef CPP_MPI
      CALL MPI_BARRIER(fmpi%mpi_comm,iErr) ! Synchronizes the RMA operations
#endif

      IF (noccbd.LE.0) CYCLE ! Note: This jump has to be after the MPI_BARRIER is called

      ! valence density in the atomic spheres
      CALL eigVecCoeffs%init(input,atoms,jspin,noccbd,noco%l_mperp)
      CALL eigVecCoeffs1%init(input,atoms,jspin,noccbd,noco%l_mperp)
      CALL eigVecCoeffsPref%init(input,atoms,jspin,noccbd,noco%l_mperp)

      IF (l_minusq) CALL eigVecCoeffs1m%init(input,atoms,jspin,noccbd,noco%l_mperp)

      DO ispin = jsp_start, jsp_end
         ! TODO: Does this spin logic hold for noco?
         CALL abcof(input,atoms,sym,cell,lapw,noccbd,usdus,noco,nococonv,ispin,&
                    eigVecCoeffs%abcof(:,0:,0,:,ispin),eigVecCoeffs%abcof(:,0:,1,:,ispin),&
                    eigVecCoeffs%ccof(-atoms%llod:,:,:,:,ispin),zMat)
         CALL abcof(input,atoms,sym,cell,lapwq,noccbd,usdus,noco,nococonv,ispin,&
                    eigVecCoeffs1%abcof(:,0:,0,:,ispin),eigVecCoeffs1%abcof(:,0:,1,:,ispin),&
                    eigVecCoeffs1%ccof(-atoms%llod:,:,:,:,ispin),zMat1)
         CALL abcof(input,atoms,sym,cell,lapw,noccbd,usdus,noco,nococonv,ispin,&
                    eigVecCoeffsPref%abcof(:,0:,0,:,ispin),eigVecCoeffsPref%abcof(:,0:,1,:,ispin),&
                    eigVecCoeffsPref%ccof(-atoms%llod:,:,:,:,ispin),zMatPref)
         IF (l_minusq) CALL abcof(input,atoms,sym,cell,lapwmq,noccbd,usdus,noco,nococonv,ispin,&
                                  eigVecCoeffs1m%abcof(:,0:,0,:,ispin),eigVecCoeffs1m%abcof(:,0:,1,:,ispin),&
                                  eigVecCoeffs1m%ccof(-atoms%llod:,:,:,:,ispin),zMat1)

         eigVecCoeffs1%abcof(:,0:,:,iDtype,ispin) = eigVecCoeffs1%abcof(:,0:,:,iDtype,ispin) + eigVecCoeffsPref%abcof(:,0:,:,iDtype,ispin)
         eigVecCoeffs1%ccof(-atoms%llod:,:,:,iDtype,ispin) = eigVecCoeffs1%ccof(-atoms%llod:,:,:,iDtype,ispin) + eigVecCoeffsPref%ccof(-atoms%llod:,:,:,iDtype,ispin)

         IF (l_minusq) THEN
            eigVecCoeffs1m%abcof(:,0:,:,iDtype,ispin) = eigVecCoeffs1m%abcof(:,0:,:,iDtype,ispin) + eigVecCoeffsPref%abcof(:,0:,:,iDtype,ispin)
            eigVecCoeffs1m%ccof(-atoms%llod:,:,:,iDtype,ispin) = eigVecCoeffs1m%ccof(-atoms%llod:,:,:,iDtype,ispin) + eigVecCoeffsPref%ccof(-atoms%llod:,:,:,iDtype,ispin)
         END IF

         IF (l_minusq) THEN
            CALL dfpt_rhomt(atoms,we,we1,noccbd,ispin,ispin,-bqpt,.TRUE.,eigVecCoeffs,eigVecCoeffs1,denCoeffs,eigVecCoeffs1m)
            CALL dfpt_rhonmt(atoms,sphhar,we,we1,noccbd,ispin,ispin,-bqpt,.TRUE.,.FALSE.,sym,eigVecCoeffs,eigVecCoeffs1,denCoeffs,eigVecCoeffs1m)
            CALL dfpt_rhomtlo(atoms,noccbd,we,we1,ispin,ispin,-bqpt,.TRUE.,eigVecCoeffs,eigVecCoeffs1,denCoeffs,eigVecCoeffs1m)
            CALL dfpt_rhonmtlo(atoms,sphhar,sym,noccbd,we,we1,eigVecCoeffs,eigVecCoeffs1,denCoeffs,ispin,ispin,.TRUE.,-bqpt,eigVecCoeffs1m)
         ELSE
            CALL dfpt_rhomt(atoms,we,we1,noccbd,ispin,ispin,bqpt,.TRUE.,eigVecCoeffs,eigVecCoeffs1,denCoeffs)
            CALL dfpt_rhonmt(atoms,sphhar,we,we1,noccbd,ispin,ispin,bqpt,.TRUE.,.FALSE.,sym,eigVecCoeffs,eigVecCoeffs1,denCoeffs)
            CALL dfpt_rhomtlo(atoms,noccbd,we,we1,ispin,ispin,bqpt,.TRUE.,eigVecCoeffs,eigVecCoeffs1,denCoeffs)
            CALL dfpt_rhonmtlo(atoms,sphhar,sym,noccbd,we,we1,eigVecCoeffs,eigVecCoeffs1,denCoeffs,ispin,ispin,.TRUE.,bqpt)
         END IF
      END DO ! end loop over ispin
      IF (noco%l_mperp) then
         IF (.NOT.l_minusq) THEN
            call timestart("denCoeffsOffdiag%calcCoefficients")
            CALL dfpt_rhomt(atoms,we,we1,noccbd,2,1,bqpt,.TRUE.,eigVecCoeffs,eigVecCoeffs1,denCoeffs,eigVecCoeffs1m)
            CALL dfpt_rhonmt(atoms,sphhar,we,we1,noccbd,2,1,bqpt,.TRUE.,.FALSE.,sym,eigVecCoeffs,eigVecCoeffs1,denCoeffs,eigVecCoeffs1m)
            CALL dfpt_rhomtlo(atoms,noccbd,we,we1,2,1,bqpt,.TRUE.,eigVecCoeffs,eigVecCoeffs1,denCoeffs,eigVecCoeffs1m)
            CALL dfpt_rhonmtlo(atoms,sphhar,sym,noccbd,we,we1,eigVecCoeffs,eigVecCoeffs1,denCoeffs,2,1,.TRUE.,bqpt,eigVecCoeffs1m)
            CALL dfpt_rhomt(atoms,we,we1,noccbd,1,2,bqpt,.TRUE.,eigVecCoeffs,eigVecCoeffs1,denCoeffs,eigVecCoeffs1m)
            CALL dfpt_rhonmt(atoms,sphhar,we,we1,noccbd,1,2,bqpt,.TRUE.,.FALSE.,sym,eigVecCoeffs,eigVecCoeffs1,denCoeffs,eigVecCoeffs1m)
            CALL dfpt_rhomtlo(atoms,noccbd,we,we1,1,2,bqpt,.TRUE.,eigVecCoeffs,eigVecCoeffs1,denCoeffs,eigVecCoeffs1m)
            CALL dfpt_rhonmtlo(atoms,sphhar,sym,noccbd,we,we1,eigVecCoeffs,eigVecCoeffs1,denCoeffs,1,2,.TRUE.,bqpt,eigVecCoeffs1m)
            call timestop("denCoeffsOffdiag%calcCoefficients")
         ELSE
            call timestart("denCoeffsOffdiag%calcCoefficients")
            CALL dfpt_rhomt(atoms,we,we1,noccbd,2,1,bqpt,.TRUE.,eigVecCoeffs,eigVecCoeffs1,denCoeffs)
            CALL dfpt_rhonmt(atoms,sphhar,we,we1,noccbd,2,1,bqpt,.TRUE.,.FALSE.,sym,eigVecCoeffs,eigVecCoeffs1,denCoeffs)
            CALL dfpt_rhomtlo(atoms,noccbd,we,we1,2,1,bqpt,.TRUE.,eigVecCoeffs,eigVecCoeffs1,denCoeffs)
            CALL dfpt_rhonmtlo(atoms,sphhar,sym,noccbd,we,we1,eigVecCoeffs,eigVecCoeffs1,denCoeffs,2,1,.TRUE.,bqpt)
            CALL dfpt_rhomt(atoms,we,we1,noccbd,1,2,bqpt,.TRUE.,eigVecCoeffs,eigVecCoeffs1,denCoeffs)
            CALL dfpt_rhonmt(atoms,sphhar,we,we1,noccbd,1,2,bqpt,.TRUE.,.FALSE.,sym,eigVecCoeffs,eigVecCoeffs1,denCoeffs)
            CALL dfpt_rhomtlo(atoms,noccbd,we,we1,1,2,bqpt,.TRUE.,eigVecCoeffs,eigVecCoeffs1,denCoeffs)
            CALL dfpt_rhonmtlo(atoms,sphhar,sym,noccbd,we,we1,eigVecCoeffs,eigVecCoeffs1,denCoeffs,1,2,.TRUE.,bqpt)
            call timestop("denCoeffsOffdiag%calcCoefficients")
         END IF
      endif

      ! valence density in the interstitial and vacuum region has to be called only once (if jspin=1) in the non-collinear case
      IF (.NOT.((jspin.EQ.2).AND.noco%l_noco)) THEN
         ! valence density in the interstitial region
         IF (.NOT.l_minusq) THEN
            CALL pwden(stars,kpts,banddosdummy ,input,fmpi,noco,nococonv,cell,atoms,sym,ikpt,&
                       jspin,lapw,noccbd,ev_list,we,eig,den,resultsdummy,f_b8_dummy,zMat,dosdummy,bqpt,lapwq,we1,zMat1,iDir)
         ELSE
            CALL pwden(stars,kpts,banddosdummy ,input,fmpi,noco,nococonv,cell,atoms,sym,ikpt,&
                       jspin,lapw,noccbd,ev_list,we,eig,den,resultsdummy,f_b8_dummy,zMat,dosdummy,bqpt,lapwq,we1,zMat1,iDir,lapwmq,zMat1m)
         END IF
         IF (input%film) THEN
            CALL vacden(vacuum,stars,input,cell,atoms,noco,nococonv,banddosdummy,&
                        we,ikpt,jspin,REAL(vTot%vac(:,1,:,:)),noccbd,ev_list,lapw,enpara%evac,den,zMat,vacdosdummy,dosdummy,lapwq,we1,zMat1)
         END IF
      END IF
   END DO ! end of k-point loop

#ifdef CPP_MPI
   DO ispin = jsp_start,jsp_end
      CALL mpi_col_den(fmpi,sphhar,atoms,stars,vacuum,input,noco,ispin,dosdummy,vacdosdummy,&
                       resultsdummy,denCoeffs,orbdummy,denCoeffsOffdiag,den)
   END DO
#endif

   IF (fmpi%irank==0) THEN
      CALL timestart("cdnmt")
      CALL cdnmt(input%jspins,input,atoms,sym,sphhar,noco,jsp_start,jsp_end,enpara,banddosdummy,&
                 vTot%mt(:,0,:,:),denCoeffs,usdus,orbdummy,denCoeffsOffdiag,den%mt,hub1inp,rhoIm=denIm%mt)
      CALL timestop("cdnmt")
   END IF

   CALL timestop("dfpt_cdnval")

END SUBROUTINE dfpt_cdnval

END MODULE m_dfpt_cdnval
