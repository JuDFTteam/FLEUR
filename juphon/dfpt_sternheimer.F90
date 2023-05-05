!--------------------------------------------------------------------------------
! Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_sternheimer
   USE m_types
   USE m_make_stars
   USE m_vgen
   USE m_dfpt_eigen_new
   USE m_dfpt_cdngen
   USE m_dfpt_vgen
   USE m_dfpt_fermie
   USE m_mix
   USE m_constants
   USE m_cdn_io
   USE m_eig66_io
   !USE m_npy

IMPLICIT NONE

CONTAINS
   SUBROUTINE dfpt_sternheimer(fi, xcpot, sphhar, stars, starsq, nococonv, qpts, fmpi, results, resultsq, enpara, hybdat, mpdata, &
                               forcetheo, rho, vTot, grRho, grVtot, grVext, grVC, iQ, iDType, iDir, dfpt_tag, eig_id, &
                               l_real, results1, dfpt_eig_id, q_eig_id, &
                               denIn1, vTot1, denIn1Im, vTot1Im, vC1, vC1Im, &
                               starsmq, resultsmq, dfpt_eigm_id, qm_eig_id, results1m, vTot1m, vTot1mIm)
      TYPE(t_fleurinput), INTENT(IN)    :: fi
      CLASS(t_xcpot),     INTENT(IN)    :: xcpot
      TYPE(t_sphhar),     INTENT(IN)    :: sphhar
      TYPE(t_stars),      INTENT(IN)    :: stars
      TYPE(t_stars),      INTENT(INOUT) :: starsq
      TYPE(t_nococonv),   INTENT(IN)    :: nococonv
      TYPE(t_kpts),       INTENT(IN)    :: qpts
      TYPE(t_mpi),        INTENT(IN)    :: fmpi
      TYPE(t_results),    INTENT(INOUT) :: results, resultsq, results1
      TYPE(t_hybdat),     INTENT(INOUT) :: hybdat
      TYPE(t_mpdata),     INTENT(INOUT) :: mpdata
      CLASS(t_forcetheo), INTENT(INOUT) :: forcetheo
      TYPE(t_enpara),     INTENT(INOUT) :: enpara
      TYPE(t_potden),     INTENT(IN)    :: rho, vTot, grRho, grVtot, grVext, grVC

      TYPE(t_potden), INTENT(INOUT) :: denIn1, vTot1, denIn1Im, vTot1Im, vC1, vC1Im

      TYPE(t_potden) :: denOut1, denOut1Im, vx, rho_loc, rho_loc0

      INTEGER, INTENT(IN) :: iQ, iDtype, iDir, eig_id, q_eig_id
      LOGICAL, INTENT(IN) :: l_real
      CHARACTER(len=20), INTENT(IN) :: dfpt_tag

      INTEGER, INTENT(OUT) :: dfpt_eig_id

      TYPE(t_stars), OPTIONAL, INTENT(INOUT) :: starsmq
      TYPE(t_results), OPTIONAL,   INTENT(INOUT) :: resultsmq, results1m
      INTEGER, OPTIONAL, INTENT(IN) :: qm_eig_id
      INTEGER, OPTIONAL, INTENT(OUT) :: dfpt_eigm_id
      TYPE(t_potden), OPTIONAL, INTENT(INOUT) :: vTot1m, vTot1mIm

#ifdef CPP_MPI
      INTEGER :: ierr
#endif

      INTEGER :: archiveType, iter, killcont(6), iterm
      REAL    :: bqpt(3), bmqpt(3)
      LOGICAL :: l_cont, l_exist, l_lastIter, l_dummy, strho, onedone, final_SH_it, l_exitus, l_minusq, l_existm

      TYPE(t_hub1data) :: hub1data
      TYPE(t_banddos)  :: banddosdummy
      TYPE(t_field)    :: field2
      TYPE(t_potden)   :: denIn1m, denIn1mIm, denOut1m, denOut1mIm

      l_minusq = PRESENT(starsmq)


      ! In this order: V1_pw_pw, T1_pw, S1_pw, V1_MT, ikGH0_MT, ikGS0_MT
      !killcont = [1,0,0,0,0,0]
      killcont = [1,1,1,1,1,1]
      !killcont = [0,0,0,0,0,1]

      CALL rho_loc%copyPotDen(rho)
      CALL rho_loc0%copyPotDen(rho)
      CALL rho_loc0%resetPotDen()
      CALL vx%copyPotDen(vTot)
      ALLOCATE(vx%pw_w, mold=vx%pw)
      vx%pw_w = vTot%pw_w

      banddosdummy = fi%banddos

      CALL make_stars(starsq, fi%sym, fi%atoms, fi%vacuum, sphhar, fi%input, fi%cell, fi%noco, fmpi, qpts%bk(:,iQ), iDtype, iDir)
      starsq%ufft = stars%ufft
      CALL denIn1%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_DEN, l_dfpt=.TRUE.)
      CALL denIn1Im%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_DEN, l_dfpt=.FALSE.)

      IF (l_minusq) THEN
         CALL make_stars(starsmq, fi%sym, fi%atoms, fi%vacuum, sphhar, fi%input, fi%cell, fi%noco, fmpi, -qpts%bk(:,iQ), iDtype, iDir)
         starsmq%ufft = stars%ufft
         CALL denIn1m%init(starsmq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_DEN, l_dfpt=.TRUE.)
         CALL denIn1mIm%init(starsmq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_DEN, l_dfpt=.FALSE.)
         INQUIRE(FILE=TRIM(dfpt_tag)//'m.hdf',EXIST=l_existm)
      END IF

      INQUIRE(FILE=TRIM(dfpt_tag)//'.hdf',EXIST=l_exist)

      archiveType = CDN_ARCHIVE_TYPE_CDN1_const
      IF (ANY(fi%noco%l_unrestrictMT)) THEN
         archiveType = CDN_ARCHIVE_TYPE_FFN_const
      ELSE IF (fi%noco%l_noco) THEN
         archiveType = CDN_ARCHIVE_TYPE_NOCO_const
      END IF

      IF (fmpi%irank == 0) THEN
         strho = .NOT.l_exist
         onedone = .NOT.strho
         final_SH_it = .FALSE.
      END IF

#ifdef CPP_MPI
      CALL MPI_BCAST(strho,1,MPI_LOGICAL,0,fmpi%mpi_comm,ierr)
      CALL MPI_BCAST(onedone,1,MPI_LOGICAL,0,fmpi%mpi_comm,ierr)
      CALL MPI_BCAST(final_SH_it,1,MPI_LOGICAL,0,fmpi%mpi_comm,ierr)
#endif

      iter = 0
      iterm = 0
      l_cont = (iter < fi%input%itmax)

      IF (fmpi%irank==0.AND.l_exist) CALL readDensity(starsq, fi%noco, fi%vacuum, fi%atoms, fi%cell, sphhar, &
                                                      fi%input, fi%sym, archiveType, CDN_INPUT_DEN_const, 0, &
                                                      results%ef, results%last_distance, l_dummy, denIn1,  &
                                                      inFilename=TRIM(dfpt_tag),denIm=denIn1Im)
      IF (fmpi%irank==0.AND.l_exist.AND.l_minusq) CALL readDensity(starsmq, fi%noco, fi%vacuum, fi%atoms, fi%cell, sphhar, &
                                                      fi%input, fi%sym, archiveType, CDN_INPUT_DEN_const, 0, &
                                                      results%ef, results%last_distance, l_dummy, denIn1m,  &
                                                      inFilename=TRIM(dfpt_tag)//'m',denIm=denIn1mIm)

      IF (fmpi%irank==0.AND..NOT.l_exist) denIn1%iter = 1

#ifdef CPP_MPI
      CALL MPI_BCAST(denIn1%iter,1,MPI_INTEGER,0,fmpi%mpi_comm,ierr)
#endif

      CALL vTot1%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
      CALL vTot1Im%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)

      CALL vC1%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
      CALL vC1Im%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)

      dfpt_eig_id = open_eig(fmpi%mpi_comm, lapw_dim_nbasfcn, fi%input%neig, fi%kpts%nkpt, fi%input%jspins, fi%noco%l_noco, &
                             .NOT.fi%INPUT%eig66(1), .FALSE., fi%noco%l_soc, fi%INPUT%eig66(1), .FALSE., fmpi%n_size)

      bqpt = qpts%bk(:, iQ)

      IF (l_minusq) THEN
         CALL vTot1m%init(starsmq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
         CALL vTot1mIm%init(starsmq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)

         dfpt_eigm_id = open_eig(fmpi%mpi_comm, lapw_dim_nbasfcn, fi%input%neig, fi%kpts%nkpt, fi%input%jspins, fi%noco%l_noco, &
                                .NOT.fi%INPUT%eig66(1), .FALSE., fi%noco%l_soc, fi%INPUT%eig66(1), .FALSE., fmpi%n_size)

         bmqpt = -qpts%bk(:, iQ)
      END IF

!#ifdef CPP_CHASE
!      CALL init_chase(fmpi, fi%input, fi%atoms, fi%kpts, fi%noco, l_real)
!#endif

      l_lastIter = .FALSE.
      scfloop: DO WHILE (l_cont)
         !IF (.NOT.strho) iter = iter + 1 !TODO: Eventually this will be right.
         IF (onedone) iter = iter + 1
         l_lastIter = l_lastIter.OR.(iter.EQ.fi%input%itmax)

         CALL timestart("Sternheimer Iteration")
         IF (fmpi%irank==0.AND.onedone) THEN
            WRITE (oUnit, FMT=8100) iter
8100        FORMAT(/, 10x, '   iter=  ', i5)
         END IF !fmpi%irank==0

#ifdef CPP_CHASE
         CALL chase_distance(results%last_distance)
#endif

         CALL denIn1%distribute(fmpi%mpi_comm)
         CALL denIn1Im%distribute(fmpi%mpi_comm)

         IF (l_minusq) THEN
            CALL denIn1m%distribute(fmpi%mpi_comm)
            CALL denIn1mIm%distribute(fmpi%mpi_comm)
         END IF

         CALL timestart("Generation of potential perturbation")
         IF (strho) THEN
            write(oUnit, *) "vExt1", iDir
            CALL dfpt_vgen(hybdat,fi%field,fi%input,xcpot,fi%atoms,sphhar,stars,fi%vacuum,fi%sym,&
                           fi%cell ,fi%sliceplot,fmpi,fi%noco,nococonv,rho_loc0,vTot,&
                           starsq,denIn1Im,vTot1,.FALSE.,vTot1Im,denIn1,iDtype,iDir,[1,1]) ! comparison is [1,0]
            IF (l_minusq) THEN
               CALL dfpt_vgen(hybdat,fi%field,fi%input,xcpot,fi%atoms,sphhar,stars,fi%vacuum,fi%sym,&
                              fi%cell ,fi%sliceplot,fmpi,fi%noco,nococonv,rho_loc0,vTot,&
                              starsmq,denIn1mIm,vTot1m,.FALSE.,vTot1mIm,denIn1m,iDtype,iDir,[1,1]) ! comparison is [1,0]
            END IF
         ELSE
            write(oUnit, *) "vEff1", iDir
            CALL dfpt_vgen(hybdat,fi%field,fi%input,xcpot,fi%atoms,sphhar,stars,fi%vacuum,fi%sym,&
                           fi%cell ,fi%sliceplot,fmpi,fi%noco,nococonv,rho_loc,vTot,&
                           starsq,denIn1Im,vTot1,.TRUE.,vTot1Im,denIn1,iDtype,iDir,[1,1]) ! comparison is [1,0]
            IF (l_minusq) THEN
               CALL dfpt_vgen(hybdat,fi%field,fi%input,xcpot,fi%atoms,sphhar,stars,fi%vacuum,fi%sym,&
                              fi%cell ,fi%sliceplot,fmpi,fi%noco,nococonv,rho_loc,vTot,&
                              starsmq,denIn1mIm,vTot1m,.TRUE.,vTot1mIm,denIn1m,iDtype,iDir,[1,1]) ! comparison is [1,0]
            END IF
         END IF

         IF (final_SH_it) THEN !!!!!!!!!!!!!!!!!!!!!!!!!!!
         !IF (.NOT.strho) THEN
            write(oUnit, *) "vC1", iDir
            CALL dfpt_vgen(hybdat,fi%field,fi%input,xcpot,fi%atoms,sphhar,stars,fi%vacuum,fi%sym,&
                           fi%cell ,fi%sliceplot,fmpi,fi%noco,nococonv,rho_loc,vTot,&
                           starsq,denIn1Im,vC1,.FALSE.,vC1Im,denIn1,iDtype,iDir,[0,0])
           !!!!!!!!!!!!!!!!!!!!!!
           !denIn1%mt(:,0:,iDtype,:) = denIn1%mt(:,0:,iDtype,:) + grRho%mt(:,0:,iDtype,:)
         END IF

         CALL timestop("Generation of potential perturbation")

#ifdef CPP_MPI
         CALL MPI_BARRIER(fmpi%mpi_comm, ierr)
#endif

         IF (strho) THEN
            vTot1%mt(:,0:,iDtype,1) = vTot1%mt(:,0:,iDtype,1) + grVext%mt(:,0:,iDtype,1)
            IF (fi%input%jspins==2) vTot1%mt(:,0:,iDtype,2) = vTot1%mt(:,0:,iDtype,2) + grVext%mt(:,0:,iDtype,1)
            CALL save_npy(TRIM(dfpt_tag)//"vExt1pw.npy",vTot1%pw)
            CALL save_npy(TRIM(dfpt_tag)//"vExt1mtre.npy",vTot1%mt)
            CALL save_npy(TRIM(dfpt_tag)//"vExt1mtim.npy",vTot1Im%mt)
            IF (l_minusq) THEN
               vTot1m%mt(:,0:,iDtype,1) = vTot1m%mt(:,0:,iDtype,1) + grVext%mt(:,0:,iDtype,1)
               IF (fi%input%jspins==2) vTot1m%mt(:,0:,iDtype,2) = vTot1m%mt(:,0:,iDtype,2) + grVext%mt(:,0:,iDtype,1)
            END IF
            !CALL timestop("Sternheimer Iteration")
            !RETURN
         ELSE
            !vC1%mt(:,0:,iDtype,:) = vC1%mt(:,0:,iDtype,:) + grVC%mt(:,0:,iDtype,:)
            !CALL save_npy(TRIM(dfpt_tag)//"vC1pw.npy",vC1%pw(:,1))
            !CALL save_npy(TRIM(dfpt_tag)//"vC1mtre.npy",vC1%mt(:,0:,1,1))
            !CALL save_npy(TRIM(dfpt_tag)//"vC1mtim.npy",vC1Im%mt(:,0:,1,1))
            !vC1%mt(:,0:,iDtype,:) = vC1%mt(:,0:,iDtype,:) - grVC%mt(:,0:,iDtype,:)
            vTot1%mt(:,0:,iDtype,:) = vTot1%mt(:,0:,iDtype,:) + grVtot%mt(:,0:,iDtype,:)
            CALL save_npy(TRIM(dfpt_tag)//"vEff1pw.npy",vTot1%pw(:,1))
            CALL save_npy(TRIM(dfpt_tag)//"vEff1mtre.npy",vTot1%mt(:,0:,1,1))
            CALL save_npy(TRIM(dfpt_tag)//"vEff1mtim.npy",vTot1Im%mt(:,0:,1,1))
            IF (l_minusq) THEN
               vTot1m%mt(:,0:,iDtype,:) = vTot1m%mt(:,0:,iDtype,:) + grVtot%mt(:,0:,iDtype,:)
            END IF
            !CALL timestop("Sternheimer Iteration")
            !RETURN
         END IF

         CALL vTot1%distribute(fmpi%mpi_comm)

         IF (l_minusq) THEN
            CALL vTot1m%distribute(fmpi%mpi_comm)
         END IF

#ifdef CPP_MPI
         CALL MPI_BARRIER(fmpi%mpi_comm, ierr)
#endif

         CALL timestart("dfpt eigen")

         !IF (.NOT. fi%input%eig66(1)) THEN
         !   CALL eigen(fi, fmpi, stars, sphhar, xcpot, forcetheo, enpara, nococonv, mpdata, &
         !              hybdat, iter, eig_id, results, rho, vTot, vx, hub1data, &
         !              bqpt=bqpt, dfpt_eig_id=dfpt_eig_id, iDir=iDir, iDtype=iDtype, &
         !              starsq=starsq, v1real=vTot1, v1imag=vTot1Im, killcont=killcont, l_real=l_real)
         !END IF
         CALL dfpt_eigen_new(fi, sphhar, results, resultsq, results1, fmpi, enpara, nococonv, starsq, vTot1, vTot1Im, &
                             vTot, rho, bqpt, eig_id, q_eig_id, dfpt_eig_id, iDir, iDtype, killcont, l_real, .NOT.final_SH_it, dfpt_tag)
         CALL timestop("dfpt eigen")

         IF (l_minusq) THEN
            CALL timestart("dfpt minus eigen")
            CALL dfpt_eigen_new(fi, sphhar, results, resultsmq, results1m, fmpi, enpara, nococonv, starsmq, vTot1m, vTot1mIm, &
                                vTot, rho, bmqpt, eig_id, qm_eig_id, dfpt_eigm_id, iDir, iDtype, killcont, l_real, .NOT.final_SH_it, dfpt_tag//"m")
            CALL timestop("dfpt minus eigen")
         END IF

#ifdef CPP_MPI
         CALL MPI_BARRIER(fmpi%mpi_comm, ierr)
#endif

         CALL timestart("Fermi energy and occupation derivative")
         IF (norm2(bqpt)<1e-8) THEN
            CALL dfpt_fermie(eig_id,dfpt_eig_id,fmpi,fi%kpts,fi%input,fi%noco,results,results1)
         ELSE
            results1%ef = 0.0
            results1%w_iks = 0.0
         END IF
         CALL timestop("Fermi energy and occupation derivative")

         IF (l_minusq) THEN
            CALL timestart("Fermi energy and occupation minus derivative")
            IF (norm2(bqpt)<1e-8) THEN
               CALL dfpt_fermie(eig_id,dfpt_eigm_id,fmpi,fi%kpts,fi%input,fi%noco,results,results1m)
            ELSE
               results1m%ef = 0.0
               results1m%w_iks = 0.0
            END IF
            CALL timestop("Fermi energy and occupation minus derivative")
         END IF

#ifdef CPP_MPI
         CALL MPI_BCAST(results1%ef, 1, MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
         CALL MPI_BCAST(results1%w_iks, SIZE(results1%w_iks), MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
#endif

         IF (l_minusq) THEN
#ifdef CPP_MPI
            CALL MPI_BCAST(results1m%ef, 1, MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
            CALL MPI_BCAST(results1m%w_iks, SIZE(results1m%w_iks), MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
#endif
         END IF

         IF (final_SH_it) THEN
            denIn1%mt(:,0:,iDtype,:) = denIn1%mt(:,0:,iDtype,:) + grRho%mt(:,0:,iDtype,:)
            CALL denIn1%distribute(fmpi%mpi_comm)

            IF (l_minusq) THEN
               denIn1m%mt(:,0:,iDtype,:) = denIn1m%mt(:,0:,iDtype,:) + grRho%mt(:,0:,iDtype,:)
               CALL denIn1m%distribute(fmpi%mpi_comm)
            END IF

            l_cont = .FALSE.
            IF (fmpi%irank==0) write(*,*) "Final Sternheimer iteration finished."
            CALL timestop("Sternheimer Iteration")
            CYCLE scfloop
         END IF

         CALL timestart("generation of new charge density (total)")
         CALL denOut1%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_DEN, l_dfpt=.TRUE.)
         CALL denOut1Im%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_DEN, l_dfpt=.FALSE.)
         denOut1%iter = denIn1%iter
         IF (.NOT.l_minusq) THEN
            CALL dfpt_cdngen(eig_id,q_eig_id,dfpt_eig_id,fmpi,fi%input,banddosdummy,fi%vacuum,&
                             fi%kpts,fi%atoms,sphhar,starsq,fi%sym,fi%gfinp,fi%hub1inp,&
                             enpara,fi%cell,fi%noco,nococonv,vTot,results,results1,&
                             archiveType,xcpot,denOut1,denOut1Im,bqpt,iDtype,iDir,l_real)
         ELSE
            CALL dfpt_cdngen(eig_id,q_eig_id,dfpt_eig_id,fmpi,fi%input,banddosdummy,fi%vacuum,&
                             fi%kpts,fi%atoms,sphhar,starsq,fi%sym,fi%gfinp,fi%hub1inp,&
                             enpara,fi%cell,fi%noco,nococonv,vTot,results,results1,&
                             archiveType,xcpot,denOut1,denOut1Im,bqpt,iDtype,iDir,l_real,&
                             qm_eig_id,dfpt_eigm_id,starsmq,results1m)
         END IF
         CALL timestop("generation of new charge density (total)")

         IF (l_minusq) THEN
            CALL timestart("generation of new charge density (total)")
            CALL denOut1m%init(starsmq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_DEN, l_dfpt=.TRUE.)
            CALL denOut1mIm%init(starsmq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_DEN, l_dfpt=.FALSE.)
            CALL dfpt_cdngen(eig_id,qm_eig_id,dfpt_eigm_id,fmpi,fi%input,banddosdummy,fi%vacuum,&
                             fi%kpts,fi%atoms,sphhar,starsmq,fi%sym,fi%gfinp,fi%hub1inp,&
                             enpara,fi%cell,fi%noco,nococonv,vTot,results,results1m,&
                             archiveType,xcpot,denOut1m,denOut1mIm,-bqpt,iDtype,iDir,l_real,&
                             q_eig_id,dfpt_eig_id,starsq,results1)
            CALL timestop("generation of new charge density (total)")
         END IF

         IF (strho) THEN
            strho = .FALSE.
            denIn1 = denOut1
            denIn1Im = denOut1Im
            denIn1%mt(:,0:,iDtype,:) = denIn1%mt(:,0:,iDtype,:) - grRho%mt(:,0:,iDtype,:)
            IF (fmpi%irank==0) write(*,*) "Starting perturbation generated."
            CALL timestop("Sternheimer Iteration")
            CALL save_npy(TRIM(dfpt_tag)//"rho1pw0.npy",denOut1%pw)
            CALL save_npy(TRIM(dfpt_tag)//"rho1mtre0.npy",denOut1%mt)
            CALL save_npy(TRIM(dfpt_tag)//"rho1mtim0.npy",denOut1Im%mt)
            IF (l_minusq) THEN
               denIn1m = denOut1m
               denIn1mIm = denOut1mIm
               CALL save_npy(TRIM(dfpt_tag)//"rho1Mpw0.npy",denOut1m%pw)
               CALL save_npy(TRIM(dfpt_tag)//"rho1Mmtre0.npy",denOut1m%mt)
               CALL save_npy(TRIM(dfpt_tag)//"rho1Mmtim0.npy",denOut1mIm%mt)
               denIn1m%mt(:,0:,iDtype,:) = denIn1m%mt(:,0:,iDtype,:) - grRho%mt(:,0:,iDtype,:)
            END IF
            !STOP
            !RETURN
            CYCLE scfloop
         END IF

         IF (.NOT.onedone) THEN
            onedone = .TRUE.
            denIn1 = denOut1
            denIn1Im = denOut1Im
            CALL save_npy(TRIM(dfpt_tag)//"_rho1pw.npy",denOut1%pw(:,1))
            CALL save_npy(TRIM(dfpt_tag)//"_rho1mtre.npy",denOut1%mt(:,0:,1,1))
            CALL save_npy(TRIM(dfpt_tag)//"_rho1mtim.npy",denOut1Im%mt(:,0:,1,1))
            denIn1%mt(:,0:,iDtype,:) = denIn1%mt(:,0:,iDtype,:) - grRho%mt(:,0:,iDtype,:)
            IF (fmpi%irank==0) write(*,*) "1st 'real' density perturbation generated."
            CALL timestop("Sternheimer Iteration")
            IF (l_minusq) THEN
               denIn1m = denOut1m
               denIn1mIm = denOut1mIm
               denIn1m%mt(:,0:,iDtype,:) = denIn1m%mt(:,0:,iDtype,:) - grRho%mt(:,0:,iDtype,:)
            END IF
            CYCLE scfloop
         END IF

         IF (final_SH_it) THEN
            denIn1 = denOut1
            denIn1Im = denOut1Im
            l_cont = .FALSE.
            CALL save_npy(TRIM(dfpt_tag)//"_fin_rho1pw.npy",denOut1%pw(:,1))
            CALL save_npy(TRIM(dfpt_tag)//"_fin_rho1mtre.npy",denOut1%mt(:,0:,1,1))
            CALL save_npy(TRIM(dfpt_tag)//"_fin_rho1mtim.npy",denOut1Im%mt(:,0:,1,1))
            IF (fmpi%irank==0) write(*,*) "Final Sternheimer iteration finished."
            CALL timestop("Sternheimer Iteration")
            IF (l_minusq) THEN
               denIn1m = denOut1m
               denIn1mIm = denOut1mIm
            END IF
            CYCLE scfloop
         END IF

         CALL denIn1%distribute(fmpi%mpi_comm)
         CALL denIn1Im%distribute(fmpi%mpi_comm)

         IF (l_minusq) THEN
            CALL denIn1m%distribute(fmpi%mpi_comm)
            CALL denIn1mIm%distribute(fmpi%mpi_comm)
         END IF

#ifdef CPP_MPI
         CALL MPI_BARRIER(fmpi%mpi_comm, ierr)
#endif

         field2 = fi%field

         ! First mixing in the 2nd "real" iteration.
         denIn1%mt(:,0:,iDtype,:) = denIn1%mt(:,0:,iDtype,:) + grRho%mt(:,0:,iDtype,:)
         CALL denIn1%distribute(fmpi%mpi_comm)

         IF (l_minusq) THEN
            denIn1m%mt(:,0:,iDtype,:) = denIn1m%mt(:,0:,iDtype,:) + grRho%mt(:,0:,iDtype,:)
            CALL denIn1m%distribute(fmpi%mpi_comm)
         END IF

#ifdef CPP_MPI
         CALL MPI_BARRIER(fmpi%mpi_comm, ierr)
#endif

         ! mix input and output densities
         CALL timestart("DFPT mixing")
         CALL mix_charge(field2, fmpi, (iter == fi%input%itmax .OR. judft_was_argument("-mix_io")), starsq, &
                         fi%atoms, sphhar, fi%vacuum, fi%input, fi%sym, fi%cell, fi%noco, nococonv, &
                         archiveType, xcpot, iter, denIn1, denOut1, results1, .FALSE., fi%sliceplot,&
                         denIn1Im, denOut1Im, dfpt_tag)
         CALL timestop("DFPT mixing")

         denIn1%mt(:,0:,iDtype,:) = denIn1%mt(:,0:,iDtype,:) - grRho%mt(:,0:,iDtype,:)

         CALL denIn1%distribute(fmpi%mpi_comm)

         IF (l_minusq) THEN
            CALL timestart("DFPT mixing")
            CALL mix_charge(field2, fmpi, (iter == fi%input%itmax .OR. judft_was_argument("-mix_io")), starsmq, &
                            fi%atoms, sphhar, fi%vacuum, fi%input, fi%sym, fi%cell, fi%noco, nococonv, &
                            archiveType, xcpot, iterm, denIn1m, denOut1m, results1m, .FALSE., fi%sliceplot,&
                            denIn1mIm, denOut1mIm, dfpt_tag//"m")
            CALL timestop("DFPT mixing")

            denIn1m%mt(:,0:,iDtype,:) = denIn1m%mt(:,0:,iDtype,:) - grRho%mt(:,0:,iDtype,:)

            CALL denIn1m%distribute(fmpi%mpi_comm)
         END IF

#ifdef CPP_MPI
         CALL MPI_BARRIER(fmpi%mpi_comm, ierr)
#endif

         IF (fmpi%irank==0) THEN
            WRITE (oUnit, FMT=8130) iter
8130        FORMAT(/, 5x, '******* it=', i3, '  is completed********', /,/)
            WRITE (*, *) "Iteration:", iter, " Distance:", results1%last_distance
         END IF ! fmpi%irank==0
         CALL timestop("Sternheimer Iteration")

#ifdef CPP_MPI
         CALL MPI_BCAST(results1%last_distance, 1, MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
         CALL MPI_BARRIER(fmpi%mpi_comm, ierr)
#endif

         IF (l_minusq) THEN
#ifdef CPP_MPI
            CALL MPI_BCAST(results1m%last_distance, 1, MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
            CALL MPI_BARRIER(fmpi%mpi_comm, ierr)
#endif
         END IF

         l_cont = l_cont .AND. (iter < fi%input%itmax)
         l_cont = l_cont .AND. ((fi%input%mindistance <= results1%last_distance))

         final_SH_it = fi%input%mindistance > results1%last_distance
         l_cont = l_cont .OR. final_SH_it ! DO one more iteration so V1, z1 and rho1 match

      END DO scfloop ! DO WHILE (l_cont)

      CALL add_usage_data("Iterations", iter)

      CALL close_eig(dfpt_eig_id)
      IF (l_minusq) THEN
         CALL close_eig(dfpt_eigm_id)
      END IF

   END SUBROUTINE dfpt_sternheimer
END MODULE m_dfpt_sternheimer
