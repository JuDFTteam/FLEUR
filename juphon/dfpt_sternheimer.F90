!--------------------------------------------------------------------------------
! Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_sternheimer
   USE m_types
   USE m_make_stars
   USE m_vgen
   USE m_eigen
   USE m_dfpt_cdngen
   USE m_dfpt_vgen
   USE m_dfpt_fermie
   USE m_mix
   USE m_constants
   USE m_cdn_io
   USE m_eig66_io

IMPLICIT NONE

CONTAINS
   SUBROUTINE dfpt_sternheimer(fi, xcpot, sphhar, stars, starsq, nococonv, qpts, fmpi, results, enpara, hybdat, mpdata, &
                               forcetheo, rho, vTot, grRho, grVtot, grVext, iQ, iDType, iDir, dfpt_tag, eig_id, &
                               l_real, results1, dfpt_eig_id, &
                               denIn1, vTot1, denIn1Im, vTot1Im, vC1, vC1Im)
      TYPE(t_fleurinput), INTENT(IN)    :: fi
      CLASS(t_xcpot),     INTENT(IN)    :: xcpot
      TYPE(t_sphhar),     INTENT(IN)    :: sphhar
      TYPE(t_stars),      INTENT(IN)    :: stars
      TYPE(t_stars),      INTENT(INOUT) :: starsq
      TYPE(t_nococonv),   INTENT(IN)    :: nococonv
      TYPE(t_kpts),       INTENT(IN)    :: qpts
      TYPE(t_mpi),        INTENT(IN)    :: fmpi
      TYPE(t_results),    INTENT(INOUT) :: results, results1
      TYPE(t_hybdat),     INTENT(INOUT) :: hybdat
      TYPE(t_mpdata),     INTENT(INOUT) :: mpdata
      CLASS(t_forcetheo), INTENT(INOUT) :: forcetheo
      TYPE(t_enpara),     INTENT(INOUT) :: enpara
      TYPE(t_potden),     INTENT(IN)    :: rho, vTot, grRho, grVtot, grVext

      TYPE(t_potden), INTENT(INOUT) :: denIn1, vTot1, denIn1Im, vTot1Im, vC1, vC1Im

      TYPE(t_potden) :: denOut1, denOut1Im, vx, rho_loc, rho_loc0

      INTEGER, INTENT(IN) :: iQ, iDtype, iDir, eig_id
      LOGICAL, INTENT(IN) :: l_real
      CHARACTER(len=20), INTENT(IN) :: dfpt_tag

      INTEGER, INTENT(OUT) :: dfpt_eig_id

#ifdef CPP_MPI
      INTEGER :: ierr
#endif

      INTEGER :: archiveType, iter, killcont(6)
      REAL    :: bqpt(3)
      LOGICAL :: l_cont, l_exist, l_lastIter, l_dummy, strho, onedone, final_SH_it, l_exitus

      TYPE(t_hub1data) :: hub1data
      TYPE(t_banddos)  :: banddosdummy
      TYPE(t_field)    :: field2

      killcont = [1,1,1,1,1,1]
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
#endif

      iter = 0
      l_cont = (iter < fi%input%itmax)

      IF (fmpi%irank==0.AND.l_exist) CALL readDensity(starsq, fi%noco, fi%vacuum, fi%atoms, fi%cell, sphhar, &
                                                      fi%input, fi%sym, archiveType, CDN_INPUT_DEN_const, 0, &
                                                      results%ef, results%last_distance, l_dummy, denIn1,  &
                                                      inFilename=TRIM(dfpt_tag),denIm=denIn1Im)

      IF (fmpi%irank==0.AND..NOT.l_exist) denIn1%iter = 1

      CALL vTot1%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
      CALL vTot1Im%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)

      CALL vC1%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
      CALL vC1Im%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)

      dfpt_eig_id = open_eig(fmpi%mpi_comm, lapw_dim_nbasfcn, fi%input%neig, fi%kpts%nkpt, fi%input%jspins, fi%noco%l_noco, &
                             .NOT.fi%INPUT%eig66(1), .FALSE., fi%noco%l_soc, fi%INPUT%eig66(1), .FALSE., fmpi%n_size)

      bqpt = qpts%bk(:, iQ)

!#ifdef CPP_CHASE
!      CALL init_chase(fmpi, fi%input, fi%atoms, fi%kpts, fi%noco, l_real)
!#endif

      l_lastIter = .FALSE.
      scfloop: DO WHILE (l_cont)
         !IF (.NOT.strho) iter = iter + 1 !TODO: Eventually this will be right.
         IF (onedone) iter = iter + 1
         l_lastIter = l_lastIter.OR.(iter.EQ.fi%input%itmax)

         CALL timestart("Iteration")
         IF (fmpi%irank==0.AND.onedone) THEN
            WRITE (oUnit, FMT=8100) iter
8100        FORMAT(/, 10x, '   iter=  ', i5)
         END IF !fmpi%irank==0

#ifdef CPP_CHASE
         CALL chase_distance(results%last_distance)
#endif

         CALL denIn1%distribute(fmpi%mpi_comm)
         CALL denIn1Im%distribute(fmpi%mpi_comm)

         CALL timestart("Generation of potential perturbation")
         IF (strho) THEN
            CALL dfpt_vgen(hybdat,fi%field,fi%input,xcpot,fi%atoms,sphhar,stars,fi%vacuum,fi%sym,&
                           fi%cell ,fi%sliceplot,fmpi,fi%noco,nococonv,rho_loc0,vTot,&
                           starsq,denIn1Im,vTot1,.FALSE.,vTot1Im,denIn1,iDtype,iDir,[1,1]) ! comparison is [1,0]
                           ! TODO: Set this back to [1,1]!
         ELSE
            CALL dfpt_vgen(hybdat,fi%field,fi%input,xcpot,fi%atoms,sphhar,stars,fi%vacuum,fi%sym,&
                           fi%cell ,fi%sliceplot,fmpi,fi%noco,nococonv,rho_loc,vTot,&
                           starsq,denIn1Im,vTot1,.TRUE.,vTot1Im,denIn1,iDtype,iDir,[1,1]) ! comparison is [1,0]
                           ! TODO: Set this back to [1,1]!
         END IF

         IF (final_SH_it) THEN
            CALL dfpt_vgen(hybdat,fi%field,fi%input,xcpot,fi%atoms,sphhar,stars,fi%vacuum,fi%sym,&
                           fi%cell ,fi%sliceplot,fmpi,fi%noco,nococonv,rho_loc,vTot,&
                           starsq,denIn1Im,vC1,.FALSE.,vC1Im,denIn1,iDtype,iDir,[0,0])
            write(*,*) "VC1 filled"
         END IF

         CALL timestop("Generation of potential perturbation")

#ifdef CPP_MPI
         CALL MPI_BARRIER(fmpi%mpi_comm, ierr)
#endif

         IF (strho) THEN
            vTot1%mt(:,0:,iDtype,:) = vTot1%mt(:,0:,iDtype,:) + grVext%mt(:,0:,iDtype,:)
         ELSE
            vTot1%mt(:,0:,iDtype,:) = vTot1%mt(:,0:,iDtype,:) + grVtot%mt(:,0:,iDtype,:)
         END IF

         IF (.NOT.strho) THEN
            killcont = [1,1,1,1,1,1]
         END IF
         CALL timestart("H1 generation (total)")

         CALL timestart("eigen")

         IF (.NOT. fi%input%eig66(1)) THEN
            CALL eigen(fi, fmpi, stars, sphhar, xcpot, forcetheo, enpara, nococonv, mpdata, &
                       hybdat, iter, eig_id, results, rho, vTot, vx, hub1data, &
                       bqpt=bqpt, dfpt_eig_id=dfpt_eig_id, iDir=iDir, iDtype=iDtype, &
                       starsq=starsq, v1real=vTot1, v1imag=vTot1Im, killcont=killcont, l_real=l_real)
         END IF
         CALL timestop("eigen")
         CALL timestop("H1 generation (total)")

#ifdef CPP_MPI
         CALL MPI_BARRIER(fmpi%mpi_comm, ierr)
#endif

         CALL timestart("Fermi energy and occupation derivative")
         CALL dfpt_fermie(eig_id,dfpt_eig_id,fmpi,fi%kpts,fi%input,fi%noco,results,results1)
         CALL timestop("Fermi energy and occupation derivative")
#ifdef CPP_MPI
         CALL MPI_BCAST(results1%ef, 1, MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
         CALL MPI_BCAST(results1%w_iks, SIZE(results1%w_iks), MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
#endif
         CALL timestart("generation of new charge density (total)")
         CALL denOut1%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_DEN, l_dfpt=.TRUE.)
         CALL denOut1Im%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_DEN, l_dfpt=.FALSE.)
         denOut1%iter = denIn1%iter
         CALL dfpt_cdngen(eig_id,dfpt_eig_id,fmpi,fi%input,banddosdummy,fi%vacuum,&
                          fi%kpts,fi%atoms,sphhar,starsq,fi%sym,fi%gfinp,fi%hub1inp,&
                          enpara,fi%cell,fi%noco,nococonv,vTot,results,results1,&
                          archiveType,xcpot,denOut1,denOut1Im,bqpt,iDtype,iDir,l_real)
         CALL timestop("generation of new charge density (total)")

         IF (strho) THEN
            strho = .FALSE.
            denIn1 = denOut1
            denIn1Im = denOut1Im
            denIn1%mt(:,0:,iDtype,:) = denIn1%mt(:,0:,iDtype,:) - grRho%mt(:,0:,iDtype,:)
            write(*,*) "Starting perturbation generated."
            CYCLE scfloop
         END IF

         IF (.NOT.onedone) THEN
            onedone = .TRUE.
            denIn1 = denOut1
            denIn1Im = denOut1Im
            denIn1%mt(:,0:,iDtype,:) = denIn1%mt(:,0:,iDtype,:) - grRho%mt(:,0:,iDtype,:)
            write(*,*) "1st 'real' density perturbation generated."
            CYCLE scfloop
         END IF

         IF (final_SH_it) THEN
            denIn1 = denOut1
            denIn1Im = denOut1Im
            l_cont = .FALSE.
            write(*,*) "Final Sternheimer iteration finished."
            CYCLE scfloop
         END IF

         field2 = fi%field
         ! First mixing in the 2nd "real" iteration.
         denIn1%mt(:,0:,iDtype,:) = denIn1%mt(:,0:,iDtype,:) + grRho%mt(:,0:,iDtype,:)
         ! mix input and output densities
         CALL mix_charge(field2, fmpi, (iter == fi%input%itmax .OR. judft_was_argument("-mix_io")), starsq, &
                         fi%atoms, sphhar, fi%vacuum, fi%input, fi%sym, fi%cell, fi%noco, nococonv, &
                         archiveType, xcpot, iter, denIn1, denOut1, results1, .FALSE., fi%sliceplot,&
                         denIn1Im, denOut1Im, dfpt_tag)

         denIn1%mt(:,0:,iDtype,:) = denIn1%mt(:,0:,iDtype,:) - grRho%mt(:,0:,iDtype,:)
         IF (fmpi%irank==0) THEN
            WRITE (oUnit, FMT=8130) iter
8130        FORMAT(/, 5x, '******* it=', i3, '  is completed********', /,/)
            WRITE (*, *) "Iteration:", iter, " Distance:", results1%last_distance
         END IF ! fmpi%irank==0
         CALL timestop("Iteration")

#ifdef CPP_MPI
         CALL MPI_BCAST(results1%last_distance, 1, MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
         CALL MPI_BARRIER(fmpi%mpi_comm, ierr)
#endif

         l_cont = l_cont .AND. (iter < fi%input%itmax)
         l_cont = l_cont .AND. ((fi%input%mindistance <= results1%last_distance))

         final_SH_it = fi%input%mindistance > results1%last_distance
         l_cont = l_cont .OR. final_SH_it ! DO one more iteration so V1, z1 and rho1 match

      END DO scfloop ! DO WHILE (l_cont)

      CALL add_usage_data("Iterations", iter)

      CALL close_eig(dfpt_eig_id)

   END SUBROUTINE dfpt_sternheimer
END MODULE m_dfpt_sternheimer
