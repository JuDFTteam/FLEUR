!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_calc_hybrid
   USE m_judft
   use m_store_load_hybrid
CONTAINS

   SUBROUTINE calc_hybrid(fi,mpdata,hybdat,fmpi,nococonv,stars,enpara,&
                          xcpot,v,iter, iterHF)
      use m_work_package
      use m_set_coul_participation
      USE m_types_hybdat
      USE m_types
      USE m_mixedbasis
      USE m_coulombmatrix
      USE m_hf_init
      USE m_hf_setup
      USE m_hsfock
      USE m_io_hybrid
      USE m_eig66_io
      use m_eig66_mpi
      use m_distribute_mpi 
      use m_create_coul_comms
      use m_eigvec_setup
      use m_distrib_vx
#ifdef CPP_MPI 
      use mpi 
#endif
#ifdef CPP_PROG_THREAD
      use m_thread_lib
#endif

      IMPLICIT NONE

      type(t_fleurinput), intent(in)    :: fi
      type(t_mpdata), intent(inout)     :: mpdata
      TYPE(t_hybdat), INTENT(INOUT)     :: hybdat
      TYPE(t_mpi), INTENT(IN)           :: fmpi
      TYPE(t_nococonv), INTENT(IN)      :: nococonv
      type(t_stars), intent(in)         :: stars
      TYPE(t_enpara), INTENT(IN)        :: enpara
      TYPE(t_xcpot_inbuild), INTENT(IN) :: xcpot
      TYPE(t_potden), INTENT(IN)        :: v
      INTEGER, INTENT(INOUT)            :: iter, iterHF

      ! local variables
      type(t_hybmpi)           :: glob_mpi, wp_mpi
      type(t_work_package)     :: work_pack(fi%input%jspins)
      INTEGER                  :: jsp, nk, err, i, wp_rank, ierr, ik
      INTEGER                  :: j_wp, n_wps, root_comm
      INTEGER                  :: max_band_pack, jq
      type(t_lapw)             :: lapw
      LOGICAL                  :: init_vex = .TRUE. !In first call we have to init v_nonlocal
      character(len=999)       :: msg
      REAL, ALLOCATABLE        :: eig_irr(:, :)
      integer, allocatable     :: vx_loc(:,:), weights(:)
      type(c_ptr)              :: threadId
      type(t_mat), allocatable :: vx_tmp(:,:)

      CALL timestart("hybrid code")

#ifdef CPP_MPI
#ifdef CPP_PROG_THREAD
      if(fmpi%l_mpi_multithreaded) call start_prog_thread(threadId)
#endif
#endif

      IF (fi%kpts%nkptf == 0) THEN
         CALL judft_error("kpoint-set of full BZ not available", &
                          hint="to generate fi%kpts in the full BZ you should specify a k-mesh in inp.xml")
      END IF

      !Check if new non-local potential shall be generated
      hybdat%l_subvxc = fi%hybinp%l_hybrid .AND. (.NOT. xcpot%is_name("exx"))
      !If this is the first iteration loop we can not calculate a new non-local potential
      !hybdat%l_calhf = (results%last_distance >= 0.0) .AND. (results%last_distance < fi%input%minDistance)
      !make sure we do at least one PBE first
      if(iter == 1 .and. iterHF == 0) hybdat%l_calhf = .False.

      IF (.NOT. hybdat%l_calhf) THEN
         hybdat%l_subvxc = hybdat%l_subvxc .AND. hybdat%l_addhf
      else
         call glob_mpi%init(fmpi%mpi_comm)
         hybdat%results%te_hfex%core = 0

         !Check if we are converged well enough to calculate a new potential
         hybdat%l_addhf = .TRUE.

         !In first iteration allocate some memory
         IF (init_vex) THEN
            call first_iteration_alloc(fi, hybdat)
            init_vex = .FALSE.
         END IF
         hybdat%l_subvxc = (hybdat%l_subvxc .AND. hybdat%l_addhf)
         IF (.NOT. ALLOCATED(hybdat%results%w_iks)) allocate(hybdat%results%w_iks(fi%input%neig, fi%kpts%nkpt, fi%input%jspins))

         iterHF = iterHF + 1

         !Delete broyd files
         CALL system("rm -f broyd*")

         !check if z-reflection trick can be used

        
         CALL timestart("Preparation for hybrid functionals")
         !construct the mixed-basis
         CALL timestart("generation of mixed basis")
         if(glob_mpi%rank == 0) write (*,*) "iterHF =    " // int2str(iterHF)
         CALL mixedbasis(fi%atoms, fi%kpts,  fi%input, fi%cell, xcpot, fi%mpinp, mpdata, fi%hybinp, hybdat,&
                        enpara, fmpi, v, iterHF)
         CALL timestop("generation of mixed basis")

         ! setup parallelization 
         n_wps = min(glob_mpi%size, fi%kpts%nkpt)
         allocate(weights(n_wps), source=0)
         do j_wp = 1, n_wps
            do ik = j_wp, fi%kpts%nkpt, n_wps
               weights(j_wp) =  weights(j_wp) + fi%kpts%eibz(ik)%nkpt
            enddo
         enddo
         call distribute_mpi(weights, glob_mpi, wp_mpi, wp_rank)
         call hybdat%set_nobd(fi)
         call hybdat%set_nbands(fi, fmpi)
         do jsp = 1,fi%input%jspins
            call work_pack(jsp)%init(fi, hybdat, mpdata, wp_mpi, jsp, wp_rank, n_wps)
         enddo

         if(.not. allocated(hybdat%zmat)) allocate(hybdat%zmat(fi%kpts%nkptf, fi%input%jspins))

         DO jsp = 1, fi%input%jspins
            DO nk = 1,fi%kpts%nkptf
!               IF (hybdat%zmat(nk, jsp)%mat%matsize2 .NE. hybdat%nbands(nk, jsp)) THEN ! This IF caused deadlocks.
                  CALL lapw%init(fi%input, fi%noco, nococonv,fi%kpts, fi%atoms, fi%sym, nk, fi%cell)
                  call eigvec_setup(hybdat%zmat(nk, jsp), fi, lapw, work_pack, fmpi, &
                                    hybdat%nbands(nk, jsp), nk, jsp, hybdat%eig_id)
!               END IF
            enddo 
         enddo
         call bcast_eigvecs(hybdat, fi, nococonv, fmpi)

         if(.not. allocated(hybdat%coul)) allocate(hybdat%coul(fi%kpts%nkpt))
         call set_coul_participation(hybdat, fi, fmpi, work_pack)
         call create_coul_comms(hybdat, fi, fmpi)

         do i =1,fi%kpts%nkpt
            if(hybdat%coul(i)%l_participate) then 
               call hybdat%coul(i)%alloc(fi, mpdata%num_radbasfn, mpdata%n_g, i, .false.)
            else 
               call hybdat%coul(i)%mini_alloc(fi)
            endif 
         enddo 

         ! use jsp=1 for coulomb work-planning
         CALL coulombmatrix(fmpi, fi, mpdata, hybdat, xcpot)
         
         do i =1,fi%kpts%nkpt
            if(hybdat%coul(i)%l_participate) then 
               call hybdat%coul(i)%mpi_bcast(fi, hybdat%coul(i)%comm, 0)
            endif
         enddo

         CALL hf_init(mpdata, fi, hybdat)
         CALL timestop("Preparation for hybrid functionals")

         call judft_comm_split(glob_mpi%comm, wp_mpi%rank, 0, root_comm)

         CALL timestart("Calculation of non-local HF potential")
         allocate(vx_loc(fi%kpts%nkpt,fi%input%jspins), source=-1)
         allocate(vx_tmp(fi%kpts%nkpt, fi%input%jspins))
         DO jsp = 1, fi%input%jspins
            CALL HF_setup(mpdata,fi, fmpi, nococonv, jsp, enpara, &
                        hybdat, v%mt(:, 0, :, :), eig_irr)

            call timestart("get max_q")
            hybdat%max_q = 0
            DO i = 1,work_pack(jsp)%k_packs(1)%size
               max_band_pack = 0 
               DO jq = 1, size(work_pack(jsp)%k_packs(i)%q_packs)
                  max_band_pack = max(max_band_pack, work_pack(jsp)%k_packs(i)%q_packs(jq)%submpi%size)
               enddo
               hybdat%max_q = hybdat%max_q + size(work_pack(jsp)%k_packs(i)%q_packs) * fi%kpts%nkptf  * max_band_pack
            END DO
#ifdef CPP_MPI
            call MPI_Allreduce(MPI_IN_PLACE, hybdat%max_q, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif
            hybdat%max_q = hybdat%max_q + 20 ! The 20 is kind of dirty. It is meant as a safety because of the MPI_BARRIER in hsfock, related to the ex_to_vx call.
            call timestop("get max_q")

            DO i = 1,work_pack(jsp)%k_packs(1)%size
               nk = work_pack(jsp)%k_packs(i)%nk
               PRINT*, 'kpoint= ', nk
               CALL lapw%init(fi%input, fi%noco, nococonv,fi%kpts, fi%atoms, fi%sym, nk, fi%cell)
               CALL hsfock(fi, work_pack(jsp)%k_packs(i), mpdata, lapw, jsp, hybdat, eig_irr, &
                           nococonv, stars, xcpot, fmpi, vx_tmp(nk, jsp))
               if(work_pack(jsp)%k_packs(i)%submpi%root()) vx_loc(nk, jsp) = fmpi%irank
            END DO

#ifdef CPP_MPI
            CALL timestart("balancing MPI_Barriers")
            DO WHILE (hybdat%max_q > 0)
               call MPI_Barrier(MPI_COMM_WORLD, ierr)
               hybdat%max_q = hybdat%max_q - 1
            END DO
            CALL timestop("balancing MPI_Barriers")
#endif

            call work_pack(jsp)%free()
         END DO
#ifdef CPP_MPI
         call timestart("MPI_Allred te_hfex%core")
         if(wp_mpi%rank == 0) call MPI_Allreduce(MPI_IN_PLACE, hybdat%results%te_hfex%core, 1, MPI_DOUBLE_PRECISION, MPI_SUM, root_comm, ierr)
         call timestop("MPI_Allred te_hfex%core")
#endif
         CALL timestop("Calculation of non-local HF potential")
#ifdef CPP_MPI
         call MPI_Allreduce(MPI_IN_PLACE, vx_loc, size(vx_loc), MPI_INTEGER, MPI_MAX, fmpi%mpi_comm, ierr)   
#endif
         call distrib_vx(fi, fmpi, nococonv, vx_loc, vx_tmp, hybdat)
         call store_hybrid_data(fi, fmpi, hybdat)

#ifdef CPP_MPI
         call timestart("Hybrid imbalance")
         call MPI_Barrier(fmpi%mpi_comm, err)
         call timestop("Hybrid imbalance")
#endif

      ENDIF

#ifdef CPP_MPI
#ifdef CPP_PROG_THREAD
      if(fmpi%l_mpi_multithreaded) call stop_prog_thread(threadId)
#endif
#endif
      CALL timestop("hybrid code")
   CONTAINS
      subroutine first_iteration_alloc(fi, hybdat)
         implicit none
         type(t_fleurinput), intent(in)    :: fi
         TYPE(t_hybdat), INTENT(INOUT)     :: hybdat

         if(allocated(hybdat%nbands)) then
            deallocate(hybdat%nbands, stat=err, errmsg=msg)
            if(err /= 0) THEN
               write (*,*) "errorcode", err
               write (*,*) "errormessage", msg
            endif
         endif

         allocate(hybdat%nbands(fi%kpts%nkptf, fi%input%jspins), source=0)

         if(allocated(hybdat%nbasm)) deallocate(hybdat%nbasm)
         allocate(hybdat%nbasm(fi%kpts%nkptf), source=0)

         if(allocated(hybdat%div_vv)) deallocate(hybdat%div_vv)
         allocate(hybdat%div_vv(fi%input%neig, fi%kpts%nkpt, fi%input%jspins), source=0.0)
      end subroutine first_iteration_alloc
   END SUBROUTINE calc_hybrid
END MODULE m_calc_hybrid
