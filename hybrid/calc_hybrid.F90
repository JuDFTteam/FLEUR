!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_calc_hybrid
   USE m_judft

CONTAINS

   SUBROUTINE calc_hybrid(eig_id,fi,mpdata,hybdat,fmpi,nococonv,stars,enpara,&
                          results,xcpot,v,iterHF)
      use m_work_package
      USE m_types_hybdat
      USE m_types
      USE m_mixedbasis
      USE m_coulombmatrix
      USE m_hf_init
      USE m_hf_setup
      USE m_hsfock
      USE m_io_hybinp
      USE m_eig66_io
      use m_eig66_mpi
#ifdef CPP_MPI 
      use mpi 
#endif

      IMPLICIT NONE

      INTEGER, INTENT(IN)               :: eig_id
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpdata), intent(inout)     :: mpdata
      TYPE(t_hybdat), INTENT(INOUT)     :: hybdat
      TYPE(t_mpi), INTENT(IN)           :: fmpi
      TYPE(t_nococonv), INTENT(IN)      :: nococonv
      type(t_stars), intent(in)         :: stars
      TYPE(t_enpara), INTENT(IN)        :: enpara
      TYPE(t_results), INTENT(INOUT)    :: results
      TYPE(t_xcpot_inbuild), INTENT(IN) :: xcpot
      TYPE(t_potden), INTENT(IN)        :: v
      INTEGER, INTENT(INOUT)            :: iterHF

      ! local variables
      type(t_hybmpi)    :: glob_mpi, wp_mpi, tmp_mpi
      type(t_work_package) :: work_pack
      INTEGER           :: jsp, nk, err, i, wp_rank, wp_size, tmp_comm
      type(t_lapw)      :: lapw
      LOGICAL           :: init_vex = .TRUE. !In first call we have to init v_nonlocal
      LOGICAL           :: l_zref
      character(len=999):: msg
      REAL, ALLOCATABLE :: eig_irr(:, :)

      CALL timestart("hybrid code")

      IF (fi%kpts%nkptf == 0) THEN
         CALL judft_error("kpoint-set of full BZ not available", &
                          hint="to generate fi%kpts in the full BZ you should specify a k-mesh in inp.xml")
      END IF

      !Check if new non-local potential shall be generated
      hybdat%l_subvxc = fi%hybinp%l_hybrid .AND. (.NOT. xcpot%is_name("exx"))
      !If this is the first iteration loop we can not calculate a new non-local potential
      hybdat%l_calhf = (results%last_distance >= 0.0) .AND. (results%last_distance < fi%input%minDistance)
      IF (.NOT. hybdat%l_calhf) THEN
         hybdat%l_subvxc = hybdat%l_subvxc .AND. hybdat%l_addhf
      else
         call glob_mpi%init(fmpi%mpi_comm)
         results%te_hfex%core = 0

         !Check if we are converged well enough to calculate a new potential
         hybdat%l_addhf = .TRUE.

         !In first iteration allocate some memory
         IF (init_vex) THEN
            call first_iteration_alloc(fi, hybdat)
            init_vex = .FALSE.
         END IF
         hybdat%l_subvxc = (hybdat%l_subvxc .AND. hybdat%l_addhf)
         IF (.NOT. ALLOCATED(results%w_iks)) allocate(results%w_iks(fi%input%neig, fi%kpts%nkpt, fi%input%jspins))

         iterHF = iterHF + 1

         !Delete broyd files
         CALL system("rm -f broyd*")

         !check if z-reflection trick can be used

         l_zref = (fi%sym%zrfs .AND. (SUM(ABS(fi%kpts%bk(3, :fi%kpts%nkpt))) < 1e-9) .AND. .NOT. fi%noco%l_noco)

         CALL timestart("Preparation for hybrid functionals")
         !construct the mixed-basis
         CALL timestart("generation of mixed basis")
         if(glob_mpi%rank == 0) write (*,*) "iterHF =    " // int2str(iterHF)
         CALL mixedbasis(fi%atoms, fi%kpts,  fi%input, fi%cell, xcpot, fi%mpinp, mpdata, fi%hybinp, hybdat,&
                        enpara, fmpi, v, iterHF)
         CALL timestop("generation of mixed basis")


         if(.not. allocated(hybdat%coul)) allocate(hybdat%coul(fi%kpts%nkpt))
         do i =1,fi%kpts%nkpt
            call hybdat%coul(i)%alloc(fi, mpdata%num_radbasfn, mpdata%n_g, i)
         enddo

         ! use jsp=1 for coulomb work-planning
         CALL coulombmatrix(fmpi, fi, mpdata, hybdat, xcpot)

         do i =1,fi%kpts%nkpt
            call hybdat%coul(i)%mpi_ibc(fi, fmpi%mpi_comm, fmpi%coulomb_owner(i))
         enddo

         CALL hf_init(eig_id, mpdata, fi, hybdat)
         CALL timestop("Preparation for hybrid functionals")

         call distrib_mpis(fi, glob_mpi, wp_mpi, wp_rank, wp_size)

         CALL timestart("Calculation of non-local HF potential")
         DO jsp = 1, fi%input%jspins
            call timestart("HF_setup")
            CALL HF_setup(mpdata,fi, fmpi, nococonv, results, jsp, enpara, &
                        hybdat, v%mt(:, 0, :, :), eig_irr)
            call timestop("HF_setup")

            call work_pack%init(fi, hybdat, wp_mpi, jsp, wp_rank, wp_size)
            
            DO i = 1,work_pack%k_packs(1)%size
               nk = work_pack%k_packs(i)%nk
               CALL lapw%init(fi%input, fi%noco, nococonv,fi%kpts, fi%atoms, fi%sym, nk, fi%cell, l_zref)
               CALL hsfock(fi, work_pack%k_packs(i), mpdata, lapw, jsp, hybdat, eig_irr, &
                           nococonv, stars, results, xcpot, fmpi)
            END DO
            call work_pack%free()
         END DO
         CALL timestop("Calculation of non-local HF potential")
#ifdef CPP_MPI
         call timestart("Hybrid imbalance")
         call MPI_Barrier(fmpi%mpi_comm, err)
         call timestop("Hybrid imbalance")
#endif

      ENDIF
      CALL timestop("hybrid code")
   CONTAINS
      subroutine first_iteration_alloc(fi, hybdat)
         implicit none
         type(t_fleurinput), intent(in)    :: fi
         TYPE(t_hybdat), INTENT(INOUT)     :: hybdat

         if(allocated(hybdat%ne_eig)) deallocate(hybdat%ne_eig)
         allocate(hybdat%ne_eig(fi%kpts%nkpt), source=0)

         if(allocated(hybdat%nbands)) then
            deallocate(hybdat%nbands, stat=err, errmsg=msg)
            if(err /= 0) THEN
               write (*,*) "errorcode", err
               write (*,*) "errormessage", msg
            endif
         endif

         allocate(hybdat%nbands(fi%kpts%nkptf), source=0)

         if(allocated(hybdat%nobd)) deallocate(hybdat%nobd)
         allocate(hybdat%nobd(fi%kpts%nkptf, fi%input%jspins), source=0)

         if(allocated(hybdat%nbasm)) deallocate(hybdat%nbasm)
         allocate(hybdat%nbasm(fi%kpts%nkptf), source=0)

         if(allocated(hybdat%div_vv)) deallocate(hybdat%div_vv)
         allocate(hybdat%div_vv(fi%input%neig, fi%kpts%nkpt, fi%input%jspins), source=0.0)
      end subroutine first_iteration_alloc

      subroutine distrib_mpis(fi, glob_mpi, wp_mpi, wp_rank, wp_size)
         USE m_types
         implicit none 
         type(t_fleurinput), intent(in)    :: fi
         type(t_hybmpi), intent(in)        :: glob_mpi
         type(t_hybmpi), intent(inout)     :: wp_mpi
         integer, intent(inout)            :: wp_rank, wp_size
   
         integer :: n_wps, i, j, cnt, j_wp, ik, idx(1), new_comm
         integer, allocatable :: nprocs(:), weights(:), color(:)
   
   
         n_wps = min(glob_mpi%size, fi%kpts%nkpt)
         allocate(nprocs(n_wps), source=0)
         allocate(weights(n_wps), source=0)
         allocate(color(glob_mpi%size), source=0)
   
         do j_wp = 1, n_wps
            do ik = j_wp, fi%kpts%nkpt, n_wps
               weights(j_wp) =  weights(j_wp) + fi%kpts%eibz(ik)%nkpt
            enddo
         enddo
   
         do i = 1,glob_mpi%size 
            idx = minloc(1.0*nprocs/weights)
            nprocs(idx(1)) = nprocs(idx(1)) + 1
         enddo
   
         cnt = 1
         do i = 1,n_wps 
            do j = 1,nprocs(i)
               color(cnt) = i - 1
               cnt = cnt + 1 
            enddo 
         enddo
   
         wp_rank = color(glob_mpi%rank+1) 
         wp_size = n_wps
         
         call judft_comm_split(glob_mpi%comm, wp_rank, glob_mpi%rank, new_comm)
   
         call wp_mpi%init(new_comm)
      end subroutine distrib_mpis
   END SUBROUTINE calc_hybrid
END MODULE m_calc_hybrid
