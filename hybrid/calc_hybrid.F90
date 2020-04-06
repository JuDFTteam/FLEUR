!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_calc_hybrid
   USE m_judft

CONTAINS

   SUBROUTINE calc_hybrid(eig_id,fi,mpdata,hybdat,mpi,nococonv,enpara,&
                          results,xcpot,v,iterHF)

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

      IMPLICIT NONE

      INTEGER, INTENT(IN)               :: eig_id
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpdata), intent(inout)     :: mpdata
      TYPE(t_hybdat), INTENT(INOUT)     :: hybdat
      TYPE(t_mpi), INTENT(IN)           :: mpi
      TYPE(t_nococonv), INTENT(IN)      :: nococonv
      TYPE(t_enpara), INTENT(IN)        :: enpara
      TYPE(t_results), INTENT(INOUT)    :: results
      TYPE(t_xcpot_inbuild), INTENT(IN) :: xcpot
      TYPE(t_potden), INTENT(IN)        :: v
      INTEGER, INTENT(INOUT)            :: iterHF

      ! local variables
      type(t_hybmpi)    :: hybmpi
      INTEGER           :: jsp, nk, err, i
      type(t_lapw)      :: lapw
      LOGICAL           :: init_vex = .TRUE. !In first call we have to init v_nonlocal
      LOGICAL           :: l_zref
      character(len=999):: msg
      REAL, ALLOCATABLE :: eig_irr(:, :)
      INTEGER, ALLOCATABLE :: my_k_list(:), k_owner(:)

      CALL timestart("hybrid code")
      call sync_eig(eig_id)

      call hybmpi%copy_mpi(mpi)
      call split_k_to_comm(fi, hybmpi, my_k_list, k_owner)

      INQUIRE (file="v_x.1", exist=hybdat%l_addhf)

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
         CALL timestop("hybrid code")
         RETURN
      ENDIF

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
      write (*,*) "iterHF = ", iterHF
      CALL mixedbasis(fi%atoms, fi%kpts,  fi%input, fi%cell, xcpot, fi%mpinp, mpdata, fi%hybinp, hybdat,&
                      enpara, mpi, v, iterHF)
      CALL timestop("generation of mixed basis")


      if(.not. allocated(hybdat%coul)) allocate(hybdat%coul(fi%kpts%nkpt))
      do i =1,fi%kpts%nkpt
         call hybdat%coul(i)%alloc(fi, mpdata%num_radbasfn, mpdata%n_g)
      enddo

      CALL coulombmatrix(mpi, fi, mpdata, hybdat, xcpot, my_k_list)

      do i =1,fi%kpts%nkpt
         call hybdat%coul(i)%mpi_ibc(fi, hybmpi, k_owner(i))
      enddo

      CALL hf_init(eig_id, mpdata, fi, hybdat)
      CALL timestop("Preparation for hybrid functionals")

      CALL timestart("Calculation of non-local HF potential")
      DO jsp = 1, fi%input%jspins
         call timestart("HF_setup")
         CALL HF_setup(mpdata,fi%hybinp, fi%input, fi%sym, fi%kpts,  fi%atoms, &
                       mpi, fi%noco, nococonv,fi%cell, fi%oneD, results, jsp, enpara, &
                       hybdat, fi%sym%invs, v%mt(:, 0, :, :), eig_irr)
         call timestop("HF_setup")

         DO i = 1,size(my_k_list)
            nk = my_k_list(i)
            CALL lapw%init(fi%input, fi%noco, nococonv,fi%kpts, fi%atoms, fi%sym, nk, fi%cell, l_zref)
            CALL hsfock(fi,nk, mpdata, lapw, jsp, hybdat, eig_irr, &
                        nococonv, results, xcpot, mpi)
         END DO
      END DO
      CALL timestop("Calculation of non-local HF potential")
#ifdef CPP_MPI
      call timestart("Hybrid imbalance")
      call MPI_Barrier(mpi%mpi_comm, err)
      call timestop("Hybrid imbalance")
#endif

      call sync_eig(eig_id)
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

      subroutine split_k_to_comm(fi, hybmpi, my_k_list, k_owner)
         implicit none

         type(t_fleurinput), intent(in)      :: fi
         type(t_hybmpi), intent(in)          :: hybmpi
         integer, allocatable, intent(inout) :: my_k_list(:)
         integer, allocatable, intent(inout)  :: k_owner(:)
         integer   :: i, irank

         if(allocated(my_k_list)) deallocate(my_k_list)
         allocate(my_k_list(0))

         if(fi%kpts%nkpt < hybmpi%size) call judft_error("not enough k-points for mpis")
         
         ! get my k-list
         do i = hybmpi%rank+1,fi%kpts%nkpt,hybmpi%size
            my_k_list = [my_k_list, i]
         enddo

         ! findout who's got the other k's
         if(.not. allocated(k_owner)) allocate(k_owner(fi%kpts%nkpt), source=-1)

         do irank = 0,hybmpi%size-1
            do i = irank+1,fi%kpts%nkpt,hybmpi%size
               k_owner(i) = irank
            enddo
         enddo

         ! sanity check
         do i = 1,size(my_k_list)
            if(k_owner(my_k_list(i)) /= hybmpi%rank) then 
               write (*,*) "my_k_list", my_k_list 
               write (*,*) "i = ", i
               write (*,*) "k_owner", k_owner
               call judft_error("I should own my own k-point")
            endif
         enddo

      end subroutine split_k_to_comm
   END SUBROUTINE calc_hybrid
END MODULE m_calc_hybrid
