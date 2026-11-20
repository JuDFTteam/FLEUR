!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt
   USE m_juDFT
   USE m_constants
   USE m_juDFT_stop, ONLY: juDFT_error
   USE m_eig66_io, ONLY: open_eig, close_eig
   USE m_dfpt_check
   USE m_dfpt_interpolation
   USE m_types_dfpt_scf
   USE m_types_phonon
   USE m_types_efield
   USE m_types_BEC
   USE m_types_bfield
   USE m_dfpt_postprocess_pot
   USE m_desym_setup
   USE m_desym_potden
   USE m_desym_eigvec
   USE m_types_desym_map
   USE m_types_lapw
   USE m_types_enpara
   USE m_types_fleurinput
   USE m_types_forcetheo
   USE m_types_hybdat
   USE m_types_kpts
   USE m_types_mpdata
   USE m_types_mpi
   USE m_types_nococonv
   USE m_types_potden
   USE m_types_misc
   USE m_types_sphhar
   USE m_types_stars
   USE m_types_sternheimerjob
   USE m_types_sym
   USE m_types_xcpot
#ifdef CPP_MPI
   USE mpi
#endif

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: dfpt

CONTAINS
   subroutine dfpt(fi, sphhar, stars, nococonv, qpts, fmpi, results, enpara, rho, vTot, vxc, eig_id, xcpot, hybdat, mpdata, forcetheo)
      !! DFPT driver. A symmetric ground state is desymmetrized in memory and the DFPT kernel runs on the setup without symmetry.
      type(t_mpi),        intent(in)    :: fmpi
      type(t_fleurinput), intent(in)    :: fi
      type(t_sphhar),     intent(in)    :: sphhar
      type(t_stars),      intent(in)    :: stars
      type(t_nococonv),   intent(in)    :: nococonv
      type(t_enpara),     intent(inout) :: enpara
      type(t_results),    intent(inout) :: results
      type(t_hybdat),     intent(inout) :: hybdat
      type(t_mpdata),     intent(inout) :: mpdata
      class(t_xcpot),     intent(in)    :: xcpot
      class(t_forcetheo), intent(inout) :: forcetheo
      type(t_kpts),       intent(in)    :: qpts
      type(t_potden),     intent(inout) :: rho
      type(t_potden),     intent(in)    :: vTot, vxc
      integer,            intent(in)    :: eig_id

      type(t_mpi)        :: fmpi_nosym
      type(t_fleurinput) :: fi_nosym
      type(t_stars)      :: stars_nosym
      type(t_sphhar)     :: sphhar_nosym
      type(t_enpara)     :: enpara_nosym
      type(t_nococonv)   :: nococonv_nosym
      type(t_desym_map)  :: map
      type(t_potden)     :: rho_nosym, vTot_nosym, vxc_nosym
      type(t_results)    :: results_nosym
      integer :: eig_id_nosym

      if (fmpi%irank == 0) call dfpt_check(fi, xcpot)

      if (fi%sym%nop == 1) then
         call dfpt_main(fi, sphhar, stars, nococonv, qpts, fmpi, results, enpara, rho, vTot, vxc, eig_id, xcpot, hybdat, mpdata, forcetheo)
         return
      end if

      call timestart("desymmetrization")
      call desym_setup(fmpi, fi, enpara, nococonv, forcetheo, fmpi_nosym, fi_nosym, stars_nosym, sphhar_nosym, enpara_nosym, nococonv_nosym, map)
      call desym_potden(fi, stars, sphhar, fi_nosym, stars_nosym, sphhar_nosym, map, rho, rho_nosym)
      call desym_potden(fi, stars, sphhar, fi_nosym, stars_nosym, sphhar_nosym, map, vTot, vTot_nosym)
      call desym_potden(fi, stars, sphhar, fi_nosym, stars_nosym, sphhar_nosym, map, vxc, vxc_nosym)
      call desym_eigvec(fi, enpara, nococonv, vTot, eig_id, results, fmpi_nosym, fi_nosym, enpara_nosym, nococonv_nosym, vTot_nosym, &
                        map, eig_id_nosym, results_nosym)
      call timestop("desymmetrization")

      call dfpt_main(fi_nosym, sphhar_nosym, stars_nosym, nococonv_nosym, fi_nosym%kpts, fmpi_nosym, results_nosym, enpara_nosym, &
                       rho_nosym, vTot_nosym, vxc_nosym, eig_id_nosym, xcpot, hybdat, mpdata, forcetheo, sym_full=fi%sym)
      call close_eig(eig_id_nosym)

   end subroutine dfpt

   SUBROUTINE dfpt_main(fi, sphhar, stars, nococonv, qpts, fmpi, results, enpara, &
                 & rho, vTot, vxc, eig_id, xcpot, hybdat, mpdata, forcetheo, sym_full)
      !! DFPT perturbations, interpolation and el-ph postprocessing for a setup without symmetry.
      !! sym_full: symmetry of the desymmetrized input, used with fi%dfpt%qvec to unfold the q-mesh.
      TYPE(t_mpi),        INTENT(IN)     :: fmpi
      TYPE(t_fleurinput), INTENT(IN)     :: fi
      TYPE(t_sphhar),     INTENT(IN)     :: sphhar
      TYPE(t_stars),      INTENT(IN)     :: stars
      TYPE(t_nococonv),   INTENT(IN)     :: nococonv
      TYPE(t_enpara),     INTENT(INOUT)  :: enpara
      TYPE(t_results),    INTENT(INOUT)  :: results
      TYPE(t_hybdat),     INTENT(INOUT)  :: hybdat
      TYPE(t_mpdata),     INTENT(INOUT)  :: mpdata

      CLASS(t_xcpot),     INTENT(IN)     :: xcpot
      CLASS(t_forcetheo), INTENT(INOUT)  :: forcetheo

      TYPE(t_kpts),       INTENT(IN)  :: qpts !Possibly replace this by fi%kpts [read correctly!]

      TYPE(t_potden),   INTENT(INOUT) :: rho
      TYPE(t_potden),   INTENT(IN)    :: vTot, vxc
      INTEGER,          INTENT(IN)    :: eig_id
      TYPE(t_sym),      INTENT(IN), OPTIONAL :: sym_full

      TYPE(t_potden)                :: grRho3(3), grVtot3(3), grVC3(3), grVext3(3)
      TYPE(t_potden)                :: grgrVext3x3(3,3)
      TYPE(t_results)               :: q_results, results1, qm_results, results1m

      INTEGER :: nspins 
      INTEGER :: q_eig_id, dfpt_eig_id, dfpt_eig_id2, qm_eig_id, dfpt_eigm_id, dfpt_eigm_id2
      LOGICAL :: l_real, l_minusq,l_gamma 

      class(t_dfpt_scf), allocatable :: phonon_obj, efield_obj, BEC_obj, bfield_obj ! we might want to have multiple scf calculation in one fleur call
      type(t_sternheimerJob)     :: sternheimerJob
#ifdef CPP_MPI
      integer :: ierr
#endif 

      INTEGER, ALLOCATABLE :: q_list(:)
      ! INTEGER, ALLOCATABLE :: sym_list(:) ! For each q: Collect, which symmetries leave q unchanged.
      
      l_real = fi%sym%invs.AND.(.NOT.fi%noco%l_soc).AND.(.NOT.fi%noco%l_noco).AND.fi%atoms%n_hia==0

      ! l_minusq is a hard false at the moment. It can be used to ignore +-q symmetries and
      ! run the Sternheimer loop etc etc for both +q and -q.
      ! This was only used for testing but may become relevant again for SOC systems with broken
      ! inversion symmetry
      l_minusq = .FALSE.

      nspins = MERGE(2, 1, fi%noco%l_noco)

      if (fmpi%n_size > 1) call juDFT_error("DFPT needs one MPI rank per k-point.", calledby="dfpt_main", &
                                            hint="Run with -pe_per_kpt 1 or choose the number of ranks as a divisor of the number of k-points.")

      ! q_results saves the eigen-info for k+q, results1 for the perturbed quantities
      CALL q_results%init(fi%input, fi%atoms, fi%kpts, fi%noco)
      CALL results1%init(fi%input, fi%atoms, fi%kpts, fi%noco)
      
      IF (l_minusq) THEN
         CALL qm_results%init(fi%input, fi%atoms, fi%kpts, fi%noco)
         CALL results1m%init(fi%input, fi%atoms, fi%kpts, fi%noco)
      END IF

      ! The eig_ids, where the stuff of k+q, the perturbed stuff and some extra dynmat stuff will be saved.
      q_eig_id = open_eig(fmpi%mpi_comm, lapw_dim_nbasfcn, fi%input%neig, fi%kpts%nkpt, fi%input%jspins, fi%noco%l_noco, &
                        .NOT.fi%INPUT%eig66(1), fi%input%l_real, fi%noco%l_soc, fi%input%eig66(1), .FALSE., fmpi%n_size)
      dfpt_eig_id = open_eig(fmpi%mpi_comm, lapw_dim_nbasfcn, fi%input%neig, fi%kpts%nkpt, fi%input%jspins, fi%noco%l_noco, &
                             .NOT.fi%INPUT%eig66(1), .FALSE., fi%noco%l_soc, fi%INPUT%eig66(1), .FALSE., fmpi%n_size)
      dfpt_eig_id2 = open_eig(fmpi%mpi_comm, lapw_dim_nbasfcn, fi%input%neig, fi%kpts%nkpt, fi%input%jspins, fi%noco%l_noco, &
                              .NOT.fi%INPUT%eig66(1), .FALSE., fi%noco%l_soc, fi%INPUT%eig66(1), .FALSE., fmpi%n_size)

      IF (l_minusq) THEN
         qm_eig_id = open_eig(fmpi%mpi_comm, lapw_dim_nbasfcn, fi%input%neig, fi%kpts%nkpt, fi%input%jspins, fi%noco%l_noco, &
                            .NOT.fi%INPUT%eig66(1), fi%input%l_real, fi%noco%l_soc, fi%input%eig66(1), .FALSE., fmpi%n_size)
         dfpt_eigm_id = open_eig(fmpi%mpi_comm, lapw_dim_nbasfcn, fi%input%neig, fi%kpts%nkpt, fi%input%jspins, fi%noco%l_noco, &
                                .NOT.fi%INPUT%eig66(1), .FALSE., fi%noco%l_soc, fi%INPUT%eig66(1), .FALSE., fmpi%n_size)
         dfpt_eigm_id2 = open_eig(fmpi%mpi_comm, lapw_dim_nbasfcn, fi%input%neig, fi%kpts%nkpt, fi%input%jspins, fi%noco%l_noco, &
                                .NOT.fi%INPUT%eig66(1), .FALSE., fi%noco%l_soc, fi%INPUT%eig66(1), .FALSE., fmpi%n_size)
      END IF

      if (fi%dfpt%l_scf) then 
         if (fi%dfpt%l_bfield) then 
            !Zeeman field
            call timestart("dfpt bfield")
            allocate(t_bfield :: bfield_obj)
            call bfield_obj%init(fi,fi%dfpt%qvec%bk)
            call sternheimerJob%init(fi,l_bfield=.true.)
            call bfield_obj%perform_scf(sternheimerJob,fi,fmpi,stars,sphhar,xcpot,forcetheo,enpara,nococonv,hybdat,fi%dfpt,rho,vTot,vxc,results,q_results,results1,eig_id,q_eig_id,dfpt_eig_id, &
                                     dfpt_eig_id2,l_minusq,qm_results,results1m,qm_eig_id,dfpt_eigm_id,dfpt_eigm_id2)
            call bfield_obj%write_outfiles(fi,fmpi,fi%dfpt)
            call timestop("dfpt bfield")
         end if 
         if (fi%dfpt%l_efield) then 
            call timestart("dfpt efield")
            ! Do a scf calculation with an electric field as the external perturbation
            allocate(t_efield_pert :: efield_obj)
            call efield_obj%init(fi,fi%dfpt%qvec_efield)
            call sternheimerJob%init(fi,l_efield=.true.)
            call efield_obj%perform_scf(sternheimerJob,fi,fmpi,stars,sphhar,xcpot,forcetheo,enpara,nococonv,hybdat,fi%dfpt,rho,vTot,vxc,results,q_results,results1,eig_id,q_eig_id,dfpt_eig_id, &
                                     dfpt_eig_id2,l_minusq,qm_results,results1m,qm_eig_id,dfpt_eigm_id,dfpt_eigm_id2)
            call efield_obj%write_outfiles(fi,fmpi,fi%dfpt)
            call timestop("dfpt efield")
         end if 
         if (fi%dfpt%l_borneffcharge) then
            call timestart("dfpt born effective charges") 
            ! Do a scf calculation for phonon for q-Vectors close to Gamma --> calculate the polar response
            allocate(t_BEC :: BEC_obj)
            call BEC_obj%init(fi,fi%dfpt%qvec_efield)
            call sternheimerJob%init(fi,l_BEC=.true.)
            call BEC_obj%perform_scf(sternheimerJob,fi,fmpi,stars,sphhar,xcpot,forcetheo,enpara,nococonv,hybdat,fi%dfpt,rho,vTot,vxc,results,q_results,results1,eig_id,q_eig_id,dfpt_eig_id, &
                                     dfpt_eig_id2,l_minusq,qm_results,results1m,qm_eig_id,dfpt_eigm_id,dfpt_eigm_id2)
            call BEC_obj%write_outfiles(fi,fmpi,fi%dfpt)
            call timestop("dfpt born effective charges")
         end if 
         if (fi%dfpt%l_phonon) then 
            allocate(t_phonon :: phonon_obj)
            call timestart("dfpt phonons")
            ! Do a scf calculation with atom displacements as the perturbation
            call phonon_obj%init(fi,fi%dfpt%qvec%bk)
            call sternheimerJob%init(fi,l_phonon=.true.)
            call phonon_obj%perform_scf(sternheimerJob,fi,fmpi,stars,sphhar,xcpot,forcetheo,enpara,nococonv,hybdat,fi%dfpt,rho,vTot,vxc,results,q_results,results1,eig_id,q_eig_id,dfpt_eig_id, &
                                      dfpt_eig_id2,l_minusq,qm_results,results1m,qm_eig_id,dfpt_eigm_id,dfpt_eigm_id2)
            call phonon_obj%write_outfiles(fi,fmpi,fi%dfpt)
            call timestop("dfpt phonons")
         end if 
      end if 

#ifdef CPP_MPI
      call MPI_BARRIER(fmpi%mpi_comm,ierr)
#endif 
      if (fi%dfpt%l_intp) then 
         ! Do post processing of converged results 
         ! Interpolate the dynamic matrix with FFT
         call timestart("dfpt interpolation")
         call dfpt_interpolation(fi,fmpi,nococonv,results,sym_full)
         call timestop("dfpt interpolation")
      end if 

      if ( fi%dfpt%l_postprocess) then 
         ! Construct the matrix element from converged potentials
         if (fi%dfpt%l_elph) then 
            call timestart("construction of el-ph matrix elements")
            call dfpt_postprocess_elph(fmpi,fi,stars,sphhar,xcpot,forcetheo,enpara,nococonv,hybdat,rho,vTot,vxc,results,eig_id,q_results,q_eig_id,l_real,sym_full)
            call timestop("construction of el-ph matrix elements")
         end if 
      end if 


      CALL close_eig(q_eig_id)
      CALL close_eig(dfpt_eig_id)
      CALL close_eig(dfpt_eig_id2)
      IF (l_minusq) THEN
         CALL close_eig(qm_eig_id)
         CALL close_eig(dfpt_eigm_id)
         CALL close_eig(dfpt_eigm_id2)
      END IF


      if (fmpi%irank==0) WRITE (oUnit,*) '------------------------------------------------------'

   END SUBROUTINE dfpt_main

END MODULE m_dfpt
