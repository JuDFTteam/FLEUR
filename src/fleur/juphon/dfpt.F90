!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt
   USE m_juDFT
   USE m_constants
   USE m_types

   IMPLICIT NONE

CONTAINS
   SUBROUTINE dfpt(fi, sphhar, stars, nococonv, qpts, fmpi, results, enpara, &
                 & rho, vTot, vxc, eig_id, xcpot, hybdat, mpdata, forcetheo)

      USE m_juDFT_stop, only : juDFT_error
      USE m_eig66_io, only : open_eig,close_eig
      USE m_dfpt_check
      USE m_dfpt_interpolation
      use m_types_dfpt
      use m_types_phonon
      use m_types_efield
      use m_types_BEC
      use m_dfpt_postprocess_pot
      


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

      TYPE(t_potden)                :: grRho3(3), grVtot3(3), grVC3(3), grVext3(3)
      TYPE(t_potden)                :: grgrVext3x3(3,3)
      TYPE(t_results)               :: q_results, results1, qm_results, results1m

      INTEGER :: nspins 
      INTEGER :: q_eig_id, dfpt_eig_id, dfpt_eig_id2, qm_eig_id, dfpt_eigm_id, dfpt_eigm_id2
      LOGICAL :: l_real, l_minusq,l_gamma 

      class(t_dfpt), allocatable :: phonon_obj , efield_obj,BEC_obj ! we might want to have multiple scf calculation in one fleur call
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
     
      if (fmpi%irank==0) call dfpt_check(fi, xcpot)

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

      if (fi%juPhon%l_scf) then 
         if (fi%juPhon%l_efield) then 
            call timestart("dfpt efield")
            ! Do a scf calculation with an electric field as the external perturbation
            allocate(t_efield_pert :: efield_obj)
            call efield_obj%init(fi,fi%juPhon%qvec_efield)
            call sternheimerJob%init(fi,l_efield=.true.)
            call efield_obj%perform_scf(sternheimerJob,fi,fmpi,stars,sphhar,xcpot,forcetheo,enpara,nococonv,hybdat,fi%juPhon,rho,vTot,vxc,results,q_results,results1,eig_id,q_eig_id,dfpt_eig_id, &
                                     dfpt_eig_id2,l_minusq,qm_results,results1m,qm_eig_id,dfpt_eigm_id,dfpt_eigm_id2)
            call efield_obj%write_outfiles(fi,fmpi,fi%juPhon)
            call timestop("dfpt efield")
         end if 
         if (fi%juPhon%l_borneffcharge) then
            call timestart("dfpt born effective charges") 
            ! Do a scf calculation for phonon for q-Vectors close to Gamma --> calculate the polar response
            allocate(t_BEC :: BEC_obj)
            call BEC_obj%init(fi,fi%juPhon%qvec_efield)
            call sternheimerJob%init(fi,l_BEC=.true.)
            call BEC_obj%perform_scf(sternheimerJob,fi,fmpi,stars,sphhar,xcpot,forcetheo,enpara,nococonv,hybdat,fi%juPhon,rho,vTot,vxc,results,q_results,results1,eig_id,q_eig_id,dfpt_eig_id, &
                                     dfpt_eig_id2,l_minusq,qm_results,results1m,qm_eig_id,dfpt_eigm_id,dfpt_eigm_id2)
            call BEC_obj%write_outfiles(fi,fmpi,fi%juPhon)
            call timestop("dfpt born effective charges")
         end if 
         if (fi%juPhon%l_phonon) then 
            allocate(t_phonon :: phonon_obj)
            call timestart("dfpt phonons")
            ! Do a scf calculation with atom displacements as the perturbation
            call phonon_obj%init(fi,fi%juPhon%qvec)
            call sternheimerJob%init(fi,l_phonon=.true.)
            call phonon_obj%perform_scf(sternheimerJob,fi,fmpi,stars,sphhar,xcpot,forcetheo,enpara,nococonv,hybdat,fi%juPhon,rho,vTot,vxc,results,q_results,results1,eig_id,q_eig_id,dfpt_eig_id, &
                                      dfpt_eig_id2,l_minusq,qm_results,results1m,qm_eig_id,dfpt_eigm_id,dfpt_eigm_id2)
            call phonon_obj%write_outfiles(fi,fmpi,fi%juPhon)
            call timestop("dfpt phonons")
         end if 
      end if 

#ifdef CPP_MPI
      call MPI_BARRIER(fmpi%mpi_comm,ierr)
#endif 
      if ( fi%juPhon%l_band .or. fi%juPhon%l_dos .or. fi%juPhon%l_intp) then 
         ! Do post processing of converged results 
         ! Interpolate the dynamic matrix with FFT
         call timestart("dfpt interpolation")
         call dfpt_interpolation(fi,fmpi,nococonv,results)
         call timestop("dfpt interpolation")
      end if 

      if ( fi%juPhon%l_elph .and. .not. fi%juPhon%l_scf) then 
         ! Construct the electron phonon matrix element from converged potentials
         call timestart("construction of el-ph matrix elements")
         call construct_elph_mat(fmpi,fi,stars,sphhar,xcpot,forcetheo,enpara,nococonv,hybdat,rho,vTot,vxc,results,eig_id,q_results,q_eig_id,l_real)
         call timestop("construction of el-ph matrix elements")
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

   END SUBROUTINE dfpt

   SUBROUTINE dfpt_desym(fmpi_nosym,fi_nosym,sphhar_nosym,stars_nosym,nococonv_nosym,enpara_nosym,results_nosym,wann_nosym,hybdat_nosym,mpdata_nosym,xcpot_nosym,forcetheo_nosym,rho_nosym,vTot_nosym,grid,inp_pref,&
                         fi,sphhar,stars,nococonv,enpara,results,rho,vTot)
      USE m_desymmetrizer
      USE m_outcdn
      USE m_plot
      USE m_fleur_init

      TYPE(t_mpi),        INTENT(INOUT) :: fmpi_nosym
      TYPE(t_fleurinput), INTENT(INOUT) :: fi_nosym
      TYPE(t_sphhar),     INTENT(INOUT) :: sphhar_nosym
      TYPE(t_stars),      INTENT(INOUT) :: stars_nosym
      TYPE(t_nococonv),   INTENT(INOUT) :: nococonv_nosym
      TYPE(t_enpara),     INTENT(INOUT) :: enpara_nosym
      TYPE(t_results),    INTENT(INOUT) :: results_nosym
      TYPE(t_wann),       INTENT(INOUT) :: wann_nosym
      TYPE(t_hybdat),     INTENT(INOUT) :: hybdat_nosym
      TYPE(t_mpdata),     INTENT(INOUT) :: mpdata_nosym

      CLASS(t_xcpot),     ALLOCATABLE, INTENT(INOUT) :: xcpot_nosym
      CLASS(t_forcetheo), ALLOCATABLE, INTENT(INOUT) :: forcetheo_nosym

      TYPE(t_potden), INTENT(INOUT):: rho_nosym, vTot_nosym
      TYPE(t_fleurinput), INTENT(IN) :: fi
      TYPE(t_sphhar),     INTENT(IN) :: sphhar
      TYPE(t_stars),      INTENT(IN) :: stars
      TYPE(t_nococonv),   INTENT(IN) :: nococonv
      TYPE(t_enpara),     INTENT(IN) :: enpara
      TYPE(t_results),    INTENT(IN) :: results
      TYPE(t_potden),   INTENT(IN) :: rho, vTot

      INTEGER, INTENT(IN) :: grid(3)

      CHARACTER(len=100), INTENT(IN) :: inp_pref

      INTEGER :: ix, iy, iz, iv_old, iflag_old, iv_new, iflag_new
      INTEGER :: iType_old, iAtom_old, iType_new, iAtom_new!, inversionOp
      REAL    :: old_point(3), new_point(3), pt_old(3), pt_new(3), xdnout_old, xdnout_new!, atom_shift(3)
      LOGICAL :: test_desym
      CALL fleur_init(fmpi_nosym, fi_nosym, sphhar_nosym, stars_nosym, nococonv_nosym, forcetheo_nosym, &
                        enpara_nosym, xcpot_nosym, results_nosym, wann_nosym, hybdat_nosym, mpdata_nosym, &
                        inp_pref)

      CALL rho_nosym%init(stars_nosym,fi_nosym%atoms,sphhar_nosym,fi_nosym%vacuum,fi_nosym%noco,fi%input%jspins,POTDEN_TYPE_DEN)
      CALL vTot_nosym%init(stars_nosym,fi_nosym%atoms,sphhar_nosym,fi_nosym%vacuum,fi_nosym%noco,fi%input%jspins,POTDEN_TYPE_POTTOT)

      ! TODO: Correctly account for such a shift in the desymmetrization.
      ! For now: Just build input, that does not necessitate a shift.
      !        inversionOp = -1
      !        symOpLoop: DO iSym = 1, sym%nop
      !           IF (ALL(sym%mrot(:,:,iSym)==invs_matrix)) THEN
      !              inversionOp = iSym
      !              EXIT symOpLoop
      !           END IF
      !        END DO symOpLoop

      !        atom_shift = 0.0
      !        IF (inversionOp.GT.0) THEN
      !           IF(ANY(ABS(sym%tau(:,inversionOp)).GT.eps7).and..not.(film.and.ABS(sym%tau(3,inversionOp))>eps7)) THEN
      !              atom_shift = 0.5*sym%tau(:,inversionOp)
      !           END IF
      !        END IF

      ALLOCATE(vTot_nosym%pw_w, mold=vTot_nosym%pw)
      vTot_nosym%pw_w = CMPLX(0.0,0.0)

      CALL desymmetrize_pw(fi%sym, stars, stars_nosym, rho%pw, rho_nosym%pw)
      CALL desymmetrize_pw(fi%sym, stars, stars_nosym, vTot%pw, vTot_nosym%pw, vTot%pw_w, vTot_nosym%pw_w)
      CALL desymmetrize_mt(fi%sym, fi_nosym%sym, fi%cell, fi%atoms, fi_nosym%atoms, sphhar, sphhar_nosym, rho%mt, rho_nosym%mt)
      CALL desymmetrize_mt(fi%sym, fi_nosym%sym, fi%cell, fi%atoms, fi_nosym%atoms, sphhar, sphhar_nosym, vTot%mt, vTot_nosym%mt)

      CALL desymmetrize_types(fi%input, fi_nosym%input, fi%atoms, fi_nosym%atoms, fi%noco, &
                              nococonv, nococonv_nosym, enpara, enpara_nosym, results, results_nosym)

      test_desym = .FALSE.
      IF (test_desym) THEN
         DO iz = 0, grid(3)-1
            DO iy = 0, grid(2)-1
               DO ix = 0, grid(1)-1
                  old_point = fi%cell%amat(:,1)*REAL(ix)/(grid(1)-1) + &
                              fi%cell%amat(:,2)*REAL(iy)/(grid(2)-1) + &
                              fi%cell%amat(:,3)*REAL(iz)/(grid(3)-1)

                  new_point = fi%cell%amat(:,1)*REAL(ix)/(grid(1)-1) + &
                              fi%cell%amat(:,2)*REAL(iy)/(grid(2)-1) + &
                              fi%cell%amat(:,3)*REAL(iz)/(grid(3)-1)! - &
                              !atom_shift

                  ! Set region specific parameters for point
                  ! Get MT sphere for point if point is in MT sphere
                  CALL getMTSphere(fi%input,fi%cell,fi%atoms,old_point,iType_old,iAtom_old,pt_old)
                  CALL getMTSphere(fi_nosym%input,fi_nosym%cell,fi_nosym%atoms,new_point,iType_new,iAtom_new,pt_new)
                  IF (iAtom_old.NE.0) THEN
                     iv_old = 0
                     iflag_old = 1
                  ELSE
                     iv_old = 0
                     iflag_old = 2
                     pt_old(:) = old_point(:)
                  END IF

                  IF (iAtom_new.NE.0) THEN
                     iv_new = 0
                     iflag_new = 1
                  ELSE
                     iv_new = 0
                     iflag_new = 2
                     pt_new(:) = new_point(:)
                  END IF

                  ! Old point:
                  CALL outcdn(pt_old,iType_old,iAtom_old,iv_old,iflag_old,1,.FALSE.,stars,&
                              fi%vacuum,sphhar,fi%atoms,fi%sym,fi%cell ,&
                              rho,xdnout_old)
                  ! New point:
                  CALL outcdn(pt_new,iType_new,iAtom_new,iv_old,iflag_old,1,.FALSE.,stars_nosym,&
                              fi_nosym%vacuum,sphhar_nosym,fi_nosym%atoms,fi_nosym%sym,fi_nosym%cell ,&
                              rho_nosym,xdnout_new)

                  WRITE(9004,*) xdnout_new-xdnout_old

                  ! Old point:
                  CALL outcdn(pt_old,iType_old,iAtom_old,iv_old,iflag_old,1,.TRUE.,stars,&
                              fi%vacuum,sphhar,fi%atoms,fi%sym,fi%cell ,&
                              vTot,xdnout_old)
                  ! New point:
                  CALL outcdn(pt_new,iType_new,iAtom_new,iv_old,iflag_old,1,.TRUE.,stars_nosym,&
                              fi_nosym%vacuum,sphhar_nosym,fi_nosym%atoms,fi_nosym%sym,fi_nosym%cell ,&
                              vTot_nosym,xdnout_new)
                  WRITE(9005,*) xdnout_new-xdnout_old
               END DO !x-loop
            END DO !y-loop
         END DO !z-loop
      END IF

   END SUBROUTINE

END MODULE m_dfpt
