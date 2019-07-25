!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_cdngen
CONTAINS

SUBROUTINE cdngen(eig_id,mpi,input,banddos,sliceplot,vacuum,&
                  dimension,kpts,atoms,sphhar,stars,sym,&
                  enpara,cell,noco,vTot,results,oneD,coreSpecInput,&
                  archiveType, xcpot,outDen,EnergyDen,gOnsite,hub1)

   !*****************************************************
   !    Charge density generator
   !    calls cdnval to generate the valence charge and the
   !    core routines for the core contribution
   !*****************************************************

   USE m_types
   USE m_constants
   USE m_juDFT
   USE m_prpqfftmap
   USE m_cdnval
   USE m_cdn_io
   USE m_wrtdop
   USE m_cdntot
   USE m_qfix
   USE m_genNewNocoInp
   USE m_xmlOutput
   USE m_magMoms
   USE m_orbMagMoms
   USE m_cdncore
   USE m_doswrite
   USE m_Ekwritesl
   USE m_banddos_io
   USE m_metagga
   USE m_unfold_band_kpts
   USE m_gfcalc
   USE m_onsite
   USE m_angles
   USE m_hubbard1_io
   USE m_denmat_dist
   USE m_triang
#ifdef CPP_MPI
   USE m_mpi_bc_potden
#endif

   IMPLICIT NONE

#ifdef CPP_MPI
   INCLUDE 'mpif.h'
#endif

   ! Type instance arguments
   TYPE(t_results),INTENT(INOUT)    :: results
   TYPE(t_mpi),INTENT(IN)           :: mpi
   TYPE(t_dimension),INTENT(IN)     :: dimension
   TYPE(t_oneD),INTENT(IN)          :: oneD
   TYPE(t_enpara),INTENT(INOUT)     :: enpara
   TYPE(t_banddos),INTENT(IN)       :: banddos
   TYPE(t_sliceplot),INTENT(IN)     :: sliceplot
   TYPE(t_input),INTENT(IN)         :: input
   TYPE(t_vacuum),INTENT(IN)        :: vacuum
   TYPE(t_noco),INTENT(INOUT)       :: noco
   TYPE(t_sym),INTENT(IN)           :: sym
   TYPE(t_stars),INTENT(IN)         :: stars
   TYPE(t_cell),INTENT(IN)          :: cell
   TYPE(t_kpts),INTENT(IN)          :: kpts
   TYPE(t_sphhar),INTENT(IN)        :: sphhar
   TYPE(t_atoms),INTENT(IN)         :: atoms
   TYPE(t_coreSpecInput),INTENT(IN) :: coreSpecInput
   TYPE(t_potden),INTENT(IN)        :: vTot
   TYPE(t_greensf),OPTIONAL,INTENT(INOUT)    :: gOnsite
   TYPE(t_hub1ham),OPTIONAL,INTENT(INOUT)    :: hub1
   CLASS(t_xcpot),INTENT(INOUT)     :: xcpot
   TYPE(t_potden),INTENT(INOUT)     :: outDen, EnergyDen

   !Scalar Arguments
   INTEGER, INTENT (IN)             :: eig_id, archiveType

   ! Local type instances
   TYPE(t_noco)          :: noco_new
   TYPE(t_regionCharges) :: regCharges, fake_regCharges
   TYPE(t_dos)           :: dos, fake_dos
   TYPE(t_moments)       :: moments, fake_moments
   TYPE(t_results)       :: fake_results
   TYPE(t_mcd)           :: mcd
   TYPE(t_slab)          :: slab
   TYPE(t_orbcomp)       :: orbcomp
   TYPE(t_cdnvalJob)     :: cdnvalJob
   TYPE(t_greensfCoeffs) :: greensfCoeffs
   TYPE(t_potden)        :: val_den, core_den


   !Local Scalars
   REAL                  :: fix, qtot, dummy,eFermiPrev
   INTEGER               :: jspin, jspmax, ierr
   INTEGER               :: dim_idx

   INTEGER               :: i_gf,m,l
   REAL                  :: n_occ
   COMPLEX               :: mmpmat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_gf,4)

#ifdef CPP_HDF
   INTEGER(HID_T)        :: banddosFile_id
#endif
   LOGICAL               :: l_error, perform_MetaGGA
   
   LOGICAL l_tria
   INTEGER ntria
   REAL    as
   INTEGER itria(3,2*kpts%nkpt)
   REAL    atr(2*kpts%nkpt)
   REAL    angle(sym%nop)



   CALL regCharges%init(input,atoms)
   CALL dos%init(input,atoms,dimension,kpts,vacuum)
   CALL moments%init(input,atoms)
   CALL mcd%init1(banddos,dimension,input,atoms,kpts)
   CALL slab%init(banddos,dimension,atoms,cell,input,kpts)
   CALL orbcomp%init(input,banddos,dimension,atoms,kpts)
   

   IF(atoms%n_gf.GT.0.AND.PRESENT(gOnsite)) THEN
      !Only calculate the greens function when needed
      CALL greensfCoeffs%init(input,lmaxU_const,atoms,noco,results%ef)
      IF(input%tria.AND.input%film) THEN
         l_tria = .true.
         CALL triang(kpts%bk,kpts%nkpt,itria,ntria,atr,as,l_tria)
         IF (sym%invs) THEN
           IF (abs(sym%nop2*as-0.5).GT.0.000001) l_tria=.false.
         ELSE
           IF (abs(sym%nop2*as-1.0).GT.0.000001) l_tria=.false.
         ENDIF
         write(*,*) as,sym%nop2,l_tria
         IF(.NOT.l_tria) CALL juDFT_warn("l_tria=F",calledby="cdngen")
      ENDIF
      IF(mpi%irank==0) THEN
         CALL gOnsite%e_contour(input,greensfCoeffs%e_bot,greensfCoeffs%e_top,results%ef)
         IF(atoms%n_hia.GT.0) hub1%mag_mom = 0.0
      ENDIF
   ENDIF

   IF(atoms%n_gf+atoms%n_u.GT.0.AND.noco%l_mperp) THEN
      CALL angles(sym,angle)
   ENDIF


   CALL outDen%init(stars,    atoms, sphhar, vacuum, noco, input%jspins, POTDEN_TYPE_DEN)
   CALL EnergyDen%init(stars, atoms, sphhar, vacuum, noco, input%jspins, POTDEN_TYPE_EnergyDen)

   IF (mpi%irank == 0) CALL openXMLElementNoAttributes('valenceDensity')

   !In a non-collinear calcuation where the off-diagonal part of the
   !density matrix in the muffin-tins is calculated, the a- and
   !b-coef. for both spins are needed at once. Thus, cdnval is only
   !called once and both spin directions are calculated in a single run.
   results%force=0.0
   jspmax = input%jspins
   IF (noco%l_mperp.OR.input%l_gfmperp) jspmax = 1
   DO jspin = 1,jspmax
      CALL cdnvalJob%init(mpi,input,kpts,noco,results,jspin,sliceplot,banddos)
      CALL cdnval(eig_id,mpi,kpts,jspin,noco,input,banddos,cell,atoms,enpara,stars,vacuum,dimension,&
                  sphhar,sym,vTot,oneD,cdnvalJob,outDen,regCharges,dos,results,moments,hub1,coreSpecInput,&
                  mcd,slab,orbcomp,greensfCoeffs,angle,ntria,as,itria,atr)
   END DO

   IF(PRESENT(gOnsite).AND.mpi%irank.EQ.0) THEN
      IF(atoms%n_gf.GT.0) THEN
         !Perform the Kramer-Kronigs-Integration
         CALL calc_onsite(atoms,input,noco,results%ef,greensfCoeffs,gOnsite,sym)
         !calculate the crystal field contribution to the local hamiltonian in LDA+Hubbard 1
         IF(atoms%n_hia.GT.0.AND.ANY(hub1%ccf(:).NE.0.0)) THEN
           CALL crystal_field(atoms,input,greensfCoeffs,hub1,vTot%mmpMat(:,:,atoms%n_u+1:atoms%n_u+atoms%n_hia,:))
         ENDIF
         CALL hybridization(3,1,gOnsite,atoms,input)
         CALL gfDOS(gOnsite,3,1,1,atoms,input)
         IF(input%jspins.EQ.2) THEN
            CALL eff_excinteraction(gOnsite,atoms,input,greensfCoeffs)
         ENDIF
         CALL occmtx(gOnsite,3,1,atoms,sym,input,mmpmat(:,:,1,:),l_write=.TRUE.)
      ENDIF
   ENDIF

   call val_den%copyPotDen(outDen)
   ! calculate kinetic energy density for MetaGGAs
   if(xcpot%exc_is_metagga()) then
      CALL calc_EnergyDen(eig_id, mpi, kpts, noco, input, banddos, cell, atoms, enpara, stars,&
                             vacuum, DIMENSION, sphhar, sym, vTot, oneD, results, EnergyDen)
   endif

   IF (mpi%irank == 0) THEN
      IF (banddos%dos.or.banddos%vacdos.or.input%cdinf) THEN
         IF (banddos%unfoldband) THEN
            eFermiPrev = 0.0
            CALL readPrevEFermi(eFermiPrev,l_error)
            CALL write_band_sc(kpts,results,eFermiPrev)
         END IF
#ifdef CPP_HDF
         CALL openBandDOSFile(banddosFile_id,input,atoms,cell,kpts,banddos)
         CALL writeBandDOSData(banddosFile_id,input,atoms,cell,kpts,results,banddos,dos,vacuum)
         CALL closeBandDOSFile(banddosFile_id)
#endif
         CALL timestart("cdngen: dos")
         CALL doswrite(eig_id,dimension,kpts,atoms,vacuum,input,banddos,sliceplot,noco,sym,cell,dos,mcd,results,slab,orbcomp,oneD)
         IF (banddos%dos.AND.(banddos%ndir == -3)) THEN
            CALL Ek_write_sl(eig_id,dimension,kpts,atoms,vacuum,input,jspmax,sym,cell,dos,slab,orbcomp,results)
         END IF
         CALL timestop("cdngen: dos")
      END IF
   END IF

   IF ((banddos%dos.OR.banddos%vacdos).AND.(banddos%ndir/=-2)) CALL juDFT_end("DOS OK",mpi%irank)
   IF (vacuum%nstm == 3) CALL juDFT_end("VACWAVE OK",mpi%irank)

   IF (mpi%irank.EQ.0) THEN
      CALL cdntot(stars,atoms,sym,vacuum,input,cell,oneD,outDen,.TRUE.,qtot,dummy)
      CALL closeXMLElement('valenceDensity')
   END IF ! mpi%irank = 0

   IF (sliceplot%slice) THEN
      IF (mpi%irank == 0) THEN
         CALL writeDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,CDN_ARCHIVE_TYPE_CDN_const,CDN_INPUT_DEN_const,&
                           0,-1.0,0.0,.FALSE.,outDen,'cdn_slice')
      END IF
      CALL juDFT_end("slice OK",mpi%irank)
   END IF

   CALL timestart("cdngen: cdncore")
   if(xcpot%exc_is_MetaGGA()) then
      CALL cdncore(mpi,dimension,oneD,input,vacuum,noco,sym,&
                   stars,cell,sphhar,atoms,vTot,outDen,moments,results, EnergyDen)
   else
      CALL cdncore(mpi,dimension,oneD,input,vacuum,noco,sym,&
                   stars,cell,sphhar,atoms,vTot,outDen,moments,results)
   endif
   call core_den%subPotDen(outDen, val_den)
   CALL timestop("cdngen: cdncore")

   CALL enpara%calcOutParams(input,atoms,vacuum,regCharges)

   IF (mpi%irank == 0) THEN
      CALL openXMLElementNoAttributes('allElectronCharges')
      CALL qfix(mpi,stars,atoms,sym,vacuum,sphhar,input,cell,oneD,outDen,noco%l_noco,.TRUE.,.true.,fix)
      CALL closeXMLElement('allElectronCharges')

      IF (input%jspins == 2) THEN
         noco_new = noco

         !Calculate and write out spin densities at the nucleus and magnetic moments in the spheres
         CALL magMoms(dimension,input,atoms,noco_new,vTot,moments)

         noco = noco_new

         !Generate and save the new nocoinp file if the directions of the local
         !moments are relaxed or a constraint B-field is calculated.
         IF (ANY(noco%l_relax(:atoms%ntype)).OR.noco%l_constr) THEN
            CALL genNewNocoInp(input,atoms,noco,noco_new)
         END IF

         IF (noco%l_soc) CALL orbMagMoms(input,atoms,noco,moments%clmom)
         
      END IF
   END IF ! mpi%irank == 0
   
   perform_MetaGGA = ALLOCATED(EnergyDen%mt) &
                   .AND. (xcpot%exc_is_MetaGGA() .or. xcpot%vx_is_MetaGGA())
   if(perform_MetaGGA) then
      call set_kinED(mpi, sphhar, atoms, sym, core_den, val_den, xcpot, &
                     input, noco, stars, cell, outDen, EnergyDen, vTot)
   endif
#ifdef CPP_MPI
   CALL MPI_BCAST(noco%l_ss,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(noco%l_mperp,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(noco%l_constr,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(noco%mix_b,1,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)

   CALL MPI_BCAST(noco%alphInit,atoms%ntype,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(noco%alph,atoms%ntype,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(noco%beta,atoms%ntype,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(noco%b_con,atoms%ntype*2,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(noco%l_relax,atoms%ntype,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
   CALL MPI_BCAST(noco%qss,3,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)

   CALL mpi_bc_potden(mpi,stars,sphhar,atoms,input,vacuum,oneD,noco,outDen)
#endif

END SUBROUTINE cdngen

END MODULE m_cdngen
