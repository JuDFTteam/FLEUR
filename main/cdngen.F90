!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_cdngen
#ifdef CPP_MPI
   use mpi
#endif
CONTAINS

SUBROUTINE cdngen(eig_id,fmpi,input,banddos,sliceplot,vacuum,&
                  kpts,atoms,sphhar,stars,sym,gfinp,hub1inp,&
                  enpara,cell,noco,nococonv,vTot,results,oneD,coreSpecInput,&
                  archiveType, xcpot,outDen,EnergyDen,greensFunction,hub1data)

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
   USE m_plot
   USE m_cdn_io
   USE m_wrtdop
   USE m_cdntot
   USE m_qfix
   USE m_genNewNocoInp
   USE m_xmlOutput
   USE m_magMoms
   USE m_magMultipoles
   USE m_orbMagMoms
   USE m_resMoms
   USE m_cdncore
   USE m_make_dos
   !USE m_Ekwritesl
   !USE m_banddos_io
   USE m_metagga
   !USE m_unfold_band_kpts
   USE m_denMultipoleExp
   USE m_greensfPostProcess
   USE m_types_greensfContourData
   USE m_angles
   USE m_types_eigdos
#ifdef CPP_MPI
   USE m_mpi_bc_potden
#endif

   IMPLICIT NONE

   ! Type instance arguments
   TYPE(t_results),INTENT(INOUT)    :: results
   TYPE(t_mpi),INTENT(IN)           :: fmpi

   TYPE(t_oneD),INTENT(IN)          :: oneD
   TYPE(t_enpara),INTENT(INOUT)     :: enpara
   TYPE(t_banddos),INTENT(IN)       :: banddos
   TYPE(t_sliceplot),INTENT(IN)     :: sliceplot
   TYPE(t_input),INTENT(IN)         :: input
   TYPE(t_vacuum),INTENT(IN)        :: vacuum
   TYPE(t_noco),INTENT(IN)          :: noco
   TYPE(t_nococonv),INTENT(INOUT)   :: nococonv
   TYPE(t_sym),INTENT(IN)           :: sym
   TYPE(t_stars),INTENT(IN)         :: stars
   TYPE(t_cell),INTENT(IN)          :: cell
   TYPE(t_kpts),INTENT(IN)          :: kpts
   TYPE(t_sphhar),INTENT(IN)        :: sphhar
   TYPE(t_atoms),INTENT(IN)         :: atoms
   TYPE(t_coreSpecInput),INTENT(IN) :: coreSpecInput
   TYPE(t_potden),INTENT(IN)        :: vTot
   TYPE(t_gfinp),INTENT(IN)         :: gfinp
   TYPE(t_hub1inp),INTENT(IN)       :: hub1inp
   TYPE(t_greensf),OPTIONAL,INTENT(INOUT)    :: greensFunction(:)
   TYPE(t_hub1data),OPTIONAL,INTENT(INOUT)    :: hub1data
   CLASS(t_xcpot),INTENT(IN)     :: xcpot
   TYPE(t_potden),INTENT(INOUT)     :: outDen, EnergyDen

   !Scalar Arguments
   INTEGER, INTENT (IN)             :: eig_id, archiveType

   ! Local type instances
   TYPE(t_regionCharges)          :: regCharges
   TYPE(t_dos),TARGET             :: dos
   TYPE(t_vacdos),TARGET          :: vacdos
   TYPE(t_moments)                :: moments
   TYPE(t_mcd),TARGET             :: mcd
   TYPE(t_slab),TARGET            :: slab
   TYPE(t_orbcomp),TARGET         :: orbcomp
   TYPE(t_jDOS),TARGET            :: jDOS
   TYPE(t_cdnvalJob)       :: cdnvalJob
   TYPE(t_greensfImagPart) :: greensfImagPart
   TYPE(t_potden)          :: val_den, core_den
   TYPE(t_greensfContourData) :: contour(gfinp%numberContours)


   !Local Scalars
   REAL                  :: fix, qtot, dummy,eFermiPrev
   INTEGER               :: jspin, ierr
   INTEGER               :: dim_idx
   INTEGER               :: i_gf,iContour

   TYPE(t_eigdos_list) :: eigdos(6)

#ifdef CPP_HDF
   INTEGER(HID_T)        :: banddosFile_id
#endif
   LOGICAL               :: l_error,Perform_metagga

   ! Initialization section
   CALL regCharges%init(input,atoms)
   CALL moments%init(fmpi,input,sphhar,atoms)
   !initalize data for DOS
   CALL dos%init(input,atoms,kpts,banddos,results%eig); eigdos(1)%p=>dos
   CALL vacdos%init(input,atoms,kpts,banddos,results%eig); eigdos(6)%p=>vacdos
   CALL mcd%init(banddos,input,atoms,kpts,results%eig);eigdos(2)%p=>mcd
   CALL slab%init(banddos,atoms,cell,input,kpts);eigdos(3)%p=>slab
   CALL orbcomp%init(input,banddos,atoms,kpts,results%eig);eigdos(4)%p=>orbcomp
   CALL jDOS%init(input,banddos,atoms,kpts);eigdos(5)%p=>jDOS

   CALL outDen%init(stars,    atoms, sphhar, vacuum, noco, input%jspins, POTDEN_TYPE_DEN)
   CALL EnergyDen%init(stars, atoms, sphhar, vacuum, noco, input%jspins, POTDEN_TYPE_EnergyDen)
   results%force=0.0


   IF(PRESENT(greensFunction).AND.gfinp%n.GT.0) THEN
      !Only calculate the greens function when needed
      IF(ANY(greensFunction(:)%l_calc)) THEN
         !calculate all contours only once
         DO iContour = 1, gfinp%numberContours
            CALL contour(iContour)%init(gfinp%contour(iContour))
            CALL contour(iContour)%eContour(gfinp%contour(iContour),results%ef,fmpi%irank)
         ENDDO
      ENDIF
      DO i_gf = 1, gfinp%n
         IF(.NOT.greensFunction(i_gf)%l_calc) CYCLE
         iContour = gfinp%elem(i_gf)%iContour
         greensFunction(i_gf)%contour = contour(iContour)
         CALL greensFunction(i_gf)%reset()
      ENDDO
      CALL greensfImagPart%init(gfinp,atoms,input,noco,ANY(greensFunction(:)%l_calc))
      IF(atoms%n_hia.GT.0 .AND. fmpi%irank==0 .AND.PRESENT(hub1data)) hub1data%mag_mom = 0.0
   ENDIF


   IF (fmpi%irank == 0) CALL openXMLElementNoAttributes('valenceDensity')

   !In a non-collinear calcuation where the off-diagonal part of the
   !density matrix in the muffin-tins is calculated, the a- and
   !b-coef. for both spins are needed at once. Thus, cdnval is only
   !called once and both spin directions are calculated in a single run.
   results%force=0.0
   DO jspin = 1,merge(1,input%jspins,noco%l_mperp.OR.banddos%l_jDOS)
      CALL cdnvalJob%init(fmpi,input,kpts,noco,results,jspin)
      IF (sliceplot%slice) CALL cdnvalJob%select_slice(sliceplot,results,input,kpts,noco,jspin)
      CALL cdnval(eig_id,fmpi,kpts,jspin,noco,nococonv,input,banddos,cell,atoms,enpara,stars,vacuum,&
                  sphhar,sym,vTot,oneD,cdnvalJob,outDen,regCharges,dos,vacdos,results,moments,gfinp,&
                  hub1inp,hub1data,coreSpecInput,mcd,slab,orbcomp,jDOS,greensfImagPart)
   END DO

   call val_den%copyPotDen(outDen)
   ! calculate kinetic energy density for MetaGGAs
   if(xcpot%exc_is_metagga()) then
      CALL calc_EnergyDen(eig_id, fmpi, kpts, noco, nococonv,input, banddos, cell, atoms, enpara, stars,&
                             vacuum,  sphhar, sym, gfinp, hub1inp, vTot, oneD, results, EnergyDen)
   endif

   IF (banddos%dos.or.banddos%band.or.input%cdinf) THEN
     IF (fmpi%irank == 0) THEN
       CALL timestart("cdngen: dos")
       call make_dos(kpts,atoms,vacuum,input,banddos,&
                      sliceplot,noco,sym,cell,results,eigdos,oneD)
       CALL timestop("cdngen: dos")
     END IF
     CALL juDFT_end("DOS OK",fmpi%irank)
   END IF

   IF ((banddos%dos.OR.banddos%vacdos)) CALL juDFT_end("DOS OK",fmpi%irank)

   CALL cdntot(stars,atoms,sym,vacuum,input,cell,oneD,outDen,.TRUE.,qtot,dummy,fmpi,.TRUE.)
   IF (fmpi%irank.EQ.0) THEN
      CALL closeXMLElement('valenceDensity')
   END IF ! fmpi%irank = 0

   IF (sliceplot%slice) THEN
      IF (fmpi%irank == 0) THEN
         IF(noco%l_alignMT) CALL juDFT_error("Relaxation of SQA and sliceplot not implemented. To perfom a sliceplot of the correct cdn deactivate realaxation.", calledby = "cdngen" )
         CALL writeDensity(stars,noco,vacuum,atoms,cell,sphhar,input,sym,oneD,CDN_ARCHIVE_TYPE_CDN_const,CDN_INPUT_DEN_const,&
                           0,-1.0,0.0,.FALSE.,outDen,'cdn_slice')
      END IF
#ifdef CPP_MPI
               CALL mpi_bc_potden(fmpi,stars,sphhar,atoms,input,vacuum,oneD,noco,outDen)
#endif
      CALL juDFT_end("slice OK",fmpi%irank)
   END IF

   !IF (sliceplot%iplot.NE.0) THEN
   !   CALL makeplots(stars, atoms, sphhar, vacuum, input, fmpi,oneD, sym, cell, noco,nococonv, outDen, PLOT_OUTDEN_Y_CORE, sliceplot)
   !END IF

   IF(PRESENT(greensFunction) .AND.gfinp%n.GT.0) THEN
      IF(greensfImagPart%l_calc) THEN
         CALL greensfPostProcess(greensFunction,greensfImagPart,atoms,gfinp,input,sym,noco,fmpi,&
                                 nococonv,vTot,enpara,hub1inp,sphhar,hub1data,results)
      ELSE
         IF(fmpi%irank.EQ.0) THEN
            WRITE(oUnit,'(/,A)') "Green's Functions are not calculated: "
            WRITE(oUnit,'(A,f12.7,TR5,A,f12.7/)') "lastDistance: ", results%last_distance,&
                                                  "minCalcDistance: ", gfinp%minCalcDistance
         ENDIF
      ENDIF
   ENDIF


   CALL timestart("cdngen: cdncore")
   if(xcpot%exc_is_MetaGGA()) then
      CALL cdncore(fmpi,oneD,input,vacuum,noco,nococonv,sym,&
                   stars,cell,sphhar,atoms,vTot,outDen,moments,results, EnergyDen)
   else
      CALL cdncore(fmpi,oneD,input,vacuum,noco,nococonv,sym,&
                   stars,cell,sphhar,atoms,vTot,outDen,moments,results)
   endif
   call core_den%subPotDen(outDen, val_den)
   CALL timestop("cdngen: cdncore")

   IF(.FALSE.) CALL denMultipoleExp(input, fmpi, atoms, sphhar, stars, sym, cell, oneD, outDen) ! There should be a switch in the inp file for this
   IF(fmpi%irank.EQ.0) THEN
      IF(input%lResMax>0) CALL resMoms(sym,input,atoms,sphhar,noco,nococonv,outDen,moments%rhoLRes) ! There should be a switch in the inp file for this
   END IF

   CALL enpara%calcOutParams(input,atoms,vacuum,regCharges)

   IF (fmpi%irank == 0) CALL openXMLElementNoAttributes('allElectronCharges')
   CALL qfix(fmpi,stars,atoms,sym,vacuum,sphhar,input,cell,oneD,outDen,noco%l_noco,.TRUE.,l_par=.TRUE.,force_fix=.TRUE.,fix=fix)
   IF (fmpi%irank == 0) THEN
      CALL closeXMLElement('allElectronCharges')

      IF (input%jspins == 2) THEN
         !Calculate and write out spin densities at the nucleus and magnetic moments in the spheres
         CALL magMoms(input,atoms,noco,nococonv,vTot,moments)
         if (sym%nop==1.and..not.input%film) call magMultipoles(sym,stars, atoms,cell, sphhar, vacuum, input, noco,nococonv,outden)
         !Generate and save the new nocoinp file if the directions of the local
         !moments are relaxed or a constraint B-field is calculated.
         IF (ANY(noco%l_relax(:atoms%ntype)).OR.noco%l_constr) THEN
          !  CALL genNewNocoInp(input,atoms,noco,noco_new)
         END IF
         IF (noco%l_soc) CALL orbMagMoms(input,atoms,noco,nococonv,moments%clmom)
      END IF
   END IF ! fmpi%irank == 0
   Perform_metagga = Allocated(Energyden%Mt) &
                   .And. (Xcpot%Exc_is_metagga() .Or. Xcpot%Vx_is_metagga())
   If(Perform_metagga) Then
     IF(noco%l_alignMT) CALL juDFT_error("Relaxation of SQA and metagga not implemented.", calledby = "cdngen" )
     CALL writeDensity(stars,noco,vacuum,atoms,cell,sphhar,input,sym,oneD,CDN_ARCHIVE_TYPE_CDN_const,CDN_INPUT_DEN_const,&
                           0,-1.0,0.0,.FALSE.,core_den,'cdnc')
   endif

#ifdef CPP_MPI
   CALL MPI_BCAST(nococonv%alph,atoms%ntype,MPI_DOUBLE_PRECISION,0,fmpi%mpi_comm,ierr)
   CALL MPI_BCAST(nococonv%beta,atoms%ntype,MPI_DOUBLE_PRECISION,0,fmpi%mpi_comm,ierr)
   CALL MPI_BCAST(nococonv%b_con,atoms%ntype*2,MPI_DOUBLE_PRECISION,0,fmpi%mpi_comm,ierr)
   CALL MPI_BCAST(nococonv%qss,3,MPI_DOUBLE_PRECISION,0,fmpi%mpi_comm,ierr)

   CALL mpi_bc_potden(fmpi,stars,sphhar,atoms,input,vacuum,oneD,noco,outDen)
#endif

END SUBROUTINE cdngen

END MODULE m_cdngen
