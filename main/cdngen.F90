!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_cdngen
#ifdef CPP_MPI
   USE mpi
#endif
CONTAINS

SUBROUTINE cdngen(eig_id,fmpi,input,banddos,sliceplot,vacuum,&
                  kpts,atoms,sphhar,stars,sym,gfinp,hub1inp,&
                  enpara,cell,noco,nococonv,vTot,results ,coreSpecInput,&
                  archiveType, xcpot,outDen,EnergyDen,greensFunction,hub1data,vxc,exc)

   !*****************************************************
   !    Charge density generator
   !    calls cdnval to generate the valence charge and the
   !    core routines for the core contribution
   !*****************************************************
   use m_types_vacdos
   use m_types_mcd
   use m_types_slab
   use m_types_orbcomp
   use m_types_jdos
   USE m_types
   USE m_constants
   USE m_juDFT
   !USE m_prpqfftmap
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
   use m_slater
   USE m_greensfPostProcess
   USE m_types_greensfContourData
   USE m_types_eigdos
   USE m_types_dos

   USE m_force_sf ! Klueppelberg (force level 3)

   IMPLICIT NONE

   ! Type instance arguments
   TYPE(t_results),INTENT(INOUT)    :: results
   TYPE(t_mpi),INTENT(IN)           :: fmpi

    
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
   TYPE(t_potden),INTENT(INOUT),OPTIONAL:: vxc, exc

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
   INTEGER               :: i_gf,iContour,n

   TYPE(t_eigdos_list),allocatable :: eigdos(:)

#ifdef CPP_HDF
   INTEGER(HID_T)        :: banddosFile_id
#endif
   LOGICAL               :: l_error,Perform_metagga

   ! Initialization section
   CALL regCharges%init(input,atoms)
   CALL moments%init(fmpi,input,sphhar,atoms)
   !initalize data for DOS
   if (noco%l_noco) results%eig(:,:,2)=results%eig(:,:,1)
   CALL dos%init(input,atoms,kpts,banddos,results%eig)
   CALL vacdos%init(input,atoms,kpts,banddos,results%eig)
   CALL mcd%init(banddos,input,atoms,kpts,results%eig)
   CALL slab%init(banddos,atoms,cell,input,kpts)
   CALL orbcomp%init(input,banddos,atoms,kpts,results%eig)
   CALL jDOS%init(input,banddos,atoms,kpts,results%eig)

   if (banddos%dos.or.banddos%band) then
     allocate(eigdos(count((/banddos%dos.or.banddos%band,banddos%vacdos,banddos%l_mcd,banddos%l_slab,banddos%l_orb,banddos%l_jDOS/))))
     n=2
     eigdos(1)%p=>dos
     if (banddos%vacdos) THEN
       eigdos(n)%p=>vacdos; n=n+1;
     endif
     if (banddos%l_mcd) THEN
       eigdos(n)%p=>mcd; n=n+1
     endif
     if (banddos%l_slab) THEN
       eigdos(n)%p=>slab; n=n+1
     endif
     if (banddos%l_orb) THEN
       eigdos(n)%p=>orbcomp; n=n+1
     endif
     if (banddos%l_jdos) eigdos(n)%p=>jDOS
   endif



   CALL outDen%init(stars,    atoms, sphhar, vacuum, noco, input%jspins, POTDEN_TYPE_DEN)
   CALL EnergyDen%init(stars, atoms, sphhar, vacuum, noco, input%jspins, POTDEN_TYPE_EnergyDen)


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
      CALL greensfImagPart%init(gfinp,atoms,input,noco,ANY(greensFunction(:)%l_calc),SIZE(fmpi%k_list))
   ENDIF

   IF(fmpi%irank==0 .AND.PRESENT(hub1data)) THEN
      hub1data%mag_mom = 0.0
      hub1data%cdn_atomic = 0.0
   ENDIF


   IF (fmpi%irank == 0) CALL openXMLElementNoAttributes('valenceDensity')

   !In a non-collinear calcuation where the off-diagonal part of the
   !density matrix in the muffin-tins is calculated, the a- and
   !b-coef. for both spins are needed at once. Thus, cdnval is only
   !called once and both spin directions are calculated in a single run.
   CALL timestart("cdngen: cdnval")
   DO jspin = 1,merge(1,input%jspins,noco%l_mperp.OR.banddos%l_jDOS)
      CALL cdnvalJob%init(fmpi,input,kpts,noco,results,jspin)
      IF (sliceplot%slice) CALL cdnvalJob%select_slice(sliceplot,results,input,kpts,noco,jspin)
      CALL cdnval(eig_id,fmpi,kpts,jspin,noco,nococonv,input,banddos,cell,atoms,enpara,stars,vacuum,&
                  sphhar,sym,vTot ,cdnvalJob,outDen,regCharges,dos,vacdos,results,moments,gfinp,&
                  hub1inp,hub1data,coreSpecInput,mcd,slab,orbcomp,jDOS,greensfImagPart)
   END DO
   CALL timestop("cdngen: cdnval")

   call val_den%copyPotDen(outDen)
   ! calculate kinetic energy density for MetaGGAs
   if(xcpot%exc_is_metagga()) then
      CALL calc_EnergyDen(eig_id, fmpi, kpts, noco, nococonv,input, banddos, cell, atoms, enpara, stars,&
                             vacuum,  sphhar, sym, gfinp, hub1inp, vTot,   results, EnergyDen)
   endif

   IF (banddos%dos.or.banddos%band.or.input%cdinf) THEN
      IF (fmpi%irank == 0) THEN
         CALL timestart("cdngen: dos")
         CALL make_dos(kpts,atoms,vacuum,input,banddos,&
                      sliceplot,noco,sym,cell,results,eigdos )
         CALL timestop("cdngen: dos")
      END IF
   END IF

   CALL cdntot(stars,atoms,sym,vacuum,input,cell ,outDen,.TRUE.,qtot,dummy,fmpi,.TRUE.)
   IF (fmpi%irank.EQ.0) THEN
      CALL closeXMLElement('valenceDensity')
   END IF ! fmpi%irank = 0

   IF(PRESENT(greensFunction) .AND.gfinp%n.GT.0) THEN
      IF(greensfImagPart%l_calc) THEN
         CALL greensfPostProcess(greensFunction,greensfImagPart,atoms,kpts,cell,gfinp,input,sym,noco,fmpi,&
                                 nococonv,vTot,enpara,hub1inp,sphhar,hub1data,results)
      ELSE
         IF(fmpi%irank.EQ.0) THEN
            WRITE(oUnit,'(/,A)') "Green's Functions are not calculated: "
            WRITE(oUnit,'(A,f12.7,TR5,A,f12.7/)') "lastDistance: ", results%last_distance,&
                                                  "minCalcDistance: ", gfinp%minCalcDistance
         ENDIF
      ENDIF
   ENDIF

   IF (banddos%vacdos.or.banddos%dos.or.banddos%band.or.input%cdinf) THEN
      CALL juDFT_end("Charge density postprocessing done.",fmpi%irank)
   END IF

   IF (sliceplot%slice) THEN
      IF (fmpi%irank == 0) THEN
         IF(any(noco%l_alignMT)) CALL juDFT_error("Relaxation of SQA and sliceplot not implemented. To perfom a sliceplot of the correct cdn deactivate realaxation.", calledby = "cdngen" )
         CALL writeDensity(stars,noco,vacuum,atoms,cell,sphhar,input,sym ,CDN_ARCHIVE_TYPE_CDN_const,CDN_INPUT_DEN_const,&
                           0,-1.0,0.0,-1.0,-1.0,.FALSE.,outDen,inFilename='cdn_slice')
      END IF
      call outDen%distribute(fmpi%mpi_comm)
      CALL juDFT_end("slice OK",fmpi%irank)
   END IF

   !IF (sliceplot%iplot.NE.0) THEN
   !   CALL makeplots(stars, atoms, sphhar, vacuum, input, fmpi , sym, cell, noco,nococonv, outDen, PLOT_OUTDEN_Y_CORE, sliceplot)
   !END IF

   CALL timestart("cdngen: cdncore")
   if(xcpot%exc_is_MetaGGA()) then
      CALL cdncore(fmpi ,input,vacuum,noco,nococonv,sym,&
                   stars,cell,sphhar,atoms,vTot,outDen,moments,results, EnergyDen)
   else
      CALL cdncore(fmpi ,input,vacuum,noco,nococonv,sym,&
                   stars,cell,sphhar,atoms,vTot,outDen,moments,results)
   endif
   call core_den%subPotDen(outDen, val_den)
   CALL timestop("cdngen: cdncore")

   IF(.FALSE.) CALL denMultipoleExp(input, fmpi, atoms, sphhar, stars, sym, cell,   outDen) ! There should be a switch in the inp file for this
   IF(fmpi%irank.EQ.0) THEN
      IF(input%lResMax>0) CALL resMoms(sym,input,atoms,sphhar,noco,nococonv,outDen,moments%rhoLRes) ! There should be a switch in the inp file for this
   END IF

   IF(atoms%n_opc>0) THEN
      do jspin=1, input%jspins
         call slater(input,jspin,atoms,vTot%mt(:,0,:,jspin),l_write=fmpi%irank==0)
      enddo
   ENDIF

   CALL enpara%calcOutParams(input,atoms,vacuum,regCharges)

   IF (fmpi%irank == 0) CALL openXMLElementNoAttributes('allElectronCharges')
   CALL qfix(fmpi,stars,atoms,sym,vacuum,sphhar,input,cell ,outDen,noco%l_noco,.TRUE.,l_par=.TRUE.,force_fix=.TRUE.,fix=fix)
   IF (fmpi%irank == 0) THEN
      CALL closeXMLElement('allElectronCharges')

      IF (input%jspins == 2) THEN
         !Calculate and write out spin densities at the nucleus and magnetic moments in the spheres
         CALL magMoms(input,atoms,noco,nococonv,vTot,moments)
         if (sym%nop==1.and..not.input%film) call magMultipoles(sym,stars, atoms,cell, sphhar, vacuum, input, noco,nococonv,outden)
         !Generate and save the new nocoinp file if the directions of the local
         !moments are relaxed or a constraint B-field is calculated.
         IF (ANY(noco%l_alignMT).OR.any(noco%l_constrained)) THEN
          !  CALL genNewNocoInp(input,atoms,noco,noco_new)
         END IF
         IF (noco%l_soc) CALL orbMagMoms(input,atoms,noco,nococonv,moments%clmom)
      END IF
   END IF ! fmpi%irank == 0
   Perform_metagga = Allocated(Energyden%Mt) &
                   .And. (Xcpot%Exc_is_metagga() .Or. Xcpot%Vx_is_metagga())
   If(Perform_metagga) Then
     IF(any(noco%l_alignMT)) CALL juDFT_error("Relaxation of SQA and metagga not implemented.", calledby = "cdngen" )
     CALL writeDensity(stars,noco,vacuum,atoms,cell,sphhar,input,sym ,CDN_ARCHIVE_TYPE_CDN_const,CDN_INPUT_DEN_const,&
                           0,-1.0,0.0,-1.0,-1.0,.FALSE.,core_den,inFilename='cdnc')
   endif

#ifdef CPP_MPI
   CALL MPI_BCAST(nococonv%alph,atoms%ntype,MPI_DOUBLE_PRECISION,0,fmpi%mpi_comm,ierr)
   CALL MPI_BCAST(nococonv%beta,atoms%ntype,MPI_DOUBLE_PRECISION,0,fmpi%mpi_comm,ierr)
   CALL MPI_BCAST(nococonv%b_con,atoms%ntype*2,MPI_DOUBLE_PRECISION,0,fmpi%mpi_comm,ierr)
   CALL MPI_BCAST(nococonv%qss,3,MPI_DOUBLE_PRECISION,0,fmpi%mpi_comm,ierr)
#endif
   CALL outDen%distribute(fmpi%mpi_comm)

   ! Klueppelberg (force level 3)
   IF (input%l_f.AND.(input%f_level.GE.3).AND.(fmpi%irank.EQ.0)) THEN
      DO jspin = 1,input%jspins ! jsp_start, jsp_end
         CALL force_sf_mt(atoms,sphhar,jspin,jspin,fmpi,vtot%mt(:,0:,:,jspin),exc%mt(:,0:,:,1),vxc%mt(:,0:,:,:),outDen%mt(:,0:,:,:),sym,cell)
      END DO
   END IF

END SUBROUTINE cdngen

END MODULE m_cdngen
