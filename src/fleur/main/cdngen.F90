!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_cdngen
#ifdef CPP_MPI
   USE mpi
#endif
   implicit none
CONTAINS

SUBROUTINE cdngen(eig_id,fmpi,input,banddos,sliceplot,vacuum,&
                  kpts,atoms,sphhar,stars,sym,juphon,gfinp,hub1inp,&
                  enpara,cell,field,noco,nococonv,vTot,results ,coreSpecInput,&
                  archiveType, xcpot,outDen,EnergyDen,core_den,greensFunction,hub1data,vxc,exc)

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
   use m_types_jointdos
   USE m_types
   USE m_constants
   USE m_juDFT
   USE m_cdnval
   USE m_plot
   USE m_cdn_io
   USE m_wrtdop
   USE m_cdntot
   USE m_qfix
   USE m_xmlOutput
   USE m_magMultipoles
   USE m_magmoments
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
   USE m_types_hyperfine

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
   TYPE(t_juphon),INTENT(IN)        :: juphon
   TYPE(t_stars),INTENT(IN)         :: stars
   TYPE(t_cell),INTENT(IN)          :: cell
   TYPE(t_field),INTENT(IN)         :: field
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
   TYPE(t_potden),INTENT(OUT),optional       :: core_den
   TYPE(t_potden),INTENT(INOUT),OPTIONAL:: vxc, exc

   !Scalar Arguments
   INTEGER, INTENT (IN)             :: eig_id, archiveType

   ! Local type instances
   TYPE(t_dos),TARGET             :: dos
   TYPE(t_vacdos),TARGET          :: vacdos
   TYPE(t_moments)                :: moments
   TYPE(t_mcd),TARGET             :: mcd
   TYPE(t_slab),TARGET            :: slab
   TYPE(t_orbcomp),TARGET         :: orbcomp
   TYPE(t_jDOS),TARGET            :: jDOS
   TYPE(t_jointDOS),TARGET       :: jointDOS
   TYPE(t_cdnvalJob)       :: cdnvalJob
   TYPE(t_greensfImagPart) :: greensfImagPart
   TYPE(t_potden)          :: val_den
   TYPE(t_greensfContourData) :: contour(gfinp%numberContours)
   TYPE(t_hyperfine)       :: hyperfine


   !Local Scalars
   REAL                  :: fix, qtot, dummy, eFermiPrev
   INTEGER               :: jspin, ierr
   INTEGER               :: dim_idx
   INTEGER               :: i_gf,iContour,n
   TYPE(t_eigdos_list),allocatable :: eigdos(:)

#ifdef CPP_HDF
   INTEGER(HID_T)        :: banddosFile_id
#endif
   LOGICAL               :: l_error,Perform_metagga

   ! Initialization section
   CALL moments%init(fmpi,input,sphhar,atoms)
   !initalize data for DOS
   if (noco%l_noco) results%eig(:,:,2)=results%eig(:,:,1)
   
   if (banddos%dos.or.banddos%band.or.input%cdinf) then
     CALL initialize_eigdos_types(eigdos, dos, jointDOS, vacdos, mcd, slab, orbcomp, jDOS, &
                                   input, atoms, kpts, banddos, noco, results, cell)
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

   CALL hyperfine%init(input, atoms)

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
                  sphhar,sym,vTot ,cdnvalJob,outDen,dos,vacdos,results,moments,gfinp,&
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
                      sliceplot,noco,nococonv,sym,cell,results,eigdos )
         CALL timestop("cdngen: dos")
      END IF
   END IF

   CALL cdntot(stars,nococonv,atoms,sym,vacuum,input,cell ,outDen,.TRUE.,qtot,dummy,fmpi,.TRUE.)
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
         CALL writeDensity(stars,noco,vacuum,atoms,cell,sphhar,input,sym,archiveType,CDN_INPUT_DEN_const,&
                           1,-1.0,0.0,-1.0,-1.0,.FALSE.,outDen,inFilename='cdn_slice')
      END IF
      call outDen%distribute(fmpi%mpi_comm)
      CALL juDFT_end("slice OK",fmpi%irank)
   END IF

   !IF (sliceplot%iplot.NE.0) THEN
   !   CALL makeplots(stars, atoms, sphhar, vacuum, input, fmpi , sym, cell, noco,nococonv, outDen, PLOT_OUTDEN_Y_CORE, sliceplot)
   !END IF

   CALL hyperfine%printValenceHyperfine(input, atoms, fmpi, moments)

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

   CALL outDen%distribute(fmpi%mpi_comm)

   IF(.FALSE.) CALL denMultipoleExp(input, fmpi, atoms, sphhar, stars, sym, cell,   outDen) ! There should be a switch in the inp file for this
   IF(fmpi%irank.EQ.0) THEN
      IF(input%lResMax>0) CALL resMoms(sym,input,atoms,sphhar,noco,nococonv,outDen,moments%rhoLRes) ! There should be a switch in the inp file for this
   END IF

   IF(atoms%n_opc>0) THEN
      do jspin=1, input%jspins
         call slater(input,jspin,atoms,vTot%mt(:,0,:,jspin),l_write=fmpi%irank==0)
      enddo
   ENDIF

  
   IF (fmpi%irank == 0) CALL openXMLElementNoAttributes('allElectronCharges')
   CALL qfix(fmpi,stars,nococonv,atoms,sym,vacuum,sphhar,input,cell,field,outDen,noco%l_noco,.TRUE.,.TRUE.,.TRUE.,fix)
   IF (fmpi%irank == 0) CALL closeXMLElement('allElectronCharges')
   IF (input%jspins == 2) THEN
      !Calculate and write out spin densities at the nucleus and magnetic moments in the spheres
      IF (fmpi%irank == 0) THEN
         CALL spinMoments(input,atoms,noco,nococonv,den=outDen,results=results)
         CALL orbMoments(input,atoms,noco,nococonv,moments)
         CALL write_output_struct_xsf(atoms,nococonv,outDen)
         if (any(noco%l_constrained).or.any(noco%l_fixedMoment)) call nococonv%update_b_cons(atoms,noco,vtot,outDen)
      END IF

      if (sym%nop==1.and..not.input%film) call magMultipoles(fmpi,sym,stars, atoms,cell, sphhar, vacuum, input, noco,nococonv,outden)
      !Generate and save the new nocoinp file if the directions of the local
      !moments are relaxed or a constraint B-field is calculated.
   END IF

   CALL hyperfine%calcPrintIsomerShifts(input,atoms,fmpi,outDen)

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
   CALL MPI_BCAST(nococonv%b_con,atoms%ntype*3,MPI_DOUBLE_PRECISION,0,fmpi%mpi_comm,ierr)
   CALL MPI_BCAST(nococonv%qss,3,MPI_DOUBLE_PRECISION,0,fmpi%mpi_comm,ierr)
#endif


   ! Klueppelberg (force level 3)
   IF (input%l_f.AND.(input%f_level.GE.3).AND.(fmpi%irank.EQ.0)) THEN
      DO jspin = 1,input%jspins ! jsp_start, jsp_end
         CALL force_sf_mt(atoms,sphhar,jspin,jspin,fmpi,vtot%mt(:,0:,:,jspin),exc%mt(:,0:,:,1),vxc%mt(:,0:,:,:),outDen%mt(:,0:,:,:),sym,cell)
      END DO
   END IF

END SUBROUTINE cdngen

SUBROUTINE write_output_struct_xsf(atoms,nococonv,outDen)

   USE m_types

   IMPLICIT NONE

   TYPE(t_atoms),INTENT(IN)      :: atoms
   TYPE(t_nococonv),INTENT(INOUT):: nococonv
   TYPE(t_potden),INTENT(IN)     :: outDen

   INTEGER, PARAMETER            :: inUnit = 97, outUnitLocal = 98
   LOGICAL                       :: l_exists
   INTEGER                       :: ios, iAtom, iType, atomCount
   INTEGER                       :: zatom
   REAL                          :: x, y, z
   REAL, ALLOCATABLE             :: mag_mom_xsf(:,:), magm_type(:,:), theta(:), phi(:)
   CHARACTER(len=1024)           :: line

   INQUIRE(file='struct.xsf', exist=l_exists)
   IF (.NOT.l_exists) RETURN

   ALLOCATE(mag_mom_xsf(3,atoms%nat), magm_type(3,atoms%ntype), theta(atoms%ntype), phi(atoms%ntype))
   CALL nococonv%avg_moments(outDen,atoms,magm_type,theta,phi)

   mag_mom_xsf = 0.0
   DO iType = 1, atoms%ntype
      DO iAtom = atoms%firstAtom(iType), atoms%firstAtom(iType) + atoms%neq(iType) - 1
         mag_mom_xsf(:,iAtom) = magm_type(:,iType)
      END DO
   END DO
   OPEN(inUnit, file='struct.xsf', status='old', action='read', iostat=ios)
   IF (ios /= 0) THEN
      DEALLOCATE(mag_mom_xsf, magm_type, theta, phi)
      RETURN
   END IF

   OPEN(outUnitLocal, file='output-struct.xsf', status='replace', action='write', iostat=ios)
   IF (ios /= 0) THEN
      CLOSE(inUnit)
      DEALLOCATE(mag_mom_xsf, magm_type, theta, phi)
      RETURN
   END IF

   DO
      READ(inUnit,'(A)',iostat=ios) line
      IF (ios /= 0) EXIT

      WRITE(outUnitLocal,'(A)') TRIM(line)
      IF (TRIM(ADJUSTL(line)) /= 'PRIMCOORD') CYCLE

      READ(inUnit,'(A)',iostat=ios) line
      IF (ios /= 0) EXIT
      WRITE(outUnitLocal,'(A)') TRIM(line)

      atomCount = atoms%nat
      DO iAtom = 1, atomCount
         READ(inUnit,'(A)',iostat=ios) line
         IF (ios /= 0) EXIT

         READ(line,*,iostat=ios) zatom, x, y, z
         IF (ios /= 0) THEN
            WRITE(outUnitLocal,'(A)') TRIM(line)
         ELSE
            WRITE(outUnitLocal,'(i4,2x,6(f0.7,1x))') zatom, x, y, z, mag_mom_xsf(:,iAtom)
         END IF
      END DO
   END DO

   CLOSE(inUnit)
   CLOSE(outUnitLocal)
   DEALLOCATE(mag_mom_xsf, magm_type, theta, phi)

END SUBROUTINE write_output_struct_xsf

SUBROUTINE initialize_eigdos_types(eigdos, dos, jointDOS, vacdos, mcd, slab, orbcomp, jDOS, &
                                    input, atoms, kpts, banddos, noco, results, cell)
   !*****************************************************
   ! Initialize all eigenvalue/DOS types and populate
   ! the eigdos pointer array
   !*****************************************************
   USE m_types_eigdos
   USE m_types_dos
   USE m_types_jointdos
   USE m_types_vacdos
   USE m_types_mcd
   USE m_types_slab
   USE m_types_orbcomp
   USE m_types_jdos
   use m_types
   
   IMPLICIT NONE
   
   ! Arguments
   TYPE(t_eigdos_list), ALLOCATABLE, INTENT(INOUT) :: eigdos(:)
   TYPE(t_dos), TARGET, INTENT(INOUT)              :: dos
   TYPE(t_jointDOS), TARGET, INTENT(INOUT)         :: jointDOS
   TYPE(t_vacdos), TARGET, INTENT(INOUT)           :: vacdos
   TYPE(t_mcd), TARGET, INTENT(INOUT)              :: mcd
   TYPE(t_slab), TARGET, INTENT(INOUT)             :: slab
   TYPE(t_orbcomp), TARGET, INTENT(INOUT)          :: orbcomp
   TYPE(t_jDOS), TARGET, INTENT(INOUT)             :: jDOS
   TYPE(t_input), INTENT(IN)                       :: input
   TYPE(t_atoms), INTENT(IN)                       :: atoms
   TYPE(t_kpts), INTENT(IN)                        :: kpts
   TYPE(t_banddos), INTENT(IN)                     :: banddos
   TYPE(t_noco), INTENT(IN)                        :: noco
   TYPE(t_results), INTENT(IN)                     :: results
   TYPE(t_cell), INTENT(IN)                        :: cell
   
   ! Local variables
   INTEGER :: n, num_types
   LOGICAL :: type_flags(7)
   
   ! Determine which types need to be initialized
   type_flags(1) = banddos%dos .OR. banddos%band .OR. input%cdinf
   type_flags(2) = banddos%l_jointDOS
   type_flags(3) = banddos%vacdos
   type_flags(4) = banddos%l_mcd
   type_flags(5) = banddos%l_slab
   type_flags(6) = banddos%l_orb
   type_flags(7) = banddos%l_jDOS
   
   ! Count number of types to allocate
   num_types = COUNT(type_flags)
   ALLOCATE(eigdos(num_types))

   
   
   ! Initialize DOS (always first)
   n = 1
   CALL dos%init(input, atoms, kpts, banddos, noco%l_noco .OR. banddos%l_jDOS, results%eig)
   eigdos(1)%p => dos
   n = 2
   
   ! Initialize optional types
   IF (banddos%l_jointDOS) THEN
      CALL jointDOS%init(input, atoms, kpts, banddos, noco%l_noco, results%eig)
      eigdos(n)%p => jointDOS
      n = n + 1
   END IF
   
   CALL vacdos%init(input, atoms, kpts, banddos, results%eig)
   IF (banddos%vacdos) THEN
      eigdos(n)%p => vacdos
      n = n + 1
   END IF
   
   IF (banddos%l_mcd) THEN
      CALL mcd%init(banddos, input, atoms, kpts, results%eig)
      eigdos(n)%p => mcd
      n = n + 1
   END IF
   
   IF (banddos%l_slab) THEN
      CALL slab%init(banddos, atoms, cell, input, kpts)
      eigdos(n)%p => slab
      n = n + 1
   END IF
   
   IF (banddos%l_orb) THEN
      CALL orbcomp%init(input, banddos, atoms, kpts, results%eig)
      eigdos(n)%p => orbcomp
      n = n + 1
   END IF
   
   IF (banddos%l_jdos) THEN
      CALL jDOS%init(input, banddos, atoms, kpts, results%eig)
      eigdos(n)%p => jDOS
   END IF
   
END SUBROUTINE initialize_eigdos_types

END MODULE m_cdngen
