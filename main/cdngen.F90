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
                  archiveType,outDen)

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
   USE m_cdnovlp
   USE m_qfix
   USE m_genNewNocoInp
   USE m_xmlOutput
   USE m_magMoms
   USE m_orbMagMoms
   USE m_cdncore
   USE m_doswrite
   USE m_Ekwritesl
   USE m_banddos_io
   USE m_unfold_band_kpts
   USE m_gOnsite
   USE m_gOnsite_radial !to be unified with m_gOnsite
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
   TYPE(t_potden),INTENT(INOUT)     :: outDen

   !Scalar Arguments
   INTEGER, INTENT (IN)             :: eig_id, archiveType

   ! Local type instances
   TYPE(t_noco)          :: noco_new
   TYPE(t_regionCharges) :: regCharges
   TYPE(t_dos)           :: dos
   TYPE(t_moments)       :: moments
   TYPE(t_mcd)           :: mcd
   TYPE(t_slab)          :: slab
   TYPE(t_orbcomp)       :: orbcomp
   TYPE(t_cdnvalJob)     :: cdnvalJob
   TYPE(t_greensf)       :: gOnsite


   !Local Scalars
   REAL                  :: fix, qtot, dummy,eFermiPrev
   INTEGER               :: jspin, jspmax, ierr
#ifdef CPP_HDF
   INTEGER(HID_T)        :: banddosFile_id
#endif
   LOGICAL               :: l_error

   CALL regCharges%init(input,atoms)
   CALL dos%init(input,atoms,dimension,kpts,vacuum)
   CALL moments%init(input,atoms)
   CALL mcd%init1(banddos,dimension,input,atoms,kpts)
   CALL slab%init(banddos,dimension,atoms,cell,input,kpts)
   CALL orbcomp%init(input,banddos,dimension,atoms,kpts)
   CALL gOnsite%init(input,atoms,kpts,dimension)

   IF (mpi%irank.EQ.0) CALL openXMLElementNoAttributes('valenceDensity')

   !In a non-collinear calcuation where the off-diagonal part of the
   !density matrix in the muffin-tins is calculated, the a- and
   !b-coef. for both spins are needed at once. Thus, cdnval is only
   !called once and both spin directions are calculated in a single run.
   jspmax = input%jspins
   IF (noco%l_mperp) jspmax = 1
   DO jspin = 1,jspmax
      CALL cdnvalJob%init(mpi,input,kpts,noco,results,jspin,sliceplot,banddos)
      CALL cdnval(eig_id,mpi,kpts,jspin,noco,input,banddos,cell,atoms,enpara,stars,vacuum,dimension,&
                  sphhar,sym,vTot,oneD,cdnvalJob,outDen,regCharges,dos,results,moments,coreSpecInput,mcd,slab,orbcomp,gOnsite)
   END DO

    IF (atoms%n_hia.GT.0) THEN
      DO jspin = 1, jspmax
         IF(input%ldahia_sphavg) THEN
            CALL calc_onsite(atoms,jspin,input%jspins,dimension%neigd,kpts%ntet,kpts%nkpt,kpts%ntetra(1:4,:),kpts%voltet(:),&
                                                   results%neig(:,jspin),results%eig(:,:,jspin),gOnsite,results%ef,sym)
         ELSE
            CALL calc_onsite_radial(atoms,enpara,vTot%mt(:,0,:,:),jspin,input%jspins,dimension%neigd,kpts%ntet,kpts%nkpt,kpts%ntetra(1:4,:),kpts%voltet(:),&
                                                   results%neig(:,jspin),results%eig(:,:,jspin),gOnsite,results%ef,sym)
         ENDIF
      ENDDO
   ENDIF

   IF (mpi%irank.EQ.0) THEN
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
         IF (banddos%dos.AND.(banddos%ndir.EQ.-3)) THEN
            CALL Ek_write_sl(eig_id,dimension,kpts,atoms,vacuum,input,jspmax,sym,cell,dos,slab,orbcomp,results)
         END IF
         CALL timestop("cdngen: dos")
      END IF
   END IF

   IF ((banddos%dos.OR.banddos%vacdos).AND.(banddos%ndir/=-2)) CALL juDFT_end("DOS OK",mpi%irank)
   IF (vacuum%nstm.EQ.3) CALL juDFT_end("VACWAVE OK",mpi%irank)

   IF (mpi%irank.EQ.0) THEN
      CALL cdntot(stars,atoms,sym,vacuum,input,cell,oneD,outDen,.TRUE.,qtot,dummy)
      CALL closeXMLElement('valenceDensity')
   END IF ! mpi%irank = 0

   IF (sliceplot%slice) THEN
      IF (mpi%irank.EQ.0) THEN
         CALL writeDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,CDN_ARCHIVE_TYPE_CDN_const,CDN_INPUT_DEN_const,&
                           0,-1.0,0.0,.FALSE.,outDen,'cdn_slice')
      END IF
      CALL juDFT_end("slice OK",mpi%irank)
   END IF

   CALL timestart("cdngen: cdncore")
   CALL cdncore(mpi,dimension,oneD,input,vacuum,noco,sym,&
                stars,cell,sphhar,atoms,vTot,outDen,moments,results)
   CALL timestop("cdngen: cdncore")

   CALL enpara%calcOutParams(input,atoms,vacuum,regCharges)

   IF (mpi%irank.EQ.0) THEN
      CALL openXMLElementNoAttributes('allElectronCharges')
      CALL qfix(stars,atoms,sym,vacuum,sphhar,input,cell,oneD,outDen,noco%l_noco,.TRUE.,.true.,fix)
      CALL closeXMLElement('allElectronCharges')

      IF (input%jspins.EQ.2) THEN
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
   END IF ! mpi%irank.EQ.0

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
