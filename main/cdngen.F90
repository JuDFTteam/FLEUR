!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_cdngen

USE m_juDFT

CONTAINS

SUBROUTINE cdngen(eig_id,mpi,input,banddos,sliceplot,vacuum,&
                  dimension,kpts,atoms,sphhar,stars,sym,obsolete,&
                  enpara,cell,noco,vTot,results,oneD,coreSpecInput,&
                  archiveType,outDen)

   !*****************************************************
   !    Charge density generator
   !    calls cdnval to generate the valence charge and the
   !    core routines for the core contribution
   !*****************************************************

   USE m_constants
   USE m_prpqfftmap
   USE m_cdnval
   USE m_cdn_io
   USE m_wrtdop
   USE m_cdntot
   USE m_cdnovlp
   USE m_qfix
   USE m_genNewNocoInp
   USE m_types
   USE m_xmlOutput
   USE m_magMoms
   USE m_orbMagMoms
   USE m_cdncore
   USE m_doswrite
   USE m_Ekwritesl
#ifdef CPP_MPI
   USE m_mpi_bc_potden
#endif

   IMPLICIT NONE

   ! Type instance arguments
   TYPE(t_results),INTENT(INOUT)    :: results
   TYPE(t_mpi),INTENT(IN)           :: mpi
   TYPE(t_dimension),INTENT(IN)     :: dimension
   TYPE(t_oneD),INTENT(IN)          :: oneD
   TYPE(t_enpara),INTENT(INOUT)     :: enpara
   TYPE(t_obsolete),INTENT(IN)      :: obsolete
   TYPE(t_banddos),INTENT(IN)       :: banddos
   TYPE(t_sliceplot),INTENT(IN)     :: sliceplot
   TYPE(t_input),INTENT(IN)         :: input
   TYPE(t_vacuum),INTENT(IN)        :: vacuum
   TYPE(t_noco),INTENT(IN)          :: noco
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
   TYPE(t_cdnvalKLoop)   :: cdnvalKLoop


   !Local Scalars
   REAL                  :: fix, qtot, dummy
   INTEGER               :: jspin, jspmax

   CALL regCharges%init(input,atoms)
   CALL dos%init(input,atoms,dimension,kpts,vacuum)
   CALL moments%init(input,atoms)

   IF (mpi%irank.EQ.0) CALL openXMLElementNoAttributes('valenceDensity')

   !In a non-collinear calcuation where the off-diagonal part of the
   !density matrix in the muffin-tins is calculated, the a- and
   !b-coef. for both spins are needed at once. Thus, cdnval is only
   !called once and both spin directions are calculated in a single run.
   jspmax = input%jspins
   IF (noco%l_mperp) jspmax = 1
   DO jspin = 1,jspmax
      CALL cdnvalKLoop%init(mpi,input,kpts,banddos,noco,results,jspin,sliceplot)
      CALL cdnval(eig_id,mpi,kpts,jspin,sliceplot,noco,input,banddos,cell,atoms,enpara,stars,vacuum,dimension,&
                  sphhar,sym,obsolete,vTot,oneD,coreSpecInput,cdnvalKLoop,outDen,regCharges,dos,results,moments,mcd,slab)
   END DO

   IF (mpi%irank.EQ.0) THEN
      IF (banddos%dos.or.banddos%vacdos.or.input%cdinf) THEN
         CALL timestart("cdngen: dos")
         CALL doswrite(eig_id,dimension,kpts,atoms,vacuum,input,banddos,sliceplot,noco,sym,cell,dos,mcd,results,slab%nsld,oneD)
         IF (banddos%dos.AND.(banddos%ndir.EQ.-3)) THEN
            CALL Ek_write_sl(eig_id,dimension,kpts,atoms,vacuum,input,jspmax,sym,cell,dos,slab)
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

   CALL cdncore(results,mpi,dimension,oneD,input,vacuum,noco,sym,&
                stars,cell,sphhar,atoms,vTot,outDen,moments)

   CALL enpara%calcOutParams(input,atoms,vacuum,regCharges)

   IF (mpi%irank.EQ.0) THEN
      CALL openXMLElementNoAttributes('allElectronCharges')
      CALL qfix(stars,atoms,sym,vacuum,sphhar,input,cell,oneD,outDen,noco%l_noco,.TRUE.,.true.,fix)
      CALL closeXMLElement('allElectronCharges')

      IF (input%jspins.EQ.2) THEN
         noco_new = noco

         !Calculate and write out spin densities at the nucleus and magnetic moments in the spheres
         CALL magMoms(dimension,input,atoms,noco_new,vTot,moments)

         !Generate and save the new nocoinp file if the directions of the local
         !moments are relaxed or a constraint B-field is calculated.
         IF (ANY(noco%l_relax(:atoms%ntype)).OR.noco%l_constr) THEN
            CALL genNewNocoInp(input,atoms,noco,noco_new)
         END IF

         IF (noco%l_soc) CALL orbMagMoms(dimension,atoms,noco,moments%clmom)
      END IF
   END IF ! mpi%irank.EQ.0

#ifdef CPP_MPI
   CALL mpi_bc_potden(mpi,stars,sphhar,atoms,input,vacuum,oneD,noco,outDen)
#endif

END SUBROUTINE cdngen

END MODULE m_cdngen
