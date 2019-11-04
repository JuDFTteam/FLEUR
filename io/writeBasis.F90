!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_writeBasis

CONTAINS

SUBROUTINE writeBasis(input,noco,kpts,atoms,sym,cell,enpara,vTot,vCoul,vx,mpi,results,eig_id,oneD,sphhar,stars,vacuum)

   USE m_types
   USE m_juDFT
#ifdef CPP_HDF   
   USE hdf5
   USE m_hdf_tools
#endif   
   USE m_genmtbasis
   USE m_pot_io
   USE m_abcof
   USE m_abcrot
   USE m_eig66_io, ONLY : read_eig

   IMPLICIT NONE
!     TYPE(t_results),INTENT(IN)    :: results
      
      TYPE(t_enpara),INTENT(IN)     :: enpara
!     TYPE(t_banddos),INTENT(IN)    :: banddos

      TYPE(t_sphhar),INTENT(IN)     :: sphhar
      TYPE(t_stars),INTENT(IN)      :: stars
      TYPE(t_vacuum),INTENT(IN)     :: vacuum

      TYPE(t_input),INTENT(IN)      :: input
      TYPE(t_noco),INTENT(IN)       :: noco
      TYPE(t_kpts),INTENT(IN)       :: kpts
      TYPE(t_atoms),INTENT(INOUT)   :: atoms
      TYPE(t_sym),INTENT(IN)        :: sym
      TYPE(t_cell),INTENT(IN)       :: cell
      TYPE(t_potden), INTENT(INOUT) :: vTot
      TYPE(t_potden), INTENT(INOUT) :: vCoul
      TYPE(t_potden), INTENT(INOUT) :: vx
      TYPE(t_mpi), INTENT(IN)       :: mpi
      TYPE(t_results), INTENT(INOUT):: results
      INTEGER, INTENT(IN)           :: eig_id
      TYPE(t_oneD), INTENT(IN)      :: oneD

      TYPE (t_usdus)                :: usdus
      TYPE(t_lapw)                  :: lapw
      TYPE (t_eigVecCoeffs)         :: eigVecCoeffs
      TYPE (t_force)                :: force
      TYPE(t_mat)                   :: zMat

#ifdef CPP_HDF   
      
      LOGICAL           :: l_exist
      CHARACTER(LEN=30) :: filename
      CHARACTER(LEN=50) :: kpt_name
      CHARACTER(LEN=30) :: jsp_name
      CHARACTER(LEN=30) :: itype_name
!     CHARACTER(LEN=30) :: l_name

   
      INTEGER(HID_T)    :: fileID
      INTEGER(HID_T)    :: metaGroupID
      INTEGER(HID_T)    :: generalGroupID
      INTEGER(HID_T)    :: cellGroupID
      INTEGER(HID_T)    :: atomsGroupID
      INTEGER(HID_T)    :: kptsGroupID
      INTEGER(HID_T)    :: kptGroupID
      INTEGER(HID_T)    :: jspGroupID
      INTEGER(HID_T)    :: itypeGroupID,  itypeSpaceID, itypeSetID
!     INTEGER(HID_T)    :: lGroupID

!     INTEGER(HID_T)    :: stringTypeID
!     INTEGER(SIZE_T)   :: stringLength

      INTEGER(HID_T)    :: bravaisMatrixSpaceID, bravaisMatrixSetID
      INTEGER(HID_T)    :: reciprocalCellSpaceID, reciprocalCellSetID

      INTEGER(HID_T)    :: atomPosSpaceID, atomPosSetID
      INTEGER(HID_T)    :: gvecSpaceID, gvecSetID
      INTEGER(HID_T)    :: eigSpaceID, eigSetID
      INTEGER(HID_T)    :: zmatSpaceID, zmatSetID
      INTEGER(HID_T)    :: atomicNumbersSpaceID, atomicNumbersSetID
      INTEGER(HID_T)    :: equivAtomsClassSpaceID, equivAtomsClassSetID

      INTEGER(HID_T)    :: kptCoordSpaceID, kptCoordSetID
      INTEGER(HID_T)    :: kptWeightSpaceID, kptWeightSetID
!     INTEGER(HID_T)    :: kptSPLabelsSpaceID, kptSPLabelsSetID
!     INTEGER(HID_T)    :: kptsSPIndicesSpaceID, kptsSPIndicesSetID
      INTEGER(HSIZE_T)  :: dims(7)

      INTEGER           :: j, iAtom, i
!     INTEGER           :: noded, nodeu

      INTEGER           :: hdfError, dimsInt(7)
      INTEGER           :: version
      REAL, ALLOCATABLE :: output(:,:,:,:),output5(:,:,:,:,:),output3(:,:,:)
!     INTEGER           :: fakeLogical
!     REAL              :: eFermiPrev
!     LOGICAL           :: l_error

      INTEGER           :: atomicNumbers(atoms%nat),ngopr_temp(atoms%nat)
      INTEGER           :: equivAtomsGroup(atoms%nat)



      REAL, ALLOCATABLE :: f(:,:,:,:),g(:,:,:,:)
      REAL, ALLOCATABLE :: flo(:,:,:)

      real, allocatable :: cof(:,:,:,:)
      integer(HSIZE_T)  :: Hdim1(4)
      INTEGER           :: lmn, na,lm,n,nn, m
      complex,parameter :: img=(0.,1.)
      REAL              :: bk(3)

      LOGICAL l_zref,l_real, link_exists
      INTEGER jsp,nk,l,itype
      INTEGER numbands, nbasfcn, ndbands !ndbands number of bands without highest (degenerate)

    !WRITE(5000,*) 'writeBasis entry'

    CALL force%init1(input,atoms)

    IF (noco%l_mperp) THEN
       ALLOCATE ( f(atoms%jmtd,2,0:atoms%lmaxd,input%jspins),g(atoms%jmtd,2,0:atoms%lmaxd,input%jspins) )
    ELSE
       ALLOCATE ( f(atoms%jmtd,2,0:atoms%lmaxd,1:input%jspins) )
       ALLOCATE ( g(atoms%jmtd,2,0:atoms%lmaxd,1:input%jspins) )
    ENDIF
    ALLOCATE (flo(atoms%jmtd,2,atoms%nlod))
    flo(:,:,:) = 0.0


    !-------------------------write potential--------------------
    IF(input%gw==1) THEN
       CALL writePotential(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,POT_ARCHIVE_TYPE_TOT_const,vTot%iter,vTot,vTot%pw_w)
       CALL writePotential(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,POT_ARCHIVE_TYPE_COUL_const,vCoul%iter,vCoul,vCoul%pw_w)
       CALL writePotential(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,POT_ARCHIVE_TYPE_X_const,vx%iter,vx,vx%pw_w)
    END IF


      l_real=sym%invs.AND..NOT.noco%l_noco
!     check if z-reflection trick can be used
      l_zref=(sym%zrfs.AND.(SUM(ABS(kpts%bk(3,:kpts%nkpt))).LT.1e-9).AND..NOT.noco%l_noco)
!     IF (mpi%n_size > 1) l_zref = .FALSE.
      version = 1
      filename = 'basis.hdf'

      INQUIRE(FILE=TRIM(ADJUSTL(filename)),EXIST=l_exist)
      IF(l_exist) THEN
         CALL system('rm '//TRIM(ADJUSTL(filename)))       
      END IF

      CALL h5fcreate_f(TRIM(ADJUSTL(filename)), H5F_ACC_TRUNC_F, fileID, hdfError, H5P_DEFAULT_F, H5P_DEFAULT_F)

      CALL h5gcreate_f(fileID, '/meta', metaGroupID, hdfError)
      CALL io_write_attint0(metaGroupID,'version',version)
      CALL h5gclose_f(metaGroupID, hdfError)

      CALL h5gcreate_f(fileID, '/general', generalGroupID, hdfError)
      CALL io_write_attint0(generalGroupID,'jspins',input%jspins)
      CALL io_write_attlog0(generalGroupID,'invs',sym%invs)
      CALL io_write_attlog0(generalGroupID,'l_soc',noco%l_soc)
      CALL io_write_attlog0(generalGroupID,'l_real',l_real)
      CALL io_write_attreal0(generalGroupID,'rkmax',input%rkmax)
      CALL h5gclose_f(generalGroupID, hdfError)      

      CALL h5gcreate_f(fileID, '/cell', cellGroupID, hdfError)

      dims(:2)=(/3,3/)
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),bravaisMatrixSpaceID,hdfError)
      CALL h5dcreate_f(cellGroupID, "amat", H5T_NATIVE_DOUBLE, bravaisMatrixSpaceID, bravaisMatrixSetID, hdfError)
      CALL h5sclose_f(bravaisMatrixSpaceID,hdfError)
      CALL io_write_real2(bravaisMatrixSetID,(/1,1/),dimsInt(:2),cell%amat)
      CALL h5dclose_f(bravaisMatrixSetID, hdfError)

      dims(:2)=(/3,3/)
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),reciprocalCellSpaceID,hdfError)
      CALL h5dcreate_f(cellGroupID, "bmat", H5T_NATIVE_DOUBLE, reciprocalCellSpaceID, reciprocalCellSetID, hdfError)
      CALL h5sclose_f(reciprocalCellSpaceID,hdfError)
      CALL io_write_real2(reciprocalCellSetID,(/1,1/),dimsInt(:2),cell%bmat)
      CALL h5dclose_f(reciprocalCellSetID, hdfError)
      
      CALL io_write_attreal0(cellGroupID,'scaleCell',input%scaleCell)
      CALL io_write_attreal0(cellGroupID,'scaleA1',input%scaleA1)
      CALL io_write_attreal0(cellGroupID,'scaleA2',input%scaleA2)
      CALL io_write_attreal0(cellGroupID,'scaleC',input%scaleC)


      CALL h5gclose_f(cellGroupID, hdfError)
      iAtom = 0
      DO iType = 1, atoms%ntype
         DO j = 1, atoms%neq(iType)
            iAtom = iAtom + 1
            atomicNumbers(iAtom) = atoms%nz(iType)
            equivAtomsGroup(iAtom) = iType
         END DO
      END DO

      CALL h5gcreate_f(fileID, '/atoms', atomsGroupID, hdfError)
      CALL io_write_attint0(atomsGroupID,'nAtoms',atoms%nat)
      CALL io_write_attint0(atomsGroupID,'nTypes',atoms%ntype)
      CALL io_write_attint0(atomsGroupID,'nlod',atoms%nlod)
      CALL io_write_attint0(atomsGroupID,'llod',atoms%llod)
      CALL io_write_attint0(atomsGroupID,'nlotot',atoms%nlotot)
      CALL io_write_attint0(atomsGroupID,'lmaxd',atoms%lmaxd)
      CALL io_write_attint0(atomsGroupID,'jmtd',atoms%jmtd)

      dims(:1)=(/atoms%ntype/)
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),atomicNumbersSpaceID,hdfError)
      CALL h5dcreate_f(atomsGroupID, "jri", H5T_NATIVE_INTEGER, atomicNumbersSpaceID, atomicNumbersSetID, hdfError)
      CALL h5sclose_f(atomicNumbersSpaceID,hdfError)
      CALL io_write_integer1(atomicNumbersSetID,(/1/),dimsInt(:1),atoms%jri)
      CALL h5dclose_f(atomicNumbersSetID, hdfError)

      dims(:1)=(/atoms%ntype/)
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),atomicNumbersSpaceID,hdfError)
      CALL h5dcreate_f(atomsGroupID, "lmax", H5T_NATIVE_INTEGER, atomicNumbersSpaceID, atomicNumbersSetID, hdfError)
      CALL h5sclose_f(atomicNumbersSpaceID,hdfError)
      CALL io_write_integer1(atomicNumbersSetID,(/1/),dimsInt(:1),atoms%lmax)
      CALL h5dclose_f(atomicNumbersSetID, hdfError)
      
      dims(:1)=(/atoms%ntype/)
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),atomicNumbersSpaceID,hdfError)
      CALL h5dcreate_f(atomsGroupID, "neq", H5T_NATIVE_INTEGER, atomicNumbersSpaceID, atomicNumbersSetID, hdfError)
      CALL h5sclose_f(atomicNumbersSpaceID,hdfError)
      CALL io_write_integer1(atomicNumbersSetID,(/1/),dimsInt(:1),atoms%neq)
      CALL h5dclose_f(atomicNumbersSetID, hdfError)

      dims(:1)=(/atoms%ntype/)
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),atomicNumbersSpaceID,hdfError)
      CALL h5dcreate_f(atomsGroupID, "nlo", H5T_NATIVE_INTEGER, atomicNumbersSpaceID, atomicNumbersSetID, hdfError)
      CALL h5sclose_f(atomicNumbersSpaceID,hdfError)
      CALL io_write_integer1(atomicNumbersSetID,(/1/),dimsInt(:1),atoms%nlo)
      CALL h5dclose_f(atomicNumbersSetID, hdfError)

      dims(:2)=(/atoms%nlod,atoms%ntype/)
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),atomPosSpaceID,hdfError)
      CALL h5dcreate_f(atomsGroupID, "llo", H5T_NATIVE_INTEGER, atomPosSpaceID, atomPosSetID, hdfError)
      CALL h5sclose_f(atomPosSpaceID,hdfError)
      CALL io_write_integer2(atomPosSetID,(/1,1/),dimsInt(:2),atoms%llo)
      CALL h5dclose_f(atomPosSetID, hdfError)

      dims(:2)=(/atoms%jmtd,atoms%ntype/)
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),atomPosSpaceID,hdfError)
      CALL h5dcreate_f(atomsGroupID, "rmsh", H5T_NATIVE_DOUBLE, atomPosSpaceID, atomPosSetID, hdfError)
      CALL h5sclose_f(atomPosSpaceID,hdfError)
      CALL io_write_real2(atomPosSetID,(/1,1/),dimsInt(:2),atoms%rmsh)
      CALL h5dclose_f(atomPosSetID, hdfError)

      dims(:2)=(/3,atoms%nat/)
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),atomPosSpaceID,hdfError)
      CALL h5dcreate_f(atomsGroupID, "ratoms", H5T_NATIVE_DOUBLE, atomPosSpaceID, atomPosSetID, hdfError)
      CALL h5sclose_f(atomPosSpaceID,hdfError)
      CALL io_write_real2(atomPosSetID,(/1,1/),dimsInt(:2),atoms%taual)
      CALL h5dclose_f(atomPosSetID, hdfError)

      dims(:1)=(/atoms%nat/)
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),atomicNumbersSpaceID,hdfError)
      CALL h5dcreate_f(atomsGroupID, "zatoms", H5T_NATIVE_INTEGER, atomicNumbersSpaceID, atomicNumbersSetID, hdfError)
      CALL h5sclose_f(atomicNumbersSpaceID,hdfError)
      CALL io_write_integer1(atomicNumbersSetID,(/1/),dimsInt(:1),atomicNumbers)
      CALL h5dclose_f(atomicNumbersSetID, hdfError)

      dims(:1)=(/atoms%ntype/)
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),atomicNumbersSpaceID,hdfError)
      CALL h5dcreate_f(atomsGroupID, "ztype", H5T_NATIVE_INTEGER, atomicNumbersSpaceID, atomicNumbersSetID, hdfError)
      CALL h5sclose_f(atomicNumbersSpaceID,hdfError)
      CALL io_write_integer1(atomicNumbersSetID,(/1/),dimsInt(:1),atoms%nz)
      CALL h5dclose_f(atomicNumbersSetID, hdfError)

      dims(:1)=(/atoms%nat/)
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),equivAtomsClassSpaceID,hdfError)
      CALL h5dcreate_f(atomsGroupID, "equivAtomsGroup", H5T_NATIVE_INTEGER, equivAtomsClassSpaceID, equivAtomsClassSetID, hdfError)
      CALL h5sclose_f(equivAtomsClassSpaceID,hdfError)
      CALL io_write_integer1(equivAtomsClassSetID,(/1/),dimsInt(:1),equivAtomsGroup)
      CALL h5dclose_f(equivAtomsClassSetID, hdfError)

      CALL h5gclose_f(atomsGroupID, hdfError)

      CALL h5gcreate_f(fileID, '/kpts', kptsGroupID, hdfError)

      CALL io_write_attint0(kptsGroupID,'nkpt',kpts%nkpt)

      dims(:2)=(/3,kpts%nkpt/)
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),kptCoordSpaceID,hdfError)
      CALL h5dcreate_f(kptsGroupID, "coordinates", H5T_NATIVE_DOUBLE, kptCoordSpaceID, kptCoordSetID, hdfError)
      CALL h5sclose_f(kptCoordSpaceID,hdfError)
      CALL io_write_real2(kptCoordSetID,(/1,1/),dimsInt(:2),kpts%bk)
      CALL h5dclose_f(kptCoordSetID, hdfError)

      dims(:1)=(/kpts%nkpt/)
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),kptWeightSpaceID,hdfError)
      CALL h5dcreate_f(kptsGroupID, "weights", H5T_NATIVE_DOUBLE, kptWeightSpaceID, kptWeightSetID, hdfError)
      CALL h5sclose_f(kptWeightSpaceID,hdfError)
      CALL io_write_real1(kptWeightSetID,(/1/),dimsInt(:1),kpts%wtkpt)
      CALL h5dclose_f(kptWeightSetID, hdfError)

      CALL h5gclose_f(kptsGroupID, hdfError)

      CALL usdus%init(atoms,input%jspins)

      DO jsp = 1,MERGE(1,input%jspins,noco%l_noco)
         write(jsp_name , '(a,i0)') '/jsp_',jsp
         CALL h5gcreate_f(fileID, TRIM(ADJUSTL(jsp_name)), jspGroupID, hdfError)
!        DO nk = mpi%n_start,kpts%nkpt,mpi%n_stride
         DO nk = 1,kpts%nkpt
            CALL lapw%init(input,noco,kpts,atoms,sym,nk,cell,l_zref)
            bk(:) = kpts%bk(:,nk)
            IF(abs(bk(1)).LT.1e-7) bk(1) = abs(bk(1))
            IF(abs(bk(2)).LT.1e-7) bk(2) = abs(bk(2))
            IF(abs(bk(3)).LT.1e-7) bk(3) = abs(bk(3))
            !write(kpt_name , '(2a,i0)') TRIM(ADJUSTL(jsp_name)),'/kpt_',nk
            write(kpt_name , '(2a,f12.10,a,f12.10,a,f12.10)') TRIM(ADJUSTL(jsp_name)),'/kpt_',bk(1),',',bk(2),',',bk(3)
            CALL h5lexists_f(fileID, TRIM(ADJUSTL(kpt_name)), link_exists, hdfError)
            IF (link_exists) CYCLE
            CALL h5gcreate_f(fileID, TRIM(ADJUSTL(kpt_name)), kptGroupID, hdfError)
!--------------------enter output gvec etc here--------------------
            CALL io_write_attint0(kptGroupID,'nv',lapw%nv(jsp))

            dims(:2)=(/3,lapw%nv(jsp)/)
            dimsInt=dims
            CALL h5screate_simple_f(2,dims(:2),gvecSpaceID,hdfError)
            CALL h5dcreate_f(kptGroupID, "gvec", H5T_NATIVE_INTEGER, gvecSpaceID, gvecSetID, hdfError)
            CALL h5sclose_f(gvecSpaceID,hdfError)
            CALL io_write_integer2(gvecSetID,(/1,1/),dimsInt(:2),lapw%gvec(:,:lapw%nv(jsp),jsp))
            CALL h5dclose_f(gvecSetID, hdfError)	
	      
            CALL h5gclose_f(kptGroupID, hdfError)   
         END DO

         DO itype = 1,atoms%ntype
	    write(itype_name , '(2a,i0)') TRIM(ADJUSTL(jsp_name)),'/itype_',itype
	    CALL h5gcreate_f(fileID, TRIM(ADJUSTL(itype_name)), itypeGroupID, hdfError)

            CALL genMTBasis(atoms,enpara,vTot,mpi,itype,jsp,usdus,f(:,:,0:,jsp),g(:,:,0:,jsp),flo)
	    dims(:3)=(/atoms%jmtd,2,atoms%lmaxd+1/)
	    dimsInt = dims
	    CALL h5screate_simple_f(3,dims(:3),itypeSpaceID,hdfError)
	    CALL h5dcreate_f(itypeGroupID, "f", H5T_NATIVE_DOUBLE, itypeSpaceID, itypeSetID, hdfError)
	    CALL h5sclose_f(itypeSpaceID,hdfError)
	    CALL io_write_real3(itypeSetID,(/1,1,1/),dimsInt(:3),f(:,:,0:,jsp))
	    CALL h5dclose_f(itypeSetID, hdfError)

	    CALL h5screate_simple_f(3,dims(:3),itypeSpaceID,hdfError)
	    CALL h5dcreate_f(itypeGroupID, "g", H5T_NATIVE_DOUBLE, itypeSpaceID, itypeSetID, hdfError)
	    CALL h5sclose_f(itypeSpaceID,hdfError)
	    CALL io_write_real3(itypeSetID,(/1,1,1/),dimsInt(:3),g(:,:,0:,jsp))
	    CALL h5dclose_f(itypeSetID, hdfError)

	    dims(:3)=(/atoms%jmtd,2,atoms%nlod/)
	    dimsInt = dims
	    CALL h5screate_simple_f(3,dims(:3),itypeSpaceID,hdfError)
	    CALL h5dcreate_f(itypeGroupID, "flo", H5T_NATIVE_DOUBLE, itypeSpaceID, itypeSetID, hdfError)
	    CALL h5sclose_f(itypeSpaceID,hdfError)
	    CALL io_write_real3(itypeSetID,(/1,1,1/),dimsInt(:3),flo(:,:,:))
	    CALL h5dclose_f(itypeSetID, hdfError)

	    dims(:1)=(/atoms%lmaxd+1/)
	    dimsInt = dims
	    CALL h5screate_simple_f(1,dims(:1),itypeSpaceID,hdfError)
	    CALL h5dcreate_f(itypeGroupID, "us", H5T_NATIVE_DOUBLE, itypeSpaceID, itypeSetID, hdfError)
	    CALL h5sclose_f(itypeSpaceID,hdfError)
	    CALL io_write_real1(itypeSetID,(/1/),dimsInt(:1),usdus%us(:,itype,jsp))
	    CALL h5dclose_f(itypeSetID, hdfError)

	    CALL h5screate_simple_f(1,dims(:1),itypeSpaceID,hdfError)
	    CALL h5dcreate_f(itypeGroupID, "dus", H5T_NATIVE_DOUBLE, itypeSpaceID, itypeSetID, hdfError)
	    CALL h5sclose_f(itypeSpaceID,hdfError)
	    CALL io_write_real1(itypeSetID,(/1/),dimsInt(:1),usdus%dus(:,itype,jsp))
	    CALL h5dclose_f(itypeSetID, hdfError)

	    CALL h5screate_simple_f(1,dims(:1),itypeSpaceID,hdfError)
	    CALL h5dcreate_f(itypeGroupID, "uds", H5T_NATIVE_DOUBLE, itypeSpaceID, itypeSetID, hdfError)
	    CALL h5sclose_f(itypeSpaceID,hdfError)
	    CALL io_write_real1(itypeSetID,(/1/),dimsInt(:1),usdus%uds(:,itype,jsp))
	    CALL h5dclose_f(itypeSetID, hdfError)

	    CALL h5screate_simple_f(1,dims(:1),itypeSpaceID,hdfError)
	    CALL h5dcreate_f(itypeGroupID, "duds", H5T_NATIVE_DOUBLE, itypeSpaceID, itypeSetID, hdfError)
	    CALL h5sclose_f(itypeSpaceID,hdfError)
	    CALL io_write_real1(itypeSetID,(/1/),dimsInt(:1),usdus%duds(:,itype,jsp))
	    CALL h5dclose_f(itypeSetID, hdfError)

	    dims(:1)=(/atoms%nlod/)
	    dimsInt = dims
	    CALL h5screate_simple_f(1,dims(:1),itypeSpaceID,hdfError)
	    CALL h5dcreate_f(itypeGroupID, "dulos", H5T_NATIVE_DOUBLE, itypeSpaceID, itypeSetID, hdfError)
	    CALL h5sclose_f(itypeSpaceID,hdfError)
	    CALL io_write_real1(itypeSetID,(/1/),dimsInt(:1),usdus%dulos(:,itype,jsp))
	    CALL h5dclose_f(itypeSetID, hdfError)

	    CALL h5screate_simple_f(1,dims(:1),itypeSpaceID,hdfError)
	    CALL h5dcreate_f(itypeGroupID, "ulos", H5T_NATIVE_DOUBLE, itypeSpaceID, itypeSetID, hdfError)
	    CALL h5sclose_f(itypeSpaceID,hdfError)
	    CALL io_write_real1(itypeSetID,(/1/),dimsInt(:1),usdus%ulos(:,itype,jsp))
	    CALL h5dclose_f(itypeSetID, hdfError)

	    dims(:1)=(/atoms%lmaxd+1/)
	    dimsInt = dims
	    CALL h5screate_simple_f(1,dims(:1),itypeSpaceID,hdfError)
	    CALL h5dcreate_f(itypeGroupID, "el0", H5T_NATIVE_DOUBLE, itypeSpaceID, itypeSetID, hdfError)
	    CALL h5sclose_f(itypeSpaceID,hdfError)
	    CALL io_write_real1(itypeSetID,(/1/),dimsInt(:1),enpara%el0(:,itype,jsp))
	    CALL h5dclose_f(itypeSetID, hdfError)

	    dims(:1)=(/atoms%nlod/)
	    dimsInt = dims
	    CALL h5screate_simple_f(1,dims(:1),itypeSpaceID,hdfError)
	    CALL h5dcreate_f(itypeGroupID, "ello0", H5T_NATIVE_DOUBLE, itypeSpaceID, itypeSetID, hdfError)
	    CALL h5sclose_f(itypeSpaceID,hdfError)
	    CALL io_write_real1(itypeSetID,(/1/),dimsInt(:1),enpara%ello0(:,itype,jsp))
	    CALL h5dclose_f(itypeSetID, hdfError)


!           write(kpt_name , '(3a,i0)') TRIM(ADJUSTL(jsp_name)),TRIM(ADJUSTL(itype_name)),'/l_',l
!	    CALL h5gcreate_f(fileID, TRIM(ADJUSTL(l_name)), lGroupID, hdfError)
!	    CALL h5gclose_f(lGroupID, hdfError)
	    CALL h5gclose_f(itypeGroupID, hdfError)
	END DO

       CALL h5gclose_f(jspGroupID, hdfError)
    END DO
    CALL h5fclose_f(fileID, hdfError)
    version = 1
    filename = 'eig_gw.hdf'


    IF(input%gw==2) THEN
      INQUIRE(FILE=TRIM(ADJUSTL(filename)),EXIST=l_exist)
      IF(l_exist) THEN
        CALL system('rm '//TRIM(ADJUSTL(filename)))       
      END IF

      CALL h5fcreate_f(TRIM(ADJUSTL(filename)), H5F_ACC_TRUNC_F, fileID, hdfError, H5P_DEFAULT_F, H5P_DEFAULT_F)

      CALL h5gcreate_f(fileID, '/meta', metaGroupID, hdfError)
      CALL io_write_attint0(metaGroupID,'version',version)
      CALL h5gclose_f(metaGroupID, hdfError)

   DO jsp = 1,MERGE(1,input%jspins,noco%l_noco)
       write(jsp_name , '(a,i0)') '/jsp_',jsp
       CALL h5gcreate_f(fileID, TRIM(ADJUSTL(jsp_name)), jspGroupID, hdfError)
!      DO nk = mpi%n_start,kpts%nkpt,mpi%n_stride
       DO nk = 1,kpts%nkpt
            CALL lapw%init(input,noco,kpts,atoms,sym,nk,cell,l_zref)
            bk(:) = kpts%bk(:,nk)
            IF(abs(bk(1)).LT.1e-7) bk(1) = abs(bk(1))
            IF(abs(bk(2)).LT.1e-7) bk(2) = abs(bk(2))
            IF(abs(bk(3)).LT.1e-7) bk(3) = abs(bk(3))
            !write(kpt_name , '(2a,i0)') TRIM(ADJUSTL(jsp_name)),'/kpt_',nk
            write(kpt_name , '(2a,f12.10,a,f12.10,a,f12.10)') TRIM(ADJUSTL(jsp_name)),'/kpt_',bk(1),',',bk(2),',',bk(3)
            CALL h5lexists_f(fileID, TRIM(ADJUSTL(kpt_name)), link_exists, hdfError)
            IF (link_exists) CYCLE
            CALL h5gcreate_f(fileID, TRIM(ADJUSTL(kpt_name)), kptGroupID, hdfError)
!--------------------abcoff, zmat, eig output here-------------------
!,results%neig(nk,jsp),results%eig(:,nk,jsp)
            numbands=results%neig(nk,jsp)
            nbasfcn = MERGE(lapw%nv(1)+lapw%nv(2)+2*atoms%nlotot,lapw%nv(1)+atoms%nlotot,noco%l_noco)
            CALL zMat%init(l_real,nbasfcn,numbands)
            CALL read_eig(eig_id,nk,jsp,zmat=zMat)
	    CALL eigVecCoeffs%init(input,atoms,noco,jsp,numbands)
            IF (input%l_f) CALL force%init2(numbands,input,atoms)
!            DO i=1,atoms%nat
!	    	ngopr_temp(i)=sym%ngopr(i)
!               sym%ngopr(i)=1
!            END DO
		CALL abcof(input,atoms,sym,cell,lapw,numbands,usdus,noco,jsp,oneD,&
		    eigVecCoeffs%acof(:,0:,:,jsp),eigVecCoeffs%bcof(:,0:,:,jsp),&
		    eigVecCoeffs%ccof(-atoms%llod:,:,:,:,jsp),zMat,results%eig(:,nk,jsp),force) 
!            DO i=1,atoms%nat
!	     	sym%ngopr(i)=ngopr_temp(i)
!            END DO
		CALL abcrot(atoms%ntype,atoms%nat,numbands,atoms%lmaxd,atoms%lmaxd*(atoms%lmaxd+2),atoms%llod,atoms%nlod,atoms%ntype,atoms%neq,&
		            numbands,atoms%lmax,atoms%nlo,atoms%llo,sym%nop,sym%ngopr,sym%mrot,sym%invsat,sym%invsatnr,cell%bmat,&
		           oneD%odi,oneD%ods,&
		           eigVecCoeffs%acof(:,0:,:,jsp),eigVecCoeffs%bcof(:,0:,:,jsp),eigVecCoeffs%ccof(-atoms%llod:,:,:,:,jsp))
!-------------------------for spex output: nbasfcn=nv(because lo info not needed) and numbands setting to numbands without highest (degenerat) state-------- 
!                nbasfcn= MERGE(lapw%nv(1)+lapw%nv(2)+2*atoms%nlotot,lapw%nv(1)+atoms%nlotot,noco%l_noco)
		ndbands=numbands-1
		DO i=(numbands-1),1,-1
		    IF (abs(results%eig(i+1,nk,jsp)-results%eig(i,nk,jsp)).LT.0.000001) THEN
		        ndbands=ndbands-1
		    ELSE
			EXIT
		    END IF
		END DO
                write(*,*)numbands,ndbands
		numbands=ndbands
!------------------------setting variables numbands and nbasfcn end -------------------
		!CALL read_eig(eig_id,nk,jsp,eig=results%eig(:,nk,jsp),zmat=zMat)
		dims(:1)=(/numbands/)
		dimsInt=dims
		CALL h5screate_simple_f(1,dims(:1),eigSpaceID,hdfError)
		CALL h5dcreate_f(kptGroupID, "eig", H5T_NATIVE_DOUBLE, eigSpaceID, eigSetID, hdfError)
		CALL h5sclose_f(eigSpaceID,hdfError)
                CALL io_write_real1(eigSetID,(/1/),dimsInt(:1),results%eig(:numbands,nk,jsp))
		CALL h5dclose_f(eigSetID, hdfError)
   
                CALL io_write_attint0(kptGroupID,'numbands',numbands)            
		IF (zMat%l_real) THEN
		      dims(:2)=(/nbasfcn,numbands/)
		      dimsInt=dims
		      CALL h5screate_simple_f(2,dims(:2),kptCoordSpaceID,hdfError)
		      CALL h5dcreate_f(kptGroupID, "pw", H5T_NATIVE_DOUBLE, kptCoordSpaceID, kptCoordSetID, hdfError)
		      CALL h5sclose_f(kptCoordSpaceID,hdfError)
            CALL io_write_real2(kptCoordSetID,(/1,1/),dimsInt(:2),zMat%data_r(:nbasfcn,:numbands))
		      CALL h5dclose_f(kptCoordSetID, hdfError)
		ELSE
                      AllOCATE(output3(2,nbasfcn,numbands))
!			SIZE(zMat%data_c,1),SIZE(zMat%data_c,2)
		      output3(1,:,:)=REAL(zMat%data_c(:,:numbands))
		      output3(2,:,:)=AIMAG(zMat%data_c(:,:numbands))
		      dims(:3)=(/2,nbasfcn,numbands/)
		      dimsInt=dims
		      CALL h5screate_simple_f(3,dims(:3),zmatSpaceID,hdfError)
		      CALL h5dcreate_f(kptGroupID, "pw", H5T_NATIVE_DOUBLE, zmatSpaceID, zmatSetID, hdfError)
		      CALL h5sclose_f(zmatSpaceID,hdfError)
		      CALL io_write_real3(zmatSetID,(/1,1,1/),dimsInt(:3),output3(:2,:nbasfcn,:numbands))
		      CALL h5dclose_f(zmatSetID, hdfError)
		      DEAllOCATE(output3)
		END IF
		!AllOCATE(output(2,numbands,atoms%lmaxd*(atoms%lmaxd+2)+1,atoms%nat))
		!output(1,:,:,:)=REAL(eigVecCoeffs%acof(:,0:,:,jsp))
		!output(2,:,:,:)=AIMAG(eigVecCoeffs%acof(:,0:,:,jsp))
		!dims(:4)=(/2,numbands,atoms%lmaxd*(atoms%lmaxd+2)+1,atoms%nat/)
		!dimsInt = dims
                !CALL h5screate_simple_f(4,dims(:4),itypeSpaceID,hdfError)
		!CALL h5dcreate_f(kptGroupID, "acof", H5T_NATIVE_DOUBLE, itypeSpaceID, itypeSetID, hdfError)
		!CALL h5sclose_f(itypeSpaceID,hdfError)
		!CALL io_write_real4(itypeSetID,(/1,1,1,1/),dimsInt(:4), output)
		!CALL h5dclose_f(itypeSetID, hdfError)
		!DEAllOCATE(output)

		!AllOCATE(output(2,numbands,atoms%lmaxd*(atoms%lmaxd+2)+1,atoms%nat))
		!output(1,:,:,:)=REAL(eigVecCoeffs%bcof(:,0:,:,jsp))
		!output(2,:,:,:)=AIMAG(eigVecCoeffs%bcof(:,0:,:,jsp))
		!dims(:4)=(/2,numbands,atoms%lmaxd*(atoms%lmaxd+2)+1,atoms%nat/)
		!dimsInt = dims
                !CALL h5screate_simple_f(4,dims(:4),itypeSpaceID,hdfError)
		!CALL h5dcreate_f(kptGroupID, "bcof", H5T_NATIVE_DOUBLE, itypeSpaceID, itypeSetID, hdfError)
		!CALL h5sclose_f(itypeSpaceID,hdfError)
		!CALL io_write_real4(itypeSetID,(/1,1,1,1/),dimsInt(:4), output)
		!CALL h5dclose_f(itypeSetID, hdfError)
		!DEAllOCATE(output)

		!AllOCATE(output5(2,atoms%llod*2+1,numbands,atoms%nlod,atoms%nat))
		!output5(1,:,:,:,:)=REAL(eigVecCoeffs%ccof(-atoms%llod:,:,:,:,jsp))
		!output5(2,:,:,:,:)=AIMAG(eigVecCoeffs%ccof(-atoms%llod:,:,:,:,jsp))
		!dims(:5)=(/2,atoms%llod*2+1,numbands,atoms%nlod,atoms%nat/)
		!dimsInt = dims
                !CALL h5screate_simple_f(5,dims(:5),itypeSpaceID,hdfError)
		!CALL h5dcreate_f(kptGroupID, "ccof", H5T_NATIVE_DOUBLE, itypeSpaceID, itypeSetID, hdfError)
		!CALL h5sclose_f(itypeSpaceID,hdfError)
		!CALL io_write_real5(itypeSetID,(/1,1,1,1,1/),dimsInt(:5), output5)
		!CALL h5dclose_f(itypeSetID, hdfError)
		!DEAllOCATE(output5)

!-------------------------output spex format-------------------
	      ! get dimensions for MT cof
	      lmn = 0
	      na  = 0
	      do n = 1,atoms%ntype
		do nn = 1,atoms%neq(n)
		  na = na + 1
		  lm = 0
		  do l = 0,atoms%lmax(n)
		    lm = lm + 2*(2*l+1)
		  enddo
		  do i = 1,atoms%nlo(n)
		    lm = lm + 2*atoms%llo(i,n)+1
		  enddo
		  lmn = max(lmn,lm)
		enddo
	      enddo
	      Hdim1 = (/ 2,lmn,na,numbands /)

		allocate ( cof(Hdim1(1),Hdim1(2),Hdim1(3),Hdim1(4)) )
		      cof = 0
		      na  = 0
		      do n = 1,atoms%ntype
			do nn = 1,atoms%neq(n)
			  na = na + 1
			  lm  = 0
			  lmn = 0
			  do l = 0,atoms%lmax(n)
			    do m = -l,l              
			      cof(1,lmn+1,na,:) = real ( eigVecCoeffs%acof(:numbands,lm,na,jsp) * img**l )
			      cof(1,lmn+2,na,:) = real ( eigVecCoeffs%bcof(:numbands,lm,na,jsp) * img**l )
			      cof(2,lmn+1,na,:) = aimag ( eigVecCoeffs%acof(:numbands,lm,na,jsp) * img**l )
			      cof(2,lmn+2,na,:) = aimag ( eigVecCoeffs%bcof(:numbands,lm,na,jsp) * img**l )
			      lm  = lm  + 1
			      lmn = lmn + 2
			      do i = 1,atoms%nlo(n)
				if(atoms%llo(i,n).eq.l) then
				  lmn = lmn + 1
				  cof(1,lmn,na,:) = real( eigVecCoeffs%ccof(m,:numbands,i,na,jsp) * img**l )
				  cof(2,lmn,na,:) = aimag( eigVecCoeffs%ccof(m,:numbands,i,na,jsp) * img**l )
				endif
			      enddo
			    enddo
			  enddo
			enddo
			if(lmn.gt.size(cof,2)) then
			  write(*,*) lmn,size(cof,2)
			endif
		      enddo

		dims(:4)=Hdim1
		dimsInt = dims
                CALL h5screate_simple_f(4,dims(:4),itypeSpaceID,hdfError)
		CALL h5dcreate_f(kptGroupID, "mt", H5T_NATIVE_DOUBLE, itypeSpaceID, itypeSetID, hdfError)
		CALL h5sclose_f(itypeSpaceID,hdfError)
     		CALL io_write_real4(itypeSetID,(/1,1,1,1/),dimsInt(:4), cof)
		CALL h5dclose_f(itypeSetID, hdfError)
		deallocate ( cof )
!-------------------------end output spex format-----------------

		CALL h5gclose_f(kptGroupID, hdfError)
            
        END DO
       CALL h5gclose_f(jspGroupID, hdfError)
   END DO
   CALL h5fclose_f(fileID, hdfError)  

   END IF
   
#else
   CALL juDFT_error("writeBasis called without HDF5! ",calledby="writeBasis")
#endif

END SUBROUTINE writeBasis

END MODULE m_writeBasis
