MODULE m_eig66_hdf
#include "juDFT_env.h"
    !*****************************************************************
    ! DESC:Module for hdf-io of eig-file
    !      To be compatible with f90 interface of HDF, use kind for vars
    !
    !      !ATTENTION before calling openeig and after calling closeeig!
    !      !the hdf library has to be initialized or finalized, respectively
    !
    !      CONTAINS the following subroutines:
    !      openeig        opens file
    !      closeeig       closes file
    !      read_keb       reads kpt, enpara and basis data
    !      read_neig      read no of eigenvalues (and eigenvalues itself)
    !      read_eig       reads eigenvectors
    !      writeeig       saves all data for kpt
    !      writesingleeig saves data for one kpt and energy
    !
    !
    !                          Daniel Wortmann, Tue Nov  512:07:522002
    !*****************************************************************
    USE m_eig66_data
#ifdef CPP_HDF
    USE hdf5
    USE m_hdf_tools
    IMPLICIT NONE

    PRIVATE
    INTEGER, PARAMETER :: one=1,two=2,three=3,zero=0
                             !to have the correct
                             !type for array constructors

#endif
    PUBLIC open_eig,close_eig
    PUBLIC read_eig
    PUBLIC write_eig!,writesingleeig,writeeigc,writebas

CONTAINS
    subroutine priv_find_data(id,d)
        INTEGER,INTENT(IN)::id
        TYPE(t_data_hdf),pointer:: d

        class(t_data),pointer   ::dp
        call eig66_find_data(dp,id)
        select type(dp)
            type is (t_data_hdf)
            d=>dp
            class default
            call judft_error("BUG: wrong datatype in eig66_hdf")
        END SELECT
    END subroutine
    !----------------------------------------------------------------------
    SUBROUTINE open_eig(id,mpi_comm,nmat,neig,nkpts,jspins,lmax,nlo,ntype,&
        &                   create,readonly,filename)

        !*****************************************************************
        !     opens hdf-file for eigenvectors+values
        !*****************************************************************
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: id,mpi_comm
        INTEGER, INTENT(IN) :: nmat,neig,nkpts,jspins,nlo,ntype,lmax
        LOGICAL, INTENT(IN) :: create,readonly
        CHARACTER(LEN=*),OPTIONAL :: filename

#ifdef CPP_HDF
      
        INTEGER         :: hdferr,access_mode
        INTEGER(HID_T)  :: creation_prp,access_prp,spaceid
        LOGICAL         :: l_exist
        INTEGER(HSIZE_T):: dims(5)
        TYPE(t_data_HDF),pointer::d
        !Set creation and access properties
#ifdef CPP_MPI
        INCLUDE 'mpif.h'
        IF (readonly) THEN
            access_prp=H5P_DEFAULT_f
            creation_prp=H5P_DEFAULT_f
        ELSE
            CALL h5pcreate_f(H5P_FILE_ACCESS_F, access_prp, hdferr)
            !      CALL h5pset_fapl_mpiposix_f(access_prp,MPI_COMM,
            !     +.false.,hdferr)
            CALL h5pset_fapl_mpio_f(access_prp, MPI_COMM, MPI_INFO_NULL,hdferr)
            creation_prp=H5P_DEFAULT_f !no special creation property
        ENDIF
#else
        access_prp=H5P_DEFAULT_f
        creation_prp=H5P_DEFAULT_f
#endif 
        IF (present(filename)) d%fname=filename
        call priv_find_data(id,d)
        !set access_flags according
        IF (readonly) THEN
            access_mode=H5F_ACC_RDONLY_F
        ELSE
            access_mode=H5F_ACC_RDWR_F
        ENDIF
        !     OPEN FILE and get D%FID's
        IF (create) THEN
            INQUIRE(FILE=trim(d%fname)//'.hdf',EXIST=l_exist)
            access_mode=H5F_ACC_TRUNC_F
            !         IF (l_exist) WRITE (*,*)'Warning: eig.hdf was overwritten'
            CALL h5fcreate_f(trim(d%fname)//'.hdf',access_Mode, d%fid, hdferr&
                &        ,creation_prp,access_prp)
            ! create dataspaces and datasets
            !   scalars
            dims(:2)=(/nkpts,jspins/)
            CALL h5screate_simple_f(2,dims(:2),spaceid,hdferr)
            CALL h5dcreate_f(d%fid, "neig", H5T_NATIVE_INTEGER, spaceid,&
                &                 d%neigsetid, hdferr)
            CALL h5dcreate_f(d%fid, "wk", H5T_NATIVE_DOUBLE,spaceid,&
                &                 d%wksetid, hdferr)
            CALL h5dcreate_f(d%fid, "nv", H5T_NATIVE_INTEGER, spaceid,&
                &                 d%nvsetid, hdferr)
            CALL h5dcreate_f(d%fid, "nmat", H5T_NATIVE_INTEGER, spaceid,&
                &                 d%nmatsetid, hdferr)
            CALL h5sclose_f(spaceid,hdferr)
            !   vectors
            dims(1:3)=(/two,nkpts,jspins/)
            CALL h5screate_simple_f(3,dims(:3),spaceid,hdferr)
            CALL h5dcreate_f(d%fid, "evac", H5T_NATIVE_DOUBLE, spaceid,&
                &                 d%evacsetid, hdferr)
            CALL h5sclose_f(spaceid,hdferr)
            dims(:3)=(/three,nkpts,jspins/)
            CALL h5screate_simple_f(3,dims(:3),spaceid,hdferr)
            CALL h5dcreate_f(d%fid, "bk", H5T_NATIVE_DOUBLE, spaceid,&
                &                 d%bksetid, hdferr)
            CALL h5sclose_f(spaceid,hdferr)
            dims(:3)=(/neig,nkpts,jspins/)
            CALL h5screate_simple_f(3,dims(:3),spaceid,hdferr)
            !     ew
            CALL h5dcreate_f(d%fid, "energy", H5T_NATIVE_DOUBLE, spaceid,&
                &                 d%energysetid, hdferr)
            CALL h5sclose_f(spaceid,hdferr)
            !     enparas
            dims(1:4)=(/lmax+1,ntype,nkpts,jspins/)
            CALL h5screate_simple_f(4,dims(1:4),spaceid,hdferr)
            CALL h5dcreate_f(d%fid, "el", H5T_NATIVE_DOUBLE, spaceid,&
                &                 d%esetid, hdferr)
            CALL h5sclose_f(spaceid,hdferr)
         
            dims(:4)=(/nlo,ntype,nkpts,jspins/)
            CALL h5screate_simple_f(4,dims(:4),spaceid,hdferr)
            CALL h5dcreate_f(d%fid, "ello", H5T_NATIVE_DOUBLE, spaceid,&
                &                 d%ellosetid, hdferr)
            CALL h5sclose_f(spaceid,hdferr)
            !     ev
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
            dims(:5)=(/one,nmat,neig,nkpts,jspins/)
#else
            dims(:5)=(/two,nmat,neig,nkpts,jspins/)
#endif
            CALL h5screate_simple_f(5,dims(:5),spaceid,hdferr)
            CALL h5dcreate_f(d%fid, "ev", H5T_NATIVE_DOUBLE, spaceid,&
                &                 d%evsetid, hdferr)
            CALL h5sclose_f(spaceid,hdferr)
            !      basis
            dims(:4)=(/nmat,three,nkpts,jspins/)
            CALL h5screate_simple_f(4,dims(:4),spaceid,hdferr)
            CALL h5dcreate_f(d%fid, "k", H5T_NATIVE_INTEGER, spaceid,&
                &                 d%ksetid, hdferr)
            CALL h5sclose_f(spaceid,hdferr)
        ELSE
            CALL h5fopen_f (trim(d%fname)//'.hdf', access_Mode, d%fid,&
                &            hdferr,access_prp)
            !get dataset-ids
            CALL h5dopen_f(d%fid, 'el', d%esetid, hdferr)
            CALL h5dopen_f(d%fid, 'evac', d%evacsetid, hdferr)
            CALL h5dopen_f(d%fid, 'ello', d%ellosetid, hdferr)
            CALL h5dopen_f(d%fid, 'bk', d%bksetid, hdferr)
            CALL h5dopen_f(d%fid, 'wk', d%wksetid, hdferr)
            CALL h5dopen_f(d%fid, 'energy', d%energysetid, hdferr)
            CALL h5dopen_f(d%fid, 'k', d%ksetid, hdferr)
            CALL h5dopen_f(d%fid, 'neig', d%neigsetid, hdferr)
            CALL h5dopen_f(d%fid, 'ev', d%evsetid, hdferr)
            CALL h5dopen_f(d%fid, 'nv', d%nvsetid, hdferr)
            CALL h5dopen_f(d%fid, 'nmat', d%nmatsetid, hdferr)
        ENDIF
        IF (.NOT.access_prp==H5P_DEFAULT_f) CALL H5Pclose_f(access_prp&
            &     ,hdferr)
#else
        CALL juDFT_error("Could not use HDF5 for IO, please recompile")
#endif
    END SUBROUTINE open_eig
    !----------------------------------------------------------------------
    SUBROUTINE close_eig(id,filename)
        !*****************************************************************
        !     closes hdf-file for eigenvectors+values
        !*****************************************************************
        IMPLICIT NONE
        INTEGER,INTENT(IN)                   :: id
        CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: filename

        INTEGER::hdferr
        TYPE(t_data_HDF),pointer::d

        !close datasets
#ifdef CPP_HDF
        call priv_find_data(id,d)

        CALL h5dclose_f(d%esetid,hdferr)
        CALL h5dclose_f(d%evacsetid,hdferr)
        CALL h5dclose_f(d%ellosetid,hdferr)
        CALL h5dclose_f(d%bksetid,hdferr)
        CALL h5dclose_f(d%wksetid,hdferr)
        CALL h5dclose_f(d%energysetid,hdferr)
        CALL h5dclose_f(d%ksetid,hdferr)
        CALL h5dclose_f(d%neigsetid,hdferr)
        CALL h5dclose_f(d%evsetid,hdferr)
        CALL h5dclose_f(d%nvsetid,hdferr)
        CALL h5dclose_f(d%nmatsetid,hdferr)
        !close file
        CALL h5fclose_f(d%fid,hdferr)
        !If a filename was given and the name is not the current filename
        IF (present(filename)) THEN
            IF (filename.NE.d%fname) THEN
                CALL system("mv "//trim(d%fname)//".hdf "//trim(filename)//".hdf")
            ENDIF
        ENDIF
        d%fname="eig"
        call eig66_remove_data(id)

#endif
    END SUBROUTINE
#ifdef CPP_HDF
    !----------------------------------------------------------------------
    SUBROUTINE priv_r_vec(d,nk,jspin,n_start,n_end,nmat,z)

        USE m_hdf_tools
        IMPLICIT NONE
        TYPE(t_data_HDF),INTENT(IN)::d
        INTEGER, INTENT(IN)  :: nk,jspin
        INTEGER, INTENT(IN)  :: n_start,n_end
        INTEGER, INTENT(OUT) :: nmat
        REAL,    INTENT(OUT) :: z(:,:)

        INTEGER i,j,neig_l

        neig_l = n_end - n_start + 1

        ! read matrix size
        CALL io_read_integer0(d%nmatsetid,(/nk,jspin/),(/1,1/),nmat)

        IF ( nmat > size(z,1) .OR. neig_l > size(z,2) ) THEN
            WRITE (6,*) nmat,size(z,1),size(z,2)
            CALL juDFT_error("eig66_hdf$read_vec",calledby ="eig66_hdf")
        ENDIF

        !read eigenvectors
        CALL io_read_real2(d%evsetid,(/1,1,n_start,nk,jspin/),&
            &                           (/1,nmat,neig_l,1,1/),&
            &                           z(:nmat,:neig_l) )

    END SUBROUTINE priv_r_vec

#endif

    SUBROUTINE write_eig(id,nk,jspin,neig,neig_total,nv,nmat,k1,k2,k3,bk,wk,&
        &                  eig,el,ello,evac,&
        &                  nlotot,kveclo,n_size,n_rank,z)

        !*****************************************************************
        !     writes all eignevecs for the nk-th kpoint
        !*****************************************************************
        IMPLICIT NONE

        INTEGER, INTENT(IN)          :: id,nk,jspin
        INTEGER, INTENT(IN),OPTIONAL :: n_size,n_rank
        REAL,    INTENT(IN),OPTIONAL :: wk
        INTEGER, INTENT(IN),OPTIONAL :: neig,nv,nmat,nlotot,neig_total
        INTEGER, INTENT(IN),OPTIONAL :: k1(:),k2(:),k3(:),kveclo(:)
        REAL,    INTENT(IN),OPTIONAL :: bk(3),eig(:),el(:,:)
        REAL,    INTENT(IN),OPTIONAL :: evac(2),ello(:,:)
        CLASS(*),INTENT(IN),OPTIONAL :: z(:,:)

        INTEGER i,j,k,nv_local,n1,n2,ne
        TYPE(t_data_HDF),pointer::d
        call priv_find_data(id,d)

#ifdef CPP_HDF
        !
        !write enparas
        !
        nv_local=huge(1)

        IF (present(el))&
            &   CALL io_write_real2(&
            &                    d%esetid,(/1,1,nk,jspin/),&
            &                    (/size(el,1),size(el,2),1,1/),el)

        IF (present(ello))&
            & CALL io_write_real2(&
            &                    d%ellosetid,(/1,1,nk,jspin/),&
            &                    (/size(ello,1),size(ello,2),1,1/),ello)

        IF (present(evac)) CALL io_write_real1(&
            &                    d%evacsetid,(/1,nk,jspin/),(/2,1,1/),evac)
        !
        !write kpts
        !

        IF (present(bk)) CALL io_write_real1(&
            &                    d%bksetid,(/1,nk,jspin/),(/3,1,1/),bk)

        IF (present(wk)) CALL io_write_real0(&
            &                    d%wksetid,(/nk,jspin/),(/1,1/),wk)
        !
        !write basis
        !

        IF (present(nv)) THEN
            nv_local=nv
            CALL io_write_integer0(d%nvsetid,(/nk,jspin/),(/1,1/),nv)
        ENDIF

        IF (present(nmat)) CALL io_write_integer0(&
            &                       d%nmatsetid,(/nk,jspin/),(/1,1/),nmat)

        IF (present(k1)) CALL io_write_integer1(&
            &              d%ksetid,(/1,1,nk,jspin/),&
            &     (/min(nv_local,size(k1)),1,1,1/),k1(:min(nv_local,size(k1))))

        IF (present(k2)) CALL io_write_integer1(&
            &              d%ksetid,(/1,2,nk,jspin/),&
            &     (/min(nv_local,size(k2)),1,1,1/),k2(:min(nv_local,size(k2))))

        IF (present(k3)) CALL io_write_integer1(&
            &              d%ksetid,(/1,3,nk,jspin/),&
            &     (/min(nv_local,size(k3)),1,1,1/),k3(:min(nv_local,size(k3))))

        IF (present(kveclo).AND.present(nlotot).AND.&
            &      present(k1).AND.present(k2).AND.present(k3)) THEN

            WRITE(*,*) kveclo,nlotot
            DO k = 1, nlotot
                CALL io_write_integer0(&
                    &        d%ksetid,(/nv+k,1,nk,jspin/),(/1,1,1,1/),k1(kveclo(k)))
                CALL io_write_integer0(&
                    &        d%ksetid,(/nv+k,2,nk,jspin/),(/1,1,1,1/),k2(kveclo(k)))
                CALL io_write_integer0(&
                    &        d%ksetid,(/nv+k,3,nk,jspin/),(/1,1,1,1/),k3(kveclo(k)))
            ENDDO
        ENDIF
        !
        !write eigenvalues
        !

        IF (present(neig_total)) THEN
            CALL io_write_integer0(d%neigsetid,(/nk,jspin/),(/1,1/),neig_total)
        ENDIF

        IF (present(n_rank).AND.present(n_size).AND.&
            &        present(eig).AND.present(neig)) THEN
            CALL io_write_real1s(&
                &                     d%energysetid,(/n_rank+1,nk,jspin/),        &
                &                     (/neig,1,1/),eig(:neig),(/n_size,1,1/))
        !write eigenvectors
        !
        ELSEIF (present(eig).AND.present(neig)) THEN
            CALL io_write_real1s(&
                &                     d%energysetid,(/1,nk,jspin/),&
                &                     (/neig,1,1/),eig(:neig),(/1,1,1/))
        ELSE
            IF (present(eig)) CALL juDFT_error("BUG in calling write_eig")
        ENDIF
        IF (present(z).AND..NOT.present(neig))&
            &    CALL juDFT_error("BUG in calling write_eig with eigenvector")

        n1=1;n2=0
        IF (present(n_size)) n1=n_size
        IF (present(n_rank)) n2=n_rank
        IF (present(z)) THEN

            SELECT TYPE(z)
                TYPE IS (real)
                   CALL io_write_real2s(&
                &                     d%evsetid,(/1,1,n2+1,nk,jspin/),&
                &           (/1,nmat,neig,1,1/),real(z(:nmat,:neig)),(/1,1,n1,1,1/))
                TYPE IS (complex)
                   CALL io_write_real2s(&
                &                     d%evsetid,(/1,1,n2+1,nk,jspin/),&
                &           (/1,nmat,neig,1,1/),real(z(:nmat,:neig)),(/1,1,n1,1,1/))
                   CALL io_write_real2s(&
                    &                     d%evsetid,(/2,1,n2+1,nk,jspin/),&
                    &           (/1,nmat,neig,1,1/),aimag(z(:nmat,:neig)),&
                    &           (/1,1,n1,1,1/))
            END SELECT
        ENDIF

#endif
    END SUBROUTINE write_eig

#ifdef CPP_HDF

    !----------------------------------------------------------------------
    SUBROUTINE priv_r_vecc(&
        &                     d,nk,jspin,n_start,n_end,&
        &                     nmat,z)

        USE m_hdf_tools
        IMPLICIT NONE
        TYPE(t_data_HDF),INTENT(IN)::d
        INTEGER, INTENT(IN)  :: nk,jspin
        INTEGER, INTENT(IN)  :: n_start,n_end
        INTEGER, INTENT(OUT) :: nmat
        COMPLEX, INTENT(OUT) :: z(:,:)

        REAL, ALLOCATABLE :: z1(:,:,:)
        INTEGER i,j,neig_l

        neig_l = n_end - n_start + 1

        ! read matrix size
        CALL io_read_integer0(&
            &                      d%nmatsetid,(/nk,jspin/),(/1,1/),&
            &                                                nmat)

        IF ( nmat > size(z,1) .OR. neig_l > size(z,2) ) THEN
            WRITE (6,*) nmat,size(z,1),size(z,2)
            CALL juDFT_error("eig66_hdf$read_vec",calledby ="eig66_hdf")
        ENDIF

        ! read eigenvectors
        ALLOCATE (z1(2,nmat,neig_l))
        CALL io_read_real3(d%evsetid,(/1,1,n_start,nk,jspin/),&
            &                      (/2,nmat,neig_l,1,1/),z1)
      
        DO i=1,neig_l
            DO j=1,nmat
                z(j,i) = cmplx( z1(1,j,i) ,z1(2,j,i) )
            ENDDO
        ENDDO

        DEALLOCATE (z1)

    END SUBROUTINE priv_r_vecc
    !-----------------------------------------------------------------------

#endif

    SUBROUTINE read_eig(id,nk,jspin,nv,nmat,k1,k2,k3,bk,wk,neig,eig,el,&
        &            ello,evac,kveclo,n_start,n_end,z)
        IMPLICIT NONE
        INTEGER, INTENT(IN)            :: id,nk,jspin
        INTEGER, INTENT(OUT),OPTIONAL  :: nv,nmat
        INTEGER, INTENT(OUT),OPTIONAL  :: neig
        REAL,    INTENT(OUT),OPTIONAL  :: eig(:)
        INTEGER, INTENT(OUT),OPTIONAL  :: k1(:),k2(:),k3(:),kveclo(:)
        REAL,    INTENT(OUT),OPTIONAL  :: evac(:),ello(:,:),el(:,:)
        REAL,    INTENT(OUT),OPTIONAL  :: bk(:),wk
        INTEGER, INTENT(IN),OPTIONAL   :: n_start,n_end
        CLASS(*),TARGET,INTENT(OUT),OPTIONAL  :: z(:,:)

#ifdef CPP_HDF
        INTEGER:: n1,n,k,k1_t,k2_t,k3_t
        TYPE(t_data_HDF),pointer::d
        call priv_find_data(id,d)


        IF (present(neig))  THEN
            CALL io_read_integer0(d%neigsetid,(/nk,jspin/),(/1,1/),neig)

            IF ( present(eig) ) THEN                           ! read eigenv
                IF ( neig > size(eig) ) THEN
                    WRITE(*,*) neig,size(eig)
                    CALL juDFT_error("eig66_hdf$readeig",calledby ="eig66_hdf")
                ENDIF
                CALL io_read_real1(d%energysetid,(/1,nk,jspin/),(/neig,1,1/),&
                    &                      eig(:neig))
            ENDIF
        ENDIF

        IF (present(k1)) THEN
            IF (.NOT.(present(k2).AND.present(k3).AND.present(kveclo)))&
                &    CALL juDFT_error("BUG1 in calling read_eig")

            CALL io_read_integer0(d%nvsetid,(/nk,jspin/),(/1,1/),n1)
            IF (n1>size(k1))  CALL juDFT_error("eig66_hdf$read_basis",&
                &     calledby="eig66_hdf")
            !read basis
            CALL io_read_integer1(d%ksetid,(/1,1,nk,jspin/),(/n1,1,1,1/),k1(:n1))
            CALL io_read_integer1(d%ksetid,(/1,2,nk,jspin/),(/n1,1,1,1/),k2(:n1))
            CALL io_read_integer1(d%ksetid,(/1,3,nk,jspin/),(/n1,1,1,1/),k3(:n1))
            DO k = 1, size(kveclo)
                CALL io_read_integer0(d%ksetid,(/n1+k,1,nk,jspin/),(/1,1,1,1/),k1_t)
                CALL io_read_integer0(d%ksetid,(/n1+k,2,nk,jspin/),(/1,1,1,1/),k2_t)
                CALL io_read_integer0(d%ksetid,(/n1+k,3,nk,jspin/),(/1,1,1,1/),k3_t)
                DO n = 1, n1
                    IF (( (k1_t == k1(n)).AND.(k2_t == k2(n)) ).AND.(k3_t == k3(n)) ) THEN
                        kveclo(k) = n
                        CYCLE
                    ENDIF
                ENDDO
            ENDDO
            IF (present(nv)) nv=n1
        ELSE
            IF (present(nv)) CALL io_read_integer0(d%nvsetid,(/nk,jspin/),(/1,1/),nv)

        ENDIF
        IF (present(nmat)) &
            & CALL io_read_integer0(d%nmatsetid,(/nk,jspin/),(/1,1/),nmat)
        IF (present(el)) CALL io_read_real2(d%esetid,(/1,1,nk,jspin/),&
            &                   (/size(el,1),size(el,2),1,1/),el(:,:))
        IF (present(ello)) CALL io_read_real2(d%ellosetid,(/1,1,nk,jspin/),&
            &                   (/size(ello,1),size(ello,2),1,1/),ello(:,:))
        IF (present(evac)) CALL io_read_real1(d%evacsetid,(/1,nk,jspin/),&
            &                 (/2,1,1/),evac)

        IF (present(bk)) CALL&
             io_read_real1(d%bksetid,(/1,nk,jspin/),(/3,1,1/),bk)
        IF (present(wk)) CALL&
             io_read_real0(d%wksetid,(/nk,jspin/),(/1,1/),wk)

        IF (present(n_start)) THEN
            IF (.NOT.present(n_end)) CALL juDFT_error("BUG3 in read_eig")
            IF (present(z)) THEN
               SELECT TYPE(z)
                 TYPE IS (real)
                       CALL priv_r_vec(d,nk,jspin,n_start,n_end,n1,z)
                 TYPE is (complex)
                       CALL priv_r_vecc(d,nk,jspin,n_start,n_end,n1,z)
               END SELECT
            ENDIF
            IF (present(nmat)) nmat=n1
        ENDIF
#endif
    END SUBROUTINE

END MODULE

