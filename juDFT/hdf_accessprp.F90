module m_hdf_accessprp
    USE hdf5
#include "juDFT_env.h"
    implicit none
    private
    !the hdf-access-properties
    integer          :: n_access_prp=0
    integer(hid_t)   :: access_prp(100)
    character(len=10) :: access_filename(100)
    public hdf_access_prp
    CONTAINS

    recursive SUBROUTINE  priv_generate_access_prp(setupfile)
      !Subroutine reads the two files ~/.gf_hdf and .gf_hdf
      !to generate access properties
      !typical file:
      !&hdf filename="embpot11",driver="mpiposix",alignment=103333
      !&hdf filename="default",driver="mpiio",hint(1)="romio_ds_read",value(1)="disable"
      IMPLICIT NONE
      character(len=*),optional::setupfile
      !generate the access_prp for later use
      character(len=10) :: filename
      character(len=9)  :: driver
      character(len=20) :: hint(10),value(10)
      logical           :: keep,gpfs,l_exist
      integer(hsize_t)  :: mem_increment
      integer(hsize_t)  :: alignment
      NAMELIST /hdf/filename,driver,hint,value,keep,gpfs,mem_increment,alignment

      integer :: n,i,hdferr,ierr,info
      character(len=128)::path
#ifdef CPP_MPI
      INCLUDE 'mpif.h'
#endif

      IF (.not.present(setupfile)) THEN
         CPP_juDFT_timestart_debug("generating access prp")
         n_access_prp=1
         call priv_generate_access_prp("default")
         call getenv("HOME",path)
         call priv_generate_access_prp(trim(path)//"/.gf_hdf")
         call priv_generate_access_prp("gf_hdf")
         CPP_juDFT_timestop_debug("generating access prp")
         return
      ENDIF

      if (.not.(setupfile=="default")) THEN
           INQUIRE(file=setupfile,exist=l_exist)
           IF (.not.l_exist) return
           OPEN(999,file=setupfile)
      ENDIF
      write(6,*) "Access properties from:",setupfile

      n=0
      readloop:DO
        filename="default"
#ifdef CPP_MPI
        driver="mpiio"
#else
        driver="default"
#endif
        hint=""
        value=""
        keep=.false.
        mem_increment=1024*1024
        alignment=-1
        gpfs=.false.

        if (.not.(setupfile=="default")) read(999,hdf,end=100)
        if (n>0.and.(setupfile=="default")) return
        !Now create the access property
        if (filename=="default") THEN
             n=1
        else
             n_access_prp=n_access_prp+1
             n=n_access_prp
             access_filename(n)=filename
        endif

        write(6,"(a,i5,8(1x,a))") "Access_prp:",n,"for:",filename,"with:",driver

        !different drivers
        if (index(driver,"default")==1) THEN
            access_prp(n)=H5P_DEFAULT_f
            cycle readloop
        endif
        call  h5pcreate_f(H5P_FILE_ACCESS_F, access_prp(n), hdferr)
        IF (index(driver,"core")==1) THEN
            CALL h5pset_fapl_core_f(access_prp(n), mem_increment, keep,hdferr)
            cycle readloop
        ENDIF
#ifdef CPP_MPI
        CALL MPI_BARRIER(MPI_COMM_WORLD,hdferr)
        IF (index(driver,"mpiio")==1) THEN
            !create info object
            CALL MPI_Info_create(info,ierr);
            i=1
            DO WHILE(i<=size(hint))
              if (len(trim(hint(i)))>1) THEN
                  write(6,*) "Hint:",hint(i),"=",value(i)
                  call MPI_Info_set(info,hint(i),value(i),hdferr)
                  i=i+1
              else
                  i=size(hint)+1
              endif
            ENDDO
            CALL h5pset_fapl_mpio_f(access_prp(n), MPI_COMM_WORLD,INFO,hdferr)
            CALL MPI_BARRIER(MPI_COMM_WORLD,hdferr)
            call mpi_info_free(info,ierr)
            if (alignment>0) CALL h5pset_alignment_f(access_prp(n), INT(0,hsize_t),alignment, hdferr)
            cycle readloop
        ENDIF
        IF (index(driver,"mpiposix")==1) THEN
            call judft_error("MPIPOSIX driver not implemented")
             !CALL h5pset_fapl_mpiposix_f(access_prp(n), MPI_COMM_WORLD,gpfs,hdferr)
            if (alignment>0) CALL h5pset_alignment_f(access_prp(n), INT(0,hsize_t),alignment, hdferr)
            cycle readloop
        ENDIF
#endif
        write(0,*) "Driver name unkown:",driver
        call judft_error("Unkown driver",calledby="gf_io2dmat")
      ENDDO readloop
100   close(999)

      END SUBROUTINE

      FUNCTION hdf_access_prp(filename)
      !return the access_prp from the list
      character(len=*),intent(in) :: filename
      INTEGER(hid_t)          :: hdf_access_prp
      INTEGER                 :: n
      CPP_juDFT_timestart_debug("getting access prp")
      !if this is the first attempt to get an access_prp, generate them
      if (n_access_prp==0) CALL priv_generate_access_prp()
      hdf_access_prp=access_prp(1)
      DO n=2,n_access_prp
          if (trim(filename)==trim(access_filename(n))) THEN
                  hdf_access_prp=access_prp(n)
                  write(6,*) "Assigned:",n," to ", filename
          ENDIF
      ENDDO
      CPP_juDFT_timestop_debug("getting access prp")
      END function

end module m_hdf_accessprp
