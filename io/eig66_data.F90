!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_eig66_data
#include "juDFT_env.h"
#ifdef CPP_HDF
    use hdf5
#endif
    implicit none

    TYPE :: t_data
       INTEGER:: io_mode
       INTEGER:: jspins,nkpts,nmat,neig,lmax,nlotot,nlo,ntype
       LOGICAL:: l_real,l_soc
    END TYPE

    TYPE,EXTENDS(t_data):: t_data_DA
        INTEGER            :: recl_vec=0,recl_wiks
        CHARACTER(LEN=20)  :: fname="eig"
        INTEGER            :: file_io_id_vec,file_io_id_wiks
    END TYPE

    TYPE,extends(t_data):: t_data_MPI
       INTEGER             :: n_size=1
       INTEGER             :: size_k,size_el,size_ello,size_eig
       INTEGER             :: eig_handle,zr_handle,zc_handle,neig_handle,w_iks_handle
       INTEGER,ALLOCATABLE :: pe_basis(:,:),slot_basis(:,:)
       INTEGER,ALLOCATABLE :: pe_ev(:,:,:),slot_ev(:,:,:)
       INTEGER             :: irank
       INTEGER,POINTER     :: neig_data(:)
       REAL,POINTER        :: eig_data(:),zr_data(:), w_iks_data(:)
       COMPLEX,POINTER     :: zc_data(:)
    END TYPE
    TYPE,EXTENDS(t_data):: t_data_hdf
#ifdef CPP_HDF
         INTEGER(HID_T) :: fid
         INTEGER(HID_T) :: neigsetid
         INTEGER(HID_T) :: energysetid,wikssetid,evsetid
         CHARACTER(LEN=20) :: fname="eig"
#endif
      END TYPE

   TYPE,EXTENDS(t_data):: t_data_mem
        INTEGER,ALLOCATABLE :: eig_int(:)
        REAL,ALLOCATABLE    :: eig_eig(:,:,:)
        REAL,ALLOCATABLE    :: eig_vecr(:,:)
        COMPLEX,ALLOCATABLE :: eig_vecc(:,:)
    END TYPE

    TYPE t_list
        INTEGER               :: id
        CLASS(t_data),POINTER :: data
        TYPE(t_list),POINTER  :: next=>null()
    END TYPE



    TYPE(t_list),POINTER :: linked_list=>null()
    private linked_list
    INTEGER, PARAMETER :: DA_mode=0,HDF_mode=1,MEM_mode=2,MPI_mode=3

    contains
    
    subroutine eig66_data_storedefault(d,jspins,nkpts,nmat,neig,lmax,nlotot,nlo,ntype,l_real,l_soc)
    CLASS(t_data)::d
    INTEGER,INTENT(IN)::jspins,nkpts,nmat,neig,lmax,nlotot,nlo,ntype
    LOGICAL,INTENT(IN):: l_real,l_soc
    d%jspins=jspins
    d%nkpts=nkpts
    d%nmat=nmat
    d%neig=neig
    d%lmax=lmax
    d%nlotot=nlotot
    d%nlo=nlo
    d%ntype=ntype
    d%l_real=l_real
    d%l_soc=l_soc
    END SUBROUTINE

    subroutine eig66_find_data(d,id,io_mode)
    IMPLICIT NONE
    INTEGER,INTENT(IN) ::id
    INTEGER,INTENT(IN),OPTIONAL :: io_mode
    CLASS(t_data),pointer::d

    TYPE(t_list),POINTER:: listpointer,lastinlist
    lastinlist=>null()

    listpointer=>linked_list
    DO WHILE(associated(listpointer))
       lastinlist=>listpointer
       if (listpointer%id==id) THEN
           d=>listpointer%data
           return
       endif
       listpointer=>listpointer%next
    enddo
    !no pointer found
    IF (present(io_mode)) THEN
       IF (.not.associated(lastinlist)) THEN
          allocate(linked_list)
          linked_list%id=id
          lastinlist=>linked_list
       ELSE
          allocate(lastinlist%next)
          lastinlist%next%id=id
          lastinlist=>lastinlist%next
       ENDIF
       SELECT CASE (io_mode)
         case (DA_MODE)
          allocate(t_data_DA::lastinlist%data)
         case (HDF_MODE)
#ifdef CPP_HDF
          allocate(t_data_HDF::lastinlist%data)
#else
          call juDFT_error("Cannot use hdf mode for IO, recompile with CPP_HDF",calledby="eig66_data")
#endif
         case (MEM_MODE)
          allocate(t_data_MEM::lastinlist%data)
         case (MPI_MODE)
          allocate(t_data_MPI::lastinlist%data)
         end select
         lastinlist%data%io_mode=io_mode
         d=>lastinlist%data
    ELSE
       call juDFT_error("BUG:Could not find data object in eig66_mpi")
    ENDIF
    END SUBROUTINE

    subroutine eig66_remove_data(id)
    INTEGER,INTENT(IN)::id

    TYPE(t_list),POINTER:: listpointer,lastpointer
    lastpointer=>null()
    listpointer=>linked_list
    loop:DO WHILE(associated(listpointer))
      IF (listpointer%id==id) THEN
        exit loop
      ENDIF
      lastpointer=>listpointer
      listpointer=>listpointer%next
    ENDDO loop

    if (.not.associated(listpointer)) call juDFT_error("BUG in eig66_data: ID not found in deleting")
    IF (associated(lastpointer)) THEN
        lastpointer%next=>listpointer%next
    ELSE
        linked_list=>listpointer%next
    ENDIF

    deallocate(listpointer)
    end subroutine

    INTEGER FUNCTION eig66_data_newid(mode)
    INTEGER,INTENT(IN)  :: mode
    TYPE(t_list),POINTER:: listpointer
    INTEGER             :: id
    CLASS(t_data),POINTER::d

    id=0
    listpointer=>linked_list
    DO WHILE(associated(listpointer))
      id=max(id,listpointer%id)
      listpointer=>listpointer%next
    ENDDO
    eig66_data_newid=id+1

    call eig66_find_data(d,id+1,mode)

    end function

    INTEGER function eig66_data_mode(id)RESULT(mode)
    INTEGER,INTENT(IN)  :: id
    TYPE(t_list),POINTER:: listpointer
    mode=-1
    listpointer=>linked_list
    DO WHILE(associated(listpointer))
      if (id==listpointer%id) THEN
           mode=listpointer%data%io_mode
           return
      ENDIF
      listpointer=>listpointer%next
    ENDDO
    END FUNCTION

end module m_eig66_data
