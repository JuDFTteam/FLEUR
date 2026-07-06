!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_writetrans
#include "juDFT_env.h"
      USE hdf5
      IMPLICIT NONE
      INTEGER(hid_t)::fid,varid
      LOGICAL       ::hdf=.false.
      ! This module writes the transmitted current to either a simple text-file or an hdf file
      PRIVATE 


      PUBLIC writetrans,gf_writetrans_hdfclose,gf_writetrans_hdfopen
      CONTAINS 

      SUBROUTINE gf_writetrans_hdf(en,nk,jspin,j)
      USE m_hdf_tools
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: en,nk,jspin
      REAL,INTENT(IN)    :: j(:)

      CALL io_WRITE(varid,(/jspin,nk,en,1/),(/1,1,1,size(j)/),j)

      END SUBROUTINE

      SUBROUTINE gf_writetrans_hdfclose()
      USE m_hdf_tools
      IMPLICIT NONE
      if (.not.hdf) return
      call io_dclose(varid)
      call io_hdfclose(fid)
      hdf=.false.
      END subroutine

      SUBROUTINE gf_writetrans_hdfopen(mpi,kpts,sym,num_layers,jspins)
      USE m_hdf_tools
      USE m_gf_types
      USE m_gf_energies
      IMPLICIT NONE
      TYPE(t_mpi),INTENT(IN) :: mpi
      TYPE(t_kpts),INTENT(IN):: kpts
      TYPE(t_sym),INTENT(IN) :: sym
      INTEGER    ,INTENT(IN) :: num_layers,jspins


      real    :: bk(2),all_bk(kpts%nkpts,sym%nop,2)
      INTEGER :: nk,nkpt(kpts%nkpts),n,i,hdferr
      INTEGER(hid_t)::access_prp
      logical :: new
      if (mpi%kl_layers(1).ne.1) return !only open for PE that work on layer==1

#ifdef CPP_HDFMPI
      call h5pcreate_f(H5P_FILE_ACCESS_F, access_prp, hdferr)
      CALL h5pset_fapl_mpiposix_f(access_prp, mpi%samelayer_subcom,.false.,hdferr)
#else
      access_prp=H5P_DEFAULT_f
#endif
      CALL h5fcreate_f('transcurrent.hdf',H5F_ACC_TRUNC_F,fid, hdferr,H5P_DEFAULT_f,access_prp)

      !Create data-sets

      !First store all information

      DO nk=1,kpts%nkpts
           nkpt(nk) = 0
           DO n = 1,sym%nop
                 bk = MATMUL(kpts%bk(1:2,nk),1.0*sym%mrot(1:2,1:2,n))
                 new = .TRUE.
                 DO i = 1,nkpt(nk)
                      IF (ALL(bk(:) == all_bk(nk,i,:))) new=.FALSE.
                 ENDDO
                 IF (new) THEN
                    nkpt(nk)=nkpt(nk)+1
                    all_bk(nk,nkpt(nk),:)=bk
                 ENDIF
           ENDDO
      ENDDO
      CALL IO_WRITE_VAR(fid,"n_kpts",nkpt)
      CALL IO_WRITE_VAR(fid,"kpts",all_bk)
      CALL IO_WRITE_VAR(fid,"energies",real(gf_allz(1)))
      !Create variable for current
      CALL io_createvar(fid,"current",H5T_NATIVE_DOUBLE,(/jspins,kpts%nkpts,gf_noen(),5*num_layers/),varid)
      hdf=.true.
      END subroutine


      !<--S:check()                                                     
      SUBROUTINE check() 
      IMPLICIT NONE 
                                !checks if file is already opened, else 
      LOGICAL::l_opened 
      INQUIRE(FILE='transcurrent',OPENED=l_opened) 
      IF (.NOT.l_opened) THEN 
         OPEN(177,FILE='transcurrent',FORM='formatted',STATUS='replace') 
      ENDIF 
      END SUBROUTINE check 
      !>                                                                
      !<--S:writetrans(en_i,nk_i,jspin,bkpts,sym,cell,nochannel,j)      
      SUBROUTINE writetrans(en_i,nk_i,jspin,bkpts,sym,cell,nochannel,j,mpi)
      USE m_gf_types,ONLY :t_sym,t_cell,t_mpi
      USE m_gf_energies,ONLY:gf_z 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)      :: en_i,nk_i,jspin,nochannel
      REAL,INTENT(IN)         :: j(:)
      REAL,INTENT(IN)         :: bkpts(:,:)
      type(t_mpi),intent(in)  :: mpi
      TYPE(t_sym),INTENT(IN)  :: sym
      TYPE(t_cell),INTENT(IN) :: cell 
      !>                                                                
      !<--local variables                                               
      REAL,ALLOCATABLE  ::j_in(:),bk(:),bk2(:,:) 
      INTEGER           ::en,nk,n_curr 
      INTEGER           ::i,n,nkpt,irank,isize
      LOGICAL           ::new 
      CHARACTER(LEN = 50) :: formatstring 
#ifdef CPP_MPI                                                          
!                                                                       
!     variables for MPI                                                 
!                                                                       
#include"cpp_double.h"                                                  
      INCLUDE 'mpif.h' 
                                    !MPI error+status                   
      INTEGER::e,s(MPI_STATUS_SIZE) 
      INTEGER::rank
#endif
      !Check if we use hdf for io
      if (hdf) THEN
         CPP_juDFT_timestart("writetrans_hdf")
         call  gf_writetrans_hdf(en_i,nk_i,jspin,j)
         CPP_juDFT_timestop("writetrans_hdf")
         return
      endif
                                                                        
      !<-- make the format string                                       
      WRITE(formatstring,"(a,i0,a)") "(i1,1x,2f10.5,f10.6,",nochannel   &
     &     ,"(1x,e14.7),99(1x,f0.4))"                                   
      !>                                                                
      !>                                                                
      !<--copy DATA to non-INTENT(IN) array                             
      n_curr=size(j) 
      ALLOCATE(j_in(n_curr)) 
      j_in=j 
      nk=nk_i 
      en=en_i 
      !>                                                                
      !<--On Parallel machine data has to be send to PE=0               
#ifdef CPP_MPI                                                          
      CALL MPI_COMM_RANK(mpi%samelayer_SUBCOM,irank,E)
      CALL MPI_COMM_SIZE(mpi%samelayer_SUBCOM,isize,E)

                                                                        
      IF (irank/=0) THEN
!     for a parallel Maschine one has to send the data to PE=0          
         CALL MPI_SEND(n_curr, 1, MPI_INTEGER, 0, 1,mpi%samelayer_SUBCOM, E)
         CALL MPI_SEND(en, 1, MPI_INTEGER, 0, 2,mpi%samelayer_SUBCOM, E)
         CALL MPI_SEND(nk, 1, MPI_INTEGER, 0, 3,mpi%samelayer_SUBCOM, E)
         CALL MPI_SEND(j_in,n_curr,CPP_MPI_REAL,                        &
     &        0,4,mpi%samelayer_SUBCOM,E)
      ELSE 
         DO rank=0,isize-1
            IF (rank/=0) THEN 
               CALL MPI_RECV(n_curr,1,MPI_INTEGER,rank,1,               &
     &              mpi%samelayer_SUBCOM,S, E)
               CALL MPI_RECV(en,1,MPI_INTEGER,rank,2,                   &
     &              mpi%samelayer_SUBCOM,S, E)
               CALL MPI_RECV(nk,1,MPI_INTEGER,rank,3,                   &
     &              mpi%samelayer_SUBCOM,S, E)
              DEALLOCATE(j_in) 
              ALLOCATE(j_in(n_curr)) 
              CALL MPI_RECV(j_in,n_curr,CPP_MPI_REAL,                   &
     &             rank,4,mpi%samelayer_SUBCOM,S, E)
           ENDIF 
#endif                                                                  
           !>                                                           
                                                                        
      !<-- Now on PE = 0 write the data                                 
      !Determine all inequivalent-kpoints                               
                   !open file                                           
      CALL check() 
      ALLOCATE(bk(2),bk2(2,sym%nop)) 
      nkpt = 0 
      DO n = 1,sym%nop 
         bk = MATMUL(bkpts(1:2,nk),1.0*sym%mrot(1:2,1:2,n)) 
         new = .TRUE. 
         DO i = 1,nkpt 
            IF (ALL(bk(:) == bk2(:,i))) new=.FALSE. 
         ENDDO 
         IF (new) THEN 
            nkpt = nkpt+1 
            bk2(:,nkpt) = bk(:) 
            !Transform to cartesian koords                              
            !bk(1:2) = MATMUL(bk(1:2),cell%bmat(1:2,1:2))               
!     ok this is PE = 0 so we WRITE the DATA                            
            WRITE(177,formatstring)                                     &
     &           jspin,bk(1),bk(2),REAL(gf_z(en,0)),(j_in(i),i = 1      &
     &           ,n_curr)                                               
         ENDIF 
      ENDDO 
      DEALLOCATE(bk,bk2) 
      !>                                                                
      !<-- END MPI-loop                                                 
#ifdef CPP_MPI                                                          
         ENDDO 
      ENDIF 
#endif                                                                  
      !>                                                                
      DEALLOCATE(j_in) 
      END SUBROUTINE writetrans 
      !>                                                                
      !<-- S: trans_close()                                             
      SUBROUTINE trans_CLOSE(jspins,curr) 
!-----------------------------------------------                        
!     Close the transcurrent file                                       
!     Additionally, read in the file again, perform                     
!     some statistics and write out the total current etc               
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      IMPLICIT NONE 
      INTEGER,INTENT(IN)     :: jspins,curr 
      !<-- Locals                                                       
                                                                        
      LOGICAL::l_opened 
      !>                                                                
      INQUIRE(FILE='transcurrent',OPENED=l_opened) 
                                !nothing to do here                     
      IF (.NOT.l_opened) RETURN 
      CLOSE(177) 
                                                                        
                                                                        
      END SUBROUTINE 
      !>                                                                
      END                                           
