!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_eig66_da
#include "juDFT_env.h"
  ! Do the IO of the eig-file in fortran direct-access
  ! The eig-file is split into two parts:
  ! eig.bas contains the basis-set information
  ! eig.vec contains the eigenvalues and the eigenvectors
  ! The record number is given by nrec=nk+(jspin-1)*nkpts
  ! each record contains:
  ! eig.bas: el,evac,ello,bkpt,wtkpt,nv,nmat
  ! eig.vec: ne,eig,z**
  !**: real or complex depending on calculation type
  USE m_eig66_data
  USE m_types
  IMPLICIT NONE

CONTAINS
  SUBROUTINE priv_find_data(id,d)
    INTEGER,INTENT(IN)            :: id
    TYPE(t_data_DA),POINTER,INTENT(out)   :: d

    CLASS(t_data),POINTER   ::dp
    CALL eig66_find_data(dp,id)
    SELECT TYPE(dp)
    TYPE is (t_data_da)
       d=>dp
       CLASS default
       CALL judft_error("BUG: wrong datatype in eig66_da")
    END SELECT
  END SUBROUTINE priv_find_data

  SUBROUTINE open_eig(id,nmat,neig,nkpts,jspins,lmax,nlo,ntype,nlotot,create,l_real,l_soc,l_dos,l_mcd,l_orb,filename,layers,nstars,ncored,nsld,nat)
    INTEGER, INTENT(IN) :: id,nmat,neig,nkpts,jspins,nlo,ntype,lmax,nlotot
    LOGICAL, INTENT(IN) :: create,l_real,l_soc
    LOGICAL,INTENT(IN),OPTIONAL ::  l_dos,l_mcd,l_orb
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: filename
    INTEGER,INTENT(IN),OPTIONAL :: layers,nstars,ncored,nsld,nat
    !locals
    LOGICAL :: l_file
    INTEGER :: i1,recl_z,recl_eig,recl_dos
    REAL    :: r1,r3(3)
    COMPLEX :: c1
    TYPE(t_data_DA),POINTER:: d

    CALL priv_find_data(id,d)

    IF (PRESENT(filename)) d%fname=filename
    CALL eig66_data_storedefault(d,jspins,nkpts,nmat,neig,lmax,nlotot,nlo,ntype,l_real,l_soc,l_dos,l_mcd,l_orb)

    !Calculate the record length

    INQUIRE(IOLENGTH=recl_eig) r1
    d%recl_wiks=recl_eig*neig
    
    recl_eig=recl_eig*(neig+2) ! add a 2 for integer 'neig'
    if (l_real.and..not.l_soc ) THEN
       INQUIRE(IOLENGTH=recl_z) r1
    else
       INQUIRE(IOLENGTH=recl_z) c1
    endif
    recl_z=recl_z*nmat*neig
    
    d%recl_vec=recl_eig+recl_z

    IF (d%l_dos) THEN
       IF (.NOT.(PRESENT(layers).AND.PRESENT(nstars).AND.PRESENT(ncored).AND.PRESENT(nsld).AND.PRESENT(nat))) &
            CALL judft_error("BUG:Could not open file for DOS-data",calledby="eig66_da")
       INQUIRE(IOLENGTH=i1) i1
       recl_dos=i1*2*neig !ksym&jsym
       INQUIRE(IOLENGTH=i1) r1
       recl_dos=recl_dos+i1*3*neig !qvac&qis
       recl_dos=recl_dos+i1*4*ntype*neig !qal
       recl_dos=recl_dos+i1*neig*2*max(1,layers) !qvlay
       IF (l_orb) THEN
          recl_dos=recl_dos+i1*2*nsld*neig !qintsl,qmtsl
          recl_dos=recl_dos+i1*24*neig*nat !qmtp,orbcomp
       ENDIF
       INQUIRE(IOLENGTH=i1) c1
       recl_dos=recl_dos+i1*nstars*neig*max(1,layers)*2 !qstars
       IF (l_mcd) recl_dos=recl_dos+i1*3*ntype*ncored*neig !mcd
    ELSE
       recl_dos=-1
    ENDIF
    d%recl_dos=recl_dos


    IF (create) THEN
       INQUIRE(file=TRIM(d%fname),opened=l_file)
       DO WHILE(l_file)
          write(*,*) "eig66_open_da:",d%fname," in use"
          d%fname=TRIM(d%fname)//"6"
          INQUIRE(file=TRIM(d%fname),opened=l_file)
       ENDDO
       d%file_io_id_vec=priv_free_uid()
       OPEN(d%file_io_id_vec,FILE=TRIM(d%fname),ACCESS='direct',FORM='unformatted',RECL=d%recl_vec,STATUS='unknown')
       d%file_io_id_wiks=priv_free_uid()
       OPEN(d%file_io_id_wiks,FILE=TRIM(d%fname)//".wiks",ACCESS='direct',FORM='unformatted',RECL=d%recl_wiks,STATUS='unknown')
       IF(d%recl_dos>0) THEN
          d%file_io_id_dos=priv_free_uid()
          OPEN(d%file_io_id_dos,FILE=TRIM(d%fname)//".dos",ACCESS='direct',FORM='unformatted',RECL=d%recl_dos,STATUS='unknown')
       ENDIF

    ELSE
       d%file_io_id_vec=priv_free_uid()
       OPEN(d%file_io_id_vec,FILE=TRIM(d%fname),ACCESS='direct',FORM='unformatted',RECL=d%recl_vec,STATUS='old')
       d%file_io_id_wiks=priv_free_uid()
       OPEN(d%file_io_id_wiks,FILE=TRIM(d%fname)//".wiks",ACCESS='direct',FORM='unformatted',RECL=d%recl_wiks,STATUS='old')
       IF(d%recl_dos>0) THEN
          d%file_io_id_dos=priv_free_uid()
          OPEN(d%file_io_id_dos,FILE=TRIM(d%fname)//".dos",ACCESS='direct',FORM='unformatted',RECL=d%recl_dos,STATUS='old')
       ENDIF
    ENDIF
  CONTAINS
    INTEGER FUNCTION priv_free_uid() RESULT(uid)
      IMPLICIT NONE
      LOGICAL::used
      used=.TRUE.
      uid=665
      DO WHILE(used)
         uid=uid+1
         INQUIRE(UNIT=uid,OPENED=used)
      END DO
    END FUNCTION priv_free_uid
  END SUBROUTINE open_eig
  SUBROUTINE close_eig(id,filename)
    INTEGER,INTENT(IN)::id
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: filename
    TYPE(t_data_DA),POINTER:: d

    CALL priv_find_data(id,d)

    CLOSE(d%file_io_id_vec)
    CLOSE(d%file_io_id_wiks)
    d%recl_vec=0
    d%recl_wiks=0

    !If a filename was given and the name is not the current filename then rename
    IF (PRESENT(filename)) THEN
       IF (filename.NE.d%fname) THEN
          CALL system("mv "//TRIM(d%fname)//" "//TRIM(filename))
       ENDIF
    ENDIF
    d%fname="eig"
    CALL eig66_remove_data(id)
  END SUBROUTINE close_eig
  SUBROUTINE read_eig(id,nk,jspin,neig,eig,w_iks,n_start,n_end,zmat)
    IMPLICIT NONE
    INTEGER, INTENT(IN)            :: id,nk,jspin
    INTEGER, INTENT(OUT),OPTIONAL  :: neig
    REAL,    INTENT(OUT),OPTIONAL  :: eig(:),w_iks(:)
    INTEGER, INTENT(IN),OPTIONAL   :: n_start,n_end
    TYPE(t_zmat),OPTIONAL  :: zmat

    !Local variables
    INTEGER:: nv_s,nmat_s,n,nrec,neig_s
    REAL   :: bkpt(3),wtkpt
    REAL,ALLOCATABLE::eig_s(:),zr_s(:,:)
    COMPLEX,ALLOCATABLE::zc_s(:,:)
    TYPE(t_data_DA),POINTER:: d



    CALL priv_find_data(id,d)
    ! check if io is performed correctly
    IF (PRESENT(n_start)) THEN
       IF (n_start/=1) &
            CALL juDFT_error("In direct access mode only all eigenstates can be read")
    ENDIF

    nrec=nk+(jspin-1)*d%nkpts
  
    IF (PRESENT(w_iks)) THEN
       print *, "R:w_iks:",nrec
        read(d%file_io_id_wiks,REC=nrec) w_iks
    ENDIF
  
    
    IF (.NOT.(PRESENT(eig).OR.PRESENT(neig).OR.PRESENT(zmat))) RETURN
    READ(d%file_io_id_vec,REC=nrec) neig_s
    IF (PRESENT(neig)) THEN
       neig=neig_s
    ENDIF
    IF (.NOT.(PRESENT(eig).OR.PRESENT(zmat))) RETURN
    ALLOCATE(eig_s(neig_s))
    IF (PRESENT(zmat)) THEN
       IF (zmat%l_real) THEN
          INQUIRE(IOLENGTH=n) neig_s,eig_s,REAL(zmat%z_r)
          IF (n>d%recl_vec) THEN
             CALL juDFT_error("BUG: Too long record")
          END IF
          READ(d%file_io_id_vec,REC=nrec) neig_s,eig_s,zmat%z_r
       ELSE
          INQUIRE(IOLENGTH=n) neig_s,eig_s,CMPLX(zmat%z_c)
          IF (n>d%recl_vec) THEN
             CALL juDFT_error("BUG: Too long record")
          END IF
          READ(d%file_io_id_vec,REC=nrec) neig_s,eig_s,zmat%z_c
       ENDIF
    ELSE
       INQUIRE(IOLENGTH=n) neig_s,eig_s
       IF (n>d%recl_vec) CALL juDFT_error("BUG: Too long record")
       READ(d%file_io_id_vec,REC=nrec) neig_s,eig_s
    ENDIF
    IF (PRESENT(eig)) eig(:min(size(eig),neig_s))=eig_s(:min(size(eig),neig_s))
   
  END SUBROUTINE read_eig

  SUBROUTINE write_eig(id,nk,jspin,neig,neig_total,eig,w_iks,n_size,n_rank,zmat)
    INTEGER, INTENT(IN)          :: id,nk,jspin
    INTEGER, INTENT(IN),OPTIONAL :: n_size,n_rank
    INTEGER, INTENT(IN),OPTIONAL :: neig,neig_total
    REAL,    INTENT(IN),OPTIONAL :: eig(:),w_iks(:)
    TYPE(t_mat),INTENT(IN),OPTIONAL :: zmat

    INTEGER:: nrec,r_len
    INTEGER:: nv_s,nmat_s
    REAL   :: bkpt(3),wtkpt
    TYPE(t_data_DA),POINTER:: d

    CALL priv_find_data(id,d)
    !This mode requires all data to be written at once!!

    IF (PRESENT(n_size).AND.PRESENT(n_rank)) THEN
       IF (n_size/=1.OR.n_rank/=0) &
            CALL juDFT_error("Direct Access IO not possible in eigenvalue parallel code")
    ENDIF
    !check record length
    !INQUIRE(iolength=r_len) nmat,el,evac,ello,bk,wk,nv,d%kvec_s,kveclo
    !if (r_len>recl_bas) call juDFT_error("BUG: too long record")

    !Now it is time for the IO :-)
    nrec=nk+(jspin-1)*d%nkpts
    IF (PRESENT(neig).AND.PRESENT(neig_total)) THEN
       IF (neig.NE.neig_total) THEN
          CALL juDFT_error("Neig and neig_total have to be equal in DA mode",calledby="eig66_da")
       ENDIF
    ENDIF
    IF (PRESENT(w_iks)) THEN
       write(d%file_io_id_wiks,REC=nrec) w_iks
    ENDIF
 

    IF (.NOT.PRESENT(eig).OR..NOT.PRESENT(neig)) RETURN
    !Now the IO of the eigenvalues/vectors
    IF (PRESENT(zmat)) THEN
       IF (zmat%l_real) THEN
          INQUIRE(IOLENGTH=r_len) neig,eig,REAL(zmat%data_r)
          IF (r_len>d%recl_vec) CALL juDFT_error("BUG: too long record")
          WRITE(d%file_io_id_vec,REC=nrec) neig,eig,REAL(zmat%data_r)
       ELSE
          INQUIRE(IOLENGTH=r_len) neig,eig(:neig),CMPLX(zmat%data_c)
          IF (r_len>d%recl_vec) CALL juDFT_error("BUG: too long record")
          WRITE(d%file_io_id_vec,REC=nrec) neig,eig(:neig),CMPLX(zmat%data_c)
       ENDIF
    ELSE
       INQUIRE(IOLENGTH=r_len) neig,eig
       IF (r_len>d%recl_vec) CALL juDFT_error("BUG: too long record")
       WRITE(d%file_io_id_vec,REC=nrec) neig,eig
    ENDIF

  END SUBROUTINE write_eig


  SUBROUTINE write_dos(id,nk,jspin,qal,qvac,qis,qvlay,qstars,ksym,jsym,mcd,qintsl,qmtsl,qmtp,orbcomp)
    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: id,nk,jspin
    REAL,INTENT(IN)              :: qal(:,:,:),qvac(:,:),qis(:),qvlay(:,:,:)
    COMPLEX,INTENT(IN)           :: qstars(:,:,:,:)
    INTEGER,INTENT(IN)           :: ksym(:),jsym(:)
    REAL,INTENT(IN),OPTIONAL  :: mcd(:,:,:)
    REAL,INTENT(IN),OPTIONAL     :: qintsl(:,:),qmtsl(:,:),qmtp(:,:),orbcomp(:,:,:)
    TYPE(t_data_DA),POINTER:: d
    INTEGER:: nrec
    CALL priv_find_data(id,d)

    nrec=nk+(jspin-1)*d%nkpts

    IF (d%l_orb.AND.PRESENT(qmtsl)) THEN
       IF (d%l_mcd) CPP_error("mcd & orbital decomposition not implemented in IO")
       WRITE(d%file_io_id_dos,REC=nrec) qal,qvac,qis,qvlay,qstars,ksym,jsym,qintsl,qmtsl,qmtp,orbcomp
    ELSEIF(d%l_mcd.AND.PRESENT(mcd)) THEN
       WRITE(d%file_io_id_dos,REC=nrec) qal,qvac,qis,qvlay,qstars,ksym,jsym,mcd
    ELSE
       WRITE(d%file_io_id_dos,REC=nrec) qal,qvac,qis,qvlay,qstars,ksym,jsym
    END IF
  END SUBROUTINE write_dos

  SUBROUTINE read_dos(id,nk,jspin,qal,qvac,qis,qvlay,qstars,ksym,jsym,mcd,qintsl,qmtsl,qmtp,orbcomp)
    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: id,nk,jspin
    REAL,INTENT(OUT)              :: qal(:,:,:),qvac(:,:),qis(:),qvlay(:,:,:)
    COMPLEX,INTENT(OUT)           :: qstars(:,:,:,:)
    INTEGER,INTENT(OUT)           :: ksym(:),jsym(:)
    REAL,INTENT(OUT),OPTIONAL  :: mcd(:,:,:)
    REAL,INTENT(OUT),OPTIONAL     :: qintsl(:,:),qmtsl(:,:),qmtp(:,:),orbcomp(:,:,:)
    TYPE(t_data_DA),POINTER:: d
    INTEGER:: nrec
    CALL priv_find_data(id,d)

    nrec=nk+(jspin-1)*d%nkpts

    IF (d%l_orb.AND.PRESENT(qmtsl)) THEN
       IF (d%l_mcd) CPP_error("mcd & orbital decomposition not implemented in IO")
       READ(d%file_io_id_dos,REC=nrec) qal,qvac,qis,qvlay,qstars,ksym,jsym,qintsl,qmtsl,qmtp,orbcomp
    ELSEIF(d%l_mcd.AND.PRESENT(mcd)) THEN
       READ(d%file_io_id_dos,REC=nrec) qal,qvac,qis,qvlay,qstars,ksym,jsym,mcd
    ELSE
       READ(d%file_io_id_dos,REC=nrec) qal,qvac,qis,qvlay,qstars,ksym,jsym
    END IF
  END SUBROUTINE read_dos


END MODULE m_eig66_da
