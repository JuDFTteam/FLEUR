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

    !Allocate the storage for the DATA always read/write
    ALLOCATE(d%el_s(0:lmax,ntype),d%ello_s(nlo,ntype),d%evac_s(2))
    ALLOCATE(d%kvec_s(nmat,3),d%kveclo_s(nlotot))
    !Calculate the record length
    INQUIRE(IOLENGTH=recl_eig) d%el_s,d%evac_s,d%ello_s,r3,r1,i1,i1,d%kvec_s,d%kveclo_s
    d%recl_bas=recl_eig
    INQUIRE(IOLENGTH=recl_eig) r1

    d%recl_wiks=recl_eig*neig
    
    print *,lmax,ntype,nlo,nlotot,nmat,neig
    
    recl_eig=recl_eig*(neig+2) ! add a 2 for integer 'neig'
    if (l_real.and..not.l_soc ) THEN
       INQUIRE(IOLENGTH=recl_z) r1
    else
       INQUIRE(IOLENGTH=recl_z) c1
    endif
    recl_z=recl_z*nmat*neig
    
    d%recl_vec=recl_eig+recl_z
    print *,l_real,l_soc
    print *,"reclen:",d%recl_vec,nmat,neig,recl_z,recl_eig

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
       d%file_io_id_bas=priv_free_uid()
       INQUIRE(file=TRIM(d%fname)//".bas",opened=l_file)
       DO WHILE(l_file)
          write(*,*) "eig66_open_da:",d%fname," in use"
          d%fname=TRIM(d%fname)//"6"
          INQUIRE(file=TRIM(d%fname)//".bas",opened=l_file)
       ENDDO
       OPEN(d%file_io_id_bas,FILE=TRIM(d%fname)//".bas",ACCESS='direct',FORM='unformatted',RECL=d%recl_bas,STATUS='unknown')
       d%file_io_id_vec=priv_free_uid()
       OPEN(d%file_io_id_vec,FILE=TRIM(d%fname)//".vec",ACCESS='direct',FORM='unformatted',RECL=d%recl_vec,STATUS='unknown')
       d%file_io_id_wiks=priv_free_uid()
       OPEN(d%file_io_id_wiks,FILE=TRIM(d%fname)//".wiks",ACCESS='direct',FORM='unformatted',RECL=d%recl_wiks,STATUS='unknown')
       IF(d%recl_dos>0) THEN
          d%file_io_id_dos=priv_free_uid()
          OPEN(d%file_io_id_dos,FILE=TRIM(d%fname)//".dos",ACCESS='direct',FORM='unformatted',RECL=d%recl_dos,STATUS='unknown')
       ENDIF

    ELSE
       d%file_io_id_bas=priv_free_uid()
       OPEN(d%file_io_id_bas,FILE=TRIM(d%fname)//".bas",ACCESS='direct',FORM='unformatted',RECL=d%recl_bas,STATUS='old')
       d%file_io_id_vec=priv_free_uid()
       OPEN(d%file_io_id_vec,FILE=TRIM(d%fname)//".vec",ACCESS='direct',FORM='unformatted',RECL=d%recl_vec,STATUS='old')
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

    DEALLOCATE(d%el_s,d%ello_s,d%evac_s,d%kvec_s,d%kveclo_s)
    CLOSE(d%file_io_id_bas)
    CLOSE(d%file_io_id_vec)
    CLOSE(d%file_io_id_wiks)
    d%recl_vec=0
    d%recl_bas=0
    d%recl_wiks=0

    !If a filename was given and the name is not the current filename then rename
    IF (PRESENT(filename)) THEN
       IF (filename.NE.d%fname) THEN
          CALL system("mv "//TRIM(d%fname)//".bas "//TRIM(filename)//".bas")
          CALL system("mv "//TRIM(d%fname)//".vec "//TRIM(filename)//".vec")
       ENDIF
    ENDIF
    d%fname="eig"
    CALL eig66_remove_data(id)
  END SUBROUTINE close_eig
  SUBROUTINE read_eig(id,nk,jspin,nv,nmat,bk,wk,neig,eig,w_iks,el,ello,evac,n_start,n_end,zmat)
    IMPLICIT NONE
    INTEGER, INTENT(IN)            :: id,nk,jspin
    INTEGER, INTENT(OUT),OPTIONAL  :: nv,nmat
    INTEGER, INTENT(OUT),OPTIONAL  :: neig
    REAL,    INTENT(OUT),OPTIONAL  :: eig(:),w_iks(:)
    REAL,    INTENT(OUT),OPTIONAL  :: evac(:),ello(:,:),el(:,:)
    REAL,    INTENT(OUT),OPTIONAL  :: bk(:),wk
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
    IF (PRESENT(el).OR.PRESENT(ello).OR.PRESENT(evac).OR.PRESENT(bk).OR.PRESENT(wk).OR.&
         PRESENT(nv).OR.PRESENT(nmat)) THEN
       !IO of basis-set information
       READ(d%file_io_id_bas,REC=nrec) nmat_s,d%el_s,d%evac_s,d%ello_s,bkpt,wtkpt,nv_s,d%kvec_s,d%kveclo_s
       IF (PRESENT(el)) el=d%el_s
       IF (PRESENT(evac)) evac=d%evac_s
       IF (PRESENT(ello)) ello=d%ello_s
       IF (PRESENT(bk)) bk=bkpt
       IF (PRESENT(wk)) wk=wtkpt
       IF (PRESENT(nv)) nv=nv_s
       IF (PRESENT(nmat)) nmat=nmat_s
    ENDIF

    IF (PRESENT(w_iks)) THEN
       print *, "R:w_iks:",nrec
        read(d%file_io_id_wiks,REC=nrec) w_iks
    ENDIF
      
    
    IF (.NOT.(PRESENT(eig).OR.PRESENT(neig).OR.PRESENT(zmat))) RETURN
    READ(d%file_io_id_vec,REC=nrec) neig_s
    IF (PRESENT(neig)) THEN
       print *,"R:",neig_s
       neig=neig_s
    ENDIF
    IF (.NOT.(PRESENT(eig).OR.PRESENT(zmat))) RETURN
    ALLOCATE(eig_s(neig_s))
    IF (PRESENT(zmat)) THEN
       IF (zmat%l_real) THEN
          INQUIRE(IOLENGTH=n) neig_s,eig_s,REAL(zmat%z_r)
          IF (n>d%recl_vec) THEN
             print *,n,d%recl_vec
             print *,size(eig_s)
             print *,size(zmat%z_r)
             CALL juDFT_error("BUG: Too long record")
          END IF
          READ(d%file_io_id_vec,REC=nrec) neig_s,eig_s,zmat%z_r
       ELSE
          INQUIRE(IOLENGTH=n) neig_s,eig_s,CMPLX(zmat%z_c)
          IF (n>d%recl_vec) THEN
             print *,n,d%recl_vec
             print *,size(eig_s)
             print *,size(zmat%z_c)
             CALL juDFT_error("BUG: Too long record")
          END IF
          READ(d%file_io_id_vec,REC=nrec) neig_s,eig_s,zmat%z_c
       ENDIF
       print *,"R:",nrec,nk,neig_s
    ELSE
       INQUIRE(IOLENGTH=n) neig_s,eig_s
       IF (n>d%recl_vec) CALL juDFT_error("BUG: Too long record")
       READ(d%file_io_id_vec,REC=nrec) neig_s,eig_s
    ENDIF
    IF (PRESENT(eig)) eig(:min(size(eig),neig_s))=eig_s(:min(size(eig),neig_s))
   
  END SUBROUTINE read_eig

  SUBROUTINE write_eig(id,nk,jspin,neig,neig_total,nv,nmat,bk,wk, &
       eig,w_iks,el,ello,evac,nlotot,n_size,n_rank,zmat)
    INTEGER, INTENT(IN)          :: id,nk,jspin
    INTEGER, INTENT(IN),OPTIONAL :: n_size,n_rank
    REAL,    INTENT(IN),OPTIONAL :: wk
    INTEGER, INTENT(IN),OPTIONAL :: neig,nv,nmat,nlotot,neig_total
    REAL,    INTENT(IN),OPTIONAL :: bk(3),eig(:),el(:,:),w_iks(:)
    REAL,    INTENT(IN),OPTIONAL :: evac(:),ello(:,:)
    TYPE(t_zmat),INTENT(IN),OPTIONAL :: zmat

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
    IF (PRESENT(nmat).AND..NOT.PRESENT(el)) THEN
       !IO of basis-set information
       READ(d%file_io_id_bas,REC=nrec,ERR=88) nmat_s,d%el_s,d%evac_s,d%ello_s,bkpt,wtkpt,nv_s,d%kvec_s
88     WRITE(d%file_io_id_bas,REC=nrec)       nmat  ,d%el_s,d%evac_s,d%ello_s,bkpt,wtkpt,nv_s,d%kvec_s
       IF (PRESENT(wk).OR.PRESENT(nv).OR.PRESENT(nlotot) &
            .OR.&
            PRESENT(bk).OR.PRESENT(ello).OR.PRESENT(evac)) THEN
          CALL juDFT_error("BUG:Direct access IO of eig-file only with all scalar data")
       ENDIF
    ELSE IF (PRESENT(el)) THEN
       IF (.NOT.(PRESENT(wk).AND.PRESENT(nv).AND.PRESENT(nmat).AND.PRESENT(nlotot) &
            .AND.&
            PRESENT(bk).AND.PRESENT(el).AND.PRESENT(ello).AND.PRESENT(evac))) THEN
          CALL juDFT_error("BUG:Direct access IO of eig-file only with all data")
       ENDIF
       IF ((SIZE(el).NE.SIZE(d%el_s)).OR.(SIZE(ello).NE.SIZE(d%ello_s).OR.(SIZE(evac).NE.SIZE(d%evac_s)))) THEN
          WRITE(*,*) SHAPE(el),SHAPE(d%el_s)
          WRITE(*,*) SHAPE(ello),SHAPE(d%ello_s)
          WRITE(*,*) SHAPE(evac),SHAPE(d%evac_s)
          CALL juDFT_error("Mismatch of sizes")
       ENDIF
       WRITE(d%file_io_id_bas,REC=nrec) nmat,el,evac,ello,bk,wk,nv,d%kvec_s
    ENDIF
    IF (PRESENT(neig).AND.PRESENT(neig_total)) THEN
       IF (neig.NE.neig_total) THEN
          CALL juDFT_error("Neig and neig_total have to be equal in DA mode",calledby="eig66_da")
       ENDIF
    ENDIF
    IF (PRESENT(w_iks)) THEN
       print *, "W:w_iks:",nrec
       write(d%file_io_id_wiks,REC=nrec) w_iks
    ENDIF
 

    IF (.NOT.PRESENT(eig).OR..NOT.PRESENT(neig)) RETURN
    !Now the IO of the eigenvalues/vectors
    IF (PRESENT(zmat)) THEN
       IF (zmat%l_real) THEN
          INQUIRE(IOLENGTH=r_len) neig,eig,REAL(zmat%z_r)
          IF (r_len>d%recl_vec) CALL juDFT_error("BUG: too long record")
          WRITE(d%file_io_id_vec,REC=nrec) neig,eig,REAL(zmat%z_r)
       ELSE
          INQUIRE(IOLENGTH=r_len) neig,eig(:neig),CMPLX(zmat%z_c)
          IF (r_len>d%recl_vec) CALL juDFT_error("BUG: too long record")
          WRITE(d%file_io_id_vec,REC=nrec) neig,eig(:neig),CMPLX(zmat%z_c)
       ENDIF
       print *,"W:",nrec,nk,neig
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
