MODULE m_eig66_mem
#include "juDFT_env.h"
  ! Do the IO of the eig-file into memory
  ! The eig-file is split into four arrays:
  ! eig_int contains the basis-set information/integers (ne)
  ! eig_eig contains the eigenvalues
  ! eig_vec contains the eigenvectors
  ! The record number is given by nrec=nk+(jspin-1)*nkpts
  USE m_eig66_data
  USE m_types
  IMPLICIT NONE
CONTAINS

  SUBROUTINE priv_find_data(id,d)
    INTEGER,INTENT(IN)::id
    TYPE(t_data_mem),POINTER,INTENT(out):: d

    CLASS(t_data),POINTER   ::dp
    CALL eig66_find_data(dp,id)
    SELECT TYPE(dp)
    TYPE is (t_data_mem)
       d=>dp
       CLASS default
       CALL judft_error("BUG: wrong datatype in eig66_mem")
    END SELECT
  END SUBROUTINE priv_find_data

  SUBROUTINE open_eig(id,nmat,neig,nkpts,jspins,lmax,nlo,ntype,l_create,l_real,l_soc,nlotot,l_noco,l_dos,l_mcd,l_orb,filename,layers,nstars,ncored,nsld,nat)
    INTEGER, INTENT(IN) :: id,nmat,neig,nkpts,jspins,nlo,ntype,lmax,nlotot
    LOGICAL, INTENT(IN) :: l_noco,l_create,l_real,l_soc
    LOGICAL,INTENT(IN),OPTIONAL::l_dos,l_mcd,l_orb
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: filename
    INTEGER,INTENT(IN),OPTIONAL :: layers,nstars,ncored,nsld,nat
    !locals
    INTEGER:: length
    TYPE(t_data_mem),POINTER:: d
    CALL priv_find_data(id,d)

    IF (ALLOCATED(d%eig_int)) THEN
       IF (.NOT.l_create) THEN
          IF (PRESENT(filename)) CALL priv_readfromfile()
          RETURN
       ENDIF
       CALL close_eig(id,.TRUE.)

    ENDIF

    CALL eig66_data_storedefault(d,jspins,nkpts,nmat,neig,lmax,nlotot,nlo,ntype,l_real,l_soc,l_dos,l_mcd,l_orb)

    !d%eig_int
    ALLOCATE(d%eig_int(jspins*nkpts))

    !d%eig_eig
    length=jspins
    IF (l_noco) length=1
    ALLOCATE(d%eig_eig(neig,2,jspins*nkpts)) !additional dimension for w_iks
    !d%eig_vec
    if (l_real.and..not.l_soc) THEN
       ALLOCATE(d%eig_vecr(nmat*neig,length*nkpts))
    else
       ALLOCATE(d%eig_vecc(nmat*neig,length*nkpts))
    endif
    length=length*nkpts
    IF (d%l_dos) THEN
       ALLOCATE(d%qal(0:3,ntype,neig,length))
       ALLOCATE(d%qvac(neig,2,length))
       ALLOCATE(d%qis(neig,length))
       ALLOCATE(d%qvlay(neig,max(layers,1),2,length))
       ALLOCATE(d%qstars(nstars,neig,max(layers,1),2,length))
       ALLOCATE(d%ksym(neig,length))
       ALLOCATE(d%jsym(neig,length))
       IF (l_mcd) ALLOCATE(d%mcd(3*ntype,ncored,neig,length))
       IF (l_orb) THEN
          ALLOCATE(d%qintsl(nsld,neig,length))
          ALLOCATE(d%qmtsl(nsld,neig,length))
          ALLOCATE(d%qmtp(neig,nat,length))
          ALLOCATE(d%orbcomp(neig,23,nat,length))
       ENDIF
    ENDIF
    IF (PRESENT(filename)) CALL priv_readfromfile()
  CONTAINS
    SUBROUTINE priv_readfromfile()
      USE m_eig66_da,ONLY:open_eig_IO=>open_eig,read_eig_IO=>read_eig,close_eig_IO=>close_eig
      INTEGER:: jspin,nk,i,ii,iii,nv,tmp_id
      REAL   :: wk,bk3(3),evac(2)
      REAL    :: eig(neig),w_iks(neig),ello(d%nlo,d%ntype),el(d%lmax,d%ntype)
      TYPE(t_zmat):: zmat

      zmat%l_real=l_real
      zmat%nbasfcn=nmat
      zmat%nbands=neig
      ALLOCATE(zmat%z_r(nmat,neig),zmat%z_c(nmat,neig))
    
      tmp_id=eig66_data_newid(DA_mode)
      IF (d%l_dos) CPP_error("Can not read DOS-data")
      CALL open_eig_IO(tmp_id,nmat,neig,nkpts,jspins,d%lmax,d%nlo,d%ntype,nlotot,.FALSE.,.FALSE.,l_real,l_soc,.FALSE.,.FALSE.,filename)
      DO jspin=1,jspins
         DO nk=1,nkpts
            CALL read_eig_IO(tmp_id,nk,jspin,i,eig,w_iks,zmat=zmat)
            !CALL write_eig(id,nk,jspin,i,i,eig,w_iks,zmat=zmat)
         ENDDO
      ENDDO
      CALL close_eig_IO(tmp_id)
    END SUBROUTINE priv_readfromfile

  END SUBROUTINE open_eig

  SUBROUTINE close_eig(id,delete,filename)
    INTEGER,INTENT(in)         :: id
    LOGICAL,INTENT(in),OPTIONAL::delete
    CHARACTER(len=*),OPTIONAL,INTENT(in)::filename
    TYPE(t_data_mem),POINTER:: d
    CALL priv_find_data(id,d)

    IF (PRESENT(filename)) CALL priv_writetofile()

    IF (PRESENT(delete)) THEN
       IF (delete) THEN
          IF (ALLOCATED(d%eig_int)) DEALLOCATE(d%eig_int)
          IF (ALLOCATED(d%eig_eig)) DEALLOCATE(d%eig_eig)
          IF (ALLOCATED(d%eig_vecr)) DEALLOCATE(d%eig_vecr)
          IF (ALLOCATED(d%eig_vecc)) DEALLOCATE(d%eig_vecc)
       ENDIF
    ENDIF
  CONTAINS
    SUBROUTINE priv_writetofile()
      USE m_eig66_DA,ONLY:open_eig_DA=>open_eig,write_eig_DA=>write_eig,close_eig_DA=>close_eig
      IMPLICIT NONE

      INTEGER:: nlotot,nk,jspin,nv,i,ii,tmp_id
      REAL   :: wk,bk3(3),evac(2)
      REAL    :: eig(SIZE(d%eig_eig,1)),w_iks(SIZE(d%eig_eig,1)),ello(d%nlo,d%ntype),el(d%lmax,d%ntype)
      TYPE(t_mat)::zmat
      zmat%l_real=d%l_real
      zmat%matsize1=d%nmat
      zmat%matsize2=SIZE(d%eig_eig,1)
      ALLOCATE(zmat%data_r(d%nmat,SIZE(d%eig_eig,1)),zmat%data_c(d%nmat,SIZE(d%eig_eig,1)))
      tmp_id=eig66_data_newid(DA_mode)
      IF (d%l_dos) CPP_error("Could not write DOS data")
      CALL open_eig_DA(tmp_id,d%nmat,d%neig,d%nkpts,d%jspins,d%lmax,d%nlo,d%ntype,d%nlotot,.FALSE.,.FALSE.,d%l_real,d%l_soc,.FALSE.,.FALSE.,filename)
      DO jspin=1,d%jspins
         DO nk=1,d%nkpts
            !TODO this code is no longer working
            STOP "BUG"
               !CALL read_eig(id,nk,jspin,nv,i,bk3,wk,ii,eig,w_iks,el,ello,evac,zmat=zmat)
               !CALL write_eig_DA(tmp_id,nk,jspin,ii,ii,nv,i,bk3,wk,eig,w_iks,el,ello,evac,nlotot,zmat=zmat)
           ENDDO
      ENDDO
      CALL close_eig_DA(tmp_id)
      CALL eig66_remove_data(id)
    END SUBROUTINE priv_writetofile
  END SUBROUTINE close_eig

  SUBROUTINE write_dos(id,nk,jspin,mcd,qintsl,qmtsl,qmtp,orbcomp)
    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: id,nk,jspin
    REAL,INTENT(IN),OPTIONAL     :: mcd(:,:,:)
    REAL,INTENT(IN),OPTIONAL     :: qintsl(:,:),qmtsl(:,:),qmtp(:,:),orbcomp(:,:,:)

    INTEGER::nrec
    TYPE(t_data_mem),POINTER:: d
    CALL priv_find_data(id,d)

    nrec=nk+(jspin-1)*d%nkpts

    IF (d%l_mcd.AND.PRESENT(mcd)) d%mcd(:,:,:,nrec)=mcd
    IF (d%l_orb.AND.PRESENT(qintsl)) THEN
       d%qintsl(:,:,nrec)=qintsl
       d%qmtsl(:,:,nrec)=qmtsl
       d%qmtp(:,:,nrec)=qmtp
       d%orbcomp(:,:,:,nrec)=orbcomp
    ENDIF
  END SUBROUTINE write_dos

  SUBROUTINE read_dos(id,nk,jspin,mcd,qintsl,qmtsl,qmtp,orbcomp)
    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: id,nk,jspin
    REAL,INTENT(OUT),OPTIONAL     :: mcd(:,:,:)
    REAL,INTENT(OUT),OPTIONAL     :: qintsl(:,:),qmtsl(:,:),qmtp(:,:),orbcomp(:,:,:)

    INTEGER::nrec
    TYPE(t_data_mem),POINTER:: d
    CALL priv_find_data(id,d)

    nrec=nk+(jspin-1)*d%nkpts

    IF (d%l_mcd.AND.PRESENT(mcd)) mcd=d%mcd(:,:,:,nrec)
    IF (d%l_orb.AND.PRESENT(qintsl)) THEN
       qintsl=d%qintsl(:,:,nrec)
       qmtsl=d%qmtsl(:,:,nrec)
       qmtp=d%qmtp(:,:,nrec)
       orbcomp=d%orbcomp(:,:,:,nrec)
    ENDIF
  END SUBROUTINE read_dos


  SUBROUTINE read_eig(id,nk,jspin,neig,eig,w_iks,n_start,n_end,zmat)
    IMPLICIT NONE
    INTEGER, INTENT(IN)            :: id,nk,jspin
    INTEGER, INTENT(OUT),OPTIONAL  :: neig
    REAL,    INTENT(OUT),OPTIONAL  :: eig(:),w_iks(:)
    INTEGER, INTENT(IN),OPTIONAL   :: n_start,n_end
    TYPE(t_zMAT),OPTIONAL  :: zmat

    INTEGER::nrec, arrayStart
    TYPE(t_data_mem),POINTER:: d
    CALL priv_find_data(id,d)

    nrec=nk+(jspin-1)*d%nkpts
    ! data from d%eig_int
    IF (PRESENT(neig)) THEN
       neig=d%eig_int(nrec)
    ENDIF
  
    !data from d%eig_eig
    IF (PRESENT(eig)) THEN
       eig=0.0
       eig=d%eig_eig(:SIZE(eig),1,nrec)
    ENDIF
    IF (PRESENT(w_iks)) THEN
       w_iks=0.0
       w_iks=d%eig_eig(:SIZE(w_iks),2,nrec)
    ENDIF
    
    !data from d%eig_vec

    arrayStart = 1
    IF(PRESENT(n_start)) THEN
       arrayStart = (n_start-1)*zMat%nbasfcn+1
    END IF

    IF (PRESENT(zmat)) THEN
      
       IF (zmat%l_real) THEN
          IF (.NOT.ALLOCATED(d%eig_vecr)) THEN
             IF (.NOT.ALLOCATED(d%eig_vecc)) CALL juDFT_error("BUG: can not read real/complex vectors from memory")
             zmat%z_r=REAL(RESHAPE(d%eig_vecc(arrayStart:arrayStart+SIZE(zmat%z_r)-1,nrec),SHAPE(zmat%z_r)))
          ELSE
             zmat%z_r=RESHAPE(d%eig_vecr(arrayStart:arrayStart+SIZE(zmat%z_r)-1,nrec),SHAPE(zmat%z_r))
          ENDIF
       ELSE !TYPE is (COMPLEX)
          IF (.NOT.ALLOCATED(d%eig_vecc)) CALL juDFT_error("BUG: can not read complex vectors from memory", calledby = "eig66_mem")
          zmat%z_c=RESHAPE(d%eig_vecc(arrayStart:arrayStart+SIZE(zmat%z_c)-1,nrec),SHAPE(zmat%z_c))
       END IF
    ENDIF
  END SUBROUTINE read_eig


  SUBROUTINE write_eig(id,nk,jspin,neig,neig_total,eig,w_iks,n_size,n_rank,zmat)
    INTEGER, INTENT(IN)          :: id,nk,jspin
    INTEGER, INTENT(IN),OPTIONAL :: n_size,n_rank
    INTEGER, INTENT(IN),OPTIONAL :: neig,neig_total
    REAL,    INTENT(IN),OPTIONAL :: eig(:),w_iks(:)
    TYPE(t_mat),INTENT(IN),OPTIONAL :: zmat
    INTEGER::nrec
    TYPE(t_data_mem),POINTER:: d
    CALL priv_find_data(id,d)

    nrec=nk+(jspin-1)*d%nkpts
    ! data from d%eig_int
    IF (PRESENT(neig)) THEN
       IF (PRESENT(neig_total)) THEN
          IF (neig.NE.neig_total) STOP "BUG in eig_mem"
          d%eig_int(nrec)=neig_total
       ELSE
          STOP "BUG2 in eig_mem"
       ENDIF
    ENDIF

  
    !data from d%eig_eig
    IF (PRESENT(eig)) THEN
       d%eig_eig(:SIZE(eig),1,nrec)=eig
    ENDIF
    IF (PRESENT(w_iks)) THEN
       d%eig_eig(:SIZE(w_iks),2,nrec)=w_iks
    ENDIF
    !data from d%eig_vec
    IF (PRESENT(zmat)) THEN
       IF (zmat%l_real) THEN
          IF (.NOT.ALLOCATED(d%eig_vecr)) THEN
             IF (.NOT.ALLOCATED(d%eig_vecc)) CALL juDFT_error("BUG: can not write complex vectors to memory")
             d%eig_vecc(:SIZE(zmat%data_r),nrec)=RESHAPE(CMPLX(zmat%data_r),(/SIZE(zmat%data_r)/)) !Type cast here
          ELSE
             d%eig_vecr(:SIZE(zmat%data_r),nrec)=RESHAPE(REAL(zmat%data_r),(/SIZE(zmat%data_r)/))
          ENDIF
       ELSE
          IF (.NOT.ALLOCATED(d%eig_vecc)) CALL juDFT_error("BUG: can not write complex vectors to memory")
          d%eig_vecc(:SIZE(zmat%data_c),nrec)=RESHAPE(zmat%data_c,(/SIZE(zmat%data_c)/))
       END IF
    ENDIF


  END SUBROUTINE write_eig


END MODULE m_eig66_mem
