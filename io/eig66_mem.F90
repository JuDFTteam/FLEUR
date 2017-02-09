MODULE m_eig66_mem
#include "juDFT_env.h"
  ! Do the IO of the eig-file into memory
  ! The eig-file is split into four arrays:
  ! eig_int contains the basis-set information/integers (nv,nmat,ne,k1,k2,k3,kveclo)
  ! eig_real contains the basis-set information/real (el,evac,ello,bkpt,wtkpt)
  ! eig_eig contains the eigenvalues
  ! eig_vec contains the eigenvectors
  ! The record number is given by nrec=nk+(jspin-1)*nkpts
  USE m_eig66_data
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
    length=3 !nv+nmat+ne
    length=length+nmat*3 !k1,k2,k3
    length=length+nlotot !kveclo

    ALLOCATE(d%eig_int(length,jspins*nkpts))

    !d%eig_real
    length=3+1+2 !bk,wk,evac
    length=length+(lmax+1)*ntype !el
    length=length+nlo*ntype  !ello
    ALLOCATE(d%eig_real(length,jspins*nkpts))
    !d%eig_eig
    length=jspins
    IF (l_noco) length=1
    ALLOCATE(d%eig_eig(neig,jspins*nkpts))
    !d%eig_vec
    if (l_real.and..not.l_soc) THEN
       print *, "Allocate real in eig66_mem"
       ALLOCATE(d%eig_vecr(nmat*neig,length*nkpts))
    else
       print *, "Allocate complex in eig66_mem"
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
      INTEGER :: k1(nmat),k2(nmat),k3(nmat),kveclo(nlotot)
      REAL    :: eig(neig),ello(d%nlo,d%ntype),el(d%lmax,d%ntype)
      REAL    :: z_r(nmat,neig)
      COMPLEX :: z_c(nmat,neig)
      tmp_id=eig66_data_newid(DA_mode)
      IF (d%l_dos) CPP_error("Can not read DOS-data")
      CALL open_eig_IO(tmp_id,nmat,neig,nkpts,jspins,d%lmax,d%nlo,d%ntype,nlotot,.FALSE.,.FALSE.,l_real,l_soc,.FALSE.,.FALSE.,filename)
      DO jspin=1,jspins
         DO nk=1,nkpts
            if (l_real) THEN
               CALL read_eig_IO(tmp_id,nk,jspin,nv,i,k1,k2,k3,bk3,wk,ii,eig,el,ello,evac,kveclo,z=z_r)
               CALL write_eig(id,nk,jspin,ii,ii,nv,i,k1,k2,k3,bk3,wk,eig,el,ello,evac,nlotot,kveclo,z=z_r)
            else
               CALL read_eig_IO(tmp_id,nk,jspin,nv,i,k1,k2,k3,bk3,wk,ii,eig,el,ello,evac,kveclo,z=z_c)
               CALL write_eig(id,nk,jspin,ii,ii,nv,i,k1,k2,k3,bk3,wk,eig,el,ello,evac,nlotot,kveclo,z=z_c)
            end if
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
          IF (ALLOCATED(d%eig_real)) DEALLOCATE(d%eig_real)
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
      INTEGER :: k1(d%nmat),k2(d%nmat),k3(d%nmat),kveclo(SIZE(d%eig_int,1)-3-3*d%nmat)
      REAL    :: eig(SIZE(d%eig_eig,1)),ello(d%nlo,d%ntype),el(d%lmax,d%ntype)
      REAL    :: z_r(d%nmat,SIZE(d%eig_eig,1))
      COMPLEX :: z_c(d%nmat,SIZE(d%eig_eig,1))
      tmp_id=eig66_data_newid(DA_mode)
      IF (d%l_dos) CPP_error("Could not write DOS data")
      CALL open_eig_DA(tmp_id,d%nmat,d%neig,d%nkpts,d%jspins,d%lmax,d%nlo,d%ntype,d%nlotot,.FALSE.,.FALSE.,d%l_real,d%l_soc,.FALSE.,.FALSE.,filename)
      DO jspin=1,d%jspins
         DO nk=1,d%nkpts
            IF (d%l_real) THEN
               CALL read_eig(id,nk,jspin,nv,i,k1,k2,k3,bk3,wk,ii,eig,el,ello,evac,kveclo,z=z_r)
               CALL write_eig_DA(tmp_id,nk,jspin,ii,ii,nv,i,k1,k2,k3,bk3,wk,eig,el,ello,evac,nlotot,kveclo,z=z_r)
            else
               CALL read_eig(id,nk,jspin,nv,i,k1,k2,k3,bk3,wk,ii,eig,el,ello,evac,kveclo,z=z_c)
               CALL write_eig_DA(tmp_id,nk,jspin,ii,ii,nv,i,k1,k2,k3,bk3,wk,eig,el,ello,evac,nlotot,kveclo,z=z_c)
            end IF
         ENDDO
      ENDDO
      CALL close_eig_DA(tmp_id)
      CALL eig66_remove_data(id)
    END SUBROUTINE priv_writetofile
  END SUBROUTINE close_eig

  SUBROUTINE write_dos(id,nk,jspin,qal,qvac,qis,qvlay,qstars,ksym,jsym,mcd,qintsl,qmtsl,qmtp,orbcomp)
    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: id,nk,jspin
    REAL,INTENT(IN)              :: qal(:,:,:),qvac(:,:),qis(:),qvlay(:,:,:)
    COMPLEX,INTENT(IN)           :: qstars(:,:,:,:)
    INTEGER,INTENT(IN)           :: ksym(:),jsym(:)
    REAL,INTENT(IN),OPTIONAL     :: mcd(:,:,:)
    REAL,INTENT(IN),OPTIONAL     :: qintsl(:,:),qmtsl(:,:),qmtp(:,:),orbcomp(:,:,:)

    INTEGER::nrec
    TYPE(t_data_mem),POINTER:: d
    CALL priv_find_data(id,d)

    nrec=nk+(jspin-1)*d%nkpts

    d%qal(:,:,:,nrec)=qal
    d%qvac(:,:,nrec)=qvac
    d%qis(:,nrec)=qis
    d%qvlay(:,:,:,nrec)=qvlay
    d%qstars(:,:,:,:,nrec)=qstars
    d%ksym(:,nrec)=ksym
    d%jsym(:,nrec)=jsym
    IF (d%l_mcd.AND.PRESENT(mcd)) d%mcd(:,:,:,nrec)=mcd
    IF (d%l_orb.AND.PRESENT(qintsl)) THEN
       d%qintsl(:,:,nrec)=qintsl
       d%qmtsl(:,:,nrec)=qmtsl
       d%qmtp(:,:,nrec)=qmtp
       d%orbcomp(:,:,:,nrec)=orbcomp
    ENDIF
  END SUBROUTINE write_dos

  SUBROUTINE read_dos(id,nk,jspin,qal,qvac,qis,qvlay,qstars,ksym,jsym,mcd,qintsl,qmtsl,qmtp,orbcomp)
    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: id,nk,jspin
    REAL,INTENT(OUT)              :: qal(:,:,:),qvac(:,:),qis(:),qvlay(:,:,:)
    COMPLEX,INTENT(OUT)           :: qstars(:,:,:,:)
    INTEGER,INTENT(OUT)           :: ksym(:),jsym(:)
    REAL,INTENT(OUT),OPTIONAL     :: mcd(:,:,:)
    REAL,INTENT(OUT),OPTIONAL     :: qintsl(:,:),qmtsl(:,:),qmtp(:,:),orbcomp(:,:,:)

    INTEGER::nrec
    TYPE(t_data_mem),POINTER:: d
    CALL priv_find_data(id,d)

    nrec=nk+(jspin-1)*d%nkpts

    qal=d%qal(:,:,:,nrec)
    qvac=d%qvac(:,:,nrec)
    qis=d%qis(:,nrec)
    qvlay=d%qvlay(:,:,:,nrec)
    qstars=d%qstars(:,:,:,:,nrec)
    ksym=d%ksym(:,nrec)
    jsym=d%jsym(:,nrec)
    IF (d%l_mcd.AND.PRESENT(mcd)) mcd=d%mcd(:,:,:,nrec)
    IF (d%l_orb.AND.PRESENT(qintsl)) THEN
       qintsl=d%qintsl(:,:,nrec)
       qmtsl=d%qmtsl(:,:,nrec)
       qmtp=d%qmtp(:,:,nrec)
       orbcomp=d%orbcomp(:,:,:,nrec)
    ENDIF
  END SUBROUTINE read_dos


  SUBROUTINE read_eig(id,nk,jspin,nv,nmat,k1,k2,k3,bk,wk,neig,eig,el,&
       ello,evac,kveclo,n_start,n_end,z)
    IMPLICIT NONE
    INTEGER, INTENT(IN)            :: id,nk,jspin
    INTEGER, INTENT(OUT),OPTIONAL  :: nv,nmat
    INTEGER, INTENT(OUT),OPTIONAL  :: neig
    REAL,    INTENT(OUT),OPTIONAL  :: eig(:)
    INTEGER, INTENT(OUT),OPTIONAL  :: k1(:),k2(:),k3(:),kveclo(:)
    REAL,    INTENT(OUT),OPTIONAL  :: evac(:),ello(:,:),el(:,:)
    REAL,    INTENT(OUT),OPTIONAL  :: bk(:),wk
    INTEGER, INTENT(IN),OPTIONAL   :: n_start,n_end
    CLASS(*),OPTIONAL  :: z(:,:)

    INTEGER::nrec
    TYPE(t_data_mem),POINTER:: d
    CALL priv_find_data(id,d)

    nrec=nk+(jspin-1)*d%nkpts
    ! data from d%eig_int
    IF (PRESENT(nv)) nv=d%eig_int(1,nrec)
    IF (PRESENT(nmat)) nmat=d%eig_int(2,nrec)
    IF (PRESENT(neig)) THEN
       neig=d%eig_int(3,nrec)
    ENDIF
    IF (PRESENT(k1)) THEN
       IF (.NOT.PRESENT(k2).OR..NOT.PRESENT(k3)) CALL juDFT_error("BUG: always read k1,k2,k3")
       k1=d%eig_int(3+1:3+size(k1),nrec)
       k2=d%eig_int(3+d%nmat+1:3+d%nmat+size(k1),nrec)
       k3=d%eig_int(3+2*d%nmat+1:3+2*d%nmat*size(k1),nrec)
    ENDIF
    IF (PRESENT(kveclo)) kveclo=d%eig_int(4+3*d%nmat:3+3*d%nmat+SIZE(kveclo),nrec)

    !data from d%eig_real
    IF (PRESENT(bk)) bk=d%eig_real(1:3,nrec)
    IF (PRESENT(wk)) wk=d%eig_real(4,nrec)
    IF (PRESENT(evac)) evac=d%eig_real(5:6,nrec)
    IF (PRESENT(el)) el=RESHAPE(d%eig_real(7:7+SIZE(el)-1,nrec),SHAPE(el))
    IF (PRESENT(ello)) ello=RESHAPE(d%eig_real(SIZE(d%eig_real,1)-SIZE(ello)+1:,nrec),SHAPE(ello))

    !data from d%eig_eig
    IF (PRESENT(eig)) THEN
       eig=0.0
       eig=d%eig_eig(:SIZE(eig),nrec)
       !print *,"R-eig:",nrec,shape(eig)
       !print*,"R-eig(data):",shape(d%eig_eig)
       !print*,"R:",eig
    ENDIF
    !data from d%eig_vec

    IF (PRESENT(z)) THEN
       SELECT TYPE(z)
       TYPE is (REAL)
          IF (.NOT.ALLOCATED(d%eig_vecr)) THEN
             IF (.NOT.ALLOCATED(d%eig_vecc)) CALL juDFT_error("BUG: can not read complex vectors from memory")
             z=REAL(RESHAPE(d%eig_vecc(:SIZE(z),nrec),SHAPE(z)))
          ELSE
             z=RESHAPE(d%eig_vecr(:SIZE(z),nrec),SHAPE(z))
          ENDIF
       TYPE is (COMPLEX)
          IF (.NOT.ALLOCATED(d%eig_vecc)) CALL juDFT_error("BUG: can not read complex vectors from memory")
          z=RESHAPE(d%eig_vecc(:SIZE(z),nrec),SHAPE(z))
       END SELECT
    ENDIF
  END SUBROUTINE read_eig


  SUBROUTINE write_eig(id,nk,jspin,neig,neig_total,nv,nmat,k1,k2,k3,bk,wk, &
       eig,el,ello,evac,                     &
       nlotot,kveclo,n_size,n_rank,z)
    INTEGER, INTENT(IN)          :: id,nk,jspin
    INTEGER, INTENT(IN),OPTIONAL :: n_size,n_rank
    REAL,    INTENT(IN),OPTIONAL :: wk
    INTEGER, INTENT(IN),OPTIONAL :: neig,neig_total,nv,nmat,nlotot
    INTEGER, INTENT(IN),OPTIONAL :: k1(:),k2(:),k3(:),kveclo(:)
    REAL,    INTENT(IN),OPTIONAL :: bk(3),eig(:),el(:,:)
    REAL,    INTENT(IN),OPTIONAL :: evac(:),ello(:,:)
    CLASS(*),INTENT(IN),OPTIONAL :: z(:,:)
    INTEGER::nrec
    TYPE(t_data_mem),POINTER:: d
    CALL priv_find_data(id,d)

    nrec=nk+(jspin-1)*d%nkpts
    ! data from d%eig_int
    IF (PRESENT(nv)) d%eig_int(1,nrec)=nv
    IF (PRESENT(nmat)) d%eig_int(2,nrec)=nmat
    IF (PRESENT(neig)) THEN
       IF (PRESENT(neig_total)) THEN
          IF (neig.NE.neig_total) STOP "BUG in eig_mem"
          d%eig_int(3,nrec)=neig_total
       ELSE
          STOP "BUG2 in eig_mem"
       ENDIF
    ENDIF

    IF (PRESENT(k1)) THEN
       IF (.NOT.PRESENT(k2).OR..NOT.PRESENT(k3)) CALL juDFT_error("BUG: always write k1,k2,k3")
       d%eig_int(3+1:3+size(k1),nrec)=k1
       d%eig_int(3+d%nmat+1:3+d%nmat+size(k1),nrec)=k2
       d%eig_int(3+2*d%nmat+1:3+2*d%nmat+size(k1),nrec)=k3
    ENDIF
    IF (PRESENT(kveclo)) d%eig_int(4+3*d%nmat:3+3*d%nmat+SIZE(kveclo),nrec)=kveclo

    !data from d%eig_real
    IF (PRESENT(bk)) d%eig_real(1:3,nrec)=bk
    IF (PRESENT(wk)) d%eig_real(4,nrec)=wk
    IF (PRESENT(evac)) d%eig_real(5:6,nrec)=evac
    IF (PRESENT(el)) d%eig_real(7:7+SIZE(el)-1,nrec)=RESHAPE(el,(/SIZE(el)/))
    IF (PRESENT(ello)) d%eig_real(SIZE(d%eig_real,1)-SIZE(ello)+1:,nrec)=RESHAPE(ello,(/SIZE(ello)/))
    !data from d%eig_eig
    IF (PRESENT(eig)) THEN
       d%eig_eig(:SIZE(eig),nrec)=eig
       !print*,"W-eig:",nrec,shape(eig)
       !print*,"W:",eig
    ENDIF
    !data from d%eig_vec
    IF (PRESENT(z)) THEN
       SELECT TYPE(z)
       TYPE IS (REAL)
          IF (.NOT.ALLOCATED(d%eig_vecr)) THEN
             IF (.NOT.ALLOCATED(d%eig_vecc)) CALL juDFT_error("BUG: can not write complex vectors to memory")
             d%eig_vecc(:SIZE(z),nrec)=RESHAPE(CMPLX(z),(/SIZE(z)/)) !Type cast here
          ELSE
             d%eig_vecr(:SIZE(z),nrec)=RESHAPE(REAL(z),(/SIZE(z)/))
          ENDIF
       TYPE IS(COMPLEX)
          IF (.NOT.ALLOCATED(d%eig_vecc)) CALL juDFT_error("BUG: can not write complex vectors to memory")
          d%eig_vecc(:SIZE(z),nrec)=RESHAPE(CMPLX(z),(/SIZE(z)/))
       END SELECT
    ENDIF


  END SUBROUTINE write_eig


END MODULE m_eig66_mem
