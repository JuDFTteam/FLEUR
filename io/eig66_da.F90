MODULE m_eig66_da
#include "juDFT_env.h"
    ! Do the IO of the eig-file in fortran direct-access
    ! The eig-file is split into two parts:
    ! eig.bas contains the basis-set information
    ! eig.vec contains the eigenvalues and the eigenvectors
    ! The record number is given by nrec=nk+(jspin-1)*nkpts
    ! each record contains:
    ! eig.bas: el,evac,ello,bkpt,wtkpt,nv,nmat,k1,k2,k3,kveclo
    ! eig.vec: ne,eig,z**
    !**: real or complex depending on calculation type
    USE m_eig66_data
    IMPLICIT NONE

CONTAINS
    subroutine priv_find_data(id,d)
        INTEGER,INTENT(IN)            :: id
        TYPE(t_data_DA),POINTER,INTENT(out)   :: d

        class(t_data),pointer   ::dp
        call eig66_find_data(dp,id)
        select type(dp)
            type is (t_data_da)
            d=>dp
            class default
            call judft_error("BUG: wrong datatype in eig66_da")
        END SELECT
    END subroutine

    SUBROUTINE open_eig(id,nmat,neig,nkpts,jspins,lmax,nlo,ntype,nlotot,create,filename)
        INTEGER, INTENT(IN) :: id,nmat,neig,nkpts,jspins,nlo,ntype,lmax,nlotot
        LOGICAL, INTENT(IN) :: create
        CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: filename
        !locals
        LOGICAL :: l_file
        INTEGER :: i1,recl_z,recl_eig
        REAL    :: r1,r3(3)
        COMPLEX :: c1
        TYPE(t_data_DA),POINTER:: d

        call priv_find_data(id,d)

        IF (present(filename)) d%fname=filename
        call eig66_data_storedefault(d,jspins,nkpts,nmat,neig,lmax,nlotot,nlo,ntype)

        !Allocate the storage for the DATA always read/write
        ALLOCATE(d%el_s(0:lmax,ntype),d%ello_s(nlo,ntype),d%evac_s(2))
        ALLOCATE(d%kvec_s(nmat,3),d%kveclo_s(nlotot))
        !Calculate the record length
        INQUIRE(IOLENGTH=recl_eig) d%el_s,d%evac_s,d%ello_s,r3,r1,i1,i1,d%kvec_s,d%kveclo_s
        d%recl_bas=recl_eig
        INQUIRE(IOLENGTH=recl_eig) r1
        recl_eig=recl_eig*(neig+2) ! add a 2 for integer 'neig'
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
            INQUIRE(IOLENGTH=recl_z) r1
#else
        INQUIRE(IOLENGTH=recl_z) c1
#endif
        recl_z=recl_z*nmat*neig

        d%recl_vec=recl_eig+recl_z


        IF (create) THEN
            d%file_io_id_bas=priv_free_uid()
            INQUIRE(file=trim(d%fname)//".bas",opened=l_file)
            DO while(l_file)
               print *,"eig66_open_da:",d%fname," in use"
               d%fname=trim(d%fname)//"6"
               INQUIRE(file=trim(d%fname)//".bas",opened=l_file)
            ENDDO
            OPEN(d%file_io_id_bas,FILE=trim(d%fname)//".bas",ACCESS='direct',FORM='unformatted',RECL=d%recl_bas,STATUS='unknown')
            d%file_io_id_vec=priv_free_uid()
            OPEN(d%file_io_id_vec,FILE=trim(d%fname)//".vec",ACCESS='direct',FORM='unformatted',RECL=d%recl_vec,STATUS='unknown')
        ELSE
            d%file_io_id_bas=priv_free_uid()
            OPEN(d%file_io_id_bas,FILE=trim(d%fname)//".bas",ACCESS='direct',FORM='unformatted',RECL=d%recl_bas,STATUS='old')
            d%file_io_id_vec=priv_free_uid()
            OPEN(d%file_io_id_vec,FILE=trim(d%fname)//".vec",ACCESS='direct',FORM='unformatted',RECL=d%recl_vec,STATUS='old')
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
        END FUNCTION
    END SUBROUTINE open_eig
    SUBROUTINE close_eig(id,filename)
        INTEGER,INTENT(IN)::id
        CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: filename
        type(t_data_DA),pointer:: d

        call priv_find_data(id,d)

        DEALLOCATE(d%el_s,d%ello_s,d%evac_s,d%kvec_s,d%kveclo_s)
        CLOSE(d%file_io_id_bas)
        CLOSE(d%file_io_id_vec)
        d%recl_vec=0
        d%recl_bas=0

        !If a filename was given and the name is not the current filename then rename
        IF (present(filename)) THEN
            IF (filename.NE.d%fname) THEN
                CALL system("mv "//trim(d%fname)//".bas "//trim(filename)//".bas")
                CALL system("mv "//trim(d%fname)//".vec "//trim(filename)//".vec")
            ENDIF
        ENDIF
        d%fname="eig"
        CALL eig66_remove_data(id)
    END SUBROUTINE close_eig
    SUBROUTINE read_eig(id,nk,jspin,nv,nmat,k1,k2,k3,bk,wk,neig,eig,el,ello,evac,kveclo,n_start,n_end,z)
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

        !Local variables
        INTEGER:: nv_s,nmat_s,n,nrec,neig_s
        REAL   :: bkpt(3),wtkpt
        REAL,ALLOCATABLE::eig_s(:),zr_s(:,:)
        COMPLEX,ALLOCATABLE::zc_s(:,:)
        REAL,POINTER :: zr(:,:)
        COMPLEX,POINTER:: zc(:,:)
        type(t_data_DA),POINTER:: d



        call priv_find_data(id,d)
        ! check if io is performed correctly
        IF (present(n_start)) THEN
            IF (n_start/=1) &
                CALL juDFT_error("In direct access mode only all eigenstates can be read")
        ENDIF

        nrec=nk+(jspin-1)*d%nkpts
        IF (present(el).OR.present(ello).OR.present(evac).OR.present(bk).OR.present(wk).OR.&
            present(nv).OR.present(nmat).OR.present(k1).OR.present(k2).OR.present(k3).OR.&
            present(kveclo)) THEN
            !IO of basis-set information
            READ(d%file_io_id_bas,REC=nrec) nmat_s,d%el_s,d%evac_s,d%ello_s,bkpt,wtkpt,nv_s,d%kvec_s,d%kveclo_s
            IF (present(el)) el=d%el_s
            IF (present(evac)) evac=d%evac_s
            IF (present(ello)) ello=d%ello_s
            IF (present(bk)) bk=bkpt
            IF (present(wk)) wk=wtkpt
            IF (present(nv)) nv=nv_s
            IF (present(nmat)) nmat=nmat_s
            IF (present(k1)) k1=d%kvec_s(:,1)
            IF (present(k2)) k2=d%kvec_s(:,2)
            IF (present(k3)) k3=d%kvec_s(:,3)
            IF (present(kveclo)) kveclo=d%kveclo_s
        ENDIF

        IF (.NOT.(present(eig).OR.present(neig).OR.present(z))) RETURN

        IF (.NOT.(present(eig).OR.present(z))) THEN
            READ(d%file_io_id_vec,REC=nrec) neig
            RETURN
        ENDIF
        IF (present(eig)) THEN
            ALLOCATE(eig_s(size(eig)))
        ENDIF
        IF (present(z).and..not.present(eig)) THEN
            ALLOCATE(eig_s(size(z,2)))
        ENDIF
        IF (present(z)) THEN
            SELECT TYPE(z)
                TYPE IS (real)
                INQUIRE(IOLENGTH=n) neig_s,eig_s,real(z)
                IF (n>d%recl_vec) CALL juDFT_error("BUG: Too long record")
                zr=>z
                READ(d%file_io_id_vec,REC=nrec) neig_s,eig_s,zr
                TYPE IS (complex)
                INQUIRE(IOLENGTH=n) neig_s,eig_s,cmplx(z)
                IF (n>d%recl_vec) CALL juDFT_error("BUG: Too long record")
                zc=>z
                READ(d%file_io_id_vec,REC=nrec) neig_s,eig_s,zc
            END SELECT
        ELSE
            INQUIRE(IOLENGTH=n) neig_s,eig_s
            IF (n>d%recl_vec) CALL juDFT_error("BUG: Too long record")
            READ(d%file_io_id_vec,REC=nrec) neig_s,eig_s
        ENDIF
        IF (present(neig)) neig=neig_s
        IF (present(eig)) eig=eig_s

    END SUBROUTINE read_eig

    SUBROUTINE write_eig(id,nk,jspin,neig,neig_total,nv,nmat,k1,k2,k3,bk,wk, &
        eig,el,ello,evac,nlotot,kveclo,n_size,n_rank,z)
        INTEGER, INTENT(IN)          :: id,nk,jspin
        INTEGER, INTENT(IN),OPTIONAL :: n_size,n_rank
        REAL,    INTENT(IN),OPTIONAL :: wk
        INTEGER, INTENT(IN),OPTIONAL :: neig,nv,nmat,nlotot,neig_total
        INTEGER, INTENT(IN),OPTIONAL :: k1(:),k2(:),k3(:),kveclo(:)
        REAL,    INTENT(IN),OPTIONAL :: bk(3),eig(:),el(:,:)
        REAL,    INTENT(IN),OPTIONAL :: evac(:),ello(:,:)
        CLASS(*),INTENT(IN),OPTIONAL :: z(:,:)

        INTEGER:: nrec,r_len
        INTEGER:: nv_s,nmat_s
        REAL   :: bkpt(3),wtkpt
        type(t_data_DA),POINTER:: d

        call priv_find_data(id,d)
        !This mode requires all data to be written at once!!

        IF (present(n_size).AND.present(n_rank)) THEN
            IF (n_size/=1.OR.n_rank/=0) &
                CALL juDFT_error("Direct Access IO not possible in eigenvalue parallel code")
        ENDIF
        !check record length
        !INQUIRE(iolength=r_len) nmat,el,evac,ello,bk,wk,nv,d%kvec_s,kveclo
        !if (r_len>recl_bas) call juDFT_error("BUG: too long record")

        !Now it is time for the IO :-)
        nrec=nk+(jspin-1)*d%nkpts
        IF (present(nmat).AND..NOT.present(el)) THEN
            !IO of basis-set information
            READ(d%file_io_id_bas,REC=nrec,ERR=88) nmat_s,d%el_s,d%evac_s,d%ello_s,bkpt,wtkpt,nv_s,d%kvec_s,d%kveclo_s
88          WRITE(d%file_io_id_bas,REC=nrec)       nmat  ,d%el_s,d%evac_s,d%ello_s,bkpt,wtkpt,nv_s,d%kvec_s,d%kveclo_s
            IF (present(wk).OR.present(nv).OR.present(nlotot) &
                .OR.present(k1).OR.present(k2).OR.present(k3).OR.present(kveclo).OR.&
                present(bk).OR.present(ello).OR.present(evac)) THEN
                CALL juDFT_error("BUG:Direct access IO of eig-file only with all scalar data")
            ENDIF
        ELSE
            IF (.NOT.(present(wk).AND.present(nv).AND.present(nmat).AND.present(nlotot) &
                .AND.present(k1).AND.present(k2).AND.present(k3).AND.present(kveclo).AND.&
                present(bk).AND.present(el).AND.present(ello).AND.present(evac))) THEN
                CALL juDFT_error("BUG:Direct access IO of eig-file only with all data")
            ENDIF
            d%kvec_s(:,1)=k1
            d%kvec_s(:,2)=k2
            d%kvec_s(:,3)=k3
            if ((size(el).ne.size(d%el_s)).or.(size(ello).ne.size(d%ello_s).or.(size(evac).ne.size(d%evac_s)))) THEN
                write(*,*) shape(el),shape(d%el_s)
                write(*,*) shape(ello),shape(d%ello_s)
                write(*,*) shape(evac),shape(d%evac_s)
                call juDFT_error("Mismatch of sizes")
            endif
            WRITE(d%file_io_id_bas,REC=nrec) nmat,el,evac,ello,bk,wk,nv,d%kvec_s,kveclo

        ENDIF
        IF (present(neig).AND.present(neig_total)) THEN
            IF (neig.NE.neig_total) THEN
                CALL juDFT_error("Neig and neig_total have to be equal in DA mode",calledby="eig66_da")
            ENDIF
        ENDIF
        IF (.NOT.present(eig).OR..NOT.present(neig)) RETURN
        !Now the IO of the eigenvalues/vectors
        IF (present(z)) THEN
            SELECT TYPE(z)
                TYPE IS (real)
                INQUIRE(IOLENGTH=r_len) neig,eig,real(z)
                IF (r_len>d%recl_vec) CALL juDFT_error("BUG: too long record")
                WRITE(d%file_io_id_vec,REC=nrec) neig,eig,real(z)
                TYPE IS (complex)
                INQUIRE(IOLENGTH=r_len) neig,eig,cmplx(z)
                IF (r_len>d%recl_vec) CALL juDFT_error("BUG: too long record")
                WRITE(d%file_io_id_vec,REC=nrec) neig,eig,cmplx(z)
            END SELECT
        ELSE
            INQUIRE(IOLENGTH=r_len) neig,eig
            IF (r_len>d%recl_vec) CALL juDFT_error("BUG: too long record")
            WRITE(d%file_io_id_vec,REC=nrec) neig,eig
        ENDIF

    END SUBROUTINE write_eig
END MODULE m_eig66_da
