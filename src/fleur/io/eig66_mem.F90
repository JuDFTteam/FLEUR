MODULE m_eig66_mem
use m_juDFT
   ! Do the IO of the eig-file into memory
   ! The eig-file is split into four arrays:
   ! eig_int contains the basis-set information/integers (ne)
   ! eig_eig contains the eigenvalues
   ! eig_vec contains the eigenvectors
   ! The record number is given by nrec=nk+(jspin-1)*nkpts
   USE m_eig66_data
   USE m_types_mat
   USE m_juDFT
   IMPLICIT NONE
CONTAINS

   SUBROUTINE priv_find_data(id, d)
      INTEGER, INTENT(IN)::id
      TYPE(t_data_mem), POINTER, INTENT(out):: d

      CLASS(t_data), POINTER   ::dp
      CALL eig66_find_data(dp, id)
      SELECT TYPE (dp)
      TYPE is (t_data_mem)
         d => dp
      CLASS default
         CALL judft_error("BUG: wrong datatype in eig66_mem")
      END SELECT
   END SUBROUTINE priv_find_data

   SUBROUTINE open_eig(id, nmat, neig, nkpts, jspins, l_create, l_real, l_soc, l_noco, l_olap, filename)
      INTEGER, INTENT(IN) :: id, nmat, neig, nkpts, jspins
      LOGICAL, INTENT(IN) :: l_noco, l_create, l_real, l_soc, l_olap
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: filename
      !locals
      INTEGER :: length, ierr
      INTEGER :: elementsize
      CHARACTER(LEN=80) errorString
      TYPE(t_data_mem), POINTER:: d
      CALL priv_find_data(id, d)

      IF (ALLOCATED(d%eig_int)) THEN
         IF (.NOT. l_create) THEN
            IF (PRESENT(filename)) CALL priv_readfromfile()
            RETURN
         ENDIF
         CALL close_eig(id, .TRUE.)

      ENDIF

      CALL eig66_data_storedefault(d, jspins, nkpts, nmat, neig, l_real, l_soc)

      !d%eig_int
      ALLOCATE (d%eig_int(jspins*nkpts))

      !d%eig_eig
      length = jspins
      IF (l_noco) length = 1
      ALLOCATE (d%eig_eig(neig, jspins*nkpts))
      !d%eig_vec
      if (l_real .and. .not. l_soc) THEN
         ALLOCATE (d%eig_vecr(nmat*neig, length*nkpts), source=0.0, STAT=ierr)
         elementsize = 8
      else
         ALLOCATE (d%eig_vecc(nmat*neig, length*nkpts), source=CMPLX(0.0,0.0), STAT=ierr)
         elementsize = 16
      endif
      IF (ierr.NE.0) THEN
         WRITE(errorString,'(a,i0,a,i0,a,i0,a,i0,a,i0,a)') "Could not allocate eigenvector array of size ", &
                                                            elementsize, " x ", nmat, " x ", neig, " x ", length, " x ", nkpts, " bytes."
         CALL juDFT_error(TRIM(ADJUSTL(errorString)), calledby = 'eig66_mem')
      END IF

      !d%olap 
      if(l_olap) then
         if (l_real .and. .not. l_soc) THEN
            ALLOCATE (d%olap_r(nmat**2, length*nkpts))
         else
            ALLOCATE (d%olap_c(nmat**2, length*nkpts))
         endif
      endif
      length = length*nkpts
      IF (PRESENT(filename)) CALL priv_readfromfile()
   CONTAINS
      SUBROUTINE priv_readfromfile()
         USE m_eig66_da, ONLY: open_eig_IO => open_eig, read_eig_IO => read_eig, close_eig_IO => close_eig
         INTEGER:: jspin, nk, i, ii, iii, nv, tmp_id
         REAL   :: wk, bk3(3), evac(2)
         REAL    :: eig(neig)
         TYPE(t_mat):: zmat

         zmat%l_real = l_real
         zmat%matsize1 = nmat
         zmat%matsize2 = neig
         ALLOCATE (zmat%data_r(nmat, neig), zmat%data_c(nmat, neig))

         tmp_id = eig66_data_newid(DA_mode)
         CALL open_eig_IO(tmp_id, nmat, neig, nkpts, jspins, .FALSE., l_real, l_soc, .false., filename)
         DO jspin = 1, jspins
            DO nk = 1, nkpts
               CALL read_eig_IO(tmp_id, nk, jspin, i, eig, zmat=zmat)
               !CALL write_eig(id,nk,jspin,i,i,eig,zmat=zmat)
            ENDDO
         ENDDO
         CALL close_eig_IO(tmp_id)
      END SUBROUTINE priv_readfromfile

   END SUBROUTINE open_eig

   SUBROUTINE close_eig(id, delete, filename)
      INTEGER, INTENT(in)         :: id
      LOGICAL, INTENT(in), OPTIONAL::delete
      CHARACTER(len=*), OPTIONAL, INTENT(in)::filename
      TYPE(t_data_mem), POINTER:: d
      CALL priv_find_data(id, d)

      IF (PRESENT(filename)) CALL priv_writetofile()

      IF (PRESENT(delete)) THEN
         IF (delete) THEN
            IF (ALLOCATED(d%eig_int))  DEALLOCATE (d%eig_int)
            IF (ALLOCATED(d%eig_eig))  DEALLOCATE (d%eig_eig)
            IF (ALLOCATED(d%eig_vecr)) DEALLOCATE (d%eig_vecr)
            IF (ALLOCATED(d%eig_vecc)) DEALLOCATE (d%eig_vecc)
            if (allocated(d%olap_r))   deallocate (d%olap_r)
            if (allocated(d%olap_c))   deallocate (d%olap_c)
         ENDIF
      ENDIF
   CONTAINS
      SUBROUTINE priv_writetofile()
         USE m_eig66_DA, ONLY: open_eig_DA => open_eig, write_eig_DA => write_eig, close_eig_DA => close_eig
         IMPLICIT NONE

         INTEGER:: nk, jspin, nv, i, ii, tmp_id
         REAL   :: wk, bk3(3), evac(2)
         REAL    :: eig(SIZE(d%eig_eig, 1))
         TYPE(t_mat)::zmat
         zmat%l_real = d%l_real
         zmat%matsize1 = d%nmat
         zmat%matsize2 = SIZE(d%eig_eig, 1)
         ALLOCATE (zmat%data_r(d%nmat, SIZE(d%eig_eig, 1)), zmat%data_c(d%nmat, SIZE(d%eig_eig, 1)))
         tmp_id = eig66_data_newid(DA_mode)
         CALL open_eig_DA(tmp_id, d%nmat, d%neig, d%nkpts, d%jspins, .FALSE., d%l_real, d%l_soc, .false., filename)
         DO jspin = 1, d%jspins
            DO nk = 1, d%nkpts
               !TODO this code is no longer working
               STOP "BUG"
               !CALL read_eig(id,nk,jspin,nv,i,bk3,wk,ii,eig,el,ello,evac,zmat=zmat)
               !CALL write_eig_DA(tmp_id,nk,jspin,ii,ii,nv,i,bk3,wk,eig,el,ello,evac,nlotot,zmat=zmat)
            ENDDO
         ENDDO
         CALL close_eig_DA(tmp_id)
         CALL eig66_remove_data(id)
      END SUBROUTINE priv_writetofile
   END SUBROUTINE close_eig

   SUBROUTINE read_eig(id, nk, jspin, neig, eig, list, zmat, smat)
      IMPLICIT NONE
      INTEGER, INTENT(IN)            :: id, nk, jspin
      INTEGER, INTENT(OUT), OPTIONAL  :: neig
      REAL, INTENT(OUT), OPTIONAL  :: eig(:)
      INTEGER, INTENT(IN), OPTIONAL   :: list(:)
      TYPE(t_mat), OPTIONAL  :: zmat, smat

      INTEGER::nrec, arrayStart, arrayStop, i
      INTEGER, ALLOCATABLE :: ind(:)
      TYPE(t_data_mem), POINTER:: d
      CALL priv_find_data(id, d)

      nrec = nk + (jspin - 1)*d%nkpts
      ! data from d%eig_int
      IF (PRESENT(neig)) THEN
         neig = d%eig_int(nrec)
      ENDIF

      !data from d%eig_eig
      IF (PRESENT(eig)) THEN
         eig = 0.0
         eig = d%eig_eig(:SIZE(eig), nrec)
      ENDIF

      !data from d%eig_vec

      IF (PRESENT(zmat)) THEN
         IF (PRESENT(list)) THEN
            ind = list
         ELSE
            ALLOCATE (ind(zmat%matsize2))
            ind = [(i, i=1, SIZE(ind))]
         END IF
         IF (zmat%l_real) THEN
            IF (.NOT. ALLOCATED(d%eig_vecr)) THEN
               IF (.NOT. ALLOCATED(d%eig_vecc)) CALL juDFT_error("BUG: can not read real/complex vectors from memory")
               DO i = 1, SIZE(ind)
                  arrayStart = (ind(i) - 1)*zMat%matsize1 + 1
                  arrayStop = ind(i)*zMat%matsize1
                  zmat%data_r(:, i) = REAL(d%eig_vecc(arrayStart:arrayStop, nrec))
               ENDDO
            ELSE
               DO i = 1, SIZE(ind)
                  arrayStart = (ind(i) - 1)*zMat%matsize1 + 1
                  arrayStop = ind(i)*zMat%matsize1
                  zmat%data_r(:, i) = d%eig_vecr(arrayStart:arrayStop, nrec)
               ENDDO
            ENDIF
         ELSE !TYPE is (COMPLEX)
            IF (.NOT. ALLOCATED(d%eig_vecc)) CALL juDFT_error("BUG: can not read complex vectors from memory", calledby="eig66_mem")
            DO i = 1, SIZE(ind)
               arrayStart = (ind(i) - 1)*zMat%matsize1 + 1
               arrayStop = ind(i)*zMat%matsize1
               zmat%data_c(:, i) = d%eig_vecc(arrayStart:arrayStop, nrec)
            END DO
         END IF
      ENDIF

      !data from d%eig_vec

      IF (PRESENT(smat)) THEN
         IF (PRESENT(list)) THEN
            ind = list
         ELSE
            ALLOCATE (ind(smat%matsize2))
            ind = [(i, i=1, SIZE(ind))]
         END IF
         IF (smat%l_real) THEN
            IF (.NOT. ALLOCATED(d%olap_r)) THEN
               IF (.NOT. ALLOCATED(d%olap_c)) CALL juDFT_error("BUG: can not read real/complex vectors from memory")
               DO i = 1, SIZE(ind)
                  arrayStart = (ind(i) - 1)*smat%matsize1 + 1
                  arrayStop = ind(i)*smat%matsize1
                  smat%data_r(:, i) = REAL(d%olap_c(arrayStart:arrayStop, nrec))
               ENDDO
            ELSE
               DO i = 1, SIZE(ind)
                  arrayStart = (ind(i) - 1)*smat%matsize1 + 1
                  arrayStop = ind(i)*smat%matsize1
                  smat%data_r(:, i) = d%olap_r(arrayStart:arrayStop, nrec)
               ENDDO
            ENDIF
         ELSE !TYPE is (COMPLEX)
            IF (.NOT. ALLOCATED(d%olap_c)) CALL juDFT_error("BUG: can not read complex vectors from memory", calledby="eig66_mem")
            DO i = 1, SIZE(ind)
               arrayStart = (ind(i) - 1)*smat%matsize1 + 1
               arrayStop = ind(i)*smat%matsize1
               smat%data_c(:, i) = d%olap_c(arrayStart:arrayStop, nrec)
            END DO
         END IF
      ENDIF
   END SUBROUTINE read_eig

   SUBROUTINE write_eig(id, nk, jspin, neig, neig_total, eig, n_size, n_rank, zmat, smat)
      INTEGER, INTENT(IN)          :: id, nk, jspin
      INTEGER, INTENT(IN), OPTIONAL :: n_size, n_rank
      INTEGER, INTENT(IN), OPTIONAL :: neig, neig_total
      REAL, INTENT(IN), OPTIONAL :: eig(:)
      TYPE(t_mat), INTENT(IN), OPTIONAL :: zmat, smat
      INTEGER::nrec
      TYPE(t_data_mem), POINTER:: d
      CALL priv_find_data(id, d)

      nrec = nk + (jspin - 1)*d%nkpts
      ! data from d%eig_int
      IF (PRESENT(neig)) THEN
         IF (PRESENT(neig_total)) THEN
            IF (neig .NE. neig_total) STOP "BUG in eig_mem"
            d%eig_int(nrec) = neig_total
         ELSE
            STOP "BUG2 in eig_mem"
         ENDIF
      ENDIF

      !data from d%eig_eig
      IF (PRESENT(eig)) THEN
         d%eig_eig(:SIZE(eig), nrec) = eig
      ENDIF
      !data from d%eig_vec
      IF (PRESENT(zmat)) THEN
         IF (zmat%l_real) THEN
            IF (.NOT. ALLOCATED(d%eig_vecr)) THEN
               IF (.NOT. ALLOCATED(d%eig_vecc)) CALL juDFT_error("BUG: can not write complex vectors to memory")
               d%eig_vecc(:SIZE(zmat%data_r), nrec) = RESHAPE(CMPLX(zmat%data_r), [SIZE(zmat%data_r)]) !Type cast here
            ELSE
               d%eig_vecr(:SIZE(zmat%data_r), nrec) = RESHAPE(REAL(zmat%data_r), [SIZE(zmat%data_r)])
            ENDIF
         ELSE
            IF (.NOT. ALLOCATED(d%eig_vecc)) CALL juDFT_error("BUG: can not write complex vectors to memory")
            d%eig_vecc(:SIZE(zmat%data_c), nrec) = RESHAPE(zmat%data_c, [SIZE(zmat%data_c)])
         END IF
      ENDIF

      IF (PRESENT(smat)) THEN
         IF (smat%l_real) THEN
            IF (.NOT. ALLOCATED(d%olap_r)) THEN
               IF (.NOT. ALLOCATED(d%olap_c)) CALL juDFT_error("BUG: can not write complex vectors to memory (olap)")
               d%olap_c(:SIZE(smat%data_r), nrec) = RESHAPE(CMPLX(smat%data_r), [SIZE(smat%data_r)]) !Type cast here
            ELSE
               d%olap_r(:SIZE(smat%data_r), nrec) = RESHAPE(REAL(smat%data_r), [SIZE(smat%data_r)])
            ENDIF
         ELSE
            IF (.NOT. ALLOCATED(d%olap_c)) CALL juDFT_error("BUG: can not write complex vectors to memory (olap)")
            d%olap_c(:SIZE(smat%data_c), nrec) = RESHAPE(smat%data_c, [SIZE(smat%data_c)])
         END IF
      ENDIF
   END SUBROUTINE write_eig
END MODULE m_eig66_mem
