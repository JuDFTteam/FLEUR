MODULE m_lh_tofrom_lm
   IMPLICIT NONE

CONTAINS
   SUBROUTINE lh_to_lm(sym,atoms, lathar, iType, flh, flm)
      USE m_types

      TYPE(t_sym),INTENT(IN)        :: sym
      TYPE(t_atoms),  INTENT(IN)    :: atoms
      TYPE(t_sphhar), INTENT(IN)    :: lathar
      INTEGER,        INTENT(IN)    :: iType
      REAL,           INTENT(IN)    :: flh(:,0:) ! (iR,iLH)
      COMPLEX,        INTENT(INOUT) :: flm(:,:) ! (iR,lm)

      INTEGER :: iAtom, iLH, ns, l, iM, m, lm, iR, lh, imem, ilm
      INTEGER, ALLOCATABLE :: llhs(:), nmem_work(:), mlh_work(:,:) 
      COMPLEX, ALLOCATABLE :: cMat(:,:), gVec(:,:), clnu_work(:,:)
      REAL, ALLOCATABLE :: fVec(:,:), flh_work(:,:)

      iAtom = SUM(atoms%neq(:iType-1)) + 1
      ns = sym%ntypsy(iAtom)

      flm = CMPLX(0.0,0.0)
      
      ALLOCATE (llhs(atoms%lmax(iType)*(atoms%lmax(iType)+2)))
      ALLOCATE (flh_work(atoms%jri(iType),atoms%lmax(iType)*(atoms%lmax(iType)+2)))
      ALLOCATE (nmem_work(atoms%lmax(iType)*(atoms%lmax(iType)+2)))
      ALLOCATE (mlh_work(lathar%memd,atoms%lmax(iType)*(atoms%lmax(iType)+2)))
      ALLOCATE (clnu_work(lathar%memd,atoms%lmax(iType)*(atoms%lmax(iType)+2)))

      llhs=0.0
      flh_work=0.0
      nmem_work=0
      mlh_work=0
      clnu_work=CMPLX(0.0,0.0)
      
      ilm=0
      DO l=0, atoms%lmax(iType)
         DO m=-l,l
            llhs(ilm)=l
            ilm=ilm+1
         END DO
      END DO 

      ilm=0
      DO lh = 0, atoms%lmax(iType)*(atoms%lmax(iType)+2)
         l=llhs(lh)
         IF (l==lathar%llh(ilm,ns)) THEN
            flh_work(:,lh)=flh(:,ilm)
            nmem_work(lh)=lathar%nmem(ilm,ns)
            DO imem=1, nmem_work(lh)
               mlh_work(imem,lh)=lathar%mlh(imem,ilm,ns)
               clnu_work(imem,lh)=lathar%mlh(imem,ilm,ns)
            END DO
            ilm=ilm+1
         END IF
      END DO

      DO l = 0, atoms%lmax(iType)
         ALLOCATE (cMat(-l:l,-l:l))
         ALLOCATE (fVec(-l:l,atoms%jri(iType)))
         ALLOCATE (gVec(-l:l,atoms%jri(iType)))

         cMat=CMPLX(0.0,0.0)

         DO M = -l, l
            lh = l*(l+1)+M
            fVec(M,:)=flh_work(:,lh)
            DO imem = 1, nmem_work(lh)
               cMat(M,mlh_work(imem,lh))=clnu_work(imem,lh)
            END DO
         END DO

         DO iR = 1, atoms%jri(iType)
            gVec(:,iR)=MATMUL(fVec(:,iR),cMat)
         END DO

         DO m = -l, l
            flm(:,l*(l+1)+1+m)=gVec(m,:)
         END DO

         DEALLOCATE (cMat,fVec,gVec)
      END DO

      DEALLOCATE (llhs,flh_work,nmem_work,mlh_work,clnu_work)

   END SUBROUTINE lh_to_lm

   SUBROUTINE lh_from_lm(sym,atoms, lathar, iType, flm, flh)
      USE m_types
      IMPLICIT NONE
      TYPE(t_sym),INTENT(IN)        :: sym
      TYPE(t_atoms),  INTENT(IN)    :: atoms
      TYPE(t_sphhar), INTENT(IN)    :: lathar
      INTEGER,        INTENT(IN)    :: iType
      COMPLEX,        INTENT(IN)    :: flm(:,:) ! (iR,lm)
      REAL,           INTENT(INOUT) :: flh(:,0:) ! (iR,iLH)


      INTEGER :: iAtom, iLH, ns, l, iM, m, lm, iR, info, lh, imem, info2, ilm
      INTEGER, ALLOCATABLE :: ipiv(:)
      COMPLEX, ALLOCATABLE :: cMat(:,:),fVec(:,:)

      EXTERNAL zgetrf, zgetrs

      iAtom = SUM(atoms%neq(:iType-1)) + 1
      ns = sym%ntypsy(iAtom)

      flh = 0.0

      DO l = 0, atoms%lmax(iType)
         ALLOCATE (cMat(-l:l,-l:l))
         ALLOCATE (fVec(-l:l,atoms%jri(iType)))
         ALLOCATE (ipiv(-l:l))

         cMat=CMPLX(0.0,0.0)

         DO M = -l, l
            lh = l*(l+1)+M
            fVec(M,:)=flm(:,lh+1)
            DO imem = 1, lathar%nmem(lh,ns)
               cMat(M,lathar%mlh(imem,lh,ns))=lathar%clnu(imem,lh,ns)
            END DO
         END DO

         CALL zgetrf(2*l+1,2*l+1,cMat,2*l+1,ipiv,info)

         DO iR = 1, atoms%jri(iType)
            CALL zgetrs('T',2*l+1,1,cMat,2*l+1,ipiv,fVec(:,iR),2*l+1,info2)
         END DO

         DO M = -l, l
            flh(:,l*(l+1)+M)=REAL(fVec(M,:))
         END DO

         DEALLOCATE (cMat,fVec,ipiv)
      END DO

   END SUBROUTINE lh_from_lm

END MODULE m_lh_tofrom_lm
