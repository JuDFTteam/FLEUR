MODULE m_hsmt_socinit
  USE m_juDFT
  IMPLICIT NONE
  TYPE t_rsoc
     REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: rsopp,rsoppd,rsopdp,rsopdpd     !(atoms%ntype,atoms%lmaxd,2,2)
     REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: rsoplop,rsoplopd,rsopdplo,rsopplo!(atoms%ntype,atoms%nlod,2,2)
     REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:) :: rsoploplop !(atoms%ntype,atoms%nlod,nlod,2,2)
  END TYPE t_rsoc


CONTAINS
  SUBROUTINE hsmt_socinit(mpi,atoms,sphhar,enpara,input,vr,noco,& !in
       rsoc,usdus,isigma) !out
    !Initialized the radial-spin-orbit elements in rsoc and the pauli matrix
    !needed in hssphn for first variation SOC
    USE m_soinit
    USE m_types
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)        :: mpi
    TYPE(t_input),INTENT(IN)      :: input
    TYPE(t_noco),INTENT(IN)       :: noco
    TYPE(t_sphhar),INTENT(IN)     :: sphhar
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_enpara),INTENT(IN)     :: enpara
    !     ..
    !     .. Scalar Arguments ..
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: vr(:,0:,:,:)!(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,DIMENSION%jspd)
    TYPE(t_usdus),INTENT(INOUT):: usdus
    TYPE(t_rsoc),INTENT(OUT) :: rsoc
    COMPLEX, INTENT (OUT)    :: isigma(2,2,3)

    !     ..
    !     .. Local Scalars ..
    INTEGER l,n
    ! for Spin-orbit...
    LOGICAL :: l_test
    LOGICAL, SAVE :: first_k = .TRUE.


    CHARACTER*3 :: chntype


    ALLOCATE(rsoc%rsopp(atoms%ntype,atoms%lmaxd,2,2),rsoc%rsoppd (atoms%ntype,atoms%lmaxd,2,2))
    ALLOCATE(rsoc%rsopdp (atoms%ntype,atoms%lmaxd,2,2),rsoc%rsopdpd(atoms%ntype,atoms%lmaxd,2,2))
    ALLOCATE(rsoc%rsoplop (atoms%ntype,atoms%nlod,2,2),rsoc%rsoplopd(atoms%ntype,atoms%nlod,2,2))
    ALLOCATE(rsoc%rsopdplo(atoms%ntype,atoms%nlod,2,2),rsoc%rsopplo (atoms%ntype,atoms%nlod,2,2))
    ALLOCATE(rsoc%rsoploplop(atoms%ntype,atoms%nlod,atoms%nlod,2,2))

    !     isigma= -i * sigma, where sigma is Pauli matrix
    isigma=CMPLX(0.0,0.0)
    isigma(1,2,1)=CMPLX(0.0,-1.0)
    isigma(2,1,1)=CMPLX(0.0,-1.0)
    isigma(1,2,2)=CMPLX(-1.0,0.0)
    isigma(2,1,2)=CMPLX(1.0,0.0)
    isigma(1,1,3)=CMPLX(0.0,-1.0)
    isigma(2,2,3)=CMPLX(0.0,1.0)

    CALL soinit(atoms,input,enpara,vr,noco%soc_opt(atoms%ntype+2),&
         rsoc%rsopp,rsoc%rsoppd,rsoc%rsopdp,rsoc%rsopdpd,usdus,&
         rsoc%rsoplop,rsoc%rsoplopd,rsoc%rsopdplo,rsoc%rsopplo,rsoc%rsoploplop)
    INQUIRE(file="socscale",exist=l_test)
    IF (l_test) THEN
       OPEN(99,file="socscale")
       READ(99,*) n
       CLOSE(99)
       WRITE(*,*) "SOC scaled by ",n,"%"
       rsoc%rsopp(:,:,:,:) = n/100.* rsoc%rsopp
       rsoc%rsopdp(:,:,:,:) =  n/100.*rsoc%rsopdp
       rsoc%rsoppd(:,:,:,:) =  n/100.*rsoc%rsoppd
       rsoc%rsopdpd(:,:,:,:) =  n/100.*rsoc%rsopdpd
       rsoc%rsoplop(:,:,:,:) =  n/100.*rsoc%rsoplop
       rsoc%rsoplopd(:,:,:,:) =  n/100.*rsoc%rsoplopd
       rsoc%rsopdplo(:,:,:,:) =  n/100.*rsoc%rsopdplo
       rsoc%rsopplo(:,:,:,:) =  n/100.* rsoc%rsopplo
       rsoc%rsoploplop(:,:,:,:,:) = n/100.*rsoc%rsoploplop
    ENDIF
    IF (noco%soc_opt(atoms%ntype+1)) THEN
       DO n= 1,atoms%ntype
          IF (.NOT. noco%soc_opt(n)) THEN 
             rsoc%rsopp(n,:,:,:) = 0.0
             rsoc%rsopdp(n,:,:,:) = 0.0
             rsoc%rsoppd(n,:,:,:) = 0.0
             rsoc%rsopdpd(n,:,:,:) = 0.0
             rsoc%rsoplop(n,:,:,:) = 0.0
             rsoc%rsoplopd(n,:,:,:) = 0.0
             rsoc%rsopdplo(n,:,:,:) = 0.0
             rsoc%rsopplo(n,:,:,:) = 0.0
             rsoc%rsoploplop(n,:,:,:,:) = 0.0
          ENDIF
       ENDDO
    ENDIF

    IF ((first_k).AND.(mpi%irank.EQ.0)) THEN
       DO n = 1,atoms%ntype
          WRITE (6,FMT=8000)
          WRITE (6,FMT=9000)
          WRITE (6,FMT=8001) (2*rsoc%rsopp(n,l,1,1),l=1,3)
          WRITE (6,FMT=8001) (2*rsoc%rsopp(n,l,2,2),l=1,3)
          WRITE (6,FMT=8001) (2*rsoc%rsopp(n,l,2,1),l=1,3)
       ENDDO
       IF (noco%soc_opt(atoms%ntype+1)) THEN
          WRITE (chntype,'(i3)') atoms%ntype
          WRITE (6,fmt='(A,2x,'//chntype//'l1)') 'SOC contribution of certain atom types included in Hamiltonian:',&
               (noco%soc_opt(n),n=1,atoms%ntype)
       ELSE
          WRITE(6,fmt='(A,1x,A)') 'SOC contribution of all atom types included in Hamiltonian.'
       ENDIF
       IF (noco%soc_opt(atoms%ntype+2)) THEN
          WRITE(6,fmt='(A)') 'SOC Hamiltonian is constructed by neglecting B_xc.'
       ENDIF
       first_k=.FALSE.
    ENDIF
8000 FORMAT (' spin - orbit parameter HR  ')
8001 FORMAT (8f8.4)
9000 FORMAT (5x,' p ',5x,' d ', 5x, ' f ')


    RETURN
  END SUBROUTINE hsmt_socinit
END MODULE m_hsmt_socinit
