!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_spnorb
  !*********************************************************************
  !     calls soinit to calculate the radial spin-orbit matrix elements:
  !     rsopp,rsopdpd,rsoppd,rsopdp
  !     and sets up the so - angular matrix elements (soangl)
  !     using the functions anglso and sgml.
  !*********************************************************************
CONTAINS
  SUBROUTINE spnorb(atoms,noco,input,mpi, enpara, vr, usdus, rsoc,l_angles,hub1)
    USE m_sorad 
    USE m_constants, only : hartree_to_ev_const
    USE m_types
    IMPLICIT NONE

    TYPE(t_mpi),INTENT(IN)      :: mpi
    TYPE(t_enpara),INTENT(IN)   :: enpara
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_usdus),INTENT(INOUT) :: usdus
    TYPE(t_rsoc),INTENT(OUT)    :: rsoc
    LOGICAL,INTENT(IN)          :: l_angles
    TYPE(t_hub1ham),OPTIONAL,INTENT(INOUT) :: hub1
    !     ..
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: vr(:,0:,:,:) !(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins)
    !     ..
    !     .. Local Scalars ..
    INTEGER is1,is2,jspin1,jspin2,l,l1,l2,m1,m2,n,i_hia
    LOGICAL, SAVE :: first_k = .TRUE.
    !     ..
  
    !Allocate space for SOC matrix elements; set to zero at the same time
    ALLOCATE(rsoc%rsopp  (atoms%ntype,atoms%lmaxd,2,2));rsoc%rsopp =0.0
    ALLOCATE(rsoc%rsoppd (atoms%ntype,atoms%lmaxd,2,2));rsoc%rsoppd=0.0
    ALLOCATE(rsoc%rsopdp (atoms%ntype,atoms%lmaxd,2,2));rsoc%rsopdp=0.0
    ALLOCATE(rsoc%rsopdpd(atoms%ntype,atoms%lmaxd,2,2));rsoc%rsopdpd=0.0
    ALLOCATE(rsoc%rsoplop (atoms%ntype,atoms%nlod,2,2));rsoc%rsoplop=0.0
    ALLOCATE(rsoc%rsoplopd(atoms%ntype,atoms%nlod,2,2));rsoc%rsoplopd=0.0
    ALLOCATE(rsoc%rsopdplo(atoms%ntype,atoms%nlod,2,2));rsoc%rsopdplo=0.0
    ALLOCATE(rsoc%rsopplo (atoms%ntype,atoms%nlod,2,2));rsoc%rsopplo=0.0
    ALLOCATE(rsoc%rsoploplop(atoms%ntype,atoms%nlod,atoms%nlod,2,2));rsoc%rsoploplop=0.0
    IF (l_angles) ALLOCATE(rsoc%soangl(atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,2,&
         atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,2))

    !Calculate radial soc-matrix elements
    DO n = 1,atoms%ntype
       CALL sorad(atoms,input,n,vr(:,0,n,:),enpara,noco%l_spav,rsoc,usdus)
    END DO
    
    !
    !Scale SOC 
    DO n= 1,atoms%ntype
       IF (ABS(noco%socscale(n)-1)>1E-5) THEN
          IF (mpi%irank==0) WRITE(6,"(a,i0,a,f10.8)") "Scaled SOC for atom ",n," by ",noco%socscale(n)
          rsoc%rsopp(n,:,:,:)    = rsoc%rsopp(n,:,:,:)*noco%socscale(n)
          rsoc%rsopdp(n,:,:,:)   = rsoc%rsopdp(n,:,:,:)*noco%socscale(n)
          rsoc%rsoppd(n,:,:,:)   = rsoc%rsoppd(n,:,:,:)*noco%socscale(n)
          rsoc%rsopdpd(n,:,:,:)  = rsoc%rsopdpd(n,:,:,:)*noco%socscale(n)
          rsoc%rsoplop(n,:,:,:)  = rsoc%rsoplop(n,:,:,:)*noco%socscale(n)
          rsoc%rsoplopd(n,:,:,:) = rsoc%rsoplopd(n,:,:,:)*noco%socscale(n)
          rsoc%rsopdplo(n,:,:,:) = rsoc%rsopdplo(n,:,:,:)*noco%socscale(n)
          rsoc%rsopplo(n,:,:,:)  = rsoc%rsopplo(n,:,:,:)*noco%socscale(n)
          rsoc%rsoploplop(n,:,:,:,:) = rsoc%rsoploplop(n,:,:,:,:)*noco%socscale(n)
       ENDIF
    ENDDO
    
    !Read in SOC-parameter for shell with hubbard 1
    IF(PRESENT(hub1).AND.mpi%irank.EQ.0) THEN
      DO i_hia = 1, atoms%n_hia
         IF(hub1%l_soc_given(i_hia)) CYCLE
         n = atoms%lda_u(atoms%n_u+i_hia)%atomType
         l = atoms%lda_u(atoms%n_u+i_hia)%l
         hub1%xi(i_hia) = (rsoc%rsopp(n,l,1,1)+rsoc%rsopp(n,l,2,2))*hartree_to_ev_const
      ENDDO
    ENDIF

    !DO some IO into out file
      IF ((first_k).AND.(mpi%irank.EQ.0)) THEN
       DO n = 1,atoms%ntype
          WRITE (6,FMT=8000)
          WRITE (6,FMT=9000)
          WRITE (6,FMT=8001) (2*rsoc%rsopp(n,l,1,1),l=1,3)
          WRITE (6,FMT=8001) (2*rsoc%rsopp(n,l,2,2),l=1,3)
          WRITE (6,FMT=8001) (2*rsoc%rsopp(n,l,2,1),l=1,3)
       ENDDO
       IF (noco%l_spav) THEN
          WRITE(6,fmt='(A)') 'SOC Hamiltonian is constructed by neglecting B_xc.'
       ENDIF
       first_k=.FALSE.
    ENDIF
8000 FORMAT (' spin - orbit parameter HR  ')
8001 FORMAT (8f8.4)
9000 FORMAT (5x,' p ',5x,' d ', 5x, ' f ')
    !

    !Calculate angular matrix elements if requested
    IF (l_angles) &
         CALL spnorb_angles(atoms,mpi,noco%theta,noco%phi,rsoc%soangl)
  END SUBROUTINE spnorb

  SUBROUTINE spnorb_angles(atoms,mpi,theta,phi,soangl)
    USE m_anglso
    USE m_sgml
    USE m_sorad 
    USE m_types
    IMPLICIT NONE
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_mpi),INTENT(IN)      :: mpi
    REAL,INTENT(IN)             :: theta,phi
    COMPLEX,INTENT(INOUT)       :: soangl(:,-atoms%lmaxd:,:,:,-atoms%lmaxd:,:)
    !     ..
    !     ..
    !     .. Local Scalars ..
    INTEGER is1,is2,jspin1,jspin2,l,l1,l2,m1,m2,n
    !     ..
    !     .. Local Arrays ..
    INTEGER ispjsp(2)
    !     ..
    !     ..
    DATA ispjsp/1,-1/

  
    IF ((ABS(theta).LT.0.00001).AND.(ABS(phi).LT.0.00001)) THEN
       !
       !       TEST for real function sgml(l1,m1,is1,l2,m2,is2)
       !
       DO l1 = 1,atoms%lmaxd
          DO l2 = 1,atoms%lmaxd
             DO jspin1 = 1,2
                DO jspin2 = 1,2
                   is1=ispjsp(jspin1)
                   is2=ispjsp(jspin2)
                   DO m1 = -l1,l1,1
                      DO m2 = -l2,l2,1
                         soangl(l1,m1,jspin1,l2,m2,jspin2) =&
                              CMPLX(sgml(l1,m1,is1,l2,m2,is2),0.0)
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       
    ELSE
       !
       !       TEST for complex function anglso(teta,phi,l1,m1,is1,l2,m2,is2)
       ! 
       DO l1 = 1,atoms%lmaxd
          DO l2 = 1,atoms%lmaxd
             DO jspin1 = 1,2
                DO jspin2 = 1,2
                   is1=ispjsp(jspin1)
                   is2=ispjsp(jspin2)
                   !
                   DO m1 = -l1,l1,1
                      DO m2 = -l2,l2,1
                         soangl(l1,m1,jspin1,l2,m2,jspin2) =&
                              anglso(theta,phi,l1,m1,is1,l2,m2,is2)
                      ENDDO
                   ENDDO
                   !
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       !
    ENDIF
    
    IF (mpi%irank.EQ.0) THEN
       WRITE (6,FMT=8002) 
       DO jspin1 = 1,2
          DO jspin2 = 1,2
             WRITE (6,FMT=*) 'd-states:is1=',jspin1,',is2=',jspin2
             WRITE (6,FMT='(7x,7i8)') (m1,m1=-3,3,1)
             WRITE (6,FMT=8003) (m2, (soangl(3,m1,jspin1,3,m2,jspin2),&
                  m1=-3,3,1),m2=-3,3,1)
          ENDDO
       ENDDO
    ENDIF
8002 FORMAT (' so - angular matrix elements')
8003 FORMAT (i8,14f8.4)

  END SUBROUTINE spnorb_angles
END MODULE m_spnorb
