!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_types_rsoc
   use m_judft
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_rsoc
 
  TYPE t_rsoc
    REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: rsopp,rsoppd,rsopdp,rsopdpd     !(atoms%ntype,atoms%lmaxd,2,2)
    REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: rsoplop,rsoplopd,rsopdplo,rsopplo!(atoms%ntype,atoms%nlod,2,2)
    REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:) :: rsoploplop !(atoms%ntype,atoms%nlod,nlod,2,2)
    COMPLEX,ALLOCATABLE,DIMENSION(:,:,:,:,:,:)::soangl
    real,allocatable :: rso(:,:,:,:,:,:) ! (icof,jcof,n_type,lmax,2,2)
    contains 
      procedure :: angles   
      procedure :: rad_matrix
      procedure :: init 
   END TYPE t_rsoc

  CONTAINS
   subroutine init(this,atoms)
      use m_types_atoms
       implicit none
     class(t_rsoc),INTENT(INOUT):: this
     class(t_atoms),INTENT(IN)   :: atoms
     ALLOCATE(this%soangl(atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,2,&
         atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,2),source=cmplx(0.0,0.0))
     allocate(this%rso(atoms%max_radial_functions,atoms%max_radial_functions,atoms%ntype,0:atoms%lmaxd,2,2),source=0.0)
   end subroutine init

  subroutine rad_matrix(rsoc,atoms,noco,nococonv,input,fmpi, enpara, vtot)
    !USE m_sorad
    USE m_constants
    USE m_types_atoms
    USE m_types_input
    USE m_types_noco
    USE m_types_nococonv
    USE m_types_potden
    USE m_types_enpara      
    USE m_types_radfun
    USE m_types_mpi
    use m_sointg
    IMPLICIT NONE

    TYPE(t_mpi),INTENT(IN)      :: fmpi
    TYPE(t_enpara),INTENT(IN)   :: enpara
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_nococonv),INTENT(IN) :: nococonv
    TYPE(t_atoms),INTENT(IN)    :: atoms
    type(t_potden),INTENT(IN)   :: vtot
    class(t_rsoc),INTENT(INOUT)    :: rsoc
  
    !     ..
    !     .. Local Scalars ..
    INTEGER:: n,i,j,l,itype,ispin,jspin,ispin1,jspin1
    LOGICAL, SAVE :: first_k = .TRUE.
    TYPE(t_radfun) :: radfun
    REAL,ALLOCATABLE:: v0(:)
    REAL,ALLOCATABLE:: vso(:,:)
    REAL:: e

    
    allocate(vso(atoms%jmtd,2))
   
    !Calculate radial soc-matrix elements
    DO itype = 1,atoms%ntype
       !
       !---> in case of jspins=1
       !
       call radfun%init(atoms, input, itype)
       call radfun%generate_radial_functions( atoms, input, enpara, fmpi, vtot, iType)
       
       !Spin averaged potential
       v0= (vtot%mt(:,0,itype,1)+vtot%mt(:,0,itype,min(2,input%jspins)))/2.
         
       DO l=0,atoms%lmax(itype)
         !
         !---> common spin-orbit integrant V   (average spin directions)
         !                                  SO
         e = (enpara%el0(l,itype,1)+enpara%el0(l,itype,min(2,input%jspins)))/2.

         CALL sointg(itype,e,vtot%mt(:,0,itype,:),v0,atoms,input,vso)
         IF (noco%l_spav) THEN
            DO i= 1,atoms%jmtd
               vso(i,1)= (vso(i,1)+vso(i,2))/2.
               vso(i,2)= vso(i,1)
            ENDDO
         ENDIF

         !                        s       s'            .s       s'
         !-->  radial integrals <u  |V  |u  > = rso(1,1), <u  |V  |u  > = rso(2,1) etc.
         !                            SO                     SO

         IF (l.GT.0) THEN ! there is no spin-orbit for s-states
            DO ispin = 1, 2
               DO jspin = 1, 2
                  ispin1=min(input%jspins,ispin)
                  jspin1=min(input%jspins,jspin)
                   
                  DO i=1,radfun%n_r(l)
                     DO j=1,radfun%n_r(l)
                        rsoc%rso(i,j,itype,l,ispin,jspin) = radso( &
                           radfun%r(:atoms%jri(itype),1,i,l,ispin1), &
                           radfun%r(:atoms%jri(itype),1,j,l,jspin1), &
                           (vso(:atoms%jri(itype),ispin1)+vso(:atoms%jri(itype),jspin1))*0.5, &
                           atoms%dx(itype), &
                           atoms%rmsh(1,itype))
                     ENDDO
                  enddo
               ENDDO
            ENDDO
         ENDIF ! l>0
      END DO
    ENDDO   

    !Scale SOC
    DO n= 1,atoms%ntype
       IF (ABS(noco%socscale(n)-1)>1E-5) THEN
          IF (fmpi%irank==0) WRITE(oUnit,"(a,i0,a,f10.8)") "Scaled SOC for atom ",n," by ",noco%socscale(n)
          rsoc%rso(:,:,n,:,:,:)    = rsoc%rso(:,:,n,:,:,:)*noco%socscale(n)
       ENDIF
    ENDDO

    !DO some IO into out file
    IF ((first_k).AND.(fmpi%irank.EQ.0)) THEN
       DO n = 1,atoms%ntype
          WRITE (oUnit,FMT=8000)
          WRITE (oUnit,FMT=9000)
          WRITE (oUnit,FMT=8001) (2*rsoc%rso(1,1,n,l,1,1),l=1,3)
          WRITE (oUnit,FMT=8001) (2*rsoc%rso(1,1,n,l,2,2),l=1,3)
          WRITE (oUnit,FMT=8001) (2*rsoc%rso(1,1,n,l,2,1),l=1,3)
       ENDDO
       IF (noco%l_spav) WRITE(oUnit,fmt='(A)') 'SOC Hamiltonian is constructed by neglecting B_xc.'
       first_k=.FALSE.
    ENDIF
8000 FORMAT (' spin - orbit parameter HR  ')
8001 FORMAT (8f8.4)
9000 FORMAT (5x,' p ',5x,' d ', 5x, ' f ')
    !
 contains 
   REAL FUNCTION radso(a,b,vso,dx,r0)
    !
    !     compute radial spin-orbit integrals
    !
    USE m_intgr, ONLY : intgr0
    IMPLICIT NONE
    !
    !     .. Scalar Arguments ..
    REAL,    INTENT (IN) :: r0,dx
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: a(:),b(:),vso(:)
    !     ..
    !     .. Local Arrays ..
    REAL q(size(a))
    !     ..
    q = a*b*vso
    CALL intgr0(q,r0,dx,size(a),radso)

    RETURN
  END FUNCTION radso
    end subroutine rad_matrix




  subroutine angles(this,atoms,fmpi,theta,phi,compo)
    USE m_constants
    USE m_anglso
    USE m_sgml
    !USE m_sorad
    USE m_types_atoms
    USE m_types_mpi
    IMPLICIT NONE
    class(t_rsoc),INTENT(INOUT):: this
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_mpi),INTENT(IN)      :: fmpi
    REAL,INTENT(IN)             :: theta,phi
    INTEGER, INTENT(IN),OPTIONAL :: compo
    !     ..
    !     ..
    !     .. Local Scalars ..
    INTEGER is1,is2,jspin1,jspin2,l,l1,l2,m1,m2,n
    !     ..
    !     .. Local Arrays ..
    INTEGER,PARAMETER:: ispjsp(2) = [1,-1]
    

    IF ((ABS(theta).LT.0.00001).AND.(ABS(phi).LT.0.00001)&
                       .AND..NOT.PRESENT(compo)) THEN
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
                         this%soangl(l1,m1,jspin1,l2,m2,jspin2) =&
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
                         this%soangl(l1,m1,jspin1,l2,m2,jspin2) =&
                           anglso(theta,phi,l1,m1,is1,l2,m2,is2,compo)
                      ENDDO
                   ENDDO
                   !
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       !
    ENDIF

    IF (fmpi%irank.EQ.0) THEN
       WRITE (oUnit,FMT=8002)
       DO jspin1 = 1,2
          DO jspin2 = 1,2
             WRITE (oUnit,FMT=*) 'd-states:is1=',jspin1,',is2=',jspin2
             WRITE (oUnit,FMT='(7x,7i8)') (m1,m1=-3,3,1)
             WRITE (oUnit,FMT=8003) (m2, (this%soangl(3,m1,jspin1,3,m2,jspin2),m1=-3,3,1),m2=-3,3,1)
          ENDDO
       ENDDO
    ENDIF
8002 FORMAT (' so - angular matrix elements')
8003 FORMAT (i8,14f8.4)

 
  end subroutine angles 

END MODULE m_types_rsoc