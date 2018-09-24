!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hsohelp
  !
  !*********************************************************************
  ! preparation of spin-orbit matrix elements: ahelp, bhelp
  ! ahelp(i,n,l,m,jspin) =Sum_(G) (conj(c(G,i,jspin)*a(G,n,l,m,jspin))
  ! bhelp - same a|->b
  ! Original version replaced by a call to abcof. Maybe not so efficient
  ! but includes now LO's and could also be used for noco
  !                                                        gb`02
  !*********************************************************************
  !
CONTAINS
  SUBROUTINE hsohelp(DIMENSION,atoms,sym,input,lapw,nsz, cell,&
       zmat,usdus, zso,noco,oneD,&
       nat_start,nat_stop,nat_l,ahelp,bhelp,chelp)
    !
    USE m_abcof_soc
    USE m_types
    IMPLICIT NONE
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER ierr(3)
#endif
    TYPE(t_dimension),INTENT(IN)   :: DIMENSION
    TYPE(t_oneD),INTENT(IN)        :: oneD
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_noco),INTENT(IN)        :: noco
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_cell),INTENT(IN)        :: cell
    TYPE(t_atoms),INTENT(IN)       :: atoms
    TYPE(t_usdus),INTENT(IN)       :: usdus
    TYPE(t_lapw),INTENT(IN)        :: lapw
    !     ..
    !     .. Scalar Arguments ..
    !     ..
    INTEGER, INTENT (IN) :: nat_start,nat_stop,nat_l
    !     .. Array Arguments ..
    INTEGER, INTENT (IN) :: nsz(DIMENSION%jspd)  
    COMPLEX, INTENT (INOUT) :: zso(:,:,:)!DIMENSION%nbasfcn,2*DIMENSION%neigd,DIMENSION%jspd)
    COMPLEX, INTENT (OUT):: ahelp(atoms%lmaxd*(atoms%lmaxd+2),nat_l,DIMENSION%neigd,input%jspins)
    COMPLEX, INTENT (OUT):: bhelp(atoms%lmaxd*(atoms%lmaxd+2),nat_l,DIMENSION%neigd,input%jspins)
    COMPLEX, INTENT (OUT):: chelp(-atoms%llod :atoms%llod, DIMENSION%neigd,atoms%nlod,nat_l,input%jspins)
    TYPE(t_mat),INTENT(IN)      :: zmat(:) ! (DIMENSION%nbasfcn,DIMENSION%neigd,DIMENSION%jspd)
    !-odim
    !+odim
    !     ..
    !     .. Locals ..
    TYPE(t_atoms)   :: atoms_local
    TYPE(t_noco)    :: noco_local
    TYPE(t_mat)     :: zMat_local
    INTEGER ispin ,l,n ,na,ie,lm,ll1,nv1(DIMENSION%jspd),m,lmd
    INTEGER, ALLOCATABLE :: g1(:,:),g2(:,:),g3(:,:)
    COMPLEX, ALLOCATABLE :: acof(:,:,:),bcof(:,:,:)
    !
    ! turn off the non-collinear part of abcof
    !
    noco_local=noco
    noco_local%l_ss   = .FALSE.
    lmd = atoms%lmaxd*(atoms%lmaxd+2)
    noco_local%qss(:) = 0.0
    atoms_local=atoms
    atoms_local%ngopr(:) = 1 ! use unrotated coeffs...
    !
    ! some praparations to match array sizes
    !
    nv1(1) = lapw%nv(1) ; nv1(DIMENSION%jspd) = lapw%nv(1)
    ALLOCATE (g1(DIMENSION%nvd,DIMENSION%jspd))
    ALLOCATE (g2(DIMENSION%nvd,DIMENSION%jspd))
    ALLOCATE (g3(DIMENSION%nvd,DIMENSION%jspd))
    g1 = 0 ; g2 = 0 ; g3 = 0
    g1(:SIZE(lapw%k1,1),1) = lapw%k1(:SIZE(lapw%k1,1),1) ; g1(:SIZE(lapw%k1,1),DIMENSION%jspd) = lapw%k1(:SIZE(lapw%k1,1),1)
    g2(:SIZE(lapw%k1,1),1) = lapw%k2(:SIZE(lapw%k1,1),1) ; g2(:SIZE(lapw%k1,1),DIMENSION%jspd) = lapw%k2(:SIZE(lapw%k1,1),1)
    g3(:SIZE(lapw%k1,1),1) = lapw%k3(:SIZE(lapw%k1,1),1) ; g3(:SIZE(lapw%k1,1),DIMENSION%jspd) = lapw%k3(:SIZE(lapw%k1,1),1)

    chelp(:,:,:,:,input%jspins) = CMPLX(0.0,0.0)

    ALLOCATE ( acof(DIMENSION%neigd,0:lmd,nat_l),bcof(DIMENSION%neigd,0:lmd,nat_l) )
    DO ispin = 1, input%jspins
       IF (zmat(1)%l_real.AND.noco%l_soc) THEN
          zso(:,1:DIMENSION%neigd,ispin) = CMPLX(zmat(ispin)%data_r(:,1:DIMENSION%neigd),0.0)
          zMat_local%l_real = .FALSE.
          zMat_local%matsize1 = zmat(1)%matsize1
          zMat_local%matsize2 = DIMENSION%neigd
          ALLOCATE(zMat_local%data_c(zmat(1)%matsize1,DIMENSION%neigd))
          zMat_local%data_c(:,:) = zso(:,1:DIMENSION%neigd,ispin)
          CALL abcof_soc(input,atoms_local,sym,cell,lapw,nsz(ispin),&
               usdus, noco_local,ispin,oneD,nat_start,nat_stop,nat_l,&
               acof,bcof,chelp(-atoms%llod:,:,:,:,ispin),zMat_local)
          DEALLOCATE(zMat_local%data_c)
          !
          !
          ! transfer (a,b)cofs to (a,b)helps used in hsoham
          !
          DO ie = 1, DIMENSION%neigd
             DO na = 1, nat_l
                DO l = 1, atoms%lmaxd
                   ll1 = l*(l+1)
                   DO m = -l,l
                      lm = ll1 + m
                      ahelp(lm,na,ie,ispin) = (acof(ie,lm,na))
                      bhelp(lm,na,ie,ispin) = (bcof(ie,lm,na))
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          chelp(:,:,:,:,ispin) = (chelp(:,:,:,:,ispin))
       ELSE
          zMat_local%l_real = zmat(1)%l_real
          zMat_local%matsize1 = zmat(1)%matsize1
          zMat_local%matsize2 = DIMENSION%neigd
          ALLOCATE(zMat_local%data_c(zmat(1)%matsize1,DIMENSION%neigd))
          zMat_local%data_c(:,:) = zmat(ispin)%data_c(:,:)
          CALL abcof_soc(input,atoms_local,sym,cell,lapw,nsz(ispin),&
               usdus, noco_local,ispin,oneD,nat_start,nat_stop,nat_l,&
               acof,bcof,chelp(-atoms%llod:,:,:,:,ispin),zMat_local)
          DEALLOCATE(zMat_local%data_c)
          !
          ! transfer (a,b)cofs to (a,b)helps used in hsoham
          !
          DO ie = 1, DIMENSION%neigd
             DO na = 1, nat_l
                DO l = 1, atoms%lmaxd
                   ll1 = l*(l+1)
                   DO m = -l,l
                      lm = ll1 + m
                      ahelp(lm,na,ie,ispin) = (acof(ie,lm,na))
                      bhelp(lm,na,ie,ispin) = (bcof(ie,lm,na))
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       !      write(54,'(6f15.8)')(((chelp(m,ie,1,na,1),m=-1,1),ie=1,5),na=1,2)
       !      write(54,'(8f15.8)')(((acof(ie,l,na),l=0,3),ie=1,5),na=1,2)
    ENDDO    ! end of spin loop (ispin)
    !
    DEALLOCATE ( acof,bcof,g1,g2,g3 )
    RETURN
  END SUBROUTINE hsohelp
END MODULE m_hsohelp
