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
  SUBROUTINE hsohelp(DIMENSION,atoms,sym,input,lapw,nsz, cell,bkpt,&
       l_real,z_r,z_c,usdus, zso,noco,oneD, kveclo, ahelp,bhelp,chelp)
    !
    USE m_abcof
    USE m_types
    IMPLICIT NONE
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
    LOGICAL,INTENT(IN) :: l_real
    !     ..
    !     .. Array Arguments ..
    INTEGER, INTENT (IN) :: nsz(DIMENSION%jspd)  
    INTEGER, INTENT (IN) :: kveclo(atoms%nlotot)
    REAL,    INTENT (IN) :: bkpt(3)  
    COMPLEX, INTENT (INOUT) :: zso(DIMENSION%nbasfcn,2*DIMENSION%neigd,DIMENSION%jspd)
    COMPLEX, INTENT (OUT):: ahelp(-atoms%lmaxd:atoms%lmaxd,atoms%lmaxd,atoms%nat,DIMENSION%neigd,DIMENSION%jspd)
    COMPLEX, INTENT (OUT):: bhelp(-atoms%lmaxd:atoms%lmaxd,atoms%lmaxd,atoms%nat,DIMENSION%neigd,DIMENSION%jspd)
    COMPLEX, INTENT (OUT):: chelp(-atoms%llod :atoms%llod, DIMENSION%neigd,atoms%nlod,atoms%nat, DIMENSION%jspd)
    REAL,INTENT(IN)      :: z_r(:,:,:) ! (DIMENSION%nbasfcn,DIMENSION%neigd,DIMENSION%jspd)
    COMPLEX,INTENT(IN)   :: z_c(:,:,:) ! (DIMENSION%nbasfcn,DIMENSION%neigd,DIMENSION%jspd)
    !-odim
    !+odim
    !     ..
    !     .. Locals ..
    TYPE(t_atoms)   :: atoms_local
    TYPE(t_noco)    :: noco_local
    TYPE(t_zMat)    :: zMat_local
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
    ALLOCATE ( g1(DIMENSION%nvd,DIMENSION%jspd),g2(DIMENSION%nvd,DIMENSION%jspd),g3(DIMENSION%nvd,DIMENSION%jspd) )
    g1(:,1) = lapw%k1(:,1) ; g1(:,DIMENSION%jspd) = lapw%k1(:,1)
    g2(:,1) = lapw%k2(:,1) ; g2(:,DIMENSION%jspd) = lapw%k2(:,1)
    g3(:,1) = lapw%k3(:,1) ; g3(:,DIMENSION%jspd) = lapw%k3(:,1)

    chelp(:,:,:,:,input%jspins) = CMPLX(0.0,0.0)

    ALLOCATE ( acof(DIMENSION%neigd,0:lmd,atoms%nat),bcof(DIMENSION%neigd,0:lmd,atoms%nat) )
    DO ispin = 1, input%jspins
       IF (l_real.AND.noco%l_soc) THEN
          zso(:,1:DIMENSION%neigd,ispin) = CMPLX(z_r(:,1:DIMENSION%neigd,ispin),0.0)
          zMat_local%l_real = .FALSE.
          zMat_local%nbasfcn = DIMENSION%nbasfcn
          zMat_local%nbands = DIMENSION%neigd
          ALLOCATE(zMat_local%z_c(DIMENSION%nbasfcn,DIMENSION%neigd))
          zMat_local%z_c(:,:) = zso(:,1:DIMENSION%neigd,ispin)
          CALL abcof(input,atoms_local,DIMENSION%neigd,sym,cell, bkpt,lapw,nsz(ispin),&
               usdus, noco_local,ispin,kveclo,oneD, acof,bcof,chelp(-atoms%llod:,:,:,:,ispin),zMat_local,.false.)
          DEALLOCATE(zMat_local%z_c)
          !
          !
          ! transfer (a,b)cofs to (a,b)helps used in hsoham
          !
          DO ie = 1, DIMENSION%neigd
             DO na = 1, atoms%nat
                DO l = 1, atoms%lmaxd
                   ll1 = l*(l+1)
                   DO m = -l,l
                      lm = ll1 + m
                      ahelp(m,l,na,ie,ispin) = (acof(ie,lm,na))
                      bhelp(m,l,na,ie,ispin) = (bcof(ie,lm,na))
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          chelp(:,:,:,:,ispin) = (chelp(:,:,:,:,ispin))
       ELSE
          zMat_local%l_real = l_real
          zMat_local%nbasfcn = DIMENSION%nbasfcn
          zMat_local%nbands = DIMENSION%neigd
          ALLOCATE(zMat_local%z_c(DIMENSION%nbasfcn,DIMENSION%neigd))
          zMat_local%z_c(:,:) = z_c(:,:,ispin)
          CALL abcof(input,atoms_local,DIMENSION%neigd,sym,cell, bkpt,lapw,nsz(ispin),&
               usdus, noco_local,ispin,kveclo,oneD, acof,bcof,chelp(-atoms%llod:,:,:,:,ispin),zMat_local,.false.)
          DEALLOCATE(zMat_local%z_c)
          !
          ! transfer (a,b)cofs to (a,b)helps used in hsoham
          !
          DO ie = 1, DIMENSION%neigd
             DO na = 1, atoms%nat
                DO l = 1, atoms%lmaxd
                   ll1 = l*(l+1)
                   DO m = -l,l
                      lm = ll1 + m
                      ahelp(m,l,na,ie,ispin) = (acof(ie,lm,na))
                      bhelp(m,l,na,ie,ispin) = (bcof(ie,lm,na))
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
