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
  SUBROUTINE hsohelp(atoms,sym,input,lapw,nsz, cell,&
       zmat,usdus, zso,noco ,nococonv,&
       nat_start,nat_stop,nat_l,ahelp,bhelp,chelp)
    !
    USE m_abcof
    
    USE m_types
#ifdef CPP_MPI
    use mpi 
#endif
    IMPLICIT NONE
#ifdef CPP_MPI
    INTEGER ierr(3)
#endif
    
     
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_noco),INTENT(IN)        :: noco
    TYPE(t_nococonv),INTENT(IN)    :: nococonv
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
    INTEGER, INTENT (IN) :: nsz(input%jspins)  
    COMPLEX, INTENT (INOUT) :: zso(:,:,:)!lapw%dim_nbasfcn(),2*input%neig,input%jspins)
    COMPLEX, INTENT (OUT):: ahelp(atoms%lmaxd*(atoms%lmaxd+2),nat_l,input%neig,input%jspins)
    COMPLEX, INTENT (OUT):: bhelp(atoms%lmaxd*(atoms%lmaxd+2),nat_l,input%neig,input%jspins)
    COMPLEX, INTENT (OUT):: chelp(-atoms%llod :atoms%llod, input%neig,atoms%nlod,nat_l,input%jspins)
    TYPE(t_mat),INTENT(IN)      :: zmat(:) ! (lapw%dim_nbasfcn(),input%neig,input%jspins)
    !-odim
    !+odim
    !     ..
    !     .. Locals ..
    INTEGER ispin ,n ,na,ie,lmd
    COMPLEX, ALLOCATABLE :: acof(:,:,:),bcof(:,:,:)
    !
    ! turn off the non-collinear part of abcof
    !
    lmd = atoms%lmaxd*(atoms%lmaxd+2)
    
    if ((nat_l)==0) return !nothing to be done here

    chelp(:,:,:,:,input%jspins) = CMPLX(0.0,0.0)
    ALLOCATE ( acof(input%neig,0:lmd,nat_l),bcof(input%neig,0:lmd,nat_l) )
    DO ispin = 1, input%jspins
          CALL abcof(input,atoms,sym,cell,lapw,nsz(ispin),&
          usdus,noco,nococonv,ispin ,&
          acof,bcof,chelp(-atoms%llod:,:,:,:,ispin),zmat(ispin),nat_start=nat_start,nat_stop=nat_stop)
          !
          ! transfer (a,b)cofs to (a,b)helps used in hsoham
          !
          DO ie = 1, input%neig
             DO na = 1, nat_l
               ahelp(:,na,ie,ispin) = acof(ie,1:lmd,na)
               bhelp(:,na,ie,ispin) = bcof(ie,1:lmd,na)
             ENDDO
          ENDDO
       
       !      write(54,'(6f15.8)')(((chelp(m,ie,1,na,1),m=-1,1),ie=1,5),na=1,2)
       !      write(54,'(8f15.8)')(((acof(ie,l,na),l=0,3),ie=1,5),na=1,2)
    ENDDO    ! end of spin loop (ispin)
    !
    DEALLOCATE ( acof,bcof )
    RETURN
  END SUBROUTINE hsohelp
END MODULE m_hsohelp
