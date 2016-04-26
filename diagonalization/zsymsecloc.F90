MODULE m_zsymsecloc
  use m_juDFT
!*******************************************************
!  Solve the generalized secular equation. 
!  For film-systems exhibiting
!  z-reflexion symmetry, the basis is transformed to
!  even and odd functions and the even-even and odd-odd 
!  blocks are diagonalized separately.
!  If local orbitals are present in a film with z-reflection,
!  locrectify is used to construct linear combinations of
!  the local orbitals that are eigenfunctions of the z-
!  reflexion operation.
!  Frank Freimuth, January 2006
!*******************************************************
CONTAINS
  SUBROUTINE zsymsecloc(jsp,input,lapw,bkpt,atoms,&
       kveclo, sym,cell, dimension,matsize,&
       nsize, jij,matind,nred, a,b, z,eig,ne)

#include"cpp_double.h"

    USE m_locrectify
    USE m_geneigprobl
    USE m_types
    IMPLICIT NONE
    TYPE(t_dimension),INTENT(IN)   :: dimension
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_jij),INTENT(IN)   :: jij
    TYPE(t_sym),INTENT(IN)   :: sym
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_atoms),INTENT(IN)   :: atoms

    TYPE(t_lapw),INTENT(IN)   :: lapw
    real,intent(in)   ::bkpt(3)
    integer,intent(in)::kveclo(atoms%nlotot)

    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN)    :: matsize,jsp 
    INTEGER, INTENT (INOUT) :: nsize
    INTEGER, INTENT (OUT)   :: ne
    INTEGER, INTENT (IN)    :: nred
    !     ..
    !     .. Array Arguments ..
    INTEGER, INTENT (IN)  :: matind(dimension%nbasfcn,2)
    REAL,    INTENT (OUT) :: eig(dimension%neigd)
#ifdef CPP_F90

#ifdef CPP_INVERSION
    REAL,  INTENT (INOUT) :: a(:),b(:)
    REAL,  INTENT (INOUT) :: z(:,:)
#else
    COMPLEX, INTENT (INOUT)::a(:),b(:)
    COMPLEX, INTENT (INOUT) :: z(:,:)
#endif

#else

#ifdef CPP_INVERSION
    REAL, ALLOCATABLE, INTENT (INOUT) :: a(:),b(:)
    REAL, ALLOCATABLE, INTENT (INOUT) :: z(:,:)
#else
    COMPLEX, ALLOCATABLE, INTENT (INOUT)::a(:),b(:)
    COMPLEX, ALLOCATABLE, INTENT (INOUT) :: z(:,:)
#endif

#endif

#ifdef CPP_INVERSION
    real locrec(atoms%nlotot,atoms%nlotot)
    real,allocatable::atemp(:,:)
    real,allocatable::btemp(:,:)
    real,allocatable::aa(:)
    real,allocatable::bb(:)
    real,allocatable::z1(:,:)
    real,allocatable::z2(:,:)
    real,allocatable::ztemp(:,:)
    real recsqrtwo,dobzero
#else
    complex locrec(atoms%nlotot,atoms%nlotot)
    complex,allocatable:: atemp(:,:)
    complex,allocatable:: btemp(:,:)
    complex,allocatable::aa(:)
    complex,allocatable::bb(:)
    complex,allocatable::z1(:,:)
    complex,allocatable::z2(:,:)
    complex,allocatable::ztemp(:,:)
    complex recsqrtwo,dobzero
#endif
    logical l_verbose
    logical kindlocrec(atoms%nlotot)
    real,allocatable::etemp1(:),etemp2(:)
    integer iind,jind,ii,info,iu,pos1,pos2,pos3,i2,j2
    integer ne1,ne2,i1,j1,i,j
    logical,allocatable:: evensort(:)
    integer evenlocs,oddlocs
    integer evenlocindex(atoms%nlotot)
    integer oddlocindex(atoms%nlotot)


    !      print*,"in zsymsecloc"
#ifndef CPP_F90
    deallocate(z)
#endif

    !******************************************
    ! l_zref=.false. => simply call eigensolver
    !******************************************
    if(.not.sym%l_zref)then
       call geneigprobl(dimension%nbasfcn, nsize,dimension%neigd,jij%l_j,a,b, z,eig,ne)

#ifndef CPP_F90
       allocate(a(dimension%nbasfcn*(dimension%nbasfcn+1)/2))
       allocate(b(dimension%nbasfcn*(dimension%nbasfcn+1)/2))
#endif
       return
       !******************************************
       ! l_zref=.true. => blockdiagonalize
       ! hamiltonian and overlap-matrix by
       ! transforming to even/odd basisfunctions
       !******************************************
    else
       inquire(file='verbose',exist=l_verbose)
       if(.not.l_verbose)inquire(file='veryverbose',exist=l_verbose)
       iu=dimension%neigd/2
#ifdef CPP_INVERSION
       recsqrtwo=1.0/sqrt(2.0)
       dobzero=0.0
#else
       recsqrtwo=cmplx(1.0/sqrt(2.0),0.0)
       dobzero=cmplx(0.0,0.0)
#endif


       if(atoms%nlotot.gt.0)then
          !*********************************************
          ! Find even/odd linear combinations of locs
          !*********************************************
          if(l_verbose)then
             print*,"find lincos of locs that are eigenfunctions of zreflection"
             print*,"apws=",lapw%nv(jsp)
             print*,"atoms%nlotot=",atoms%nlotot
             print*,"basis-size=",nsize
          endif
          IF(nsize/=(lapw%nv(jsp)+atoms%nlotot)) &
               CALL juDFT_error("nsize /= lapw%nv + atoms%nlotot" ,calledby ="zsymsecloc")
          call locrectify(jsp,input,lapw,bkpt,atoms, kveclo, sym,cell,&
                    locrec,kindlocrec,evenlocs,oddlocs, evenlocindex,oddlocindex)


          !*********************************************
          ! Perform basis-transformation of Hamiltonian
          ! and Overlap matrix. The new basis is the basis
          ! of even / odd (with respect to z-reflection)
          ! local orbitals.
          !*********************************************
          allocate(atemp(lapw%nv(jsp),atoms%nlotot))
          allocate(btemp(atoms%nlotot,lapw%nv(jsp)))
          do iind=1,atoms%nlotot
             pos1=((lapw%nv(jsp)+iind)*(lapw%nv(jsp)+iind-1))/2
             atemp(1:lapw%nv(jsp),iind)=a(pos1+1:pos1+lapw%nv(jsp))
          enddo
#ifdef CPP_INVERSION
          call CPP_BLAS_sgemm('T','T',atoms%nlotot,lapw%nv(jsp),atoms%nlotot,real(1.0),&
                   locrec,atoms%nlotot,atemp,lapw%nv(jsp),real(0.0),btemp,atoms%nlotot)
#else
          call CPP_BLAS_cgemm('C','T',atoms%nlotot,lapw%nv(jsp),atoms%nlotot,cmplx(1.0,0.0),&
                   locrec,atoms%nlotot,atemp,lapw%nv(jsp),cmplx(0.0,0.0),btemp,atoms%nlotot)
#endif
          atemp(:,:) = transpose(btemp)
          do iind=1,atoms%nlotot
             pos1=((lapw%nv(jsp)+iind)*(lapw%nv(jsp)+iind-1))/2
             a(pos1+1:pos1+lapw%nv(jsp))=atemp(1:lapw%nv(jsp),iind)
          enddo

          do iind=1,atoms%nlotot
             pos1=((lapw%nv(jsp)+iind)*(lapw%nv(jsp)+iind-1))/2
             atemp(1:lapw%nv(jsp),iind)=b(pos1+1:pos1+lapw%nv(jsp))
          enddo
#ifdef CPP_INVERSION
          call CPP_BLAS_sgemm('T','T',atoms%nlotot,lapw%nv(jsp),atoms%nlotot,real(1.0),&
                   locrec,atoms%nlotot,atemp,lapw%nv(jsp),real(0.0),btemp,atoms%nlotot)
#else
          call CPP_BLAS_cgemm('C','T',atoms%nlotot,lapw%nv(jsp),atoms%nlotot,cmplx(1.0,0.0),&
                   locrec,atoms%nlotot,atemp,lapw%nv(jsp),cmplx(0.0,0.0),btemp,atoms%nlotot)
#endif
          atemp(:,:) = transpose(btemp)
          do iind=1,atoms%nlotot
             pos1=((lapw%nv(jsp)+iind)*(lapw%nv(jsp)+iind-1))/2
             b(pos1+1:pos1+lapw%nv(jsp))=atemp(1:lapw%nv(jsp),iind)
          enddo
          deallocate(atemp)
          deallocate(btemp)

          allocate(atemp(atoms%nlotot,atoms%nlotot))
          allocate(btemp(atoms%nlotot,atoms%nlotot))
          do iind=1,atoms%nlotot
             pos1=((lapw%nv(jsp)+iind)*(lapw%nv(jsp)+iind-1))/2
             do jind=1,iind-1
                atemp(iind,jind)=a(pos1+lapw%nv(jsp)+jind)
#ifdef CPP_INVERSION
                atemp(jind,iind)=a(pos1+lapw%nv(jsp)+jind)
#else
                atemp(jind,iind)=conjg(a(pos1+lapw%nv(jsp)+jind))
#endif
             enddo
             atemp(iind,iind)=a(pos1+lapw%nv(jsp)+iind)
          enddo
#ifdef CPP_INVERSION
          call CPP_BLAS_sgemm('T','N',atoms%nlotot,atoms%nlotot,atoms%nlotot,real(1.0),&
                    locrec,atoms%nlotot,atemp,atoms%nlotot,real(0.0),btemp,atoms%nlotot)
          call CPP_BLAS_sgemm('N','N',atoms%nlotot,atoms%nlotot,atoms%nlotot,real(1.0),&
                    btemp,atoms%nlotot,locrec,atoms%nlotot,real(0.0),atemp,atoms%nlotot)
#else
          call CPP_BLAS_cgemm('C','N',atoms%nlotot,atoms%nlotot,atoms%nlotot,cmplx(1.0,0.0),&
                    locrec,atoms%nlotot,atemp,atoms%nlotot,cmplx(0.0,0.0),btemp,atoms%nlotot)
          call CPP_BLAS_cgemm('N','N',atoms%nlotot,atoms%nlotot,atoms%nlotot,cmplx(1.0,0.0),&
                    btemp,atoms%nlotot,locrec,atoms%nlotot,cmplx(0.0,0.0),atemp,atoms%nlotot)
#endif
          do iind=1,atoms%nlotot
             pos1=((lapw%nv(jsp)+iind)*(lapw%nv(jsp)+iind-1))/2
             do jind=1,iind
                a(pos1+lapw%nv(jsp)+jind)=atemp(iind,jind)
             enddo
          enddo

          do iind=1,atoms%nlotot
             pos1=((lapw%nv(jsp)+iind)*(lapw%nv(jsp)+iind-1))/2
             do jind=1,iind-1
                atemp(iind,jind)=b(pos1+lapw%nv(jsp)+jind)
#ifdef CPP_INVERSION
                atemp(jind,iind)=b(pos1+lapw%nv(jsp)+jind)
#else
                atemp(jind,iind)=conjg(b(pos1+lapw%nv(jsp)+jind))
#endif
             enddo
             atemp(iind,iind)=b(pos1+lapw%nv(jsp)+iind)
          enddo
#ifdef CPP_INVERSION
          call CPP_BLAS_sgemm('T','N',atoms%nlotot,atoms%nlotot,atoms%nlotot,real(1.0),&
                    locrec,atoms%nlotot,atemp,atoms%nlotot,real(0.0),btemp,atoms%nlotot)
          call CPP_BLAS_sgemm('N','N',atoms%nlotot,atoms%nlotot,atoms%nlotot,real(1.0),&
                    btemp,atoms%nlotot,locrec,atoms%nlotot,real(0.0),atemp,atoms%nlotot)
#else
          call CPP_BLAS_cgemm('C','N',atoms%nlotot,atoms%nlotot,atoms%nlotot,cmplx(1.0,0.0),&
                    locrec,atoms%nlotot,atemp,atoms%nlotot,cmplx(0.0,0.0),btemp,atoms%nlotot)
          call CPP_BLAS_cgemm('N','N',atoms%nlotot,atoms%nlotot,atoms%nlotot,cmplx(1.0,0.0),&
                    btemp,atoms%nlotot,locrec,atoms%nlotot,cmplx(0.0,0.0),atemp,atoms%nlotot)
#endif
          do iind=1,atoms%nlotot
             pos1=((lapw%nv(jsp)+iind)*(lapw%nv(jsp)+iind-1))/2
             do jind=1,iind
                b(pos1+lapw%nv(jsp)+jind)=atemp(iind,jind)
             enddo
          enddo

          deallocate(atemp)
          deallocate(btemp)

       else
          evenlocs=0
          oddlocs=0
       endif ! atoms%nlotot.gt.0

       !*********************************************
       ! Test matind.
       !*********************************************
       ii=0
       jind=0
       pos1=0
       do iind=1,nred
          IF(matind(iind,1)<matind(iind,2))  CALL juDFT_error("mat1mat2" ,calledby ="zsymsecloc")
          IF(matind(iind,1)<pos1) CALL juDFT_error("matpos1",calledby ="zsymsecloc")
          pos1=matind(iind,1)
          if(matind(iind,1).ne.matind(iind,2))then
             ii=ii+1
          else
             jind=jind+1
          endif
       enddo
       IF(2*ii+jind/=lapw%nv(jsp)) CALL juDFT_error("matind",calledby="zsymsecloc")

       !*****************************************************************
       !Transform into representation with even-even- and odd-odd-blocks.
       !First step: Transform the lapw-lapw-part.
       !*****************************************************************
       allocate(    aa(  ((nred+evenlocs)*(nred+evenlocs+1))/2   )  )
       allocate(    bb(  ((nred+evenlocs)*(nred+evenlocs+1))/2   )  )
       i2=0
       j2=0
       DO i=1,nred
          DO j=1,i
             i1=(matind(i,1)-1)*matind(i,1)/2+matind(j,1)
             j1=(matind(i,1)-1)*matind(i,1)/2+matind(j,2)
             i2=i2+1
             aa(i2)=a(i1)+a(j1)
             bb(i2)=b(i1)+b(j1)
             IF ((matind(i,1).NE.matind(i,2)).AND. (matind(j,1).NE.matind(j,2))) THEN
                j2=j2+1
                a(j2)=a(i1)-a(j1)
                b(j2)=b(i1)-b(j1)
             ENDIF
          ENDDO
       ENDDO

       if(atoms%nlotot.gt.0)then
          !******************************************************
          !the lapw-lo- and lo-lo-parts for the even-even-block.
          !******************************************************
          allocate(atemp(lapw%nv(jsp),atoms%nlotot))
          allocate(btemp(nred,evenlocs))
          !hamiltonian
          do iind=1,atoms%nlotot
             pos1=((lapw%nv(jsp)+iind)*(lapw%nv(jsp)+iind-1))/2
             atemp(1:lapw%nv(jsp),iind)=a(pos1+1:pos1+lapw%nv(jsp))
          enddo
          do jind=1,evenlocs
             do iind=1,nred
                btemp(iind,jind)=(atemp(matind(iind,1),evenlocindex(jind))&
                     +atemp(matind(iind,2),evenlocindex(jind)))*recsqrtwo
             enddo
          enddo
          do iind=1,evenlocs
             pos1=((nred+iind)*(nred+iind-1))/2
             aa(pos1+1:pos1+nred)=btemp(1:nred,iind)
          enddo
          !metric
          do iind=1,atoms%nlotot
             pos1=((lapw%nv(jsp)+iind)*(lapw%nv(jsp)+iind-1))/2
             atemp(1:lapw%nv(jsp),iind)=b(pos1+1:pos1+lapw%nv(jsp))
          enddo
          do jind=1,evenlocs
             do iind=1,nred
                btemp(iind,jind)=(atemp(matind(iind,1),evenlocindex(jind))&
                     +atemp(matind(iind,2),evenlocindex(jind)))*recsqrtwo
             enddo
          enddo
          do iind=1,evenlocs
             pos1=((nred+iind)*(nred+iind-1))/2
             bb(pos1+1:pos1+nred)=btemp(1:nred,iind)
          enddo
          deallocate(btemp)
          deallocate(atemp)
          !hamiltonian
          do iind=1,evenlocs
             do jind=1,iind
                pos1=((lapw%nv(jsp)+evenlocindex(iind))*(lapw%nv(jsp)+evenlocindex(iind)-1))/2
                pos2=((nred+iind)*(nred+iind-1))/2
                aa(pos2+nred+jind)=a(pos1+lapw%nv(jsp)+evenlocindex(jind))
             enddo
          enddo
          !metric
          do iind=1,evenlocs
             do jind=1,iind
                pos1=((lapw%nv(jsp)+evenlocindex(iind))*(lapw%nv(jsp)+evenlocindex(iind)-1))/2
                pos2=((nred+iind)*(nred+iind-1))/2
                bb(pos2+nred+jind)=b(pos1+lapw%nv(jsp)+evenlocindex(jind))
             enddo
          enddo


          !******************************************************
          !the lapw-lo- and lo-lo-parts for the odd-odd-block
          !******************************************************
          allocate(atemp(lapw%nv(jsp),atoms%nlotot))
          allocate(btemp(lapw%nv(jsp)-nred,oddlocs))
          !hamiltonian
          do iind=1,atoms%nlotot
             pos1=((lapw%nv(jsp)+iind)*(lapw%nv(jsp)+iind-1))/2
             atemp(1:lapw%nv(jsp),iind)=a(pos1+1:pos1+lapw%nv(jsp))
          enddo
          do jind=1,oddlocs
             ii=0
             do iind=1,nred
                if(matind(iind,1).ne.matind(iind,2))then
                   ii=ii+1
                   btemp(ii,jind)=(atemp(matind(iind,1),oddlocindex(jind))&
                        -atemp(matind(iind,2),oddlocindex(jind)))*recsqrtwo
                endif
             enddo
          enddo
          do iind=1,oddlocs
             pos1=((lapw%nv(jsp)-nred+iind)*(lapw%nv(jsp)-nred+iind-1))/2
             a(pos1+1:pos1+lapw%nv(jsp)-nred)=btemp(1:lapw%nv(jsp)-nred,iind)
          enddo
          !metric
          do iind=1,atoms%nlotot
             pos1=((lapw%nv(jsp)+iind)*(lapw%nv(jsp)+iind-1))/2
             atemp(1:lapw%nv(jsp),iind)=b(pos1+1:pos1+lapw%nv(jsp))
          enddo
          do jind=1,oddlocs
             ii=0
             do iind=1,nred
                if(matind(iind,1).ne.matind(iind,2))then
                   ii=ii+1
                   btemp(ii,jind)=(atemp(matind(iind,1),oddlocindex(jind))&
                        -atemp(matind(iind,2),oddlocindex(jind)))*recsqrtwo
                endif
             enddo
          enddo
          do iind=1,oddlocs
             pos1=((lapw%nv(jsp)-nred+iind)*(lapw%nv(jsp)-nred+iind-1))/2
             b(pos1+1:pos1+lapw%nv(jsp)-nred)=btemp(1:lapw%nv(jsp)-nred,iind)
          enddo
          deallocate(btemp)
          deallocate(atemp)
          !hamiltonian
          do iind=1,oddlocs
             do jind=1,iind
                pos1=((lapw%nv(jsp)+oddlocindex(iind))*(lapw%nv(jsp)+oddlocindex(iind)-1))/2
                pos2=((lapw%nv(jsp)-nred+iind)*(lapw%nv(jsp)-nred+iind-1))/2
                a(pos2+lapw%nv(jsp)-nred+jind)=a(pos1+lapw%nv(jsp)+oddlocindex(jind))
             enddo
          enddo
          !metric
          do iind=1,oddlocs
             do jind=1,iind
                pos1=((lapw%nv(jsp)+oddlocindex(iind))*(lapw%nv(jsp)+oddlocindex(iind)-1))/2
                pos2=((lapw%nv(jsp)-nred+iind)*(lapw%nv(jsp)-nred+iind-1))/2
                b(pos2+lapw%nv(jsp)-nred+jind)=b(pos1+lapw%nv(jsp)+oddlocindex(jind))
             enddo
          enddo



       endif !atoms%nlotot.gt.0

       !******************************************************
       !Solve the eigenvalue problem for the odd-odd block
       !******************************************************

       allocate(etemp2(iu+1))
       call geneigprobl(lapw%nv(jsp)-nred+oddlocs, lapw%nv(jsp)-nred+oddlocs,iu,jij%l_j,a,b, z2,etemp2,ne2)

       if(l_verbose)then   
          print*,"odd block diagonalized"
       endif

       !***********************************************************
       !Solve the eigenvalue problem for the even-even block
       !***********************************************************
       allocate(etemp1(iu+1))
       call geneigprobl(nred+evenlocs, nred+evenlocs,iu,jij%l_j,aa,bb, z1,etemp1,ne1)
       if(l_verbose)then
          print*,"even block diagonalized"
       endif


       ne=ne1+ne2
       !********************************************************************
       !  Recover eigenvectors of original eigenvalue problem.
       !  Sort eigenvalues and eigenvectors according to increasing eigenvalue.
       !  etemp1 holds eigenvalues of even block.
       !  etemp2 holds eigenvalues of odd block.
       !  z1 holds eigenvectors of even block.
       !  z2 holds eigenvectors of odd block.
       !********************************************************************
#ifndef CPP_F90
       allocate(z(dimension%nbasfcn,dimension%neigd))
#endif
       allocate(evensort(ne))
       etemp1(ne1+1)=99.9e9
       etemp2(ne2+1)=99.9e9
       jind=1
       iind=1
       !evensort(ii)=.true.  => eigenvalue ii belongs to even spectrum
       !evensort(ii)=.false. => eigenvalue ii belongs to odd spectrum
       do ii=1,ne
          if(etemp1(iind).lt.etemp2(jind)) then
             evensort(ii)=.true.
             iind=iind+1
          else
             evensort(ii)=.false.
             jind=jind+1
          endif
       enddo
       iind=1 !position in the even-arrays
       jind=1 !position in the oneD%odd-arrays
       do ii=1,ne
          if(evensort(ii))then
             eig(ii)=etemp1(iind)
             z(1:lapw%nv(jsp)+atoms%nlotot,ii)=dobzero
             !Recover the eigenvectors of the original problem for the even block
             do i=1,nred
                z(matind(i,1),ii)=z1(i,iind)*recsqrtwo
                z(matind(i,2),ii)=z(matind(i,2),ii)+z1(i,iind)*recsqrtwo
             enddo !i
             if(atoms%nlotot.gt.0) z(lapw%nv(jsp)+1:lapw%nv(jsp)+atoms%nlotot,ii)=dobzero
             do pos1=1,evenlocs
                do i=1,atoms%nlotot
#ifdef CPP_INVERSION
                   z(lapw%nv(jsp)+i,ii)=z(lapw%nv(jsp)+i,ii)+ z1(nred+pos1,iind)*locrec(i,evenlocindex(pos1))
#else
                   z(lapw%nv(jsp)+i,ii)=z(lapw%nv(jsp)+i,ii)+ z1(nred+pos1,iind)*conjg(locrec(i,evenlocindex(pos1)))
#endif
                enddo
             enddo
             iind=iind+1
          else
             !Recover the eigenvectors of the original problem for the odd block
             eig(ii)=etemp2(jind)
             j1=0
             do i=1,nred
                if(matind(i,1).ne.matind(i,2))then
                   j1=j1+1
                   z(matind(i,1),ii)=z2(j1,jind)*recsqrtwo
                   z(matind(i,2),ii)=-z2(j1,jind)*recsqrtwo
                else
                   z(matind(i,1),ii)=dobzero
                endif
             enddo !i
             if(atoms%nlotot.gt.0) z(lapw%nv(jsp)+1:lapw%nv(jsp)+atoms%nlotot,ii)=dobzero
             do pos1=1,oddlocs
                do i=1,atoms%nlotot
#ifdef CPP_INVERSION
                   z(lapw%nv(jsp)+i,ii)=z(lapw%nv(jsp)+i,ii)+ z2(lapw%nv(jsp)-nred+pos1,jind)*locrec(i,oddlocindex(pos1))
#else
                   z(lapw%nv(jsp)+i,ii)=z(lapw%nv(jsp)+i,ii)+ z2(lapw%nv(jsp)-nred+pos1,jind)*conjg(locrec(i,oddlocindex(pos1)))
#endif
                enddo
             enddo
             jind=jind+1
          endif !evensort
       enddo !ii

#ifndef CPP_F90
       allocate(a(dimension%nbasfcn*(dimension%nbasfcn+1)/2))
       allocate(b(dimension%nbasfcn*(dimension%nbasfcn+1)/2))
#endif
    endif !sym%l_zref

    deallocate ( z1,z2,etemp1,etemp2,evensort )

  END SUBROUTINE zsymsecloc
      END MODULE
