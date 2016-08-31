!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_projmethod
      use m_juDFT
      contains
      subroutine wann_projmethod(
     >               fullnkpts,
     >               l_projmethod,l_bestproj,
     >               ikpt,nwfs,nslibd,amn,eig,
     <               psiw,hwfr)
c**********************************************************************
c     Construct Wannier functions from the projections of the
c     Bloch functions onto the trial orbitals. These projections
c     are provided as input to this subroutine in the amn-matrix.
c     The Wannier functions are orthonormalized. The orthonormalization
c     requires the construction of square root of the overlap 
c     matrix (omn). Two different algorithms for the calculation of the
c     square root are available (projmethod/bestproj).
c
c     YM && FF
c**********************************************************************
      use m_types
      implicit none
      integer, intent(in)    :: fullnkpts
      logical, intent(in)    :: l_projmethod
      logical, intent(in)    :: l_bestproj
      integer, intent(in)    :: ikpt
      integer, intent(in)    :: nwfs
      integer, intent(in)    :: nslibd
      complex, intent(in)    :: amn(:,:,:)
      real,    intent(in)    :: eig(:)

      complex, intent(inout) :: psiw(:,:,:)
      complex, intent(inout) :: hwfr(:,:)

      real                   :: ei(nwfs)
      complex, allocatable   :: work(:)
      integer, allocatable   :: iwork(:)
      real, allocatable      :: rwork(:)
      integer                :: info,lwork,lrwork,liwork,lee
      integer                :: nwf,nwfp,i,n,lda,j
      complex                :: omn(nwfs,nwfs),vec(nwfs,nwfs)
      complex                :: smn(nwfs,nwfs)
      complex, allocatable   :: amn_copy(:,:)
      character              :: jobu*1,jobv*1
      character              :: jobz*1,uplo*1
      complex, allocatable   :: sunit(:,:)
      real,    allocatable   :: sdiag(:)
      complex, allocatable   :: umn(:,:),vmn(:,:)
      integer                :: ldu,ldv,m
      complex                :: hwfk(nwfs,nwfs)

      if(l_projmethod) then
         omn=cmplx(0.0,0.0)
c..      now we determine the Omn matrix, PRB 71,125119, Eq.(13)
         do nwf = 1,nwfs
          do nwfp = 1,nwfs
           do i = 1,nslibd
              omn(nwf,nwfp) = omn(nwf,nwfp) +
     +               conjg(amn(i,nwf,ikpt))*amn(i,nwfp,ikpt)
           enddo
          enddo
         enddo

c        lwork = 2*nwfs + nwfs*nwfs
c        lrwork = 1 + 4*nwfs + 5*nwfs*nwfs
c        liwork = 2 + 6*nwfs
         lee = log( dble(nwfs) )/log(2.d0) + 1
         lwork = 1 + 5*nwfs + 2*nwfs*lee + 3*(nwfs**2)
         lrwork = 1 + 4*nwfs + 2*nwfs*lee + 3*(nwfs**2)
         liwork = 2 + 5*nwfs +1

         allocate (work(lwork),rwork(lrwork),iwork(liwork))

         jobz = 'V' ; uplo = 'L' ; n = nwfs ; lda = nwfs

         call zheevd(jobz,uplo,n,vec,lda,ei,work,lwork,
     &             rwork,lrwork,iwork,liwork,info)

         if (info.lt.0) write (6,*)
     &               'ith argument had an illegal value ',info
         IF (info>0)  CALL juDFT_error("not converged diagonalization"
     +        ,calledby ="wann_projmethod")

         do i = 1,nwfs
          if (ei(i).le.(-1.e-6)) then
           print*,'eigenvalue is less than zero'
          elseif (abs(ei(i)).le.1.e-16) then
           print*,'eigenvalue is very close to zero:',ei(i)
          endif
         enddo

c..      constructing the square root of the O matrix, Eq.(14)

         smn(:,:) = cmplx(0.,0.)

         do i = 1,nwfs
          do j = 1,nwfs
           do n = 1,nwfs
            smn(i,j) = smn(i,j) +
     +         conjg(vec(i,n))*(1./sqrt(ei(n)))*vec(j,n)
           enddo
          enddo
         enddo

         deallocate (work,rwork,iwork)
c..      now constructing the overlaps of the wavefunctions
c..      with the wannier functions in reciprocal space, Eq.(16)

         psiw(:,:,ikpt) = cmplx(0.,0.)
         do i = 1,nslibd
          do j = 1,nwfs
           do n = 1,nwfs
               psiw(i,j,ikpt) = psiw(i,j,ikpt) +
     +                       smn(j,n)*amn(i,n,ikpt)
           enddo
          enddo
         enddo
      endif !wann%l_projmethod

      if(l_bestproj) then

c..      psiw matrix has the meaning of the Umn rotational matrix
c..      applied to the wavefunctions at each point in the BZone
c..      here it is calculated differently,
c..      according to the PRB 65,035109 (Eq.23)
c..
c..      the singe value decomposition of the Amn matrix is found:
c..                      A = U*SS*V
c..      Then the psiw matrix, psiw = Amn*Smn is given by:
c..                     psiw = U*1*V

         allocate (amn_copy(nslibd,nwfs),sunit(nslibd,nwfs))
         allocate (umn(nslibd,nslibd),vmn(nwfs,nwfs))
         allocate (sdiag(max(nslibd,nwfs)))

         lwork = max( 3*min(nslibd,nwfs) + max(nslibd,nwfs),
     &                        5*min(nslibd,nwfs) )

         allocate ( work(lwork),rwork(max(1,5*min(nslibd,nwfs))) )
!
!        Compute the eigenvalues and eigenvectors.
!
         jobu = 'A' ; jobv = 'A'
         lda = nslibd ; ldu = nslibd ; ldv = nwfs
!
!        The input matrix is destroyed by the routine.  Since we need to keep
!        it around, we only pass a copy to the routine.
!
         amn_copy(1:nslibd,1:nwfs) = amn(1:nslibd,1:nwfs,ikpt)

         call zgesvd(jobu,jobv,nslibd,nwfs,amn_copy,lda,
     >             sdiag,umn,ldu,vmn,ldv,work,lwork,rwork,info)

         if (info /= 0) then
            write (*,*) 'LAPACK routine ZGESVD returned a nonzero'
            write (*,*) 'value of the error flag, INFO =',info
         end if
!
!        Make the MxN identical Sunit matrix
!
         sunit(1:nslibd,1:nwfs) = cmplx(0.,0.)
         do i = 1,min(nslibd,nwfs)
            sunit(i,i) = cmplx(1.,0.)
         end do


!        and finally the psiw matrix

         psiw(:,:,ikpt) = cmplx(0.,0.)
         do i = 1,nslibd
          do j = 1,nwfs
           do n = 1,nslibd
            do m = 1,nwfs
               psiw(i,j,ikpt) = psiw(i,j,ikpt) +
     +                     umn(i,n)*sunit(n,m)*vmn(m,j)
            enddo
           enddo
          enddo
         enddo

         deallocate (work,rwork,amn_copy,sunit,vmn,umn,sdiag)

         write (*,*) 'The Psiw matrix was calculated via SVD'

      endif   !wann%l_bestproj



c..constructing the WF-hamiltonian in reciprocal space Eq.(23)

      hwfk(:,:) = cmplx(0.,0.)

      do i = 1,nwfs
         do j = 1,nwfs
            do n = 1,nslibd
               hwfk(i,j) = hwfk(i,j) +
     +            eig(n)*psiw(n,i,ikpt)*conjg(psiw(n,j,ikpt))

            enddo
         enddo
      enddo

c..   adding up the k-point reciprocal hamiltonians Eq.(20)
c..   to get the hopping elements inside the same unit cell

      hwfr=hwfr+hwfk/fullnkpts


c..   now we diagonalize the reciprocal hamiltonian for the
c..   purpose of testing

      do nwf = 1,nwfs
       do nwfp = 1,nwfs
        vec(nwf,nwfp) = hwfk(nwf,nwfp)
       enddo
      enddo

      lee = log( dble(nwfs) )/log(2.d0) + 1
      lwork = 1 + 5*nwfs + 2*nwfs*lee + 3*(nwfs**2)
      lrwork = 1 + 4*nwfs + 2*nwfs*lee + 3*(nwfs**2)
      liwork = 2 + 5*nwfs +1

      allocate (work(lwork),rwork(lrwork),iwork(liwork))

      jobz = 'V' ; uplo = 'L' ; n = nwfs ; lda = nwfs

      call zheevd(jobz,uplo,n,vec,lda,ei,work,lwork,
     &             rwork,lrwork,iwork,liwork,info)

      if (info.lt.0) write (6,*)
     &               'ith argument had an illegal value ',info
      IF (info>0)  CALL juDFT_error("not converged diagonalization"
     +     ,calledby ="wann_projmethod")

      do i = 1,nwfs
         write(6,*) ei(i)
      enddo

      deallocate (work,rwork,iwork)


      end subroutine wann_projmethod
      end module m_wann_projmethod
