!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter GrÃ¼nberg Institut, Forschungszentrum JÃ¼lich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_kptsrotate
      use m_judft
c****************************************
c     Rotate the wave function.
c     Frank Freimuth
c****************************************
      contains
      subroutine wann_kptsrotate(
     >               atoms,
     >               invsat,
     >               l_noco,l_soc,
     >               jspin_in,
     >               oper,nop,mrot,nvd,
     >               shiftkpt,
     >               tau,
     >               lapw,
!     x               bkpt,k1,k2,k3,
     x               zMat,nsfactor)

      USE m_types
      use m_constants
      use m_inv3

      implicit none
      TYPE(t_atoms),INTENT(IN)             :: atoms
      INTEGER, INTENT (IN)  :: invsat(:)
      logical,intent(in)    :: l_noco
      logical,intent(in)    :: l_soc
      integer,intent(in)    :: jspin_in
      integer,intent(in)    :: oper
      integer,intent(in)    :: nop
      integer,intent(in)    :: mrot(3,3,nop)
      integer,intent(in)    :: nvd
      integer,intent(in)    :: shiftkpt(3)
      real,intent(in)       :: tau(3,nop)
!      real,intent(inout)    :: bkpt(3)
!      integer,intent(inout) :: k1(:,:),k2(:,:),k3(:,:) !nvd,jspd
      TYPE(t_lapw),intent(inout)  :: lapw

      
      TYPE(t_mat), INTENT (INOUT) :: zMat !z(nbasfcn,noccbd) !can be real/complex

      complex,intent(out)   :: nsfactor !phase of non-symmorphic ops

      real    :: bkrot(3),tpi,arg
      integer :: rotmat(3,3),determ,iatom,itype
      integer :: j1,j2,j3,k,j,at_ind,invoper,nkvec
      real    :: shiftnonsymm(3)
      real    :: phase
      integer :: natom,ntyp,lm,m,l,lo,nn,jspin
      integer :: absoper,jj,jsp_start,jsp_end
      integer :: testmat(3,3),nbasf

      call timestart("wann_kptsrotate")
      tpi=2.0*pimach()

      absoper=abs(oper)

      call inv3(mrot(:,:,absoper),rotmat,determ)
      shiftnonsymm(:)=matmul(rotmat,tau(:,absoper))

c      testmat=matmul(mrot(:,:,absoper),rotmat)
c      print*,testmat
c      testmat=matmul(rotmat,mrot(:,:,absoper))
c      print*,testmat

      IF(.NOT.zMat%l_real) THEN
         if(oper.lt.0)then
            zMat%data_c = CONJG(zMat%data_c)
            shiftnonsymm=-1.0*shiftnonsymm
         endif
      END IF

      if(l_noco) then
         jsp_start=1
         jsp_end=2
      else
         jsp_start=jspin_in
         jsp_end=jspin_in
      endif



      do jspin=jsp_start,jsp_end
          if(.not.(l_noco.and.(jspin.eq.2)))then
            bkrot(:)=0.0
            do k=1,3
              bkrot(:)=bkrot(:)+mrot(k,:,absoper)*lapw%bkpt(k)
            enddo
            lapw%bkpt(:)=bkrot(:)
          endif
          do j=1,lapw%nv(jspin)  !rotate reciprocal vector
               j1=mrot(1,1,absoper)*lapw%k1(j,jspin)+
     +            mrot(2,1,absoper)*lapw%k2(j,jspin)+
     +            mrot(3,1,absoper)*lapw%k3(j,jspin)
               j2=mrot(1,2,absoper)*lapw%k1(j,jspin)+
     +            mrot(2,2,absoper)*lapw%k2(j,jspin)+
     +            mrot(3,2,absoper)*lapw%k3(j,jspin)
               j3=mrot(1,3,absoper)*lapw%k1(j,jspin)+
     +            mrot(2,3,absoper)*lapw%k2(j,jspin)+
     +            mrot(3,3,absoper)*lapw%k3(j,jspin)
               lapw%k1(j,jspin)=j1
               lapw%k2(j,jspin)=j2
               lapw%k3(j,jspin)=j3
          enddo 
      enddo !jspin  

      if(oper.lt.0)then !time-inversion symmetry
         lapw%k1   = -lapw%k1
         lapw%k2   = -lapw%k2
         lapw%k3   = -lapw%k3
         lapw%bkpt = -lapw%bkpt
      endif

      do jspin=jsp_start,jsp_end
        jj=0 
        if(l_noco.and.(jspin.eq.2))then
           jj=lapw%nv(1)+atoms%nlotot
        endif
        do j=1,lapw%nv(jspin)
         phase = lapw%k1(j,jspin) * shiftnonsymm(1)
     +         + lapw%k2(j,jspin) * shiftnonsymm(2)
     +         + lapw%k3(j,jspin) * shiftnonsymm(3)
         phase = tpi*phase
         phase = cos(phase)
         IF(zMat%l_real) THEN
            zMat%data_r(j+jj,:)  = phase * zMat%data_r(j+jj,:) 
         ELSE
            zMat%data_c(j+jj,:)  = phase * zMat%data_c(j+jj,:) 
         END IF
        enddo    
        jj=jj+lapw%nv(jspin)
c$$$        do ilo=1,nlotot
c$$$           call judft_error("BUG in wann_kptsrotate:LOs missing")
c$$$         !j=kveclo(ilo)
c$$$         phase = lapw%k1(j,jspin) * shiftnonsymm(1)
c$$$     +         + lapw%k2(j,jspin) * shiftnonsymm(2)
c$$$     +         + lapw%k3(j,jspin) * shiftnonsymm(3)
c$$$         phase = tpi*phase
c$$$         phase = cos(phase)
c$$$         IF(zMat%l_real) THEN
c$$$            zMat%data_r(jj+ilo,:)  = phase * zMat%data_r(jj+ilo,:) 
c$$$         ELSE
c$$$            zMat%data_c(jj+ilo,:)  = phase * zMat%data_c(jj+ilo,:) 
c$$$         END IF
c$$$  enddo
        do iatom=1,atoms%nat
          iType = atoms%itype(iAtom)
          do lo=1,atoms%nlo(itype)
             DO nkvec = 1, lapw%nkvec(lo,iAtom)
!                   iLAPW = lapw%kvec(nkvec,lo,iAtom)
                   nbasf=lapw%nv(jspin)+lapw%index_lo(lo,iatom)+nkvec
                   IF(zMat%l_real) THEN
            zMat%data_r(nbasf,:)  = phase * zMat%data_r(nbasf,:) 
                   ELSE
            zMat%data_c(nbasf,:)  = phase * zMat%data_c(nbasf,:) 
                   END IF                   
             END DO
          enddo !lo   
        enddo !iatom
        
      enddo  


      lapw%bkpt(:)   = lapw%bkpt(:)   - shiftkpt(:)
      lapw%k1(:,:)   = lapw%k1(:,:)   + shiftkpt(1)
      lapw%k2(:,:)   = lapw%k2(:,:)   + shiftkpt(2)
      lapw%k3(:,:)   = lapw%k3(:,:)   + shiftkpt(3)

      write(oUnit,*)"in wann_kptsrotate:"
      write(oUnit,*) "bkpt=",lapw%bkpt

      arg = tpi*(
     +        lapw%bkpt(1)*shiftnonsymm(1)+
     +        lapw%bkpt(2)*shiftnonsymm(2)+
     +        lapw%bkpt(3)*shiftnonsymm(3)  )
      
      nsfactor = cmplx(cos(arg),sin(arg))
      call timestop("wann_kptsrotate")

      end subroutine wann_kptsrotate
      end module m_wann_kptsrotate
