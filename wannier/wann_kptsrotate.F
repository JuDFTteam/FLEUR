      module m_wann_kptsrotate
c****************************************
c     Rotate the wave function.
c     Frank Freimuth
c****************************************
      contains
      subroutine wann_kptsrotate(
     >               natd,nlod,llod,
     >               ntypd,nlo,llo,invsat,
     >               l_noco,l_soc,
     >               ntype,neq,nlotot,
     >               kveclo,jspin_in,
     >               oper,nop,mrot,nvd,
     >               nv,shiftkpt,
     >               tau,
     x               bkpt,k1,k2,k3,
     x               z,nsfactor)
      use m_inv3
      use m_constants,only:pimach
      implicit none
      integer,intent(in)    :: natd
      integer,intent(in)    :: nlod
      integer,intent(in)    :: llod
      integer,intent(in)    :: ntypd
      INTEGER, INTENT (IN)  :: nlo(ntypd)
      integer, intent(in)   :: llo(nlod,ntypd)
      INTEGER, INTENT (IN)  :: invsat(natd)
      logical,intent(in)    :: l_noco
      logical,intent(in)    :: l_soc
      integer,intent(in)    :: ntype
      integer,intent(in)    :: neq(ntype)
      integer,intent(in)    :: nlotot
      integer,intent(in)    :: kveclo(nlotot)
      integer,intent(in)    :: jspin_in
      integer,intent(in)    :: oper
      integer,intent(in)    :: nop
      integer,intent(in)    :: mrot(3,3,nop)
      integer,intent(in)    :: nvd,nv(:)
      integer,intent(in)    :: shiftkpt(3)
      real,intent(in)       :: tau(3,nop)
      real,intent(inout)    :: bkpt(3)
      integer,intent(inout) :: k1(:,:),k2(:,:),k3(:,:) !nvd,jspd
#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      complex,intent(inout) :: z(:,:)
#else
      real,intent(inout)    :: z(:,:)
#endif
      complex,intent(out)   :: nsfactor !phase of non-symmorphic ops

      real    :: bkrot(3),tpi,arg
      integer :: rotmat(3,3),determ,ilo
      integer :: j1,j2,j3,k,j,at_ind,invoper
      real    :: shiftnonsymm(3)
      real    :: phase
      INTEGER :: kvec(2*(2*llod+1),nlod,natd)
      integer :: natom,ntyp,lm,m,l,lo,nn,jspin
      integer :: absoper,jj,jsp_start,jsp_end
      integer :: testmat(3,3)

      tpi=2.0*pimach()

      absoper=abs(oper)

      call inv3(mrot(:,:,absoper),rotmat,determ)
      shiftnonsymm(:)=matmul(rotmat,tau(:,absoper))

c      testmat=matmul(mrot(:,:,absoper),rotmat)
c      print*,testmat
c      testmat=matmul(rotmat,mrot(:,:,absoper))
c      print*,testmat


#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      if(oper.lt.0)then
         z=conjg(z) ! time-inversion symmetry
         shiftnonsymm=-1.0*shiftnonsymm
      endif
#endif

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
              bkrot(:)=bkrot(:)+mrot(k,:,absoper)*bkpt(k)
            enddo
            bkpt(:)=bkrot(:)
          endif
          do j=1,nv(jspin)  !rotate reciprocal vector
               j1=mrot(1,1,absoper)*k1(j,jspin)+
     +            mrot(2,1,absoper)*k2(j,jspin)+
     +            mrot(3,1,absoper)*k3(j,jspin)
               j2=mrot(1,2,absoper)*k1(j,jspin)+
     +            mrot(2,2,absoper)*k2(j,jspin)+
     +            mrot(3,2,absoper)*k3(j,jspin)
               j3=mrot(1,3,absoper)*k1(j,jspin)+
     +            mrot(2,3,absoper)*k2(j,jspin)+
     +            mrot(3,3,absoper)*k3(j,jspin)
               k1(j,jspin)=j1
               k2(j,jspin)=j2
               k3(j,jspin)=j3
          enddo 
      enddo !jspin  

      if(oper.lt.0)then !time-inversion symmetry
         k1   = -k1
         k2   = -k2
         k3   = -k3
         bkpt = -bkpt
      endif

      do jspin=jsp_start,jsp_end
        jj=0 
        if(l_noco.and.(jspin.eq.2))then
           jj=nv(1)+nlotot
        endif
        do j=1,nv(jspin)
         phase = k1(j,jspin) * shiftnonsymm(1)
     +         + k2(j,jspin) * shiftnonsymm(2)
     +         + k3(j,jspin) * shiftnonsymm(3)
         phase = tpi*phase
         phase = cos(phase)
         z(j+jj,:)  = phase * z(j+jj,:) 
        enddo    
        jj=jj+nv(jspin)
        do ilo=1,nlotot
         j=kveclo(ilo)
         phase = k1(j,jspin) * shiftnonsymm(1)
     +         + k2(j,jspin) * shiftnonsymm(2)
     +         + k3(j,jspin) * shiftnonsymm(3)
         phase = tpi*phase
         phase = cos(phase)
         z(jj+ilo,:)  = phase * z(jj+ilo,:) 
        enddo
      enddo  


      bkpt(:)   = bkpt(:)   - shiftkpt(:)
      k1(:,:)   = k1(:,:)   + shiftkpt(1)
      k2(:,:)   = k2(:,:)   + shiftkpt(2)
      k3(:,:)   = k3(:,:)   + shiftkpt(3)

      write(6,*)"in wann_kptsrotate:"
      write(6,*) "bkpt=",bkpt

      arg = tpi*(
     +        bkpt(1)*shiftnonsymm(1)+
     +        bkpt(2)*shiftnonsymm(2)+
     +        bkpt(3)*shiftnonsymm(3)  )
      
      nsfactor = cmplx(cos(arg),sin(arg))




      end subroutine wann_kptsrotate
      end module m_wann_kptsrotate
