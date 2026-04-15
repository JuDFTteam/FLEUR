      module m_wann_torque
            use m_judft
      private
      public :: wann_torque
      contains
      subroutine wann_torque(
     >               atoms,
     >               noccbd,
     >               acof,bcof,ccof,
     >               us,dus,duds,uds,
     >               ulos,dulos,
     >               atomlist_num,atomlist, 
     >               l_soc,nbnd,
     &               torque)
c******************************************************
c     Calculate matrix elements of the surface 
c     current operator.
c     Frank Freimuth
c******************************************************
      use m_types
      use m_constants, only:pimach

      implicit none

      TYPE(t_atoms), INTENT(IN) :: atoms

      integer, intent(in) :: noccbd

      complex, intent(in) :: acof(:,0:,:,:)
      complex, intent(in) :: bcof(:,0:,:,:)
      complex, intent(in) :: ccof(-atoms%llod:,:,:,:,:)

      real,    intent(in) :: us(0:,:,:)
      real,    intent(in) :: dus(0:,:,:)
      real,    intent(in) :: duds(0:,:,:)
      real,    intent(in) :: uds(0:,:,:)

      real,    intent(in) :: ulos(:,:,:)
      real,    intent(in) :: dulos(:,:,:)

      integer, intent(in) :: atomlist_num
      integer, intent(in) :: atomlist(:)

      logical, intent(in) :: l_soc
      integer, intent(in) :: nbnd

      complex, intent(out) :: torque(:,:,:,:)

      integer             :: at,n,itype,nt,minat,maxat
      integer             :: m,l,lm,j,i,jspin,lmd
      complex,allocatable :: surfcof(:,:,:,:)
      complex,allocatable :: surfder(:,:,:,:)
      complex,allocatable :: currmat(:,:,:,:)
      real                :: bmat(3,3)
      integer             :: tlp1,ilo,iat,jj
      integer(kind=8)     :: varint
      real                :: sqrtpi
      integer             :: lp,lpp,lppp,lmp,lmin,lmaxx
      integer             :: spin1,spin2,mp,mppp,lmppp,mpp,lmpp
      complex,parameter   :: ci=(0.0,1.0)
      integer             :: dir,list_ind
      logical             :: l_inthelist


      call timestart("wann_torque")
      lmd = (atoms%lmaxd+1)**2 - 1
      allocate ( surfcof(noccbd,0:lmd,atoms%nat,2) )
      allocate ( surfder(noccbd,0:lmd,atoms%nat,2) )  

      surfcof=0.0
      surfder=0.0
      do jspin=1,2
       nt=0
       DO itype=1,atoms%ntype
        minat=nt+1
        maxat=nt+atoms%neq(itype)
        do at=minat,maxat  
         do l=0,atoms%lmax(itype)  
          do m=-l,l  
           lm=l*(l+1)+m
           do j=1,noccbd
            surfcof(j,lm,at,jspin)=
     &        ( acof(j,lm,at,jspin)* us(l,itype,jspin)
     &                    +
     &          bcof(j,lm,at,jspin)* uds(l,itype,jspin) )*
     &              atoms%rmt(itype)*ci**l

            surfder(j,lm,at,jspin)=
     &        ( acof(j,lm,at,jspin)*dus(l,itype,jspin)
     &                    +
     &          bcof(j,lm,at,jspin)*duds(l,itype,jspin) )*
     &              atoms%rmt(itype)*ci**l

           enddo !j
          enddo !m 
         enddo !l 
         do ilo=1,atoms%nlo(itype)
          l=atoms%llo(ilo,itype)
          do m=-l,l
           lm=l*(l+1)+m
           do j=1,noccbd
c              if(.false.)then
            surfcof(j,lm,at,jspin)=surfcof(j,lm,at,jspin)+
     &        ( ccof(m,j,ilo,at,jspin)* ulos(ilo,itype,jspin) )*
     &              atoms%rmt(itype)*ci**l

            surfder(j,lm,at,jspin)=surfder(j,lm,at,jspin)+
     &        ( ccof(m,j,ilo,at,jspin)*dulos(ilo,itype,jspin) )*
     &              atoms%rmt(itype)*ci**l             
c            endif
           enddo !j
          enddo !m   
         enddo !ilo
        enddo !at
        nt=nt+atoms%neq(itype)
       enddo !itype
      enddo !jspin

      allocate( currmat(noccbd,noccbd,2,2) )
      at=0
      do itype=1,atoms%ntype
       do iat=1,atoms%neq(itype)
        at=at+1
        l_inthelist=.false.
        do list_ind=1,atomlist_num
          if(atomlist(list_ind).eq.at)then
             l_inthelist=.true.
             exit
          endif
        enddo !list_ind
        if(.not.l_inthelist)cycle
        currmat=cmplx(0.0,0.0)
        do spin2=1,2
         do spin1=1,2
          do l=0,atoms%lmax(itype)
           do m=-l,l
            lm=l*(l+1)+m
            do j=1,noccbd
             do jj=1,noccbd  
              currmat(jj,j,spin1,spin2)=
     &           currmat(jj,j,spin1,spin2)+
     &             conjg(surfcof(jj,lm,at,spin1))*
     &             surfder(j,lm,at,spin2)
             enddo !jj
            enddo !j
           enddo !m
          enddo !l
         enddo !spin1
        enddo !spin2

c-------x-component
c-------Pauli-x:
c-------( 0  1 )
c-------( 1  0 )
        if(l_soc)then
         do j=1,noccbd
          do jj=1,noccbd
           torque(1,jj,j,list_ind)=
     &            currmat(jj,j,2,1)+
     &            currmat(jj,j,1,2)
          enddo
         enddo
        else
         do j=1,noccbd
          do jj=1,noccbd
           torque(1,jj+nbnd,j,list_ind)=
     &            currmat(jj,j,2,1)
           torque(1,jj,j+nbnd,list_ind)=
     &            currmat(jj,j,1,2)
          enddo
         enddo
        endif

c-------y-component
c-------Pauli-y:
c-------( 0 -i )
c-------( i  0 )
        if(l_soc)then
         do j=1,noccbd
          do jj=1,noccbd
           torque(2,jj,j,list_ind)=
     &              ci*currmat(jj,j,2,1)-
     &              ci*currmat(jj,j,1,2)
          enddo
         enddo
        else
         do j=1,noccbd
          do jj=1,noccbd
           torque(2,jj+nbnd,j,list_ind)=
     &              ci*currmat(jj,j,2,1)
           torque(2,jj,j+nbnd,list_ind)=
     &             -ci*currmat(jj,j,1,2)

          enddo
         enddo
        endif

c-------z-component
c-------Pauli-z:
c-------( 1  0 )
c-------( 0 -1 )
        if(l_soc)then
         do j=1,noccbd
          do jj=1,noccbd
           torque(3,jj,j,list_ind)=
     &              currmat(jj,j,1,1)-
     &              currmat(jj,j,2,2)
          enddo
         enddo
        else
         do j=1,noccbd
          do jj=1,noccbd
           torque(3,jj,j,list_ind)=
     &              currmat(jj,j,1,1)

           torque(3,jj+nbnd,j+nbnd,list_ind)=
     &             -currmat(jj,j,2,2)
          enddo
         enddo
        endif
       enddo !iat
      enddo !itype

c      write(*,*)torque

c      stop 'ok till here'
      call timestop("wann_torque")
      end subroutine wann_torque
      
      end module m_wann_torque