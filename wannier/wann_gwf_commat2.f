c*************************c
c  routine to set up the  c
c  composite matrices as  c
c  input for wannier90    c
c*************************c
      module m_wann_gwf_commat2
      USE m_juDFT
      implicit none
      contains

      subroutine wann_gwf_commat2(nkpts,nntot,bpt,nqpts,nntot_q,
     >                           bpt_q,gb,gb_q,latt_const_q,
     >                           l_unformatted,l_amn,l_mmn)
      use m_wann_gwf_tools
      use m_constants, only : pimach

      implicit none
      
      ! input parameters
      logical,intent(in) :: l_unformatted,l_amn,l_mmn
      integer,intent(in) :: nkpts,nntot,nqpts,nntot_q
      integer,intent(in) :: bpt(nntot,nkpts),bpt_q(nntot_q,nqpts)
      integer,intent(in) :: gb(3,nntot,nkpts),gb_q(3,nntot_q,nqpts)
      real,   intent(in) :: latt_const_q

      ! allocatable arrays
      integer,allocatable :: gb_kq(:,:,:)
      integer,allocatable :: bpt_kq(:,:)
      integer,allocatable :: g1(:,:,:),g2(:,:)
      integer,allocatable :: b1(:,:),b2(:)
      complex,allocatable :: tmp_a(:,:,:),tmp_m(:,:,:,:)
      complex,allocatable :: mmnk(:,:,:,:),mmnq(:,:,:,:)
      complex,allocatable :: mk(:,:,:),mq(:,:,:)
      complex,allocatable :: mmn(:,:,:,:)
      character(len=20),allocatable :: feig(:),famn(:)
      character(len=20),allocatable :: fmmn(:),fmmn2(:)

      ! further variables
      integer :: nbnd,nwfs,kb,qb,nn_kq,kqb
      integer :: i,j,ikqpt,ikqpt_b,nk,nq,g(3)
      integer :: nbnd_arti,nqpts_arti,nntot_arti
      integer :: nkqpts,nntot_kq,nn,ikqpt_help
      integer :: dump,nn_arti,nkp,k,b,q,nwf,n,is,ie
      real :: mmn_r,mmn_i,tau,tpi,a_arti,b_arti,eig
      complex :: mmn_arti,amn_arti,amn
      logical :: l_exist,l_miss,l_fe,l_fa,l_fm,l_fm2,l_proj
      character(len=20) :: fq

      integer :: ncpu

      call timestart("wann_gwf_commat2")

      ncpu = 8

      write(*,*)'create HDWFs input for w90...'

      nkqpts=nkpts*nqpts
      tpi = 2.0*pimach()
      a_arti = latt_const_q
      b_arti = 0.98*latt_const_q
      l_fe   = .true.
      l_fa   = .true.
      l_fm   = .true.
      l_fm2  = .true.
      l_exist= .true.

! compute auxiliary overlaps and projections
      tau = tpi/real(nqpts)/a_arti 
      mmn_arti = 2.0*tpi*tpi*sin(tau*b_arti/2.0)
      mmn_arti = mmn_arti/(tpi*tpi-tau*tau*b_arti*b_arti)
      mmn_arti = mmn_arti/(tau*b_arti)

      amn_arti = cmplx(1.0,0.0)

      open(777,file='ma_aux')
      write(777,'(2(f24.18))')real(mmn_arti),aimag(mmn_arti)
      write(777,'(2(f24.18))')real(amn_arti),aimag(amn_arti)
      close(777)

! are all necessary files present?
      allocate(feig(nqpts),famn(nqpts),fmmn(nqpts),fmmn2(nqpts))
      l_miss=.false.
      do q=1,nqpts
         write(fq,'(i4.4)')q
         feig(q) = 'WF1_'//trim(fq)//'.eig'    ! energies
         famn(q) = 'WF1_'//trim(fq)//'.amn'    ! projections
         fmmn(q) = 'WF1_'//trim(fq)//'.mmn'    ! k-overlaps
         fmmn2(q)= 'param_'//trim(fq)//'.mmn'  ! q-overlaps

         inquire(file=feig(q) ,exist=l_fe )
         if(l_amn) inquire(file=famn(q) ,exist=l_fa )
         if(l_mmn) inquire(file=fmmn(q) ,exist=l_fm )
         if(l_mmn) inquire(file=fmmn2(q),exist=l_fm2)

         if(.not.l_fe ) write(*,*)'missing: ',feig(q)
         if(.not.l_fa ) write(*,*)'missing: ',famn(q)
         if(.not.l_fm ) write(*,*)'missing: ',fmmn(q)
         if(.not.l_fm2) write(*,*)'missing: ',fmmn2(q)
         if(.not.(l_fe.and.l_fa.and.l_fm.and.l_fm2))l_miss=.true.
      enddo
      if(l_mmn) inquire(file='bkqpts',exist=l_exist)
      if(.not.l_exist) write(*,*)'missing: bkqpts'
      inquire(file='proj',exist=l_proj)
      if(.not.l_proj) write(*,*)'missing: proj'
      if(l_miss.or.(.not.l_exist).or.(.not.l_proj))
     >   call juDFT_error("missing file(s) for HDWFs")

! get number of bands and wfs from proj
      open(405,file='proj',status='old')
      read(405,*)nwfs,nbnd
      close(405)
      write(*,*)'nbnd=',nbnd
      write(*,*)'nwfs=',nwfs


c*****************************c
c     COMPOSITE .EIG FILE     c
c*****************************c
      open(305,file='WF1_gwf.eig')
      do q=1,nqpts
         open(405,file=feig(q))
         do k=1,nkpts
            ikqpt=get_index_kq(k,q,nkpts)
            do i=1,nbnd
               read(405,*)b,nk,eig
               write(305,'(2i12,f19.13)')b,ikqpt,eig
            enddo
         enddo
         close(405)
      enddo
      close(305)


c*****************************c
c     COMPOSITE .AMN FILE     c
c*****************************c
      if(.not.l_amn) goto 100 ! skip amn part

      if(l_unformatted) then
         allocate(tmp_a(nbnd,nwfs,nkqpts))
         do q=1,nqpts
            is = get_index_kq(1,q,nkpts)
            ie = get_index_kq(nkpts,q,nkpts)
            open(405,file=famn(q),form='unformatted')
            read(405)b,k,nwf
            read(405)tmp_a(:,:,is:ie)
            close(405)
         enddo

         open(306,file='WF1_gwf.amn',form='unformatted')
         write(306)nbnd,nkqpts,nwfs
         write(306)tmp_a*amn_arti
         close(306)
         deallocate(tmp_a)
      else
         open(306,file='WF1_gwf.amn')
         write(306,*)'Projections for HDWFs'
         write(306,'(i5,i7,i5)')nbnd,nkqpts,nwfs
         
         do q=1,nqpts
           open(405,file=famn(q))
           read(405,*)!title
           read(405,*)b,k,nwf
           do k=1,nkpts
            ikqpt=get_index_kq(k,q,nkpts)
            do nwf=1,nwfs
             do i=1,nbnd
              read(405,*)b,n,nk,mmn_r,mmn_i
              amn = cmplx(mmn_r,mmn_i)
              amn = amn*amn_arti
              write(306,'(i5,i5,i7,3x,2f18.12)')
     >                     b,n,ikqpt,real(amn),aimag(amn)
             enddo
            enddo
           enddo
         enddo
         close(306)
      endif!l_unformatted

 100  continue


c*****************************c
c     COMPOSITE .MMN FILE     c
c*****************************c
      if(.not.l_mmn) goto 200 ! skip mmn part

      open (202,file='bkqpts',form='formatted',status='old')
      rewind (202)
      read (202,'(i4)') nntot_kq
      write (*,*) 'nntot_kq=',nntot_kq
      allocate ( gb_kq(4,nntot_kq,nkqpts),
     &           bpt_kq(nntot_kq,nkqpts))
      do ikqpt=1,nkqpts
       do nn=1,nntot_kq
         read (202,'(2i6,3x,4i4)')
     &    ikqpt_help,bpt_kq(nn,ikqpt),(gb_kq(i,nn,ikqpt),i=1,4)
       enddo
      enddo
      close (202)

      if(l_unformatted) then
         allocate(mmnk(nbnd,nbnd,nntot,nkpts))
         allocate(mmnq(nbnd,nbnd,nntot_q,nkpts))
         allocate(g1(3,nntot,nkpts),g2(3,nntot_q))
         allocate(b1(nntot,nkpts),b2(nntot_q))
         allocate(mmn(nbnd,nbnd,nntot_kq,nkqpts))

         do q=1,nqpts
           is = get_index_kq(1,q,nkpts)
           ie = get_index_kq(nkpts,q,nkpts)
         
           open(405,file=fmmn(q),form='unformatted')
           read(405)b,k,nn
           read(405)b1,g1
           read(405)mmnk
           close(405)

           open(405,file=fmmn2(q),form='unformatted')
           read(405)b,k,nn
           read(405)b2,g2
           read(405)mmnq
           close(405)
           mmnq = mmnq*conjg(mmn_arti)

           do k=1,nkpts
             ikqpt=get_index_kq(k,q,nkpts)           
             do nn=1,nntot_kq
               kqb = bpt_kq(nn,ikqpt)
               kb=get_index_k(kqb,nkpts)
               qb=get_index_q(kqb,nkpts)
               if(q.eq.qb) then ! k-overlap
                  nk=get_index_nn_k(bpt(:,k),nntot,kb,
     >                         gb_kq(:,nn,ikqpt),
     >                         gb(1:3,1:nntot,k))
                  mmn(:,:,nn,ikqpt)=mmnk(:,:,nk,k)              
               elseif(k.eq.kb) then ! q-overlap
                  nq=get_index_nn_q(bpt_q(:,q),nntot_q,qb,
     >                         gb_kq(:,nn,ikqpt),
     >                         gb_q(1:3,1:nntot_q,q))
                  mmn(:,:,nn,ikqpt)=mmnq(:,:,nq,k)
               else ! otherwise problem
                  call juDFT_error("problem in overlap mmn gwf",
     >                        calledby="wann_gwf_commat")
               endif
             enddo!nn
           enddo!k
         enddo!q
         deallocate(b1,b2,g1,g2,mmnk,mmnq)

         open(307,file='WF1_gwf.mmn',form='unformatted')
         write(307)nbnd,nkqpts,nntot_kq
         write(307)bpt_kq,gb_kq
         write(307)mmn
         close(307)
         deallocate(mmn)
      else
         allocate(mk(nbnd,nbnd,nntot))
         allocate(mq(nbnd,nbnd,nntot_q))
         open(307,file='WF1_gwf.mmn')
         write(307,*)'Overlaps for HDWFs'
         write(307,'(i5,i7,i5)')nbnd,nkqpts,nntot_kq

         do q=1,nqpts       
           open(405,file=fmmn(q))
           read(405,*)!title
           read(405,*)b,k,nn
           open(406,file=fmmn2(q))
           read(406,*)!title
           read(406,*)b,k,nn
           do k=1,nkpts
             ikqpt=get_index_kq(k,q,nkpts)           

             do nn=1,nntot ! read k-overlaps
               kb = bpt(nn,k)
               read(405,*)nkp,b,g(1:3)
               do i=1,nbnd
                  do j=1,nbnd
                     read(405,*)mmn_r,mmn_i
                     mk(j,i,nn)=cmplx(mmn_r,mmn_i)
                  enddo
               enddo               
             enddo  

             do nn=1,nntot_q ! read q-overlaps
               qb = bpt_q(nn,q)
               read(406,*)nkp,b,g(1:3)
               do i=1,nbnd
                  do j=1,nbnd
                     read(406,*)mmn_r,mmn_i
                     mq(j,i,nn)=cmplx(mmn_r,mmn_i)
                  enddo
               enddo                              
             enddo
             mq = mq*conjg(mmn_arti)   

             do nn=1,nntot_kq ! write composite overlaps
               kqb = bpt_kq(nn,ikqpt)
               kb=get_index_k(kqb,nkpts)
               qb=get_index_q(kqb,nkpts)
               write(307,'(i7,i7,3x,4i4)')ikqpt,kqb,gb_kq(1:4,nn,ikqpt)
               if(q.eq.qb) then ! k-overlap
                  nk=get_index_nn_k(bpt(:,k),nntot,kb,
     >                         gb_kq(:,nn,ikqpt),
     >                         gb(1:3,1:nntot,k))
                  do i=1,nbnd
                     do j=1,nbnd
                        write(307,'(2f24.18)')
     >                      real(mk(j,i,nk)),aimag(mk(j,i,nk))
                     enddo
                  enddo
               elseif(k.eq.kb) then ! q-overlap
                  nq=get_index_nn_q(bpt_q(:,q),nntot_q,qb,
     >                         gb_kq(:,nn,ikqpt),
     >                         gb_q(1:3,1:nntot_q,q))
                  do i=1,nbnd
                     do j=1,nbnd
                        write(307,'(2f24.18)')
     >                      real(mq(j,i,nq)),aimag(mq(j,i,nq))
                     enddo
                  enddo
               else ! otherwise problem
                  call juDFT_error("problem in overlap mmn gwf",
     >                        calledby="wann_gwf_commat")
               endif
             enddo!nn
           enddo!k
           close(406)
           close(405)
         enddo!q           
         close(307)
         deallocate(mk,mq)
      endif!l_unformatted

      deallocate(gb_kq,bpt_kq)
 200  continue
      deallocate(feig,famn,fmmn,fmmn2)


      call timestop("wann_gwf_commat2")
      end subroutine wann_gwf_commat2  


      subroutine split_array(numpoints,n,ncpu,counts,displs)
      implicit none
      integer, intent(in) :: numpoints,ncpu,n
      integer, intent(out) :: counts,displs

      integer :: ratio, remainder, i

      call timestart("split_array")

      ratio = numpoints / ncpu
      remainder = mod(numpoints,ncpu)

      if(n < remainder) then
         counts = ratio+1
         displs = i*(ratio+1)
      else
         counts = ratio
         displs = remainder*(ratio+1)+(i-remainder)*ratio
      endif   
      
      call timestop("split_array")
      end subroutine split_array

      end module m_wann_gwf_commat2
