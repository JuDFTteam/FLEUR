      module m_wann_uHu_commat
      USE m_juDFT
      implicit none
      contains

      subroutine wann_uHu_commat(nkpts,nntot,bpt,nqpts,nntot_q,
     >                           bpt_q,gb,gb_q,latt_const_q,
     >                           l_unformatted,l_dim,nparampts,
     >                           param_vec)
      use m_wann_gwf_tools
      use m_wann_gwf_auxovlp
      use m_constants, only : pimach

      implicit none
      
      ! input parameters
      logical,intent(in) :: l_unformatted,l_dim(3)
      integer,intent(in) :: nkpts,nntot,nqpts,nntot_q,nparampts
      integer,intent(in) :: bpt(nntot,nkpts),bpt_q(nntot_q,nqpts)
      integer,intent(in) :: gb(3,nntot,nkpts),gb_q(3,nntot_q,nqpts)
      real,   intent(in) :: latt_const_q,param_vec(3,nparampts)

      ! allocatable arrays
      integer,allocatable :: gb_kq(:,:,:)
      integer,allocatable :: bpt_kq(:,:)
      integer,allocatable :: g1(:,:,:),g2(:,:)
      integer,allocatable :: b1(:,:),b2(:)
      complex,allocatable :: tmp_a(:,:,:),tmp_m(:,:,:,:)
      complex,allocatable :: uHuk(:,:,:,:,:)
      complex,allocatable :: mk(:,:,:,:)
      complex,allocatable :: uHu(:,:,:,:,:)
      character(len=20),allocatable :: fuHu(:),fuHu2(:)

      ! further variables
      integer :: nbnd,nwfs,kb,qb,kb2,qb2,nn_kq,kqb,kqb2,arr_len,shift(3)
      integer :: i,j,ikqpt,ikqpt_b,nk,nq,nk2,nq2,g(3)
      integer :: nbnd_arti,nqpts_arti,nntot_arti
      integer :: nkqpts,nntot_kq,nn,nn2,ikqpt_help
      integer :: dump,nn_arti,nkp,k,b,q,nwf,n,is,ie
      real :: mmn_r,mmn_i,tau,tpi,a_arti,b_arti,eig
      complex :: ovaux
      logical :: l_exist,l_miss,l_fe,l_fa,l_fm,l_fm2,l_proj
      character(len=20) :: fq
      character(len=32) :: fmt,fmt2

      write(*,*)'create uHu for HDWF case...'

      call get_dimension(l_dim,arr_len)
      call get_shift(l_dim,shift)
      
      write(*,*)'dimension:',arr_len
      write(*,*)'x?,y?,z? :',l_dim(1:3)
      if(arr_len.le.3) call juDFT_error("dimension<4",
     >                                  calledby='wann_gwf_commat')

      nkqpts=nkpts*nqpts
      tpi = 2.0*pimach()
      l_fm   = .true.
      l_fm2  = .true.
      l_exist= .true.

! are all necessary files present?
      allocate(fuHu(nqpts),fuHu2(nqpts))
      l_miss=.false.
      do q=1,nqpts
         write(fq,'(i4.4)')q
         fuHu(q) = 'WF1_'//trim(fq)//'.uHu'    ! k-overlaps
         fuHu2(q)= 'param_'//trim(fq)//'.uHu'  ! q-overlaps

         inquire(file=fuHu(q) ,exist=l_fm )
         l_fm2 = .true.
         !inquire(file=fuHu2(q),exist=l_fm2)

         if(.not.l_fm ) write(*,*)'missing: ',fuHu(q)
         if(.not.l_fm2) write(*,*)'missing: ',fuHu2(q)
         if(.not.(l_fm.and.l_fm2))l_miss=.true.
      enddo
      inquire(file='bkqpts',exist=l_exist)
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
c     COMPOSITE .MMN FILE     c
c*****************************c

      write(fmt,'(a,i1,a)')'(2i6,3x,',arr_len,'i4)'
      write(fmt2,'(a,i1,a)')'(i7,i7,3x,',arr_len,'i4)'

      open (202,file='bkqpts',form='formatted',status='old')
      rewind (202)
      read (202,'(i4)') nntot_kq
      write (*,*) 'nntot_kq=',nntot_kq
      allocate ( gb_kq(arr_len,nntot_kq,nkqpts),
     &           bpt_kq(nntot_kq,nkqpts))
      do ikqpt=1,nkqpts
       do nn=1,nntot_kq
          read (202,fmt)
     &      ikqpt_help,bpt_kq(nn,ikqpt),(gb_kq(i,nn,ikqpt),i=1,arr_len)
       enddo
      enddo
      close (202)

      if(l_unformatted) then
         allocate(uHuk(nbnd,nbnd,nntot,nntot,nkpts))
         allocate(g1(3,nntot,nkpts),g2(3,nntot_q))
         allocate(b1(nntot,nkpts),b2(nntot_q))
         allocate(uHu(nbnd,nbnd,nntot_kq,nntot_kq,nkqpts))

         do q=1,nqpts
           is = get_index_kq(1,q,nkpts)
           ie = get_index_kq(nkpts,q,nkpts)
         
           open(405,file=fuHu(q),form='unformatted')
           read(405)b,k,nn
           read(405)b1,g1
           read(405)uHuk
           close(405)
c           open(405,file=fuHu2(q),form='unformatted')
c           read(405)b,k,nn
c           read(405)b2,g2
c           read(405)mmnq
c           close(405)
c           mmnq = mmnq*conjg(mmn_arti)

           do k=1,nkpts
             ikqpt=get_index_kq(k,q,nkpts)           
             do nn=1,nntot_kq
               kqb = bpt_kq(nn,ikqpt)
               kb=get_index_k(kqb,nkpts)
               qb=get_index_q(kqb,nkpts)

               do nn2=1,nntot_kq
               kqb2= bpt_kq(nn2,ikqpt)
               kb2=get_index_k(kqb2,nkpts)
               qb2=get_index_q(kqb2,nkpts)

               CALL wann_gwf_auxovlp(param_vec(:,qb),param_vec(:,qb2),
     >                               latt_const_q,ovaux)
               if(ovaux.ne.1.0) write(*,*)'ovaux',ovaux

               if((q.eq.qb) .and. (q.eq.qb2)) then ! k neighbors
                  nk=get_index_nn_k(bpt(:,k),nntot,kb,
     >                         gb_kq(:,nn,ikqpt),
     >                         gb(1:3,1:nntot,k),arr_len)
                  nk2=get_index_nn_k(bpt(:,k),nntot,kb2,
     >                         gb_kq(:,nn2,ikqpt),
     >                         gb(1:3,1:nntot,k),arr_len)
                  uHu(:,:,nn2,nn,ikqpt)=uHuk(:,:,nk2,nk,k)*ovaux
c                  write(*,'(a,5i8)')'<k|k>',k,kb,kb2,nk,nk2
               elseif((k.eq.kb) .and. (k.eq.kb2)) then ! q neighbors
                  nq=get_index_nn_q(bpt_q(:,q),nntot_q,qb,
     >                         gb_kq(:,nn,ikqpt),
     >                         gb_q(1:3,1:nntot_q,q),
     >                         arr_len,shift,l_dim)
                  nq2=get_index_nn_q(bpt_q(:,q),nntot_q,qb2,
     >                         gb_kq(:,nn2,ikqpt),
     >                         gb_q(1:3,1:nntot_q,q),
     >                         arr_len,shift,l_dim)
                  uHu(:,:,nn2,nn,ikqpt)=cmplx(0.,0.)
c                  write(*,'(a,5i8)')'<q|q>',q,qb,qb2,nq,nq2
               else ! otherwise problem
                  if((k.eq.kb).and.(k.ne.kb2).and.
     >               (q.ne.qb).and.(q.eq.qb2)) then
                  nk2=get_index_nn_k(bpt(:,k),nntot,kb2,
     >                         gb_kq(:,nn2,ikqpt),
     >                         gb(1:3,1:nntot,k),arr_len)
                  nq=get_index_nn_q(bpt_q(:,q),nntot_q,qb,
     >                         gb_kq(:,nn,ikqpt),
     >                         gb_q(1:3,1:nntot_q,q),
     >                         arr_len,shift,l_dim)
                  uHu(:,:,nn2,nn,ikqpt)=cmplx(0.,0.)
c                  write(*,'(a,8i8)')'<q|k>',q,k,qb,qb2,kb,kb2,nq,nk2
                  elseif((k.ne.kb).and.(k.eq.kb2).and.
     >                   (q.eq.qb).and.(q.ne.qb2)) then
                  nk=get_index_nn_k(bpt(:,k),nntot,kb,
     >                         gb_kq(:,nn,ikqpt),
     >                         gb(1:3,1:nntot,k),arr_len)
                  nq2=get_index_nn_q(bpt_q(:,q),nntot_q,qb2,
     >                         gb_kq(:,nn2,ikqpt),
     >                         gb_q(1:3,1:nntot_q,q),
     >                         arr_len,shift,l_dim)
                  uHu(:,:,nn2,nn,ikqpt)=cmplx(0.,0.)
c                  write(*,'(a,8i8)')'<k|q>',q,k,qb,qb2,kb,kb2,nk,nq2
                  else
                     call juDFT_error("problem k,q neighbors",
     >                                calledby="wann_uHu_commat")
                  endif
               endif
               enddo

             enddo!nn
           enddo!k
         enddo!q
         deallocate(b1,b2,g1,g2,uHuk)

         open(307,file='WF1_gwf.uHu',form='unformatted')
         write(307)nbnd,nkqpts,nntot_kq
         write(307)bpt_kq,gb_kq
         write(307)uHu
         close(307)
         deallocate(uHu)
      else
         allocate(mk(nbnd,nbnd,nntot,nntot))
         open(307,file='WF1_gwf.uHu')
         write(307,*)'uHu for HDWFs'
         write(307,'(i5,i7,i5)')nbnd,nkqpts,nntot_kq

         do q=1,nqpts       
           open(405,file=fuHu(q))
           read(405,*)!title
           read(405,*)b,k,nn
           do k=1,1!nkpts
             ikqpt=get_index_kq(k,q,nkpts)           

             mk = cmplx(0.,0.)
             do nn=1,nntot ! read k-overlaps
               kb = bpt(nn,k)
               do nn2=1,nntot
               kb2= bpt(nn2,k)
               read(405,*)nkp,b,i!g(1:3)
               do i=1,nbnd
                  do j=1,nbnd
                     read(405,*)mmn_r,mmn_i
                     mk(j,i,nn2,nn)=cmplx(mmn_r,mmn_i)
                  enddo
               enddo              
               enddo 
             enddo  

             do nn=1,nntot_kq ! write composite overlaps
               kqb = bpt_kq(nn,ikqpt)
               kb=get_index_k(kqb,nkpts)
               qb=get_index_q(kqb,nkpts)

               do nn2=1,nntot_kq
               kqb2 = bpt_kq(nn2,ikqpt)
               kb2=get_index_k(kqb2,nkpts)
               qb2=get_index_q(kqb2,nkpts)

               write(307,'(3i7)')ikqpt,kqb,kqb2!gb_kq(1:arr_len,nn,ikqpt)
               if((q.eq.qb).and.(q.eq.qb2)) then ! k-overlap
                  !write(307,*)'<k|k>'
                  nk=get_index_nn_k(bpt(:,k),nntot,kb,
     >                         gb_kq(:,nn,ikqpt),
     >                         gb(1:3,1:nntot,k),arr_len)
                  nk2=get_index_nn_k(bpt(:,k),nntot,kb2,
     >                         gb_kq(:,nn2,ikqpt),
     >                         gb(1:3,1:nntot,k),arr_len)
                  do i=1,nbnd
                     do j=1,nbnd
                        write(307,'(2f24.18)')
     >                      real(mk(j,i,nk2,nk)),aimag(mk(j,i,nk2,nk))
                     enddo
                  enddo
               elseif((k.eq.kb).and.(k.eq.kb2)) then ! q-overlap

c                  write(*,*)'<q|q>',qb,qb2

                  nq=get_index_nn_q(bpt_q(:,q),nntot_q,qb,
     >                         gb_kq(:,nn,ikqpt),
     >                         gb_q(1:3,1:nntot_q,q),
     >                         arr_len,shift,l_dim)
                  nq2=get_index_nn_q(bpt_q(:,q),nntot_q,qb2,
     >                         gb_kq(:,nn2,ikqpt),
     >                         gb_q(1:3,1:nntot_q,q),
     >                         arr_len,shift,l_dim)

c                  CALL wann_gwf_auxovlp( param_vec(:,qb)+gb_q(:,nq,q),
c     >              param_vec(:,qb2)+gb_q(:,nq2,q),latt_const_q,ovaux)

                  do i=1,nbnd
                     do j=1,nbnd
                        write(307,'(2f24.18)') 0.0,0.0
c     >                      real(mq(j,i,nq)),aimag(mq(j,i,nq))
                     enddo
                  enddo

               else 

                  if((k.eq.kb).and.(k.ne.kb2).and.
     >               (q.ne.qb).and.(q.eq.qb2)) then
c                  write(*,*)'<q|k>',qb,qb2
                  nk2=get_index_nn_k(bpt(:,k),nntot,kb2,
     >                         gb_kq(:,nn2,ikqpt),
     >                         gb(1:3,1:nntot,k),arr_len)
                  nq=get_index_nn_q(bpt_q(:,q),nntot_q,qb,
     >                         gb_kq(:,nn,ikqpt),
     >                         gb_q(1:3,1:nntot_q,q),
     >                         arr_len,shift,l_dim)

c                  CALL wann_gwf_auxovlp( param_vec(:,qb)+gb_q(:,nq,q),
c     >              param_vec(:,qb2),latt_const_q,ovaux)

                  do i=1,nbnd
                     do j=1,nbnd
                        write(307,'(2f24.18)') 0.0,0.0
c     >                      real(mq(j,i,nq)),aimag(mq(j,i,nq))
                     enddo
                  enddo
                  
                  elseif((k.ne.kb).and.(k.eq.kb2).and.
     >                   (q.eq.qb).and.(q.ne.qb2)) then
c                  write(*,*)'<k|q>',qb,qb2
                  nk=get_index_nn_k(bpt(:,k),nntot,kb,
     >                         gb_kq(:,nn,ikqpt),
     >                         gb(1:3,1:nntot,k),arr_len)
                  nq2=get_index_nn_q(bpt_q(:,q),nntot_q,qb2,
     >                         gb_kq(:,nn2,ikqpt),
     >                         gb_q(1:3,1:nntot_q,q),
     >                         arr_len,shift,l_dim)

c                  CALL wann_gwf_auxovlp( param_vec(:,qb),
c     >              param_vec(:,qb2)+gb_q(:,nq2,q),latt_const_q,ovaux)

                  do i=1,nbnd
                     do j=1,nbnd
                        write(307,'(2f24.18)') 0.0,0.0
c     >                      real(mq(j,i,nq)),aimag(mq(j,i,nq))
                     enddo
                  enddo

                  else
                     call juDFT_error("problem k,q neighbors",
     >                                calledby="wann_uHu_commat")
                  endif
               endif
               enddo

             enddo!nn
           enddo!k
           close(405)
         enddo!q           
         close(307)
         deallocate(mk)
      endif!l_unformatted

      deallocate(gb_kq,bpt_kq)
      deallocate(fuHu)
      deallocate(fuHu2)

      end subroutine wann_uHu_commat    
      end module m_wann_uHu_commat
