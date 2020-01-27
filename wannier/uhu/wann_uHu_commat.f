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
      logical,allocatable :: l_fm(:),l_fm2(:)
      integer,allocatable :: gb_kq(:,:,:)
      integer,allocatable :: bpt_kq(:,:)
      integer,allocatable :: g1(:,:,:),g2(:,:)
      integer,allocatable :: b1(:,:),b2(:)
      complex,allocatable :: tmp_a(:,:,:),tmp_m(:,:,:,:)
      complex,allocatable :: uHukk(:,:,:,:,:)
      complex,allocatable :: uHukq(:,:,:,:,:)
      complex,allocatable :: mkk(:,:,:,:)
      complex,allocatable :: mkq(:,:,:,:)
      complex,allocatable :: uHu(:,:,:)
      character(len=20),allocatable :: fuHu(:),fuHu2(:)

      ! further variables
      integer :: nbnd,nwfs,kb,qb,kb2,qb2,nn_kq,kqb,kqb2,arr_len,shift(3)
      integer :: i,j,ikqpt,ikqpt_b,nk,nq,nk2,nq2,g(3)
      integer :: nbnd_arti,nqpts_arti,nntot_arti,nnkq_len,nnkq_ind
      integer :: nkqpts,nntot_kq,nn,nn2,ikqpt_help,testnn
      integer :: dump,nn_arti,nkp,k,b,q,nwf,n,is,ie,countkq,countqk
      real :: mmn_r,mmn_i,tau,tpi,a_arti,b_arti,eig
      complex :: ovaux
      logical :: l_exist,l_miss,l_fe,l_fa,l_proj
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
      a_arti = latt_const_q
      b_arti = 0.98*latt_const_q
      tau = tpi/real(nqpts)/a_arti
      ovaux = 2.0*tpi*tpi*sin(tau*b_arti/2.0)
      ovaux = ovaux/(tpi*tpi-tau*tau*b_arti*b_arti)
      ovaux = ovaux/(tau*b_arti)
      write(*,*)'ov(0)=',ovaux

      allocate(fuHu(nqpts),fuHu2(nqpts),l_fm(nqpts),l_fm2(nqpts))

! are all necessary files present?
      l_miss=.false.
      do q=1,nqpts
         write(fq,'(i4.4)')q
         fuHu(q) = 'WF1_'//trim(fq)//'.uHu'    ! kk-overlaps
         fuHu2(q)= 'WF1_'//trim(fq)//'.uHu_kq'  ! kq-overlaps

         inquire(file=fuHu(q) ,exist=l_fm(q) )
         inquire(file=fuHu2(q),exist=l_fm2(q))

         if(.not.l_fm(q) ) write(*,*)'missing: ',fuHu(q)
         if(.not.l_fm2(q)) write(*,*)'missing: ',fuHu2(q)
         if(.not.(l_fm(q).or.l_fm2(q))) l_miss=.true.
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
         nnkq_len = (nntot_kq+1)*nntot_kq/2
         allocate(uHukk(nbnd,nbnd,nntot,nntot,nkpts))
         allocate(uHukq(nbnd,nbnd,nntot_q,nntot,nkpts))
         allocate(g1(3,nntot,nkpts))
         allocate(b1(nntot,nkpts))
c         allocate(uHu(nbnd,nbnd,nntot_kq,nntot_kq))!,nkqpts))
         allocate(uHu(nbnd,nbnd,nnkq_len))

         open(307,file='WF1_gwf.uHu',form='unformatted')
         write(307)nbnd,nkqpts,nntot_kq,nnkq_len
         write(307)bpt_kq,gb_kq

         !open(853,file='debug_commat')

         countkq=0
         countqk=0

         do q=1,nqpts
           is = get_index_kq(1,q,nkpts)
           ie = get_index_kq(nkpts,q,nkpts)
         
           if(l_fm(q)) then
            open(405,file=fuHu(q),form='unformatted')
            read(405)b,k,nn
            read(405)b1,g1
            read(405)uHukk
            close(405,status='delete')
           else
            uHukk=cmplx(0.,0.)
            write(*,*)'k-k overlaps set to zero for q=',q
           endif

           if(l_fm2(q)) then
            open(405,file=fuHu2(q),form='unformatted')
            read(405)b,k,nn
            read(405)b1,g1
            read(405)uHukq
            close(405,status='delete')
           else
            uHukq=cmplx(0.,0.)
            write(*,*)'k-q overlaps set to zero for q=',q
           endif

c           mmnq = mmnq*conjg(mmn_arti)

           do k=1,nkpts
             uHu = cmplx(0.,0.)
             ikqpt=get_index_kq(k,q,nkpts)
             nnkq_ind=0 
             do nn=1,nntot_kq
               kqb = bpt_kq(nn,ikqpt)
               kb=get_index_k(kqb,nkpts)
               qb=get_index_q(kqb,nkpts)

               nk=get_index_nn_k(bpt(:,k),nntot,kb,
     >                         gb_kq(:,nn,ikqpt),
     >                         gb(1:3,1:nntot,k),arr_len)
               nq=get_index_nn_q(bpt_q(:,q),nntot_q,qb,
     >                         gb_kq(:,nn,ikqpt),
     >                         gb_q(1:3,1:nntot_q,q),
     >                         arr_len,shift,l_dim)

               do nn2=1,nn!nntot_kq
               nnkq_ind = nnkq_ind+1
               kqb2= bpt_kq(nn2,ikqpt)
               kb2=get_index_k(kqb2,nkpts)
               qb2=get_index_q(kqb2,nkpts)

               nk2=get_index_nn_k(bpt(:,k),nntot,kb2,
     >                         gb_kq(:,nn2,ikqpt),
     >                         gb(1:3,1:nntot,k),arr_len)
               nq2=get_index_nn_q(bpt_q(:,q),nntot_q,qb2,
     >                         gb_kq(:,nn2,ikqpt),
     >                         gb_q(1:3,1:nntot_q,q),
     >                         arr_len,shift,l_dim)
               testnn = nk/abs(nk)+nk2/abs(nk2)+nq/abs(nq)+nq2/abs(nq2)
               if(testnn.ne.0) stop 'testnn.ne.0'
c               write(*,*)'kq',ikqpt,nn,nn2
c               write(*,*)'k ',k,nk,nk2
c               write(*,*)'q ',q,nq,nq2

               !write(853,'(a,1x,5(5i))')'kq,nkq',ikqpt,kqb,kqb2,nn,nn2
               !write(853,'(a,1x,5(5i))')'k,nk',k,kb,kb2,nk,nk2
               !write(853,'(a,1x,5(5i))')'q,nq',q,qb,qb2,nq,nq2
               !write(853,*)'nnkq_ind',nnkq_ind,nnkq_len

c               CALL wann_gwf_auxovlp(param_vec(:,qb),param_vec(:,qb2),
c     >                               latt_const_q,ovaux)
c               if(ovaux.ne.1.0) write(*,*)'ovaux',ovaux

               if((nq.lt.0) .and. (nq2.lt.0)) then ! kk
                  uHu(:,:,nnkq_ind)=uHukk(:,:,nk2,nk,k)!*ovaux
c                  write(*,'(a,5i8)')'<k|k>'
                  !write(853,*)'<k|k>'
               elseif((nk.lt.0) .and. (nk2.lt.0)) then ! qq
                  uHu(:,:,nnkq_ind)=cmplx(0.,0.)
c                  write(*,'(a,5i8)')'<q|q>'
                  !write(853,*)'<q|q>'
               else
                  if((nk.lt.0) .and. (nk2.ge.0) .and.
     >               (nq.ge.0) .and. (nq2.lt.0) ) then
                     uHu(:,:,nnkq_ind)
     >                  = transpose(conjg(uHukq(:,:,nq,nk2,k)))
     >                    *conjg(ovaux)
                  !write(853,*)'<q|k>'
c                     write(*,'(a,8i8)')'<q|k>'
c                  write(*,*)nn,nn2
c                  write(*,*)'<q|k>',q,bpt_q(nq,q),k,bpt(nk2,k)
                  if(q.eq.bpt_q(nq,q)) stop 'q.eq.bpt_q'
                  countqk=countqk+1
                  elseif((nk.ge.0) .and. (nk2.lt.0) .and. 
     >                   (nq.lt.0) .and. (nq2.ge.0) ) then
                     uHu(:,:,nnkq_ind)
     >                  = uHukq(:,:,nq2,nk,k) * ovaux
                  !write(853,*)'<k|q>'
c                     write(*,'(a,8i8)')'<k|q>'
c                  write(*,*)nn,nn2
c                  write(*,*)'<k|q>',k,bpt(nk,k),q,bpt_q(nq2,q)
                  if(q.eq.bpt_q(nq2,q)) stop 'q.eq.bpt_q'
                  countkq=countkq+1
                  else
                     call juDFT_error("problem k,q neighbors",
     >                                calledby="wann_uHu_commat")
                  endif
               endif
               enddo

             enddo!nn
             write(307)ikqpt,uHu
           enddo!k
         enddo!q
         deallocate(b1,g1,uHukk,uHukq)

         !close(853)
         close(307)
         deallocate(uHu)

         write(*,*)'countkq',countkq
         write(*,*)'countqk',countqk

      else
         allocate(mkk(nbnd,nbnd,nntot,nntot))
         allocate(mkq(nbnd,nbnd,nntot_q,nntot))
         open(307,file='WF1_gwf.uHu')
         write(307,*)'uHu for HDWFs'
         write(307,'(i5,i7,i5)')nbnd,nkqpts,nntot_kq

         do q=1,nqpts      
 
           if(l_fm(q)) then
            open(405,file=fuHu(q))
            read(405,*)!title
            read(405,*)b,k,nn
           endif
           if(l_fm2(q)) then
            open(505,file=fuHu2(q))
            read(505,*)!title
            read(505,*)b,k,nn
           endif

           do k=1,nkpts
             ikqpt=get_index_kq(k,q,nkpts)           

             mkk = cmplx(0.,0.)
             if(l_fm(q)) then
              do nn=1,nntot
               do nn2=1,nntot
                read(405,*)nkp,b,i!g(1:3)
                do i=1,nbnd
                   do j=1,nbnd
                      read(405,*)mmn_r,mmn_i
                      mkk(j,i,nn2,nn)=cmplx(mmn_r,mmn_i)
                   enddo
                enddo              
               enddo 
              enddo
             endif

             mkq = cmplx(0.,0.)
             if(l_fm2(q)) then
              do nn=1,nntot
               do nn2=1,nntot_q
                read(505,*)nkp,b,i!g(1:3)
                do i=1,nbnd
                   do j=1,nbnd
                      read(505,*)mmn_r,mmn_i
                      mkq(j,i,nn2,nn)=cmplx(mmn_r,mmn_i)
                   enddo
                enddo              
               enddo 
              enddo
              mkq = mkq*ovaux
             endif

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
     >                      real(mkk(j,i,nk2,nk)),aimag(mkk(j,i,nk2,nk))
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
                        write(307,'(2f24.18)')
     >                    real(mkq(i,j,nq,nk2)),-aimag(mkq(i,j,nq,nk2))
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
                        write(307,'(2f24.18)')
     >                    real(mkq(j,i,nq2,nk)),aimag(mkq(j,i,nq2,nk))
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
           close(505)
         enddo!q           
         close(307)
         deallocate(mkk,mkq)
      endif!l_unformatted

      deallocate(gb_kq,bpt_kq)
      deallocate(fuHu,l_fm)
      deallocate(fuHu2,l_fm2)

      end subroutine wann_uHu_commat    
      end module m_wann_uHu_commat
