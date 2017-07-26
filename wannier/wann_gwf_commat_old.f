      module m_wann_gwf_commat
      USE m_fleurenv
      implicit none

      contains

      subroutine wann_gwf_commat(nbnd,nwfs,nkpts,nntot,bpt,
     >                     nqpts,nntot_q,bpt_q,gb,gb_q,
     >                     mmnk,mmnq,amn,eig,latt_const_q)
      use m_wann_gwf_tools
      use m_wann_gwf_write_mmnk
      use m_wann_write_amn
      use m_constants, only : pimach

      implicit none
      integer,intent(in) :: nbnd,nkpts,nntot,nqpts,nntot_q,nwfs
      integer,intent(in) :: bpt(nntot,nkpts),bpt_q(nntot_q,nqpts)
      integer,intent(in) :: gb(3,nntot,nkpts),gb_q(3,nntot_q,nqpts)
      complex,intent(in) :: mmnk(nbnd,nbnd,nntot,nkpts,nqpts)
      complex,intent(inout) :: mmnq(nbnd,nbnd,nntot_q,nqpts,nkpts)
      complex,intent(inout) :: amn(nbnd,nwfs,nkpts,nqpts)
      real,   intent(in) :: eig(nbnd,nkpts,nqpts)
      real,   intent(in) :: latt_const_q

      integer :: iqpt,iqpt_b,i,j,ikqpt,ikqpt_b,iqpt_help
      integer :: ikpt,ikpt_b,nbnd_arti,nqpts_arti,nntot_arti
      integer :: nkqpts,nntot_kq,nn,ikqpt_help
      integer,allocatable :: gb_kq(:,:,:),bpt_kq(:,:)
      integer,allocatable :: bpt_arti(:,:)
      integer :: dump,nn_arti
      logical :: l_file,l_debug
      complex :: temp_mmn
      complex,allocatable :: mmn_arti(:,:),amn_arti(:)
      complex,allocatable :: mmn_comp(:,:,:,:)
      complex,allocatable :: amn_comp(:,:,:)
      real,allocatable :: eig_comp(:,:)
      real :: mmn_r,mmn_i,tau,tpi,a_arti,b_arti
      character(len=20) :: filename

      l_debug = .false.
      tpi = 2.0*pimach()
      a_arti = latt_const_q
      b_arti = 0.98*latt_const_q

! read in overlaps to artifical potential
      allocate(mmn_arti(nntot_q,nqpts))
      allocate(bpt_arti(nntot_q,nqpts))
      mmn_arti = cmplx(0.,0.)

      tau = tpi/real(nqpts)/a_arti 
      temp_mmn = 2.0*tpi*tpi*sin(tau*b_arti/2.0)
      temp_mmn = temp_mmn/(tpi*tpi-tau*tau*b_arti*b_arti)
      temp_mmn = temp_mmn/(tau*b_arti)

!      temp_mmn=1.0

      write(*,*) temp_mmn

c      inquire(file='WF1_arti.mmn',exist=l_file)
c      if(.not.l_file) 
c     >  call fleur_err("provide WF1_arti.mmn",calledby ="wann_gwf_mmkb")
c
c      open(200,file='WF1_arti.mmn',status='old',action='read')
c      read(200,*)!header
c      read(200,*)nbnd_arti,nqpts_arti,nntot_arti
c
c      if(nbnd_arti.ne.1 .or. nqpts_arti.ne.nqpts 
c     >                  .or. nntot_arti.ne.nntot_q) then
c         close(200)
c         call fleur_err("check format WF1_arti.mmn",
c     >                  calledby="wann_gwf_mmkb")
c      endif
c
c      do iqpt=1,nqpts
c         do iqpt_b=1,nntot_q
c            read(200,*)iqpt_help,bpt_arti(iqpt_b,iqpt),dump,dump,dump
c            read(200,*)mmn_r,mmn_i
c            mmn_arti(iqpt_b,iqpt) = cmplx(mmn_r,-mmn_i)   
c         enddo
c      enddo
c      close(200)

! read in projections to artificial potential
      allocate(amn_arti(nqpts))
      amn_arti = cmplx(0.,0.)

      amn_arti = 1.0

c      inquire(file='WF1_arti.amn',exist=l_file)
c      if(.not.l_file) 
c     >  call fleur_err("provide WF1_arti.amn",calledby ="wann_gwf_mmkb")
c
c      open(200,file='WF1_arti.amn',status='old',action='read')
c      read(200,*)!header
c      read(200,*)nbnd_arti,nqpts_arti,nntot_arti
c
c      if(nbnd_arti.ne.1 .or. nqpts_arti.ne.nqpts 
c     >                  .or. nntot_arti.ne.1) then
c         close(200)
c         call fleur_err("check format WF1_arti.amn",
c     >                  calledby="wann_gwf_mmkb")
c      endif

c      do iqpt=1,nqpts
c         read(200,*)nbnd_arti,nbnd_arti,iqpt_help,mmn_r,mmn_i
c         if(iqpt_help.ne.iqpt) call fleur_err("iqpt_help.ne.iqpt",
c     >                                        calledby="wann_gwf_mmkb")
c         amn_arti(iqpt) = cmplx(mmn_r,mmn_i) 
c      enddo
c      close(200)

c      call fleur_end("stop in wann_gwf_mmkb")


! read in bkqpts
       nkqpts=nkpts*nqpts
       inquire (file='bkqpts',exist=l_file)
       if (.not.l_file)  CALL fleur_err("need bkqpts for l_gwf"
     +      ,calledby ="wann_gwf_mmkb")
       open (202,file='bkqpts',form='formatted',status='old')
       rewind (202)
       read (202,'(i4)') nntot_kq
       write (*,*) 'nntot_kq=',nntot_kq
       allocate ( gb_kq(4,nntot_kq,nkqpts),
     &            bpt_kq(nntot_kq,nkqpts))
       do ikqpt=1,nkqpts
        do nn=1,nntot_kq
         read (202,'(2i6,3x,4i4)')
     &     ikqpt_help,bpt_kq(nn,ikqpt),(gb_kq(i,nn,ikqpt),i=1,4)
         if (ikqpt/=ikqpt_help)  CALL fleur_err("ikqpt.ne.ikqpt_help"
     +        ,calledby ="wann_gwf_mmkb")
         if (bpt_kq(nn,ikqpt)>nkqpts)
     &        CALL fleur_err("bpt_kq.gt.nkqpts",
     >                       calledby ="wann_gwf_mmkb")
c         write(*,'(2i6,3x,4i4)')ikqpt_help,bpt_kq(nn,ikqpt),
c     >                          (gb_kq(i,nn,ikqpt),i=1,4)
        enddo
       enddo
       close (202)

       !call fleur_end("stop")

! construct composite mmn and amn matrices and eig
      allocate(mmn_comp(nbnd,nbnd,nntot_kq,nkqpts))
      allocate(amn_comp(nbnd,nwfs,nkqpts))
      allocate(eig_comp(nbnd,nkqpts))
      mmn_comp=cmplx(0.,0.)
      amn_comp=cmplx(0.,0.)
      eig_comp=cmplx(0.,0.)

      do ikqpt=1,nkqpts
         ikpt=get_index_k(ikqpt,nkpts)
         iqpt=get_index_q(ikqpt,nkpts)
         
         ! composite amn
         amn_comp(:,:,ikqpt) = amn(:,:,ikpt,iqpt)*amn_arti(iqpt)
         
         ! composite eig
         eig_comp(:,ikqpt) = eig(:,ikpt,iqpt)

         ! composite mmn
         do ikqpt_b=1,nntot_kq
            ikpt_b=get_index_k(bpt_kq(ikqpt_b,ikqpt),nkpts)
            iqpt_b=get_index_q(bpt_kq(ikqpt_b,ikqpt),nkpts)

c            write(*,*)'kq',ikqpt,'kqb',ikqpt_b
c            write(*,*)'k',ikpt,'kb',ikpt_b,'q',iqpt,'qb',iqpt_b

            if(iqpt.eq.iqpt_b) then !overlap in k
               nn=get_index_nn_k(bpt(:,ikpt),nntot,ikpt_b,
     >                         gb_kq(:,ikqpt_b,ikqpt),
     >                         gb(1:3,1:nntot,ikpt))
c               write(*,*)'k neighbor nr.',nn
               mmn_comp(:,:,ikqpt_b,ikqpt) = mmnk(:,:,nn,ikpt,iqpt)
            elseif(ikpt.eq.ikpt_b) then !overlap in q
               nn=get_index_nn_q(bpt_q(:,iqpt),nntot_q,iqpt_b,
     >                           gb_kq(:,ikqpt_b,ikqpt),
     >                           gb_q(1:3,1:nntot_q,iqpt))
c               write(*,*)'q neighbor nr.',nn
c               do nn_arti=1,nntot_q
c                  if(bpt_arti(nn_arti,iqpt).eq.bpt_q(nn,iqpt)) exit
c                  if(nn_arti.eq.nntot_q) call fleur_err(
c     >                   "nn_arti not found",calledby="wann_gwf_mmkb")
c               enddo
c               write(*,*)'q nn',nn,'arti nn',nn_arti
               mmn_comp(:,:,ikqpt_b,ikqpt) = mmnq(:,:,nn,iqpt,ikpt)
     >                                     * temp_mmn
c     >                                     * mmn_arti(nn_arti,iqpt) 
            else !problem
               call fleur_err("overlap mmn gwf",
     >                        calledby="wann_gwf_mmkb")
            endif
            !write(*,*)

         enddo
      enddo

! write composite mmn
      call wann_gwf_write_mmnk("WF1_gwf",nkqpts,nntot_kq,nbnd,bpt_kq,
     >                         gb_kq,mmn_comp)

! write composite amn
      call wann_write_amn(.true.,"WF1_gwf.amn",                  
     >           'Overlaps of the wavefunct. with the trial orbitals',
     >           nbnd,nkqpts,nwfs,1,1,.false.,.true.,amn_comp)

! write composite eig
      open(300,file="WF1_gwf.eig",form="formatted",
     >     status="unknown",action="write")
      do ikqpt=1,nkqpts
         do i=1,nbnd
            write(300,'(2i12,f19.13)')i,ikqpt,eig_comp(i,ikqpt)
         enddo
      enddo
      close(300)

! further output for testing and debugging
      if(l_debug) then
      do ikpt=1,nkpts

         write(filename,'("WF1_q_",i4.4)')ikpt 
         
         ! write eig_q file
         open(300,file=trim(filename)//".eig",
     >        status='unknown',action='write')
         do iqpt=1,nqpts
            do i=1,nbnd
               write(300,*)i,iqpt,eig(i,ikpt,iqpt)
            enddo
         enddo
         close(300)

         ! write mmnq
         call wann_gwf_write_mmnk(filename,nqpts,nntot_q,nbnd,bpt_q,
     >                    gb_q,mmnq(:,:,:,:,ikpt))

         ! write composite amn
         call wann_write_amn(.true.,trim(filename)//".amn",                  
     >           'Overlaps of the wavefunct. with the trial orbitals',
     >           nbnd,nqpts,nwfs,1,1,.false.,.true.,amn(:,:,ikpt,:))
      enddo
      endif!debugging
     
      deallocate(mmn_arti,mmn_comp)
      deallocate(amn_comp,eig_comp)
      deallocate(gb_kq,bpt_kq)

      end subroutine wann_gwf_commat    


      end module m_wann_gwf_commat
