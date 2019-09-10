      MODULE m_wann_convert_fleur_w90
      contains
      subroutine wann_convert_fleur_w90(jspins_in,l_nocosoc,wann)
      use m_types
      use m_juDFT
      IMPLICIT NONE
      integer,intent(in) :: jspins_in
      logical,intent(in) :: l_nocosoc
      type(t_wann), intent(in) :: wann

      logical :: l_readunf,l_writeunf,l_conjg
      logical :: l_spinmat,l_paulimat
      character(len=20) :: filenamewrite
      character(len=15) :: writeform
      character(len=20) :: filenameread(3)
      character(len=15) :: readform
      integer :: filestoread,fileidx
      integer :: spn_in,nb_tmp,num_bands,nkp_tmp
      integer :: num_kpts,ierr,ik,counter,m,n,compo,s
      character(len=60) :: header
      complex,allocatable :: oper_temp(:,:)
      complex,allocatable :: oper_o(:,:,:,:)
      integer :: num_compos,jspins,nbnd,fullnkpts,nwfs,nkp,i,j
      integer :: num_bands1,num_bands2
      complex,allocatable :: matrix6(:,:,:,:,:,:)
      real :: a,b
      real :: s_real,s_img
      complex,parameter   :: ci=(0.0,1.0)
      integer :: dummy1,dummy2,dummy3,dummy4,dummy5,dummy6
      integer :: spin1,spin2,spinmat_dims,ikpt,dir,ii,jj
      real :: conversionfactor
      real,parameter      :: hart=27.21138505
      integer :: jj_tmp,ii_tmp,j_tmp,i_tmp,dir_tmp,ikpt_tmp,test_tmp
      integer :: write_bands,firstband
      integer :: map3(3)
      
      do ik=1,3
         map3(ik)=ik
      enddo   

      if(l_nocosoc)then
         jspins=1
      else
         jspins=min(jspins_in,2)
      endif

! Determine number of kpoints.
      spn_in=916
      open(spn_in,file='WF1.amn',form='formatted')
      read(spn_in,*)
      read(spn_in,*)nbnd,num_kpts,nwfs
      close(spn_in)


! For hsomtxvec we need the total number of bands.
      num_bands1=wann%band_max(1)-wann%band_min(1)+1
      num_bands2=wann%band_max(2)-wann%band_min(2)+1

      if(l_nocosoc)then
        num_bands=num_bands1
      else  !assume that we use socinterpol
         num_bands=num_bands1+num_bands2
      endif   


      if(nbnd.ne.num_bands1) then
        write(*,*)"num_bands1=",num_bands1
        write(*,*)"nbnd=",nbnd
        call juDFT_error('discrepancy convert fleur_w90')
      endif   
      
      ! Defaults
      filestoread=1
      l_conjg=.false.
      l_spinmat=.false.
      l_paulimat=.false.
      num_compos=3
      write_bands=num_bands
      if(wann%l_mmn0_unf_to_spn_unf)then
         l_readunf=.true.
         readform='unformatted'
         l_writeunf=.true.
         writeform='unformatted'
         filenameread(1)='updown.mmn0_unf'
         filestoread=3
         if(l_nocosoc)then
             filenameread(2)='WF1.socmmn0_unf'
             filenameread(3)='WF2.socmmn0_unf'
         else
             filenameread(2)='WF1.mmn0_unf'
             filenameread(3)='WF2.mmn0_unf'
         endif
         filenamewrite='WF1.spn'
         l_conjg=.true.
         conversionfactor=1.0
         l_paulimat=.true.
         
      elseif(wann%l_mmn0_to_spn_unf)then
         l_readunf=.false.
         readform='formatted'
         l_writeunf=.true.
         writeform='unformatted'
         filenameread(1)='updown.mmn0'
         filestoread=3
         if(l_nocosoc)then
            filenameread(2)='WF1.socmmn0'
            filenameread(3)='WF2.socmmn0'
         else
            filenameread(2)='WF1.mmn0'
            filenameread(3)='WF2.mmn0'
         endif
         filenamewrite='WF1.spn'
         l_conjg=.true.
         conversionfactor=1.0
         l_paulimat=.true.
         
      elseif(wann%l_mmn0_to_spn)then
         l_readunf=.false.
         readform='formatted'
         l_writeunf=.false.
         writeform='formatted'
         filenameread(1)='updown.mmn0'
         filestoread=3
         if(l_nocosoc)then
            filenameread(2)='WF1.socmmn0'
            filenameread(3)='WF2.socmmn0'
         else
            filenameread(2)='WF1.mmn0'
            filenameread(3)='WF2.mmn0'
         endif
         filenamewrite='WF1.spn'
         l_conjg=.true.
         conversionfactor=1.0
         l_paulimat=.true.

        
      elseif(wann%l_mmn0_to_spn2)then ! when socinterpolation is used
         l_readunf=.false.
         readform='formatted'
         l_writeunf=.false.
         writeform='formatted'
         filenameread(1)='updown.mmn0'
         filestoread=3
         if(l_nocosoc)then
            filenameread(2)='WF1.socmmn0'
            filenameread(3)='WF2.socmmn0'
         else
            filenameread(2)='WF1.mmn0'
            filenameread(3)='WF2.mmn0'
         endif
         filenamewrite='WF1.spn'
         l_conjg=.true.
         conversionfactor=1.0
         l_paulimat=.true.
         write_bands=num_bands/2



      elseif(wann%l_mmn0_unf_to_spn)then
         l_readunf=.true.
         readform='unformatted'
         l_writeunf=.false.
         writeform='formatted'
         filenameread(1)='updown.mmn0_unf'
         filestoread=3
         if(l_nocosoc)then
            filenameread(2)='WF1.socmmn0_unf'
            filenameread(3)='WF2.socmmn0_unf'
         else
            filenameread(2)='WF1.mmn0_unf'
            filenameread(3)='WF2.mmn0_unf'
         endif
         filenamewrite='WF1.spn'
         l_conjg=.true.
         conversionfactor=1.0
         l_paulimat=.true.

        
      elseif(wann%l_perpmag_unf_to_tor_unf)then
         l_readunf=.true.
         readform='unformatted'
         l_writeunf=.true.
         writeform='unformatted'
         filestoread=1
         filenameread(1)='updown.perpmag_unf'
         filenamewrite='WF1.tor'
         conversionfactor=hart
         map3(1)=2
         map3(2)=1
        
      elseif(wann%l_perpmag_to_tor_unf)then
         l_readunf=.false.
         readform='formatted'
         l_writeunf=.true.
         writeform='unformatted'
         filestoread=1
         filenameread(1)='updown.perpmag'
         filenamewrite='WF1.tor'
         conversionfactor=hart
         map3(1)=2
         map3(2)=1
        
      elseif(wann%l_perpmag_to_tor)then
         l_readunf=.false.
         readform='formatted'
         l_writeunf=.false.
         writeform='formatted'
         filestoread=1
         filenameread(1)='updown.perpmag'
         filenamewrite='WF1.tor'
         conversionfactor=hart
         map3(1)=2
         map3(2)=1
         
      elseif(wann%l_perpmag_unf_to_tor)then
         l_readunf=.true.
         readform='unformatted'
         l_writeunf=.false.
         writeform='formatted'
         filestoread=1
         filenameread(1)='updown.perpmag_unf'
         filenamewrite='WF1.tor'
         conversionfactor=hart
         map3(1)=2
         map3(2)=1
        
      elseif(wann%l_hsomtxvec_unf_to_lmpzsoc_unf)then
         l_readunf=.true.
         readform='unformatted'
         l_writeunf=.true.
         writeform='unformatted'
         filestoread=1
         filenameread(1)='WF1.hsomtxvec_unf'
         filenamewrite='WF1.lmpzsoc'
         l_spinmat=.true.
         l_conjg=.true.
         spinmat_dims=3
         conversionfactor=hart         
!         write_bands=num_bands1
         
         if(l_nocosoc) call juDFT_error('noco_or_soc and hsomtxvec')      
      elseif(wann%l_hsomtxvec_to_lmpzsoc_unf)then
         l_readunf=.false.
         readform='formatted'
         l_writeunf=.true.
         writeform='unformatted'
         filestoread=1
         spinmat_dims=3
         filenameread(1)='WF1.hsomtxvec'
         filenamewrite='WF1.lmpzsoc'
         l_spinmat=.true.
         l_conjg=.true.
         conversionfactor=hart
!         write_bands=num_bands1
        
         if(l_nocosoc) call juDFT_error('noco_or_soc and hsomtxvec')   
      elseif(wann%l_hsomtxvec_to_lmpzsoc)then
         l_readunf=.false.
         readform='formatted'
         l_writeunf=.false.
         writeform='formatted'
         filestoread=1
         filenameread(1)='WF1.hsomtxvec'
         filenamewrite='WF1.lmpzsoc'
         l_spinmat=.true.
         spinmat_dims=3
         l_conjg=.true.
         conversionfactor=hart
!         write_bands=num_bands1
        
         if(l_nocosoc) call juDFT_error('noco_or_soc and hsomtxvec')   
      elseif(wann%l_hsomtxvec_unf_to_lmpzsoc)then
         l_readunf=.true.
         readform='unformatted'
         l_writeunf=.false.
         writeform='formatted'
         filestoread=1
         filenameread(1)='WF1.hsomtxvec_unf'
         filenamewrite='WF1.lmpzsoc'
         l_spinmat=.true.
         spinmat_dims=3
         l_conjg=.true.
         conversionfactor=hart
!         write_bands=num_bands1
        
         if(l_nocosoc) call juDFT_error('noco_or_soc and hsomtxvec')   
      elseif(wann%l_hsomtx_unf_to_hsoc_unf)then
         l_readunf=.true.
         readform='unformatted'
         l_writeunf=.true.
         writeform='unformatted'
         filestoread=1
         filenameread(1)='WF1.hsomtx_unf'
         filenamewrite='WF1.hsoc'
         l_spinmat=.true.
         num_compos=1
           spinmat_dims=1
         l_conjg=.true.
         conversionfactor=hart         
!       
            
      elseif(wann%l_hsomtx_to_hsoc_unf)then
         l_readunf=.false.
         readform='formatted'
         l_writeunf=.true.
         writeform='unformatted'
         filestoread=1
         filenameread(1)='WF1.hsomtx'
         filenamewrite='WF1.hsoc'
         l_spinmat=.true.
         num_compos=1
           spinmat_dims=1
         l_conjg=.true.
         conversionfactor=hart
!         
      elseif(wann%l_hsomtx_to_hsoc)then
         l_readunf=.false.
         readform='formatted'
         l_writeunf=.false.
         writeform='formatted'
         filestoread=1
         filenameread(1)='WF1.hsomtx'
         filenamewrite='WF1.hsoc'
         l_spinmat=.true.
         num_compos=1
           spinmat_dims=1
         l_conjg=.true.
         conversionfactor=hart
!       
      elseif(wann%l_hsomtx_unf_to_hsoc)then
         l_readunf=.true.
         readform='unformatted'
         l_writeunf=.false.
         writeform='formatted'
         filestoread=1
         num_compos=1
         filenameread(1)='WF1.hsomtx_unf'
         filenamewrite='WF1.hsoc'
         l_spinmat=.true.
           spinmat_dims=1
         l_conjg=.true.
         conversionfactor=hart
!       
         
         
      endif   

!---- read data in fleur-format
      spn_in=916
    
      do fileidx=1,filestoread
         write(*,*)"fileidx=",fileidx
         write(*,*)"filenameread(fileidx)=",filenameread(fileidx)
      open(spn_in,file=trim(filenameread(fileidx)),form=trim(readform))
      if(l_readunf)then
         header=trim(filenameread(fileidx))
        if(l_spinmat)then
           
           allocate( matrix6(2,2,num_bands1,num_bands1,spinmat_dims,num_kpts) )
           
           if(spinmat_dims==1)then
            read(spn_in)dummy1,dummy2,dummy3,dummy4,dummy5
           else
            read(spn_in)dummy1,dummy2,dummy3,dummy4,dummy5,dummy6
           endif
           
           do nkp=1,num_kpts
            read(spn_in)matrix6(:,:,:,:,:,nkp)
           enddo 
        else
	  read(spn_in)nbnd,fullnkpts,nwfs
          num_bands=nbnd
          if(.not.l_nocosoc)num_bands=2.0*nbnd
          num_kpts=fullnkpts
          if(.not.allocated(oper_o))then
          allocate(oper_o(num_bands,num_bands,num_kpts,num_compos) )
          oper_o=cmplx(0.0,0.0)
          endif
          write(*,*)"before read(spn_in):nbnd,fullnkpts,fileidx=",nbnd,fullnkpts,fileidx
          if(.not.l_nocosoc  .and.  fileidx==1  )then
            read(spn_in)oper_o(1:nbnd,1+nbnd:2*nbnd,1:fullnkpts,fileidx)
          else
	        read(spn_in)oper_o(1:nbnd,1:nbnd,1:fullnkpts,fileidx)
	      endif
        endif
      else
        if(l_spinmat)then
           spin1=2
           spin2=2
           
!          num_bands=nbnd
           allocate( matrix6(2,2,num_bands,num_bands,spinmat_dims,num_kpts) )
           num_bands1=nbnd
           num_bands2=nbnd
           read(spn_in,*)header
           read(spn_in,*)
           do ikpt=1,num_kpts
            do dir=1,spinmat_dims  
             do i = 1,num_bands2
              do j = 1,num_bands1
               do ii=1,spin1
                do jj=1,spin2
                 if(spinmat_dims==3)then
       read (spn_in,'(6i5,3x,2f18.12)') jj_tmp,ii_tmp,j_tmp,i_tmp,dir_tmp,ikpt_tmp,   a,    b
                test_tmp=abs(jj_tmp-jj)+abs(ii_tmp-ii)+abs(j_tmp-j)+abs(i_tmp-i)+abs(dir_tmp-dir)+abs(ikpt_tmp-ikpt)
                  if(test_tmp.ne.0)call juDFT_error('convert: test_tmp')
                  matrix6(jj,ii,j,i,dir,ikpt) = cmplx(a,b)
                 else
                 read (spn_in,'(5i5,3x,2f18.12)') jj_tmp,ii_tmp,j_tmp,i_tmp,ikpt_tmp,   a,    b
                test_tmp=abs(jj_tmp-jj)+abs(ii_tmp-ii)+abs(j_tmp-j)+abs(i_tmp-i)+abs(ikpt_tmp-ikpt)
                  if(test_tmp.ne.0)call juDFT_error('convert: test_tmp')
                  matrix6(jj,ii,j,i,dir,ikpt) = cmplx(a,b)
                 endif
                enddo !jj 
               enddo !ii
              enddo !j
             enddo !i
            enddo !dir 
           enddo !ikpt
        else
	   read(spn_in,*)header
	   read(spn_in,*)nbnd,fullnkpts,nwfs

           num_kpts=fullnkpts
           if(l_nocosoc)then
             num_bands=nbnd
             firstband=0
           else
             num_bands=2*nbnd
             if(fileidx.eq.1)then
               firstband=nbnd
             else
               firstband=0
             endif
           endif  

           if(.not.allocated(oper_o))then
           allocate(oper_o(num_bands,num_bands,num_kpts,num_compos) )
           endif
           

           do nkp=1,num_kpts
            do i=1,nbnd
             do j=1,nbnd
              read(spn_in,*)dummy1,dummy2,dummy3,a,b
!             paulimat(1,j,i+num_bands*(jspins-1),nkp)=cmplx(a,b)

               oper_o(j,i+firstband,nkp,fileidx)=cmplx(a,b)
             enddo !j
            enddo !i
           enddo !nkp
       endif
      endif
      close(spn_in)
      enddo !filestoread


      write(*,*)"before spinmat:l_spinmat=",l_spinmat
      if(l_spinmat)then


          write(*,*)"matrix6=",matrix6(1,1,2,1,2,1)

!         num_bands1=num_bands
!         num_bands2=num_bands
!         num_bands=2*num_bands
!         write(*,*)"num_bands1,num_bands2,num_bands=",num_bands1,
!     &                   num_bands2,num_bands
         allocate( oper_o(num_bands,num_bands,num_kpts,spinmat_dims) )
         if(.false.)then
!   Old variant with num_bands1==num_bands2:
           if(num_bands1.ne.num_bands2) call juDFT_error('convert: num_bands1.ne.num_bands2')
           do ikpt=1,fullnkpts
            do dir=1,spinmat_dims  
             do i = 1,num_bands2
              do j = 1,num_bands2
               do ii=1,spin1
                do jj=1,spin2
         oper_o(j+(jj-1)*num_bands,i+(ii-1)*num_bands,dir,ikpt)=matrix6(jj,ii,j,i,dir,ikpt)
                enddo !jj 
               enddo !ii
              enddo !j
             enddo !i
            enddo !dir 
           enddo !ikpt
        else   
          if(num_bands1.ne.num_bands2) call juDFT_error('convert: num_bands1.ne.num_bands2')

          if(l_nocosoc)then !compute the full SOC-matrix
           do ikpt=1,num_kpts
            do dir=1,spinmat_dims   
             do i = 1,num_bands1
              do j = 1,num_bands1
                 oper_o(j,i,ikpt,dir)=matrix6(1,1,j,i,dir,ikpt)+ &
                                    matrix6(1,2,j,i,dir,ikpt)+ &
                                    matrix6(2,1,j,i,dir,ikpt)+ &
                                    matrix6(2,2,j,i,dir,ikpt)
              enddo !j
             enddo !i
            enddo !dir
           enddo !ikpt
          else !rewrite the SOC-matrix for the purpose of socinterpol


          write(*,*)"before first loop, fullnkpts,num_dims,num_bands1=",num_kpts,spinmat_dims,num_bands1
           do ikpt=1,num_kpts
            do dir=1,spinmat_dims  

             do i = 1,num_bands1
              do j = 1,num_bands1
                   oper_o(j,i,ikpt,dir)=matrix6(1,1,j,i,dir,ikpt)
!                   write(*,*)"j,i,ikpt,dir,matrix6=",j,i,ikpt,dir,matrix6(1,1,j,i,dir,ikpt)
              enddo !j
             enddo !i

             do i = 1,num_bands1
              do j = 1,num_bands2
          oper_o(j+num_bands1,i,ikpt,dir)=matrix6(2,1,j,i,dir,ikpt)
              enddo !j
             enddo !i

             do i = 1,num_bands2
              do j = 1,num_bands1
         oper_o(j,i+num_bands2,ikpt,dir)=matrix6(1,2,j,i,dir,ikpt)
              enddo !j
             enddo !i

             do i = 1,num_bands2
              do j = 1,num_bands2
                   oper_o(j+num_bands1,i+num_bands1,ikpt,dir)=matrix6(2,2,j,i,dir,ikpt)
              enddo !j
             enddo !i

            enddo !dir 
           enddo !ikpt
          endif !nocosoc?
        endif   

!        write(*,*)"oper_o=",oper_o(2,1,1,2)

      endif

      do nkp=1,num_kpts
!          paulimat(2,:,:,nkp)=ci*paulimat(1,:,:,nkp)
!          paulimat(1,:,:,nkp)=paulimat(1,:,:,nkp)+transpose(conjg( paulimat(1,:,:,nkp) ))
!          paulimat(2,:,:,nkp)=paulimat(2,:,:,nkp)+transpose(conjg( paulimat(2,:,:,nkp) ))

         if(l_paulimat)then
!            write(*,*)"convert paulimat: num_bands=",num_bands
            if(l_nocosoc)then
             do i=1,num_bands1
              do j=1,num_bands1
                 oper_o(j,i,nkp,3)=oper_o(j,i,nkp,2)-oper_o(j,i,nkp,3)
              enddo   
             enddo
            else
             do i=1,num_bands1
              do j=1,num_bands1
                 oper_o(j,i,nkp,3)=oper_o(j,i,nkp,2)
              enddo   
             enddo

!             write(*,*)"paulimat conversion:"
!             write(*,*)"num_bands1=",num_bands1

             do i=1,num_bands1
              do j=1,num_bands1
                 oper_o(j+num_bands1,i+num_bands1,nkp,3)=-1.0*oper_o(j,i,nkp,3)
              enddo   
             enddo

!             do i=1,num_bands1
!                write(*,*)"i,oper_o=",i,oper_o(i,i,nkp,3)
!             enddo !i

             do i=1,num_bands1
              do j=1,num_bands1
                 oper_o(j,i+num_bands1,nkp,1)=oper_o(j,i,nkp,1)
              enddo   
             enddo

             do i=1,num_bands1
              do j=1,num_bands1
                 oper_o(j,i,nkp,1)=cmplx(0.0,0.0)
              enddo   
             enddo

            endif  
         endif


!         l_conjg=.false.
         if(l_conjg)then
!           write(*,*)"complex conjugation: num_bands=",num_bands
          do i=1,num_bands
           do j=1,num_bands
              oper_o(j,i,nkp,:)=conjg(oper_o(j,i,nkp,:))
           enddo   
          enddo
         endif   

         if(.not.l_spinmat)then
         if(l_conjg)then
            oper_o(:,:,nkp,2)=-ci*oper_o(:,:,nkp,1)
         else
            oper_o(:,:,nkp,2)=ci*oper_o(:,:,nkp,1)
         endif
         oper_o(:,:,nkp,1)=      oper_o(:,:,nkp,1)+transpose(conjg(oper_o(:,:,nkp,1)))
         oper_o(:,:,nkp,2)=      oper_o(:,:,nkp,2)+transpose(conjg(oper_o(:,:,nkp,2)))
         endif

      enddo	


!---- write data in w90-format
      spn_in=916
      nb_tmp=write_bands
      nkp_tmp=num_kpts
      open(spn_in,file=trim(filenamewrite),form=trim(writeform))
      if(l_writeunf)then
        write (spn_in) header
        write (spn_in) nb_tmp, nkp_tmp

        allocate (oper_temp(num_compos, (write_bands*(write_bands + 1))/2), stat=ierr)
        if (ierr /= 0) call juDFT_error('Error in allocating oper_temp',calledby='wann_convert_fleur_w90')

        do ik = 1, num_kpts

          counter = 0
          do m = 1, write_bands
            do n = 1, m
              counter = counter + 1
              do compo=1,num_compos
                oper_temp(compo, counter)=oper_o(n, m, ik, map3(compo))
!              oper_o(n, m, ik, compo) = oper_temp(compo, counter)
!              oper_o(m, n, ik, compo) = conjg(oper_temp(compo, counter))
              enddo !compo
            end do !n
          end do !m

          oper_temp=oper_temp*conversionfactor

          write (spn_in) ((oper_temp(s, m), s=1, num_compos), m=1, (write_bands*(write_bands + 1))/2)
        end do !ik

      else
        write (spn_in, *) header
        write (spn_in, *) nb_tmp, nkp_tmp

        do ik = 1, num_kpts
          do m = 1, write_bands
            do n = 1, m
              do compo=1,num_compos
                  s_real=real(oper_o(n, m, ik, map3(compo)))*conversionfactor
                  s_img=imag(oper_o(n, m, ik, map3(compo)))*conversionfactor
                  write (spn_in, *) s_real, s_img   !,compo,n, m, ik 
!               oper_o(n, m, ik, compo) = cmplx(s_real, s_img, dp)
              ! Read upper-triangular part, now build the rest
!              oper_o(m, n, ik, compo) = conjg(oper_o(n, m, ik, compo))
              enddo 
            end do
          end do
        enddo

      endif  
      close(spn_in)

      if(wann%l_mmn0_to_spn2)then
         open(spn_in,file="WF2.spn",form=trim(writeform))
         
         write (spn_in, *) header
         write (spn_in, *) nb_tmp, nkp_tmp

         do ik = 1, num_kpts
          do m = 1+write_bands,2*write_bands
            do n = 1+write_bands, m
              do compo=1,num_compos
                  s_real=real(oper_o(n, m, ik, compo))*conversionfactor
                  s_img=imag(oper_o(n, m, ik, compo))*conversionfactor
                  write (spn_in, *) s_real, s_img
!               oper_o(n, m, ik, compo) = cmplx(s_real, s_img, dp)
              ! Read upper-triangular part, now build the rest
!              oper_o(m, n, ik, compo) = conjg(oper_o(n, m, ik, compo))
              enddo 
            end do
          end do
         enddo

         close(spn_in)
      endif !  wann%l_mmn0_to_spn2

      end subroutine wann_convert_fleur_w90
      END MODULE m_wann_convert_fleur_w90