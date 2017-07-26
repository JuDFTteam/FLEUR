      module m_wann_get_qpts
      USE m_fleurenv
      contains
      subroutine wann_get_qpts(
     >               l_bzsym,film,l_onedimens,l_readqpts,
     <               nqpts,qpoints,param_file)
c********************************************************
c     Read in the q-points from qpts file.
c 
c   
c********************************************************
      implicit none
      logical,intent(in)  :: l_bzsym,film
      logical,intent(in)  :: l_onedimens,l_readqpts
      integer,intent(out) :: nqpts
      real,intent(inout)  :: qpoints(:,:)
      character(len=20),intent(in) :: param_file

      real             :: scale
      !integer          :: at,j
      integer          :: iter!,len,num_wann,num_bands,nn,i
      logical          :: l_file

      if(l_bzsym)then
         inquire(file='w90qpts',exist=l_file)
         IF(.NOT.l_file) CALL fleur_err("where is w90qpts?",calledby
     +        ="wann_get_qpts")
         open(987,file='w90qpts',status='old',form='formatted')
         read(987,*)nqpts, scale
         write(6,*)"wann_get_qpts: nqpts=",nqpts
         if(l_readqpts)then
            IF(SIZE(qpoints,1)/=3) CALL fleur_err("wann_get_qpts: 1"
     +           ,calledby ="wann_get_qpts")
            IF(SIZE(qpoints,2)/=nqpts) CALL fleur_err("wann_get_qpts: 2"
     +           ,calledby ="wann_get_qpts")
            do iter=1,nqpts
               read(987,*)qpoints(:,iter)
            enddo
         endif   
      else
         inquire(file=param_file,exist=l_file)
         IF(.NOT.l_file) CALL fleur_err(
     >         "where is "//trim(param_file)//"?",calledby
     +        ="wann_get_qpts")
         open(987,file=param_file,status='old',form='formatted')
         read(987,*)nqpts,scale
         write(6,*)"wann_get_qpts: nqpts=",nqpts
         if(l_readqpts)then
            IF(SIZE(qpoints,1)/=3) CALL fleur_err("wann_get_qpts: 1"
     +           ,calledby ="wann_get_qpts")
            IF(SIZE(qpoints,2)/=nqpts) CALL fleur_err("wann_get_qpts: 2"
     +           ,calledby ="wann_get_qpts")
            do iter=1,nqpts
               read(987,*)qpoints(:,iter)
            enddo
         endif   
      endif

      close(987)

      if(l_readqpts)then
         qpoints=qpoints/scale !* 2.0 !2xBZ
         if(film.and..not.l_onedimens)then 
            qpoints(3,:)=0.0
         endif   
         do iter=1,nqpts
            write(6,*)qpoints(:,iter)
         enddo
      endif

      end subroutine wann_get_qpts
      end module m_wann_get_qpts
