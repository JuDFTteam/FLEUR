      module m_wann_mmnk_symm
      use m_juDFT
      private
      public:: wann_mmnk_symm
      contains
c******************************************************************
c     Find out minimal set of k-point-pairs that have to be
c     calculated; map symmetry-related k-point-pairs to this
c     minimal set.
c     Frank Freimuth
c******************************************************************      
      subroutine wann_mmnk_symm(
     >               fullnkpts,nntot,bpt,gb,l_bzsym,
     >               irreduc,mapkoper,l_p0,film,nop,
     >               invtab,mrot,l_onedimens,tau,
     <               pair_to_do,maptopair,kdiff)

      implicit none
      integer,intent(in) :: nop
      integer,intent(in) :: mrot(3,3,nop)
      integer,intent(in) :: invtab(nop)
      integer,intent(in) :: fullnkpts
      integer,intent(in) :: nntot
      integer,intent(in) :: bpt(nntot,fullnkpts)
      integer,intent(in) :: gb(3,nntot,fullnkpts)
      logical,intent(in) :: l_bzsym
      integer,intent(in) :: irreduc(fullnkpts)
      integer,intent(in) :: mapkoper(fullnkpts)
      logical,intent(in) :: l_p0,film
      real,intent(in)    :: tau(3,nop)

      integer,intent(out):: pair_to_do(fullnkpts,nntot)
      integer,intent(out):: maptopair(3,fullnkpts,nntot)
      real,intent(out)   :: kdiff(3,nntot)

      integer :: ikpt,ikpt_k,kptibz,ikpt_b
      integer :: num_pair,kptibz_b
      integer :: index(fullnkpts,nntot)
      integer :: oper,oper_b,num_rot,num_conj,repo,repo_b
      integer :: k,kx,ky,repkpt,repkpt_b,repkpt_bb
      integer :: sign,sign_b,ngis,ngis_b
      logical :: startloop
      logical :: l_file,l_onedimens
      real    :: kpoints(3,fullnkpts)
      real    :: kdiffvec(3)
      integer :: multtab(nop,nop)
      integer :: fullnkpts_tmp,kr
      real    :: scale
      real    :: brot(3)
      logical :: l_testnosymm,l_nosymm1,l_nosymm2

      index(:,:)=0
      num_pair=0
      num_rot=0
      num_conj=0
      pair_to_do(:,:)=0

      inquire(file='testnosymm',exist=l_testnosymm)

c-----Test for nonsymmorphic space groups
      if(l_bzsym)then
       do ikpt=1,fullnkpts
         oper=abs(mapkoper(ikpt))
         if( any( abs(tau(:,oper)).gt.1.e-6 ) ) l_testnosymm=.true.
       enddo
      endif

      if(l_bzsym)then
         call close_pt(nop,mrot,multtab)
      endif

      do 10 ikpt = 1,fullnkpts  ! loop by k-points starts
        l_nosymm1=.false. 
        kptibz=ikpt 
        if(l_bzsym) then
           kptibz=irreduc(ikpt)
           oper=mapkoper(ikpt)
           if(oper.lt.0)then
              oper=-oper
              sign=-1
           else
              sign=1
           endif
           if( 
     &          any( abs(tau(:,oper)).gt.1.e-6 ) 
     &         )l_nosymm1=.true.
        endif

        do 15 ikpt_b = 1,nntot
         l_nosymm2=.false.  
         if(index(ikpt,ikpt_b).eq.1)cycle
         kptibz_b=bpt(ikpt_b,ikpt)
         if(l_bzsym) then
            oper_b=mapkoper(kptibz_b)
            kptibz_b=irreduc(kptibz_b)
            if(oper_b.lt.0)then
               oper_b=-oper_b
               sign_b=-1
            else
               sign_b=1
            endif
            if( 
     &          any( abs(tau(:,oper_b)).gt.1.e-6 ) 
     &         )l_nosymm2=.true.
         endif

         if(l_testnosymm)goto 33
         if(l_nosymm1.or.l_nosymm2)goto 33

c***************************************************************
c..the conjugation selection rule, which speeds up significantly
c***************************************************************
        do ikpt_k = 1,nntot
         if((bpt(ikpt_k,bpt(ikpt_b,ikpt)).eq.ikpt).and.
     &       gb(1,ikpt_b,ikpt).eq.(-gb(1,ikpt_k,bpt(ikpt_b,ikpt))).and.
     &       gb(2,ikpt_b,ikpt).eq.(-gb(2,ikpt_k,bpt(ikpt_b,ikpt))).and.
     &       gb(3,ikpt_b,ikpt).eq.(-gb(3,ikpt_k,bpt(ikpt_b,ikpt))).and.
     &       (index(bpt(ikpt_b,ikpt),ikpt_k).eq.1))then
            index(ikpt,ikpt_b)=1
            maptopair(1,ikpt,ikpt_b)=bpt(ikpt_b,ikpt)
            maptopair(2,ikpt,ikpt_b)=ikpt_k
            maptopair(3,ikpt,ikpt_b)=1
c            print*,"conjugation"
            num_conj=num_conj+1
           goto 15  
         endif
        enddo !ikpt_k
c****************************************************************
c     check whether k-point pairs can be mapped onto each other
c         by rotation
c****************************************************************
        if(l_bzsym)then
c         if(all(gb(:,ikpt_b,ikpt).eq.0))then
          do k=1,fullnkpts
           if(irreduc(k).eq.kptibz)then
             repkpt=k
             repo=mapkoper(k)
             if(repo.lt.0)then
                repo=-repo
                ngis=-1
             else
                ngis=1
             endif

             do kx=1,fullnkpts
              if(irreduc(kx).eq.kptibz_b)then
               repkpt_bb=kx
               repo_b=mapkoper(kx)
               if(repo_b.lt.0)then
                  repo_b=-repo_b
                  ngis_b=-1
               else
                  ngis_b=1
               endif 
               do ky=1,nntot
                if(bpt(ky,repkpt).eq.repkpt_bb)then
                 repkpt_b=ky
                 if (index(repkpt,repkpt_b).eq.1)then
                  if(.not.all(gb(:,ikpt_b,ikpt).eq.0))then
                   brot(:)=0.0
                   do kr=1,3
                     brot(:)=brot(:)
     &                 -sign_b*mrot(kr,:,oper)*gb(kr,ikpt_b,ikpt)
     &                 +ngis_b*mrot(kr,:,oper_b)*gb(kr,repkpt_b,repkpt)
                   enddo
                   if( any(   abs(brot).gt.1e-6       )   )cycle
                  endif
                  if(sign*ngis*multtab(invtab(oper),repo).eq.
     &               sign_b*ngis_b*multtab(invtab(oper_b),repo_b))then  
                    maptopair(1,ikpt,ikpt_b)=repkpt
                    maptopair(2,ikpt,ikpt_b)=repkpt_b
                    maptopair(3,ikpt,ikpt_b)=2+(1-ngis*sign)/2
                    index(ikpt,ikpt_b)=1
                    num_rot=num_rot+1
                    goto 15
                  endif
                 endif                   
                endif   
               enddo    
              endif   
             enddo 
           endif   
          enddo   
c        endif !gb=0   
        endif

 33     continue

        index(ikpt,ikpt_b)=1
        num_pair=num_pair+1
        pair_to_do(ikpt,ikpt_b)=num_pair
        
15    continue !loop over nearest neighbor k-points
10    continue ! end of cycle by the k-points  

      if(l_p0)then
      write(6,*)"pairs to calculate: ",num_pair
      write(6,*)"maps by conjugation: ",num_conj
      write(6,*)"maps by rotation:", num_rot
      write(6,*)"num_pair+num_rot+num_conj:",num_pair+num_conj+num_rot
      write(6,*)"fullnkpts*nntot:", fullnkpts*nntot
      endif !l_p0

c*****************************************************************
c     determine difference vectors that occur on the k-mesh
c*****************************************************************      
      if (l_bzsym) then
         l_file=.false.
         inquire(file='w90kpts',exist=l_file)
         IF(.NOT.l_file)  CALL juDFT_error
     +        ("w90kpts not found, needed if bzsym",calledby
     +        ="wann_mmnk_symm")
         open(412,file='w90kpts',form='formatted')
         read(412,*)fullnkpts_tmp,scale
         do k=1,fullnkpts
               read(412,*)kpoints(:,k)
         enddo   
         kpoints=kpoints/scale
         close(412)
      else   
            open(412,file='kpts',form='formatted')
            read(412,*)fullnkpts_tmp,scale
            do k=1,fullnkpts
               read(412,*)kpoints(:,k)
            enddo   
            kpoints(:,:)=kpoints/scale
            if (film.and..not.l_onedimens) kpoints(3,:)=0.0
            close(412)
      endif
      if(l_p0)then
         print*,"vectors combining nearest neighbor k-points:"
      endif   
      ky=1
      do k=1,fullnkpts
         do kx=1,nntot
            kdiffvec=kpoints(:,bpt(kx,k))+gb(:,kx,k)-kpoints(:,k)
            do ikpt=1,ky-1
               if(all(abs(kdiff(:,ikpt)-kdiffvec).le.0.0001))goto 200
            enddo 
            IF(ky>nntot)  CALL juDFT_error("problem in wann_mmnk_symm"
     +           ,calledby ="wann_mmnk_symm")
            kdiff(:,ky)=kdiffvec(:)
            if(l_p0)then
               print*,ky,k,kx,kdiff(:,ky)
            endif

            ky=ky+1
 200        continue

         enddo
      enddo   
      end subroutine
      SUBROUTINE close_pt(
     >                    nops,mrot,
     <                    mtable)

      IMPLICIT NONE

      INTEGER, INTENT (IN)  :: nops,mrot(3,3,nops)
      INTEGER, INTENT (OUT) :: mtable(nops,nops)   ! table(i,j) = {R_i|0}{R_j|0}

      INTEGER              :: i,j,k,mp(3,3),map(nops)

!---> loop over all operations
      DO j=1,nops

         map(1:nops) = 0

!--->    multiply {R_j|0}{R_i|0}
         DO i=1,nops
            mp = matmul( mrot(:,:,j) , mrot(:,:,i) )

!--->       determine which operation this is
            DO k = 1, nops
              IF ( all( mp(:,:)==mrot(:,:,k) ) ) THEN
                 IF ( map(i) .eq. 0 ) THEN
                    map(i) = k
                 ELSE
                    WRITE (6,'(" Symmetry error : multiple ops")')
                    CALL juDFT_error("close_pt: Multiple ops (Bravais)"
     +                   ,calledby ="wann_mmnk_symm")
                 ENDIF
              ENDIF
            ENDDO

            IF (map(i).eq.0) THEN
               WRITE (6,'(" Group not closed (Bravais lattice)")')
               WRITE (6,'(" operation j=",i2,"  map=",12i4,:/,
     &                  (21x,12i4))')  j, map(1:nops)
               CALL juDFT_error("close_pt: Not closed",calledby
     +              ="wann_mmnk_symm")
            ENDIF
         ENDDo
         mtable(j,1:nops) = map(1:nops)
      ENDDO

      END SUBROUTINE close_pt

      end module m_wann_mmnk_symm
