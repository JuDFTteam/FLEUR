PROGRAM inpgen
!----------------------------------------------------------------------------+
!   Set up a FLEUR inp-file from basic input data; for use and docu please   !
!   refer to inpgen.html (or see http://www.flapw.de/docs/inpgen.html)       !
!                                                                            |
!   The program is based on the input conventions of the FLAIR-code, so that !
!   some compatibility is ensured. The symmetry generator was written by     ! 
!   M.Weinert and implemented in the FLAIR-code by G.Schneider.              !
!                                                                    gb`02   |
!----------------------------------------------------------------------------+
      use m_juDFT
      USE m_structinput
      USE m_crystal
      USE m_socorssdw
      USE m_rwsymfile
      USE m_setinp
      USE m_writestruct
      USE m_xsf_io, ONLY : xsf_write_atoms
      USE m_types
      IMPLICIT NONE
    
      INTEGER natmax,nop48,nline,natin,ngen,i,j
      INTEGER nops,no3,no2,na,numSpecies
      INTEGER infh,errfh,bfh,warnfh,symfh,dbgfh,outfh,dispfh
      LOGICAL cal_symm,checkinp,newSpecies,noangles
      LOGICAL cartesian,oldfleur,l_hyb  ,inistop
      REAL    aa
 
      REAL a1(3),a2(3),a3(3),scale(3),factor(3)
      INTEGER, ALLOCATABLE :: mmrot(:,:,:)
      REAL,    ALLOCATABLE :: ttr(:, :),atompos(:, :),atomid(:) 
      REAL,    ALLOCATABLE :: idlist(:)
      INTEGER, ALLOCATABLE ::  ntyrep(:)              ! these variables are allocated with
      INTEGER, ALLOCATABLE :: natype(:),natrep(:),natmap(:)  ! or  'nat'
      INTEGER, ALLOCATABLE :: speciesRepAtomType(:),atomTypeSpecies(:)
     
      INTEGER, PARAMETER :: xl_buffer=16384              ! maximum length of read record
      CHARACTER(len=xl_buffer) :: buffer

      CHARACTER(len=80):: title
      CHARACTER(len=7) :: symfn
      CHARACTER(len=4) :: dispfn

            TYPE(t_input)    :: input
          TYPE(t_atoms)    :: atoms
          TYPE(t_cell)     :: cell
          TYPE(t_sym)      :: sym
          TYPE(t_noco)     :: noco
          TYPE(t_vacuum)   :: vacuum
       
      nop48 = 48
      natmax = 9999
      ngen = 0
      infh = 5
      errfh = 6 ; warnfh = 6 ; dbgfh = 6 ; outfh = 6
      bfh = 93
      symfh = 94
      symfn = 'sym    '
      dispfh = 97
      dispfn='disp'
      nline = 0

      input%l_inpXML = .FALSE.

      ALLOCATE ( mmrot(3,3,nop48), ttr(3,nop48) )
      ALLOCATE ( atompos(3,natmax),atomid(natmax) )

!      OPEN (5,file='inp2',form='formatted',status='old')
      OPEN (6,file='out',form='formatted',status='unknown')

      CALL struct_input(&
     &                  infh,errfh,bfh,warnfh,symfh,symfn,&
     &                  natmax,nop48,&
     &                  nline,xl_buffer,buffer,&
     &                  title,input%film,cal_symm,checkinp,sym%symor,&
     &                  cartesian,oldfleur,a1,a2,a3,vacuum%dvac,aa,scale,noangles,&
     &                 factor,natin,atomid,atompos,ngen,mmrot,ttr,&
     &                  l_hyb,noco%l_soc,noco%l_ss,noco%theta,noco%phi,noco%qss,inistop)!keep

!      CLOSE (5)

      IF (.not.input%film) vacuum%dvac=a3(3)
      WRITE (6,*)
      WRITE (6,*) title
      WRITE (6,*) 'film=',input%film,'cartesian=',cartesian
      WRITE (6,*) 'checkinp=',checkinp,'symor=',sym%symor
      WRITE (6,*)
      WRITE (6,'(a5,3f10.5)') 'a1 = ',a1(:)
      WRITE (6,'(a5,3f10.5)') 'a2 = ',a2(:)
      WRITE (6,'(a5,3f10.5)') 'a3 = ',a3(:)
      WRITE (6,*)
      WRITE (6,'(2(a5,f10.5))') 'dvac=',vacuum%dvac,' aa =',aa
      WRITE (6,'(a8,3f10.5)') 'scale = ',scale(:)
      WRITE (6,*)
      WRITE (6,'(a6,i3,a6,999i5)') 'natin=',natin,' Z = ',&
     &                             (nint(atomid(i)),i=1,abs(natin))
      WRITE (6,*) 'positions: '
      WRITE (6,'(3(3x,f10.5))') ((atompos(j,i),j=1,3),i=1,abs(natin))
      WRITE (6,*)
      WRITE (6,*) 'generators: ',ngen,'(excluding identity)'
      DO i = 2, ngen+1
         WRITE (6,*) i
         WRITE (6,'(3i5,f8.3)') (mmrot(1,j,i),j=1,3),ttr(1,i)
         WRITE (6,'(3i5,f8.3)') (mmrot(2,j,i),j=1,3),ttr(2,i)
         WRITE (6,'(3i5,f8.3)') (mmrot(3,j,i),j=1,3),ttr(3,i)
      ENDDO
      IF (noco%l_soc) WRITE(6,'(a4,2f10.5)') 'soc:',noco%theta,noco%phi
      IF (noco%l_ss)  WRITE(6,'(a4,3f10.5)') 'qss:',noco%qss(:)
!
! --> generate symmetry from input (atomic positions, generators or whatever)
!     
      CALL crystal(&
     &             dbgfh,errfh,outfh,dispfh,dispfn,&
     &             cal_symm,cartesian,sym%symor,input%film,&
     &             natin,natmax,nop48,&
     &             atomid,atompos,a1,a2,a3,aa,scale,noangles,&
     &             sym%invs,sym%zrfs,sym%invs2,sym%nop,sym%nop2,&
     &             ngen,mmrot,ttr,atoms%ntype,atoms%nat,nops,&
     &             atoms%neq,ntyrep,atoms%zatom,natype,natrep,natmap,&
     &             sym%mrot,sym%tau,atoms%pos,cell%amat,cell%bmat,cell%omtil)

      IF (noco%l_ss.OR.noco%l_soc)  THEN
         CALL soc_or_ssdw(&
     &                    noco%l_soc,noco%l_ss,noco%theta,noco%phi,noco%qss,cell%amat,&
     &                    sym%mrot,sym%tau,sym%nop,sym%nop2,atoms%nat,atomid,atompos,&
     &                    mmrot,ttr,no3,no2,atoms%ntype,atoms%neq,natmap,&
     &                    ntyrep,natype,natrep,atoms%zatom,atoms%pos)
         sym%nop = no3 ; sym%nop2 = no2
         sym%mrot(:,:,1:sym%nop) = mmrot(:,:,1:sym%nop)
         sym%tau(:,1:sym%nop) = ttr(:,1:sym%nop)
      ENDIF
      DEALLOCATE ( mmrot, ttr, atompos )

      ALLOCATE ( atoms%taual(3,atoms%nat),idlist(atoms%ntype) ) 
      WRITE (6,*)
      WRITE (6,'(a6,i3,a6,i3)') 'atoms%ntype=',atoms%ntype,' atoms%nat= ',atoms%nat
      na = 0
      DO i = 1, atoms%ntype
        WRITE (6,'(a3,i3,a2,i3,a6,i3)') ' Z(',i,')=',nint(atoms%zatom(i)),&
     &                                             ' atoms%neq= ',atoms%neq(i)
        DO j = 1, atoms%neq(i)
           WRITE (6,'(3f10.6,10x,i7)')&
     &           atoms%pos(:,natmap(na+j)),natmap(na+j)
           atoms%taual(:,na+j) = atoms%pos(:,natmap(na+j))      ! reorder coordinates
           idlist(i)    = atomid(natmap(na+j))      !     and atomic id's
        ENDDO
        na = na + atoms%neq(i)
      ENDDO
      DO i=1,atoms%nat
        atoms%pos(:,i) = matmul( cell%amat , atoms%taual(:,i) )
      ENDDO

!
! --> write a file 'sym.out' with accepted symmetry operations
!
      nops = sym%nop
      symfn = 'sym.out'
      IF (.not.input%film) sym%nop2=sym%nop
      CALL rw_symfile(&
     &                'W',symfh,symfn,nops,cell%bmat,&
     &                 sym%mrot,sym%tau,sym%nop,sym%nop2,sym%symor)

      ALLOCATE (atomTypeSpecies(atoms%ntype))
      ALLOCATE (speciesRepAtomType(atoms%nat))
      numSpecies = 0
      speciesRepAtomType = -1
      atomTypeSpecies = -1
      DO i = 1, atoms%nat
         newSpecies = .TRUE.
         DO j = 1, i-1
            IF(atomid(i).EQ.atomid(j)) THEN
               newSpecies = .FALSE.
               atomTypeSpecies(natype(i)) = atomTypeSpecies(natype(j))
               EXIT
            END IF
         END DO
         IF(newSpecies) THEN
            numSpecies = numSpecies + 1
            speciesRepAtomType(numSpecies) = natype(i)
            atomTypeSpecies(natype(i)) = numSpecies
         END IF
      END DO

!
! --> set defaults for FLEUR inp-file
!
      ALLOCATE ( atoms%rmt(atoms%ntype) )
      atoms%nlod=9  ! This fixed dimensioning might have to be made more dynamical!
      CALL set_inp(&
     &             infh,nline,xl_buffer,buffer,l_hyb,&
     &             atoms,sym,cell,title,idlist,&
     &             input,vacuum,noco,&
     &             atomTypeSpecies,speciesRepAtomType,&
     &             a1,a2,a3)

      DEALLOCATE (atomTypeSpecies,speciesRepAtomType)
      DEALLOCATE ( ntyrep, natype, natrep, atomid )

!
! --> Structure in povray or xsf-format
!
      IF (.false.) THEN
         CALL write_struct(&
     &                  atoms%ntype,atoms%nat,atoms%neq,&
     &                  atoms%rmt,atoms%pos,natmap,cell%amat)!keep
      ELSE 
         OPEN (55,file="struct.xsf")
         CALL xsf_WRITE_atoms(&
     &                        55,atoms,input%film,.false.,cell%amat)
         CLOSE (55)
      ENDIF

      DEALLOCATE (vacuum%izlay)
      DEALLOCATE ( atoms%taual,sym%mrot,sym%tau,atoms%neq,atoms%zatom,atoms%rmt,natmap,atoms%pos,idlist )

      IF (inistop)  CALL juDFT_end("Symmetry done")

      END 
