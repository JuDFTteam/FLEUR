!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

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
      USE m_inpgen_help
      USE m_constants
      IMPLICIT NONE
    
      INTEGER natmax,nop48,nline,natin,ngen,i,j,bfh
      INTEGER nops,no3,no2,na,numSpecies,i_c,element
      INTEGER infh,errfh,warnfh,symfh,dbgfh,outfh,dispfh
      LOGICAL cal_symm,checkinp,newSpecies
      LOGICAL cartesian,oldfleur,l_hyb  ,inistop
      REAL    aa
 
      REAL a1(3),a2(3),a3(3),scale(3),factor(3)
      INTEGER              :: elementNumSpecies(0:104)
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
      CHARACTER(LEN=8) :: tempNumberString
      CHARACTER(len=20), ALLOCATABLE :: atomLabel(:)

      TYPE(t_input)    :: input
      TYPE(t_atoms)    :: atoms
      TYPE(t_cell)     :: cell
      TYPE(t_sym)      :: sym
      TYPE(t_noco)     :: noco
      TYPE(t_vacuum)   :: vacuum

      !Start program and greet user
      CALL inpgen_help()

      INQUIRE(file='inp.xml',exist=l_inpxml)
      IF (l_inpxml.AND..NOT.judft_was_argument("-inp.xml")) CALL judft_error("inp.xml exists and can not be overwritten")
      
      IF (judft_was_argument("-inp")) THEN
         STOP "not yet"
         !CALL read_old_input()
         full_input=.TRUE.
      ELSEIF (judft_was_argument("-inp.xml")) THEN
         STOP "not yet"
         !CALL r_inpXML()
         full_input=.TRUE.
      ELSEIF(judft_was_argument("-f")) THEN
         !read the input
         
         CALL read_inpgen_input(atom_pos,atom_id,atom_label,amat,


         film,symor,atomid,atompos,atomlabel,amat,dvac,noco)
         full_input=.FALSE.
      ELSE
         CALL judft_error("You should either specify -inp,-inp.xml or -f command line options. Check -h if unsure")
      ENDIF
      IF (.NOT.full_input) THEN
         !First we determine the spacegoup and map the atoms to groups
         CALL make_crystal(film,symor,atomid,atompos,atomlabel,amat,dvac,noco,&
              cell,sym,atoms)          
         
         !All atom related parameters are set here. Note that some parameters might
         !have been set in the read_input call before by adding defaults to the atompar module
         CALL make_atomic_defaults(input,vacuum,cell,oneD,atoms)
         
         !Set all defaults that have not been specified before or can not be specified in inpgen
         CALL make_defaults(atoms,vacuum,input,stars,sliceplot,forcetheo,banddos,&
              cell,sym,xcpot,noco,oneD,hybrid,kpts)
      ENDIF
      !
      ! k-points can also be modified here
      !
      call make_kpoints()
      !
      !Now the IO-section
      !
      IF (.NOT.l_inpxml) THEN
         !the inp.xml file
         l_explicit=judft_was_argument("-explicit")
         CALL dump_FleurInputSchema()
         CALL w_inpxml(&
              atoms,vacuum,input,stars,sliceplot,forcetheo,banddos,&
              cell,sym,xcpot,noco,oneD,hybrid,kpts,&
              div,l_gamma,& !should be in kpts!?
              namex,relcor,dtild_opt,name_opt,&!?should be somewhere...
              .false.,"inp.xml",&
              l_explicit,enpara)
         !the sym.xml file
         CALL sym%writeXML(0,"sym.xml")
      ENDIF
      CALL kpts%writeXML(0,"kpts.xml")
      
      ! Structure in  xsf-format
      OPEN (55,file="struct.xsf")
      CALL xsf_WRITE_atoms(55,atoms,input%film,.FALSE.,cell%amat)
      CLOSE (55)

      
      CALL juDFT_end("All done",1)

    END PROGRAM inpgen
