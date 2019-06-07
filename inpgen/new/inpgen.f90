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
      
      CALL inpgen_help()


      CALL read_input()

      CALL make_crystal(film, symor,atomid,atompos,atomlabel,amat,dvac,noco,&
           cell,sym,atoms)

          

      CALL make_atomic_defaults(cell,atoms)

      !Now the IO-section
      
      CALL w_inpxml()
      CALL write_sym()
!
! --> Structure in  xsf-format
!
      OPEN (55,file="struct.xsf")
      CALL xsf_WRITE_atoms(55,atoms,input%film,.FALSE.,cell%amat)
      CLOSE (55)

      
      CALL juDFT_end("All done",1)

    END PROGRAM inpgen
