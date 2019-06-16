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
  USE m_juDFT
  USE m_inpgen_help
  USE m_read_inpgen_input
  USE m_make_crystal
  USE m_make_atomic_defaults
  USE m_make_defaults
  USE m_make_kpoints
  USE m_winpxml
  USE m_xsf_io
  USE m_types_input
  USE m_types_atoms
  USE m_types_cell
  USE m_types_sym
  USE m_types_noco
  USE m_types_vacuum
  USE m_types_banddos
  USE m_types_hybrid
  USE m_types_xcpot_inbuild_nofunction
  USE m_types_forcetheo
  USE m_types_kpts
  USE m_types_enpara
  USE m_types_oneD
  USE m_types_sliceplot
  USE m_types_stars
  
      IMPLICIT NONE
    
      REAL,    ALLOCATABLE :: atompos(:, :),atomid(:) 
      CHARACTER(len=20), ALLOCATABLE :: atomLabel(:)
      LOGICAL               :: l_fullinput,l_explicit,l_inpxml,l_include(4)
      
      TYPE(t_input)    :: input
      TYPE(t_atoms)    :: atoms
      TYPE(t_cell)     :: cell
      TYPE(t_sym)      :: sym
      TYPE(t_noco)     :: noco
      TYPE(t_vacuum)   :: vacuum
      TYPE(t_banddos)  :: banddos
      TYPE(t_hybrid)   :: hybrid
      TYPE(t_xcpot_inbuild_nf)::xcpot
      TYPE(t_enpara)   :: enpara
      TYPE(t_forcetheo):: forcetheo
      TYPE(t_kpts)     :: kpts
      TYPE(t_oned)     :: oned
      TYPE(t_sliceplot):: sliceplot
      TYPE(t_stars)    :: stars
      
      CHARACTER(len=40):: kpts_str
      kpts_str=""
      
      !Start program and greet user
      CALL inpgen_help()
      l_explicit=judft_was_argument("-explicit")
      
      OPEN(6,file='out')

      INQUIRE(file='inp.xml',exist=l_inpxml)
      IF (l_inpxml.AND..NOT.(judft_was_argument("-inp.xml").or.judft_was_argument("-overwrite")))&
           CALL judft_error("inp.xml exists and can not be overwritten")
      
      IF (judft_was_argument("-inp")) THEN
         STOP "not yet"
         !CALL read_old_input()
         l_fullinput=.TRUE.
      ELSEIF (judft_was_argument("-inp.xml")) THEN
         call read_old_inp(input,DIMENSION,atoms,sphhar,cell,stars,sym,noco,vacuum,forcetheo,&
              sliceplot,banddos,enpara,xcpot,kpts,hybrid, oneD)
         l_fullinput=.TRUE.
      ELSEIF(judft_was_argument("-f")) THEN
         !read the input
         
         CALL read_inpgen_input(atompos,atomid,atomlabel,kpts_str,&
              input,sym,noco,vacuum,stars,xcpot,cell,hybrid)
         l_fullinput=.FALSE.
      ELSE
         CALL judft_error("You should either specify -inp,-inp.xml or -f command line options. Check -h if unsure")
      ENDIF
      IF (.NOT.l_fullinput) THEN
         !First we determine the spacegoup and map the atoms to groups
         CALL make_crystal(input%film,atomid,atompos,atomlabel,vacuum%dvac,noco,&
              cell,sym,atoms)          
         
         !All atom related parameters are set here. Note that some parameters might
         !have been set in the read_input call before by adding defaults to the atompar module
         CALL make_atomic_defaults(input,vacuum,cell,oneD,atoms,enpara)
         
         !Set all defaults that have not been specified before or can not be specified in inpgen
         CALL make_defaults(atoms,sym,cell,vacuum,input,stars,&
                   xcpot,noco,hybrid)
      ENDIF
      !
      ! k-points can also be modified here
      !
      call make_kpoints(kpts,cell,sym,input%film,noco%l_ss.or.noco%l_soc,kpts_str)
      !
      !Now the IO-section
      !
      IF (.NOT.l_inpxml.or.judft_was_argument("-overwrite")) THEN
         call determine_includes(l_include)
         !the inp.xml file
         !CALL dump_FleurInputSchema()
         CALL w_inpxml(&
              atoms,vacuum,input,stars,sliceplot,forcetheo,banddos,&
              cell,sym,xcpot,noco,oneD,hybrid,kpts,enpara,&
              5,l_explicit,l_include,"inp.xml")
         if (.not.l_include) CALL sym%print_XML(99,"sym.xml")
      ENDIF
      if (.not.l_include) CALL kpts%print_XML(99,"kpts.xml")
      
      ! Structure in  xsf-format
      OPEN (55,file="struct.xsf")
      CALL xsf_WRITE_atoms(55,atoms,input%film,.FALSE.,cell%amat)
      CLOSE (55)
      CLOSE(6)
      
      CALL juDFT_end("All done")

    contains
      subroutine determine_includes(l_include)
        logical,intent(out)::l_include(4)  !kpts,operations,species,position

        l_include=[.false.,.false.,.true.,.true.]

        IF (judft_was_argument("-inc")) THEN
           str=judft_string_for_argument("-inc")
           
           do while(len_trim(str)>0) then
              if (str(1:1)=='-') then
                 incl=.false.
                 str=str(2:)
              else
                 incl=.true.
                 if (str(1:1)=='+') str=str(2:)
              endif
              select case(str(1:1))
              case ('k','K')
                 l_include(1)=incl
              case ('o','O')
                 l_include(2)=incl
              case ('s','S')
                 l_include(3)=incl
              case ('p','P')
                 l_include(4)=incl
              case ('a','A')
                 l_include(:)=incl
              end select
              if (index(str,"'")>0) then
                 str=str(index(str,"'")+1:)
              else
                 str=""
              end if
           end do
        endif
        
        IF (LEN_TRIM(str)>1) CALL judft_error("Do not specify k-points in file and on command line")
        str=judft_string_for_argument("-k")
     END IF
      
    END PROGRAM inpgen
