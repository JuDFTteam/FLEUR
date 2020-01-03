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
  USE m_types_hybinp
  USE m_types_xcpot_inbuild_nofunction
  USE m_types_forcetheo
  USE m_types_kpts
  USE m_types_enpara
  USE m_types_oneD
  USE m_types_sliceplot
  USE m_types_stars
  use m_read_old_inp
  use m_fleurinput_read_xml
  USE m_types_mpinp

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
      TYPE(t_mpinp)    :: mpinp
      TYPE(t_hybinp)   :: hybinp
      TYPE(t_xcpot_inbuild_nf)::xcpot
      TYPE(t_enpara)   :: enpara
      TYPE(t_forcetheo):: forcetheo
      TYPE(t_kpts)     :: kpts
      TYPE(t_oned)     :: oned
      TYPE(t_sliceplot):: sliceplot
      TYPE(t_stars)    :: stars

      CHARACTER(len=40):: kpts_str
      LOGICAL          :: l_exist
      INTEGER          :: idum

      INTERFACE
       FUNCTION dropDefaultEConfig() BIND(C, name="dropDefaultEconfig")
         USE iso_c_binding
         INTEGER(c_int) dropDefaultEConfig
       END FUNCTION dropDefaultEConfig
      END INTERFACE

      kpts_str=""

      !Start program and greet user
      CALL inpgen_help()
      l_explicit=judft_was_argument("-explicit")

      INQUIRE(file='default.econfig',exist=l_exist)
      IF (.NOT.l_exist) idum=dropDefaultEconfig()

      OPEN(6,file='out')

      INQUIRE(file='inp.xml',exist=l_inpxml)
      IF (l_inpxml.AND..NOT.(judft_was_argument("-inp.xml").or.judft_was_argument("-overwrite")))&
           CALL judft_error("inp.xml exists and can not be overwritten")

      IF (judft_was_argument("-inp")) THEN
         call read_old_inp(input,atoms,cell,stars,sym,noco,vacuum,forcetheo,&
              sliceplot,banddos,enpara,xcpot,kpts,hybinp, oneD)
         l_fullinput=.TRUE.
      ELSEIF (judft_was_argument("-inp.xml")) THEN
         !not yet
         call Fleurinput_read_xml(cell,sym,atoms,input,noco,vacuum,&
         sliceplot=Sliceplot,banddos=Banddos,hybinp=hybinp,oned=Oned,xcpot=Xcpot,kpts=Kpts)
         Call Cell%Init(Dot_product(Atoms%Volmts(:),Atoms%Neq(:)))
         call atoms%init(cell)
         Call Sym%Init(Cell,Input%Film)
         l_fullinput=.TRUE.
      ELSEIF(judft_was_argument("-f")) THEN
         !read the input

         CALL read_inpgen_input(atompos,atomid,atomlabel,kpts_str,&
              input,sym,noco,vacuum,stars,xcpot,cell,hybinp)
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
                   xcpot,noco,mpinp,hybinp)
      ENDIF
      !
      ! k-points can also be modified here
      !
      call make_kpoints(kpts,cell,sym,hybinp,input%film,noco%l_ss.or.noco%l_soc,kpts_str)
      !
      !Now the IO-section
      !
      IF (.NOT.l_inpxml.or.judft_was_argument("-overwrite")) THEN
         call determine_includes(l_include)
         !the inp.xml file
         !CALL dump_FleurInputSchema()
         CALL w_inpxml(&
              atoms,vacuum,input,stars,sliceplot,forcetheo,banddos,&
              cell,sym,xcpot,noco,oneD,mpinp,hybinp,kpts,enpara,&
              l_explicit,l_include,"inp.xml")
         if (.not.l_include(1)) CALL sym%print_XML(99,"sym.xml")
      ENDIF
      if (.not.l_include(2)) THEN
        inquire(file="kpts.xml",exist=l_include(2))
        if (l_include(2)) call system("rm kpts.xml")
        CALL kpts%print_XML(99,"kpts.xml")
      endif
      ! Structure in  xsf-format
      OPEN (55,file="struct.xsf")
      CALL xsf_WRITE_atoms(55,atoms,input%film,.FALSE.,cell%amat)
      CLOSE (55)
      CLOSE(6)

      CALL juDFT_end("All done")

    CONTAINS
      SUBROUTINE determine_includes(l_include)
        LOGICAL,INTENT(out)::l_include(4)  !kpts,operations,species,position

        CHARACTER(len=100)::str=''
        LOGICAL           ::incl

        l_include=[.FALSE.,.FALSE.,.TRUE.,.TRUE.]

        IF (judft_was_argument("-inc")) THEN
           str=judft_string_for_argument("-inc")

           DO WHILE(LEN_TRIM(str)>0)
              IF (str(1:1)=='-') THEN
                 incl=.FALSE.
                 str=str(2:)
              ELSE
                 incl=.TRUE.
                 IF (str(1:1)=='+') str=str(2:)
              ENDIF
              SELECT CASE(str(1:1))
              CASE ('k','K')
                 l_include(1)=incl
              CASE ('o','O')
                 l_include(2)=incl
              CASE ('s','S')
                 l_include(3)=incl
              CASE ('p','P')
                 l_include(4)=incl
              CASE ('a','A')
                 l_include(:)=incl
              END SELECT
              IF (INDEX(str,"'")>0) THEN
                 str=str(INDEX(str,"'")+1:)
              ELSE
                 str=""
              END IF
           END DO
        ENDIF

      END SUBROUTINE determine_includes

    END PROGRAM inpgen
