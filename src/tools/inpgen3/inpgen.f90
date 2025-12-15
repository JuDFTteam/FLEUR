!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

PROGRAM inpgen
!--------------------------------------------------------------------------------+
!   Set up a FLEUR inp.xml-file from basic input data; for use and docu please   !
!   refer to http://www.flapw.de/docs/inpgen.html)                               !
!   Two modes of operation are possible:                                         !
!   1) Generating a full inp.xml from a simple input file (use -f option)        !
!   2) Adding a k-point set to an existing inp.xml (use -kpt option)             !
!--------------------------------------------------------------------------------+
     USE m_juDFT
     USE m_inpgen_help
     USE m_inpgen_version
     USE m_fleur_dropxmlschema
     USE m_inpgen_make_inp
     USE m_xsf_io
     use m_inpgen_make_kpts 
      IMPLICIT NONE

      LOGICAL               :: l_fullinput,l_explicit(3),l_inpxml,l_include(4)

      LOGICAL               :: l_nosym,l_overwrite,l_noInpgenComment
      LOGICAL               :: l_dropXMLSchema,l_noksym
      CHARACTER(len=100)    :: kpt_string,profileName,kptsPath
      INTEGER               :: simple_inp_file_id
      INTEGER               :: oUnit,inpgenIUnit,inpOldUnit
      CHARACTER(len=100)    :: filename
      CHARACTER(len=200)    :: line
      LOGICAL               :: l_exist
      INTEGER               :: ios
   
      CALL judft_init(oUnit,.FALSE.)
      CALL timestart("inpgen")
      OPEN(oUnit,file='out')

      !Start program and greet user
      CALL inpgen_help()
      call inpgen_version()
      
      call process_arguments(l_nosym,l_explicit,filename,kpt_string,profileName,&
           kptsPath,l_overwrite,l_noInpgenComment,l_dropXMLSchema,l_noksym,l_include)
      if (l_dropXMLSchema) call fleur_dropxmlschema()

      if(len_trim(filename)>0) then
         ! mode: generate inp.xml from input file
         INQUIRE(file=filename,exist=l_exist)
         IF (.NOT.l_exist) CALL judft_error("Input file specified is not readable")
         INQUIRE(file='inp.xml',exist=l_inpxml)
         IF (l_inpxml .AND. .NOT.l_overwrite)&
              CALL judft_error("inp.xml exists and can not be overwritten")

         if (.not.l_noInpgenComment) THEN
            simple_inp_file_id=97
            open(simple_inp_file_id,file=filename,status='old',action='read')
            call make_inp_xml(simple_inp_file_id, "inp.xml",  .true., &
                              profileName, l_include, l_explicit,l_noksym)

            OPEN (inpgenIUnit,file=TRIM(filename),action="read")
            OPEN (inpOldUnit, file="inp.xml", action="write", status='old', access='append')
            WRITE(inpOldUnit,'(a)') ''
            WRITE(inpOldUnit,'(a)') '<!--'
            WRITE(inpOldUnit,'(a)') 'Command line when calling inpgen (only for documentation purposes):'
            CALL GET_COMMAND(line)
            WRITE(inpOldUnit,'(a)') TRIM(line)
            WRITE(inpOldUnit,'(a)') ''
            WRITE(inpOldUnit,'(a)') 'Initial (original) inpgen input (only for documentation purposes):'
            ios = 0
            DO WHILE(ios==0)
               READ(inpgenIUnit,'(a)',iostat=ios) line
               IF (ios.EQ.0) WRITE(inpOldUnit,'(a)') TRIM(line)
            END DO
            WRITE(inpOldUnit,'(a)') '-->'
            CLOSE (inpOldUnit)
            CLOSE (inpgenIUnit)
         endif
         ! Structure in  xsf-format
         !OPEN (55,file="struct.xsf")
         !CALL xsf_WRITE_atoms(55,atoms,input%film,cell%amat)
         !CLOSE (55)
      
      else
         !mode: add k-point set
         call add_kpoints(kpt_string, kptsPath, l_noksym)
      endif   
      CLOSE(oUnit)
      CALL timestop("inpgen")
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

         
      subroutine process_arguments(l_nosym,l_explicit,filename,kpt_string,profileName,kptsPath,l_overwrite,l_noInpgenComment,l_dropXMLSchema,l_noksym,l_include)
      implicit none
      logical, intent(out) :: l_nosym,l_overwrite,l_noInpgenComment,l_dropXMLSchema,l_noksym
      character(len=100), intent(out) :: filename, kpt_string, profileName, kptsPath
      logical, intent(out) :: l_include(4),l_explicit(3) 

      l_nosym=judft_was_argument("-nosym")
      l_explicit(1)=judft_was_argument("-explicit")
      l_explicit(2)=judft_was_argument("-noco")
      l_explicit(3)=judft_was_argument("-greensf")
      l_overwrite=judft_was_argument("-overwrite")
      l_noInpgenComment=juDFT_was_argument("-noInpgenComment")
      l_dropXMLSchema=juDFT_was_argument("-dropXMLSchema")
      l_noksym=judft_was_argument("-noKsym")

      IF (judft_was_argument("-profile")) THEN 
         profileName=TRIM(ADJUSTL(judft_string_for_argument("-profile"))) 
      ELSE
        profileName='default'
        print *, "================= WARNING ============ WARNING ================="
        print *, "No profile name provided, using 'default' settings"
        print *, "You are strongly advised to provide a profile name using the -profile option."
        print *, "=================================================================="
      ENDIF
      kptsPath = 'default'
      filename = ''
      kpt_string = ''
      if (juDFT_was_argument("-f")) filename=juDFT_string_for_argument("-f")         
      if (juDFT_was_argument("-kpt")) kpt_string=juDFT_string_for_argument("-kpt")
      IF (judft_was_argument("-kptsPath")) kptsPath = judft_string_for_argument("-kptsPath")
            
      call determine_includes(l_include)

      if(len_trim(filename)==0) then
         if (len_trim(kpt_string)==0) then
            CALL judft_error("No operating mode selected, use either the -f or the -kpt option")
         else
            Write(*,*) "Mode selection: adding k-point set"
         endif
      else
         if (len_trim(kpt_string)>0) then
            CALL judft_error("Two operating mode selected, use either the -f or the -kpt option")
         else
            Write(*,*) "Mode selection: generating inp.xml from input file"
         endif
      endif
      end subroutine process_arguments

    END PROGRAM inpgen
