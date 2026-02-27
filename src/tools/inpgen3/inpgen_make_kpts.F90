!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------

MODULE m_inpgen_make_kpts
   IMPLICIT NONE
CONTAINS
   subroutine add_kpoints(kpts_str,kpts_path,l_onlyIdentitySym)
      use m_types_cell
      use m_types_sym
      use m_types_hybinp
      use m_types_input
      use m_types_kpts
      use m_types_noco
      use m_fleurinput_read_xml
      use m_make_kpoints
      use m_constants
      implicit none
      character(len=*), intent(inout) :: kpts_str,kpts_path
      logical, intent(in) :: l_onlyIdentitySym
      type(t_kpts):: kpts
      type(t_kpts), allocatable :: kpts_old(:)
      type(t_cell) :: cell
      type(t_sym) :: sym
      type(t_hybinp) :: hybinp
      type(t_input) :: input
      type(t_noco) :: noco
      character(len=40) :: kpts_name
      
  
      integer :: bz_integration,i,j
      logical :: l_new

      print *,"READING EXISTING inp.xml TO DETERMINE OLD K-POINT SET NAMES"
      call Fleurinput_read_xml(0,"",cell,sym,input=input,noco=noco,kptsArray=kpts_old,hybinp=hybinp)
      call cell%init(0.0)
      call sym%init(cell,input%film)


      !determine name of new kpoint set
      DO i=1,99 
         write(kpts_name, '(A,I0)') 'set_',i
         l_new=.true.
         DO j=1,size(kpts_old)
            IF (trim(kpts_name) == trim(kpts_old(j)%kptsName)) THEN
               l_new=.false.
               EXIT
            END IF
         END DO
         IF (l_new)  EXIT !this is a new name
      END DO
      print *,"NEW K-POINT SET NAME: ",trim(kpts_name)

      bz_integration=BZINT_METHOD_HIST !just a default here, will be overwritten in make_kpoints
      call make_kpoints(kpts, cell, sym, input%film, noco%l_soc, &
                       bz_integration, .false., kpts_str, kpts_name, &
                       kpts_path, l_onlyIdentitySym)
      print *,"ADDING NEW K-POINT SET TO inp.xml"
      call add_kpts_set(kpts)                 

   end subroutine add_kpoints


   subroutine add_kpts_set(kpts)
      use m_types_kpts
      implicit none
      type(t_kpts), intent(in) :: kpts
      character(len=500) :: line
      integer :: ierr, start_index, end_index
      logical :: l_written

      l_written=.false.

      !first open inp.xml,replace selected kpoint set and try to write there
      call copy_file('inp.xml','inp_old.xml')

      OPEN(UNIT=99,FILE='inp_old.xml',STATUS='OLD',FORM='FORMATTED',ACTION="READ",IOSTAT=ierr)
      OPEN(UNIT=100,FILE='inp.xml',STATUS='REPLACE',FORM='FORMATTED',ACTION="WRITE",IOSTAT=ierr)
      !Loop over full inp.xml file and find the k-points section
      do 
         READ(99, '(a)', IOSTAT=ierr) line
         IF (ierr /= 0) EXIT !done with file
         ! Check if the line contains something interesting
         IF (index(line, 'kPointListSelection') > 0) THEN
            ! Found the k-points selection tag
            start_index=index(line,'listName')
            start_index=index(line(start_index:),'"')+start_index-1
            end_index=index(line(start_index+1:),'"')+start_index
            write(100, '(3a)') line(1:start_index),trim(kpts%kptsName),trim(line(end_index:))
         elseif (index(line,"<kPointLists>")>0.and..not.l_written) then 
            write(100, '(a)') line
            call kpts%print_XML(100)
            l_written=.true.
         else
            write(100, '(a)') trim(line)
         end if
      end do 
      close(99)
      close(100)

      !if not written yet, open kpts.xml and write there
      if (.not.l_written) then
         call copy_file("kpts.xml","kpts_old.xml")
         OPEN(UNIT=99,FILE='kpts_old.xml',STATUS='OLD',FORM='FORMATTED',ACTION="READ",IOSTAT=ierr)
         OPEN(UNIT=100,FILE='kpts.xml',STATUS='REPLACE',FORM='FORMATTED',ACTION="WRITE",IOSTAT=ierr)
         !Loop over full kpts.xml file and find the k-points section
         do 
            READ(99, '(a)', IOSTAT=ierr) line
            IF (ierr /= 0) EXIT !done with file
            ! Check if the line contains kpointlists tag
            IF (index(line,"<kPointLists>")>0.and..not.l_written) then 
               write(100, '(a)') line
               call kpts%print_XML(100)
               l_written=.true.
            else
               write(100, '(a)') trim(line)
            end if
         end do 
      endif
      close(99)
      close(100)

      contains 
      subroutine copy_file(sourcefile,targetfile)
         character(len=*), intent(in) :: sourcefile,targetfile
         integer :: ierr
         character(len=500) :: line

         OPEN(UNIT= 99, FILE=sourcefile, STATUS='OLD', FORM='FORMATTED', ACTION="READ", IOSTAT=ierr)
         OPEN(UNIT= 100, FILE=targetfile, STATUS='REPLACE', FORM='FORMATTED', ACTION="WRITE", IOSTAT=ierr)
         do
            READ(99,'(a)',IOSTAT=ierr) line
            IF (ierr /= 0) EXIT
            WRITE(100,'(a)') line
         end do
         close(99)
         close(100)
      end subroutine copy_file

   end subroutine add_kpts_set


END MODULE m_inpgen_make_kpts