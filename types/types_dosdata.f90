!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_dosdata
  USE m_juDFT
  IMPLICIT NONE


  PRIVATE
  PUBLIC:: t_eigdesc

  TYPE:: t_dos_data
    !each eigenvalue might be described by weights
    CHARACTER(len=20),ALLOCATABLE:: weight_names(:)!This must be allocated in init of derived type
    real,ALLOCATABLE             :: e(:)
    real,allocatable             :: dos(:,:,:)
  CONTAINS
    PROCEDURE :: init
    PROCEDURE :: smooth
    PROCEDURE :: write
  END TYPE


CONTAINS

  subroutine write(dos_data,filename)
    class(t_dos_data),INTENT(INOUT):: dos_data
    character(len=*),OPTIONAL,intent(in)::filename

    integer:: jspin,i,ind
    character(len=100)::file

    DO jspin=1,size(dos_data%dos,2)
      if (present(filename)) THEN
        write(file,"(a,a,i0)") trim(adjustl(filename)),".",jspin
      ELSE
        write(file,"(a,i0)") "DOS.",jspin
      ENDIF
      open(999,file)
      write(999,"(999a21)") "#energy",dos_data%weight_names
      DO i=1,size(dos_data%e)
        write(999,"(999(f12.8,1x))") dos_data%e(i),dos_data%dos(i,jspin,:)
      ENDDO
      close(999)
    ENDDO
  END subroutine

subroutine smooth(dos_data,banddos)
  use m_types_banddos
  class(t_dos_data),INTENT(INOUT):: dos_data
  type(t_banddos),INTENT(IN)     :: banddos

  integer :: jspin,i

  if (banddos%sig_dos==0.0) RETURN

  IF ( NINT((maxval(dos_data%e) - minval(dos_data%e))/(banddos%sig_dos*hartree_to_ev_const)) > size(dos_data%e) ) THEN
    WRITE(oUnit,*) 'sig_dos too small for DOS smoothing:'
    WRITE(oUnit,*) 'Reduce energy window or enlarge banddos%sig_dos!'
    WRITE(oUnit,*) 'For now: no smoothing done'
    return
  ENDIF

  DO i=1,size(dos_data%dos,3)
    DO jspin=1,size(dos_data%dos,2)
      CALL smooth(dos_data%e,dos%data%dos(:,jspin,i),banddos%sig_dos*hartree_to_ev_const,size(dos_data%e))
  ENDDO
END subroutine


subroutine init(dos_data,ned,jspins,emin,emax,efermi,eigdesc)
  class(t_dos_data),INTENT(OUT):: dos_data
  integer,intent(in)           :: jspins,ned
  REAL,INTENT(IN)              :: emin,emax,efermi
  class(t_eigdesc),INTENT(IN)  :: eigdesc(:) !These are e.g. of t_dos,t_orbcomp...

  integer:: ind
  ind=0
  DO i=1,size(eigdesc)
    ind=ind+eigdesc(i)%get_num_weights()
  ENDDO
  ALLOCATE(dos_data%weight_names(ind),dos_data%dos(ned,jspins,ind))
  ind=0
  DO i=1,size(eigdesc)
    DO j=1,eigdesc(i)%get_num_weights()
      ind=ind+1
      dos%weight_names(ind)=eigdesc(i)%get_weight_names(j)
    ENDDO
  ENDDO

  dos_data%e=generate_grid(ned,emin,emax,efermi)

end subroutine

function generate_grid(ned,emin_arg,emax_arg,efermi_arg)result(e)
  INTEGER,INTENT(IN) :: ned
  REAL,INTENT(IN)    :: emin_arg,emax_arg,efermi_arg
  REAL,ALLOCATABLE   :: e(:)

  REAL :: emin,emax,efermi,de
  INTEGER :: i
  ! scale energies
  emin =emin_arg*hartree_to_ev_const
  emax =emax_arg*hartree_to_ev_const
  efermi = efermiarg*hartree_to_ev_const
  !
  !     create energy grid
  !
  emax = max(emin,emax) - efermi
  emin = min(emax,emin) - efermi
  de = (emax-emin)/(ned-1)

  ALLOCATE(e(ned+1))
  DO i=1,ned+1
    e(i) = emin + (i-1)*de
  ENDDO

end function

END MODULE m_types_dosdata
