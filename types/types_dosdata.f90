!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_dosdata
  USE m_juDFT
  use m_types
  use m_constants
  IMPLICIT NONE

  PRIVATE
  PUBLIC:: t_dos_data
  TYPE:: t_dos_data
    !each eigenvalue might be described by weights
    CHARACTER(len=20),ALLOCATABLE:: weight_name(:)!This must be allocated in init of derived type
    real,ALLOCATABLE             :: e(:)
    real,allocatable             :: dos(:,:,:)
  CONTAINS
    PROCEDURE :: init
    PROCEDURE :: smooth=>dosdata_smooth
    PROCEDURE :: write
  END TYPE


CONTAINS

  

subroutine dosdata_smooth(dos_data,banddos)
  use m_smooth
  use m_types_banddos
  class(t_dos_data),INTENT(INOUT):: dos_data
  type(t_banddos),INTENT(IN)     :: banddos

  integer :: jspin,i

  if (banddos%sig_dos==0.0) RETURN

  IF ( NINT((maxval(dos_data%e) - minval(dos_data%e))/banddos%sig_dos) > size(dos_data%e) ) THEN
    WRITE(oUnit,*) 'sig_dos too small for DOS smoothing:'
    WRITE(oUnit,*) 'Reduce energy window or enlarge banddos%sig_dos!'
    WRITE(oUnit,*) 'For now: no smoothing done'
    return
  ENDIF

  DO i=1,size(dos_data%dos,3)
    DO jspin=1,size(dos_data%dos,2)
      CALL smooth(dos_data%e,dos_data%dos(:,jspin,i),banddos%sig_dos,size(dos_data%e))
    ENDDO
  ENDDO
END subroutine


subroutine init(dos_data,banddos,jspins,efermi,eigdesc)
  use m_types_eigdesc
  class(t_dos_data),INTENT(OUT):: dos_data
  type(t_banddos),INTENT(IN)   :: banddos
  integer,intent(in)           :: jspins
  REAL,INTENT(IN)              :: efermi
  class(t_eigdesc_list),INTENT(IN)  :: eigdesc(:) !These are e.g. of t_dos,t_orbcomp...

  integer:: ind,i,j
  ind=0
  DO i=1,size(eigdesc)
    ind=ind+eigdesc(i)%p%get_num_weights()
  ENDDO
  ALLOCATE(dos_data%weight_name(ind),dos_data%dos(banddos%ndos_points,jspins,ind))
  ind=0
  DO i=1,size(eigdesc)
    DO j=1,eigdesc(i)%p%get_num_weights()
      ind=ind+1
      dos_data%weight_name(ind)=eigdesc(i)%p%get_weight_name(j)
    ENDDO
  ENDDO

  dos_data%e=generate_grid(banddos%ndos_points,banddos%e1_dos,banddos%e2_dos,efermi)

end subroutine

function generate_grid(ned,emin_arg,emax_arg,efermi_arg)result(e)
  INTEGER,INTENT(IN) :: ned
  REAL,INTENT(IN)    :: emin_arg,emax_arg,efermi_arg
  REAL,ALLOCATABLE   :: e(:)

  REAL :: emin,emax,efermi,de
  INTEGER :: i
  ! scale energies
  emin =emin_arg
  emax =emax_arg
  efermi = efermi_arg
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
