!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_eigdos
  USE m_juDFT
  use m_constants
  IMPLICIT NONE
  PRIVATE
  PUBLIC:: t_eigdos,t_eigdos_list,t_eigdos_make_dos

  TYPE,abstract :: t_eigdos
    CHARACTER(len=20)            :: name_of_dos="DOS"
    !each eigenvalue might be described by weights
    REAL,ALLOCATABLE             :: eig(:,:,:)
    REAL,ALLOCATABLE             :: dos_grid(:) !This is the grid the DOS part uses internally (FOR IO use the routine below)
    REAL,ALLOCATABLE             :: dos(:,:,:) !(grid,spins,weights)
  CONTAINS
    procedure(get_weight_name),DEFERRED    :: get_weight_name
    procedure(get_num_weights),DEFERRED    :: get_num_weights
    procedure(get_weight_eig),DEFERRED     :: get_weight_eig !should be overwritten in derived type
    procedure          :: get_dims
    procedure          :: get_eig
    procedure          :: get_dos_grid
    procedure          :: make_dos=>t_eigdos_make_dos
    procedure          :: smooth=>dosdata_smooth
    procedure          :: write_raw
    procedure          :: write_dos
  END TYPE

  type::t_eigdos_list
    class(t_eigdos),POINTER:: p
  end type

  INTERFACE
    function get_weight_eig(this,id)
      import t_eigdos
      class(t_eigdos),intent(in):: this
      INTEGER,intent(in)         :: id
      real,allocatable:: get_weight_eig(:,:,:)
    end function
  END interface
  INTERFACE
    integer function get_num_weights(this)
      import t_eigdos
      class(t_eigdos),intent(in):: this
    end function
  end interface
  INTERFACE
    character(len=20) function get_weight_name(this,id)
      import t_eigdos
      class(t_eigdos),intent(in):: this
      INTEGER,intent(in)         :: id
    end function
  end interface
CONTAINS


  function get_dims(this)
    CLASS(t_eigdos),INTENT(IN)::this
    INTEGER :: get_dims(2)
    get_dims(1)=size(this%eig,1)
    get_dims(2)=size(this%eig,3)
  END function
  function get_eig(this,id)
    CLASS(t_eigdos),INTENT(IN):: this
    INTEGER,INTENT(IN)     :: id
    real,allocatable:: get_eig(:,:,:)
    get_eig=this%eig
  END function

  function get_dos_grid(this)
    CLASS(t_eigdos),INTENT(IN):: this
    real,allocatable:: get_dos_grid(:)
    get_dos_grid=this%dos_grid
  END function

subroutine dosdata_smooth(eigdos,banddos)
  use m_smooth
  use m_types_banddos
  class(t_eigdos),INTENT(INOUT)  :: eigdos
  type(t_banddos),INTENT(IN)     :: banddos

  integer :: jspin,i
  real,allocatable :: dos_grid(:)

  if (.not.allocated(eigdos%dos)) return
  if (size(eigdos%dos)==0) return
  if (size(eigdos%dos)<1) return
  if (banddos%sig_dos==0.0) RETURN
  dos_grid=eigdos%get_dos_grid()

  IF ( NINT((maxval(dos_grid) - minval(dos_grid))/banddos%sig_dos) > size(dos_grid) ) THEN
    WRITE(oUnit,*) 'sig_dos too small for DOS smoothing:'
    WRITE(oUnit,*) 'Reduce energy window or enlarge banddos%sig_dos!'
    WRITE(oUnit,*) 'For now: no smoothing done'
    return
  ENDIF

  DO i=1,size(eigdos%dos,3)
    DO jspin=1,size(eigdos%dos,2)
      CALL smooth(dos_grid,eigdos%dos(:,jspin,i),banddos%sig_dos,size(eigdos%dos_grid))
    ENDDO
  ENDDO
END subroutine

subroutine write_dos(eigdos,filename)
    class(t_eigdos),INTENT(INOUT):: eigdos
    character(len=*),OPTIONAL,intent(in)::filename

    integer:: jspin,i,ind,id
    character(len=100)::file
    real,allocatable:: dos_grid(:)

    if (.not.allocated(eigdos%dos)) return
    if (size(eigdos%dos)==0) return

    DO jspin=1,size(eigdos%dos,2)
      if (present(filename)) THEN
        write(file,"(a,a,i0)") trim(adjustl(filename)),".",jspin
      ELSE
        write(file,"(a,a,i0)") trim(eigdos%name_of_dos),".",jspin
      ENDIF
      open(999,file=file)
      write(999,"(999a21)") "#energy",(eigdos%get_weight_name(id),id=1,eigdos%get_num_weights())
      write(*,"(999a21)") file,(eigdos%get_weight_name(id),id=1,eigdos%get_num_weights())
      dos_grid=eigdos%get_dos_grid()
      DO i=1,size(dos_grid)
        write(999,"(999(e20.8,1x))") dos_grid(i)*hartree_to_ev_const,(eigdos%dos(i,jspin,id),id=1,eigdos%get_num_weights())
      ENDDO
      close(999)
      write(*,*) "done:",file
    ENDDO
  END subroutine

  subroutine t_eigdos_make_dos(eigdos,kpts,input,banddos,efermi)
    use m_types_banddos
    use m_types_input
    use m_dosbin
    use m_dostetra
    use m_types_kpts

    class(t_eigdos),intent(inout):: eigdos
    type(t_banddos),intent(in)   :: banddos
    type(t_kpts),intent(in)      :: kpts
    type(t_input),intent(in)     :: input
    real,intent(in)              :: efermi

    integer ::n,dims(2)
    real :: emin,emax

    emin=min(banddos%e1_dos,banddos%e2_dos)-efermi
    emax=max(banddos%e1_dos,banddos%e2_dos)-efermi
    if (allocated(eigdos%dos)) return
    dims=eigdos%get_dims()
    allocate(eigdos%dos(banddos%ndos_points,dims(2),eigdos%get_num_weights()))
    !Generate DOS grid
    if (allocated(eigdos%dos_grid)) deallocate(eigdos%dos_grid)
    allocate(eigdos%dos_grid(banddos%ndos_points))
    DO n=1,banddos%ndos_points
      eigdos%dos_grid(n)=emin+(emax-emin)/(banddos%ndos_points-1.0)*(n-1.0)
    ENDDO

    DO n=1,eigdos%get_num_weights()
      print *,eigdos%name_of_dos,n
      if (kpts%ntet==0) then
        call dos_bin(input%jspins,kpts%wtkpt,eigdos%dos_grid,eigdos%get_eig(n),eigdos%get_weight_eig(n),eigdos%dos(:,:,n))
      ELSE
        CALL dostetra(kpts,input,eigdos%dos_grid,eigdos%get_eig(n),eigdos%get_weight_eig(n),eigdos%dos(:,:,n))
      endif
    end do
  END subroutine



  subroutine write_raw(this,id)
    class(t_eigdos),INTENT(IN):: this
    INTEGER,INTENT(IN)        :: id


  end subroutine



END MODULE m_types_eigdos
