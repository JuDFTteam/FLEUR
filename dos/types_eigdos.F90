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
    CHARACTER(len=20)            :: name_of_dos="unnamed"
    !each eigenvalue might be described by weights
    REAL,ALLOCATABLE             :: eig(:,:,:)
    REAL,ALLOCATABLE             :: dos_grid(:) !This is the grid the DOS part uses internally (FOR IO use the routine below)
    REAL,ALLOCATABLE             :: dos(:,:,:) !(grid,spins,weights)
  CONTAINS
    procedure(get_weight_name),DEFERRED    :: get_weight_name
    procedure(get_num_weights),DEFERRED    :: get_num_weights
    procedure(get_weight_eig),DEFERRED     :: get_weight_eig !should be overwritten in derived type
    procedure          :: get_spins
    procedure          :: get_neig
    procedure          :: get_eig
    procedure          :: get_dos_grid
    procedure          :: make_dos=>t_eigdos_make_dos
    procedure          :: smooth=>dosdata_smooth
    procedure          :: write_raw   !should be implemented later to allow eig66 functionality
    procedure          :: write_dos
    procedure          :: write_band
    procedure          :: write_EVData
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

  function get_neig(this)
    CLASS(t_eigdos),INTENT(IN)::this
    integer,allocatable :: get_neig(:,:)
    real,allocatable::ev(:,:,:)
    integer ::k,j

    ev=this%get_eig()
    allocate(get_neig(size(ev,2),size(ev,3)))
    DO j=1,this%get_spins()
      DO k=1,size(ev,2)
        get_neig(k,j)=count(ev(:,k,j)<1E99)
      ENDDO
    ENDDO
  end function

  pure integer function get_spins(this)
    CLASS(t_eigdos),INTENT(IN)::this
    get_spins=size(this%eig,3)
  END function
  function get_eig(this)
    CLASS(t_eigdos),INTENT(IN):: this
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
    DO jspin=1,eigdos%get_spins()
      CALL smooth(dos_grid,eigdos%dos(:,jspin,i),banddos%sig_dos,size(eigdos%dos_grid))
    ENDDO
  ENDDO
END subroutine

subroutine write_dos(eigdos,hdf_id)
#ifdef CPP_HDF
    use HDF5
    use m_banddos_io
#endif
    class(t_eigdos),INTENT(INOUT):: eigdos
#ifdef CPP_HDF
    integer(HID_T),intent(in) ::hdf_id
#else
    integer,       intent(in) ::hdf_id
#endif
    integer:: jspin,i,ind,id, n
    character(len=100)::filename
    real,allocatable:: dos_grid(:)
    LOGICAL l_printTextDOS

    l_printTextDOS = .TRUE.

#ifdef CPP_HDF
    DO n=1,eigdos%get_num_weights()
      print *, "writedos:",n,eigdos%get_num_weights()
      call writedosData(hdf_ID,eigdos%name_of_dos,eigdos%get_dos_grid(),eigdos%get_weight_name(n),eigdos%dos(:,:,n))
    enddo
    IF(eigdos%get_num_weights().GT.40) THEN
       WRITE(*,*) 'Number of weights in ', TRIM(ADJUSTL(eigdos%name_of_dos)),' DOS too large for simple text output.'
       WRITE(*,*) 'Output only in banddos.hdf file.'
       l_printTextDOS = .FALSE.
    END IF
#endif

    IF (.NOT.l_printTextDOS) RETURN
    if (.not.allocated(eigdos%dos)) return
    if (size(eigdos%dos)==0) return
    DO jspin=1,eigdos%get_spins()
      write(filename,"(a,a,i0)") trim(eigdos%name_of_dos),".",jspin
      open(999,file=filename)
      write(999,"(999a21)") "#energy",(eigdos%get_weight_name(id),id=1,eigdos%get_num_weights())
      write(*,"(999a21)") filename,(eigdos%get_weight_name(id),id=1,eigdos%get_num_weights())
      dos_grid=eigdos%get_dos_grid()
      DO i=1,size(dos_grid)
        write(999,"(999(e20.8,1x))") dos_grid(i)*hartree_to_ev_const,(eigdos%dos(i,jspin,id)/hartree_to_ev_const,id=1,eigdos%get_num_weights())
      ENDDO
      close(999)
      write(*,*) "done:",filename
    ENDDO
  END subroutine

  subroutine write_band(eigdos,kpts,title,cell,hdf_id,efermi,banddos)
    use m_types_kpts
    use m_types_cell
    use m_gnuplot_BS
    use m_types_banddos
#ifdef CPP_HDF
     use HDF5
     use m_banddos_io
#endif
    class(t_eigdos),INTENT(INOUT):: eigdos
    type(t_kpts),intent(in)      :: kpts
    type(t_banddos),INTENT(IN)   :: banddos
    type(t_cell),intent(in)      :: cell
    CHARACTER(LEN=*), INTENT(IN) :: title
    real,intent(in)              :: efermi
#ifdef CPP_HDF
    integer(HID_T),intent(in) ::hdf_id
    INTEGER::n

#else
    integer,intent(in):: hdf_id !not used
#endif

    integer:: jspin,i,k
    real,allocatable  :: ev(:,:,:),kx(:)
    real              :: vkr(3),vkr_prev(3)
    character(len=100)::file
    allocate(kx(kpts%nkpt))
#ifdef CPP_HDF
    DO n=1,eigdos%get_num_weights()
      call writebandData(hdf_id,kpts,eigdos%name_of_dos,eigdos%get_weight_name(n),eigdos%get_weight_eig(n),eigdos%get_eig())
    enddo
#endif
    if (eigdos%name_of_dos.ne."Local") then
#ifndef CPP_HDF
      print *,"WARNING, only very basic BS written without HDF"
#endif
      return
    endif
    DO jspin=1,eigdos%get_spins()
      write(file,"(a,i0)") "bands.",jspin
      open(18,file=file)
      ev=eigdos%get_eig()
      kx(1) = 0.0
      vkr_prev=matmul(kpts%bk(:,1),cell%bmat)
      DO k = 2, kpts%nkpt
        vkr=matmul(kpts%bk(:,k),cell%bmat)
        kx(k)=kx(k-1)+ sqrt(dot_product(vkr-vkr_prev,vkr-vkr_prev))
        vkr_prev=vkr
      ENDDO
      DO i = 1, minval(eigdos%get_neig())
        DO k = 1, kpts%nkpt
          write(18,'(999f15.9)') kx(k),(ev(i,k,jspin)-efermi)*hartree_to_ev_const
          !-eFermiCorrection
        ENDDO
      ENDDO
      CLOSE (18)
    enddo
    call gnuplot_bs(kpts,title,cell,eigdos%get_spins())
    IF (banddos%unfoldband) call write_gnu_sc(banddos,kpts,title,cell,eigdos%get_spins())
  end subroutine

  subroutine write_EVData(eigdos,hdf_id)
#ifdef CPP_HDF
     use HDF5
     use m_banddos_io
#endif
     class(t_eigdos),INTENT(INOUT):: eigdos
#ifdef CPP_HDF
     integer(HID_T),intent(in) ::hdf_id
     INTEGER::n
#else
     integer,intent(in):: hdf_id !not used
#endif

#ifdef CPP_HDF
     DO n = 1, eigdos%get_num_weights()
        CALL writeEVData(hdf_id,eigdos%name_of_dos,eigdos%get_weight_name(n),eigdos%get_eig(),eigdos%get_weight_eig(n))
     END DO
#endif
  end subroutine

  subroutine t_eigdos_make_dos(eigdos,kpts,input,banddos,efermi)
    use m_types_banddos
    use m_types_input
    use m_dosbin
    use m_ptdos
    use m_tetra_dos
    use m_dostetra
    use m_types_kpts

    class(t_eigdos),intent(inout):: eigdos
    type(t_banddos),intent(in)   :: banddos
    type(t_kpts),intent(in)      :: kpts
    type(t_input),intent(in)     :: input
    real,intent(in)              :: efermi

    integer ::n
    real :: emin,emax

    emin=min(banddos%e1_dos,banddos%e2_dos)-efermi
    emax=max(banddos%e1_dos,banddos%e2_dos)-efermi
    if (allocated(eigdos%dos)) return

    allocate(eigdos%dos(banddos%ndos_points,eigdos%get_spins(),eigdos%get_num_weights()))
    !Generate DOS grid
    if (allocated(eigdos%dos_grid)) deallocate(eigdos%dos_grid)
    allocate(eigdos%dos_grid(banddos%ndos_points))
    DO n=1,banddos%ndos_points
      eigdos%dos_grid(n)=emin+(emax-emin)/(banddos%ndos_points-1.0)*(n-1.0)
    ENDDO

    DO n=1,eigdos%get_num_weights()
      print *,eigdos%name_of_dos,n,eigdos%get_num_weights()
      SELECT CASE(input%bz_integration)

      CASE(BZINT_METHOD_HIST, BZINT_METHOD_GAUSS)
        call dos_bin(input%jspins,kpts%wtkpt,eigdos%dos_grid,eigdos%get_eig(),eigdos%get_weight_eig(n),eigdos%dos(:,:,n),efermi)
      CASE(BZINT_METHOD_TRIA)
        IF(input%film) THEN
          CALL ptdos(input%jspins,kpts,eigdos%dos_grid,eigdos%get_eig(),eigdos%get_weight_eig(n),eigdos%dos(:,:,n),efermi)
        ELSE
          CALL tetra_dos(input%jspins,kpts,eigdos%dos_grid,eigdos%get_neig(),eigdos%get_eig(),&
                         eigdos%get_weight_eig(n),eigdos%dos(:,:,n),efermi)
        ENDIF
      CASE(BZINT_METHOD_TETRA)
        CALL dostetra(kpts,input,eigdos%dos_grid,eigdos%get_eig(),eigdos%get_weight_eig(n),eigdos%dos(:,:,n),efermi)
      END SELECT

    end do
  END subroutine



  subroutine write_raw(this,id)
    class(t_eigdos),INTENT(IN):: this
    INTEGER,INTENT(IN)        :: id


  end subroutine



END MODULE m_types_eigdos
