!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_kpts
  
  TYPE t_kpts
     INTEGER               :: nkpt=0
     INTEGER               :: ntet=0
     LOGICAL               :: l_gamma=.false.
     !(3,nkpt) k-vectors internal units
     REAL,ALLOCATABLE      :: bk(:,:)
     !(nkpts) weights
     REAL,ALLOCATABLE      :: wtkpt(:)
     INTEGER               :: nkptf=0   !<k-vectors in full BZ
     INTEGER               :: nkpt3(3)=[0,0,0]
     REAL   ,ALLOCATABLE   :: bkf(:,:)
     INTEGER,ALLOCATABLE   :: bkp(:)
     INTEGER,ALLOCATABLE   :: bksym(:)
     INTEGER                       :: numSpecialPoints=0
     INTEGER, ALLOCATABLE          :: specialPointIndices(:)
     CHARACTER(LEN=50),ALLOCATABLE :: specialPointNames(:)
     REAL   ,ALLOCATABLE           :: specialPoints(:,:)
     INTEGER,ALLOCATABLE           :: ntetra(:,:)
     REAL   ,ALLOCATABLE           :: voltet(:)
     REAL   ,ALLOCATABLE           :: sc_list(:,:) !list for all information about folding of bandstructure (need for unfoldBandKPTS)((k(x,y,z),K(x,y,z),m(g1,g2,g3)),(nkpt),k_original(x,y,z))
   contains
     procedure :: init_defaults
     procedure :: init_by_density
     procedure :: init_by_number
     procedure :: init_by_grid
     procedure :: init_special
     procedure :: read_xml
     procedure :: add_special_line
  ENDTYPE t_kpts

contains
  subroutine init_defaults(kpts,film,area,omtil,nop2,nop)
    class(t_kpts),intent(out):: kpts
    logical,intent(in)       :: film
    real,intent(in)          :: area,omtil
    integer,intent(in)       :: nop2,nop

    integer:: nkpt
    IF (film) THEN
       nkpt = MAX(nint((3600/area)/nop2),1)
    ELSE
       nkpt = MAX(nint((216000/omtil)/nop),1)
    ENDIF
    call kpts%init_by_number(nkpt)
  end subroutine init_defaults

    subroutine init_by_density(kpts,density,bmat)
    class(t_kpts),intent(out):: kpts
    real,intent(in)          :: density,bmat(3,3)

    real    :: length
    integer :: n,grid(3)

    do n=1,3
       length=sqrt(dot_product(bmat(n,:),bmat(n,:)))  !TODO why not bmat(:,n)???
       grid(n)=ceiling(density*length)
    end do
    call kpts%init_by_grid(grid)
    call kpts%init_grid(kpts
  end subroutine init_by_density

  
END MODULE m_types_kpts
