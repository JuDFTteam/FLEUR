!--------------------------------------------------------------------------------
! Copyright (c) 2020 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_mcd
  use m_judft
  use m_types_eigdos
  implicit none
  PRIVATE
  public t_mcd
  TYPE,extends(t_eigdos):: t_mcd
    REAL                 :: emcd_lo, emcd_up, maxE_mcd
    INTEGER              :: numInvolvedAtomTypes
    INTEGER, ALLOCATABLE :: atomTypeIndices(:)
    INTEGER, ALLOCATABLE :: ncore(:)
    INTEGER,ALLOCATABLE  :: n_dos_to_type(:)
    REAL,    ALLOCATABLE :: e_mcd(:,:,:)
    REAL,    ALLOCATABLE :: mcd(:,:,:,:,:)
    COMPLEX, ALLOCATABLE :: m_mcd(:,:,:,:)
    REAL,    ALLOCATABLE :: mcd_grid(:)
  CONTAINS
    PROCEDURE,PASS :: init => mcd_init
    procedure      :: get_dos_grid
    PROCEDURE      :: make_dos
    PROCEDURE      :: get_weight_eig
    PROCEDURE      :: get_num_weights
    PROCEDURE      :: get_weight_name
    END TYPE t_mcd
contains

function get_dos_grid(this)
  class(t_mcd),intent(in):: this
  real,allocatable:: get_dos_grid(:)

  INTEGER :: ind,ntype,nc
  REAL:: e_core
  get_dos_grid=this%mcd_grid
end function

subroutine make_dos(eigdos,kpts,input,banddos,efermi)
    use m_types_banddos
    use m_types_input
    use m_types_kpts

    class(t_mcd),intent(inout)   :: eigdos
    type(t_banddos),intent(in)   :: banddos
    type(t_input),intent(in)     :: input
    type(t_kpts),intent(in)      :: kpts
    real,intent(in)              :: efermi

    integer ::n,i,ind,ntype,nc,k,l,jspin
    real    :: e1,e2,e_lo,e_up,fac
    real,allocatable:: dos(:,:,:), totDos(:,:,:)
    if (allocated(eigdos%dos)) return
    !Call the routine of the parent-class
    call t_eigdos_make_dos(eigdos,kpts,input,banddos,efermi)

    !Only unoccupied states
    DO n=1,size(eigdos%dos_grid)
      if (eigdos%dos_grid(n)<0.0)  eigdos%dos(n,:,:)=0.0
    enddo

    !Map the values to MCD grid

    e_lo =  minval(eigdos%e_mcd)-efermi-maxval(eigdos%dos_grid) - 3.0*banddos%sig_dos
    e_up =  eigdos%maxE_mcd-efermi + 3.0*banddos%sig_dos
    ALLOCATE(eigdos%mcd_grid(size(eigdos%dos_grid)))
    DO i=1,size(eigdos%dos_grid)
      eigdos%mcd_grid(i)=e_lo+(i-1)*(e_up-e_lo)/(size(eigdos%mcd_grid)-1)
    ENDDO

    allocate(dos,mold=eigdos%dos)
    ALLOCATE(totDOS(size(eigdos%mcd_grid),size(eigdos%e_mcd,2),size(eigdos%ncore)*3))
    dos=0.0
    ind=0
    DO ntype=1,size(eigdos%ncore)
      DO nc=1,eigdos%ncore(ntype)
        DO k = 1,3
          ind=ind+1
          DO jspin=1,size(eigdos%e_mcd,2)
            DO i=1,size(eigdos%dos_grid)-1
              e1=-1*eigdos%dos_grid(i)-efermi+eigdos%e_mcd(ntype,jspin,nc)
              e2=-1*eigdos%dos_grid(i+1)-efermi+eigdos%e_mcd(ntype,jspin,nc)
              DO l=1,size(eigdos%mcd_grid)
                IF ((e2.LE.eigdos%mcd_grid(l)).AND. (e1.GT.eigdos%mcd_grid(l))) THEN
                  fac = (eigdos%mcd_grid(l)-e1)/(e2-e1)
                  dos(l,jspin,ind) = dos(l,jspin,ind)+ eigdos%dos(i,jspin,ind)*(1.-fac) + fac * eigdos%dos(i+1,jspin,ind)
                  totDos(l,jspin,3*(ntype-1)+k) = totDos(l,jspin,3*(ntype-1)+k) + eigdos%dos(i,jspin,ind)*(1.-fac) + fac * eigdos%dos(i+1,jspin,ind)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    DO k = 1, size(eigdos%ncore)*3
       dos(:,:,ind+k) = totDos(:,:,k)
    END DO
    DEALLOCATE(totDos)
    
    eigdos%dos=dos
  end subroutine

function get_weight_eig(this,id)
  class(t_mcd),intent(in):: this
  INTEGER,intent(in)         :: id
  real,allocatable:: get_weight_eig(:,:,:)

  INTEGER :: ind,ntype,nc, typeIndex, offset, n_dos, mcdIndex
  ind=0
  DO ntype=1,size(this%ncore)
    DO nc=1,this%ncore(ntype)
      ind=ind+1
      if (ind==id) get_weight_eig=this%mcd(3*(ntype-1)+1,nc,:,:,:)
      ind=ind+1
      if (ind==id) get_weight_eig=this%mcd(3*(ntype-1)+2,nc,:,:,:)
      ind=ind+1
      if (ind==id) get_weight_eig=this%mcd(3*(ntype-1)+3,nc,:,:,:)
      IF(ind>id) return
    ENDDO
  ENDDO

  typeIndex = (id - 3*sum(this%ncore) - 1) / 3 + 1
  offset = mod(id,3)

  IF (offset.EQ.0) offset = 3
  IF ((id.GT.3*sum(this%ncore)).AND.(typeIndex.LE.this%numInvolvedAtomTypes)) THEN
     ALLOCATE(get_weight_eig(SIZE(this%mcd,3),SIZE(this%mcd,4),SIZE(this%mcd,5)))
     get_weight_eig = 0.0
     DO ntype=1,size(this%ncore)
        DO nc=1,this%ncore(ntype)
           mcdIndex = 3*(ntype-1)+offset
           n_dos = (mcdIndex-1)/3+1
           IF (this%atomTypeIndices(this%n_dos_to_type(n_dos)).EQ.typeIndex) THEN
              get_weight_eig(:,:,:) = get_weight_eig(:,:,:) + this%mcd(mcdIndex,nc,:,:,:)
           END IF
        END DO
     END DO
  END IF
  IF (ind>id) CALL judft_error("Types_mcd: data not found")

END function

integer function get_num_weights(this)
  class(t_mcd),intent(in):: this
  
  get_num_weights=3*sum(this%ncore) + 3*this%numInvolvedAtomTypes
end function

  character(len=20) function get_weight_name(this,id)
    class(t_mcd),intent(in):: this
    INTEGER,intent(in)         :: id

    character(len=3):: c
    INTEGER :: ind,n_dos,nc,n, typeIndex, iType
    select case(mod(id,3))
    case(1)
      c="pos"
    case(2)
      c="neg"
    case(0)
      c="cir"
    end select
    ind=0
    DO n=1,size(this%mcd,1)
      n_dos=(n-1)/3+1
      DO nc=1,this%ncore(n_dos)
        ind=ind+1
        if (ind==id) THEN
          write(get_weight_name,"(a,i0,a,i0,a)") "At",this%n_dos_to_type(n_dos),"NC",nc,c
          RETURN
        ENDIF
      ENDDO
    ENDDO
    typeIndex = (id - 3*sum(this%ncore) - 1) / 3 + 1
    IF ((typeIndex.GT.0).AND.(typeIndex.LE.this%numInvolvedAtomTypes)) THEN
       DO iType = 1, SIZE(this%atomTypeIndices)
          IF(typeIndex.EQ.this%atomTypeIndices(iType)) THEN
             write(get_weight_name,"(a,i0,a)") "At",typeIndex,c
             RETURN
          END IF
       END DO
    END IF
    IF(ind>id) then
       CALL judft_error("Types_mcd: data not found")
    END IF
  end function



SUBROUTINE mcd_init(thisMCD,banddos,input,atoms,kpts,eig)
  USE m_types_setup
  USE m_types_kpts

   IMPLICIT NONE

   CLASS(t_mcd),          INTENT(INOUT) :: thisMCD
   TYPE(t_banddos),       INTENT(IN)    :: banddos

   TYPE(t_input),         INTENT(IN)    :: input
   TYPE(t_atoms),         INTENT(IN)    :: atoms
   TYPE(t_kpts),          INTENT(IN)    :: kpts
   real,INTENT(IN)                      :: eig(:,:,:)

   integer :: ntype !no of types for which MCD is calculated
   INTEGER :: numInvolvedAtomTypes, i
   LOGICAL :: involvedAtomTypes(atoms%ntype)
   
   ALLOCATE (thisMCD%atomTypeIndices(atoms%ntype))
   thisMCD%atomTypeIndices = 0
   involvedAtomTypes = .FALSE.
   thisMCD%n_dos_to_type=banddos%dos_typelist
   ntype=size(banddos%dos_typelist)
   DO i = 1, ntype
      involvedAtomTypes(thisMCD%n_dos_to_type(i)) = .TRUE.
   END DO
   numInvolvedAtomTypes = 0
   DO i = 1, atoms%ntype
      IF (involvedAtomTypes(i)) numInvolvedAtomTypes = numInvolvedAtomTypes + 1
      thisMCD%atomTypeIndices(i) = numInvolvedAtomTypes
   END DO
   thisMCD%numInvolvedAtomTypes = numInvolvedAtomTypes
   
   thisMCD%name_of_dos="MCD"
   ALLOCATE (thisMCD%ncore(ntype))
   ALLOCATE (thisMCD%e_mcd(ntype,input%jspins,29))
   IF (banddos%l_mcd) THEN
      thisMCD%emcd_lo = banddos%e_mcd_lo
      thisMCD%emcd_up = banddos%e_mcd_up
      ALLOCATE (thisMCD%m_mcd(29,(3+1)**2,3*ntype,2))
      ALLOCATE (thisMCD%mcd(3*ntype,29,input%neig,kpts%nkpt,input%jspins) )
      IF (.NOT.banddos%dos) WRITE (*,*) 'For mcd-spectra set banddos%dos=T!'
   ELSE
      ALLOCATE(thisMCD%dos(0,0,0)) !indicated no DOS should be calculated
      ALLOCATE (thisMCD%m_mcd(1,1,1,1))
      ALLOCATE (thisMCD%mcd(1,1,1,1,input%jspins))
   ENDIF

   thisMCD%maxE_mcd = -1000000.0
   thisMCD%ncore = 0
   thisMCD%e_mcd = 0.0
   thisMCD%mcd = 0.0
   thisMCD%m_mcd = CMPLX(0.0,0.0)

   thisMCD%eig=eig




END SUBROUTINE mcd_init
end module m_types_mcd
