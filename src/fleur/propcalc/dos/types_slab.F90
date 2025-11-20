!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_slab
  use m_judft
  use m_types_eigdos
  implicit none
  PRIVATE
  public t_slab
   TYPE,extends(t_eigdos):: t_slab
      INTEGER              :: nsld, nsl

      INTEGER, ALLOCATABLE :: nmtsl(:,:)
      INTEGER, ALLOCATABLE :: nslat(:,:)
      REAL,    ALLOCATABLE :: zsl(:,:)
      REAL,    ALLOCATABLE :: volsl(:)
      REAL,    ALLOCATABLE :: volintsl(:)
      REAL,    ALLOCATABLE :: qintsl(:,:,:,:)
      REAL,    ALLOCATABLE :: qmtsl(:,:,:,:)

      CONTAINS
         PROCEDURE,PASS :: init => slab_init
         PROCEDURE      :: get_num_weights
         PROCEDURE      :: get_weight_eig
         PROCEDURE      :: get_weight_name
         PROCEDURE      :: postprocessing
         PROCEDURE,PASS :: calc_mt_slab
         PROCEDURE,PASS :: calc_int_slab 
   END TYPE t_slab
CONTAINS
  subroutine postprocessing(this, noco,nococonv, banddos,alldos, ef)
      use m_types_atoms
      use m_types_noco
      use m_types_nococonv
      use m_types_banddos
      class(t_slab), intent(inout):: this
      TYPE(t_noco), INTENT(IN)    :: noco
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_banddos), INTENT(IN)    :: banddos
      class(t_eigdos_list), intent(in), optional :: alldos(:)
      real, intent(in), optional :: ef

      return !currently no postprocessing needed for slab
   end subroutine postprocessing 
  integer function get_num_weights(this)
    class(t_slab),intent(in):: this
    get_num_weights=2*this%nsl
  END function

  character(len=20) function get_weight_name(this,id)
    class(t_slab),intent(in):: this
    INTEGER,intent(in)         :: id

    INTEGER :: ind,n
    ind=0
    DO n=1,this%nsl
      ind=ind+1
      if (ind==id) write(get_weight_name,"(a,i0)") "SLAB(INT):",n
      ind=ind+1
      if (ind==id) write(get_weight_name,"(a,i0)") "SLAB(MT):",n
      IF(ind>id) return
    ENDDO
  end function

  function get_weight_eig(this,id)
    class(t_slab),intent(in):: this
    INTEGER,intent(in)      :: id
    real,allocatable:: get_weight_eig(:,:,:)


    INTEGER :: ind,n
    ind=0
    DO n=1,this%nsl
      ind=ind+1
      if (ind==id) get_weight_eig=this%qintsl(n,:,:,:)
      ind=ind+1
      if (ind==id) get_weight_eig=this%qmtsl(n,:,:,:)
    ENDDO
  end function



  SUBROUTINE slab_init(thisSlab,banddos,atoms,cell,input,kpts)
   USE m_types_setup
   USE m_types_kpts
   USE m_slabdim
   USE m_slabgeom

   IMPLICIT NONE

   CLASS(t_slab),      INTENT(INOUT) :: thisSlab
   TYPE(t_banddos),    INTENT(IN)    :: banddos

   TYPE(t_atoms),      INTENT(IN)    :: atoms
   TYPE(t_cell),       INTENT(IN)    :: cell
   TYPE(t_input),      INTENT(IN)    :: input
   TYPE(t_kpts),       INTENT(IN)    :: kpts

   INTEGER :: nsld
   thisSlab%l_initialized = .TRUE.
   thisSlab%name_of_dos="SLAB"
   nsld=1
   IF (banddos%l_slab.AND.banddos%dos) THEN
      CALL slab_dim(atoms, nsld)
      ALLOCATE (thisSlab%nmtsl(atoms%ntype,nsld))
      ALLOCATE (thisSlab%nslat(atoms%nat,nsld))
      ALLOCATE (thisSlab%zsl(2,nsld))
      ALLOCATE (thisSlab%volsl(nsld))
      ALLOCATE (thisSlab%volintsl(nsld))
      ALLOCATE (thisSlab%qintsl(nsld,input%neig,kpts%nkpt,input%jspins))
      ALLOCATE (thisSlab%qmtsl(nsld,input%neig,kpts%nkpt,input%jspins))
      CALL slabgeom(atoms,cell,nsld,thisSlab%nsl,thisSlab%zsl,thisSlab%nmtsl,&
                    thisSlab%nslat,thisSlab%volsl,thisSlab%volintsl)
   ELSE
     allocate(thisSlab%dos(0,0,0))
      ALLOCATE (thisSlab%nmtsl(1,1))
      ALLOCATE (thisSlab%nslat(1,1))
      ALLOCATE (thisSlab%zsl(1,1))
      ALLOCATE (thisSlab%volsl(1))
      ALLOCATE (thisSlab%volintsl(1))
      ALLOCATE (thisSlab%qintsl(1,1,1,input%jspins))
      ALLOCATE (thisSlab%qmtsl(1,1,1,input%jspins))
   END IF
   thisSlab%nsld = nsld

   thisSlab%nmtsl = 0
   thisSlab%nslat = 0
   thisSlab%zsl = 0.0
   thisSlab%volsl = 0.0
   thisSlab%volintsl = 0.0
   thisSlab%qintsl = 0.0
   thisSlab%qmtsl = 0.0

END SUBROUTINE slab_init

  !***********************************************************************
  ! Calculates the mt-spheres contribution to the layer charge for states
  !  {En} at the current k-point.
  !                                      Yury Koroteev 2003
  !                     from eparas.F  by  Philipp Kurz 99/04
  !
  !***********************************************************************
  !
  SUBROUTINE calc_mt_slab(slab,itype,jsp,ikpt,atoms,ev_list,ne,abc,radfun)
    USE m_types_setup
    USE m_types_abc
    USE m_types_radfun
    IMPLICIT NONE
    TYPE(t_atoms),INTENT(IN)        :: atoms
    TYPE(t_abc),INTENT(IN)          :: abc
    CLASS(t_slab), INTENT(INOUT)     :: slab
    TYPE(t_radfun), INTENT(IN)      :: radfun
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: jsp,itype
    INTEGER, INTENT (IN) :: ne,ikpt 

    INTEGER, INTENT (IN) :: ev_list(:)

    !     ..
    !     .. Local Scalars ..
    INTEGER:: i,l,natom,m,j,jj
    INTEGER:: lm,ll1,nl
    REAL :: sabd
    

    if (.not. slab%l_initialized) return
    
    DO i = 1,ne  
      sabd = 0.0  !sum over l,m,natoms,j,jj
      DO l = 0,atoms%lmax(itype)
         ll1 = l* (l+1)
         DO m = -l,l
            lm = ll1 + m
            DO natom = 1, atoms%neq(itype)
            DO j = 1, radfun%n_r(l)
               DO jj = 1, radfun%n_r(l)
                  sabd = sabd + abc%cof(i, lm, j, natom)*CONJG(abc%cof(i, lm, jj, natom))*radfun%integral(j, jj, l, jsp, jsp)
               END DO
            END DO
            ENDDO
         enddo
      enddo
      !Map to slabs
       DO nl = 1,slab%nsl
          slab%qmtsl(nl,ev_list(i),ikpt,jsp) = sabd/atoms%neq(itype)*slab%nmtsl(itype,nl)
       ENDDO
    ENDDO
    
  END SUBROUTINE calc_mt_slab

 SUBROUTINE calc_int_slab(slab,isp,ikpt,stars,atoms,sym,cell,ne,ev_list,lapw,zMat)
    !     *******************************************************
    !     calculate the charge of the En(k) state
    !     in the interstitial region of each leyer
    !                                             Yu.M. Koroteev
    !             From pwden_old.F and pwint.F by  c.l.fu
    !     *******************************************************

    USE m_pwintsl
    USE m_types
    IMPLICIT NONE

    TYPE(t_lapw),INTENT(IN)   :: lapw
     
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_stars),INTENT(IN)  :: stars
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_mat),INTENT(IN)    :: zMat
    CLASS(t_slab),INTENT(INOUT):: slab
    !
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ne,isp,ikpt

    INTEGER, INTENT (IN) :: ev_list(ne)

    !     ..
    !     .. Local Scalars ..
    REAL q1,zsl1,zsl2,qi,volsli,volintsli
    INTEGER i ,indp,ix1,iy1,iz1,j,n,ns,ind
    COMPLEX x,phase,phasep
    !     ..
    !     .. Local Arrays ..
    COMPLEX, ALLOCATABLE :: stfunint(:,:),z_z(:)
    if (.not. slab%l_initialized) return
    !     ..
    !
    !     calculate the star function expansion coefficients of
    !     the plane wave charge density for each En(k)
    !
    !     ----> g=0 star
    !
    ALLOCATE ( stfunint(stars%ng3,slab%nsl), z_z(stars%ng3) )
    !
    !  -----> calculate the integrals of star functions over
    !                     the layer interstitial
    !
    DO i = 1,slab%nsl
       zsl1 = slab%zsl(1,i)
       zsl2 = slab%zsl(2,i)
       volsli = slab%volsl(i)
       volintsli = slab%volintsl(i)
       DO j = 1,stars%ng3
          CALL pwint_sl(stars,atoms,sym,zsl1,zsl2,&
                        volsli,volintsli,cell,slab%nmtsl(1,i),stars%kv3(1,j),x)
          stfunint(j,i) =  x*stars%nstr(j)
       ENDDO  ! over 3D stars
    ENDDO     ! over banddos%layers
    !
    ! Here, I reordered the stuff to save memory
    !
    DO  n = 1,ne
       z_z(:) = CMPLX(0.0,0.0)
       q1 = 0.0
       IF (zmat%l_real) THEN
          DO  i = 1,lapw%nv(isp)
             q1 = q1 + zMat%data_r(i,n)*zMat%data_r(i,n)
          ENDDO
       ELSE
          DO  i = 1,lapw%nv(isp)
             q1 = q1 + REAL(zMat%data_c(i,n)*CONJG(zMat%data_c(i,n)))
          ENDDO
       ENDIF
       z_z(1) = q1/cell%omtil
       !
       !     ----> g.ne.0 stars
       !
       DO  i = 1,lapw%nv(isp)
          DO  j = 1,i-1
             ix1 = lapw%k1(j,isp) - lapw%k1(i,isp)
             iy1 = lapw%k2(j,isp) - lapw%k2(i,isp)
             iz1 = lapw%k3(j,isp) - lapw%k3(i,isp)
             IF (iabs(ix1).GT.stars%mx1) CYCLE
             IF (iabs(iy1).GT.stars%mx2) CYCLE
             IF (iabs(iz1).GT.stars%mx3) CYCLE
             ind = stars%ig(ix1,iy1,iz1)
             indp = stars%ig(-ix1,-iy1,-iz1)
             IF (ind.EQ.0 .OR. indp.EQ.0) CYCLE
             phase = stars%rgphs(ix1,iy1,iz1)/ (stars%nstr(ind)*cell%omtil)
             phasep = stars%rgphs(-ix1,-iy1,-iz1)/ (stars%nstr(indp)*cell%omtil)
             IF (zmat%l_real) THEN
                z_z(ind)  = z_z(ind)  + zMat%data_r(j,n)*zMat%data_r(i,n)*REAL(phase)
                z_z(indp) = z_z(indp) + zMat%data_r(i,n)*zMat%data_r(j,n)*REAL(phasep)
             ELSE
                z_z(ind) = z_z(ind) +zMat%data_c(j,n)*CONJG(zMat%data_c(i,n))*phase
                z_z(indp)= z_z(indp)+zMat%data_c(i,n)*CONJG(zMat%data_c(j,n))*phasep
             ENDIF
          ENDDO
       ENDDO
       ! ----> calculate a charge in the layer interstitial region of the film
       !
       DO i = 1,slab%nsl
          qi = 0.0
          DO j = 1,stars%ng3
             qi = qi + z_z(j)*stfunint(j,i)
          ENDDO
          slab%qintsl(i,ev_list(n),ikpt,isp) = qi
       ENDDO    ! over banddos%layers

    ENDDO ! over states

    DEALLOCATE ( stfunint, z_z )

  END SUBROUTINE calc_int_slab


end module m_types_slab
