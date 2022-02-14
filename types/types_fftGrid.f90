!--------------------------------------------------------------------------------
! Copyright (c) 2020 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_fftGrid
   use m_constants
   TYPE t_fftGrid

      INTEGER :: extend(3) = [-1, -1, -1]
      INTEGER :: dimensions(3) = [-1, -1, -1]
      INTEGER :: gridLength = -1
      COMPLEX, ALLOCATABLE :: grid(:)

   CONTAINS

      PROCEDURE :: init => t_fftGrid_init
      PROCEDURE :: free => t_fftGrid_free
      PROCEDURE :: putFieldOnGrid
      PROCEDURE :: takeFieldFromGrid
      PROCEDURE :: getRealPartOfGrid
      PROCEDURE :: putStateOnGrid
      procedure :: put_state_on_external_grid
      PROCEDURE :: putRealStateOnGrid
      PROCEDURE :: putComplexStateOnGrid
      PROCEDURE :: fillStateIndexArray
      PROCEDURE :: fillFieldSphereIndexArray
      PROCEDURE :: getElement
      procedure :: g2fft => map_g_to_fft_grid
      PROCEDURE :: perform_fft
   END TYPE t_fftGrid

   PUBLIC :: t_fftGrid, put_real_on_external_grid, put_cmplx_on_external_grid

CONTAINS
subroutine perform_fft(grid,forward)
   use m_types_fft
   implicit none 
   CLASS(t_fftGrid), INTENT(INOUT) :: grid
   LOGICAL,INTENT(IN)              :: forward
   
   type(t_fft) :: fft

   if (size(grid%grid) .ne. product(grid%dimensions)) call juDFT_error('array bounds are inconsistent', calledby='perform_fft')
   
   call fft%init(grid%dimensions, forward)
   call fft%exec(grid%grid)
   call fft%free()
end subroutine perform_fft

function map_g_to_fft_grid(grid, g_in) result(g_idx)
   implicit none 
   CLASS(t_fftGrid), INTENT(IN)    :: grid
   integer, intent(in)             :: g_in(3)
   integer                         :: g_idx

   integer ::  shifted_g(3)

   ! the fft-grid starts at g=0, not -g_max
   ! therefore all negative g's need to be shifted
   
   shifted_g = merge(g_in + grid%dimensions, g_in, g_in < 0)

   ! map it to 1d
   g_idx = shifted_g(1) &
         + shifted_g(2) * grid%dimensions(1) &
         + shifted_g(3) * grid%dimensions(1) * grid%dimensions(2)
 end function map_g_to_fft_grid

   SUBROUTINE t_fftGrid_init(this, cell, sym, gCutoff,dims)
      USE m_constants
      USE m_boxdim
      USE m_spgrot
      USE m_ifft
      USE m_types_cell
      USE m_types_sym
      IMPLICIT NONE
      CLASS(t_fftGrid), INTENT(INOUT) :: this
      TYPE(t_cell), INTENT(IN)    :: cell
      TYPE(t_sym), INTENT(IN)    :: sym
      REAL, INTENT(IN)    :: gCutoff
      INTEGER,INTENT(IN),OPTIONAL :: dims(3)

      INTEGER, ALLOCATABLE :: ig(:, :, :)

      INTEGER :: k1, k2, k3, i
      INTEGER :: mxx(3), kVec(3), kRot(3, sym%nop), inv_du(sym%nop), tempDim(3)
      REAL    :: gCutoffSquared, gSquared
      REAL    :: arltv(3), g(3)

      this%extend     = calc_extend(cell, sym, gCutoff)
      this%dimensions = 2*this%extend + 1

      do i = 1,3
         this%dimensions(i) = next_optimal_fft_size(this%dimensions(i))
      enddo
      if (present(dims)) this%dimensions=dims
      this%gridLength = product(this%dimensions)
   
      IF (ALLOCATED(this%grid)) DEALLOCATE (this%grid)
      ALLOCATE (this%grid(0:this%gridLength - 1), source=CMPLX_NOT_INITALIZED)
   END SUBROUTINE t_fftGrid_init

   SUBROUTINE putFieldOnGrid(this, stars, cell,field, gCutoff, firstderiv,secondderiv)
      USE m_types_stars
      USE m_types_cell
      IMPLICIT NONE
      CLASS(t_fftGrid), INTENT(INOUT) :: this
      TYPE(t_stars), INTENT(IN)       :: stars
      TYPE(t_cell),INTENT(IN)         :: cell
      COMPLEX, INTENT(IN)             :: field(:) ! length is stars%ng3
      REAL, OPTIONAL, INTENT(IN)      :: gCutoff
      real,optional,intent(in)        :: firstderiv(3),secondderiv(3)

      INTEGER :: x, y, z, iStar
      INTEGER :: xGrid, yGrid, zGrid, layerDim
      REAL    :: gCutoffInternal
      COMPLEX :: fct

      gCutoffInternal = 1.0e99
      IF (PRESENT(gCutoff)) gCutoffInternal = gCutoff

      layerDim = this%dimensions(1)*this%dimensions(2)

      this%grid(:) = CMPLX(0.0, 0.0)

      DO z = -stars%mx3, stars%mx3
         zGrid = MODULO(z, this%dimensions(3))
         DO y = -stars%mx2, stars%mx2
            yGrid = MODULO(y, this%dimensions(2))
            DO x = -stars%mx1, stars%mx1
               iStar = stars%ig(x, y, z)
               fct=stars%rgphs(x, y, z)
               if (present(firstderiv)) THEN
                  fct=fct*cmplx(0.0,-1*dot_product(firstderiv,matmul(real([x,y,z]),cell%bmat)))
                  if (present(secondderiv)) fct=fct*cmplx(0.0,-1*dot_product(secondderiv,matmul(real([x,y,z]),cell%bmat)))
               endif   
               IF (iStar .EQ. 0) CYCLE
               IF (stars%sk3(iStar) .GT. gCutoffInternal) CYCLE
               xGrid = MODULO(x, this%dimensions(1))
               this%grid(xGrid + this%dimensions(1)*yGrid + layerDim*zGrid) = field(iStar)*fct
            END DO
         END DO
      END DO

   END SUBROUTINE putFieldOnGrid

   SUBROUTINE takeFieldFromGrid(this, stars, field, gCutoff)
      USE m_types_stars
      IMPLICIT NONE
      CLASS(t_fftGrid), INTENT(IN)    :: this
      TYPE(t_stars), INTENT(IN)    :: stars
      COMPLEX, INTENT(INOUT) :: field(:)
      REAL, OPTIONAL, INTENT(IN)    :: gCutoff

      INTEGER :: x, y, z, iStar
      INTEGER :: xGrid, yGrid, zGrid, layerDim
      REAL    :: elementWeight, gCutoffInternal

      gCutoffInternal = 1.0e99
      IF (PRESENT(gCutoff)) gCutoffInternal = gCutoff

      field(:) = CMPLX(0.0, 0.0)
      layerDim = this%dimensions(1)*this%dimensions(2)

      DO z = -stars%mx3, stars%mx3
         zGrid = MODULO(z, this%dimensions(3))
         DO y = -stars%mx2, stars%mx2
            yGrid = MODULO(y, this%dimensions(2))
            DO x = -stars%mx1, stars%mx1
               iStar = stars%ig(x, y, z)
               IF (iStar .EQ. 0) CYCLE
               IF (stars%sk3(iStar) .GT. gCutoffInternal) CYCLE
               xGrid = MODULO(x, this%dimensions(1))
               field(iStar) = field(iStar) + this%grid(xGrid + this%dimensions(1)*yGrid + layerDim*zGrid) * CONJG(stars%rgphs(x, y, z))
            END DO
         END DO
      END DO

      elementWeight = 1.0/(this%dimensions(1)*this%dimensions(2)*this%dimensions(3))

      field(:) = elementWeight*field(:)/stars%nstr(:)

   END SUBROUTINE takeFieldFromGrid

   SUBROUTINE putStateOnGrid(this, lapw, iSpin, zMat, iState)
      USE m_types_lapw
      USE m_types_mat
      IMPLICIT NONE
      CLASS(t_fftGrid), INTENT(INOUT) :: this
      TYPE(t_lapw), INTENT(IN)    :: lapw
      TYPE(t_mat), INTENT(IN)    :: zMat
      INTEGER, INTENT(IN)    :: iSpin
      INTEGER, INTENT(IN)    :: iState

      IF (zMat%l_real) THEN
         CALL putRealStateOnGrid(this, lapw, iSpin, zMat%data_r(:, iState))
      ELSE
         CALL putComplexStateOnGrid(this, lapw, iSpin, zMat%data_c(:, iState))
      END IF
   END SUBROUTINE putStateOnGrid

   SUBROUTINE put_state_on_external_grid(this, lapw, iSpin, zMat, iState, ext_grid, l_gpu)
      USE m_types_lapw
      USE m_types_mat
      IMPLICIT NONE
      CLASS(t_fftGrid), INTENT(INOUT) :: this
      TYPE(t_lapw), INTENT(IN)    :: lapw
      TYPE(t_mat), INTENT(IN)    :: zMat
      INTEGER, INTENT(IN)    :: iSpin
      INTEGER, INTENT(IN)    :: iState
      complex, intent(inout) :: ext_grid(0:)
      logical, intent(in), optional :: l_gpu

      if (zMat%l_real) then
         call put_real_on_external_grid(this, lapw, ispin, zMat%data_r(:, iState), ext_grid, l_gpu)
      else
         call put_cmplx_on_external_grid(this, lapw, ispin, zMat%data_c(:, iState), ext_grid, l_gpu)
      endif
   end subroutine put_state_on_external_grid

   SUBROUTINE putRealStateOnGrid(this, lapw, iSpin, state)
      USE m_types_lapw
      IMPLICIT NONE
      CLASS(t_fftGrid), INTENT(INOUT) :: this
      TYPE(t_lapw), INTENT(IN)    :: lapw
      REAL, INTENT(IN)    :: state(:)
      INTEGER, INTENT(IN)    :: iSpin

      call put_real_on_external_grid(this, lapw, ispin, state, this%grid)
   END SUBROUTINE putRealStateOnGrid

   subroutine put_real_on_external_grid(this, lapw, ispin, state, ext_grid, l_gpu)   
     USE m_types_lapw
      IMPLICIT NONE
      CLASS(t_fftGrid), INTENT(INOUT) :: this
      TYPE(t_lapw), INTENT(IN)    :: lapw
      REAL, INTENT(IN)       :: state(:)
      complex, intent(inout) :: ext_grid(0:)
      INTEGER, INTENT(IN)    :: iSpin
      logical, intent(in), optional :: l_gpu

      logical :: use_gpu
      INTEGER :: xGrid, yGrid, zGrid, layerDim, iLAPW

      layerDim = this%dimensions(1)*this%dimensions(2)

      if(present(l_gpu)) then 
         use_gpu = l_gpu 
      else
         use_gpu = .False. 
      endif 

      if(use_gpu) then
         !$acc kernels
         ext_grid = cmplx_0
         !$acc end kernels

         !$acc parallel loop default(none) private(xGrid, yGrid, zGrid) &
         !$acc present(lapw, lapw%nv, lapw%gvec, this, this%dimensions, ext_grid, state) &
         !$acc copyin(layerDim, iSpin)
         DO iLAPW = 1, lapw%nv(iSpin)
            xGrid = MODULO(lapw%gvec(1, iLAPW, iSpin), this%dimensions(1))
            yGrid = MODULO(lapw%gvec(2, iLAPW, iSpin), this%dimensions(2))
            zGrid = MODULO(lapw%gvec(3, iLAPW, iSpin), this%dimensions(3))
            ext_grid(xGrid + this%dimensions(1)*yGrid + layerDim*zGrid) = state(iLAPW)
         END DO
         !$acc end parallel loop
      else
         ext_grid = cmplx_0

         DO iLAPW = 1, lapw%nv(iSpin)
            xGrid = MODULO(lapw%gvec(1, iLAPW, iSpin), this%dimensions(1))
            yGrid = MODULO(lapw%gvec(2, iLAPW, iSpin), this%dimensions(2))
            zGrid = MODULO(lapw%gvec(3, iLAPW, iSpin), this%dimensions(3))
            ext_grid(xGrid + this%dimensions(1)*yGrid + layerDim*zGrid) = state(iLAPW)
         END DO
      endif
   end subroutine put_real_on_external_grid

   SUBROUTINE putComplexStateOnGrid(this, lapw, iSpin, state)
      USE m_types_lapw
      IMPLICIT NONE
      CLASS(t_fftGrid), INTENT(INOUT) :: this
      TYPE(t_lapw), INTENT(IN)    :: lapw
      COMPLEX, INTENT(IN)    :: state(:)
      INTEGER, INTENT(IN)    :: iSpin

      call put_cmplx_on_external_grid(this, lapw, ispin, state, this%grid)
   END SUBROUTINE putComplexStateOnGrid

   SUBROUTINE put_cmplx_on_external_grid(this, lapw, iSpin, state, ext_grid, l_gpu)
      USE m_types_lapw
      use m_judft
      IMPLICIT NONE
      CLASS(t_fftGrid), INTENT(INOUT) :: this
      TYPE(t_lapw), INTENT(IN)    :: lapw
      COMPLEX, INTENT(IN)    :: state(:)
      complex, intent(inout) :: ext_grid(0:)
      INTEGER, INTENT(IN)    :: iSpin
      logical, intent(in), optional :: l_gpu

      logical :: use_gpu
      INTEGER :: xGrid, yGrid, zGrid, iLAPW

      if(present(l_gpu)) then 
         use_gpu = l_gpu 
      else
         use_gpu = .False. 
      endif 

      if(use_gpu) then
         !$acc kernels
         ext_grid = cmplx_0
         !$acc end kernels

         !$acc parallel loop default(none) private(xGrid, yGrid, zGrid) &
         !$acc present(lapw, lapw%nv, lapw%gvec, this, this%dimensions, ext_grid, state, iSpin) 
         DO iLAPW = 1, lapw%nv(iSpin)
            xGrid = MODULO(lapw%gvec(1, iLAPW, iSpin), this%dimensions(1))
            yGrid = MODULO(lapw%gvec(2, iLAPW, iSpin), this%dimensions(2))
            zGrid = MODULO(lapw%gvec(3, iLAPW, iSpin), this%dimensions(3))
            ext_grid(xGrid + this%dimensions(1)*yGrid + (this%dimensions(1)*this%dimensions(2))*zGrid) = state(iLAPW)
         END DO
         !$acc end parallel loop
      else
         ext_grid = cmplx_0

         DO iLAPW = 1, lapw%nv(iSpin)
            xGrid = MODULO(lapw%gvec(1, iLAPW, iSpin), this%dimensions(1))
            yGrid = MODULO(lapw%gvec(2, iLAPW, iSpin), this%dimensions(2))
            zGrid = MODULO(lapw%gvec(3, iLAPW, iSpin), this%dimensions(3))
            ext_grid(xGrid + this%dimensions(1)*yGrid + (this%dimensions(1)*this%dimensions(2))*zGrid) = state(iLAPW)
         END DO
      endif
   end SUBROUTINE put_cmplx_on_external_grid

   SUBROUTINE fillStateIndexArray(this, lapw, ispin, indexArray)
      USE m_types_lapw
      IMPLICIT NONE
      CLASS(t_fftGrid), INTENT(INOUT) :: this
      TYPE(t_lapw), INTENT(IN)        :: lapw
      INTEGER, INTENT(IN)             :: iSpin
      INTEGER, INTENT(INOUT)          :: indexArray(lapw%nv(ispin))

      INTEGER :: xGrid, yGrid, zGrid, layerDim, iLAPW

      layerDim = this%dimensions(1)*this%dimensions(2)

      DO iLAPW = 1, lapw%nv(iSpin)
         xGrid = MODULO(lapw%gvec(1, iLAPW, iSpin), this%dimensions(1))
         yGrid = MODULO(lapw%gvec(2, iLAPW, iSpin), this%dimensions(2))
         zGrid = MODULO(lapw%gvec(3, iLAPW, iSpin), this%dimensions(3))
         indexArray(iLAPW) = xGrid + this%dimensions(1)*yGrid + layerDim*zGrid
      END DO
   END SUBROUTINE fillStateIndexArray

   SUBROUTINE fillFieldSphereIndexArray(this, stars, gCutoff, indexArray)
      USE m_types_stars
      IMPLICIT NONE
      CLASS(t_fftGrid), INTENT(IN)        :: this
      TYPE(t_stars), INTENT(IN)           :: stars
      REAL, INTENT(IN)                    :: gCutoff
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: indexArray(:)

      INTEGER :: x, y, z, iStar
      INTEGER :: xGrid, yGrid, zGrid, layerDim, tempArrayIndex
      INTEGER :: tempArray((2*stars%mx1+1)*(2*stars%mx2+1)*(2*stars%mx3+1))

      layerDim = this%dimensions(1)*this%dimensions(2)

      tempArrayIndex = 0
      DO z = -stars%mx3, stars%mx3
         zGrid = MODULO(z, this%dimensions(3))
         DO y = -stars%mx2, stars%mx2
            yGrid = MODULO(y, this%dimensions(2))
            DO x = -stars%mx1, stars%mx1
               iStar = stars%ig(x, y, z)
               IF (iStar .EQ. 0) CYCLE
               IF (stars%sk3(iStar) .GT. gCutoff) CYCLE
               xGrid = MODULO(x, this%dimensions(1))
               tempArrayIndex = tempArrayIndex + 1
               tempArray(tempArrayIndex) = xGrid + this%dimensions(1)*yGrid + layerDim*zGrid
            END DO
         END DO
      END DO

      IF(ALLOCATED(indexArray)) DEALLOCATE (indexArray)
      ALLOCATE(indexArray(tempArrayIndex))

      indexArray(1:tempArrayIndex) = tempArray(1:tempArrayIndex)

   END SUBROUTINE fillFieldSphereIndexArray

   COMPLEX FUNCTION getElement(this, x, y, z)
      IMPLICIT NONE
      CLASS(t_fftGrid), INTENT(IN)    :: this
      INTEGER, INTENT(IN)    :: x, y, z

      INTEGER :: xGrid, yGrid, zGrid, layerDim

      layerDim = this%dimensions(1)*this%dimensions(2)

      xGrid = MODULO(x, this%dimensions(1))
      yGrid = MODULO(y, this%dimensions(2))
      zGrid = MODULO(z, this%dimensions(3))

      getElement = this%grid(xGrid + this%dimensions(1)*yGrid + layerDim*zGrid)

   END FUNCTION getElement

   SUBROUTINE getRealPartOfGrid(this, realGrid)
      IMPLICIT NONE
      CLASS(t_fftGrid), INTENT(IN)    :: this
      REAL, INTENT(INOUT) :: realGrid(:)

      realGrid(:) = REAL(this%grid(:))

   END SUBROUTINE getRealPartOfGrid

   subroutine t_fftGrid_free(fftGrid)
      implicit none
      CLASS(t_fftGrid), INTENT(INOUT)    :: fftGrid

      fftGrid%extend     = -1
      fftGrid%dimensions = -1
      fftGrid%gridLength = -1

      if(allocated(fftGrid%grid)) deallocate(fftGrid%grid)
   end subroutine t_fftGrid_free

   function calc_extend(cell, sym, gCutoff) result(mxx)
      USE m_constants
      USE m_boxdim
      USE m_spgrot
      USE m_ifft
      USE m_types_cell
      USE m_types_sym
      IMPLICIT NONE
      
      TYPE(t_cell), INTENT(IN)  :: cell
      TYPE(t_sym), INTENT(IN)   :: sym
      REAL, INTENT(IN)          :: gCutoff

      INTEGER, ALLOCATABLE :: ig(:, :, :)

      INTEGER :: k1, k2, k3, i, j
      INTEGER :: mxx(3), kVec(3), kRot(3, sym%nop), inv_du(sym%nop), tempDim(3)
      REAL    :: gCutoffSquared, gSquared
      REAL    :: arltv(3), g(3)

      CALL boxdim(cell%bmat, arltv(1), arltv(2), arltv(3))

      tempDim(:) = INT(gCutoff/arltv(:)) + 1

      DO i = 1, sym%nop
         inv_du(i) = i ! dummy array for spgrot
      END DO

      ALLOCATE (ig(-tempDim(1):tempDim(1), &
                   -tempDim(2):tempDim(2), &
                   -tempDim(3):tempDim(3)), source=0)

      mxx(:) = 0
      gCutoffSquared = gCutoff*gCutoff
      DO k1 = tempDim(1), -tempDim(1), -1
         kVec(1) = k1
         DO k2 = tempDim(2), -tempDim(2), -1
            kVec(2) = k2
            DO k3 = tempDim(3), -tempDim(3), -1
               IF (ig(k1, k2, k3) .NE. 0) CYCLE

               kVec(3) = k3

               DO i = 1, 3
                  g(i) = DOT_PRODUCT(cell%bmat(:, i), kVec(:))  ! loop to be replaced by MATMUL call later.
               END DO

               gSquared = g(1)**2 + g(2)**2 + g(3)**2

               IF (gSquared .LE. gCutoffSquared) THEN
                  CALL spgrot(sym%nop, .true., sym%mrot, sym%tau, inv_du, kVec, kRot)
                  DO i = 1, sym%nop
                     do j = 1,3 
                        mxx(j) = max(mxx(j), kRot(j,i))
                     enddo
                     ig(kRot(1, i), kRot(2, i), kRot(3, i)) = 1
                  END DO
               END IF
            END DO
         END DO
      END DO
   END function calc_extend

   function calc_fft_dim(cell, sym, gCutoff) result(dims)
      USE m_ifft
      USE m_types_cell
      USE m_types_sym
      implicit none
      TYPE(t_cell), INTENT(IN)  :: cell
      TYPE(t_sym), INTENT(IN)   :: sym
      REAL, INTENT(IN)          :: gCutoff 
      integer :: dims(3)
      integer :: i

      dims = 2* calc_extend(cell, sym, gCutoff) + 1

      do i = 1,3
         dims(i) = next_optimal_fft_size(dims(i))
      enddo
   end function calc_fft_dim

END MODULE m_types_fftGrid
