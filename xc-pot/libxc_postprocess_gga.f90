!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_libxc_postprocess_gga
CONTAINS

   SUBROUTINE libxc_postprocess_gga_mt(xcpot,atoms,sphhar,n,v_xc,grad, atom_num)
      USE m_mt_tofrom_grid
      USE m_types
      USE m_npy
      use m_judft_string

      IMPLICIT NONE
      CLASS(t_xcpot),INTENT(IN)   :: xcpot
      TYPE(t_atoms),INTENT(IN)    :: atoms
      TYPE(t_sphhar),INTENT(IN)   :: sphhar
      INTEGER,INTENT(IN)          :: n
      REAL,INTENT(INOUT)          :: v_xc(:,:)
      TYPE(t_gradients),INTENT(IN):: grad
      INTEGER, OPTIONAL           :: atom_num

      INTEGER :: nsp,n_sigma,i
      REAL,ALLOCATABLE:: vsigma(:,:),vsigma_mt(:,:,:)
      TYPE(t_gradients)::grad_vsigma
      character(len=:), allocatable :: fname

      n_sigma=MERGE(1,3,SIZE(v_xc,2)==1) !Number of contracted gradients in libxc 1 for non-spin-polarized, 3 otherwise
      nsp=SIZE(v_xc,1) !no of points
      ALLOCATE(vsigma(nsp,n_sigma),vsigma_mt(atoms%jri(n),0:sphhar%nlhd,n_sigma))
      vsigma_mt=0.0
      vsigma=TRANSPOSE(grad%vsigma) !create a (nsp,n_sigma) matrix
      CALL mt_from_grid(atoms,sphhar,n,n_sigma,vsigma,vsigma_mt)
      DO i=1,atoms%jri(n)
         vsigma_mt(i,:,:)=vsigma_mt(i,:,:)*atoms%rmsh(i,n)**2
      ENDDO
      ALLOCATE(grad_vsigma%gr(3,nsp,n_sigma))
      CALL mt_to_grid(xcpot,n_sigma,atoms,sphhar,vsigma_mt,n,grad=grad_vsigma)

      fname = merge("mt=" //int2str(atom_num) // "_lapl.npy","mt_lapl.npy", present(atom_num))
      call save_npy(fname, grad%laplace)

      fname = merge("mt="// int2str(atom_num) // "_grad.npy", "mt_grad.npy", present(atom_num))
      call save_npy(fname, grad%gr) 
      CALL libxc_postprocess_gga(transpose(grad%vsigma),grad,grad_vsigma,v_xc)
   END SUBROUTINE libxc_postprocess_gga_mt

   SUBROUTINE libxc_postprocess_gga_pw(xcpot,stars,cell,v_xc,grad)
      USE m_pw_tofrom_grid
      USE m_types
      USE m_npy

      IMPLICIT NONE
      CLASS(t_xcpot),INTENT(IN)   :: xcpot
      TYPE(t_stars),INTENT(IN)    :: stars
      TYPE(t_cell),INTENT(IN)     :: cell
      REAL,INTENT(INOUT)          :: v_xc(:,:)
      TYPE(t_gradients),INTENT(IN):: grad

      COMPLEX,ALLOCATABLE:: vsigma_g(:,:)
      REAL,ALLOCATABLE:: vsigma(:,:)
      TYPE(t_gradients)::grad_vsigma
      INTEGER :: nsp,n_sigma

      nsp=SIZE(v_xc,1) !no of points
      n_sigma=MERGE(1,3,SIZE(v_xc,2)==1) !See in _mt routine
      ALLOCATE(vsigma_g(stars%ng3,n_sigma),vsigma(nsp,n_sigma)); vsigma_g=0.0
      vsigma=TRANSPOSE(grad%vsigma) !create a (nsp,n_sigma) matrix
      CALL pw_from_grid(xcpot,stars,.FALSE.,vsigma,vsigma_g)
      !vsigma_g(:,1)=vsigma_g(:,1)*stars%nstr(:)
      ALLOCATE(grad_vsigma%gr(3,nsp,n_sigma))
      CALL pw_to_grid(xcpot,n_sigma,.false.,stars,cell,vsigma_g,grad_vsigma)

      call save_npy("pw_lapl.npy", grad%laplace)
      call save_npy("pw_grad.npy", grad%gr) 
      CALL libxc_postprocess_gga(transpose(grad%vsigma),grad,grad_vsigma,v_xc)

   END SUBROUTINE libxc_postprocess_gga_pw

   SUBROUTINE libxc_postprocess_gga(vsigma,grad,grad_vsigma,v_xc)
      USE m_types
      IMPLICIT NONE
      REAL,INTENT(IN)             :: vsigma(:,:)
      TYPE(t_gradients),INTENT(IN):: grad,grad_vsigma
      REAL,INTENT(INOUT)          :: v_xc(:,:)
      INTEGER:: i
      IF (SIZE(v_xc,2)==1) THEN !Single spin
         DO i=1,SIZE(v_xc,1) !loop over points
            v_xc(i,1)=v_xc(i,1)-2*dot_PRODUCT(grad_vsigma%gr(:,i,1),grad%gr(:,i,1))-2*vsigma(i,1)*grad%laplace(i,1)
         ENDDO
      ELSE  !two spins
         DO i=1,SIZE(v_xc,1) !loop over points
            v_xc(i,1)=v_xc(i,1)-2*dot_PRODUCT(grad_vsigma%gr(:,i,1),grad%gr(:,i,1))-2*vsigma(i,1)*grad%laplace(i,1)-&
                       dot_PRODUCT(grad_vsigma%gr(:,i,2),grad%gr(:,i,2))-vsigma(i,2)*grad%laplace(i,2)
            v_xc(i,2)=v_xc(i,2)-2*dot_PRODUCT(grad_vsigma%gr(:,i,3),grad%gr(:,i,2))-2*vsigma(i,3)*grad%laplace(i,2)-&
                       dot_PRODUCT(grad_vsigma%gr(:,i,2),grad%gr(:,i,1))-vsigma(i,2)*grad%laplace(i,1)
         ENDDO
      END IF

   END SUBROUTINE libxc_postprocess_gga

END MODULE m_libxc_postprocess_gga
