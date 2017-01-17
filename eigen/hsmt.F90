MODULE m_hsmt
  USE m_juDFT
  IMPLICIT NONE
CONTAINS
  SUBROUTINE hsmt(DIMENSION,atoms,sphhar,sym,enpara,SUB_COMM,n_size,n_rank,&
       jsp,input,mpi,lmaxb,gwc,noco,cell,&
       lapw, bkpt, vr,vs_mmp, oneD,usdus, kveclo,tlmplm,l_real,hamOvlp)
    !=====================================================================
    !     redesign of old hssphn into hsmt
    !           splitted into several parts:
    !                ->hsmt_init_soc
    !                ->hsmt_sph (spherical part)
    !                ->hsmt_extra (LOs and LDA+U)
    !                ->hsmt_nonsp (non-spherical part)
    !
    !                     Daniel Wortmann 2014
    !=====================================================================
    !*********************************************************************
    !     updates the hamiltonian and overlap matrices with the
    !     contributions from the spheres, both spherical and non-
    !     spherical.
    !                m. weinert  1986
    !     modified for inversion symmetry (real version), rp.p nov. 92
    !
    !     Since this non-collinear version allows only one definite
    !     direction of magnetic field in each muffin-tin, no matrix
    !     potential has to be used inside the MT's. However, the boundary
    !     conditions change, because the local (spin-) frame is different
    !     from the global frame. As a result, the usual spin-up and -down
    !     Hamitonian- and overlapp-matrix elements have to be
    !     premultiplied by a phasefactor, depending on the direction of
    !     the magentic field of each MT, to form the elements of the full
    !     complex matrices.
    !
    !     Philipp Kurz 98/01/27
    !*********************************************************************
    !------------------------------------------------------------------+
    ! Note for ev-parallelization:                                     |
    !                             unlike in hsint and hsvac, we do not |
    ! move through the H-matrix block by block (i.e. first the up/up,  |
    ! then the down/down, and finally the up/down spin blocks) but     |
    ! go  through one block only and update all blocks simultaniously. |
    ! Therefore, virtually we move through the whole matrix (up to     |
    ! nv(1)+nv(2)) and project back to the first block.                |
    !                                                   gb-00          |
    !------------------------------------------------------------------+
    !******** ABBREVIATIONS ***********************************************
    !     aa       : hamitonian matrix
    !     bb       : overlapp  matrix
    !     alph,beta: Euler angles of the local magnetic field direction of
    !                each atom (-type).
    !     chi      : Pauli spinors of the local spin-coordinate-frames
    !     chinn    : prefactors to be multiplied with the Hamiltonian- and
    !                overlappmatrix elements
    !**********************************************************************
    !
    !******** Spin-orbit interaction *************************************
    !   When l_soc=true & l_noco=true & l_ss=false, spin-orbit interaction
    !   is added in the first variation
    !   The formula by Youn et al. J.Comp.Phys. 172, 387 (2001) is used
    !   (note the wrong sign in the paper):
    !   <G|Hso|G'>= -i sigma . (G x G') sum_l (2l+1)/4pi Rso P'_l(G.G')
    !
    !   Jussi Enkovaara 2004, Juelich
    !********************************************************************
#include"cpp_double.h"
    USE m_hsmt_socinit
    USE m_hsmt_nonsph
    USE m_hsmt_sph
    USE m_hsmt_extra
    USE m_types
    USE m_hsmt_fjgj
    USE m_hsmt_overlap
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)        :: mpi
    TYPE(t_dimension),INTENT(IN)  :: DIMENSION
    TYPE(t_oneD),INTENT(IN)       :: oneD
    TYPE(t_input),INTENT(IN)      :: input
    TYPE(t_noco),INTENT(IN)       :: noco
    TYPE(t_sym),INTENT(IN)        :: sym
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_sphhar),INTENT(IN)     :: sphhar
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_enpara),INTENT(IN)     :: enpara
    TYPE(t_lapw),INTENT(INOUT)    :: lapw !nlotot,nv_tot&nmat can be updated
    TYPE(t_hamOvlp),INTENT(INOUT) :: hamOvlp
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: SUB_COMM,n_size,n_rank 
    INTEGER, INTENT (IN) :: jsp  
    INTEGER, INTENT (IN) :: lmaxb
    INTEGER, INTENT (IN) :: gwc

    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: bkpt(3) 
    REAL,    INTENT (IN) :: vr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,DIMENSION%jspd)
    COMPLEX,INTENT(IN):: vs_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb,atoms%n_u,input%jspins)
    TYPE(t_usdus),INTENT(INOUT)  :: usdus

    INTEGER, INTENT (OUT) :: kveclo(atoms%nlotot)

    LOGICAL,INTENT(IN)    :: l_real

    TYPE(t_tlmplm),INTENT(INOUT) :: tlmplm


    !     ..
    !     .. Local Scalars ..
    INTEGER :: k,i,hlpmsize,isp,jsp_start,jsp_end,nintsp,iintsp,nc,ab_dim
    LOGICAL :: l_socfirst
    !     ..
    !     .. Local Arrays ..
    REAL                 :: v(3)
    COMPLEX              :: isigma(2,2,3)
    REAL, ALLOCATABLE    :: fj(:,:,:,:),gj(:,:,:,:)
    REAL, ALLOCATABLE    :: gk(:,:,:),vk(:,:,:)
    REAL, PARAMETER      :: eps = 1.0e-30

    TYPE(t_rsoc):: rsoc

#if 1==2
    REAL, ALLOCATABLE    :: fj_test(:,:,:,:),gj_test(:,:,:,:)
    COMPLEX,ALLOCATABLE  :: bb_test(:)
#endif

    !
    CALL timestart("hsmt init")
    l_socfirst= noco%l_soc .AND. noco%l_noco .AND. (.NOT. noco%l_ss)
    IF (l_socfirst) CALL hsmt_socinit(mpi,atoms,sphhar,enpara,input,vr,noco,& !in
         rsoc,usdus,isigma) !out

#ifdef CPP_MPI
    !
    ! determine size of close-packed matrix
    !
    nc = 0
    k = 0
    DO  i=1+n_rank, lapw%nv(1)+atoms%nlotot, n_size
       nc = nc + 1
       k = k + n_size*(nc-1) + n_rank + 1
    ENDDO
    hlpmsize = k
#else
    hlpmsize = (DIMENSION%nvd+atoms%nlotot)*(DIMENSION%nvd+atoms%nlotot+1)/2
#endif


    !!---> pk non-collinear
    !---> loop over the local spins and interstitial spins
    nintsp = 1
    IF (noco%l_noco) THEN
       jsp_start = 1
       jsp_end   = 2
       !--->    in a spin-spiral calculations the a- and b-coeff and other
       !--->    quantities depend on the interstitial spin. therefore, an
       !--->    additional loop over the interstitial spin is needed in some
       !--->    places.
       IF (noco%l_ss) nintsp = 2
    ELSE
       jsp_start = jsp
       jsp_end   = jsp
    ENDIF
    !Set up the k+G+qss vectors
    ab_dim = 1
    IF (noco%l_ss) ab_dim = 2
    ALLOCATE(vk(DIMENSION%nvd,3,ab_dim),gk(DIMENSION%nvd,3,ab_dim))

    DO iintsp = 1,nintsp
       DO k = 1,lapw%nv(iintsp)
          IF (iintsp ==1) THEN
             v=bkpt+(/lapw%k1(k,iintsp),lapw%k2(k,iintsp),lapw%k3(k,iintsp)/)-noco%qss/2
          ELSE
             v=bkpt+(/lapw%k1(k,iintsp),lapw%k2(k,iintsp),lapw%k3(k,iintsp)/)+noco%qss/2
          ENDIF
          vk(k,:,iintsp) = v
          gk(k,:,iintsp) = MATMUL(TRANSPOSE(cell%bmat),v)/MAX (lapw%rk(k,iintsp),eps)
       END DO
    ENDDO
    !! the fj and gj arrays are constructed in hssphn_sph and used later
    IF (noco%l_constr.OR.l_socfirst) THEN
       ALLOCATE ( fj(DIMENSION%nvd,0:atoms%lmaxd,atoms%ntype,2),gj(DIMENSION%nvd,0:atoms%lmaxd,atoms%ntype,2))
    ELSE
       ALLOCATE(fj(DIMENSION%nvd,0:atoms%lmaxd,atoms%ntype,ab_dim))
       ALLOCATE(gj(DIMENSION%nvd,0:atoms%lmaxd,atoms%ntype,ab_dim))
    ENDIF
    CALL timestop("hsmt init")

    DO isp = jsp_start,jsp_end
#if 1==2
       ALLOCATE(bb_test(size(bb)))
       bb_test=0.0
       bb_test=bb
       fj=0.0
       gj=0.0
#endif
       CALL timestart("hsmt spherical")
       CALL hsmt_sph(sym,DIMENSION,atoms,SUB_COMM,n_size,n_rank,sphhar,isp,ab_dim,&
            input,hlpmsize,noco,l_socfirst,cell,nintsp,lapw,enpara%el0,usdus,&
            vr,gk,rsoc,isigma,fj,gj,l_real,hamOvlp)
       CALL timestop("hsmt spherical")
#if 1==2
       ALLOCATE(fj_test(size(fj,1),0:size(fj,2)-1,size(fj,3),size(fj,4)))
       ALLOCATE(gj_test(size(fj,1),0:size(fj,2)-1,size(fj,3),size(fj,4)))
       fj_test=0.0
       gj_test=0.0
       CALL timestart("hsmt_fjgj")
       CALL hsmt_fjgj(atoms,isp,noco,l_socfirst,cell,nintsp, lapw,usdus,fj_test,gj_test)
       CALL timestop("hsmt_fjgj")
       print *,"fj diff",maxval(abs(fj-fj_test))
       print *,"gj diff",maxval(abs(gj-gj_test))
       deallocate(fj_test,gj_test)
       bb=0.0
       bb_test=0.0
       CALL timestart("hsmt_overlap")
       CALL  hsmt_overlap(input,atoms,n_size,n_rank,isp,l_socfirst,hlpmsize,ab_dim,&
       noco,cell,nintsp, lapw,usdus,gk,fj,gj,bb_test)
       CALL timestop("hsmt_overlap")
       CALL timestart("hsmt_overlap_zherk")
       CALL  hsmt_overlap_zherk(input,sym,atoms,n_size,n_rank,isp,l_socfirst,hlpmsize,ab_dim,&
       noco,cell,nintsp, lapw,usdus,gk,vk,fj,gj,bb)
       CALL timestop("hsmt_overlap_zherk")
       print *,"bb_diff",maxval(abs(bb-bb_test))
       
       deallocate(bb_test)
#endif
       IF (.NOT.input%secvar) THEN
          CALL timestart("hsmt extra")
          IF (ANY(atoms%nlo>0).OR.ANY(atoms%lda_u%l.GE.0)) &
               CALL hsmt_extra(DIMENSION,atoms,sym,isp,n_size,n_rank,input,nintsp,sub_comm,&
               hlpmsize,lmaxb,gwc,noco,l_socfirst,lapw,cell,enpara%el0,&
               fj,gj,gk,vk,tlmplm,usdus, vs_mmp,oneD,& !in
               kveclo,l_real,hamOvlp%a_r,hamOvlp%b_r,hamOvlp%a_c,hamOvlp%b_c) !out/in
          CALL timestop("hsmt extra")
          CALL timestart("hsmt non-spherical")
          CALL hsmt_nonsph(DIMENSION,atoms,sym,SUB_COMM,n_size,n_rank,input,isp,nintsp,&
               hlpmsize,noco,l_socfirst,lapw,cell,tlmplm,fj,gj,gk,vk,oneD,l_real,hamOvlp%a_r,hamOvlp%a_c)

          CALL timestop("hsmt non-spherical")
       ENDIF
    ENDDO
#if 1==2
       print *,"hsmt checkmode"
       call judft_end("test",0)
#endif
    RETURN
  END SUBROUTINE hsmt
END MODULE m_hsmt
