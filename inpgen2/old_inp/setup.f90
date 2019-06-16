      MODULE m_setup
      USE m_juDFT
      CONTAINS
        SUBROUTINE setup(mpi,atoms,kpts,DIMENSION,sphhar,&
             obsolete,sym,stars,oneD, input,noco,vacuum,cell,xcpot, sliceplot,enpara)
          !stripped down version 

          
          !
          !----------------------------------------
          ! this routine is called by: fleur.F
          !
          ! setup --+
          !         +-- spg2set
          !         +-- local_sym -+- ptsym
          !         |              +- lhcal -+- gaussp
          !         |                        +- gtest
          !         |                        +- ylm4
          !         +-- strgn1 -+- boxdim
          !         |           +- sort
          !         |           +- dset
          !         |           +- spgrot
          !         |           +- unor2or
          !         +-- mapatom -+- dotset
          !         |            +- dotirl
          !         +-- inpeig -- gkptwgt
          !         +-- gaunt2 -- grule
          !         +-- prp_qfft -+- boxdim
          !         |             +- ifft235
          !         +-- prp_xcfft -+- boxdim
          !         |              +- ifft235
          !         +-- stepf
          !         +-- convn
          !         +-- efield
          !----------------------------------------

          !
          USE m_types
          USE m_localsym
          USE m_rwsymfile
          USE m_spg2set
          USE m_dwigner
          USE m_strgn
          USE m_stepf
          USE m_cdn_io
          USE m_mapatom
          USE m_convn
          USE m_prpqfft
          USE m_prpxcfft
          USE m_inpeig
          USE m_efield
          USE m_ylm
          !-odim
          USE m_od_mapatom
          USE m_od_chisym
          USE m_od_strgn1
          !+odim
          IMPLICIT NONE
          !     ..
          !     .. Scalars Arguments ..
          TYPE(t_mpi),INTENT(IN)         :: mpi
          TYPE(t_atoms),INTENT(INOUT)    :: atoms
          TYPE(t_kpts),INTENT(INOUT)     :: kpts
          TYPE(t_dimension),INTENT(INOUT):: DIMENSION
          TYPE(t_sphhar),INTENT(INOUT)   :: sphhar
          TYPE(t_obsolete),INTENT(INOUT) :: obsolete
          TYPE(t_sym),INTENT(INOUT)      :: sym
          TYPE(t_stars),INTENT(INOUT)    :: stars
          TYPE(t_oneD),INTENT(INOUT)     :: oneD
          TYPE(t_input),INTENT(INOUT)    :: input
          TYPE(t_noco),INTENT(INOUT)     :: noco
          TYPE(t_vacuum),INTENT(INOUT)   :: vacuum
          TYPE(t_cell),INTENT(INOUT)     :: cell
          CLASS(t_xcpot),INTENT(INOUT)   :: xcpot
          TYPE(t_sliceplot),INTENT(INOUT):: sliceplot
          TYPE(t_enpara),INTENT(INOUT)   :: enpara
          !     ..
          !     .. Local Scalars ..
          REAL       :: rkmaxx
          INTEGER    :: ntp1,ii,i,j,n1,n2,na,np1,n
          INTEGER, ALLOCATABLE :: lmx1(:), nq1(:), nlhtp1(:)
          !
        IF ( mpi%irank == 0 ) THEN
          IF (sym%namgrp.EQ.'any ') THEN
             CALL rw_symfile('R',94,'sym.out',sym%nop,cell%bmat, sym%mrot,sym%tau,sym%nop,sym%nop2,sym%symor)
          ELSE
             CALL spg2set(sym%nop,sym%zrfs,sym%invs,sym%namgrp,cell%latnam, sym%mrot,sym%tau,sym%nop2,sym%symor)
          ENDIF
          IF (input%film.AND..NOT.sym%symor) CALL juDFT_warn("Films&Symor",hint&
               &     ="Films should be symmorphic",calledby ='setup')


          !
          CALL inpeig(atoms,cell,input,oneD%odd%d1,kpts,enpara)
          !
          !
        ENDIF ! (mpi%irank == 0)

        END SUBROUTINE setup
      END MODULE m_setup
