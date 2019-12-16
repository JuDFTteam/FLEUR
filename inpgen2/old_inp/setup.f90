      MODULE m_setup
      USE m_juDFT
      CONTAINS
        SUBROUTINE setup(atoms,kpts,&
             sym,oneD, input,cell,enpara,latnam,namgrp)
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
          USE m_types_atoms
          USE m_types_kpts
          USE m_types_sym
          USE m_types_oned
          USE m_types_input
          USE m_types_cell
          USE m_types_enpara
     
          !USE m_localsym
          USE m_rwsymfile
          USE m_spg2set
          !USE m_dwigner
          !USE m_strgn
          !USE m_mapatom
          !USE m_convn
          USE m_inpeig
          !USE m_ylm
          !-odim
          !USE m_od_mapatom
          !USE m_od_chisym
          !USE m_od_strgn1
          !+odim
          IMPLICIT NONE
          !     ..
          !     .. Scalars Arguments ..
          TYPE(t_atoms),INTENT(INOUT)    :: atoms
          TYPE(t_kpts),INTENT(INOUT)     :: kpts
          TYPE(t_sym),INTENT(INOUT)      :: sym
          TYPE(t_oneD),INTENT(INOUT)     :: oneD
          TYPE(t_input),INTENT(INOUT)    :: input
          TYPE(t_cell),INTENT(INOUT)     :: cell
          TYPE(t_enpara),INTENT(INOUT)   :: enpara
          CHARACTER(len=*),intent(in)    :: latnam,namgrp
          !     ..
          !     .. Local Scalars ..
          REAL       :: rkmaxx
          INTEGER    :: ntp1,ii,i,j,n1,n2,na,np1,n
          INTEGER, ALLOCATABLE :: lmx1(:), nq1(:), nlhtp1(:)
          !
          IF (namgrp.EQ.'any ') THEN
             CALL rw_symfile('R',94,'sym.out',sym%nop,cell%bmat, sym%mrot,sym%tau,sym%nop,sym%nop2,sym%symor)
          ELSE
             CALL spg2set(sym%nop,sym%zrfs,sym%invs,namgrp,latnam, sym%mrot,sym%tau,sym%nop2,sym%symor)
          ENDIF
          IF (input%film.AND..NOT.sym%symor) CALL juDFT_warn("Films&Symor",hint&
               &     ="Films should be symmorphic",calledby ='setup')


          !
          CALL inpeig(atoms,cell,input,oneD%odd%d1,kpts,enpara,latnam=latnam)
          !
          !
      
        END SUBROUTINE setup
      END MODULE m_setup
