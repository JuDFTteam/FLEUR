      MODULE m_setup
      USE m_juDFT
      CONTAINS
        SUBROUTINE setup(mpi,atoms,kpts,DIMENSION,sphhar,&
                  obsolete,sym,stars,oneD, input,noco,vacuum,cell,xcpot, sliceplot,enpara,l_opti)
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
          LOGICAL, INTENT (IN) :: l_opti  
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
          CALL ylmnorm_init(atoms%lmaxd)
          IF (.NOT.oneD%odd%d1) THEN
             CALL local_sym(&
                  atoms%lmaxd,atoms%lmax,sym%nop,sym%mrot,sym%tau,&
                  atoms%nat,atoms%ntype,atoms%neq,cell%amat,cell%bmat,atoms%taual,&
                  sphhar%nlhd,sphhar%memd,sphhar%ntypsd,.FALSE.,&
                  atoms%nlhtyp,atoms%ntypsy,sphhar%nlh,sphhar%llh,sphhar%nmem,sphhar%mlh,sphhar%clnu)
             sym%nsymt = sphhar%ntypsd

             oneD%mrot1(:,:,:) = sym%mrot(:,:,:)
             oneD%tau1(:,:) = sym%tau(:,:)
          ELSEIF (oneD%odd%d1) THEN
             CALL od_chisym(oneD%odd,oneD%mrot1,oneD%tau1,sym%zrfs,sym%invs,sym%invs2,cell%amat)
             ntp1 = atoms%nat
             ALLOCATE (nq1(ntp1),lmx1(ntp1),nlhtp1(ntp1))
             ii = 1
             DO i = 1,atoms%ntype
                DO j = 1,atoms%neq(i)
                   nq1(ii) = 1
                   lmx1(ii) = atoms%lmax(i)
                   ii = ii + 1
                END DO
             END DO
             CALL local_sym(&
                  atoms%lmaxd,lmx1,sym%nop,sym%mrot,sym%tau,&
                  atoms%nat,ntp1,nq1,cell%amat,cell%bmat,atoms%taual,&
                  sphhar%nlhd,sphhar%memd,sphhar%ntypsd,.FALSE.,&
                  nlhtp1,atoms%ntypsy,sphhar%nlh,sphhar%llh,sphhar%nmem,sphhar%mlh,sphhar%clnu)

             sym%nsymt = sphhar%ntypsd
             ii = 1
             DO i = 1,atoms%ntype
                atoms%nlhtyp(i) = nlhtp1(ii)
                ii = ii + atoms%neq(i)
             END DO
             DEALLOCATE (lmx1,nlhtp1)
          END IF
          !+odim
          IF (atoms%n_u.GT.0) THEN
             CALL d_wigner(sym%nop,sym%mrot,cell%bmat,3, sym%d_wgn)
          ENDIF
          !
          !+odim
          IF (.NOT.oneD%odd%d1) THEN
             CALL mapatom(sym,atoms, cell,input, noco)
             oneD%ngopr1(1:atoms%nat) = atoms%ngopr(1:atoms%nat)
             !        DEALLOCATE ( nq1 )
          ELSE
             !-odim
             CALL juDFT_error("The oneD version is broken here. Compare call to mapatom with old version")
             CALL mapatom(sym,atoms, cell,input, noco)
             CALL od_mapatom(oneD,atoms,sym,cell)
          END IF

          ! Store structure data

          CALL storeStructureIfNew(input,stars, atoms, cell, vacuum, oneD, sym,mpi,sphhar,noco)

          !+odim
          IF (input%film.OR.(sym%namgrp.NE.'any ')) THEN
             CALL strgn1(stars,sym,atoms, vacuum,sphhar, input,cell,xcpot)
             !-odim
             IF (oneD%odd%d1) THEN
                CALL od_strgn1(xcpot,cell,sym,oneD)
             END IF
             !+odim
          ELSE
             CALL strgn2(stars,sym,atoms, vacuum,sphhar, input,cell,xcpot)
          ENDIF
          !
          !IF (.NOT.l_opti) THEN
             CALL inpeig(atoms,cell,input,oneD%odd%d1,kpts,enpara)
          !ENDIF
          !
          !-----> prepare dimensions for charge density fft-box in pwden.f
          !
          CALL  prp_qfft(stars, cell,noco, input)

          !
          !-----> prepare dimensions for xc fft-box in visxc(g).f
          !
          
         
          CALL  prp_xcfft(stars,input, cell, xcpot)
          !
        ENDIF ! (mpi%irank == 0)
          CALL stepf(sym,stars,atoms,oneD, input,cell, vacuum,mpi)
          IF (.NOT.sliceplot%iplot) THEN
             IF ( mpi%irank == 0 ) THEN
                CALL convn(DIMENSION,atoms,stars)

                !--->    set up electric field parameters (if needed) 
                ! CALL e_field(atoms, DIMENSION, stars, sym, vacuum, cell, input,field)
             ENDIF
          ENDIF

        END SUBROUTINE setup
      END MODULE m_setup
