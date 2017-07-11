          IF (    (xcpot%icorr.EQ.icorr_hf ) .OR. (xcpot%icorr.EQ.icorr_pbe0)&
               &    .OR.(xcpot%icorr.EQ.icorr_exx) .OR. (xcpot%icorr.EQ.icorr_hse)&
               &    .OR.(xcpot%icorr.EQ.icorr_vhse) ) THEN
             IF (input%film .OR. oneD%odi%d1)&
                  &    CALL juDFT_error("2D film and 1D calculations not implemented"&
                  &                 //"for HF/EXX/PBE0/HSE", calledby ="fleur",&
                  &                 hint="Use a supercell or a different functional")

             IF( ANY( atoms%l_geo  ) )&
                  &     CALL juDFT_error("Forces not implemented for HF/PBE0/HSE ",&
                  &                    calledby ="fleur")

             IF (.NOT. obsolete%pot8) STOP 'Choose pot8=T'
             !calculate whole Brilloun zone
             CALL gen_bz(kpts,sym)
             CALL gen_map(&
                  &          atoms,sym,oneD,hybrid)
             !
             ! calculate d_wgn
             !
             ALLOCATE (hybrid%d_wgn2(-atoms%lmaxd:atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,0:atoms%lmaxd,sym%nsym))
             CALL d_wigner(sym%nop,sym%mrot,cell%bmat,atoms%lmaxd,hybrid%d_wgn2(:,:,1:,:sym%nop))
             hybrid%d_wgn2(:,:,0,:) = 1

             DO isym = sym%nop+1,sym%nsym
                iisym = isym - sym%nop
                DO l = 0,atoms%lmaxd
                   DO m2 = -l,l
                      DO m1 = -l,-1
                         cdum                  = hybrid%d_wgn2( m1,m2,l,iisym)
                         hybrid%d_wgn2( m1,m2,l,isym) = hybrid%d_wgn2(-m1,m2,l,iisym)*(-1)**m1
                         hybrid%d_wgn2(-m1,m2,l,isym) = cdum                  *(-1)**m1
                      END DO
                      hybrid%d_wgn2(0,m2,l,isym) = hybrid%d_wgn2(0,m2,l,iisym)
                   END DO
                END DO
             END DO
