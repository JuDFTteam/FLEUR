!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_vecforlo
  USE m_juDFT
  !----------------------------------------------------------------------------
  ! For a given atom (na) set up 2*llo(lo)+1 k-vectors for each lo on this atom.
  ! if it is the first of two inversion-related atoms, set up twice as much.
  !
  ! nkvec(nlod,1) = number of k-vectors that were set up
  ! kvec(2*(2*llod+1),nlod) = index of these k-vectors. Stored as kveclo on the
  !                           eig-file, for later use in charge-density part.
  !----------------------------------------------------------------------------
CONTAINS
  SUBROUTINE vec_for_lo(atoms,nintsp,sym,na,&
       n,np,noco, lapw,cell, gk,vk, nkvec,kvec)
    USE m_constants,ONLY: tpi_const,fpi_const
    USE m_orthoglo
    USE m_ylm

    USE m_types
    IMPLICIT NONE
    TYPE(t_noco),INTENT(IN)   :: noco
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_lapw),INTENT(IN)   :: lapw
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: nintsp ,na,n,np 
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: gk(:,:,:)!(dimension%nvd,3,nintsp)
    REAL,    INTENT (IN) :: vk(:,:,:)!(dimension%nvd,3,nintsp)
    INTEGER, INTENT (OUT):: kvec(2*(2*atoms%llod+1),atoms%nlod),nkvec(atoms%nlod,nintsp)
    !     ..
    !     .. Local Scalars ..
    COMPLEX term1 
    REAL th,con1
    INTEGER l,lo ,mind,ll1,lm,iintsp,k,nkmin,ntyp,lmp,m
    LOGICAL linind,enough,l_lo1
    !     ..
    !     .. Local Arrays ..
    REAL qssbti(3),bmrot(3,3),v(3),vmult(3)
    REAL :: gkrot(SIZE(gk,1),3,nintsp)
    REAL :: rph(SIZE(gk,1),nintsp)
    REAL :: cph(SIZE(gk,1),nintsp)
    COMPLEX ylm( (atoms%lmaxd+1)**2 )
    COMPLEX cwork(-2*atoms%llod:2*atoms%llod+1,2*(2*atoms%llod+1),atoms%nlod ,nintsp)
    !     ..
    !     .. Data statements ..
    REAL, PARAMETER :: eps = 1.0E-30
    REAL, PARAMETER :: linindq = 1.0e-4

    con1=fpi_const/SQRT(cell%omtil)
    ntyp = n
    DO iintsp = 1,nintsp
       IF (iintsp.EQ.1) THEN
          qssbti = - noco%qss/2
       ELSE
          qssbti = + noco%qss/2
       ENDIF

       !--->    set up phase factors
       DO k = 1,lapw%nv(iintsp)
          th= tpi_const*DOT_PRODUCT((/lapw%k1(k,iintsp),lapw%k2(k,iintsp),lapw%k3(k,iintsp)/)+qssbti,atoms%taual(:,na))
          rph(k,iintsp) = COS(th)
          cph(k,iintsp) = -SIN(th)
       END DO

       IF (np.EQ.1) THEN
          gkrot(:,:,iintsp)=gk(:,:,iintsp)
       ELSE
          bmrot=MATMUL(1.*sym%mrot(:,:,np),cell%bmat)
          DO k = 1,lapw%nv(iintsp)
             !-->           apply the rotation that brings this atom into the
             !-->           representative (this is the definition of ngopr(na))
             !-->           and transform to cartesian coordinates
             v(:) = vk(k,:,iintsp)
             gkrot(k,:,iintsp) = MATMUL(v,bmrot)
          END DO
       END IF
       !--->   end loop over interstitial spin
    ENDDO

    nkvec(:,:) = 0
    cwork(:,:,:,:) = CMPLX(0.0,0.0)
    enough=.FALSE.
    DO k = 1,MIN(lapw%nv(1),lapw%nv(nintsp))
       IF (ANY(lapw%rk(k,:nintsp)).LT.eps) CYCLE
       IF (.NOT.enough) THEN
          DO iintsp = 1,nintsp

             !-->        generate spherical harmonics
             vmult(:) =  gkrot(k,:,iintsp)
             CALL ylm4(atoms%lnonsph(ntyp),vmult, ylm)
                enough = .TRUE.
                term1 = con1* ((atoms%rmt(ntyp)**2)/2)* CMPLX(rph(k,iintsp),cph(k,iintsp))
                DO lo = 1,atoms%nlo(ntyp)
                   IF (atoms%invsat(na).EQ.0) THEN
                      IF ((nkvec(lo,iintsp)).LT. (2*atoms%llo(lo,ntyp)+1)) THEN
                         enough = .FALSE.
                         nkvec(lo,iintsp) = nkvec(lo,iintsp) + 1
                         l = atoms%llo(lo,ntyp)
                         ll1 = l*(l+1) + 1
                         DO m = -l,l
                            lm = ll1 + m
                            cwork(m,nkvec(lo,iintsp),lo,iintsp) = term1*ylm(lm)
                         END DO
                         CALL orthoglo(&
                              sym%invs,atoms,nkvec(lo,iintsp),lo,l,linindq,.FALSE., cwork(-2*atoms%llod,1,1,iintsp),linind)
                         IF (linind) THEN
                            kvec(nkvec(lo,iintsp),lo) = k
                         ELSE
                            nkvec(lo,iintsp) = nkvec(lo,iintsp) - 1
                         ENDIF
                      ENDIF
                   ELSE
                      IF ((atoms%invsat(na).EQ.1) .OR. (atoms%invsat(na).EQ.2)) THEN
                         IF (nkvec(lo,iintsp).LT.2*(2*atoms%llo(lo,ntyp)+1)) THEN
                            enough = .FALSE.
                            nkvec(lo,iintsp) = nkvec(lo,iintsp) + 1
                            l = atoms%llo(lo,ntyp)
                            ll1 = l*(l+1) + 1
                            DO m = -l,l
                               lm = ll1 + m
                               mind = -l + m
                               cwork(mind,nkvec(lo,iintsp),lo,iintsp) = term1*ylm(lm)
                               mind = l + 1 + m
                               lmp = ll1 - m
                               cwork(mind,nkvec(lo,iintsp),lo,iintsp) = ((-1)** (l+m))*CONJG(term1*ylm(lmp))
                            END DO
                            CALL orthoglo(&
                                 sym%invs,atoms,nkvec(lo,iintsp),lo,l,linindq,.TRUE., cwork(-2*atoms%llod,1,1,iintsp),linind)
                            IF (linind) THEN
                               kvec(nkvec(lo,iintsp),lo) = k
                               !                          write(*,*) nkvec(lo,iintsp),k,' <- '
                            ELSE
                               nkvec(lo,iintsp) = nkvec(lo,iintsp) - 1
                            END IF
                         END IF
                      END IF
                   END IF
                END DO
                IF ((k.EQ.lapw%nv(iintsp)) .AND. (.NOT.enough)) THEN
                   WRITE (6,FMT=*) 'abccoflo did not find enough linearly independent'
                   WRITE (6,FMT=*) 'clo coefficient-vectors. the linear independence'
                   WRITE (6,FMT=*) 'quality, linindq, is set: ',linindq
                   WRITE (6,FMT=*) 'this value might be to large.'
                   WRITE(*,*) na,k,lapw%nv 
                   CALL juDFT_error("not enough lin. indep. clo-vectors" ,calledby ="vec_for_lo")
                END IF
             ! -- >        end of abccoflo-part           
          ENDDO
       ENDIF

       ! -->    check whether we have already enough k-vecs
       enough=.TRUE.
       DO lo = 1,atoms%nlo(ntyp)
          IF (nkvec(lo,1).EQ.nkvec(lo,nintsp)) THEN   ! k-vec accepted by both spin channels
             IF (atoms%invsat(na).EQ.0) THEN
                IF ( nkvec(lo,1).LT.(2*atoms%llo(lo,ntyp)+1) ) THEN 
                   enough=.FALSE.
                ENDIF
             ELSE
                IF ( nkvec(lo,1).LT.(2*(2*atoms%llo(lo,ntyp)+1))  ) THEN
                   enough=.FALSE.
                ENDIF
             ENDIF
          ELSE
             nkmin = MIN(nkvec(lo,1),nkvec(lo,nintsp)) ! try another k-vec
             nkvec(lo,1) = nkmin ; nkvec(lo,nintsp) = nkmin
             enough=.FALSE.
          ENDIF
       ENDDO
       IF ( enough ) THEN
          !           DO  lo = 1,nlo(ntyp)
          !             DO l = 1, nkvec(lo,1)
          !              write(*,*) lo,l,kvec(l,lo)
          !             ENDDO
          !           ENDDO
          RETURN
       ENDIF
    ENDDO

  END SUBROUTINE vec_for_lo
      END MODULE m_vecforlo
