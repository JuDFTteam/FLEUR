      MODULE m_loddop
      USE m_juDFT
      CONTAINS
        SUBROUTINE loddop(&
             &                  stars,vacuum,atoms,sphhar,&
             &                  input,sym,&
             &                  nu,&
             &                  it,fr,fpw,fz,fzxy)
          !     ***********************************************************
          !     reload formatted density or potential   c.l.fu
          !     ***********************************************************
          USE m_types
          IMPLICIT NONE
          !
          !     .. Scalar Arguments ..
          TYPE(t_stars),INTENT(IN)  :: stars
          TYPE(t_vacuum),INTENT(IN) :: vacuum
          TYPE(t_atoms),INTENT(IN)  :: atoms
          TYPE(t_sphhar),INTENT(IN) :: sphhar
          TYPE(t_input),INTENT(IN)  :: input
          TYPE(t_sym),INTENT(IN)    :: sym

          INTEGER, INTENT (IN) :: nu  
          INTEGER, INTENT (OUT):: it
          !     ..
          !     .. Array Arguments ..
          COMPLEX, INTENT (OUT):: fpw(:,:),fzxy(:,:,:,:)!(stars%ng3,input%jspins),fzxy(vacuum%nmzxyd,stars%ng2-1,2,input%jspins)
          REAL,    INTENT (OUT):: fr(:,0:,:,:)!(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins)
          REAL,    INTENT (OUT):: fz(:,:,:)!(vacuum%nmzd,2,input%jspins)
          CHARACTER(len=8) :: dop,iop,name(10)
          !     ..
          !     .. Local Scalars ..
          REAL delzn,dxn,rmtn,z1n,dummy
          INTEGER i,ivac,ivdummy,j,jrin,jsp,jspdum,k,lh,lhdummy,n,ndum,nlhn,&
               &        nmzn,nmzxyn,nn,nq2n,nq3n,ntydum,n_diff,na
          CHARACTER(len=2) namaux
          !     ..
          !     .. Local Arrays ..
          REAL, ALLOCATABLE :: fpwr(:,:),fzxyr(:,:,:,:)

          CHARACTER(len=8) space(10)
          !     ..
          !     .. Intrinsic Functions ..
          INTRINSIC cmplx
          !     ..
          !     .. Data statements ..
          DATA space/10*'        '/
          !     ..

          fr = 0 ; fzxy = 0 ; fr = 0 ; fz = 0

          IF (sym%invs) ALLOCATE ( fpwr(stars%ng3,SIZE(fpw,2)) )
          IF (sym%invs2) ALLOCATE ( fzxyr(vacuum%nmzxyd,stars%ng2-1,2,SIZE(fzxy,4)) )

          name = space
          READ (nu,END=200,ERR=200) name
          !      WRITE (*,FMT=8000) name
          ! 8000 FORMAT (' loddop title:',10a8)
          READ (nu,END=200,ERR=200) iop,dop,it
          DO  jsp = 1,SIZE(fr,4)
             READ (nu,END=200,ERR=200) jspdum
             READ (nu,END=200,ERR=200) nn
             IF (nn/=atoms%ntype) CALL juDFT_error("nn.NE.ntype",calledby =&
                  &       "loddop" )

             na = 1
             DO  n = 1,nn
                READ (nu,END=200,ERR=200) namaux,ndum,jrin,rmtn,dxn
                READ (nu,END=200,ERR=200) ntydum,nlhn
                !+gu
                IF ( nlhn.GT.sphhar%nlh(sym%ntypsy(na)) ) THEN
                   WRITE (*,*) 'nlh (',nlhn,') set to (',sphhar%nlh(sym%ntypsy(na)),')'
                   n_diff = nlhn - sphhar%nlh(sym%ntypsy(na))
                   nlhn = sphhar%nlh(sym%ntypsy(na))
                ELSE
                   n_diff = 0 
                ENDIF
                !-gu
                DO  lh = 0,nlhn
                   READ (nu,END=200,ERR=200) lhdummy
                   READ (nu,END=200,ERR=200) (fr(i,lh,n,jsp),i=1,jrin)
                ENDDO
                IF (nlhn.LT.sphhar%nlh(sym%ntypsy(na))) THEN
                   DO lh = nlhn + 1,sphhar%nlh(sym%ntypsy(na))
                      DO i = 1,atoms%jri(n)
                         fr(i,lh,n,jsp) = 0.
                      ENDDO
                   ENDDO
                ELSE
                   DO lh = 1, n_diff
                      READ (nu,END=200,ERR=200) lhdummy
                      READ (nu,END=200,ERR=200) dummy
                   ENDDO
                ENDIF

                na = na + atoms%neq(n)
             ENDDO
             IF (jsp<=SIZE(fpw,2)) THEN
                READ (nu,END=200,ERR=200) nq3n
                !+gu
                IF (nq3n.GT.stars%ng3) THEN
                   WRITE (*,*) 'nq3n (',nq3n,') reduced to nq3 (',stars%ng3,')'
                   nq3n = stars%ng3
                ENDIF
                !-gu
                IF (sym%invs) THEN
                   READ (nu,END=200,ERR=200) (fpwr(k,jsp),k=1,nq3n)
                   fpw(:nq3n,jsp) = CMPLX(fpwr(:nq3n,jsp),0.)
                   
                ELSE
                   READ (nu,END=200,ERR=200) (fpw(k,jsp),k=1,nq3n)
                END IF
                IF (nq3n.LT.stars%ng3) THEN
                   fpw(nq3n+1:,jsp) = (0.,0.)
                END IF
             ENDIF
             IF (input%film) THEN
                IF (jsp<=SIZE(fz,3)) THEN
                   DO  ivac = 1,vacuum%nvac
                      READ (nu,END=200,ERR=200) ivdummy
                      READ (nu,END=200,ERR=200) nmzn,z1n,delzn
                      READ (nu,END=200,ERR=200) (fz(i,ivac,jsp),i=1,nmzn)
                      IF (vacuum%nvac.EQ.1) THEN
                         DO i=1,nmzn
                            fz(i,2,jsp)=fz(i,1,jsp)
                         ENDDO
                      ENDIF
                      IF (jsp<=SIZE(fzxy,4)) THEN
                         READ (nu,END=200,ERR=200) nq2n,nmzxyn
                         !+gu
                         IF (nq2n.GT.stars%ng2) THEN
                            WRITE (*,*) 'nq2n (',nq2n,') reduced to nq2 (',stars%ng2,')'
                            n_diff = nq2n - stars%ng2
                            nq2n = stars%ng2
                         ELSE
                            n_diff = 0
                         ENDIF
                         !-gu
                         DO  k = 2,nq2n
                            IF (sym%invs2) THEN
                               READ (nu,END=200,ERR=200) &
                                    &                              (fzxyr(j,k-1,ivac,jsp),j=1,nmzxyn)
                               fzxy(:nmzxyn,k-1,ivac,jsp) = CMPLX(fzxyr(:nmzxyn,k-1,ivac,&
                                    &                                         jsp),0.)
                            ELSE
                               READ (nu,END=200,ERR=200)  &
                                    &                               (fzxy(j,k-1,ivac,jsp),j=1,nmzxyn)
                            END IF
                            IF (vacuum%nvac.EQ.1) THEN
                               IF (sym%invs) THEN
                                  DO j = 1,nmzxyn
                                     fzxy(j,k-1,2,jsp) = CONJG(fzxy(j,k-1,1,jsp))
                                  ENDDO
                               ELSE
                                  DO j = 1,nmzxyn
                                     fzxy(j,k-1,2,jsp) = fzxy(j,k-1,1,jsp)
                                  ENDDO
                               ENDIF
                            ENDIF
                         ENDDO
                         !+gu
                         DO k = 1,n_diff
                            READ (nu,END=200,ERR=200) dummy
                         ENDDO
                         !-gu
                         IF (nq2n.LT.stars%ng2) THEN
                            fzxy(:nmzxyn,nq2n:,ivac,jsp) = (0.,0.)
                         END IF
                      ENDIF
                   ENDDO
                END IF
             ENDIF
          ENDDO
          !
          IF (sym%invs) DEALLOCATE (fpwr)
          IF (sym%invs2) DEALLOCATE ( fzxyr )
          RETURN

200       WRITE (6,*) 'error reading dop nr.',nu
          IF (nu /= 98)  CALL juDFT_error("error reading d/p-file!",calledby&
          &     ="loddop")

        END SUBROUTINE loddop
      END MODULE m_loddop
