      MODULE m_checkdop
      CONTAINS
        SUBROUTINE checkdop(&
             &                    p,np,n,na,ivac,iflag,jsp,&
             &                    DIMENSION,atoms,sphhar,stars,sym,&
             &                    vacuum,cell,oneD,potden)
          ! ************************************************************
          !     subroutines checks the continuity of coulomb           *
          !     potential or valence charge density                    *
          !     across the m.t. and vacuum boundaries                  *
          !                                     c.l.fu                 *
          !     (unifies vcheck and cdnchk      g.b.                   *
          ! YM:  this routine doesn't really work in the vacuum in 1D case yet
          ! ************************************************************

          USE m_juDFT
          USE m_starf, ONLY : starf2,starf3
          USE m_angle
          USE m_ylm
          USE m_types
          USE m_constants
          USE m_fitchk
          IMPLICIT NONE
          !     ..
          !     .. Scalar Arguments ..
          TYPE(t_dimension),INTENT(IN) :: dimension
          type(t_sphhar),intent(in)    :: sphhar      
          TYPE(t_stars),INTENT(IN)     :: stars
          TYPE(t_atoms),INTENT(IN)     :: atoms
          TYPE(t_sym),INTENT(IN)       :: sym
          TYPE(t_vacuum),INTENT(IN)    :: vacuum
          TYPE(t_oneD),INTENT(IN)      :: oneD
          TYPE(t_cell),INTENT(IN)      :: cell
          TYPE(t_potden),INTENT(IN)    :: potden

          INTEGER, INTENT (IN) :: iflag,ivac,n,na,np,jsp
          !-odim
          !+odim
          !     .. Array Arguments ..
          REAL,    INTENT (IN) :: p(:,:)!(3,DIMENSION%nspd)
          !     ..
          !     .. Local Scalars ..
          REAL av,dms,rms,s,ir2,help,phi
          INTEGER i,j,k,lh,mem,nd,lm,ll1,nopa ,gz,m
          COMPLEX ic
          LOGICAL l_cdn
          !     ..
          !     .. Local Arrays ..
          COMPLEX sf2(stars%ng2),sf3(stars%ng3),ylm( (atoms%lmaxd+1)**2 )
          REAL rcc(3),v1(SIZE(p,2)),v2(SIZE(p,2)),x(3),ri(3)

          l_cdn = .FALSE. ! By default we assume that the input is a potential.
          IF (potden%potdenType.LE.0) CALL juDFT_error('unknown potden type', calledby='checkdop')
          IF (potden%potdenType.GT.1000) l_cdn = .TRUE. ! potdenTypes > 1000 are reserved for densities
          

          !     ..
          !     ..
#ifdef  __TOS_BGQ__
          RETURN
#endif
          ic = CMPLX(0.,1.)
          IF (.NOT.iflag.LT.0) THEN
             !     ---> Interstitial part
             DO j = 1,np
                IF (.NOT.oneD%odi%d1) THEN
                   !CALL cotra1(p(:,j),rcc,cell%bmat)
                   rcc=MATMUL(cell%bmat,p(:,j))/tpi_const
                   CALL starf3(&
                        &               sym%nop,stars%ng3,sym%symor,stars%kv3,sym%mrot,sym%tau,p(:,j),sym%invtab,&
                        &               sf3)!keep
                ENDIF
                !
                IF (oneD%odi%d1) THEN
                   !CALL cotra1(p(:,j),rcc,cell%bmat)
                   rcc=MATMUL(cell%bmat,p(:,j))/tpi_const
                   CALL starf3(&
                        &               sym%nop,stars%ng3,sym%symor,stars%kv3,sym%mrot,sym%tau,rcc,sym%invtab,&
                        &               sf3)!keep
                ENDIF
                v1(j) = 0.0
                DO k = 1,stars%ng3
                   v1(j) = v1(j) + REAL(potden%pw(k,jsp)*sf3(k))*stars%nstr(k)
                ENDDO
             ENDDO
             !     ---> vacuum part
             IF (l_cdn) THEN
                WRITE (6,FMT=9000) ivac
             ELSE
                WRITE (6,FMT=8000) ivac
             ENDIF
             DO  j = 1,np
                IF (.NOT.oneD%odi%d1) THEN
                   CALL starf2(&
                        &           sym%nop2,stars%ng2,stars%kv2,sym%mrot,sym%symor,sym%tau,p(1,j),sym%invtab,&
                        &           sf2)!keep
                   v2(j) = potden%vacz(1,ivac,jsp)
                   DO  k = 2,stars%ng2
                      v2(j) = v2(j) + REAL(potden%vacxy(1,k-1,ivac,jsp)*sf2(k))*stars%nstr2(k)
                   ENDDO
                ELSE
                   !-odim
                   v2(j) = potden%vacz(1,ivac,jsp)
                   phi = angle(p(1,j),p(2,j))
                   DO  k = 2,oneD%odi%nq2
                      m = oneD%odi%kv(2,k)
                      gz = oneD%odi%kv(1,k)
                      v2(j) = v2(j) + REAL(potden%vacxy(1,k-1,ivac,jsp)*&
                           &           EXP(ic*m*phi)*EXP(ic*cell%bmat(3,3)*gz*p(3,j)))*oneD%odi%nst2(k)
                   ENDDO
                   !+odim
                END IF
                IF (oneD%odi%d1) THEN
                   !CALL cotra1(p(:,j),rcc,cell%bmat)
                   rcc=MATMUL(cell%bmat,p(:,j))/tpi_const

                   WRITE (6,FMT=8020)  rcc,(p(i,j),i=1,3),v1(j),v2(j)
                ELSE
                   !CALL cotra0(p(1,j),rcc,cell%amat)
                   rcc=MATMUL(cell%amat,p(:,j))
                   WRITE (6,FMT=8020) (p(i,j),i=1,3),rcc,v1(j),v2(j)
                ENDIF
             ENDDO
             CALL fitchk(v1(:np),v2(:np),av,rms,dms)
             WRITE (6,FMT=8030) av,rms,dms
             RETURN
          ENDIF
          !      ----> interstitial part
          DO j = 1,np
             !CALL cotra1(p(1,j),rcc,cell%bmat)
             rcc=MATMUL(cell%bmat,p(:,j))/tpi_const

             CALL starf3(&
                  &               sym%nop,stars%ng3,sym%symor,stars%kv3,sym%mrot,sym%tau,rcc,sym%invtab,&
                  &               sf3)!keep
             !
             v1(j) = 0.0
             DO k = 1,stars%ng3
                v1(j) = v1(j) + REAL(potden%pw(k,jsp)*sf3(k))*stars%nstr(k)
             ENDDO
          ENDDO
          !     ----> m.t. part
          IF (l_cdn) THEN
             WRITE (6,FMT=9010) n
          ELSE
             WRITE (6,FMT=8010) n
          ENDIF
          ir2 = 1.0
          IF (l_cdn) ir2 = 1.0 / ( atoms%rmt(n)*atoms%rmt(n) )
          nd = sym%ntypsy(na)
          nopa = sym%ngopr(na)
          IF (oneD%odi%d1) THEN
             nopa = oneD%ods%ngopr(na)
             nopa = oneD%ods%invtab(nopa)
          END IF
          DO j = 1,np
             DO i = 1,3
                x(i) = p(i,j) - atoms%pos(i,na)
             ENDDO
             ! new
             IF (nopa.NE.1) THEN
                !CALL cotra1(x,rcc,cell%bmat)  ! switch to internal units
                rcc=MATMUL(cell%bmat,x)/tpi_const

                DO i = 1,3               ! rotate into representative
                   ri(i) = 0.
                   DO k = 1,3
                      IF (oneD%odi%d1) THEN
                         ri(i) = ri(i) + oneD%ods%mrot(i,k,nopa)*rcc(k)
                      ELSE
                         ri(i) = ri(i) + sym%mrot(i,k,nopa)*rcc(k)
                      END IF
                   ENDDO
                ENDDO
                !CALL cotra0(ri,x,cell%amat)    !switch back to cartesian units
                x=MATMUL(cell%amat,ri)

             END IF
             ! new
             CALL ylm4(&
                  &             atoms%lmax(n),x,&
                  &             ylm)
             help = 0.0
             DO lh = 0,sphhar%nlh(nd)
                s = 0.0
                ll1 = sphhar%llh(lh,nd) * ( sphhar%llh(lh,nd) + 1 ) + 1
                DO mem = 1,sphhar%nmem(lh,nd)
                   lm = ll1 + sphhar%mlh(mem,lh,nd)
                   s = s + REAL( sphhar%clnu(mem,lh,nd)* ylm(lm) )
                ENDDO
                help = help + potden%mt(atoms%jri(n),lh,n,jsp) * s
             ENDDO
             v2(j) = help * ir2 
             IF (j.LE.8) THEN
                !CALL cotra1(p(1,j),rcc,cell%bmat)
                rcc=MATMUL(cell%bmat,p(:,j))/tpi_const

                WRITE (6,FMT=8020) rcc, (p(i,j),i=1,3),v1(j),v2(j)
             END IF
          ENDDO
          CALL fitchk(v1(:np),v2(:np),av,rms,dms)
          WRITE (6,FMT=8030) av,rms,dms
8000      FORMAT (/,'    int.-vac. boundary (potential): ivac=',i2,/,t10,&
               &       'int-coord',t36,'cart-coord',t57,' inter. ',t69,' vacuum ')
8010      FORMAT (/,'    int.-m.t. boundary (potential): atom type=',i2,/,&
               &       t10,'int-coord',t36,'cart-coord',t57,' inter. ',t69,&
               &       '  m. t. ')
8020      FORMAT (1x,2 (3f8.3,2x),2f12.6)
8030      FORMAT (/,10x,'average value = ',f10.6,/,t11,'rms,dmx=',2f7.3,&
               &       ' per cent')
9000      FORMAT (/,'    int.-vac. boundary (density): ivac=',i2,/,t10,&
               &       'int-coord',t36,'cart-coord',t57,' inter. ',t69,' vacuum ')
9010      FORMAT (/,'    int.-m.t. boundary (density): atom type=',i2,/,t10,&
               &       'int-coord',t36,'cart-coord',t57,' inter. ',t69,'  m. t. ')
          RETURN
        END SUBROUTINE checkdop
      END MODULE m_checkdop
