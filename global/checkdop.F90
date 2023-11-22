      MODULE m_checkdop
      CONTAINS
        SUBROUTINE checkdop(p,np,n,na,ivac,iflag,jsp,atoms,sphhar,stars,sym,&
                            vacuum,cell,potden,potdenIm)
          ! ************************************************************
          !     subroutines checks the continuity of coulomb           *
          !     potential or valence charge density                    *
          !     across the m.t. and vacuum boundaries                  *
          !                                     c.l.fu                 *
          !     (unifies vcheck and cdnchk      g.b.                   *
          ! YM:  this routine doesn't really work in the vacuum in 1D case yet
          ! ************************************************************

          USE m_types
          USE m_constants
          USE m_juDFT
          USE m_starf, ONLY : starf2,starf3
          USE m_angle
          USE m_ylm
          USE m_fitchk

          IMPLICIT NONE

          !     ..
          !     .. Scalar Arguments ..
          
          type(t_sphhar),intent(in)    :: sphhar
          TYPE(t_stars),INTENT(IN)     :: stars
          TYPE(t_atoms),INTENT(IN)     :: atoms
          TYPE(t_sym),INTENT(IN)       :: sym
          TYPE(t_vacuum),INTENT(IN)    :: vacuum
          TYPE(t_cell),INTENT(IN)      :: cell
          TYPE(t_potden),INTENT(IN)    :: potden
          TYPE(t_potden),INTENT(IN), OPTIONAL    :: potdenIm
          INTEGER, INTENT (IN) :: iflag,ivac,n,na,np,jsp
          !-odim
          !+odim
          !     .. Array Arguments ..
          REAL,    INTENT (IN) :: p(:,:)!(3,DIMENSION%nspd)
          !     ..
          !     .. Local Scalars ..
          REAL av,dms,rms,s,ir2,help,phi,helpIm
          INTEGER i,j,k,lh,mem,nd,lm,ll1,nopa ,gz,m
          COMPLEX ic
          LOGICAL l_cdn, l_dfpt
          !     ..
          !     .. Local Arrays ..
          COMPLEX sf2(stars%ng2),sf3(stars%ng3),ylm( (atoms%lmaxd+1)**2 )
          REAL rcc(3),v1(SIZE(p,2)),v2(SIZE(p,2)),x(3),ri(3), v1Im(SIZE(p,2)), v2Im(SIZE(p,2))
          l_dfpt = .FALSE.
          l_dfpt = PRESENT(potdenIm)
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
                rcc=MATMUL(cell%bmat,p(:,j))/tpi_const
                IF (l_dfpt) THEN 
                  CALL starf3(sym%nop,stars%ng3,sym%symor,stars%kv3,sym%mrot,sym%tau,p(:,j),sym%invtab,sf3,stars%center) ! Stars shifted by q
                ELSE
                  CALL starf3(sym%nop,stars%ng3,sym%symor,stars%kv3,sym%mrot,sym%tau,p(:,j),sym%invtab,sf3) 
                END IF
               !
                v1(j) = 0.0
                v1Im(j)=0.0
                If (l_dfpt) v1Im(j) = 0.0 
                DO k = 1,stars%ng3
                   v1(j) = v1(j) + REAL(potden%pw(k,jsp)*sf3(k))*stars%nstr(k)
                   IF (l_dfpt) v1Im(j) = v1Im(j) + AIMAG(potden%pw(k,jsp)*sf3(k))*stars%nstr(k)
                ENDDO
             ENDDO
             !     ---> vacuum part
             IF (l_cdn) THEN
               IF (l_dfpt) THEN
                  WRITE (oUnit,FMT=9005) ivac
                ELSE 
                  WRITE (oUnit,FMT=9000) ivac
                END IF
             ELSE
                IF (l_dfpt) THEN
                  WRITE (oUnit,FMT=8005) ivac
                ELSE 
                  WRITE (oUnit,FMT=8000) ivac
                END IF
             ENDIF
             DO  j = 1,np
                   IF (l_dfpt) THEN
                        CALL starf2(sym%nop2,stars%ng2,stars%kv2,sym%mrot,sym%symor,sym%tau,p(1:3,j),sym%invtab,sf2,stars%center)
                   ELSE
                        CALL starf2(sym%nop2,stars%ng2,stars%kv2,sym%mrot,sym%symor,sym%tau,p(1:3,j),sym%invtab,sf2)
                   END IF 
                   v2(j)=0.0
                   v2Im(j)=0.0
                   IF ((norm2(stars%center)<1e-8)) THEN 
                     v2(j) = REAL(potden%vac(1,1,ivac,jsp))
                     IF (l_dfpt) v2Im(j)=AIMAG(potden%vac(1,1,ivac,jsp)) !(l_dfpt Gamma Point)
                     DO  k = 2,stars%ng2
                       v2(j) = v2(j) + REAL(potden%vac(1,k,ivac,jsp)*sf2(k))*stars%nstr2(k)
                       IF (l_dfpt) v2Im(j) = v2Im(j) + AIMAG(potden%vac(1,k,ivac,jsp)*sf2(k))*stars%nstr2(k)
                     END DO 
                   ELSE ! (l_dfpt)
                     DO  k = 1,stars%ng2
                        v2(j) = v2(j) + REAL(potden%vac(1,k,ivac,jsp)*sf2(k))*stars%nstr2(k)
                        v2Im(j) = v2Im(j) + AIMAG(potden%vac(1,k,ivac,jsp)*sf2(k))*stars%nstr2(k)
                     ENDDO
                  END IF 
                   rcc=MATMUL(cell%amat,p(:,j))
                   IF (l_dfpt) THEN 
                     WRITE (oUnit,FMT=8025) (p(i,j),i=1,3),rcc, v1(j),v1Im(j),v2(j),v2Im(j)
                   ELSE
                     WRITE (oUnit,FMT=8020) (p(i,j),i=1,3),rcc, v1(j),v2(j)
                   END IF
             END DO
             CALL fitchk(v1(:np),v2(:np),av,rms,dms)
             WRITE (oUnit,FMT=8030) av,rms,dms
             IF (l_dfpt) THEN
               CALL fitchk(v1Im(:np),v2Im(:np),av,rms,dms)
               WRITE (oUnit,*) "Imaginary Part:"
               WRITE (oUnit,FMT=8030) av,rms,dms
             END IF 
             RETURN
          ENDIF
          !      ----> interstitial part
          DO j = 1,np
             rcc=MATMUL(cell%bmat,p(:,j))/tpi_const
             IF (l_dfpt) THEN
               CALL starf3(sym%nop,stars%ng3,sym%symor,stars%kv3,sym%mrot,sym%tau,rcc,sym%invtab,sf3,stars%center)
             ELSE
               CALL starf3(sym%nop,stars%ng3,sym%symor,stars%kv3,sym%mrot,sym%tau,rcc,sym%invtab,sf3)
             END IF 
             !
             v1(j) = 0.0
             IF (l_dfpt) v1Im(j) = 0.0 
             DO k = 1,stars%ng3
                v1(j) = v1(j) + REAL(potden%pw(k,jsp)*sf3(k))*stars%nstr(k)
                IF (l_dfpt) v1Im(j) = v1Im(j) + AIMAG(potden%pw(k,jsp)*sf3(k))*stars%nstr(k)
             ENDDO
          ENDDO
          !     ----> m.t. part
          IF (l_cdn) THEN
            IF (l_dfpt) THEN
               WRITE (oUnit,FMT=9015) n
             ELSE
               WRITE (oUnit,FMT=9010) n
             END IF 
          ELSE
             IF (l_dfpt) THEN
               WRITE (oUnit,FMT=8015) n
             ELSE
               WRITE (oUnit,FMT=8010) n
             END IF 
          ENDIF
          ir2 = 1.0
          IF (l_cdn) ir2 = 1.0 / ( atoms%rmt(n)*atoms%rmt(n) )
          nd = sym%ntypsy(na)
          nopa = sym%ngopr(na)
          DO j = 1,np
             DO i = 1,3
                x(i) = p(i,j) - atoms%pos(i,na)
             ENDDO
             ! new
             IF (nopa.NE.1) THEN
                rcc=MATMUL(cell%bmat,x)/tpi_const

                DO i = 1,3               ! rotate into representative
                   ri(i) = 0.
                   DO k = 1,3
                         ri(i) = ri(i) + sym%mrot(i,k,nopa)*rcc(k)
                   ENDDO
                ENDDO
                x=MATMUL(cell%amat,ri)

             END IF
             ! new
             CALL ylm4(atoms%lmax(n),x,ylm)
             help = 0.0
             helpIm= 0.0 
             DO lh = 0,sphhar%nlh(nd)
                s = 0.0
                ll1 = sphhar%llh(lh,nd) * ( sphhar%llh(lh,nd) + 1 ) + 1
                DO mem = 1,sphhar%nmem(lh,nd)
                   lm = ll1 + sphhar%mlh(mem,lh,nd)
                   s = s + REAL( sphhar%clnu(mem,lh,nd)* ylm(lm) )
                ENDDO
                help = help + potden%mt(atoms%jri(n),lh,n,jsp) * s
                IF (l_dfpt) helpIm=helpIm + potdenIm%mt(atoms%jri(n),lh,n,jsp) * s
             ENDDO
             v2(j) = help * ir2
             IF (l_dfpt) v2Im(j) = helpIm * ir2
             IF (j.LE.8) THEN
                rcc=MATMUL(cell%bmat,p(:,j))/tpi_const

                IF (l_dfpt) THEN 
                  WRITE (oUnit,FMT=8025) rcc, (p(i,j),i=1,3),v1(j),v1Im(j),v2(j),v2Im(j)
                ELSE
                  WRITE (oUnit,FMT=8020) rcc, (p(i,j),i=1,3),v1(j),v2(j)
                END IF
             END IF
         END DO 
          CALL fitchk(v1(:np),v2(:np),av,rms,dms)
          WRITE (oUnit,FMT=8030) av,rms,dms
          IF (l_dfpt) THEN
               CALL fitchk(v1Im(:np),v2Im(:np),av,rms,dms)
               WRITE (oUnit,*) "Imaginary Part:"
               WRITE (oUnit,FMT=8030) av,rms,dms
          END IF 
8000      FORMAT (/,'    int.-vac. boundary (potential): ivac=',i2,/,t10,&
               &       'int-coord',t36,'cart-coord',t57,' inter. ',t69,' vacuum ')
8005      FORMAT (/,'    int.-vac. boundary (potential): ivac=',i2,/,t10,&
               &       'int-coord',t36,'cart-coord',t64,' inter. ',t87,' vacuum ')
8010      FORMAT (/,'    int.-m.t. boundary (potential): atom type=',i2,/,&
               &       t10,'int-coord',t36,'cart-coord',t57,' inter. ',t69,&
               &       '  m. t. ')
8015      FORMAT (/,'    int.-m.t. boundary (potential): atom type=',i2,/,&
               &       t10,'int-coord',t36,'cart-coord',t64,' inter. ',t87,&
               &       '  m. t. ')
8020      FORMAT (1x,2 (3f8.3,2x),2f12.6)
8025      FORMAT (1x,2 (3f8.3,2x),4f12.6) !;_dfpt
8030      FORMAT (/,10x,'average value = ',f12.6,/,t11,'rms,dmx=',2f10.3,&
               &       ' per cent')
9000      FORMAT (/,'    int.-vac. boundary (density): ivac=',i2,/,t10,&
               &       'int-coord',t36,'cart-coord',t57,' inter. ',t69,' vacuum ')
9005      FORMAT (/,'    int.-vac. boundary (density): ivac=',i2,/,t10,&
               &       'int-coord',t36,'cart-coord',t64,' inter. ',t87,' vacuum ')
9010      FORMAT (/,'    int.-m.t. boundary (density): atom type=',i2,/,t10,&
               &       'int-coord',t36,'cart-coord',t57,' inter. ',t69,'  m. t. ')
9015      FORMAT (/,'    int.-m.t. boundary (density): atom type=',i2,/,t10,&
               &       'int-coord',t36,'cart-coord',t64,' inter. ',t87,'  m. t. ')
          RETURN
        END SUBROUTINE checkdop
      END MODULE m_checkdop
