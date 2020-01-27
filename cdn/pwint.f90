      MODULE m_pwint
!     ******************************************************************
!     calculate the integral of a star function over the interstial    *
!     region              c.l.fu                                       *
!     ******************************************************************
      CONTAINS
      SUBROUTINE pwint(&
     &                 stars,atoms,sym,oneD,&
     &                 cell,ng,&
     &                 x)

      USE m_spgrot
      USE m_od_cylbes
      use m_juDFT
      USE m_types 
      USE m_constants
      IMPLICIT NONE
!     ..
!     .. Scalar Arguments ..
      TYPE(t_stars),INTENT(IN) :: stars
      TYPE(t_atoms),INTENT(IN) :: atoms
      TYPE(t_sym),INTENT(IN)   :: sym
      TYPE(t_oneD),INTENT(IN)  :: oneD
      TYPE(t_cell),INTENT(IN)  :: cell
      INTEGER,INTENT(IN)       :: ng
      COMPLEX, INTENT (OUT):: x
!     ..
!     .. Array Arguments ..
!-odim
!+odim
!     ..
!     .. Local Scalars ..
      COMPLEX s1,sfs
      REAL arg,g,s,srmt,gr,fJ
      INTEGER ig2d,ig3d,n,nn,na,ii
!     ..
!     .. Local Arrays ..
      COMPLEX ph(sym%nop)
      INTEGER kr(3,sym%nop)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC cmplx,cos,exp,sin
!     ..
      ig3d = stars%ig(stars%kv3(1,ng),stars%kv3(2,ng),stars%kv3(3,ng))
      IF (ig3d.EQ.0) THEN
         x = (0.,0.)
         RETURN
      END IF
      IF (ig3d.EQ.1) THEN
         x = cmplx(cell%volint,0.0)
         RETURN
      ELSE

         IF (oneD%odi%d1) THEN
            IF (stars%kv3(3,ng).EQ.0) THEN
               g = (stars%kv3(1,ng)*cell%bmat(1,1) + stars%kv3(2,ng)*cell%bmat(2,1))**2 +&
     &             (stars%kv3(1,ng)*cell%bmat(1,2) + stars%kv3(2,ng)*cell%bmat(2,2))**2
               gr = sqrt(g)
               g  = gr*cell%z1
               CALL od_cylbes(1,g,fJ)
               x = cmplx(2*cell%vol*fJ/g,0.0)
            ELSE
               x = (0.0,0.0)
            END IF
          ELSE
             ig2d = stars%ig2(ig3d)
             IF (ig2d.EQ.1) THEN
                g = stars%kv3(3,ng)*cell%bmat(3,3)*cell%z1
                x = cmplx(cell%vol*sin(g)/g,0.0)
             ELSE
                x = (0.0,0.0)
             END IF
          END IF

      END IF
!     -----> sphere contributions
      s = stars%sk3(ig3d)
      na = 1

      IF (.NOT.oneD%odi%d1) THEN
         CALL spgrot(&
     &           sym%nop,sym%symor,sym%mrot,sym%tau,sym%invtab,&
     &           stars%kv3,&
     &           kr,ph)
          DO  n = 1,atoms%ntype
            srmt = s*atoms%rmt(n)
            sfs = (0.0,0.0)
            DO  nn = 1,sym%nop
               arg = tpi_const * dot_product(real(kr(:,nn)),atoms%taual(:,na))
               sfs = sfs + cmplx(cos(arg),sin(arg))*ph(nn)
            ENDDO
            sfs = sfs/sym%nop
!     -----3*ji(gr)/gr term
            s1 = 3.* (sin(srmt)/srmt-cos(srmt))/ (srmt*srmt)
            x = x - atoms%volmts(n)*atoms%neq(n)*s1*sfs
            na = na + atoms%neq(n)
         ENDDO
      ELSE
!-odim
         DO 21 n = 1,atoms%ntype
            DO ii = 1,atoms%neq(n)
               srmt = s*atoms%rmt(n)
               CALL spgrot(&
     &              sym%nop,sym%symor,sym%mrot,sym%tau,sym%invtab,&
     &              stars%kv3,&
     &              kr,ph)
               sfs = (0.0,0.0)
               DO nn = 1,sym%nop
                 arg = tpi_const * dot_product(real(kr(:,nn)),atoms%taual(:,na))
                 sfs = sfs + cmplx(cos(arg),sin(arg))*ph(nn)
               ENDDO
               sfs = sfs/sym%nop
!     -----3*ji(gr)/gr term
               s1 = 3.* (sin(srmt)/srmt-cos(srmt))/ (srmt*srmt)
               x = x - atoms%volmts(n)*s1*sfs
               na = na + 1
            END DO
 21      CONTINUE
!+odim
      ENDIF

      END SUBROUTINE pwint
      SUBROUTINE pwint_all(&
     &                 stars,atoms,sym,oneD,&
     &                 cell,x_start,x_end,&
     &                 x)

      USE m_spgrot
      USE m_od_cylbes
      use m_juDFT
      USE m_types 
      USE m_constants
      IMPLICIT NONE
!     ..

      TYPE(t_stars),INTENT(IN) :: stars
      TYPE(t_atoms),INTENT(IN) :: atoms
      TYPE(t_sym),INTENT(IN)   :: sym
      TYPE(t_oneD),INTENT(IN)  :: oneD
      TYPE(t_cell),INTENT(IN)  :: cell
      INTEGER, INTENT (IN) :: x_start,x_end
      COMPLEX, INTENT (OUT):: x(:)
!     ..
!-odim
!+odim
!     ..
!     .. Local Scalars ..
      COMPLEX s1,sfs
      REAL arg,g,s,srmt,gr,fJ
      INTEGER ig2d,ig3d,n,nn,na,ii,ng
!     ..
!     .. Local Arrays ..
      COMPLEX ph(sym%nop)
      INTEGER kr(3,sym%nop)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC cmplx,cos,exp,sin
!     ..
      
!$OMP PARALLEL DO default(shared)  &
!$OMP PRIVATE(ng,ig3d,g,gr,fj,ig2d,s,na,kr,ph,n)&
!$OMP PRIVATE(srmt,nn,sfs,arg,s1,ii)
      starloop:DO ng=x_start,x_end
      ig3d = stars%ig(stars%kv3(1,ng),stars%kv3(2,ng),stars%kv3(3,ng))
      IF (ig3d.EQ.0) THEN
         x(ng) = (0.,0.)
         cycle starloop
      END IF

      IF (ig3d.EQ.1) THEN
         x(ng) = cmplx(cell%volint,0.0)
         cycle starloop
      ELSE

         IF (oneD%odi%d1) THEN
            IF (stars%kv3(3,ng).EQ.0) THEN
               g = (stars%kv3(1,ng)*cell%bmat(1,1) + stars%kv3(2,ng)*cell%bmat(2,1))**2 +&
     &             (stars%kv3(1,ng)*cell%bmat(1,2) + stars%kv3(2,ng)*cell%bmat(2,2))**2
               gr = sqrt(g)
               g  = gr*cell%z1
               CALL od_cylbes(1,g,fJ)
               x(ng) = cmplx(2*cell%vol*fJ/g,0.0)
            ELSE
               x(ng) = (0.0,0.0)
            END IF
          ELSE
             ig2d = stars%ig2(ig3d)
             IF (ig2d.EQ.1) THEN
                g = stars%kv3(3,ng)*cell%bmat(3,3)*cell%z1
                x(ng) = cmplx(cell%vol*sin(g)/g,0.0)
             ELSE
                x(ng) = (0.0,0.0)
             END IF
          END IF

      END IF
!     -----> sphere contributions
      s = stars%sk3(ig3d)
      na = 1

      IF (.NOT.oneD%odi%d1) THEN
         CALL spgrot(&
     &           sym%nop,sym%symor,sym%mrot,sym%tau,sym%invtab,&
     &           stars%kv3(:,ng),&
     &           kr,ph)
          DO  n = 1,atoms%ntype
            srmt = s*atoms%rmt(n)
            sfs = (0.0,0.0)
            DO  nn = 1,sym%nop
             arg = tpi_const* (kr(1,nn)*atoms%taual(1,na)+kr(2,nn)*atoms%taual(2,na)+&
     &              kr(3,nn)*atoms%taual(3,na))
               sfs = sfs + exp(cmplx(0.0,arg))*ph(nn)
            ENDDO
            sfs = sfs/sym%nop
!     -----3*ji(gr)/gr term
            s1 = 3.* (sin(srmt)/srmt-cos(srmt))/ (srmt*srmt)
            x(ng) = x(ng) - atoms%volmts(n)*atoms%neq(n)*s1*sfs
            na = na + atoms%neq(n)
         ENDDO
      ELSE
!-odim
         DO 21 n = 1,atoms%ntype
            DO ii = 1,atoms%neq(n)
               srmt = s*atoms%rmt(n)
               CALL spgrot(&
     &              sym%nop,sym%symor,sym%mrot,sym%tau,sym%invtab,&
     &              stars%kv3(:,ng),&
     &              kr,ph)
               sfs = (0.0,0.0)
               DO 11 nn = 1,sym%nop
                  arg = tpi_const* (kr(1,nn)*atoms%taual(1,na)+&
     &                 kr(2,nn)*atoms%taual(2,na)+&
     &                 kr(3,nn)*atoms%taual(3,na))
                  sfs = sfs + exp(cmplx(0.0,arg))*ph(nn)
 11            CONTINUE
               sfs = sfs/sym%nop
!     -----3*ji(gr)/gr term
               s1 = 3.* (sin(srmt)/srmt-cos(srmt))/ (srmt*srmt)
               x(ng) = x(ng) - atoms%volmts(n)*s1*sfs
               na = na + 1
            END DO
 21      CONTINUE
!+odim
      ENDIF
      ENDDO starloop
!$OMP end parallel do

      END SUBROUTINE pwint_all
      END MODULE m_pwint
