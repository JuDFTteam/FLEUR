      MODULE m_gkptwgt
      CONTAINS
        SUBROUTINE gkptwgt(&
             &                   kpts,cell)
          !                                                             
          !     this subroutine generates the weight factor of a k-point 
          !     in the irreducible wedge of a 2d-brillouin zone         
          !     ver are vertex coordinates in counter clockwise order and 
          !     in units of pi/a(1) and pi/a(2)                              
          !     it is assumed that each k-point has the same 'area'. D.S.Wang     
          !     
          !     changed by                        Stefan Bl"ugel, IFF, Jan.96
          !                                                           
          USE m_constants
          USE m_types
          USE m_juDFT
          IMPLICIT NONE
          TYPE(t_cell),INTENT(IN)   :: cell
          TYPE(t_kpts),INTENT(INOUT):: kpts

          !                                                        
          !     .. was an Argument
          REAL    :: wt
          !     ..
          !     .. Array Arguments ..
          !     ..
          !     .. Local Scalars ..
          REAL cross,dot,eps,x1,x2,y1,y2
          REAL s1,s2
          INTEGER i,j,ikpt,nver
          !     ..
          !     .. Local Arrays ..
          REAL ver(2,4),dummy(2,kpts%nkpt)
          !     ..
          !     .. Intrinsic Functions ..
          INTRINSIC abs,atan2
          !     ..
          !     .. Data Statements ..
          DATA ver/0.e0,0.e0,1.e0,0.e0,1.e0,1.e0,0.e0,1.e0/
          DATA eps/1.e-6/
          !     ..

          !      nver = 3
          IF ( cell%latnam.EQ.'squ' ) THEN
             nver = 3
          ELSEIF ( cell%latnam.EQ.'p-r' .OR. cell%latnam.EQ.'c-r' ) THEN
             nver = 4
          ELSEIF ( cell%latnam.EQ.'hex' ) THEN
             nver = 3
             ver(2,3) = 1./3.
          ELSEIF ( cell%latnam.EQ.'hx3' .OR. cell%latnam.EQ. 'obl' ) THEN
             CALL juDFT_error("weights for hx3 or obl not defined" ,calledby&
                  &        ="gkptwgt")
          ENDIF
          !                                                          
          !     transform from internal coordinates to xy-coordinates
          !
          !                            changed by shz Feb.96
          DO 10 ikpt = 1 , kpts%nkpt
             kpts%wtkpt(ikpt) = 0
             dummy(1,ikpt)=kpts%bk(1,ikpt)
             dummy(2,ikpt)=kpts%bk(2,ikpt)
             s1 = 0.0
             s2 = 0.0
             DO i = 1,2 
                s1 = s1+cell%bmat(i,1)*kpts%bk(i,ikpt)
                s2 = s2+cell%bmat(i,2)*kpts%bk(i,ikpt)
             ENDDO
             IF (cell%latnam.EQ.'hex') THEN
                kpts%bk(1,ikpt) = s1*cell%amat(2,2)/tpi_const
                kpts%bk(2,ikpt) = s2*cell%amat(1,1)/pi_const
             ELSE
                kpts%bk(1,ikpt) = s1*cell%amat(1,1)/pi_const
                kpts%bk(2,ikpt) = s2*cell%amat(2,2)/pi_const
             ENDIF
10        ENDDO

          DO 20 ikpt = 1 , kpts%nkpt
             DO 30 i = 1,nver
                x1 = ( ver(1,i)-kpts%bk(1,ikpt) ) / cell%amat(1,1)
                y1 = ( ver(2,i)-kpts%bk(2,ikpt) ) / cell%amat(2,2)
                j  = i + 1
                IF ( j.GT.nver ) j = 1
                x2 = ( ver(1,j)-kpts%bk(1,ikpt) ) / cell%amat(1,1)
                y2 = ( ver(2,j)-kpts%bk(2,ikpt) ) / cell%amat(2,2)
                dot = x1*x2 + y1*y2
                cross = x1*y2 - y1*x2
                IF ( ABS(cross).GE.eps ) THEN
                   kpts%wtkpt(ikpt) = kpts%wtkpt(ikpt) + ATAN2(cross,dot)
                ENDIF
30           ENDDO
20        ENDDO
          !
          DO ikpt = 1 , kpts%nkpt
             kpts%wtkpt(ikpt) = kpts%wtkpt(ikpt) /tpi_const
          ENDDO
          !   
          wt = 0.0
          DO ikpt = 1,kpts%nkpt
             wt = wt + kpts%wtkpt(ikpt)
             kpts%bk(1,ikpt)=dummy(1,ikpt)
             kpts%bk(2,ikpt)=dummy(2,ikpt)
          ENDDO

          RETURN
        END SUBROUTINE gkptwgt
      END
