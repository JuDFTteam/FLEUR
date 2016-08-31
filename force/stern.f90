MODULE m_stern
  !     **************************************************************
  !     returns star of recipocal space vector g
  !     called by force_a8 - APW+LO package
  !     *************************************************************
CONTAINS
  SUBROUTINE stern(sym,cell,g, nst,stg,taup,gl,rstg)

    USE m_constants, ONLY : tpi_const
    USE m_types
    IMPLICIT NONE
    TYPE(t_sym),INTENT(IN)   :: sym
    TYPE(t_cell),INTENT(IN)   :: cell
    !     ..
    !     .. Arguments
    INTEGER, INTENT (IN)  :: g(3) 

    INTEGER, INTENT (OUT) :: nst,stg(3,sym%nop)
    REAL,    INTENT (OUT) :: gl,rstg(3,sym%nop)
    COMPLEX, INTENT (OUT) :: taup(sym%nop)
    !     ..
    !     .. Local Variables
    INTEGER               :: i,m,j,k,l,ind(sym%nop)
    REAL                  :: tk,s,rg(3)

    ind(1:sym%nop)  = 0
    taup(1:sym%nop) = 0.0
    nst = 0                                                             

    rg(:) = REAL( g(:) )
    gl = SQRT( DOT_PRODUCT(rg,MATMUL(rg,cell%bbmat)))

    i_loop:DO i = 1,sym%nop
       tk=0.                                                          
       DO j=1,3                                                     
          tk=tk+sym%tau(j,i)*g(j)*tpi_const                                     
          k=0                                                         
          DO l=1,3
             k=sym%mrot(l,j,i)*g(l)+k                                       
          ENDDO
          stg(j,i)=k                                                  
       ENDDO
       IF (nst.NE.0) THEN                                              
          DO m = 1,nst                                                   
               IF (ALL(stg(:,m)==stg(:,i))) THEN
                  ind(m)=ind(m)+1
                  taup(m)=taup(m) + CMPLX(COS(tk),SIN(tk))
                  CYCLE i_loop
               ENDIF
            ENDDO
         ENDIF 
         nst=nst+1
         stg(:,nst)=stg(:,i)
         DO j = 1,3
            rstg(j,nst) = DOT_PRODUCT(stg(:,nst),cell%bmat(:,j))
         ENDDO
         ind(nst)  = 1
         taup(nst) = CMPLX(COS(tk),SIN(tk))
      ENDDO i_loop                                                        
                              
      taup(:nst)=taup(:nst)/ind(:nst)                                            

      RETURN
   END SUBROUTINE stern
END MODULE m_stern
