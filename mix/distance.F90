!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_distance
contains
  SUBROUTINE distance(irank,vol,jspins,fsm,sm,iter,outden,results)
    use m_types
    use m_types_mixvector
    use m_xmlOutput
 
    implicit none
    integer,intent(in)             :: irank,jspins,iter
    real,intent(in)                :: vol
    type(t_mixvector),INTENT(IN)   :: fsm,sm
    TYPE(t_potden),INTENT(IN)      :: outden
    TYPE(t_results),INTENT(INOUT)  :: results
    
    integer         ::js
    REAL            :: dist(6) !1:up,2:down,3:spinoff,4:total,5:magnet,6:noco
    type(t_mixvector)::fmMet,fsm_mag
    character(len=100)::attributes(2)
    
    CALL fmMet%alloc()
    if (input%jspins==2) CALL fsm_mag%alloc()
 
    ! calculate Magnetisation-difference
    CALL fsm_mag%from_density(outden,swapspin=.true.)
    fsm_mag=fsm_mag-sm(it)

    ! Apply metric w to fsm and store in fmMet:  w |fsm>
    fmMet=fsm%apply_metric()
  
    dist(:) = 0.0
    DO js = 1,input%jspins
       dist(js) = fsm%multiply_dot_mask(fmMet,(/.true.,.true.,.true.,.false./),js)
    END DO
    dist(6) = fsm%multiply_dot_mask(fmMet,(/.true.,.true.,.true.,.false./),3)
    IF (input%jspins.EQ.2) then
       dist(3) = fsm_mag%multiply_dot_mask(fmMet,(/.true.,.true.,.true.,.false./),1)+&
            fsm_mag%multiply_dot_mask(fmMet,(/.true.,.true.,.true.,.false./),2)
       dist(4) = dist(1) + dist(2) + 2.0e0*dist(3)
       dist(5) = dist(1) + dist(2) - 2.0e0*dist(3)
    endif
    results%last_distance=maxval(1000*SQRT(ABS(dist/vol)))
    if (irank>1) return

    !calculate the distance of charge densities for each spin
    CALL openXMLElement('densityConvergence',(/'units'/),(/'me/bohr^3'/))
    
    DO js = 1,input%jspins         
       attributes = ''
       WRITE(attributes(1),'(i0)') js
       WRITE(attributes(2),'(f20.10)') 1000*SQRT(ABS(dist(js)/vol))
       CALL writeXMLElementForm('chargeDensity',(/'spin    ','distance'/),attributes,reshape((/4,8,1,20/),(/2,2/)))
       WRITE ( 6,FMT=7900) js,iter,1000*SQRT(ABS(dist(js)/vol))
    END DO
    
    IF (abs(dist(6))>1E-15) WRITE (6,FMT=7900) 3,iter,1000*SQRT(ABS(dist(6)/vol))
    
    !calculate the distance of total charge and spin density
    !|rho/m(o) - rho/m(i)| = |rh1(o) -rh1(i)|+ |rh2(o) -rh2(i)| +/_
    !                        +/_2<rh2(o) -rh2(i)|rh1(o) -rh1(i)>
    IF (input%jspins.EQ.2) THEN
       CALL writeXMLElementFormPoly('overallChargeDensity',(/'distance'/),&
            (/1000*SQRT(ABS(dist(4)/vol))/),reshape((/10,20/),(/1,2/)))
       CALL writeXMLElementFormPoly('spinDensity',(/'distance'/),&
            (/1000*SQRT(ABS(dist(5)/vol))/),reshape((/19,20/),(/1,2/)))
       WRITE ( 6,FMT=8000) iter,1000*SQRT(ABS(dist(4)/vol))
       WRITE ( 6,FMT=8010) iter,1000*SQRT(ABS(dist(5)/vol))
       
       !dist/vol should always be >= 0 ,
       !but for dist=0 numerically you might obtain dist/vol < 0
       !(e.g. when calculating non-magnetic systems with jspins=2).
    END IF
    CALL closeXMLElement('densityConvergence')


7900  FORMAT (/,'---->    distance of charge densities for spin ',i2,'                 it=',i5,':',f13.6,' me/bohr**3')
8000 FORMAT (/,'---->    distance of charge densities for it=',i5,':', f13.6,' me/bohr**3')
8010 FORMAT (/,'---->    distance of spin densities for it=',i5,':', f13.6,' me/bohr**3')
8020 FORMAT (4d25.14)
8030  FORMAT (10i10)
  end SUBROUTINE distance
end module m_distance
