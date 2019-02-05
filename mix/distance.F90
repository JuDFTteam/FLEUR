!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_distance
contains
  SUBROUTINE distance(irank,fsm)

    
    REAL            :: dist(6) !1:up,2:down,3:total,4:
 


    CALL fsm_up%alloc()
    IF (input%jspin
    fsm_up=
    
    

    dist(:) = 0.0
    ! Apply metric w to fsm and store in fmMet:  w |fsm>
 
      !calculate the distance of charge densities for each spin
      IF(hybrid%l_calhf) THEN
         CALL openXMLElement('densityConvergence',(/'units  ','comment'/),(/'me/bohr^3','HF       '/))
      ELSE
         CALL openXMLElement('densityConvergence',(/'units'/),(/'me/bohr^3'/))
      END IF

      DO js = 1,input%jspins
         dist(js) = dot_product(fsm%vec(nmaph*(js-1)+1:nmaph*(js-1)+nmaph),fmMet%vec(nmaph*(js-1)+1:nmaph*(js-1)+nmaph))
         
         attributes = ''
         WRITE(attributes(1),'(i0)') js
         WRITE(attributes(2),'(f20.10)') 1000*SQRT(ABS(dist(js)/cell%vol))
         CALL writeXMLElementForm('chargeDensity',(/'spin    ','distance'/),attributes,reshape((/4,8,1,20/),(/2,2/)))
         IF( hybrid%l_calhf ) THEN
            WRITE ( 6,FMT=7901) js,inDen%iter,1000*SQRT(ABS(dist(js)/cell%vol))
         ELSE
            WRITE ( 6,FMT=7900) js,inDen%iter,1000*SQRT(ABS(dist(js)/cell%vol))
         END IF
      END DO
      IF (noco%l_noco) dist(6) = dot_product(fsm%vec(nmaph*2+1:nmaph*2+nmap-2*nmaph),fmMet%vec(nmaph*2+1:nmaph*2+nmap-2*nmaph))
      IF (noco%l_noco) WRITE (6,FMT=7900) 3,inDen%iter,1000*SQRT(ABS(dist(6)/cell%vol))

      !calculate the distance of total charge and spin density
      !|rho/m(o) - rho/m(i)| = |rh1(o) -rh1(i)|+ |rh2(o) -rh2(i)| +/_
      !                        +/_2<rh2(o) -rh2(i)|rh1(o) -rh1(i)>
      IF (input%jspins.EQ.2) THEN
         dist(3) = fsm.dot.fmMet
         dist(4) = dist(1) + dist(2) + 2.0e0*dist(3)
         dist(5) = dist(1) + dist(2) - 2.0e0*dist(3)
         CALL writeXMLElementFormPoly('overallChargeDensity',(/'distance'/),&
                                      (/1000*SQRT(ABS(dist(4)/cell%vol))/),reshape((/10,20/),(/1,2/)))
         CALL writeXMLElementFormPoly('spinDensity',(/'distance'/),&
                                      (/1000*SQRT(ABS(dist(5)/cell%vol))/),reshape((/19,20/),(/1,2/)))
         IF( hybrid%l_calhf ) THEN
            WRITE ( 6,FMT=8001) inDen%iter,1000*SQRT(ABS(dist(4)/cell%vol))
            WRITE ( 6,FMT=8011) inDen%iter,1000*SQRT(ABS(dist(5)/cell%vol))
         ELSE
            WRITE ( 6,FMT=8000) inDen%iter,1000*SQRT(ABS(dist(4)/cell%vol))
            WRITE ( 6,FMT=8010) inDen%iter,1000*SQRT(ABS(dist(5)/cell%vol))
         END IF

         !dist/vol should always be >= 0 ,
         !but for dist=0 numerically you might obtain dist/vol < 0
         !(e.g. when calculating non-magnetic systems with jspins=2).
      END IF
      results%last_distance=maxval(1000*SQRT(ABS(dist/cell%vol)))
      DEALLOCATE (smMet,fmMet)
      CALL closeXMLElement('densityConvergence')
