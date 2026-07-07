!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_current4 
      USE m_constants, ONLY: oUnit
      IMPLICIT NONE
      PRIVATE 
      PUBLIC::gf_current4 
      LOGICAL :: new=.false.,firstcall=.true.

      type t_curr
         integer:: num_currents,totalnum_currents
         logical,allocatable:: bardeen(:),landauer(:)
         integer,allocatable:: layer(:)
         real,allocatable   :: imag(:)
      endtype

      type(t_curr)::curr


      CONTAINS 

    subroutine gf_current4_read(num_layers)
      use m_juDFT
      integer,intent(in) :: num_layers
      integer:: n,layer
      logical:: landauer,bardeen
      real   :: imag

      NAMELIST /current/layer,landauer,bardeen,imag
      !count the entries of the gf_current4 file
      open(99,file="gf_current4")
      n=0
      DO
         read(99,current,end=100)
         n=n+1
      ENDDO
  100 if (n==0) call juDFT_error("Incorrect gf_current4 file")
      rewind(99)

      !now we can read
      allocate(curr%bardeen(n),curr%landauer(n),curr%imag(n),curr%layer(n))
      curr%num_currents=n
      curr%totalnum_currents=0
      DO n=1,curr%num_currents
          landauer=.true.;bardeen=.false.;imag=0;layer=1
          read(99,current)
          curr%landauer(n)=landauer
          curr%bardeen(n)=bardeen
          if ((landauer.and.bardeen).or..not.(landauer.or.bardeen)) call juDFT_error("You must specify either landauer or bardeen",calledby="gf_current4")
          if (landauer) curr%totalnum_currents= curr%totalnum_currents+2
          if (bardeen)  curr%totalnum_currents= curr%totalnum_currents+3
          curr%layer(n)=layer
          if (layer<1.or.layer>num_layers) call juDFT_error("Impossible choice of layer index",calledby="gf_current4")
          curr%imag(n)=imag

          write(oUnit,*) n,curr%totalnum_currents,curr%layer(n),curr%imag(n)

       ENDDO
       close(99)
      END subroutine gf_current4_read

      SUBROUTINE gf_current4(                                           &
     &           layers,l_noco,en,nk,jspin,lapw,lapw_gf,                        &
     &           bkpts,sym,cell,gfinp,mpi)
!************************************************                       
!     Calculate the current on the basis of                             
!     the embedding potentials of the same planes.                      
!     This is done for all dividing planes of                           
!     the composite system.                                             
!     Frank Freimuth, November 2007                                     
!************************************************                       
      USE m_gf_iotmat 
      USE m_gf_types 
      USE m_gf_embedding 
      USE m_gf_propaemb 
      USE m_gf_writetrans,ONLY:writetrans 
      USE m_gf_propagate_embpot 
      USE m_gf_math 
                                                                        
      IMPLICIT NONE 
      TYPE(t_layers),INTENT(IN)::layers 
      type(t_gfmpi),intent(in)   :: mpi
      LOGICAL,INTENT(IN)       :: l_noco 
      INTEGER,INTENT(IN)::en 
      INTEGER,INTENT(IN)::nk 
      INTEGER,INTENT(IN)::jspin 
      TYPE(t_lapw),INTENT(IN)::lapw 
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      REAL,INTENT(IN)::bkpts(:,:) 
      TYPE(t_sym),INTENT(IN)::sym 
      TYPE(t_cell),INTENT(IN)::cell 
      TYPE(t_embinp),INTENT(IN)::gfinp 

      if (firstcall) THEN
              inquire(file="gf_current4",exist=new)
              firstcall=.false.
              if (new) call gf_current4_read(layers%num_layers)
      ENDIF
      if (.not.new) THEN
         call gf_current4_old( &
     &           layers,l_noco,en,nk,jspin,lapw,lapw_gf,                        &
     &           bkpts,sym,cell,gfinp,mpi)
         return
      endif
      !
      call gf_current4_new( &
     &           layers,l_noco,en,nk,jspin,lapw,lapw_gf,                        &
     &           bkpts,sym,cell,gfinp,mpi)

      END subroutine gf_current4

      SUBROUTINE gf_current4_new(                                           &
     &           layers,l_noco,en,nk,jspin,lapw,lapw_gf,                        &
     &           bkpts,sym,cell,gfinp,mpi)
!************************************************
!     Calculate the current on the basis of
!     the embedding potentials of the same planes.
!     This is done for all dividing planes of
!     the composite system.
!     Frank Freimuth, November 2007
!************************************************
      USE m_gf_iotmat
      USE m_gf_types
      USE m_gf_embedding
      USE m_gf_propaemb
      USE m_gf_writetrans,ONLY:writetrans
      USE m_gf_propagate_embpot
      USE m_gf_math

      IMPLICIT NONE
      TYPE(t_layers),INTENT(IN)::layers
      type(t_gfmpi),intent(in)   :: mpi
      LOGICAL,INTENT(IN)       :: l_noco
      INTEGER,INTENT(IN)::en
      INTEGER,INTENT(IN)::nk
      INTEGER,INTENT(IN)::jspin
      TYPE(t_lapw),INTENT(IN)::lapw
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      REAL,INTENT(IN)::bkpts(:,:)
      TYPE(t_sym),INTENT(IN)::sym
      TYPE(t_cell),INTENT(IN)::cell
      TYPE(t_embinp),INTENT(IN)::gfinp
                                                                        
      COMPLEX,ALLOCATABLE:: g1(:,:),g2(:,:) 
      INTEGER            :: layer,n,c
      REAL               :: j(curr%totalnum_currents)

      ALLOCATE( g1(lapw_gf%nv2_tot,lapw_gf%nv2_tot) )
      ALLOCATE( g2(lapw_gf%nv2_tot,lapw_gf%nv2_tot) )

      n = 1
      DO c=1,curr%num_currents
         layer=curr%layer(c)
         CALL gf_getemb2(g1,1,layer+1,en,nk,jspin,lapw,lapw_gf)
         CALL gf_getemb2(g2,2,layer,  en,nk,jspin,lapw,lapw_gf)
         IF (curr%landauer(c)) THEN
            CALL gf_landauer1plane(l_noco,                                 &
     &        g1,                                                       &
     &        g2,                                                       &
     &        j(n:n+1),curr%imag(c))
             n=n+2
         ENDIF
         IF (curr%bardeen(c)) THEN
             CALL gf_bardeen1plane(g1,g2,j(n:n+2))
             n = n+3
         ENDIF
      ENDDO
      CALL writetrans(en,nk,jspin,bkpts,sym,cell,n,j,mpi)

      END SUBROUTINE gf_current4_new


      SUBROUTINE gf_current4_old(                                           &
     &           layers,l_noco,en,nk,jspin,lapw,lapw_gf,                        &
     &           bkpts,sym,cell,gfinp,mpi)
!************************************************
!     Calculate the current on the basis of
!     the embedding potentials of the same planes.
!     This is done for all dividing planes of
!     the composite system.
!     Frank Freimuth, November 2007
!************************************************
      USE m_gf_iotmat
      USE m_gf_types
      USE m_gf_embedding
      USE m_gf_propaemb
      USE m_gf_writetrans,ONLY:writetrans
      USE m_gf_propagate_embpot
      USE m_gf_math

      IMPLICIT NONE
      TYPE(t_layers),INTENT(IN)::layers
      type(t_gfmpi),intent(in)   :: mpi
      LOGICAL,INTENT(IN)       :: l_noco
      INTEGER,INTENT(IN)::en
      INTEGER,INTENT(IN)::nk
      INTEGER,INTENT(IN)::jspin
      TYPE(t_lapw),INTENT(IN)::lapw
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      REAL,INTENT(IN)::bkpts(:,:)
      TYPE(t_sym),INTENT(IN)::sym
      TYPE(t_cell),INTENT(IN)::cell
      TYPE(t_embinp),INTENT(IN)::gfinp

      COMPLEX,ALLOCATABLE:: g1(:,:),g2(:,:)
      INTEGER            :: layer,n
      REAL               :: j(5*layers%num_layers-5)

      ALLOCATE( g1(lapw_gf%nv2_tot,lapw_gf%nv2_tot) )
      ALLOCATE( g2(lapw_gf%nv2_tot,lapw_gf%nv2_tot) )

      n = 1
      DO layer = 1,layers%num_layers-1
         CALL gf_getemb2(g1,1,layer+1,en,nk,jspin,lapw,lapw_gf)
         CALL gf_getemb2(g2,2,layer,  en,nk,jspin,lapw,lapw_gf)
         CALL gf_landauer1plane(l_noco,                                 &
     &        g1,                                                       &
     &        g2,                                                       &
     &        j(n:n+1),0.0)
         CALL gf_bardeen1plane(g1,g2,j(n+2:n+4))
         n = n+5
      ENDDO
      CALL writetrans(en,nk,jspin,bkpts,sym,cell,n,j,mpi)

      END SUBROUTINE gf_current4_old
      !<-- S: gf_Landauer1Plane                                         
      SUBROUTINE gf_Landauer1Plane(l_noco,                              &
     &                             g1,g2,                               &
     &                             j,imaginary)
!c********************************************************************* 
!c     subroutine to calculate the current from two embedding           
!c     potentials on the same plane                                     
!c                                                                      
!c                                      Daniel Wortmann                 
!c********************************************************************* 
      USE m_gf_math 
      IMPLICIT NONE 
      LOGICAL,INTENT(IN)  :: l_noco 
      COMPLEX,INTENT(IN)  :: g1(:,:),g2(:,:) 
      REAL,INTENT(OUT)    :: j(2) 
      real,intent(in)     :: imaginary
                                                                        
      INTEGER             :: n 
      COMPLEX             :: G(SIZE(g1,1),SIZE(g1,1)) 
      COMPLEX             :: A(SIZE(g1,1),SIZE(g1,1)),B(SIZE(g1,1)      &
     &     ,SIZE(g1,1))                                                 
                                                                        

      IF(.FALSE.)THEN
         G = mat_inverse(G1+G2)
         j(1) = 2.0*REAL(trace(MATMUL(MATMUL(imag2d(G1),G),             &
     &     matmul(imag2d(g2)                                            &
     &     ,TRANSPOSE(CONJG(G))))))                                     
                                                                        
      ELSE 
        ! longer version                                                
        ! Re(G1 G12 G2 G21^* - G1 G12 G2^* G21^*)                       
        G = mat_inverse(G1+G2+cmplx(0.0,imaginary))
        A=matmul(matmul(g1,g),g2)
        B=matmul(matmul(g1,g),transpose(conjg(g2)))
        if (abs(imaginary)>epsilon(1.0)) G = mat_inverse(G1+G2-cmplx(0.0,imaginary))
        A=matmul(A,transpose(conjg(g))) 
        B=matmul(B,transpose(conjg(g))) 
                                                                        
        j(1) =-2.*REAL(trace(A-B)) 
      ENDIF 
                                                                        
      IF (l_noco) THEN 
         n = SIZE(g1,1)/2 
         j(2) = 2.0*REAL(trace(MATMUL(MATMUL(imag2d(G1(:n,:n)),G(:n,:n))&
     &        ,MATMUL(imag2d(g2(:n,:n)),TRANSPOSE(CONJG(G(:n,:n)))))))  
      ELSE 
         j(2) = 0.0 
      ENDIF 
                                                                        
      END SUBROUTINE gf_landauer1plane 
      !>                                                                
      !<-- S: gf_Bardeen1Plane                                          
      SUBROUTINE gf_Bardeen1Plane(                                      &
     &                             g1,g2,                               &
     &                             j)                                   
!c********************************************************************* 
!c     subroutine to calculate the current from two embedding           
!c     potentials on the same plane                                     
!c                                                                      
!c                                      Daniel Wortmann                 
!c********************************************************************* 
      USE m_gf_math 
      IMPLICIT NONE 
      COMPLEX,INTENT(IN)  :: g1(:,:),g2(:,:) 
      REAL,INTENT(OUT)    :: j(3) 
                                                                        
      COMPLEX             :: A(SIZE(g1,1),SIZE(g1,1)),g_1(SIZE(g1,1)    &
     &     ,SIZE(g1,1)),g_2(SIZE(g1,1),SIZE(g1,1)),sigma(SIZE(g1,1)     &
     &     ,SIZE(g1,1))                                                 
      !calculate conductance by using the left side for both electrodes 
                                !make sigma real                        
      Sigma = 2*(g1-IMAG2d(g1)) 
      G_1   = mat_inverse(G1+G2-imag2d(g2)) 
      A   = MATMUL(sigma,imag2d(g_1)) 
      j(1) = 0.5*trace(MATMUL(a,a)) 
      !calculate conductance by using the right side for both electrodes
      sigma = 2*(g2-IMAG2d(g2)) 
      G_2   = mat_inverse(G1+G2-imag2d(g1)) 
      A     = MATMUL(sigma,imag2d(g_2)) 
      j(2)  = 0.5*trace(MATMUL(a,a)) 
      !calculate conductance by using both sides                        
      sigma = (g2-IMAG2d(g2)+g1-IMAG2d(g1)) 
      j(3)  = 0.5*trace(MATMUL(sigma,MATMUL(imag2d(g_1),MATMUL(sigma    &
     &     ,imag2d(g_2)))))                                             
      !WRITE(*,*) trace(sigma),trace(g_2)
      END SUBROUTINE gf_bardeen1plane 
      !>                                                                
                                                                        
      END                                           
