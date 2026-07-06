!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_CBS 
      use m_juDFT
      USE hdf5
      IMPLICIT NONE
      INTEGER(hid_t)::fid,re_k_var,im_k_var,spinp_var,curr_var
      PRIVATE 
      PUBLIC gf_CBS,outCBS_openfile,outCBS_closefile
      CONTAINS 
                                                                        
      !<--S: gf_CBS(nk,Nv2,En,                                          
                                                                        
      SUBROUTINE gf_CBS(layers,layer,nk,bk,En,jspin,gfinp,lapw,         &
     &     cell,T,mrot,mpi,T1,T2)                                       
!********************************************************************** 
!     * This SUBROUTINE takes the Transfer-Matrix and calculates        
!     *the CBS by calculating Eigenvectors and Eigenvalues              
!     *                                                                 
!     *                           Daniel Wortmann, Tokyo, 2001          
! *                                                                     
!********************************************************************** 
#include "juDFT_env.h"
      USE m_gf_types 
      USE m_gf_math, ONLY  : eigenvalues,tpi_const 
      USE m_gf_io2dmat 
      USE m_gf_embedding 
      IMPLICIT NONE 
                                                                        
      !<--Arguments                                                     
                                                                        
      TYPE(t_layers),INTENT(IN)  :: layers 
      INTEGER,INTENT(IN)         :: layer 
      INTEGER, INTENT(IN)       :: jspin, En,nk 
      REAL, INTENT(IN)          :: bk(2) 
      INTEGER, INTENT(IN)       :: mrot(:,:,:) 
      TYPE(t_lapw),INTENT(IN)   :: lapw 
      TYPE(t_mpi)               :: mpi 
      TYPE(t_cell),INTENT(IN)   :: cell 
      TYPE(t_gfinp), INTENT(IN) :: gfinp 
      COMPLEX, INTENT(INOUT)    :: T(2*Lapw%nv2_Tot,2*Lapw%nv2_Tot) 
      COMPLEX, INTENT(IN)       :: T1(2*Lapw%nv2_Tot,2*Lapw%nv2_Tot)    &
     &     ,T2(2*Lapw%nv2_Tot,2*Lapw%nv2_Tot)                           
                                                                        
      !>                                                                
      !<-- Locals                                                       
                                                                        
      INTEGER              :: n,state,n2,nn,nloc(2*lapw%nv2_tot)
      REAL                 :: bz 
                                                 !put into input-data la
      LOGICAL, PARAMETER   :: L_STATES = .FALSE. 
                                                 !Symmetry info         
      INTEGER ,ALLOCATABLE :: syms(:,:) 
                                                 !Eigenvalues+vecs      
      COMPLEX, ALLOCATABLE :: ew(:,:),ev(:,:) 
                                                 !Current of EV         
      REAL, ALLOCATABLE    :: curr(:),spinp(:,:)
      COMPLEX,ALLOCATABLE  :: logderiv(:,:)
                                                                    !cma
      COMPLEX              :: T_temp(2*Lapw%nv2_Tot,2*Lapw%nv2_Tot) 
      LOGICAL              :: l_noembtest 
      COMPLEX              :: vector(lapw%nv2_tot,2*lapw%nv2_tot) 
      INTEGER              :: layer_ind 
      complex              :: cnorm
                                                                        
      !>                                                                
      !<-- Allocate storage,initialize                                  
                                                                        
      ALLOCATE(ew(2*Lapw%nv2_Tot,2),spinp(0:3,2*lapw%nv2_tot),curr(2*Lapw%nv2_Tot),ev(2           &
     &     *Lapw%nv2_Tot,2*Lapw%nv2_Tot),logderiv(4,2*lapw%nv2_tot))
                                                                        
!      bz = tpi_const/(2*cell%Z1-gfinp%d)                               
       bz=tpi_const/sum(layers%c)
                                                                        
      !>                                                                
                                                                        
      !<--now calculate Eigenvectors and eigenvalues                    
                                                                        
      t_temp = T 
      CPP_juDFT_timestart("CBS:diagonalization")
      CALL eigenvalues(T_temp,ew(:,1),ev) 
      CPP_juDFT_timestop("CBS:diagonalization")
      ew(:,2) = CMPLX(0,bz/tpi_const)*LOG(ew(:,1)) 
!      do n=1,2*lapw%nv2_tot                                            
!         if(abs(ew(n,1)).lt.1.0e-10)                                   
!      enddo                                                            
      !>                                                                
      CPP_juDFT_timestart("CBS:current")
      CALL priv_current(en,nk,mpi,lapw,gfinp,ew,ev,curr) 
      CPP_juDFT_timestop("CBS:current")
                                                                        
                                                                        
      !<--For Debugging: write to stdout&STOP                                                                    
      IF ( L_STATES ) THEN 
         WRITE (*,*) '2D-Eigenstates:' 
         DO state = 1, Lapw%nv2_Tot 
            WRITE (*,*) 'State:', state, ew(state,2) 
            DO n = 1, Lapw%nv2_Tot 
               WRITE (*,*) lapw%KP%k1p(n,jspin),lapw%kp%k2p(n,jspin),   &
     &              REAL(Ev(n,state))                                   
            ENDDO 
            WRITE (*,*) 'Derivative' 
            DO n = 1, Lapw%nv2_Tot 
               WRITE (*,*)lapw%KP%k1p(n,jspin),lapw%kp%k2p(n,jspin),    &
     &              REAL(Ev(n+Lapw%nv2_Tot,state))                      
            ENDDO 
         ENDDO 
         CALL juDFT_error('CBS: 2d-states') 
      ENDIF 
                                                                        
      !>                                                                
                                                                        
      !<-- Debugging: output of statistics of the eigen-states          
                                                                        
#ifdef CPP_DEBUG                                                        
      DO n=1,2*lapw%nv2_tot 
         WRITE(16,*) 'State:',n 
         WRITE(16,*) ew(n,2),curr(n) 
         WRITE(16,*) maxval(abs(ev(:,n))),minval(abs(ev(:,n))) 
         WRITE(16,*) sum(ev(:,n)) 
      ENDDO 
#endif                                                                  
                                                                        
      !>                                                                
                                                                        
      !<--Test the Symmetry of the states                               
                                                                        
       ALLOCATE(syms(size(ev,2),size(mrot,3)+1)) 
       syms=0 
       IF (lapw%nv2(1) == lapw%nv2_tot)                                 &
     &      CALL symcheck(jspin,bk,cell%amat,ev,ew(:,2),curr,mrot,lapw  &
     &      ,syms)                                                      
                                                                        
      !>                                                                
                                                                        
      !<-- output of CBS (eigenvalues)                                  
      CPP_juDFT_timestart("writing of CBS")
      DO n=1,size(ew,1)
        cnorm=dot_product(ev(:lapw%nv2(1),n),ev(:lapw%nv2(1),n))
        logderiv(1,n)=ev(1,n)
        logderiv(2,n)=ev(lapw%nv2_tot+1,n)
        !logderiv(1,n)=dot_product(ev(lapw%nv2(1)+1:lapw%nv2_tot+lapw%nv2(1),n),  &
     !&      ev(:lapw%nv2(1),n))/cnorm
        if (Lapw%nv2_Tot>lapw%nv2(1)) then
           !this is a noco calculation
           logderiv(3,n)=ev(lapw%nv2(1)+1,n)
           logderiv(4,n)=ev(lapw%nv2(1)+lapw%nv2_tot+1,n)
           !cnorm=dot_product(ev(:lapw%nv2_tot,n),ev(:lapw%nv2_tot,n))
           !logderiv(2,n)=dot_product(ev(lapw%nv2_tot+1:2*lapw%nv2_tot,n),       &
     !&      ev(:lapw%nv2_tot,n))/cnorm
        endif
      enddo


      !Generate spin-projection
      if (Lapw%nv2_Tot>lapw%nv2(1)) then
           !this is a noco calculation
           n2=lapw%nv2(1)
           DO n=1,size(ew,1)
                 spinp(0,n)=dot_product(ev(:n2,n),ev(:n2,n))/dot_product(ev(:2*n2,n),ev(:2*n2,n))
                 spinp(1,n)=ev(1,n)*conjg(ev(n2+1,n))+conjg(ev(1,n))*ev(n2+1,n)
                 spinp(2,n)=cmplx(0.,-1.)*(ev(1,n)*conjg(ev(n2+1,n))-conjg(ev(1,n))*ev(n2+1,n))
                 spinp(3,n)=ev(1,n)*conjg(ev(1,n))-ev(n2+1,n)*conjg(ev(n2+1,n))
                 spinp(1:3,n)=spinp(1:3,n)/sqrt(dot_product(spinp(1:3,n),spinp(1:3,n)))
           enddo
      else
          spinp=0.0
      endif
                                                                        
      !set crit. for second complex :: line                             
      IF (gfinp%CBS_bz>0) bz=gfinp%CBS_bz 
      if (.not.gfinp%l_hdfio) THEN
               CALL outCBS(mpi,nk,en,lapw%nv2_tot,ew(:,2),bz,bk,syms,curr,spinp,logderiv)
      else
               CALL outCBS_HDF(nk,en,jspin,ew(:,2),syms,curr,spinp)
      endif
      DEALLOCATE(syms) 
      CPP_juDFT_timestop("writing of CBS")
      !>                                                                
                                                                        
      !<-- check symmetry of Eigenvalues                                
                                                                        
      IF (mpi%pe0) THEN 
         WRITE(6,'(a,e10.5,a,e10.5,a,f10.5)') 'CBS-eigenvalues max:'    &
     &        ,MAXVAL(ABS(ew(:,1))),'min:',MINVAL(ABS(ew(:,1)))         &
     &        ,'Product:',ABS(ew(MINLOC(ABS(ew(:,1))),1)                &
     &        *ew(MAXLOC(ABS(ew(:,1))),1))                              
      ENDIF 
                                                                        
      !>                                                                
      !<-- calculate embedding potential!                               
                                                                        
      CPP_juDFT_timestart("CBS:gen of Sigma")
      INQUIRE(FILE='noembtest',EXIST=l_noembtest) 
      IF(.NOT.l_noembtest)THEN 
       CALL gf_generateEmbPot(en,nk,jspin,mpi%pe0,ev,ew(:,2),T1         &
     &     ,T2,curr,lapw,gfinp,layer)                                   
      ENDIF 
      CPP_juDFT_timestop("CBS:gen of Sigma")
      !>                                                                
                                                                        
      CALL priv_embtest(jspin,nk,en,lapw%nv2_tot,lapw) 
      DEALLOCATE(ew,curr,ev) 
      END SUBROUTINE gf_CBS 
                                                                        
      !>                                                                
                                                                        
      !<-- S:priv_current(mpi,lapw,gfinp,ew,ev,curr)                    
      recursive SUBROUTINE priv_current(en,nk,mpi,lapw,gfinp,ew,ev,curr,eps_non_bloch_in)
!-----------------------------------------------                        
!                                                                       
!           (last modified:09-10-19) D. Wortmann                        
!-----------------------------------------------                        
      USE m_gf_energies,ONLY   : gf_Z 
      USE m_gf_curvy2dprojector 
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)       :: en,nk 
      TYPE(t_mpi),INTENT(IN)   :: mpi 
      TYPE(t_lapw),INTENT(IN ) :: lapw 
      TYPE(t_gfinp),INTENT(IN) :: gfinp 
      COMPLEX,INTENT(IN)       :: ev(:,:) 
      COMPLEX,INTENT(IN)       :: ew(:,:) 
      REAL   ,INTENT(OUT)      :: curr(:) 
      real, intent(in),optional :: eps_non_bloch_in
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: n 
      LOGICAL             :: l_exist 
      real                :: eps_non_bloch
      integer             :: n_evanescent,n_decay,n_bloch,n_blochin
      !>
      if (present(eps_non_bloch_in)) then
         eps_non_bloch=eps_non_bloch_in
      else
         eps_non_bloch=gfinp%eps_non_bloch
      endif
      !<-- calculate the current for each state!                        
#ifdef CPP_CUOVLP                                                       
      IF(.NOT.ALLOCATED(basisoverlaps))THEN 
           CALL gf_curvy2dprojector(cell,lapw,.FALSE.) 
         ENDIF 
         vector = matmul(basisoverlaps(1:lapw%nv2_tot,1:lapw%nv2_tot,1),&
     &      Ev(lapw%nv2_tot+1:2*lapw%nv2_tot,1:2*lapw%nv2_tot))         
         DO n = 1, 2*Lapw%nv2_Tot 
           curr(n) = AIMAG(  DOT_PRODUCT(  Ev(1:Lapw%nv2_Tot,n)         &
     &        ,vector(1:lapw%nv2_tot,n)     ))                          
         ENDDO 
#else                                                                   
         DO n = 1, 2*Lapw%nv2_Tot 
           curr(n) = AIMAG(  DOT_PRODUCT(  Ev(1:Lapw%nv2_Tot,n)         &
     &        ,Ev(Lapw%nv2_Tot+1:2*Lapw%nv2_Tot,n)     ))               
         ENDDO 
#endif                                                                  
      !>                                                                
      !<--check IF the evanescent waves DO not carry current            
                                                                        
      ! & set the sign of the 'current' of these states to indicate the 
      ! direction of decay
      n_decay=0
      n_evanescent=0
      n_blochin=0
      DO n = 1, 2*Lapw%nv2_Tot 
         IF ( ABS(ABS(ew(n,1))-1)>eps_non_bloch ) THEN
            n_evanescent=n_evanescent+1
            IF  (ABS(ew(n,1))>1) n_decay=n_decay+1
            IF (( ABS(curr(n))>gfinp%eps_current .AND.                  &
     &           AIMAG(gf_Z(en,0))<1E-10).AND.mpi%pe0) THEN             
               WRITE(6,*)'WARNING! Evanescent with current detected!' 
               WRITE (6,*) '|Lambda|=', ABS(ew(n,1)) 
               WRITE (6,*) 'curr    =', curr(n) 
            ENDIF 
!     set sign of current to indicate direction of decay                
            IF (ABS(curr(n))>EPSILON(1.0)) THEN 
               curr(n) =- ABS(curr(n))*SIGN(1.0,ABS(ew(n,1))-1) 
            ELSE 
               curr(n) =- SIGN(TINY(1.0),ABS(ew(n,1))-1) 
            ENDIF 
         ELSE !Its a Bloch state!!
            if (curr(n)>0) n_blochin=n_blochin+1
         ENDIF
      ENDDO 
      !>                                                                
                                                                        
      IF (COUNT(curr>0)/=Lapw%nv2_Tot) THEN
         !Decomposition of states failed, first try again
         n_bloch=count(ABS(ABS(ew(:2*Lapw%nv2_Tot,1))-1)<eps_non_bloch)
         if (eps_non_bloch<gfinp%eps_non_bloch.and.eps_non_bloch>5*epsilon(1.0)) then
                write(6,*) "Decomposition of States failed,retry with smaller eps_non_bloch:",eps_non_bloch/2.0
                write(6,*) "Bloch states:",n_bloch," 'in':",n_blochin
                write(6,*) "Evanescent:",n_evanescent," decaying:",n_decay
                call priv_current(en,nk,mpi,lapw,gfinp,ew,ev,curr,eps_non_bloch/2.0)
                return
         endif
         if (eps_non_bloch<gfinp%eps_non_bloch*1000) then
                if (eps_non_bloch<=5*epsilon(1.0)) eps_non_bloch=gfinp%eps_non_bloch
                write(6,*) "Decomposition of States failed,retry with larger eps_non_bloch:",eps_non_bloch*2.0
                write(6,*) "Bloch states:",n_bloch," 'in':",n_blochin
                write(6,*) "Evanescent:",n_evanescent," decaying:",n_decay
                call priv_current(en,nk,mpi,lapw,gfinp,ew,ev,curr,eps_non_bloch*2.0)
                return
         endif
         if (n_blochin*2.ne.n_bloch) THEN
                WRITE(*,*) "CBS failed to decompose into in/out states"
                write(*,*) "Bloch states:",n_bloch," 'in':",n_blochin
                WRITE(*,*) "nk   en   n   EW (Re+Im) Direction  Abs(ew)-1"
                DO n = 1, 2*Lapw%nv2_Tot
                WRITE(*,"(3i5,10(e15.4,1x))") nk,en,n,ew(n,1),curr(n)       &
     &           ,SIGN(1.0,curr(n)),ABS(ew(n,1))-1
                ENDDO
                call juDFT_error("CBS_decomposition failed")
         else
             WRITE(*,*) "CBS failed to decompose into in/out states (evanescent)"
             write(*,*) "evanescent:",n_evanescent," decaying:",n_decay
             WRITE(*,*) "nk   en   n   EW (Re+Im) Direction  Abs(ew)-1"
             DO n = 1, 2*Lapw%nv2_Tot
               WRITE(*,"(3i5,10(e15.4,1x))") nk,en,n,ew(n,1),curr(n)       &
     &           ,SIGN(1.0,curr(n)),ABS(ew(n,1))-1                      
             ENDDO
             call juDFT_error("CBS_decomposition failed")
         endif
         !INQUIRE(FILE ="curr_fix",EXIST = l_exist)
         !IF(.NOT.l_exist) CALL juDFT_error("CBS decomposition failed")
         !n = COUNT(curr<0)-Lapw%nv2_Tot
         !IF (n>0) THEN
         !   curr(:n) =-curr(:n)
         !ELSE
         !   curr(Lapw%nv2_Tot+n:)=-curr(Lapw%nv2_Tot+n:)
         !ENDIF
      ENDIF 
                                                                        
      END SUBROUTINE 
      !>                                                                


      SUBROUTINE outCBS_HDF(nk,en,jspin,ew,syms,curr,spinp)
      !<-- S: outCBS(mpi,nk,en,lapw%nv2_tot,ew,bz,bk,syms)              
      USE m_hdf_tools,ONLY:io_WRITE
      USE m_hdf_accessprp
      IMPLICIT NONE
      !<--Arguments
      INTEGER,INTENT(IN)    :: nk,en,jspin
      INTEGER,INTENT(IN)    :: syms(:,:)
      COMPLEX,INTENT(IN)    :: ew(:)
      REAL   ,INTENT(IN)    :: curr(:),spinp(0:,:)


      CALL io_write(re_k_var,(/1,en,nk,jspin/),(/size(ew),1,1,1/),real(ew))
      CALL io_write(im_k_var,(/1,en,nk,jspin/),(/size(ew),1,1,1/),aimag(ew))

      CALL io_write(curr_var,(/1,en,nk,jspin/),(/size(ew),1,1,1/),curr)
      CALL io_write(spinp_var,(/1,1,en,nk,jspin/),(/4,size(ew),1,1,1/),spinp)
      !CALL io_write(sym_var,(/1,1,en,nk,jspin/),(/size(ew),size(syms,2),1,1,1/),syms)

      END SUBROUTINE

      SUBROUTINE outCBS_openfile(jspin,jspins,kpts,c,nv2d,l_hdfio)
      USE m_gf_types
      USE m_hdf_tools
      USE m_gf_energies
      USE m_gf_math, ONLY  : tpi_const
      USE m_hdf_accessprp
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: jspin,jspins,nv2d
      REAL,INTENT(IN)    :: c
      LOGICAL,INTENT(IN) :: l_hdfio
      TYPE(t_kpts),INTENT(IN):: kpts

      integer :: hdferr
      if (l_hdfio) THEN
        if (jspin>1) return !open the file only once

        CALL h5fcreate_f("CBS.hdf",H5F_ACC_TRUNC_F,fid,hdferr,H5P_DEFAULT_F,hdf_access_prp("CBS.hdf"))
        CALL io_createvar(fid,"Re_k",H5T_NATIVE_DOUBLE,(/2*nv2d,gf_noen(),kpts%nkpts,jspins/),re_k_var)
        CALL io_createvar(fid,"Im_k",H5T_NATIVE_DOUBLE,(/2*nv2d,gf_noen(),kpts%nkpts,jspins/),Im_k_var)
        CALL io_createvar(fid,"Curr",H5T_NATIVE_DOUBLE,(/2*nv2d,gf_noen(),kpts%nkpts,jspins/),curr_var)
        CALL io_createvar(fid,"Spinp",H5T_NATIVE_DOUBLE,(/4,2*nv2d,gf_noen(),kpts%nkpts,jspins/),spinp_var)

        CALL io_write_var(fid,"kpts",kpts%bk(:2,:))
        CALL io_write_var(fid,"energies",real(gf_allz(0)))
        CALL io_write_var(fid,"bz",(/tpi_const/c/))
      else

        OPEN(92,FILE ='CBS_eigenvalues.'//CHAR(jspin+48),FORM      ='formatted')
        OPEN(93,FILE ='CBS_bloch.'//CHAR(jspin+48),FORM      ='formatted')
      endif

      END subroutine

      SUBROUTINE outCBS_closefile(close)
      use m_hdf_tools
      implicit none
      LOGICAL,INTENT(IN) :: close
#ifdef CPP_MPI
      integer:: hdferr
      if (.not.close) return
      CALL io_dclose(re_k_var,hdferr)
      CALL io_dclose(im_k_var,hdferr)
      CALL io_dclose(curr_var,hdferr)
      CALL io_dclose(spinp_var,hdferr)
      CALL io_hdfclose(fid,hdferr)
#else
      CLOSE (92)
      CLOSE(93)
#endif
      END subroutine
      SUBROUTINE outCBS(mpi,nk,en,nv2,ew,bz,bk_in,syms,curr_in,spinp_in,logderiv_in)
!*****************************************************************      
!     This SUBROUTINE writes the CBS into formated files                
!     the DATA is written in the following FORMAT:                      
!     k_z,REAL(z),axis,nkp,Re(kz),Im(kz),ABS(kz)                        
!     WHERE axis indicates on which axis the k-point is located         
!                                                                       
!                       D. Wortmann, Tokyo, 2001                        
!****************************************************************       
      USE m_gf_energies,ONLY  :gf_Z 
      USE m_gf_types,ONLY   :t_mpi 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)    :: nv2,nk,en 
      TYPE(t_mpi),INTENT(IN) :: mpi 
      REAL,INTENT(IN)       ::bz,bk_in(2) 
      INTEGER,INTENT(IN)    ::syms(:,:) 
      COMPLEX,INTENT(IN)    ::ew(:) 
      REAL   ,INTENT(IN)    :: curr_in(:),spinp_in(0:,:)
      COMPLEX,INTENT(IN)    :: logderiv_in(:,:)
      !>                                                                
      !<-- Locals                                                       
      INTEGER            ::n,nkp,nnv2 
      COMPLEX,ALLOCATABLE ::kz(:) 
      INTEGER,ALLOCATABLE::nsyms(:,:) 
      REAL   ,ALLOCATABLE :: curr(:),spinp(:,:)
      COMPLEX,ALLOCATABLE :: logderiv(:,:)
      COMPLEX            ::zz 
      REAL               ::ikz,realz,rkz,akz,eps,bk(2) 
      LOGICAL            :: l_bloch 
      !>                                                                
      !<-- some MPI Stuff                                               
#ifdef CPP_MPI                                                          
#include"cpp_double.h"                                                  
      INCLUDE 'mpif.h' 
      INTEGER::rank 
                                    !MPI error+status                   
      INTEGER::e,s(MPI_STATUS_SIZE) 
#endif                                                                  
      !>                                                                
      !<--copy Data to non-INTENT(IN) local arrays                      
      ALLOCATE(kz(1:2*nv2),nsyms(2*nv2,size(syms,2))) 
      ALLOCATE(curr(1:2*nv2),spinp(0:3,2*nv2),logderiv(4,2*nv2))
      nnv2 = nv2 
      zz=gf_z(en,0) 
      kz(1:2*nv2)=ew 
      nsyms(1:2*nv2,:)=syms(1:2*nv2,:) 
      curr(1:2*nv2) = curr_in(1:2*nv2) 
      spinp=spinp_in(0:,:2*nv2)
      logderiv=logderiv_in(:,:2*nv2)
      nkp=nk 
      bk=bk_in 
      !>                                                                
                                                                        
#ifdef CPP_MPI                                                          
      IF (mpi%irank/=0) THEN 
         !<--Send data to PE=0                                          
         CALL MPI_SEND(nkp, 1, MPI_INTEGER, 0, 1,MPI_COMM_WORLD, E) 
         CALL MPI_SEND(zz,1,CPP_MPI_COMPLEX,0,2,MPI_COMM_WORLD,E) 
         CALL MPI_SEND(nnv2,1, MPI_INTEGER, 0, 4,MPI_COMM_WORLD, E) 
         CALL MPI_SEND(bk,2, CPP_MPI_REAL,                              &
     &        0, 5,MPI_COMM_WORLD,E)                                    
         CALL MPI_SEND(kz, 2*nnv2, CPP_MPI_COMPLEX,                     &
     &        0, 3,MPI_COMM_WORLD, E)                                   
         CALL mpi_send(curr,2*nnv2,CPP_MPI_REAL,                        &
     &        0,7,MPI_COMM_WORLD,E)                                     
         CALL mpi_send(spinp,4*2*nnv2,CPP_MPI_REAL,                        &
     &        0,8,MPI_COMM_WORLD,E)
     CALL mpi_send(logderiv,4*2*nnv2,CPP_MPI_COMPLEX,                        &
     &        0,9,MPI_COMM_WORLD,E)
         CALL MPI_SEND(nsyms,2*nnv2*size(syms,2),MPI_INTEGER, 0, 6,     &
     &        MPI_COMM_WORLD, E)                                        
         !>                                                             
      ELSE 
         DO rank=0,mpi%isize-1 
            !<--On PE0 recieve Data                                     
            IF (rank/=0) THEN 
               CALL MPI_RECV(nkp,1,MPI_INTEGER,rank,1,                  &
     &              MPI_COMM_WORLD,S, E)                                
              CALL MPI_RECV(zz,1,CPP_MPI_COMPLEX,                       &
     &             rank,2, MPI_COMM_WORLD,S, E)                         
              CALL MPI_RECV(nnv2,1,MPI_INTEGER,rank,4,                  &
     &             MPI_COMM_WORLD,S, E)                                 
              CALL MPI_RECV(bk,2,CPP_MPI_REAL,                          &
     &             rank,5, MPI_COMM_WORLD,S, E)                         

              DEALLOCATE(kz,nsyms) 
              ALLOCATE(kz(1:2*nnv2),nsyms(2*nnv2,SIZE(syms,2))) 
              DEALLOCATE(curr,spinp,logderiv)
              ALLOCATE(curr(1:2*nnv2),spinp(0:3,2*nnv2),logderiv(4,2*nnv2))
              CALL MPI_RECV(kz,2*nnv2,CPP_MPI_COMPLEX,                  &
     &             rank,3,MPI_COMM_WORLD,S, E)                          
              CALL MPI_RECV(curr,2*nnv2,CPP_MPI_REAL,                   &
     &             rank,7,MPI_COMM_WORLD,S,E)                           
              CALL MPI_RECV(spinp,4*2*nnv2,CPP_MPI_REAL,                   &
     &             rank,8,MPI_COMM_WORLD,S,E)
              CALL MPI_RECV(logderiv,4*2*nnv2,CPP_MPI_COMPLEX,                   &
     &             rank,9,MPI_COMM_WORLD,S,E)
              CALL MPI_RECV(nsyms,2*nnv2*SIZE(syms,2),MPI_INTEGER,rank, &
     &             6, MPI_COMM_WORLD,S, E)                              
            ENDIF 
            !>                                                          
#endif                                                                  
                                                                        
            !<--OK this is PE=0 so we WRITE the DATA                    
                                                                        
            INQUIRE(UNIT = 93,OPENED = l_bloch) 
                                                                        
            eps = 1.E-3 
            realz=REAL(zz) 
            DO n=1,2*nnv2 
!               rkz=ABS(REAL(kz(n)))                                    
!               ikz=ABS(AIMAG(kz(n)))                                   
!               akz=ABS(kz(n))                                          
               rkz=REAL(kz(n)) 
               ikz=AIMAG(kz(n)) 
               akz=ABS(kz(n)) 
               IF (ABS(rkz)<=bz/2+eps) THEN 
                  IF ( ABS(ikz)<eps ) THEN 
                     !IF (rkz.LE.bz/2)                                  
                     WRITE (92,999)ABS(rkz),realz,1,nkp,rkz,ikz,akz,bk  &
     &                    ,spinp(0:,n),logderiv(:,n),nsyms(n,:)
                     IF (l_bloch) WRITE (93,"(f0.8,1x,i3,2(1x,f0.8))")  &
     &                    realz,nkp,rkz,curr(n)                         
                     ! plot COMPLEX :: bands                            
                  ELSEIF ( ABS(rkz)<eps ) THEN 
                     !along Re=0                                        
                     WRITE (92,999) -ABS(ikz),realz,2,nkp,rkz,ikz,akz,bk&
     &                    ,spinp(0:,n),logderiv(:,n),nsyms(n,:)
                ELSEIF ( ABS(ABS(rkz)-bz/2)<eps ) THEN 
                     !along Re=bz/2                                     
                     WRITE (92,999)bz/2+ABS(ikz),realz,3,nkp,rkz,ikz,akz&
     &                    ,bk,spinp(0:,n),logderiv(:,n),nsyms(n,:)
                  ELSE 
                     !general COMPLEX :: k                              
                     WRITE (92,999) ABS(rkz),realz,0,nkp,rkz,ikz,akz,bk &
     &                    ,spinp(0:,n),logderiv(:,n),nsyms(n,:)
                     WRITE (92,999) -ABS(ikz),realz,0,nkp,rkz,ikz,akz,bk&
     &                    ,spinp(0:,n),logderiv(:,n),nsyms(n,:)
                  ENDIF 
               ENDIF 
                                                                        
               !>                                                       
         ENDDO 
                                                                        
#ifdef CPP_MPI                                                          
      ENDDO 
      ENDIF 
#endif                                                                  
      DEALLOCATE(KZ,nsyms) 
  999 FORMAT(f12.6,f12.6,2i5,9(f12.6,1x),'LD: ',8(e15.8,1x),i3,550(i2,1x))
      END SUBROUTINE outCBS 
                                                                        
      !>                                                                
                                                                        
      !<--S: symcheck                                                   
                                                                        
      SUBROUTINE symcheck(jspin,bk,amat,ev,ew,curr,mrot,lapw,sym) 
!*****************************************************************      
!                                                                       
! Calculate the irreducible representations of the electronic states    
!                                                                       
! Double groups do not work yet...                                      
!                                                                       
! Jussi Enkovaara Jlich 2004                                            
!****************************************************************       
      USE m_gf_types,ONLY:t_lapw 
      USE m_gf_math,  ONLY:mat_inverse 
      USE m_fleur_interface 
                                                                        
      IMPLICIT NONE 
      !<--Argument                                                      
      INTEGER,INTENT(IN)       :: jspin 
      COMPLEX,INTENT(IN)       :: ev(:,:),ew(:) 
      REAL, INTENT(IN)         :: curr(:),bk(2) 
      INTEGER,INTENT(IN)       :: mrot(:,:,:) 
      TYPE(t_lapw),INTENT(IN)  :: lapw 
      INTEGER,INTENT(OUT)      ::sym(:,:) 
                                                                        
      !>                                                                
      !<-- Locals                                                       
                                                                        
      REAL :: bk3(3),amat(3,3) 
      INTEGER::nv2,nop,i,j,k 
      INTEGER :: d 
      REAL    ::kv(2),kvtest(2) 
      REAL,PARAMETER::eps=1E-5 
      INTEGER :: nclass,n,n2,c 
      INTEGER :: mtmpinv(3,3) 
      INTEGER :: mrot_k(3,3,10) 
      COMPLEX, SAVE, ALLOCATABLE :: char_table(:,:) 
      COMPLEX :: c_table(9,9) 
      COMPLEX, ALLOCATABLE :: chars(:,:) 
      COMPLEX, ALLOCATABLE :: csum(:,:,:),overlap(:,:) 
      COMPLEX :: ortho2 
      INTEGER, ALLOCATABLE :: gmap(:,:),irrep(:) 
      LOGICAL, ALLOCATABLE :: symdone(:) 
      REAL :: degthre1,degthre2 
      INTEGER :: n1,ndeg,nstates 
      INTEGER, ALLOCATABLE :: deg(:) 
      CHARACTER(LEN=7) :: grpname
      CHARACTER(LEN=5) :: irrname(9)
                                                                        
      LOGICAL, SAVE :: char_written 
      INTEGER, SAVE :: nclass_old=0,nirr 
      CHARACTER(LEN=5), SAVE :: grpname_old='What' 
                                                                        
      LOGICAL :: soc 
      COMPLEX, ALLOCATABLE :: su(:,:,:) 
                                                                        
      !>                                                                
      nv2=size(ev,1)/2 
      soc = (nv2/=lapw%nv2(1)) 
      nstates = 2*nv2 
      nop=size(mrot,3) 
      sym=0 
                                                                        
      ALLOCATE(deg(nstates)) 
      ALLOCATE(symdone(nstates)) 
                                                                        
      !<--Determine the group of k                                      
                                                                        
      bk3(1:2)=bk 
      bk3(3)=0.12345 
      IF (soc) THEN 
         nv2=nv2/2 
         ALLOCATE(su(2,2,2*nop)) 
         CALL fleur_grp_k(mrot,mrot_k,amat,bk3,nclass,nirr,c_table,     &
     &     grpname,irrname,su)                                          
      ELSE 
         CALL fleur_grp_k(mrot,mrot_k,amat,bk3,nclass,nirr,c_table,     &
     &     grpname,irrname)                                             
      ENDIF 
                                                                        
      !>                                                                
                                                                        
      IF (ALLOCATED(char_table)) THEN 
         IF (SIZE(char_table,1)/=nclass) THEN 
            DEALLOCATE(char_table) 
            ALLOCATE(char_table(nclass,nclass)) 
         ENDIF 
      ELSE 
         ALLOCATE(char_table(nclass,nclass)) 
      ENDIF 
                                                                        
      char_table=c_table(1:nclass,1:nclass) 
                                                                        
      IF ((grpname_old/=grpname).OR.(nclass_old/=nclass)) THEN 
         char_written=.FALSE. 
         grpname_old=grpname 
         nclass_old=nclass 
      ENDIF 
      IF (soc) THEN 
         ALLOCATE(chars(8*nv2,nclass)) 
         ALLOCATE(irrep(8*nv2)) 
      ELSE 
         ALLOCATE(chars(2*nv2,nclass)) 
         ALLOCATE(gmap(nv2,nclass)) 
         ALLOCATE(irrep(2*nv2)) 
         ALLOCATE(csum(nstates,nstates,nclass),overlap(nstates,nstates)) 
      ENDIF 
      irrep=0 
      chars=0.0 
                                                                        
      !<--map the g-vectors related by inv(rot)                         
                                                                        
      gmap=0 
      DO c=1,nclass 
         mtmpinv = mat_inverse(mrot_k(:,:,c)*1.0) 
         iloop: DO i=1,nv2 
            kv(1)=lapw%kp%k1p(i,jspin) 
            kv(2)=lapw%kp%k2p(i,jspin) 
            kv=kv+bk 
            kvtest=matmul(kv,mtmpinv(1:2,1:2))-bk 
            DO j=1,nv2 
               kv(1) = lapw%kp%k1p(j,jspin) 
               kv(2) = lapw%kp%k2p(j,jspin) 
               IF (ABS(kvtest(1)-kv(1))<eps.AND.                     &
     &           ABS(kvtest(2)-kv(2))<eps) THEN                      
                  gmap(i,c)=j 
                  CYCLE iloop 
               ENDIF 
            ENDDO 
            WRITE(24,*) 'Sym. analysis is not done, cannot find         &
     &         rotated kv for' , i,kvtest(1),kvtest(2)                  
            RETURN 
         ENDDO iloop 
      ENDDO 
                                                                        
      !>                                                                
                                                                        
                                                                        
      !<--calculate the characters and irreducible representations      
                                                                        
      symdone=.FALSE. 
      degthre1=0.005 
      degthre2=0.005 
                                                                        
      stateloop: DO i=1,nstates 
         IF (symdone(i)) CYCLE stateloop 
         ndeg=0 
         deg=0 
! Find degenerate states                                                
         DO n=1,nstates 
            IF (                                                        &
     &  ((ABS(REAL(ew(i))-REAL(ew(n)))<ABS(REAL(ew(i)))*degthre1).OR.&
     &   (ABS(REAL(ew(i))-REAL(ew(n)))<degthre2)).AND.               &
     &           ((ABS(AIMAG(ew(i))-AIMAG(ew(n)))<ABS(AIMAG(ew(i)))     &
     &           *degthre1).OR.(ABS(AIMAG(ew(i))-AIMAG(ew(n)))<degthre2)&
     &           )) THEN                                                
               ndeg=ndeg+1 
               deg(ndeg)=n 
            ENDIF 
         ENDDO 
! Calculate the representation matrices. In general the states are not o
! so we need also the overlap                                           
         csum=0.0 
         DO c=1,nclass 
            overlap=0.0 
            DO n1=1,ndeg 
               DO n2=1,ndeg 
                  DO k=1,nv2 
                     IF (soc) THEN 
                       csum(n1,n2,c)=csum(n1,n2,c)+conjg(ev(k,deg(n1)))*&
     &                       (su(1,1,c)*ev(gmap(k,c),deg(n2))+          &
     &                       su(1,2,c)*ev(gmap(k,c)+nv2,deg(n2)))+      &
     &                       conjg(ev(k+nv2,deg(n1)))*                  &
     &                       (su(2,1,c)*ev(gmap(k,c),deg(n2))+          &
     &                       su(2,2,c)*ev(gmap(k,c)+nv2,deg(n2)))       
                       overlap(n1,n2)=overlap(n1,n2)+                   &
     &                       conjg(ev(k,deg(n1)))*ev(k,deg(n2))+        &
     &                       conjg(ev(k+nv2,deg(n1)))*ev(k+nv2,deg(n2)) 
                    ELSE 
                       csum(n1,n2,c)=csum(n1,n2,c)+conjg(ev(k,deg(n1)))*&
     &                      ev(gmap(k,c),deg(n2))                       
                       overlap(n1,n2)=overlap(n1,n2)+                   &
     &                      conjg(ev(k,deg(n1)))*ev(k,deg(n2))          
                    ENDIF 
                 ENDDO 
              ENDDO 
           ENDDO 
           overlap(1:ndeg,1:ndeg)=                                      &
     &          mat_inverse(overlap(1:ndeg,1:ndeg))                     
           csum(1:ndeg,1:ndeg,c)=matmul(overlap(1:ndeg,1:ndeg),         &
     &          csum(1:ndeg,1:ndeg,c))                                  
        ENDDO 
                                                                        
! We might have taken degenerate states which are not degenerate due sym
        DO n1=1,ndeg 
           chars(deg(n1),:)=0.0 
           DO n2=1,ndeg 
              IF (ANY(ABS(csum(n1,n2,:))>0.05)) THEN 
                 chars(deg(n1),:)=chars(deg(n1),:)+csum(n2,n2,:) 
              ENDIF 
           ENDDO 
           sym(deg(n1),2:nclass+1)=NINT(REAL(chars(deg(n1),:))) 
           symdone(deg(n1))=.TRUE. 
        ENDDO 
                                                                        
! determine the irreducible presentation                                
         irrloop: DO n1=1,ndeg 
         DO c=1,nclass 
            IF (ALL(ABS(chars(deg(n1),1:nclass)-                        &
     &           char_table(c,1:nclass))<0.001)) THEN                   
               irrep(deg(n1))=c 
               CYCLE irrloop 
! The group was unknown                                                 
            ELSE IF (ALL(ABS(char_table(c,1:nclass))<0.001)) THEN 
               char_table(c,:)=chars(deg(n1),:) 
               irrep(deg(n1))=c 
               CYCLE irrloop 
            ENDIF 
         ENDDO 
      ENDDO irrloop 
                                                                        
      ENDDO stateloop 
                                                                        
      !>                                                                
                                                                        
      IF (.NOT.char_written) THEN 
         WRITE(24,122) bk 
         WRITE(24,*) 'Group is ' ,grpname 
         DO c=1,nclass 
           WRITE(24,123) c,irrname(c),(REAL(char_table(c,n)),n=1,nclass) 
         ENDDO 
         char_written=.TRUE. 
      ENDIF 
  122 FORMAT('Character :: table for k ',2f8.4) 
  123 FORMAT(i3,1x,a5,1x,20f7.3) 
                                                                        
      sym(:,1)=irrep 
                                                                        
      !<--Some tests about the orthogonality of states                  
                                                                        
                                                                        
      WRITE(24,124) 
      WRITE(24,125) 
      DO i=1,nstates 
         DO j=1,i 
            ortho2=0.0 
            DO k=1,2*nv2 
               ortho2=ortho2+conjg(ev(k,i))*ev(k,j) 
              ENDDO 
              IF (ABS(ortho2)>0.001) THEN 
               IF ((irrep(i)/=irrep(j)).AND.irrep(i)>0.AND.irrep(j)  &
     &                >0) THEN                                       
                  WRITE(24,126) i,j,irrep(i),irrep(j),ortho2 
               ENDIF 
              ENDIF 
           ENDDO 
        ENDDO 
                                                                        
  124   FORMAT('Orthogonality of states belonging to different reps.') 
  125   FORMAT('  i',1x,'  j',2x,'irr(i)',2x,'irr(j)',9x,'<i|j>') 
  126   FORMAT(i3,1x,i3,2x,i3,5x,i3,5x,2f9.5) 
                                                                        
        !>                                                              
                                                                        
        DEALLOCATE(csum,overlap) 
        DEALLOCATE(chars) 
        DEALLOCATE(gmap) 
        DEALLOCATE(irrep) 
        DEALLOCATE(deg) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S: priv_embtest                                              
                                                                        
      SUBROUTINE priv_embtest(jspin,nk,en,nv2,lapw) 
!******************************************                             
!    Different tests on the embedding potential....                     
!     Test J.Inglesfield's Idea to calculate eigenvectors               
!     of the Im part of the Embedding potential                         
!     D. Wortmann                                                       
!******************************************                             
      USE m_gf_io2dmat 
      USE m_gf_math 
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN) :: jspin,en,nk,nv2 
      TYPE(t_lapw),INTENT(IN) :: lapw 
      !>                                                                
      !<-- Locals                                                       
                                                                        
      COMPLEX :: mat(nv2,nv2) 
      COMPLEX :: ev(nv2,nv2),ew(nv2) 
      INTEGER :: n 
      LOGICAL :: l_exist 
      !>                                                                
      !Test J.Inglesfield's channel functions....                       
      INQUIRE(FILE ="jingle",EXIST= l_exist) 
      IF (l_exist) THEN 
         IF(.NOT.gf_read2dmat(IO2D_EMB,2,2,en,nk,jspin,lapw,mat)) RETURN 
         mat = AIMAG(mat) 
         CALL eigenvalues(mat,ew,ev) 
         OPEN(99,FILE ="jingle",STATUS ="old",POSITION ="append",FORM="formatted")
         DO n = 1,nv2 
            IF (ABS(ew(n))>1E-3) THEN 
               WRITE(99,'(3(i4,1x),2(f10.4,1x))') jspin,nk,en,ew(n) 
            ENDIF 
         ENDDO 
         CLOSE(99) 
      ENDIF 
      !Calculate Sigma*Sigma^*^-1                                       
      INQUIRE(FILE ="SigmaSigma",EXIST = l_exist) 
      IF (l_exist) THEN 
         IF(.NOT.gf_read2dmat(IO2D_EMB,2,2,en,nk,jspin,lapw,mat)) RETURN 
         mat = MATMUL(CONJG(mat),mat_inverse(mat)) 
         CALL eigenvalues(mat,ew,ev) 
         OPEN(99,FILE ="SigmaSigma",STATUS ="old",POSITION ="append"    &
     &        ,FORM   ="formatted")                                     
         DO n = 1,nv2 
            WRITE(99,'(3(i4,1x),2(f10.4,1x))') jspin,nk,en,ew(n) 
         ENDDO 
         CLOSE(99) 
      ENDIF 
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      END                                           
