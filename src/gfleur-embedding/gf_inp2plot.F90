!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_inp2plot 
      use m_juDFT 
      USE m_gf_types 
      USE m_gf_embdesc 
      IMPLICIT NONE
      PRIVATE 
                                        !Max no of atoms in plot        
      INTEGER,PARAMETER :: MAXAT = 1000 
      CHARACTER*2       :: namat(0:103) 
      DATA namat/'VA',' H','HE','LI','BE',' B',' C',' N',' O',' F','NE',&
     &     'NA','MG','AL','SI',' P',' S','CL','AR',' K','CA','SC','TI', &
     &     ' V','CR','MN','FE','CO','NI','CU','ZN','GA','GE','AS','SE', &
     &     'BR','KR','RB','SR',' Y','ZR','NB','MO','TC','RU','RH','PD', &
     &     'AG','CD','IN','SN','SB','TE',' J','XE','CS','BA','LA','CE', &
     &     'PR','ND','PM','SM','EU','GD','TB','DY','HO','ER','TM','YB', &
     &     'LU','HF','TA',' W','RE','OS','IR','PT','AU','HG','TL','PB', &
     &     'BI','PO','AT','RN','FR','RA','AC','TH','PA',' U','NP','PU', &
     &     'AM','CM','BK','CF','ES','FM','MD','NO','LW'/                
      PUBLIC gf_inp2plot 
      CONTAINS 
      !<-- S: gf_inp2plot                                               
      SUBROUTINE gf_inp2plot(atoms,cell,gfinp) 
!-----------------------------------------------                        
!     Use atomic positions to generate imput files for                  
!      visualisation tools                                              
!           (last modified: 05-01-10) D. Wortmann                       
!-----------------------------------------------                        
                                                                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_atoms),INTENT(IN)  :: atoms 
      TYPE(t_cell),INTENT(IN)   :: cell 
      TYPE(t_embinp),INTENT(IN)  :: gfinp 
      !>                                                                
      !<--Locals                                                        
                                                                        
      TYPE(t_embdesc)           :: descR,descL 
      INTEGER                   :: nx,ny,nl,nr 
      REAL                      :: atpos(3,MAXAT),at_rmt(MAXAT) 
      INTEGER                   :: at_z(MAXAT),at_type(MAXAT),natoms 
                                                                        
      !>                                                                
                                        !Do nothing if not asked for    
      IF (.NOT.gfinp%l_inp2plot) RETURN 
                                      !output                           
                                                                        
               !Interactive IO starts now                               
      CLOSE(6) 
      WRITE (*,*) "Inp2Plot routine" 
      WRITE (*,*) "Specify the 2D-Unit cell" 
      WRITE (*,*) "E.g. 1 1 for p(1x1) or 3 4 for p(3x4)" 
      READ (*,*) nx,ny 
                                                                        
      !<-- Read the embedding-plane descriptors                         
      CALL gf_readdescriptor(1,1,cell%amat,descL) 
      CALL gf_readdescriptor(1,2,cell%amat,descR) 
                                                                        
      IF (descL%valid) THEN 
         WRITE(*,*) "How many principle layers of the "                 &
     &        //"left electrode should be added?"                       
         READ(*,*) nl 
      ELSE 
         nl = 0 
      ENDIF 
      IF (descL%valid) THEN 
         WRITE(*,*) "How many principle layers of the "                 &
     &        //"right electrode should be added?"                      
         READ(*,*) nr 
      ELSE 
         nr = 0 
      ENDIF 
      !>                                                                
                                                                        
      !<-- Generate all atomic positions                                
      CALL priv_makeatomspos(atoms,cell,descL,descR,nx,ny,nl,nr         &
     &     ,atpos,at_z,at_TYPE,at_rmt,natoms)                           
      WRITE(*,*) "No of atoms found:",natoms 
      !>                                                                
                                                                        
      !<-- Write Output files                                           
      CALL priv_xyz(atpos(:,:natoms),at_z(:natoms)) 
      CALL priv_pdb(atpos(:,:natoms),at_z(:natoms)) 
      CALL priv_pov(atpos(:,:natoms),at_z(:natoms),at_rmt(:natoms)) 
      !>                                                                
      CALL juDFT_error("gf_inp2plot done") 
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S: priv_xyz(atpos,at_z)                                      
      SUBROUTINE priv_xyz(atpos,at_z) 
!-----------------------------------------------                        
!   write a xyz-file                                                    
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      REAL   ,INTENT(IN)     :: atpos(:,:) 
      INTEGER,INTENT(IN)     :: at_z(:) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: n 
                                                                        
      !>                                                                
      OPEN(99,FILE ='inp.xyz',FORM='formatted') 
      WRITE(99,'(i3)') size(at_z) 
      WRITE(99,*) 'Crystall structure from GF-FLEUR' 
      DO n = 1,SIZE(at_z) 
         WRITE(99,'(a4,3f12.7)') namat(at_z(n)),atpos(:,n) 
      ENDDO 
      CLOSE(99) 
      END SUBROUTINE 
      !>                                                                
      !<-- S: priv_pdb(atpos,at_z)                                      
      SUBROUTINE priv_pdb(atpos,at_z) 
!-----------------------------------------------                        
!   write a pdb-file                                                    
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      REAL   ,INTENT(IN)     :: atpos(:,:) 
      INTEGER,INTENT(IN)     :: at_z(:) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: n 
      !>                                                                
                                                                        
      OPEN(99,FILE ='inp.pdb',FORM='formatted') 
      WRITE(99,'(i3)') size(at_z) 
      WRITE(99,'(a6,2x,a60)') "TITLE ",'Crystall structure from GFFLEUR' 
      DO n = 1,SIZE(at_z) 
         WRITE(99,'(a4,3f12.7)') namat(at_z(n)),atpos(:,n) 
         WRITE(99,FMT=98) "ATOM  ",n,namat(at_z(n))//'      ',' '       &
     &                    ,'   ','A',0,' ',atpos(:,n),0.0,0.0,          &
     &                    ' ',namat(at_z(n))                            
      ENDDO 
      CLOSE(99) 
   98 FORMAT(a6,i5,a4,a1,a3,a1,i4,a1,f8.3,f8.3,f8.3,f6.2,f6.2,a4,a2) 
      END SUBROUTINE 
      !>                                                                
      !<-- S: priv_pov(atpos,at_z,at_rmt)                               
      SUBROUTINE priv_pov(atpos,at_z,at_rmt) 
!-----------------------------------------------                        
!      write a povray file                                              
!           (last modified: 2004-00-00) D. Wortmann                     
!           (last modified: 2006-09-00) A. Hanuschkin                   
!-----------------------------------------------                        
      USE m_gf_povrayColors 
      IMPLICIT NONE 
      !<--Arguments                                                     
      REAL   ,INTENT(IN)     :: atpos(:,:),at_rmt(:) 
      INTEGER,INTENT(IN)     :: at_z(:) 
      !>                                                                
      !<-- Locals                                                       
      REAL                :: col(3,103) 
      INTEGER             :: nn 
      !>                                                                
      col=getColors() 
                                                                        
      !<-- WRITE HEADER                                                 
      OPEN (99,FILE ='struct.pov',FORM='formatted',STATUS='unknown') 
      WRITE (99,*) '#include "colors.inc"' 
      WRITE (99,*) '#include "shapes.inc"' 
      WRITE (99,*) 'global_settings { max_trace_level 20 ' 
      WRITE (99,*) '                  assumed_gamma 2.2 }' 
      WRITE (99,*) 'light_source { <60,-10,-10>' 
      WRITE (99,*) '              color rgb <2.5,2.5,2.5> }' 
      WRITE (99,*) 'camera { location <90,10,10>' 
      WRITE (99,*) '         look_at <0.0,0.0,0.0> angle 20 }' 
      WRITE (99,*) 'background {color White}' 
      !Colors                                                           
      WRITE (99,*) ' #declare R1 =  pigment{ color Black } ' 
      WRITE (99,*) ' #declare Rd =  0.05 ;' 
      DO nn = 1,size(at_z) 
        IF (nn<10) THEN 
          WRITE (99,1005) nn,col(:,at_z(nn)+1) 
        ELSEIF ((nn>9).AND.(nn<100)) THEN 
          WRITE (99,1006) nn,col(:,at_z(nn)+1) 
        ELSE 
          WRITE (99,1007) nn,col(:,at_z(nn)+1) 
        ENDIF 
      ENDDO 
                                                                        
      DO nn = size(at_z)+1,103 
        IF (nn<10) THEN 
          WRITE (99,1005) nn,col(:,size(at_z)+1) 
        ELSEIF ((nn>9).AND.(nn<100)) THEN 
          WRITE (99,1006) nn,col(:,size(at_z)+1) 
        ELSE 
          WRITE (99,1007) nn,col(:,size(at_z)+1) 
        ENDIF 
      ENDDO 
                                                                        
                                                                        
 1005 FORMAT ('#declare Acol',i1.0,'= color rgb <',3f4.1,'>;') 
 1006 FORMAT ('#declare Acol',i2.0,'= color rgb <',3f4.1,'>;') 
 1007 FORMAT ('#declare Acol',i3.0,'= color rgb <',3f4.1,'>;') 
      !>                                                                
      !<--Assign colors to atoms                                        
      DO nn = 1,9 
         WRITE (99,1010) nn,nn 
      ENDDO 
 1010 FORMAT ('#declare Atom',i1.0,' = pigment { Acol',i1,' }') 
      DO nn = 10,99 
         WRITE (99,1011) nn,nn 
      ENDDO 
 1011 FORMAT ('#declare Atom',i2.0,' = pigment { Acol',i2,' }') 
      DO nn = 100,103 
         WRITE (99,1012) nn,nn 
      ENDDO 
 1012 FORMAT ('#declare Atom',i3.0,' = pigment { Acol',i3,' }') 
                                                                        
                                                                        
      WRITE (99,1015) 
 1015 FORMAT (/,'#declare Ascale = 0.80 ;') 
                                                                        
      DO nn = 1,size(at_rmt) 
        IF (nn<10) THEN 
         WRITE (99,1030) nn,at_rmt(nn) 
        ELSEIF ((nn>9).AND.(nn<100)) THEN 
         WRITE (99,1031) nn,at_rmt(nn) 
        ELSE 
         WRITE (99,1032) nn,at_rmt(nn) 
        ENDIF 
      ENDDO 
 1030 FORMAT ('#declare Asize',i1.0,' = ',f6.3,'*Ascale ;') 
 1031 FORMAT ('#declare Asize',i2.0,' = ',f6.3,'*Ascale ;') 
 1032 FORMAT ('#declare Asize',i3.0,' = ',f6.3,'*Ascale ;') 
                                                                        
                                                                        
                                                                        
                                                                        
      !>                                                                
      !<-- output the atomic positions                                  
                                                                        
      DO nn = 1, size(atpos,2) 
        IF (nn<10) THEN 
         WRITE (99,1020)                                                &
     &           atpos(:,nn),nn,nn,nn,at_z(nn)                          
        ELSEIF ((nn>9).AND.(nn<100)) THEN 
         WRITE (99,1021)                                                &
     &           atpos(:,nn),nn,nn,nn,at_z(nn)                          
        ELSE 
         WRITE (99,1022)                                                &
     &           atpos(:,nn),nn,nn,nn,at_z(nn)                          
        ENDIF 
      ENDDO 
 1020 FORMAT('sphere { <',3(f8.3,','),'>, Asize',i1.0,' texture { Atom',&
     &       i1.0,' } } // ',2i4)                                       
 1021 FORMAT('sphere { <',3(f8.3,','),'>, Asize',i2.0,' texture { Atom',&
     &       i2.0,' } } // ',2i4)                                       
 1022 FORMAT('sphere { <',3(f8.3,','),'>, Asize',i3.0,' texture { Atom',&
     &       i2.0,' } } // ',2i4)                                       
                                                                        
                                                                        
      !>                                                                
      CLOSE (99) 
      END SUBROUTINE 
      !>                                                                
      !<-- S: priv_makeatomspos(atoms,cell,descL,descR,nx,ny,nl,nr,atpos
      SUBROUTINE priv_makeatomspos(atoms,cell,descL,descR,nx,ny,nl      &
     &     ,nr,atpos,at_z,at_TYPE,at_rmt,natoms)                        
!-----------------------------------------------                        
!     generates a list of atomic positions                              
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_shiftatoms 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_atoms),INTENT(IN)    :: atoms 
      TYPE(t_cell),INTENT(IN)     :: cell 
      TYPE(t_embdesc),INTENT(IN)  :: descL,descR 
      INTEGER,INTENT(IN)          :: nx,ny,nl,nr 
      REAL   ,INTENT(OUT)         :: atpos(3,MAXAT),at_rmt(MAXAT) 
      INTEGER,INTENT(OUT)         :: at_z(MAXAT),at_type(MAXAT) 
      INTEGER                     :: natoms 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: na,n,nn 
                                                                        
      !>                                                                
      at_type=0 
      !<-- First get Atoms within unit-cell                             
      natoms=1 
      na=1 
      DO n=1,atoms%ntype 
         DO nn=1,atoms%neq(n) 
            atpos(:,natoms) = atoms%pos(:,na) 
            at_z(natoms) =int(atoms%zatom(n)) 
            at_rmt(natoms)=atoms%rmt(n) 
            natoms=natoms+1 
            IF (natoms>maxat) CALL juDFT_error("natoms>maxat in inp2plot") 
            na=na+1 
         ENDDO 
      ENDDO 
      at_type(:natoms)=1 
      !>                                                                
      !<-- Add atoms from left side                                     
                                                                        
      IF (descL%valid.AND.nl>0) THEN 
         DO n = 1,nl 
            DO nn = 1,descL%all_atoms 
               atpos(:,natoms) = (/descL%dvec(1),descL%dvec(2),cell%z1  &
     &              -descL%dist_aux/)+descL%atoms_pos(:,nn,3)+(n-1)     &
     &              *descL%dvec                                         
               atpos(3,natoms) = -1.* atpos(3,natoms) 
    !           atpos(:,natoms) = (/descL%dvec(1),descL%dvec(2),-cell%z1
    ! +              +descL%dist_aux/)-descL%atoms_pos(:,nn,3)-(n-1)    
    ! +              *descL%dvec                                        
               at_z(natoms)    = descL%atoms_z(nn,3) 
               at_rmt(natoms)  = descL%atoms_rmt(nn,3) 
               at_type(natoms)  = 2 
               natoms=natoms+1 
               IF (natoms>maxat) CALL juDFT_error                         &
     &              ("natoms>maxat in inp2plot")                        
            ENDDO 
         ENDDO 
      ENDIF 
                                                                        
      !>                                                                
      !<--Add atoms from right side                                     
      IF (descr%valid.AND.nr>0) THEN 
         DO n = 1,nr 
            DO nn = 1,descr%all_atoms 
               atpos(:,natoms) = (/descR%dvec(1),descR%dvec(2),cell%z1  &
     &              -descR%dist_aux/)+descR%atoms_pos(:,nn,3)+(n-1)     &
     &              *descR%dvec                                         
               at_z(natoms)    = descR%atoms_z(nn,3) 
               at_rmt(natoms)  = descR%atoms_rmt(nn,3) 
               at_type(natoms)  = 3 
               natoms=natoms+1 
               IF (natoms>maxat) CALL juDFT_error                         &
     &              ("natoms>maxat in inp2plot")                        
            ENDDO 
         ENDDO 
      ENDIF 
      !>                                                                
      natoms=natoms-1 
      !<-- move all atoms to positive (x,y)-Coordinates                 
      CALL gf_shiftatoms(atpos(:2,:natoms),cell%amat(1:2,1:2)) 
      !>                                                                
                                                                        
      !<-- Write atom info to STDOUT                                    
      WRITE(*,*) "Barrier-unit cell from",-1.*cell%z1+descL%dist_aux    &
     &     ," to ",cell%z1-descR%dist_aux                               
      WRITE(*,*) "Aux volumes:",descL%dist_aux,descR%dist_aux 
      WRITE(*,*) "Atoms in plot (in single unit cell):" 
      DO n = 1,natoms 
         WRITE(*,'(i3,1x,a2,1x,3(f10.5,1x),i2)') n,namat(at_z(n))       &
     &        ,atpos(:,n),at_TYPE(n)                                    
      ENDDO 
      !>                                                                
      !<-- Add more atoms for in-plane supercell                        
      IF (nx*ny*natoms>maxat) CALL juDFT_error("natoms>maxat in inp2plot") 
      DO n = 2,nx 
         DO nn = 1,natoms 
            atpos(:2,natoms*(n-1)+nn) = atpos(:2,nn)                    &
     &           +cell%amat(1:2,1)*(n-1)                                
         ENDDO 
         atpos(3,natoms*(n-1)+1:natoms*n) = atpos(3,:natoms) 
         at_z(natoms*(n-1)+1:natoms*n) = at_z(:natoms) 
         at_type(natoms*(n-1)+1:natoms*n) = at_type(:natoms) 
         at_rmt(natoms*(n-1)+1:natoms*n) = at_rmt(:natoms) 
      ENDDO 
      natoms = natoms*nx 
      DO n = 2,ny 
         DO nn = 1,natoms 
            atpos(:2,natoms*(n-1)+nn) = atpos(:2,nn)                    &
     &           +cell%amat(1:2,2)*(n-1)                                
         ENDDO 
         atpos(3,natoms*(n-1)+1:natoms*n) = atpos(3,:natoms) 
         at_type(natoms*(n-1)+1:natoms*n) = at_type(:natoms) 
         at_z(natoms*(n-1)+1:natoms*n) = at_z(:natoms) 
         at_rmt(natoms*(n-1)+1:natoms*n) = at_rmt(:natoms) 
      ENDDO 
      natoms = natoms*ny 
                                                                        
      !>                                                                
      END SUBROUTINE 
      !>                                                                
      END                                           
