      MODULE m_vdWfleur
!
! implements Grimmes D2 and D3 based on routines from V. Caciuc (`19)
!
      CONTAINS

      SUBROUTINE vdW_fleur(
     >                     ntype,natd,neq,film,amat,bmat,pos,zatom,
     >                     nop,ngopr,invtab,mrot,
     <                     e_vdW,f_vdW)

      USE m_types, ONLY: sc_data,atom_data
      USE DFT_D2,  ONLY: driver_DFT_D2
      USE DFT_D3,  ONLY: driver_DFT_D3
      USE m_cotra, ONLY: cotra1,cotra0
      USE m_genernewcell

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: ntype,natd,nop
      LOGICAL, INTENT (IN) :: film
      INTEGER, INTENT (IN) :: neq(ntype),ngopr(natd)
      INTEGER, INTENT (IN) :: mrot(3,3,nop),invtab(nop)
      REAL,    INTENT (IN) :: zatom(ntype),pos(3,natd)
      REAL,    INTENT (IN) :: amat(3,3), bmat(3,3)
      REAL,    INTENT (OUT) :: e_vdW,f_vdW(3,ntype)

      INTEGER  :: NrAtomType,nsize,max_cyc,iop
      INTEGER  :: NrAtoms,i,j,i1,i2,na,ifxyz,irep(3),nat_big
      CHARACTER(LEN=60) :: pfxyz
      REAL :: start,finish,toler,delta
      REAL :: test(3,8),brmin(3),brmax(3),force_i(3),f_rot(3)
      LOGICAL l_D2,l_in_au,l_exist
      TYPE(  sc_data) :: sc_old
      TYPE(atom_data) :: atom,atom_new
      REAL, ALLOCATABLE :: ener(:),force_vdW(:,:),force_max(:,:)

      max_cyc = 20              ! max. tries for bigger supercells
      toler = 0.0005            ! energy tolerance required (eV)
      ALLOCATE ( ener(max_cyc) )

      l_D2 = .false.
      l_in_au = .true. ! .false.
      ifxyz = 869
      pfxyz = "SuperCell.dat"

      sc_old%alat = 1.0
      sc_old%bravais(:,:) = amat(:,:) !  transpose(amat(:,:))
!      sc_old%bravais(:,:) = transpose(amat(:,:))
      NrAtomType=ntype
      ALLOCATE( atom%nr_atom_type(NrAtomType) )                ! neq(ntype)
      atom%nr_atom_type(:) = neq(:)
      NrAtoms = natd
      ALLOCATE( atom%atomic_number(NrAtoms) )                  ! like zatom(natd) not (ntype)
      ALLOCATE( atom%coord_bravais(3,NrAtoms) )                ! taual(3,natd)
      ALLOCATE( atom%coord_cart(3,NrAtoms) )                   ! pos(3,natd)
      ALLOCATE( force_vdW(3,NrAtoms),force_max(NrAtoms,max_cyc) )

      na = 0
      DO i1=1,ntype
        DO i2 = 1,neq(i1)
          na = na + 1
          atom%atomic_number(na) = NINT(zatom(i1))
          CALL cotra1(pos(:,na),atom%coord_bravais(:,na),bmat)
          atom%coord_cart(:,na)= pos(:,na)
        ENDDO
      ENDDO

      l_exist=.false.
      INQUIRE (file=pfxyz,exist=l_exist)       ! check, if convergence test was done
       
      IF (l_exist) THEN
        OPEN (ifxyz,file=pfxyz)  
        READ (ifxyz,'(4i10)',END=99,ERR=99) nat_big,irep(:)
        WRITE(6,*) 'Taking info from SuperCell.dat',irep(:)
        max_cyc = 1                ! just calculate once for parameters found in file
 99     CLOSE (ifxyz)
      ENDIF

      IF (max_cyc > 1) THEN ! determine a supercell size where the vdW energy converges

        test(:,1) = 0.0
        test(:,2) = sc_old%bravais(1,:) ! tr!
        test(:,3) = sc_old%bravais(2,:)
        test(:,4) = sc_old%bravais(3,:)
        test(:,5) = test(:,2) + test(:,3)
        test(:,6) = test(:,3) + test(:,4)
        test(:,7) = test(:,4) + test(:,2)
        test(:,8) = test(:,5) + test(:,4)
      
        brmin(:) = minval(test(:,1:8),2) 
        brmax(:) = maxval(test(:,1:8),2) 
     
        irep(:) = 45.0 / (brmax(:)-brmin(:))
        IF (film) irep(3) = 0
        
      ENDIF

      cyc: DO nsize = 1, max_cyc
!      print'(A)',"#-----------------------------------------------#"
!      print'(A)',"#--------- generate supercell for vdW  ---------#"
!      print'(A)',"#-----------------------------------------------#"
      CALL gener_new_cell(atom,sc_old,irep,ifxyz,pfxyz,
     <                    atom_new)
!
      IF (l_D2) THEN
        print'(A)',"#-----------------------------------------------#"
        print'(A)',"#------------- using DFT-D2 method -------------#"
        print'(A)',"#-----------------------------------------------#"
        CALL cpu_time(start)
        CALL driver_DFT_D2(    atom%coord_cart,    atom%atomic_number, 
     >                     atom_new%coord_cart,atom_new%atomic_number)
        CALL cpu_time(finish)
        PRINT '("Time for DFT_D2 = ",f12.3," seconds.")',finish-start
      ENDIF
!
! note that in the driver_DFT_D3 atom%coord_cart and atom_new%coord_cart
!   are modified from Angstrom to au [they have INTENT(INOUT)]
!
!      print'(A)',"#-----------------------------------------------#"
!      print'(A)',"#------------- using DFT-D3 method -------------#"
!      print'(A)',"#-----------------------------------------------#"
      CALL cpu_time(start)
      CALL driver_DFT_D3(    atom%coord_cart,    atom%atomic_number, 
     >                   atom_new%coord_cart,atom_new%atomic_number,
     >                   l_in_au,
     <                   ener(nsize),force_vdW)

      force_max(:,nsize) = sqrt(force_vdW(1,:)**2 + force_vdW(2,:)**2
     +                                            + force_vdW(3,:)**2)

      IF (max_cyc == 1) THEN ! data taken from file
        delta = 0.0
        EXIT cyc
      ELSE
        IF (nsize > 1) THEN
          delta = ener(nsize)-ener(nsize-1)
          WRITE(6,*) 'Delta = ',delta,' eV'
          DO i1=1,NrAtoms
            WRITE (6,'(a15,i4,a3,3f15.9)') 'Delta vdW force',i1,' : ',
     &                       force_max(i1,nsize)-force_max(i1,nsize-1)
          ENDDO
          IF (abs(delta) < toler) EXIT cyc
        ENDIF
      ENDIF
 
      CALL cpu_time(finish)
!      print '("Time for DFT_D3 = ",f12.3," seconds.")',finish-start

        irep(1:2) = irep(1:2) + 1
        IF (.not.film) irep(3) =  irep(3) + 1
      ENDDO cyc

      IF (abs(delta) > toler) THEN
        WRITE (6,*) 'vdW did not converge with cell size!'
        e_vdW = 0.0
      ELSE
        CALL cpu_time(finish)
        WRITE (6,'(a13,3i5,a4,f12.3,a5)') 'vdW converged',irep(:),
     &                                 ' in ',finish-start,' sec.'
        e_vdW = ener(nsize)/27.21138386
!        WRITE (*,*) 'vdW energy :',e_vdW,' htr'
        WRITE ( 6,8060) e_vdW
        WRITE (16,8060) e_vdW
        DO i1=1,NrAtoms
          WRITE (6,'(a15,i4,a3,6f15.9)') 'vdW force atom ',i1,' : ',
     &                         force_vdW(:,i1),atom%coord_cart(:,i1)
        ENDDO
      ENDIF
 8060 FORMAT (/,/,' ----> vdW (D3)  energy=',t40,f20.10,' htr')

      na = 0
      DO i1=1,NrAtomType
          f_rot(:) = 0.0
          DO i2 = 1,neq(i1)   ! here symmetrization should be done
            na = na + 1
            iop = ngopr(na) ! invtab(ngopr(na))
            CALL cotra1(force_vdW(:,na),force_i(:),bmat)
            DO i = 1,3
              DO j = 1,3
                f_rot(i) = f_rot(i) + mrot(i,j,iop) * force_i(j)
              ENDDO
            ENDDO
          ENDDO
          CALL cotra0(f_rot(:),f_vdW(:,i1),amat)
          f_vdW(:,i1) = f_vdW(:,i1)/neq(i1)
!          write(*,*) 's f:',f_vdW(:,i1)
      ENDDO
!
      END SUBROUTINE vdW_fleur
!
      END MODULE m_vdWfleur
