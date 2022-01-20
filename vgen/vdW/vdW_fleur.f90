MODULE m_vdWfleur_grimme
    !
    !   implements Grimmes D2 and D3 based on routines from V. Caciuc (`19)
    ! 
    
    !
    ! Types for vdW-forces
    !
    PRIVATE
    TYPE atom_data
    INTEGER, DIMENSION(:),  POINTER :: nr_atom_type
    INTEGER, DIMENSION(:),  POINTER :: atomic_number
    REAL,    DIMENSION(:,:),POINTER :: coord_bravais
    REAL,    DIMENSION(:,:),POINTER :: coord_cart
    END TYPE atom_data
    
    PUBLIC :: vdw_fleur_grimme

    CONTAINS
    SUBROUTINE vdW_fleur_grimme(atoms,sym,cell,film,e_vdW,f_vdW)
        USE m_constants,only:oUnit,tpi_const
        USE m_types,only: t_atoms,t_cell,t_sym,t_xcpot
        USE DFT_D2,  ONLY: driver_DFT_D2
        USE DFT_D3,  ONLY: driver_DFT_D3
       
        IMPLICIT NONE
        TYPE(t_atoms),INTENT(IN) :: atoms
        TYPE(t_cell),INTENT(IN)  :: cell
        TYPE(t_sym),INTENT(IN)   :: sym

        LOGICAL, INTENT (IN) :: film
        REAL,    INTENT (OUT) :: e_vdW,f_vdW(:,:)
        
        INTEGER  :: NrAtomType,nsize,max_cyc,iop
        INTEGER  :: NrAtoms,i,j,i1,i2,na
        INTEGER,SAVE:: irep(3)=0 !will be determined at first call and reused later
        REAL :: start,finish,toler,delta
        REAL :: test(3,8),brmin(3),brmax(3),force_i(3),f_rot(3)
        LOGICAL l_D2,l_in_au
        TYPE(atom_data) :: atom,atom_new
        REAL, ALLOCATABLE :: ener(:),force_vdW(:,:),force_max(:,:)
        
                   
        max_cyc=merge(20,1,all(irep==0)) ! max. tries for bigger supercells
        toler = 0.0005            ! energy tolerance required (eV)
        ALLOCATE ( ener(max_cyc) )
        
        l_D2 = .false.
        l_in_au = .true. ! .false.
        
        NrAtomType=atoms%ntype
        ALLOCATE( atom%nr_atom_type(NrAtomType) )                ! neq(atoms%ntype)
        atom%nr_atom_type(:) = atoms%neq(:)
        NrAtoms = atoms%nat
        ALLOCATE( atom%atomic_number(NrAtoms) )                  ! like zatom(natd) not (atoms%ntype)
        ALLOCATE( atom%coord_bravais(3,NrAtoms) )                ! taual(3,natd)
        ALLOCATE( atom%coord_cart(3,NrAtoms) )                   ! pos(3,natd)
        ALLOCATE( force_vdW(3,NrAtoms),force_max(NrAtoms,max_cyc) )
        
        na = 0
        DO i1=1,atoms%ntype
            DO i2 = 1,atoms%neq(i1)
                na = na + 1
                atom%atomic_number(na) = NINT(atoms%zatom(i1))
                atom%coord_bravais(:,na)=matmul(cell%atoms%pos(:,na))/tpi_const
                atom%coord_cart(:,na)= atoms%pos(:,na)
            ENDDO
        ENDDO
        

        
        IF (max_cyc > 1) THEN ! determine a supercell size where the vdW energy converges
            
            test(:,1) = 0.0
            test(:,2) = cell%amat(1,:) ! tr!
            test(:,3) = cell%amat(2,:)
            test(:,4) = cell%amat(3,:)
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
            ! generate supercell for vdW  
        CALL gener_new_cell(atom,cell,irep,atom_new)
        !
        IF (l_D2) THEN
             ! using DFT-D2 method
            CALL driver_DFT_D2(atom%coord_cart,atom%atomic_number,atom_new%coord_cart,atom_new%atomic_number)
        ENDIF
        !
        ! note that in the driver_DFT_D3 atom%coord_cart and atom_new%coord_cart
        !   are modified from Angstrom to au [they have INTENT(INOUT)]
        ! using DFT-D3 method 
        CALL driver_DFT_D3( atom%coord_cart,atom%atomic_number,atom_new%coord_cart,atom_new%atomic_number,l_in_au,ener(nsize),force_vdW)
        
        force_max(:,nsize) = sqrt(force_vdW(1,:)**2 + force_vdW(2,:)**2+ force_vdW(3,:)**2)
        
        IF (max_cyc == 1) THEN ! data taken from file
            delta = 0.0
            EXIT cyc
        ELSE
            IF (nsize > 1) THEN
                delta = ener(nsize)-ener(nsize-1)
                WRITE(oUnit,*) 'Delta = ',delta,' eV'
                DO i1=1,NrAtoms
                    WRITE (oUnit,'(a15,i4,a3,3f15.9)') 'Delta vdW force',i1,' : ',force_max(i1,nsize)-force_max(i1,nsize-1)
                ENDDO
                IF (abs(delta) < toler) EXIT cyc
            ENDIF
        ENDIF
        
        irep(1:2) = irep(1:2) + 1
        IF (.not.film) irep(3) =  irep(3) + 1
    ENDDO cyc
    
    IF (abs(delta) > toler) THEN
        WRITE (oUnit,*) 'vdW did not converge with cell size!'
        e_vdW = 0.0
    ELSE
        WRITE (oUnit,'(a13,3i5,a4,f12.3,a5)') 'vdW converged',irep(:)
        e_vdW = ener(nsize)/27.21138386
        WRITE ( oUnit,8060) e_vdW
        DO i1=1,NrAtoms
            WRITE (oUnit,'(a15,i4,a3,6f15.9)') 'vdW force atom ',i1,' : ',force_vdW(:,i1),atom%coord_cart(:,i1)
        ENDDO
    ENDIF
    8060 FORMAT (/,/,' ----> vdW (D3)  energy=',t40,f20.10,' htr')
    
    na = 0
    DO i1=1,NrAtomType
        f_rot=0.0
        DO i2 = 1,atoms%neq(i1)   ! here symmetrization should be done
            na = na + 1
            iop = sym%ngopr(na) ! invtab(ngopr(na))
            force_i=matmul(cell%bmat,force_vdW(:,na))/tpi_const
            f_rot=f_rot+matmul(real(sym%mrot(:,:,iop)),force_i)
        ENDDO
        f_vdW(:,i1)=matmul(cell%amat,f_rot)/atoms%neq(i1)
    ENDDO
    !
END SUBROUTINE vdW_fleur_grimme
!

SUBROUTINE gener_new_cell(atom,cell,irep,atom_new)
    
    USE m_types,only: t_cell
        
    IMPLICIT NONE
    
    INTEGER, intent(IN)         :: irep(3)
    TYPE(t_cell),intent(IN)     :: cell
    TYPE(atom_data),intent(IN)  :: atom
    TYPE(atom_data),intent(OUT) :: atom_new
    
    
    INTEGER :: irepeat1,irepeat2,irepeat3
    INTEGER :: i1,i2,j1,j2,j3
    INTEGER :: iatom,iline
    INTEGER :: NrNewAtoms
    !
    irepeat1=irep(1)
    irepeat2=irep(2)
    irepeat3=irep(3)

    !
    NrNewAtoms=(2*irepeat1+1)*(2*irepeat2+1)*(2*irepeat3+1)*sum(atom%nr_atom_type)
   
    !
    allocate(atom_new%atomic_number(1:NrNewAtoms))
    allocate(atom_new%coord_cart(3,NrNewAtoms))
    allocate(atom_new%coord_bravais(3,NrNewAtoms))
    !
    !
    iatom=0
    iline=0
    do i1=1,size(atom%nr_atom_type)	! how many different types of atoms
        do i2=1,atom%nr_atom_type(i1)	! how many atoms for each atom type
            ! iline is used to count the number of atoms in the initial unit cell
            iline=iline+1
            do j1=-irepeat1,irepeat1
                do j2=-irepeat2,irepeat2
                    do j3=-irepeat3,irepeat3
                        iatom=iatom+1
                        atom_new%atomic_number(iatom)=atom%atomic_number(iline)
                        atom_new%coord_bravais(:,iatom)=atom%coord_bravais(:,iline)+(/real(j1),real(j2),real(j3)/)
                        atom_new%coord_cart(:,iatom)=matmul(cell%amat,atom_new%coord_bravais(:,iatom))
                    enddo
                enddo
            enddo
        enddo
    enddo
    !
    !
END SUBROUTINE gener_new_cell

END MODULE m_vdWfleur_grimme   
