      MODULE m_genernewcell
!
      CONTAINS
      SUBROUTINE gener_new_cell(
     >                          atom,sc_old,irep,ifxyz,pfxyz,
     <                          atom_new)

      USE m_types, ONLY: sc_data,atom_data

      IMPLICIT NONE

      INTEGER, intent(IN)         :: ifxyz
      INTEGER, intent(IN)         :: irep(3)
      TYPE(  sc_data),intent(IN)  :: sc_old
      TYPE(atom_data),intent(IN)  :: atom
      TYPE(atom_data),intent(OUT) :: atom_new
      CHARACTER(LEN=60),intent(IN):: pfxyz


      INTEGER :: irepeat1,irepeat2,irepeat3
      INTEGER :: i1,i2,j1,j2,j3
      INTEGER :: iatom,iline
      INTEGER :: NrNewAtoms
!
      irepeat1=irep(1)
      irepeat2=irep(2)
      irepeat3=irep(3)
!
!      print'(3(A11,I3))'," irepeat1= ",irepeat1, 
!     +                   " irepeat2= ",irepeat2, 
!     +                   " irepeat3= ",irepeat3
!
      NrNewAtoms=0
!
      do i1=1,size(atom%nr_atom_type)
        do i2=1,atom%nr_atom_type(i1)
	  do j1=-irepeat1,irepeat1
	    do j2=-irepeat2,irepeat2
	      do j3=-irepeat3,irepeat3
	        NrNewAtoms=NrNewAtoms+1
	      enddo
	    enddo
	  enddo
	enddo
      enddo
!
!      print*,'NrNewAtoms: ',NrNewAtoms
!
      allocate(atom_new%atomic_number(1:NrNewAtoms))
      allocate(atom_new%coord_cart(3,NrNewAtoms))
      allocate(atom_new%coord_bravais(3,NrNewAtoms))
!
      OPEN(ifxyz,FILE=pfxyz)
      WRITE(ifxyz,'(4i10)') NrNewAtoms,irep(:)
      WRITE(ifxyz,*)
!
      iatom=0
      iline=0
      do i1=1,size(atom%nr_atom_type)	! how many different types of atoms
        do i2=1,atom%nr_atom_type(i1)	! how many atoms for each atom type
!
! iline is used to count the number of atoms in the initial unit cell
!
	  iline=iline+1
!	  print'(I5,3F20.14)',atom%atomic_number(iline), 
!     +                      atom%coord_bravais(1:3,iline)
!         print'(5x,3F20.14)',atom%coord_cart(1:3,iline)

	  do j1=-irepeat1,irepeat1
	    do j2=-irepeat2,irepeat2
	      do j3=-irepeat3,irepeat3
	        iatom=iatom+1
		atom_new%atomic_number(iatom)=atom%atomic_number(iline)
		atom_new%coord_bravais(1,iatom)= 
     +           atom%coord_bravais(1,iline)+real(j1)
		atom_new%coord_bravais(2,iatom)= 
     +           atom%coord_bravais(2,iline)+real(j2)
		atom_new%coord_bravais(3,iatom)= 
     +           atom%coord_bravais(3,iline)+real(j3) 
		atom_new%coord_cart(1,iatom)=              
     +           atom_new%coord_bravais(1,iatom)*sc_old%bravais(1,1)+ 
     +           atom_new%coord_bravais(2,iatom)*sc_old%bravais(1,2)+ 
     +           atom_new%coord_bravais(3,iatom)*sc_old%bravais(1,3)
	        atom_new%coord_cart(2,iatom)=              
     +           atom_new%coord_bravais(1,iatom)*sc_old%bravais(2,1)+ 
     +           atom_new%coord_bravais(2,iatom)*sc_old%bravais(2,2)+ 
     +           atom_new%coord_bravais(3,iatom)*sc_old%bravais(2,3)
	        atom_new%coord_cart(3,iatom)=              
     +           atom_new%coord_bravais(1,iatom)*sc_old%bravais(3,1)+ 
     +           atom_new%coord_bravais(2,iatom)*sc_old%bravais(3,2)+ 
     +           atom_new%coord_bravais(3,iatom)*sc_old%bravais(3,3)
                write(ifxyz,'(I2,2x,3F20.14)') 
     +           atom_new%atomic_number(iatom),
     +           atom_new%coord_cart(1:3,iatom)
	      enddo
	    enddo
	  enddo
	enddo
      enddo
!
      close(ifxyz)
!
      END SUBROUTINE gener_new_cell

      END MODULE m_genernewcell
