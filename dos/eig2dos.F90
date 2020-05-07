module m_eig2dos
  implicit none
CONTAINS
  function generate_grid(emin_arg,emax_arg,efermi_arg)result(e)
    REAL,INTENT(IN):: emin_arg,emax_arg,efermi_arg
    REAL,ALLOCATABLE :: e(:)

    INTEGER, PARAMETER ::  ned = 1301
    REAL :: emin,emax,efermi,de
    INTEGER :: i
! scale energies
      emin =emin_arg*hartree_to_ev_const
      emax =emax_arg*hartree_to_ev_const
      efermi = efermiarg*hartree_to_ev_const
!
!     create energy grid
!
      emax = max(emin,emax) - efermi
      emin = min(emax,emin) - efermi
      de = (emax-emin)/(ned-1)

      ALLOCATE(e(ned+1))
      DO i=1,ned+1
         e(i) = emin + (i-1)*de
      ENDDO

  end function

  subroutine transform_to_dos(eigdesc,banddos,kpts,ev,e,dos_data)
    class(t_eigdesc),INTENT(IN):: eigdesc(:) !These are e.g. of t_dos,t_orbcomp...
    real,intent(in):: e(:) ! The grid
    real,intent(in):: ev(:,:,:) !eigenvalues
    type(t_dos_data),INTENT(OUT):: dos_data(:) !the dos generated for each weight in eigdesc

    DO i=1,size(eigdesc)
      DO n=1,eigdesc(n)%num_weights()
        select case(banddos%dos_mode)
        case ("hist")
          CALL dos_bin(e,results%neig,kpts%wtkpt,ev,eigdesc(i)%get_weights(n), dos_data(i)%dos(:,:,n))
        case ("tetra")
          IF ( input%film ) THEN
            CALL ptdos..
          ELSE
            CALL tetra_dos...
          ENDIF

  end subroutine
end module
