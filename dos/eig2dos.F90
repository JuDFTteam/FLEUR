module m_eig2dos
  implicit none
CONTAINS


  subroutine transform_to_dos(eigdesc,banddos,kpts,ev,e,dos_data)
    class(t_eigdesc),INTENT(IN)    :: eigdesc(:) !These are e.g. of t_dos,t_orbcomp...
    TYPE(t_banddos),INTENT(IN)     :: banddos
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_vacuum),INTENT(IN)      :: vacuum
    TYPE(t_noco),INTENT(IN)        :: noco
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_cell),INTENT(IN)        :: cell
    TYPE(t_results),INTENT(IN)     :: results
    TYPE(t_dos),INTENT(IN)         :: dos
    TYPE(t_kpts),INTENT(IN)        :: kpts
    TYPE(t_atoms),INTENT(IN)       :: atoms
    real,intent(in):: ev(:,:,:) !eigenvalues

    type(t_dos_data),INTENT(INOUT):: dos_data !the dos generated for each weight in all eigdesc
    
    ind=0
    DO i=1,size(eigdesc)
      DO n=1,eigdesc(n)%num_weights()
        ind=ind+1
        select case(banddos%dos_mode)
        case ("hist")
          CALL dos_bin(dos_data%e,results%neig,kpts%wtkpt,ev,eigdesc(i)%get_weights(n), dos_data%dos(:,:,ind))
        !case ("tetra")
        !  IF ( input%film ) THEN
        !    CALL ptdos..
        !  ELSE
        !    CALL tetra_dos...
        !  ENDIF
        end select

  end subroutine
end module
