MODULE m_fleurinput_mpi_bc
  USE m_types_fleurinput
  IMPLICIT NONE
CONTAINS
  SUBROUTINE fleurinput_mpi_bc(cell,sym,atoms,input,noco,vacuum,field,&
       sliceplot,banddos,mpinp,hybinp ,coreSpecInput,wann,&
       xcpot,forcetheo_data,kpts,enparaXML,gfinp,hub1inp,mpi_comm,juPhon,rank)
    USE m_types_xml


    TYPE(t_cell),INTENT(INOUT)::cell
    TYPE(t_sym),INTENT(INOUT)::sym
    TYPE(t_atoms),INTENT(INOUT)::atoms
    TYPE(t_input),INTENT(INOUT)::input
    TYPE(t_noco),INTENT(INOUT)::noco
    TYPE(t_vacuum),INTENT(INOUT)::vacuum
    TYPE(t_field),INTENT(INOUT)::field
    TYPE(t_sliceplot),INTENT(INOUT)::sliceplot
    TYPE(t_banddos),INTENT(INOUT)::banddos
    TYPE(t_mpinp), INTENT(INOUT):: mpinp
    TYPE(t_hybinp),INTENT(INOUT)::hybinp
     
    TYPE(t_coreSpecInput),INTENT(INOUT)::coreSpecInput
    TYPE(t_wann),INTENT(INOUT)::wann
    CLASS(t_xcpot),ALLOCATABLE,INTENT(INOUT)::xcpot
    TYPE(t_forcetheo_data),INTENT(INOUT)::forcetheo_data
    TYPE(t_enparaXML),INTENT(INOUT)::enparaXML
    TYPE(t_kpts),INTENT(INOUT)::kpts
    TYPE(t_gfinp),INTENT(INOUT)::gfinp
    TYPE(t_hub1inp),INTENT(INOUT)::hub1inp
    TYPE(t_juPhon),INTENT(INOUT)::juPhon
    INTEGER,INTENT(IN)::mpi_comm
    INTEGER,OPTIONAL,INTENT(IN)::rank


    CALL cell%mpi_bc(mpi_comm,rank)
    CALL sym%mpi_bc(mpi_comm,rank)
    CALL atoms%mpi_bc(mpi_comm,rank)
    CALL input%mpi_bc(mpi_comm,rank)
    CALL noco%mpi_bc(mpi_comm,rank)
    CALL vacuum%mpi_bc(mpi_comm,rank)
    CALL field%mpi_bc(mpi_comm,rank)
    CALL sliceplot%mpi_bc(mpi_comm,rank)
    CALL banddos%mpi_bc(mpi_comm,rank)
    CALL hybinp%mpi_bc(mpi_comm,rank)
    CALL mpinp%mpi_bc(mpi_comm, rank)
    CALL coreSpecInput%mpi_bc(mpi_comm,rank)
    CALL wann%mpi_bc(mpi_comm,rank)
    CALL forcetheo_data%mpi_bc(mpi_comm,rank)
    CALL enparaXML%mpi_bc(mpi_comm,rank)
    CALL kpts%mpi_bc(mpi_comm,rank)
    CALL gfinp%mpi_bc(mpi_comm,rank)
    CALL hub1inp%mpi_bc(mpi_comm,rank)
    CALL xcpot%mpi_bc(mpi_comm,rank)
    CALL juPhon%mpi_bc(mpi_comm,rank)

  END SUBROUTINE fleurinput_mpi_bc
END MODULE m_fleurinput_mpi_bc
