MODULE m_fleurinput_mpi_bc
  USE m_types_fleurinput
  IMPLICIT NONE
CONTAINS
  SUBROUTINE fleurinput_mpi_bc(cell,sym,atoms,input,noco,vacuum,field,&
       sliceplot,banddos,hybrid,oneD,coreSpecInput,wann,&
       xcpot,forcetheo_data,kpts,enparaXML,mpi_comm,rank)
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
    TYPE(t_hybrid),INTENT(INOUT)::hybrid
    TYPE(t_oneD),INTENT(INOUT)::oneD
    TYPE(t_coreSpecInput),INTENT(INOUT)::coreSpecInput
    TYPE(t_wann),INTENT(INOUT)::wann
    CLASS(t_xcpot),INTENT(INOUT)::xcpot
    TYPE(t_forcetheo_data),INTENT(INOUT)::forcetheo_data
    TYPE(t_enparaXML),INTENT(INOUT)::enparaXML
    TYPE(t_kpts),INTENT(INOUT)::kpts
    INTEGER,INTENT(IN)::mpi_comm
    INTEGER,INTENT(IN),OPTIONAL::rank

    
    CALL cell%mpi_bc(mpi_comm,rank)
    CALL sym%mpi_bc(mpi_comm,rank)
    CALL atoms%mpi_bc(mpi_comm,rank)
    CALL input%mpi_bc(mpi_comm,rank)
    CALL noco%mpi_bc(mpi_comm,rank)
    CALL vacuum%mpi_bc(mpi_comm,rank)
    CALL field%mpi_bc(mpi_comm,rank)
    CALL sliceplot%mpi_bc(mpi_comm,rank)
    CALL banddos%mpi_bc(mpi_comm,rank)
    CALL hybrid%mpi_bc(mpi_comm,rank)
    CALL oneD%mpi_bc(mpi_comm,rank)
    CALL coreSpecInput%mpi_bc(mpi_comm,rank)
    CALL wann%mpi_bc(mpi_comm,rank)
    CALL xcpot%mpi_bc(mpi_comm,rank)
    CALL forcetheo_data%mpi_bc(mpi_comm,rank)
    CALL enparaXML%mpi_bc(mpi_comm,rank)
    CALL kpts%mpi_bc(mpi_comm,rank)

  END SUBROUTINE fleurinput_mpi_bc
END MODULE m_fleurinput_mpi_bc
