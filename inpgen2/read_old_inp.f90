MODULE m_read_old_inp
  IMPLICIT NONE
CONTAINS
  SUBROUTINE read_old_inp(input,DIMENSION,atoms,sphhar,cell,stars,sym,noco,vacuum,forcetheo,&
       sliceplot,banddos,enpara,xcpot,kpts,hybrid,&
       oneD)
    USE m_types_input
    USE m_types_atoms
    USE m_types_cell
    USE m_types_sym
    USE m_types_noco
    USE m_types_vacuum
    USE m_types_banddos
    USE m_types_hybrid
    USE m_types_xcpot_inbuild_nofunction
    USE m_types_forcetheo
    USE m_types_kpts
    USE m_types_enpara
    USE m_types_oneD
    USE m_types_sliceplot
    USE m_types_stars
    USE m_fleur_init_old



    TYPE(t_input),INTENT(OUT):: input
    TYPE(t_atoms),INTENT(OUT):: atoms
    TYPE(t_cell) ,INTENT(OUT):: cell
    TYPE(t_sym)  ,INTENT(OUT):: sym
    TYPE(t_noco) ,INTENT(OUT):: noco
    TYPE(t_vacuum) ,INTENT(OUT):: vacuum
    TYPE(t_banddos),INTENT(OUT):: banddos
    TYPE(t_hybrid) ,INTENT(OUT):: hybrid
    TYPE(t_xcpot_inbuild_nf),INTENT(OUT)::xcpot
    TYPE(t_oned) ,INTENT(OUT):: oned
    TYPE(t_sliceplot),INTENT(OUT):: sliceplot
    TYPE(t_stars),INTENT(OUT):: stars
    TYPE(t_enpara) ,INTENT(OUT):: enpara

    TYPE(t_forcetheo):: forcetheo
    TYPE(t_kpts) ,INTENT(OUT):: kpts


    TYPE(t_coreSpecInput)::coreSpecInput


    CALL fleur_init_old(&
         input,DIMENSION,atoms,sphhar,cell,stars,sym,noco,vacuum,forcetheo,&
         sliceplot,banddos,enpara,xcpot,kpts,hybrid,&
         oneD,coreSpecInput)

  END SUBROUTINE read_old_inp
 END MODULE m_read_old_inp
