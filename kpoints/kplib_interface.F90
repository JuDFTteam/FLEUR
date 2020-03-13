MODULE m_kpts_kplib
  IMPLICIT NONE

  CONTAINS
  SUBROUTINE kpts_kplib(cell,sym,kpts,minDistance)
    use iso_c_binding
    USE m_types_cell
    USE m_types_sym
    USE m_types_kpts
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_kpts),INTENT(INOUT):: kpts
    REAL,INTENT(in),OPTIONAL  :: minDistance
  INTERFACE
     FUNCTION kplib_generate(primVectorsArray,latticePointOperatorsArray,numOperators,isConventionalHexagonal,Gamma,minDistance,minSize)BIND(C,name="kplib_generate")
       use iso_c_binding
       REAL(c_double)::primVectorsArray(3,3),minDistance
       INTEGER(c_int)::latticePointOperatorsArray(3,3,*)
       INTEGER(c_int):: numOperators,minSize
       LOGICAL(c_bool):: isConventionalHexagonal,Gamma
       INTEGER(c_int):: kplib_generate
     END FUNCTION kplib_generate
  END INTERFACE
  INTERFACE
     FUNCTION kplib_setkpoints(kpoints,weights)BIND(C,name="kplib_setkpoints")
       use iso_c_binding
       REAL(c_double)::kpoints(3,*),weights(*)
       INTEGER(c_int):: kplib_setkpoints
     END FUNCTION kplib_setkpoints
  END INTERFACE

   

    INTEGER:: i
    REAL   :: mind
    i=kpts%nkpt
    mind=MERGE(minDistance,1E-10,PRESENT(minDistance))
    WRITE(kpts%comment,"(a,i0,a,f10.4,a,l1)") "kplib: minSize=",i," mindist=",mind," gamma=",kpts%l_gamma
    kpts%nkpt=kplib_generate(cell%amat,sym%mrot,sym%nop,LOGICAL(.FALSE.,kind=c_bool),LOGICAL(kpts%l_gamma,kind=c_bool),mind,i)

    ALLOCATE(kpts%bk(3,kpts%nkpt))
    ALLOCATE(kpts%wtkpt(kpts%nkpt))

    i=kplib_setkpoints(kpts%bk,kpts%wtkpt)

    
  END SUBROUTINE kpts_kplib
END MODULE m_kpts_kplib
