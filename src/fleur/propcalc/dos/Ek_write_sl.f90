MODULE m_Ekwritesl
  use m_juDFT
CONTAINS
  SUBROUTINE Ek_write_sl(eig_id,kpts,atoms,vacuum,input,jspin,sym,cell,dos,slab,orbcomp,results)
    !-----------------------------------------------------------------
    !-- now write E(k) for all kpts if on T3E
    !-- now read data from tmp_dos and write of E(k) in  ek_orbcomp
    !-----------------------------------------------------------------
    USE m_types
    IMPLICIT NONE
    
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_vacuum),INTENT(IN)      :: vacuum
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_cell),INTENT(IN)        :: cell
    TYPE(t_dos),INTENT(IN)         :: dos
    TYPE(t_kpts),INTENT(IN)        :: kpts
    TYPE(t_atoms),INTENT(IN)       :: atoms
    TYPE(t_slab),INTENT(IN)        :: slab
    TYPE(t_orbcomp),INTENT(IN)     :: orbcomp
    TYPE(t_results),INTENT(IN)     :: results
    !	..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: eig_id
    INTEGER, INTENT (IN) :: jspin
    !     ..
    !     .. Local Scalars
    INTEGER :: nbands,ikpt,kspin,j,i,n,it ,na,iband,mt,l
    INTEGER :: ivac,m
    REAL    :: wk
    !     ..
    !     .. Local Arrays
    INTEGER  norb(23),iqsl(slab%nsld),iqvacpc(2)
    REAL     qvact(2)
    REAL, ALLOCATABLE :: eig(:)
    CHARACTER (len=2) :: chntype
    CHARACTER (len=99) :: chform
    !     ..
    IF (slab%nsl.GT.slab%nsld)  THEN
       CALL juDFT_error("nsl.GT.nsld",calledby="Ek_write_sl")
    ENDIF
    ALLOCATE(eig(input%neig))
    !  --->     open files for a bandstucture with an orbital composition
    !  --->     in the case of the film geometry
    !
    IF (jspin.EQ.1)  OPEN (130,file='ek_orco_11') 
    IF (jspin.EQ.2)  OPEN (130,file='ek_orco_12')
    !
    ! ----->       write bandstructure to ek_orbcomp - file
    ! 
    WRITE (chntype,'(i2)') slab%nsl
    chform = "('E',i3,'= ',f10.4,4x,'vac ( banddos%layers ) vac = ',i3,' ('&
         &        ,"//chntype//"(i3,2x),')',i3))"
    WRITE (130,FMT=901) 
    WRITE (130,FMT=902) 
    WRITE (130,FMT=901) 
    WRITE (130,FMT=903) slab%nsl,vacuum%nvac,kpts%nkpt
    WRITE (130,FMT=904) atoms%ntype,(atoms%neq(n),n=1,atoms%ntype) 
    WRITE (130,FMT=805)  
    DO j=1,slab%nsl
       WRITE (130,FMT=806) j,(slab%nslat(i,j),i=1,atoms%nat)
    ENDDO
    DO kspin = 1,input%jspins
       WRITE (130,FMT=907)  kspin,input%jspins
       !============================================================== 
901    FORMAT (5X,'--------------------------------------')
902    FORMAT (5X,'-------- E(k) for a input%film  ------------')
903    FORMAT (5X,' nsl =',i3,'   vacuum%nvac =',i2,'   kpts%nkpt =',i4,/)
904    FORMAT (5X,' atoms%ntype = ',i3,' atoms%neq(n) = ',50i3)
907    FORMAT (/,5X,' kspin = ',i4,' input%jspins = ',i4)
805    FORMAT (5X,'  nsl   nslat(1:nate,nsli)  ')
806    FORMAT (5X,51i4)
       !==============================================================
       DO ikpt=1,kpts%nkpt

          WRITE (130,FMT=8000) (kpts%bk(i,ikpt),i=1,3)
8000      FORMAT (/,3x,'  k =',3f10.5,/)
          !
          DO iband = 1,results%neig(ikpt,kspin)
             qvact = 0.0
             DO ivac = 1,vacuum%nvac
                qvact(ivac) = dos%qvac(iband,ivac,ikpt,kspin)
             ENDDO
             IF (vacuum%nvac==1)    qvact(2) = qvact(1)
             iqvacpc(:) = nint(qvact(:)*100.0)
             DO j = 1,slab%nsl
                iqsl(j) = nint((slab%qintsl(j,iband,ikpt,kspin) + slab%qmtsl(j,iband,ikpt,kspin))*100.0) 
             ENDDO
             WRITE(130,FMT=chform) iband,results%eig(iband,ikpt,kspin),iqvacpc(2),(iqsl(l),l=1,slab%nsl),iqvacpc(1)
             WRITE(130,FMT=9) 
             WRITE(130,FMT=8)
             WRITE(130,FMT=9) 
             DO n = 1,slab%nsl
                mt=0 
                DO  it=1,atoms%ntype
                   DO  m=1,atoms%neq(it)
                      mt=mt+1	
                      na = slab%nslat(mt,n) 
                      IF (na.EQ.1) THEN
                         DO  j=1,23
                            norb(j) = nint ( orbcomp%comp(iband,j,mt,ikpt,kspin) )
                         ENDDO
                         WRITE (130,FMT=5) n,it,m,(norb(l),l=1,23),orbcomp%qmtp(iband,mt,ikpt,kspin)
                      ENDIF
                   ENDDO
                enddo
             ENDDO              ! over ( n = 1,nsl ) 
             WRITE(130,FMT=9) 
          ENDDO           ! over ( iband = 1,results%neig(ikpt,kspin) ) 
       ENDDO        ! over ( ikpt=1,kpts%nkpt )
    ENDDO	  ! over ( kspin = 1,input%jspins )  
    CLOSE (130)
    !
    ! 8040 FORMAT ('E',i3,'= ',f10.4,4x,'vac | layers | vac = ',i3,' | ',50(i3,2x),' | ',i3)
8   FORMAT('|lyr,tp,at| S | Px  Py  Pz | Dxy  Dyz  Dzx  Dx-y Dz2 |',&
         &  ' Fx3  Fy3  Fz3  Fx2y Fy2z Fz2x Fxyz| Fz2x Fz2y Fz3  Fxyz Fx2z',&
         &  ' Fx3  Fy3 |  mt  |') 
5   FORMAT('|',i3,',',i2,',',i2,'|',i3,'|',3(i3,1x),'|',&
         &        5(1x,i3,1x),'|',&
         &        7(1x,i3,1x),'|',7(1x,i3,1x),'|',f6.1,'|')
9   FORMAT(133('-'))
    !
    DEALLOCATE ( eig )

  END SUBROUTINE Ek_write_sl
END MODULE m_Ekwritesl
