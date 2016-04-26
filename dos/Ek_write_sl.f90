MODULE m_Ekwritesl
  use m_juDFT
CONTAINS
  SUBROUTINE Ek_write_sl(&
       &                       dimension,kpts,atoms,vacuum,nsld,&
       &                       input,jspin,&
       &                       sym,cell,&
       &                       nsl,nslat)
    !-----------------------------------------------------------------
    !-- now write E(k) for all kpts if on T3E
    !-- now read data from tmp_dos and write of E(k) in  ek_orbcomp
    !-----------------------------------------------------------------
    USE m_types
    IMPLICIT NONE
    TYPE(t_dimension),INTENT(IN)   :: dimension
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_vacuum),INTENT(IN)      :: vacuum
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_cell),INTENT(IN)        :: cell
    TYPE(t_kpts),INTENT(IN)        :: kpts
    TYPE(t_atoms),INTENT(IN)       :: atoms
    !	..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: nsld
    INTEGER, INTENT (IN) :: nsl ,jspin 
    !     ..
    !     .. Array Arguments ..
    INTEGER, INTENT (IN) :: nslat(atoms%natd,nsld)
    !     ..
    !     .. Local Scalars
    INTEGER :: nbands,ikpt,kspin,j,i,n,it ,na,iband,mt,l
    INTEGER :: ivac,m
    REAL    :: wk
    !     ..
    !     .. Local Arrays
    INTEGER  norb(23),iqsl(nsld),iqvacpc(2)
    REAL     bkpt(3),qvact(2)
    REAL, ALLOCATABLE :: eig(:),qvac(:,:,:,:),orbcomp(:,:,:,:,:)
    REAL, ALLOCATABLE :: qintsl(:,:,:,:),qmtsl(:,:,:,:),qmtp(:,:,:,:)
    CHARACTER (len=2) :: chntype
    CHARACTER (len=99) :: chform
    !     ..
    !     .. Intrinsic Functions
    INTRINSIC  nint 
    !
    IF (nsl.GT.nsld)  THEN
       CALL juDFT_error("nsl.GT.nsld",calledby="Ek_write_sl")
    ENDIF
    ALLOCATE(eig(dimension%neigd),orbcomp(dimension%neigd,23,atoms%natd,kpts%nkptd,dimension%jspd))
    ALLOCATE(qvac(dimension%neigd,2,kpts%nkptd,dimension%jspd),qintsl(nsld,dimension%neigd,kpts%nkptd,dimension%jspd))
    ALLOCATE(qmtsl(nsld,dimension%neigd,kpts%nkptd,dimension%jspd),qmtp(dimension%neigd,atoms%natd,kpts%nkptd,dimension%jspd))
    !
    !  --->     open files for a bandstucture with an orbital composition
    !  --->     in the case of the film geometry
    !
    IF (jspin.EQ.1)  OPEN (130,file='ek_orco_11') 
    IF (jspin.EQ.2)  OPEN (130,file='ek_orco_12')
    !
    ! ----->       write bandstructure to ek_orbcomp - file
    ! 
    WRITE (chntype,'(i2)') nsl
    chform = "('E',i3,'= ',f10.4,4x,'vac ( vacuum%layers ) vac = ',i3,' ('&
         &        ,"//chntype//"(i3,2x),')',i3))"
    WRITE (130,FMT=901) 
    WRITE (130,FMT=902) 
    WRITE (130,FMT=901) 
    WRITE (130,FMT=903) nsl,vacuum%nvac,kpts%nkpt
    WRITE (130,FMT=904) atoms%ntype,(atoms%neq(n),n=1,atoms%ntype) 
    WRITE (130,FMT=805)  
    DO j=1,nsl
       WRITE (130,FMT=806) j,(nslat(i,j),i=1,atoms%natd)
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
          !                
          READ (129,rec=kpts%nkpt*(kspin-1)+ikpt)  bkpt,wk,nbands,eig,&
               &            qvac(:,:,ikpt,kspin),qintsl(:,:,ikpt,kspin),&
               &            qmtsl(:,:,ikpt,kspin),&
               &            orbcomp(:,:,:,ikpt,kspin),qmtp(:,:,ikpt,kspin)
          !            write(*,*) kspin,nkpt,qmtp(1,:,ikpt,kspin)
          !
          WRITE (130,FMT=8000) (bkpt(i),i=1,3)
8000      FORMAT (/,3x,'  k =',3f10.5,/)
          !
          DO iband = 1,nbands
             qvact = 0.0
             DO ivac = 1,vacuum%nvac
                qvact(ivac) = qvac(iband,ivac,ikpt,kspin)
             ENDDO
             IF (sym%invs .OR. sym%zrfs)    qvact(2) = qvact(1)
             iqvacpc(:) = nint(qvact(:)*100.0)
             DO j = 1,nsl
                iqsl(j) = nint( ( qintsl(j,iband,ikpt,kspin) + &
                     &                               qmtsl(j,iband,ikpt,kspin) )*100.0 ) 
             ENDDO
             WRITE (130,FMT=chform) iband,eig(iband),iqvacpc(2),&
                  &                               (iqsl(l),l=1,nsl),iqvacpc(1)
             WRITE(130,FMT=9) 
             WRITE(130,FMT=8)
             WRITE(130,FMT=9) 
             DO n = 1,nsl
                mt=0 
                DO  it=1,atoms%ntype
                   DO  m=1,atoms%neq(it)
                      mt=mt+1	
                      na = nslat(mt,n) 
                      IF (na.EQ.1) THEN
                         DO  j=1,23
                            norb(j) = &
                                 &                nint ( orbcomp(iband,j,mt,ikpt,kspin) )
                         ENDDO
                         WRITE (130,FMT=5) n,it,m,&
                              &	   		                  (norb(l),l=1,23),&
                              &                                    qmtp(iband,mt,ikpt,kspin)
                      ENDIF
                   ENDDO
                enddo
             ENDDO              ! over ( n = 1,nsl ) 
             WRITE(130,FMT=9) 
          ENDDO           ! over ( iband = 1,nbands ) 
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
    DEALLOCATE ( eig,qvac,orbcomp,qintsl,qmtsl,qmtp )

  END SUBROUTINE Ek_write_sl
END MODULE m_Ekwritesl
