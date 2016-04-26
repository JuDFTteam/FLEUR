      module m_wann_socmat
      contains
      subroutine wann_socmat(
     >               lmaxd,ntypd,nlod,natd,neigd,
     >               llod,jmtd,jspd,nlhd,soc_opt,neq,
     >               ntype,theta_in,phi_in,jspins,irank,
     >               jri,lmax,dx,rmsh,el,ello,nlo,llo,
     >               l_dulo,ulo_der,vr,
     >               acof,bcof,chelp,
     <               hsomtx)
c***********************************************************************
c     Calculate matrix elements of the spin-orbit interaction for those
c     Bloch functions out of which the Wannier functions are 
c     constructed. From these matrix elements the spin-orbit Hamiltonian
c     in the basis set of Wannier functions may be constructed.
c     
c     Frank Freimuth
c***********************************************************************
      USE m_spnorb 
      USE m_hsoham
      implicit none
      integer, intent(in) :: lmaxd
      integer, intent(in) :: ntypd
      integer, intent(in) :: nlod
      integer, intent(in) :: natd
      integer, intent(in) :: neigd

      integer, intent(in) :: llod
      integer, intent(in) :: jmtd
      integer, intent(in) :: jspd
      integer, intent(in) :: nlhd
      LOGICAL, INTENT (IN):: soc_opt(:) !(ntypd+2)
      integer, intent(in) :: neq(:) !ntype

      integer, intent(in) :: ntype
      real,    intent(in) :: theta_in
      real,    intent(in) :: phi_in
      integer, intent(in) :: jspins
      integer, intent(in) :: irank

      integer, intent(in) :: jri(:)  ! ntypd
      integer, intent(in) :: lmax(:) ! ntypd
      real,    intent(in) :: dx(:)   ! ntypd
      real,    intent(in) :: rmsh(:,:)  ! jmtd,ntypd
      real,    intent(in) :: el(0:,:,:) ! 0:lmaxd,ntypd,max(2,jspd)
      real,    intent(in) :: ello(:,:,:)! nlod,ntypd,max(2,jspd)
      integer, intent(in) :: nlo(:)     ! ntypd
      integer, intent(in) :: llo(:,:)   ! nlod,ntypd

      logical, intent(in) :: l_dulo(:,:)  ! l_dulo(nlod,ntypd)
      integer, intent(in) :: ulo_der(:,:) ! ulo_der(nlod,ntypd)
      real,    intent(in) :: vr(:,0:,:,:)    ! jmtd,0:nlhd,ntypd,jspd   
  
      complex, intent(in) :: acof(:,0:,:,:) !acof(noccbd,0:lmd,natd,jspd)
      complex, intent(in) :: bcof(:,0:,:,:) !bcof(noccbd,0:lmd,natd,jspd)
      complex, intent(in) :: chelp(-llod:,:,:,:,:) !chelp(-llod:llod,neigd,nlod,natd,jspd)

      complex, intent(out):: hsomtx(:,:,:,:) !(2,2,neigd,neigd)

      integer :: n,l,nwdd,nw,ispin,ie,na,ll1,m,lm,i,nsz(2)
      real    :: s(3),r2
      logical :: l_all
      CHARACTER*3 chntype
      real    :: theta,phi,pi

      real,allocatable :: ddn(:,:,:) ! 0:lmaxd,ntypd,jspd
      real,allocatable :: us(:,:,:)  ! 0:lmaxd,ntypd,jspd
      real,allocatable :: dus(:,:,:) ! 0:lmaxd,ntypd,jspd
      real,allocatable :: uds(:,:,:) ! 0:lmaxd,ntypd,jspd
      real,allocatable :: duds(:,:,:)! 0:lmaxd,ntypd,jspd

      real,allocatable :: ulos(:,:,:)  ! nlod,ntypd,jspd
      real,allocatable :: dulos(:,:,:) ! nlod,ntypd,jspd
      real,allocatable :: uulon(:,:,:) ! nlod,ntypd,jspd
      real,allocatable :: dulon(:,:,:) ! nlod,ntypd,jspd  

      COMPLEX,ALLOCATABLE :: ahelp(:,:,:,:,:)
      complex,allocatable :: bhelp(:,:,:,:,:)

      REAL,    ALLOCATABLE :: rsopdp(:,:,:,:),rsopdpd(:,:,:,:)
      REAL,    ALLOCATABLE :: rsopp(:,:,:,:),rsoppd(:,:,:,:)
      REAL,    ALLOCATABLE :: rsoplop(:,:,:,:)
      REAL,    ALLOCATABLE :: rsoplopd(:,:,:,:),rsopdplo(:,:,:,:)
      REAL,    ALLOCATABLE :: rsopplo(:,:,:,:),rsoploplop(:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: soangl(:,:,:,:,:,:)

      nwdd=1
      nw=1
      nsz=neigd

      IF (.true.) THEN
        theta= theta_In
        phi= phi_In
      ELSE
        pi= 4.*ATAN(1.)
        theta= -theta_In
        phi=   phi_In+pi
        ! now the definition of rotation matrices
        ! is equivalent to the def in the noco-routines
      ENDIF

      ALLOCATE ( ahelp(-lmaxd:lmaxd,lmaxd,natd,neigd,jspd) )
      ALLOCATE ( bhelp(-lmaxd:lmaxd,lmaxd,natd,neigd,jspd) )

      do ispin=1,jspd
       DO ie = 1, neigd
        DO na = 1, natd
         DO l = 1, lmaxd
          ll1 = l*(l+1)
          DO m = -l,l
             lm = ll1 + m
             ahelp(m,l,na,ie,ispin) = (acof(ie,lm,na,ispin))
             bhelp(m,l,na,ie,ispin) = (bcof(ie,lm,na,ispin))
          ENDDO !m
         ENDDO !l
        ENDDO !na
       ENDDO !ie
      enddo !ispin

      ALLOCATE( rsopdp(ntypd,lmaxd,2,2),rsopdpd(ntypd,lmaxd,2,2),
     +          rsopp(ntypd,lmaxd,2,2),rsoppd(ntypd,lmaxd,2,2),
     +          rsoplop(ntypd,nlod,2,2),rsoplopd(ntypd,nlod,2,2),
     +          rsopdplo(ntypd,nlod,2,2),rsopplo(ntypd,nlod,2,2),
     +          rsoploplop(ntypd,nlod,nlod,2,2),
     +          soangl(lmaxd,-lmaxd:lmaxd,2,lmaxd,-lmaxd:lmaxd,2) )

      ALLOCATE( us(0:lmaxd,ntypd,jspd), dus(0:lmaxd,ntypd,jspd),
     +          uds(0:lmaxd,ntypd,jspd),duds(0:lmaxd,ntypd,jspd),
     +          ddn(0:lmaxd,ntypd,jspd),
     +          ulos(nlod,ntypd,jspd),dulos(nlod,ntypd,jspd),
     +          uulon(nlod,ntypd,jspd),dulon(nlod,ntypd,jspd))

      soangl(:,:,:,:,:,:) = cmplx(0.0,0.0)
      CALL spnorb(
     >         ntypd,lmaxd,jmtd,jspd,nwdd,nlhd,nlod,
     >         theta,phi,jspins,ntype,nw,irank,
     >         jri,lmax,dx,rmsh,el,ello,nlo,llo,l_dulo,ulo_der,vr,
     >         soc_opt(ntype+2), 
     <         rsopp,rsoppd,rsopdp,rsopdpd,ddn,
     <         rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop,
     <         us,dus,uds,duds,ulos,dulos,uulon,dulon,soangl)

      l_all=.false.
      IF (irank==0) THEN
           IF ( soc_opt(ntype+1).or.l_all ) THEN
             IF (l_all) THEN
               WRITE (6,fmt='(A)') 'Only SOC contribution of certain'
     &          //' atom types included in Hamiltonian.'
             ELSE 
               WRITE (chntype,'(i3)') ntype
               WRITE (6,fmt='(A,2x,'//chntype//'l1)') 'SOC contributi'
     &          //'on of certain atom types included in Hamiltonian:',
     &          (soc_opt(n),n=1,ntype)
             ENDIF
           ELSE
             WRITE(6,fmt='(A,1x,A)') 'SOC contribution of all atom'//
     &        ' types inculded in Hamiltonian.'
           ENDIF
           IF (soc_opt(ntype+2)) THEN
             WRITE(6,fmt='(A)')
     &        'SOC Hamiltonian is constructed by neglecting B_xc.'
           ENDIF
      ENDIF
      soangl(:,:,:,:,:,:) = conjg(soangl(:,:,:,:,:,:))

      CALL hsoham(
     >    jspd,neigd,natd,lmaxd,ntypd,nlod,llod,
     >    ntype,soc_opt,jspins,nsz,lmax,neq,nlo,llo,chelp,
     >    rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop,
     >    ahelp,bhelp,rsopp,rsoppd,rsopdp,rsopdpd,soangl,
     <    hsomtx)
 
      end subroutine wann_socmat
      end module m_wann_socmat
