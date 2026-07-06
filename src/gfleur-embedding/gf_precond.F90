!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_gf_precond
    use m_juDFT
    implicit none
    contains
    subroutine gf_precond(layer,ld,gmpi,mix,jspins,d_qpw,d_rho)
    use m_fleur_psqpw
    use m_qpwtonmt
    use m_gf_types
    use m_gf_stepsanaly
    use m_gf_precond_stepweight
    implicit none
    integer,intent(in)        :: layer
    type(t_gflayer),intent(in),target :: ld
    type(t_gfmpi),intent(in)  :: gmpi
    TYPE(t_gfmix),intent(in)  :: mix
    integer,intent(in)        :: jspins
    complex,intent(inout)     :: d_qpw(:,:)
    real,intent(inout)        :: d_rho(:,:,:,:)

    integer :: jspin,n
    real    :: rho4(size(d_rho,1),size(d_rho,2),size(d_rho,3),1)
    complex :: psq(size(d_qpw,1))
    type(t_potden) :: dden
    type(t_stars),pointer :: stars
    type(t_atoms),pointer :: atoms
    type(t_sphhar),pointer :: sphhar
    type(t_cell),pointer :: cell
    type(t_sym),pointer :: sym
    stars=>ld%stars; atoms=>ld%fi%atoms; sphhar=>ld%sphhar
    cell=>ld%fi%cell; sym=>ld%fi%sym
    if (btest(mix%precond,6)) then
             call gf_diagonalize_step(layer,stars)
             call gf_remove_variation(layer,stars,1E-3,d_qpw)
    endif


    DO jspin = 1,jspins


         call priv_plot(stars,d_qpw(:,jspin),1,layer,"Charge difference")

         call dden%init(ld%stars,ld%fi%atoms,ld%sphhar,ld%fi%vacuum,          &
             ld%fi%noco,1,1001)
         dden%pw(:,1)=d_qpw(:,jspin)
         dden%mt(:,:,:,1)=d_rho(:,:,:,jspin)
         call fleur_psqpw(gmpi%fmpi,ld,dden,1,psq,no_core=.true.)

         call priv_plot(stars,psq,2,layer,"Pseudo_Charge")

         if (btest(mix%precond,2)) then
             !restrict maximal change of total charge
             psq(1) = max(min(real(psq(1)),abs(mix%g0max)),-abs(mix%g0max))
             call priv_plot(stars,psq,5,layer,"After restriction")
         endif

         if (btest(mix%precond,5)) then
             !restrict maximal change of total charge
             psq(1) = psq(1)*mix%g0scale
             call priv_plot(stars,psq,5,layer,"After restriction")
         endif

         if (btest(mix%precond,3)) then
            ! apply step function
            CALL gf_initstepsanaly(stars,0)
            call gf_gspaceconvolve(layer,stars,0.0,psq,d_qpw(:,jspin))
            call priv_plot(stars,d_qpw(:,jspin),3,layer,"After Step function")
         else
            d_qpw(:,jspin)=psq
         endif

        if (btest(mix%precond,4)) then
            ! apply simple step function
            call priv_simplestep(stars,cell,d_qpw(:,jspin))
            call priv_plot(stars,d_qpw(:,jspin),4,layer,"After Simple Step function")
         endif

         if (btest(mix%precond,1)) then
            !Apply Kerker scheme
            DO n = 2,stars%ng3
                d_qpw(n,jspin) = d_qpw(n,jspin)/(1.0+mix%k_kerker**2/stars%sk3(n)**2)
            ENDDO
            call priv_plot(stars,d_qpw(:,jspin),5,layer,"After Kerker")
         endif




         rho4 = 0.0
         CALL qpw_to_nmt(sphhar,atoms,stars,sym,cell,gmpi%fmpi,1,4,     &
             d_qpw(:,jspin)-psq,rho4)
         d_rho(:,:,:,jspin)=d_rho(:,:,:,jspin)+rho4(:,:,:,1)
    enddo

    end subroutine gf_precond

    SUBROUTINE priv_plot(stars,qpw,num,layer,text)
    USE m_gf_fft_singleton
    use m_gf_types
    IMPLICIT NONE
    TYPE(t_stars),INTENT(IN)    :: stars
    COMPLEX,INTENT(IN)          :: qpw(:)
    integer,intent(in)          :: layer,num
    character(len=*),intent(in) :: text

    complex   :: rho(stars%mx3*2+1)
    INTEGER   :: n,nn,ind
    rho=0.0
    DO n=-stars%mx3,stars%mx3
       if (n<0) then
         nn=n+2*stars%mx3+2
       else
         nn=n+1
       endif
       ind=stars%ig(0,0,n)
       if (ind==0) cycle
       rho(nn)=qpw(ind)
    enddo
    rho=fft(rho,inv=.true.)

    write(88,*) text
    DO n=1,size(rho)
       write(88,"(3i4,2(1x,f0.10))") layer,num,n,rho(n)
    enddo
    end subroutine


    SUBROUTINE priv_simplestep(stars,cell,qpw)
    USE m_gf_fft_singleton
    use m_gf_types
    IMPLICIT NONE
    TYPE(t_stars),INTENT(IN)    :: stars
    type(t_cell),intent(in)     :: cell
    COMPLEX,INTENT(INOUT)       :: qpw(:)

    complex   :: rho(stars%mx3*2+1)
    INTEGER   :: n,nn,ind,nz1,nz2,n2
    integer,parameter:: n2_maxstar=6

    nz1=cell%z1/cell%amat(3,3)*(stars%mx3*2+1)+1
    nz2=(stars%mx3*2+1)-nz1


    DO n2=1,min(stars%ng2,n2_maxstar)
        rho=0.0
        DO n=-stars%mx3,stars%mx3
            if (n<0) then
                 nn=n+2*stars%mx3+2
            else
                 nn=n+1
            endif
            ind=stars%ig(stars%kv2(1,n2),stars%kv2(2,n2),n)
            if (ind==0) cycle
            rho(nn)=qpw(ind)
        enddo
        rho=fft(rho,inv=.true.)
        DO n=nz1,nz2
           rho(n)=0.0
        enddo
        rho=fft(rho)/size(rho)
        DO n=-stars%mx3,stars%mx3
            if (n<0) then
                 nn=n+2*stars%mx3+2
            else
                 nn=n+1
            endif
            ind=stars%ig(stars%kv2(1,n2),stars%kv2(2,n2),n)
            if (ind==0) cycle
            qpw(ind)=rho(nn)
        enddo
    enddo

  end subroutine
end module m_gf_precond
