!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_gf_vacuum_charge
    use m_juDFT
    implicit none
    private
    integer,parameter      :: nmz=250,nmzxy=100
    real,parameter         :: delz=0.1
!Saved variable
    complex,allocatable,save :: qz(:,:)
    complex,allocatable,save :: qxy(:,:,:)

    public gf_vacuum_makecharge,gf_vacuum_writecharge
    contains
    recursive subroutine gf_vacuum_makecharge(fmpi,vacuum,input,nococonv,stars,lapw,lapw_gf,sym,cell,kpts,enpara,layers,l_noco,jspin,jspins,nk)
    ! subroutine that generates a new vacuum charge for all energies and the current spin&kpoint
    use m_gf_iodop
    use m_gf_types
    use m_gf_embedding
    use m_gf_energies
    use m_gf_math
    use m_gf_vacuum_hs
    use m_vacuz
    use m_vacudz
    USE m_constants, ONLY: pi_const, oUnit

    implicit none
    !Arguments
    type(t_mpi),intent(in)   :: fmpi
    type(t_vacuum),intent(in):: vacuum
    type(t_input),intent(in) :: input
    type(t_nococonv),intent(in):: nococonv
    type(t_stars),intent(in) :: stars
    type(t_cell),intent(in)  :: cell
    type(t_lapw),intent(in)  :: lapw
    type(t_lapw_gf),intent(in) :: lapw_gf
    type(t_sym),intent(in)   :: sym
    type(t_enpara),intent(in):: enpara
    type(t_kpts),intent(in)  :: kpts
    type(t_layers),intent(in) :: layers
    integer,intent(in)       :: jspins,jspin,nk
    logical,intent(in)       :: l_noco

    !Local variables
    integer  :: en,n1,n2,n,nv2,ng,n3d_star,nn
    real     :: v(3),ev,vz(nmz,jspins)
    complex  :: scale,vxy(nmzxy,stars%ng2-1,1)
    complex  :: h_mat(lapw_gf%nv2(jspin)*2,lapw_gf%nv2(jspin)*2)
    complex  :: s_mat(lapw_gf%nv2(jspin)*2,lapw_gf%nv2(jspin)*2)
    complex  :: g(lapw_gf%nv2(jspin)*2,lapw_gf%nv2(jspin)*2)
    complex  :: sigma(lapw_gf%nv2(jspin),lapw_gf%nv2(jspin))
    complex  :: g1(nmz),g2(nmz)

    real,allocatable    :: u(:,:,:),ud(:,:,:)
    real,allocatable    :: uz(:,:),udz(:,:),duz(:,:),dudz(:,:),ddnv(:,:)

    real :: norm,evac

    norm=pi_const*cell%area

    if (jspins==1) norm=norm/2

    nv2=lapw_gf%nv2(jspin)

    if (l_noco) then
        call juDFT_warn("Generation of vacuum charge is spin-diagonal in noco-case")
        call gf_vacuum_makecharge(fmpi,vacuum,input,nococonv,stars,lapw,lapw_gf,sym,cell,kpts,enpara,layers,.false.,1,jspins,nk)
        call gf_vacuum_makecharge(fmpi,vacuum,input,nococonv,stars,lapw,lapw_gf,sym,cell,kpts,enpara,layers,.false.,2,jspins,nk)
        return
    endif


    !read the vacuum potential
    call gf_iodop_readvacuum(GF_POTFILE,vxy,vz)

    !Generate z-dependent basis functions
    ALLOCATE(u(nmz,lapw_gf%nv2(jspin),jspins),ud(nmz,lapw_gf%nv2(jspin),jspins))
    allocate(uz(nv2,jspins),udz(nv2,jspins),duz(nv2,jspins),dudz(nv2,jspins),ddnv(nv2,jspins))


    if (nk==1.and.jspin==1) then
        allocate(qz(nmz,jspins),qxy(nmzxy,stars%ng2-1,jspins))
        qz=0.0
        qxy=0.0
    endif

        evac=enpara%evac(1,jspin)
        !evac=gf_z(en,layers%num_layers)
        if (evac>vz(nmz,jspin)-0.01) then
            call juDFT_warn("Vacuum level below highest energy point")
            evac=vz(nmz,jspin)-0.01
        endif
        DO ng = 1,lapw_gf%nv2(jspin)
            v(1:2)=kpts%bk(1:2,nk)+(/lapw_gf%k1p(ng,jspin),lapw_gf%k2p(ng,jspin)/)
            v(3) = 0.0
            ev = evac - 0.5*dot_product(v,matmul(transpose(cell%bbmat),v))
            CALL vacuz(ev,vz(:,jspin),vz(nmz,jspin),nmz,delz,uz(ng,jspin),duz(ng,jspin),u(:,ng,jspin))
            CALL vacudz(ev,vz(:,jspin),vz(nmz,jspin),nmz,delz,udz(ng,jspin),dudz(ng,jspin),ddnv(ng,jspin),  &
    &                  ud(:,ng,jspin),duz(ng,jspin),u(:,ng,jspin))
!--->       make sure the solutions satisfy the wronksian
            scale = 2.0/ (udz(ng,jspin)*duz(ng,jspin)-dudz(ng,jspin)*uz(ng,jspin))
            udz(ng,jspin) = scale*udz(ng,jspin)
            dudz(ng,jspin) = scale*dudz(ng,jspin)
            ddnv(ng,jspin) = scale*ddnv(ng,jspin)
            ud(:,ng,jspin) = scale*ud(:,ng,jspin)
        end do
        !Generate Hamiltonian & overlapp matrix
        Call gf_vacuum_getHS(fmpi,vacuum,input,nococonv,stars,lapw,lapw_gf,jspin,jspins,nk,l_noco,sym,kpts,evac,cell,h_mat,s_mat,l_sym=.false.)



    DO en=1,gf_noen()
        g=0.0
        !call juDFT_error("Not implemeted")
        call gf_GETEMB2(sigma,1,layers%num_layers+1,en,nk,jspin,lapw,lapw_gf)
        sigma=sigma
        g=gf_z(en,layers%num_layers)*s_mat-h_mat
        !g=conjg(g)
            DO n1=1,nv2
                DO n2=1,nv2
                    g(n1,n2)         =g(n1,n2)        +0.5*sigma(n1,n2)*uz(n1,jspin)*uz(n2,jspin)
                    g(n1+nv2,n2)     =g(n1+nv2,n2)    +0.5*sigma(n1,n2)*udz(n1,jspin)*uz(n2,jspin)
                    g(n1,n2+nv2)     =g(n1,nv2+n2)    +0.5*sigma(n1,n2)*uz(n1,jspin)*udz(n2,jspin)
                    g(n1+nv2,n2+nv2) =g(n1+nv2,n2+nv2)+0.5*sigma(n1,n2)*udz(n1,jspin)*udz(n2,jspin)
                enddo
            enddo
            !g=0.5*(g+transpose(conjg(g)))


            g(:nv2,:nv2)=mat_inverse(g(:nv2,:nv2))
            !g=mat_inverse(g)

            DO n=1,stars%ng2
                g1=0;g2=0
                DO n1=1,nv2
                    DO n2=1,nv2
                        n3d_star=stars%ig(lapw_gf%k1p(n1,jspin)-lapw_gf%k1p(n2,jspin),lapw_gf%k2p(n1,jspin)-lapw_gf%k2p(n2,jspin),0)
                        if (n3d_star/=0) then
                           if (stars%ig2(n3d_star)==n)then
                              g1=g1+g(n1,n2)*u(:,n1,jspin)*u(:,n2,jspin)
              !                g1=g1+g(n1+nv2,n2+nv2)*ud(:,n1,jspin)*ud(:,n2,jspin)
              !                g1=g1+g(n1,n2+nv2)*u(:,n1,jspin)*ud(:,n2,jspin)
              !                g1=g1+g(n1+nv2,n2)*ud(:,n1,jspin)*u(:,n2,jspin)
                              !Debug
                           endif
                        endif
                        n3d_star=stars%ig(lapw_gf%k1p(n2,jspin)-lapw_gf%k1p(n1,jspin),lapw_gf%k2p(n2,jspin)-lapw_gf%k2p(n1,jspin),0)
                        if (n3d_star/=0) then
                           if (stars%ig2(n3d_star)==n) then
                               g2=g2+g(n1,n2)*u(:,n1,jspin)*u(:,n2,jspin)
             !                  g2=g2+g(n1+nv2,n2+nv2)*ud(:,n1,jspin)*ud(:,n2,jspin)
             !                  g2=g2+g(n1,n2+nv2)*u(:,n1,jspin)*ud(:,n2,jspin)
             !                  g2=g2+g(n1+nv2,n2)*ud(:,n1,jspin)*u(:,n2,jspin)

                           endif
                        endif
                    enddo
                enddo
                if (n==1) then
                    qz(:,jspin)=qz(:,jspin)+kpts%wtkpt(nk)*gf_weightz(en)*(g1-conjg(g2))*cmplx(0,0.5)/norm/stars%nstr2(n)
                else
                    qxy(:,n-1,jspin)=qxy(:,n-1,jspin)+kpts%wtkpt(nk)*gf_weightz(en)*cmplx(0,0.5)*(g1(:nmzxy)-cmplx(g2(:nmzxy)))/norm/stars%nstr2(n)
                endif
            enddo
        enddo

    end subroutine

    subroutine gf_vacuum_writecharge(mpi)
    !write out the new vacuum charge
    use m_gf_types
    use m_gf_iodop
    implicit none
    TYPE(t_gfmpi),intent(in)   :: mpi
    call gf_iodop_writevacuum(GF_CDNFILE,qxy,real(qz),mpi%self_subcom)
    deallocate(qz,qxy)
    end subroutine
end module m_gf_vacuum_charge
