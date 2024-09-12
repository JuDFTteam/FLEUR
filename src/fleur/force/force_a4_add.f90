!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_force_a4_add
   USE m_juDFT
   USE m_types
   !------------------------------------------------------------------------------
   ! Adds the force contribution from core-tails in addition to formula A4 of Yu, 
   ! Singh & Krakauer minus a surface term that was included conveniently in
   ! force_a4.f90. The contributions are calculated in cdnovlp.F90.
   ! 
   ! It evaluates the term int [nabla rho_core]V d^3r for the core density of each
   ! atom in each region but its own Muffin Tin (i.e. the interstitial and the
   ! other MT).
   ! The corresponding formula is extended over the unit cell.
   ! Klueppelberg Sep'12 (force level 1)
   !------------------------------------------------------------------------------

   IMPLICIT NONE

   REAL   , ALLOCATABLE, PUBLIC, SAVE :: force_a4_mt(:,:,:)
   COMPLEX, ALLOCATABLE, PUBLIC, SAVE :: force_a4_is(:,:,:)

CONTAINS
   SUBROUTINE alloc_fa4_arrays(atoms,input)
      ! This subroutine allocates the arrays filled in cdnovlp.F90
      ! so their content can be provided in totale.f90, where the addition takes place

      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_atoms), INTENT(IN) :: atoms

      CALL timestart("force_a4_alloc")

      ALLOCATE(force_a4_mt(3,atoms%ntype,input%jspins), &
               force_a4_is(3,atoms%ntype,input%jspins))

      CALL timestop("force_a4_alloc")

   END SUBROUTINE alloc_fa4_arrays

   SUBROUTINE force_a4_add(atoms,input,results)
      ! This subroutine adds the coretail contribution to the atomic forces.

      USE m_constants

      TYPE(t_input),   INTENT(IN)    :: input
      TYPE(t_atoms),   INTENT(IN)    :: atoms
      TYPE(t_results), INTENT(INOUT) :: results

      INTEGER :: jsp,n,dir

      COMPLEX :: force_a4_mt_loc(3,atoms%ntype,input%jspins)
   
      CALL timestart("force_a4_add")

      ! The following film part needs to be implemented into cdnovlp, since
      ! the coretail force contribution has been relocated to cdnovlp.F
      !       CALL cpu_time(time1) TODO: Readd this film part.
      ! !**********************************************************************
      ! !     in case of film calculations: calculate integral over vacuum
      ! !**********************************************************************
      !       force_a4_1d = czero
      !       force_a4_2d = czero
      !       IF (film.AND..not.odi%d1) THEN
      !       DO i = 1,max(nmzxy,nmz)
      !         z(i) = i*delz
      !       END DO ! i
      !       DO jsp = 1,jspins
      !       DO ivac = 1,nvac
      !       sumG = czero
      !       integrand2 = czero
      !       DO k = 2,nq3
      !         k1 = kv3(1,k)
      !         k2 = kv3(2,k)
      !         k3 = kv3(3,k)
      !         ind2 = ig2(k)
      ! !         IF (ind2.EQ.1) CYCLE
      ! 
      !         CALL spgrot(
      !      >              nop,symor,tpi,mrot,tau,invtab,
      !      >              kv3(:,k),
      !      <              kr3d,phas3d)
      ! 
      !         CALL spgrot(
      !      >              nop2,symor,tpi,mrot,tau,invtab,
      !      >              (/k1,k2,0/),
      !      <              kr2d,phas2d)
      ! 
      !         carg = ci * k3 * bmat(3,3) ! k3 in external coordinates
      ! 
      !         pot = 0.0
      !         IF (ind2.EQ.1) THEN
      !           pot(1:nmz  ) = vz( 1:nmz         ,ivac,jsp)
      !         ELSE
      !           pot(1:nmzxy) = vxy(1:nmzxy,ind2-1,ivac,jsp)
      !         END IF
      ! 
      ! 
      !         DO j = 1,nop
      !           kp = (k-1)*nop + j
      !         DO jp = 1,nop2
      !           IF (kr3d(1,j).NE.kr2d(1,jp)) CYCLE
      !           IF (kr3d(2,j).NE.kr2d(2,jp)) CYCLE
      !           sumG(:,kp) = sumG(:,kp) +exp(carg*z(:))*pot(:)*phas2d(jp)/nop2
      !         END DO ! jp 2d operations
      !         END DO ! j operations
      ! 
      !         IF (k.EQ.nq3) THEN
      !         DO n = 1,ntype
      !           DO j = 1,nop
      !             kp = (k-1)*nop + j
      !             DO dir = 1,3
      !               integrand2(dir,1:mshd,n) = integrand2(dir,1:mshd,n)
      !      +                                 + dqpwcatom(dir,kp,n,jsp)
      !      *                                 * sumG(1:mshd,kp)
      !             END DO ! directions
      !           END DO ! j operations
      !         END DO ! n types of atoms
      !         END IF
      ! 
      !       END DO ! k stars
      ! 
      !       DO n = 1,ntype
      !       DO dir = 1,3
      !        CALL qsf(delz, real(integrand2(dir,:,n)),integral2(1,dir),mshd,0)
      !        CALL qsf(delz,aimag(integrand2(dir,:,n)),integral2(2,dir),mshd,0)
      !         integral2(:,dir) = integral2(:,dir) * area
      !         force_a4_2d(dir,n,jsp) = force_a4_2d(dir,n,jsp)
      !      +                         + integral2(1,dir) + ci* integral2(2,dir)
      !       END DO ! dir ections
      !       END DO ! n types of atoms
      ! 
      !       END DO ! ivac vacuum regions (1 or 2)
      !       END DO ! jsp spins
      !       END IF
      !       CALL cpu_time(time2)
      !       WRITE (*,*) 'time to calculate lower dimensional contributions:',
      !      + time2-time1

      ! Add the results to the total force.

      force_a4_mt_loc=force_a4_mt*CMPLX(1.0,0.0)
      DO jsp = 1, input%jspins
         DO n = 1, atoms%ntype
            IF (atoms%l_geo(n)) THEN

               results%force(:,n,jsp) = results%force(:,n,jsp) + real(force_a4_is(:,n,jsp)) + real(force_a4_mt(:,n,jsp))
               ! + real(force_a4_2d(:,n,jsp)) ! is calculated only if film.and..not.odi%d1, otherwise 0
               ! + real(force_a4_1d(:,n,jsp)) ! is calculated only if odi%d1, otherwise 0

               ! Write force contributions to out file.
               WRITE (oUnit,FMT=8010) n
               WRITE (oUnit,FMT=8020) ((force_a4_is(dir,n,jsp)),dir=1,3) ! 8020
               WRITE (oUnit,FMT=8015) n
               WRITE (oUnit,FMT=8070) ((force_a4_mt_loc(dir,n,jsp)),dir=1,3) ! 8070
8010           FORMAT (' FORCES: IS ADDITION TO EQUATION A4 FOR ATOM TYPE',i4)
8015           FORMAT (' FORCES: MT ADDITION TO EQUATION A4 FOR ATOM TYPE',i4)
!8025           FORMAT (' FORCES: VACUUM ADD. TO EQUATION A4 FOR ATOM TYPE',i4)
!8020           FORMAT (' FX_A4=',2f19.15,' FY_A4=',2f19.15,' FZ_A4=',2f19.15)
8020           FORMAT (' FX_IS=',2f10.6,' FY_IS=',2f10.6,' FZ_IS=',2f10.6)
8070           FORMAT (' FX_MT=',2f10.6,' FY_MT=',2f10.6,' FZ_MT=',2f10.6)

            END IF ! atoms%l_geo(n)

         END DO ! n types of atoms
      END DO ! jsp spins
    
      DEALLOCATE ( force_a4_mt, force_a4_is )
      CALL timestop("force_a4_add")

   END SUBROUTINE force_a4_add

END MODULE m_force_a4_add
