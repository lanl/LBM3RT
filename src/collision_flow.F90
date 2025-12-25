module collision_flow_module

  contains

      subroutine collisionf(fi,rho,ux,uy,uz,walls)

      use equilf_module
      
      use comvar_module
      
      implicit none
      
      integer :: i,j,k,s,istep

      real*8,dimension(0:18)::feq,ffb,tmp
      
      real*8, dimension(mx,my+2,mz+2,0:18) :: fi
      real*8, dimension(mx,my+2,mz+2) :: rho,ux,uy,uz
      logical,dimension(mx,my+2,mz+2) :: walls
      


      do i = 1,mx
        do j = 1,my+2
          do k = 1,mz+2
            if(.not.walls(i,j,k)) then

              call equilf(feq,rho(i,j,k),ux(i,j,k),uy(i,j,k),uz(i,j,k))
              do s = 0,18 
                tmp(s) = cix(s)*ux(i,j,k) + ciy(s)*uy(i,j,k) + ciz(s)*uz(i,j,k)
                ffb(s) = wts(s)*((3.0d0*(cix(s)-ux(i,j,k)) + 9.0d0*tmp(s)*cix(s))*ffx)

              enddo

              fi(i,j,k,:)=fi(i,j,k,:)-1.d0/tau*(fi(i,j,k,:)-feq(:)) + ffb(:)*(1.0d0-0.5d0/tau)
            endif

          enddo
        enddo
      enddo


      ! do i = 1,mx
      !   do j = 1,my+2
      !     do k = 1,mz+2
      !       if(.not.walls(i,j,k)) then

      !         call equilf(feq,rho(i,j,k),ux(i,j,k),uy(i,j,k),uz(i,j,k))
      !         do s = 0,18 
                
      !           ffb(s) = wts(s)*cix(s)*ffx*3.0d0

      !         enddo

      !         fi(i,j,k,:)=fi(i,j,k,:)-1.d0/tau*(fi(i,j,k,:)-feq(:)) + ffb(:)*(1.0d0-0.5d0/tau)
      !       endif

      !     enddo
      !   enddo
      ! enddo
    
      end subroutine collisionf
end module collision_flow_module