module boundaryf_module

contains

      subroutine boundaryf(fi,rho,ux,uy,uz,rhol,rhor,uin,walls)
      
      use comvar_module
      
      implicit none
      
      integer :: i,j,k,ii,s
      real*8 :: nyx_inj,nzx_inj,nyx_out,nzx_out,rhol,rhor,uin

      real*8, dimension(1:mx,1:my+2,1:mz+2,0:18) :: fi
      real*8, dimension(1:mx,1:my+2,1:mz+2) :: rho,ux,uy,uz
      logical, dimension(1:mx,1:my+2,1:mz+2) :: walls
      
!      goto 100


      !Const velocity BC at inlet
      ! ii = 1
      ! uy(ii,:,:) = 0.0d0
      ! uz(ii,:,:) = 0.0d0
      ! do j=1,my
      !   do k=1,mz+2
      !     if(.not.walls(ii,j,k)) then
      !       ux(ii,j,k) = uin
      !       rho(ii,j,k) = 1.0d0/(1.0d0-uin) * (fi(3,ii,j,k) + fi(4,ii,j,k) + fi(5,ii,j,k) + fi(6,ii,j,k) &
      !                   + fi(15,ii,j,k) + fi(16,ii,j,k) + fi(17,ii,j,k) + fi(18,ii,j,k) + fi(0,ii,j,k) &
      !                   + 2.0d0*(fi(2,ii,j,k) + fi(8,ii,j,k) + fi(10,ii,j,k) + fi(12,ii,j,k) + fi(14,ii,j,k)))


      ! !Const pressure BC at inlet
      ii = 1
      do j=1,my+2
        do k=1,mz+2
          if(.not.walls(ii,j,k)) then
            rho(ii,j,k) = rhol
            uy(ii,j,k) = uy(ii+1,j,k)!0.0d0
            uz(ii,j,k) = uz(ii+1,j,k) !0.0d0            
            ux(ii,j,k) = 1.0d0 - 1.0d0/rhol * (fi(ii,j,k,3) + fi(ii,j,k,4) + fi(ii,j,k,5) + fi(ii,j,k,6) &
                        + fi(ii,j,k,15) + fi(ii,j,k,16) + fi(ii,j,k,17) + fi(ii,j,k,18) + fi(ii,j,k,0) &
                        + 2.0d0*(fi(ii,j,k,2) + fi(ii,j,k,8) + fi(ii,j,k,10) + fi(ii,j,k,12) + fi(ii,j,k,14)))

            !common part
            nyx_inj = 0.5d0 * (fi(ii,j,k,3) + fi(ii,j,k,15) + fi(ii,j,k,17) &
                              - fi(ii,j,k,4) - fi(ii,j,k,16) - fi(ii,j,k,18)) - 1.0d0/3.0d0*rho(ii,j,k)*uy(ii,j,k)
            nzx_inj = 0.5d0 * (fi(ii,j,k,5) + fi(ii,j,k,16) + fi(ii,j,k,15) &
                              - fi(ii,j,k,6) - fi(ii,j,k,17) - fi(ii,j,k,18)) - 1.0d0/3.0d0*rho(ii,j,k)*uz(ii,j,k)
            
            fi(ii,j,k,1) = fi(ii,j,k,2) + rho(ii,j,k)/3.0d0*ux(ii,j,k)
            fi(ii,j,k,9) = fi(ii,j,k,8) + rho(ii,j,k)/6.0d0*(ux(ii,j,k) - uy(ii,j,k)) + nyx_inj
            fi(ii,j,k,7) = fi(ii,j,k,10) + rho(ii,j,k)/6.0d0*(ux(ii,j,k) + uy(ii,j,k)) - nyx_inj
            fi(ii,j,k,11) = fi(ii,j,k,14) + rho(ii,j,k)/6.0d0*(ux(ii,j,k) + uz(ii,j,k)) - nzx_inj
            fi(ii,j,k,13) = fi(ii,j,k,12) + rho(ii,j,k)/6.0d0*(ux(ii,j,k) - uz(ii,j,k)) + nzx_inj
          else
            ux(ii,j,k) = 0.0d0
 
          endif
        enddo
      enddo



      ii = mx
      do j=1,my+2
        do k=1,mz+2
          if (.not.walls(ii,j,k)) then
            rho(ii,j,k) = rhor
            uy(ii,j,k) = uy(ii-1,j,k)
            uz(ii,j,k) = uz(ii-1,j,k)
            ux(ii,j,k) = - 1.0d0 + 1.0d0/rho(ii,j,k) * (fi(ii,j,k,3) + fi(ii,j,k,4) + fi(ii,j,k,5) + fi(ii,j,k,6) &
                        + fi(ii,j,k,15) + fi(ii,j,k,16) + fi(ii,j,k,17) + fi(ii,j,k,18) + fi(ii,j,k,0) &
                        + 2.0d0*(fi(ii,j,k,1) + fi(ii,j,k,9) + fi(ii,j,k,7) + fi(ii,j,k,11) + fi(ii,j,k,13)))
            nyx_out = 0.5d0 * (fi(ii,j,k,3) + fi(ii,j,k,15) + fi(ii,j,k,17) &
                              - fi(ii,j,k,4) - fi(ii,j,k,16) - fi(ii,j,k,18)) - 1.0d0/3.0d0*rho(ii,j,k)*uy(ii,j,k)
            nzx_out = 0.5d0 * (fi(ii,j,k,5) + fi(ii,j,k,16) + fi(ii,j,k,15) &
                              - fi(ii,j,k,6) - fi(ii,j,k,17) - fi(ii,j,k,18)) - 1.0d0/3.0d0*rho(ii,j,k)*uz(ii,j,k)
            fi(ii,j,k,2) = fi(ii,j,k,1) + rho(ii,j,k)/3.0d0*(-ux(ii,j,k))
            fi(ii,j,k,8) = fi(ii,j,k,9) + rho(ii,j,k)/6.0d0*(-ux(ii,j,k) + uy(ii,j,k)) - nyx_out
            fi(ii,j,k,10) = fi(ii,j,k,7) + rho(ii,j,k)/6.0d0*(-ux(ii,j,k) - uy(ii,j,k)) + nyx_out
            fi(ii,j,k,14) = fi(ii,j,k,11) + rho(ii,j,k)/6.0d0*(-ux(ii,j,k) - uz(ii,j,k)) + nzx_out
            fi(ii,j,k,12) = fi(ii,j,k,13) + rho(ii,j,k)/6.0d0*(-ux(ii,j,k) + uz(ii,j,k)) - nzx_out
          else
            ux(ii,j,k) = 0.0d0
            uy(ii,j,k) = 0.0d0
            uz(ii,j,k) = 0.0d0  
          endif
        enddo
      enddo


      end subroutine boundaryf


      subroutine boundaryf_z(fi,rho,ux,uy,uz,rhol,rhor,uin,walls,npz)
      
      use comvar_module
      use mpi
      
      implicit none
      
      integer :: i,j,k,kk,s,npz
      real*8 :: nyz_inj,nxz_inj,nyz_out,nxz_out,rhol,rhor,uin

      real*8, dimension(1:mx,1:my+2,1:mz+2,0:18) :: fi
      real*8, dimension(1:mx,1:my+2,1:mz+2) :: rho,ux,uy,uz
      logical, dimension(1:mx,1:my+2,1:mz+2) :: walls  

      if(mpicoords(3)==0) then
      kk = 2 
      do j=1,my+2
        do i=1,mx
          if(.not.walls(i,j,kk)) then
            rho(i,j,kk) = rhol
            uy(i,j,kk) = uy(i,j,kk+1)!0.0d0
            ux(i,j,kk) = ux(i,j,kk+1) !0.0d0            
            uz(i,j,kk) = 1.0d0 - 1.0d0/rhol * (fi(i,j,kk,1) + fi(i,j,kk,2) + fi(i,j,kk,3) + fi(i,j,kk,4) &
                        + fi(i,j,kk,7) + fi(i,j,kk,8) + fi(i,j,kk,9) + fi(i,j,kk,10) + fi(i,j,kk,0) &
                        + 2.0d0*(fi(i,j,kk,6) + fi(i,j,kk,13) + fi(i,j,kk,14) + fi(i,j,kk,17) + fi(i,j,kk,18)))

            !common part
            nyz_inj = 0.5d0 * (fi(i,j,kk,3) + fi(i,j,kk,7) + fi(i,j,kk,8) &
                              - fi(i,j,kk,4) - fi(i,j,kk,9) - fi(i,j,kk,10)) - 1.0d0/3.0d0*rho(i,j,kk)*uy(i,j,kk)
            nxz_inj = 0.5d0 * (fi(i,j,kk,1) + fi(i,j,kk,7) + fi(i,j,kk,9) &
                              - fi(i,j,kk,2) - fi(i,j,kk,8) - fi(i,j,kk,10)) - 1.0d0/3.0d0*rho(i,j,kk)*ux(i,j,kk)
            
            fi(i,j,kk,5) = fi(i,j,kk,6) + rho(i,j,kk)/3.0d0*uz(i,j,kk)
            fi(i,j,kk,11) = fi(i,j,kk,14) + rho(i,j,kk)/6.0d0*(uz(i,j,kk) + ux(i,j,kk)) - nxz_inj
            fi(i,j,kk,12) = fi(i,j,kk,13) + rho(i,j,kk)/6.0d0*(uz(i,j,kk) - ux(i,j,kk)) + nxz_inj
            fi(i,j,kk,15) = fi(i,j,kk,18) + rho(i,j,kk)/6.0d0*(uz(i,j,kk) + uy(i,j,kk)) - nyz_inj
            fi(i,j,kk,16) = fi(i,j,kk,17) + rho(i,j,kk)/6.0d0*(uz(i,j,kk) - uy(i,j,kk)) + nyz_inj
          else
            ux(i,j,kk) = 0.0d0
 
          endif
        enddo
      enddo
      endif


      if(mpicoords(3)==npz-1) then
      kk = mz+1
      do j=1,my+2
        do i=1,mx
          if (.not.walls(i,j,kk)) then
            rho(i,j,kk) = rhor
            uy(i,j,kk) = uy(i,j,kk-1)
            ux(i,j,kk) = ux(i,j,kk-1)
            uz(i,j,kk) = - 1.0d0 + 1.0d0/rho(i,j,kk) * (fi(i,j,kk,1) + fi(i,j,kk,2) + fi(i,j,kk,3) + fi(i,j,kk,4) &
                        + fi(i,j,kk,7) + fi(i,j,kk,8) + fi(i,j,kk,9) + fi(i,j,kk,10) + fi(i,j,kk,0) &
                        + 2.0d0*(fi(i,j,kk,5) + fi(i,j,kk,11) + fi(i,j,kk,12) + fi(i,j,kk,15) + fi(i,j,kk,16)))

            nyz_out = 0.5d0 * (fi(i,j,kk,3) + fi(i,j,kk,7) + fi(i,j,kk,8) &
                              - fi(i,j,kk,4) - fi(i,j,kk,9) - fi(i,j,kk,10)) - 1.0d0/3.0d0*rho(i,j,kk)*uy(i,j,kk)
            nxz_out = 0.5d0 * (fi(i,j,kk,1) + fi(i,j,kk,7) + fi(i,j,kk,9) &
                              - fi(i,j,kk,2) - fi(i,j,kk,8) - fi(i,j,kk,10)) - 1.0d0/3.0d0*rho(i,j,kk)*ux(i,j,kk)

            fi(i,j,kk,6) = fi(i,j,kk,5) + rho(i,j,kk)/3.0d0*(-uz(i,j,kk))
            fi(i,j,kk,14) = fi(i,j,kk,11) + rho(i,j,kk)/6.0d0*(-uz(i,j,kk) - ux(i,j,kk)) + nxz_out
            fi(i,j,kk,13) = fi(i,j,kk,12) + rho(i,j,kk)/6.0d0*(-uz(i,j,kk) + ux(i,j,kk)) - nxz_out
            fi(i,j,kk,18) = fi(i,j,kk,15) + rho(i,j,kk)/6.0d0*(-uz(i,j,kk) - uy(i,j,kk)) + nyz_out
            fi(i,j,kk,17) = fi(i,j,kk,16) + rho(i,j,kk)/6.0d0*(-uz(i,j,kk) + uy(i,j,kk)) - nyz_out
          else
            ux(i,j,kk) = 0.0d0
            uy(i,j,kk) = 0.0d0
            uz(i,j,kk) = 0.0d0  
          endif
        enddo
      enddo
      endif

      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      end subroutine boundaryf_z

      subroutine boundaryf_z_incompressible(fi,rho,ux,uy,uz,rhol,rhor,uin,walls,npz)
      use comvar_module
      use mpi

      implicit none

      integer :: i,j,k,kk,s,npz
      real*8 :: nyz_inj,nxz_inj,nyz_out,nxz_out,rhol,rhor,uin

      real*8, dimension(1:mx,1:my+2,1:mz+2,0:18) :: fi
      real*8, dimension(1:mx,1:my+2,1:mz+2) :: rho,ux,uy,uz
      logical, dimension(1:mx,1:my+2,1:mz+2) :: walls

      if(mpicoords(3)==0) then
      kk = 2
      do j=1,my+2
        do i=1,mx
          if(.not.walls(i,j,kk)) then
            rho(i,j,kk) = rhol
            uy(i,j,kk) = uy(i,j,kk+1)!0.0d0
            ux(i,j,kk) = ux(i,j,kk+1) !0.0d0            
            uz(i,j,kk) = rhol - (fi(i,j,kk,1) + fi(i,j,kk,2) + fi(i,j,kk,3) + fi(i,j,kk,4) &
                        + fi(i,j,kk,7) + fi(i,j,kk,8) + fi(i,j,kk,9) + fi(i,j,kk,10) + fi(i,j,kk,0) &
                        + 2.0d0*(fi(i,j,kk,6) + fi(i,j,kk,13) + fi(i,j,kk,14) + fi(i,j,kk,17) + fi(i,j,kk,18)))

            !common part
            nyz_inj = 0.5d0 * (fi(i,j,kk,3) + fi(i,j,kk,7) + fi(i,j,kk,8) &
                              - fi(i,j,kk,4) - fi(i,j,kk,9) - fi(i,j,kk,10)) - 1.0d0/3.0d0*uy(i,j,kk)
            nxz_inj = 0.5d0 * (fi(i,j,kk,1) + fi(i,j,kk,7) + fi(i,j,kk,9) &
                              - fi(i,j,kk,2) - fi(i,j,kk,8) - fi(i,j,kk,10)) - 1.0d0/3.0d0*ux(i,j,kk)

            fi(i,j,kk,5) = fi(i,j,kk,6) + 1.d0/3.0d0*uz(i,j,kk)
            fi(i,j,kk,11) = fi(i,j,kk,14) + 1.d0/6.0d0*(uz(i,j,kk) + ux(i,j,kk)) - nxz_inj
            fi(i,j,kk,12) = fi(i,j,kk,13) + 1.d0/6.0d0*(uz(i,j,kk) - ux(i,j,kk)) + nxz_inj
            fi(i,j,kk,15) = fi(i,j,kk,18) + 1.d0/6.0d0*(uz(i,j,kk) + uy(i,j,kk)) - nyz_inj
            fi(i,j,kk,16) = fi(i,j,kk,17) + 1.d0/6.0d0*(uz(i,j,kk) - uy(i,j,kk)) + nyz_inj
          else
            uy(i,j,kk) = 0.0d0
            ux(i,j,kk) = 0.0d0                  
            uz(i,j,kk) = 0.0d0

          endif
        enddo
      enddo
      endif


      if(mpicoords(3)==npz-1) then
      kk = mz+1
      do j=1,my+2
        do i=1,mx
          if (.not.walls(i,j,kk)) then
            rho(i,j,kk) = rhor
            uy(i,j,kk) = uy(i,j,kk-1)
            ux(i,j,kk) = ux(i,j,kk-1)
            uz(i,j,kk) = - rhor +  (fi(i,j,kk,1) + fi(i,j,kk,2) + fi(i,j,kk,3) + fi(i,j,kk,4) &
                        + fi(i,j,kk,7) + fi(i,j,kk,8) + fi(i,j,kk,9) + fi(i,j,kk,10) + fi(i,j,kk,0) &
                        + 2.0d0*(fi(i,j,kk,5) + fi(i,j,kk,11) + fi(i,j,kk,12) + fi(i,j,kk,15) + fi(i,j,kk,16)))

            nyz_out = 0.5d0 * (fi(i,j,kk,3) + fi(i,j,kk,7) + fi(i,j,kk,8) &
                              - fi(i,j,kk,4) - fi(i,j,kk,9) - fi(i,j,kk,10)) - 1.0d0/3.0d0*uy(i,j,kk)
            nxz_out = 0.5d0 * (fi(i,j,kk,1) + fi(i,j,kk,7) + fi(i,j,kk,9) &
                              - fi(i,j,kk,2) - fi(i,j,kk,8) - fi(i,j,kk,10)) - 1.0d0/3.0d0*ux(i,j,kk)

            fi(i,j,kk,6) = fi(i,j,kk,5) + 1.d0/3.0d0*(-uz(i,j,kk))
            fi(i,j,kk,14) = fi(i,j,kk,11) + 1.d0/6.0d0*(-uz(i,j,kk) - ux(i,j,kk)) + nxz_out
            fi(i,j,kk,13) = fi(i,j,kk,12) + 1.d0/6.0d0*(-uz(i,j,kk) + ux(i,j,kk)) - nxz_out
            fi(i,j,kk,18) = fi(i,j,kk,15) + 1.d0/6.0d0*(-uz(i,j,kk) - uy(i,j,kk)) + nyz_out
            fi(i,j,kk,17) = fi(i,j,kk,16) + 1.d0/6.0d0*(-uz(i,j,kk) + uy(i,j,kk)) - nyz_out
          else
            ux(i,j,kk) = 0.0d0
            uy(i,j,kk) = 0.0d0
            uz(i,j,kk) = 0.0d0
          endif
        enddo
      enddo
      endif

      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      end subroutine boundaryf_z_incompressible    
end module boundaryf_module
