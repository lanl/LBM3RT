module initialize_flow_module

contains

      subroutine initialize_flow(fi,rho_lbm,ux,uy,uz,walls,wall,nkin,rhol,rhor,myrank,npy,npz)
      
      use comvar_module
      
      implicit none
      
      integer :: i,j,k
      real*8 :: rhol,rhor,rhol_mpi,rhor_mpi
      integer :: myrank,npy,npz,nkin
      real*8, dimension(mx,my+2,mz+2,0:18) :: fi
      real*8, dimension(mx,my+2,mz+2) :: rho_lbm,ux,uy,uz
      logical, dimension(mx,0:my+3,0:mz+3) :: walls
      logical, dimension(nkin,mx,0:my+3,0:mz+3) :: wall
      
      fi=0.0d0

      ! do i=1,mx
      !   rho(i,:,:)=rhol-(rhol-rhor)*(i-1)/(mx-1)
      ! enddo
      
      rhor_mpi = (rhol-rhor)/npz*(npz-mpicoords(3)-1) + rhor
      rhol_mpi = rhor_mpi + (rhol-rhor)/npz
      do k = 1,mz+2
        do j = 1,my+2
          do i = 1,mx
            if(.not.walls(i,j,k)) then
              rho_lbm(i,j,k) = rhol_mpi-(rhol_mpi-rhor_mpi)*(k-1)/(mz+1)
            endif
          enddo
        enddo
      enddo

      ux=0.d0
      uy=0.d0
      uz=0.d0
      
      walls=.false. 
      wall = .false.
      !wall(2,1,:,:) = .true.
      !wall(2,mx,:,:) = .true.
      ! if(myrank==0) then
      !   wall(2,:,12:41,12:41)=.true.
      ! endif

      ! if (mpicoords(2)==2.and.(mpicoords(3)==4.or.mpicoords(3)==5)) then
      !   wall(2,21:30,:,:)=.true.
      ! endif
      ! if (mpicoords(2)==2.and.mpicoords(3)==2) then
      !   wall(2,21:30,:,:)=.true.
      ! endif
      ! if ((mpicoords(2)>=3.and.mpicoords(2)<=6).and.(mpicoords(3)>=3.and.mpicoords(3)<=6)) then
      !   wall(2,16:35,:,:)=.true.
      ! endif
      ! if ((mpicoords(2)>=4.and.mpicoords(2)<=5).and.(mpicoords(3)>=10.and.mpicoords(3)<=14)) then
      !   wall(2,21:30,:,:)=.true.
      ! endif
      ! if ((mpicoords(2)>=10.and.mpicoords(2)<=14).and.(mpicoords(3)>=10.and.mpicoords(3)<=14)) then
      !   wall(2,21:30,:,:)=.true.
      ! endif
      ! if ((mpicoords(2)>=4.and.mpicoords(2)<=5).and.(mpicoords(3)>=20.and.mpicoords(3)<=29)) then
      !   wall(2,21:30,:,:)=.true.
      ! endif
      ! if (mpicoords(2)==0) then
      !   wall(2,21:30,:,:)=.true.
      ! endif
      !if ((mpicoords(2)>=5.and.mpicoords(2)<=19).and.(mpicoords(3)>=5.and.mpicoords(3)<=19)) then
      !  wall(2,11:40,:,:)=.true.
      !endif
      ! if ((mpicoords(2)>=11.and.mpicoords(2)<=13).and.(mpicoords(3)>=11.and.mpicoords(3)<=13)) then
      !   wall(2,23:28,:,:)=.true.
      ! endif      
      !if ((mpicoords(2)>=1.and.mpicoords(2)<=3).and.(mpicoords(3)>=1.and.mpicoords(3)<=3)) then
      !  wall(2,11:40,:,:)=.true.
      !endif
      ! wall(2,7:14,8:15,8:15) = .true.
      ! if ((mpicoords(2)>=5.and.mpicoords(2)<=19).and.(mpicoords(3)>=5.and.mpicoords(3)<=19)) then
      !   wall(2,11:40,:,:)=.true.
      ! endif
      ! if (mpicoords(2)==2.and.mpicoords(3)==2) then
      !   wall(2,9:12,:,:)=.true.
      ! endif
      ! if (myrank==0) then
      !   write(*,'(''average density on rank 0 = '',1pe14.6)') sum(fi)/(mx*(my+2)*(mz+2))
      ! endif




      end subroutine initialize_flow
      
end module initialize_flow_module

