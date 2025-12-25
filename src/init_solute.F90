module initialize_solute_module

  contains

      subroutine initialize_solute(cj,walls,ncomp)
      

      use comvar_module
      
      implicit none
      
      integer :: i,j,k,p,ncomp
      real*8, dimension(1:ncomp,mx,my+2,mz+2) :: cj
      logical, dimension(mx,0:my+3,0:mz+3) :: walls



      do i=1,mx
        do j=1,my+2
          do k = 1,mz+2
            !if(walls(i,j,k)==.false.) then
              cj(:,i,j,k)=cjini(1:ncomp)
            !else
            !  cj(:,i,j,k)=1.d-20
            !endif
          enddo
        enddo
      enddo


      ! cj(1,25:26,40:41,10:11)=1.0d0
      ! if (mpicoords(2)==1.and.mpicoords(3)==1) then
      !   cj(1,25:26,20:21,5:6)=1.0d0
      ! endif

      ! if ((mpicoords(2)>=5.and.mpicoords(2)<=19).and.(mpicoords(3)>=5.and.mpicoords(3)<=19)) then
      !   do i=1,mx
      !     do j=1,my+2
      !       do k = 1,mz+2
      !         cj(:,i,j,k)=cjin(1:m)
      !       enddo
      !     enddo
      !   enddo
      ! else
      !   do i=1,mx
      !     do j=1,my+2
      !       do k = 1,mz+2
      !         cj(:,i,j,k)=cjini(1:m)
      !       enddo
      !     enddo
      !   enddo
      ! endif
      
      return
      end subroutine initialize_solute





      ! subroutine initialize_solute_restart(cj,ci,psi_lbm)
      
      ! use equilg_module
      ! use con_module

      ! use ptran_global_module
      
      ! use comvar_module
      
      ! implicit none
      
      ! integer :: i,j,k,p
   
      ! real*8, dimension(1:m,mx,my,mz+2) :: cj,psi_lbm    
      ! real*8, dimension(1:n,mx,my,mz+2) :: ci      

      

      ! call con_ci(ci,cj)

      ! psi_lbm=cj
      ! do j=1,m
      ! 	do i=1,n
      ! 	  psi_lbm(j,:,:,:)=psi_lbm(j,:,:,:)+shom(j,i)*ci(i,:,:,:)
      ! 	enddo
      ! enddo
      
      ! return
      ! end subroutine initialize_solute_restart

end module initialize_solute_module
