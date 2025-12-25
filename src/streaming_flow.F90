module streaming_flow_module

contains

  subroutine streamingf(fi,walls,myrank)
  
  use comvar_module
  use mpi
  
  implicit none
  
  integer :: s,i,j,k,myrank
  real*8, dimension(mx,my+2,mz+2,0:18) :: fi,fi_tmp
  logical, dimension(1:mx,1:my+2,1:mz+2) :: walls


    do s = 1, 18
      if (cixi(s) .ne. 0) then
        fi(:,:,:,s) = cshift(fi(:,:,:,s), shift=-cixi(s), dim=1)
      endif

      if (ciyi(s) .ne. 0) then
        fi(:,:,:,s) = cshift(fi(:,:,:,s), shift=-ciyi(s), dim=2)
      endif

      if (cizi(s) .ne. 0) then
        fi(:,:,:,s) = cshift(fi(:,:,:,s), shift=-cizi(s), dim=3)
      endif
    enddo


    do i=1,mx
      do j=1,my+2
        do k=1,mz+2
          if(walls(i,j,k)) then
            do s = 1,18
              fi_tmp(i,j,k,ind_ci_opp(s)) = fi(i,j,k,s)
            enddo

            do s = 1,18
              fi(i,j,k,s) = fi_tmp(i,j,k,s)
            enddo

          endif
        enddo
      enddo
    enddo


  end subroutine streamingf

end module streaming_flow_module
