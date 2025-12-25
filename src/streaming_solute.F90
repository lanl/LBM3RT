module streaming_solute_module

contains


  subroutine streamingg(gi,ncomp)
  
  use comvar_module
  use mpi
  implicit none
  

  integer :: s,ncomp
  real*8, dimension(1:ncomp,1:mx,1:my+2,1:mz+2,0:6) :: gi


    do s = 1, 6 
      if (cjxi(s) .ne. 0) then
        gi(:,:,:,:,s) = cshift(gi(:,:,:,:,s), shift=-cjxi(s), dim=2)
      endif
      if (cjyi(s) .ne. 0) then
        gi(:,:,:,:,s) = cshift(gi(:,:,:,:,s), shift=-cjyi(s), dim=3)
      endif
      if (cjzi(s) .ne. 0) then
        gi(:,:,:,:,s) = cshift(gi(:,:,:,:,s), shift=-cjzi(s), dim=4)
      endif

      call MPI_Barrier(MPI_COMM_WORLD, ierr)

    enddo

  end subroutine streamingg

end module streaming_solute_module