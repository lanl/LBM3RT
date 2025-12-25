module setupMPI_module

contains
  subroutine setupMPI(nkin,ncomp,commsize,myrank,npy,npz)
  use comvar_module
  use mpi
  implicit none
  integer :: myrank,commsize,npy,npz
  integer :: i,j,nkin,ncomp


  rreorder=.false.
!
  Tperiodic(1) = .true.
  Tperiodic(2) = .true.
  Tperiodic(3) = .true.
!
  prgrid(1) = 1
  prgrid(2) = npy
  prgrid(3) = npz


  call MPI_cart_create(MPI_COMM_WORLD, mpid, prgrid, &
                          Tperiodic,rreorder,lbecomm,ierr)

  call MPI_comm_rank(lbecomm, myrank, ierr)

  call MPI_cart_coords(lbecomm, myrank, mpid, &
                        mpicoords, ierr)



  call mpi_barrier(lbecomm,ierr)

  mx = lx
  my = ly / npy
  mz = lz / npz
  ny_local(:) = ly / npy
  nz_local(:) = lz / npz
  my_local(:) = ly / npy
  mz_local(:) = lz / npz
  zrem = mod(lz,npz)
  yrem = mod(ly,npy)

  if(zrem.ne.0) then
    do i = 1, npz
      if((i<=zrem)) then
        nz_local(i) = nz_local(i) + 1
      endif
      mz_local(i) = nz_local(i)
    enddo
  endif


   do i = npz, 2, -1
      do j = 1,i-1
         nz_local(i) = nz_local(i) + nz_local(j) !nz_local: end-point position
      enddo
   enddo

  if(yrem.ne.0) then
    do i = 1, npy
      if((i<=yrem)) then
        ny_local(i) = ny_local(i) + 1
      endif
      my_local(i) = ny_local(i)
    enddo
  endif

   do i = npy, 2, -1 
      do j = 1,i-1
         ny_local(i) = ny_local(i) + ny_local(j) !ny_local: end-point position
      enddo
   enddo


  if(zrem.ne.0) then
    do i = 1, npz
      if((i<=zrem).and.(mpicoords(3)==i-1)) then
        mz = mz + 1
      endif
    enddo
  endif

  if(yrem.ne.0) then
    do i = 1, npy
      if((i<=yrem).and.(mpicoords(2)==i-1)) then
        my = my + 1
      endif
    enddo
  endif



  call MPI_type_contiguous((my+2)*(mx), MPI_DOUBLE_PRECISION, xyplane_db1, ierr)
  call MPI_type_commit(xyplane_db1,ierr)!fi

  call mpi_type_vector(mz+2, mx, (my + 2)*mx, MPI_DOUBLE_PRECISION, xzplane_db1, ierr)
  call mpi_type_commit(xzplane_db1, ierr)

  call MPI_type_contiguous((my+2)*(mx)*ncomp, MPI_DOUBLE_PRECISION, xyplane_db2, ierr)
  call MPI_type_commit(xyplane_db2,ierr)!gi,cj

  call mpi_type_vector(mz+2, ncomp * mx, (my + 2)*mx*ncomp, MPI_DOUBLE_PRECISION, xzplane_db2, ierr)
  call mpi_type_commit(xzplane_db2, ierr)

  call MPI_type_contiguous((my+2)*(mx)*ncomp, MPI_DOUBLE_PRECISION, xyplane_db2x, ierr)
  call MPI_type_commit(xyplane_db2x,ierr)!psi

  call mpi_type_vector(mz+2, ncomp * mx, (my + 2)*mx*ncomp, MPI_DOUBLE_PRECISION, xzplane_db2x, ierr)
  call mpi_type_commit(xzplane_db2x, ierr)

  call MPI_type_contiguous((my+2)*(mx), MPI_DOUBLE_PRECISION, xyplane_db3, ierr)
  call MPI_type_commit(xyplane_db3,ierr)!3d tensors

  call mpi_type_vector(mz+2, mx, (my + 2)*mx, MPI_DOUBLE_PRECISION, xzplane_db3, ierr)
  call mpi_type_commit(xzplane_db3, ierr)

  call MPI_type_contiguous((my+4)*(mx)*nkin, MPI_DOUBLE_PRECISION, xyplane_db5, ierr)
  call MPI_type_commit(xyplane_db5,ierr)!bnd

  call mpi_type_vector(mz+4, nkin * mx,  (my + 4)*mx*nkin, MPI_DOUBLE_PRECISION, xzplane_db5, ierr)
  call mpi_type_commit(xzplane_db5, ierr)

  call mpi_type_vector(7, ncomp * mx * (my+2), (mz+2)*(my + 2)*mx*ncomp, MPI_DOUBLE_PRECISION, xyplane_db6, ierr)
  call mpi_type_commit(xyplane_db6, ierr)!gi_full

  call mpi_type_vector(7*(mz+2), ncomp * mx, (my + 2)*mx*ncomp, MPI_DOUBLE_PRECISION, xzplane_db6, ierr)
  call mpi_type_commit(xzplane_db6, ierr)


  call mpi_type_vector(19,mx * (my+2), (mz+2)*(my + 2)*mx, MPI_DOUBLE_PRECISION, xyplane_db7, ierr)
  call mpi_type_commit(xyplane_db7, ierr)!fi_full

  call mpi_type_vector(19*(mz+2), mx, (my + 2)*mx, MPI_DOUBLE_PRECISION, xzplane_db7, ierr)
  call mpi_type_commit(xzplane_db7, ierr)

  call MPI_type_contiguous((my+4)*(mx)*nkin*2, MPI_LOGICAL, xyplane_logical, ierr)
  call MPI_type_commit(xyplane_logical,ierr)!wall

  call mpi_type_vector(mz+4, nkin * mx * 2,  (my + 4)*mx*nkin, MPI_LOGICAL, xzplane_logical, ierr)
  call mpi_type_commit(xzplane_logical, ierr)





  ! Shift in each dimension to get neighbors
  call MPI_cart_shift(lbecomm, 0, 1, left(2), right(2), ierr)
  call MPI_cart_shift(lbecomm, 1, 1, rear(2), front(2), ierr)  
  call MPI_cart_shift(lbecomm, 2, 1, down(2), up(2), ierr)


  ! Determine diagonal ranks (combinations of shifts)
  call MPI_Cart_rank(lbecomm, [mpicoords(1), mpicoords(2)+1, mpicoords(3)-1], front_down, ierr)
  call MPI_Cart_rank(lbecomm, [mpicoords(1), mpicoords(2)+1, mpicoords(3)+1], front_up, ierr)
  call MPI_Cart_rank(lbecomm, [mpicoords(1), mpicoords(2)-1, mpicoords(3)-1], rear_down, ierr)
  call MPI_Cart_rank(lbecomm, [mpicoords(1), mpicoords(2)-1, mpicoords(3)+1], rear_up, ierr)


  end subroutine setupMPI




end module setupMPI_module
