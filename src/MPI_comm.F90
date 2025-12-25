module mpi_comm_module

  contains
  subroutine mpi_fi(fi)

  use comvar_module
  use mpi
  !mpi parameters
  integer      :: status_up(MPI_STATUS_SIZE)
  integer      :: reqs_up(5)
  integer      :: status_down(MPI_STATUS_SIZE)
  integer      :: reqs_down(5)
  integer      :: status_front(MPI_STATUS_SIZE)
  integer      :: reqs_front(5)
  integer      :: status_rear(MPI_STATUS_SIZE)
  integer      :: reqs_rear(5)
  integer      :: tag
  real*8, dimension(mx,my+2,mz+2,0:18) :: fi   



  tag = 005
  call mpi_isend(fi(1,1,mz+1,5), 1, xyplane_db1, up(2), tag, lbecomm, reqs_up(1), ierr)

  tag = 011
  call mpi_isend(fi(1,1,mz+1,11), 1, xyplane_db1, up(2), tag, lbecomm, reqs_up(2), ierr)

  tag = 012
  call mpi_isend(fi(1,1,mz+1,12), 1, xyplane_db1, up(2), tag, lbecomm, reqs_up(3), ierr)

  tag = 015
  call mpi_isend(fi(1,1,mz+1,15), 1, xyplane_db1, up(2), tag, lbecomm, reqs_up(4), ierr)

  tag = 016
  call mpi_isend(fi(1,1,mz+1,16), 1, xyplane_db1, up(2), tag, lbecomm, reqs_up(5), ierr)


  tag = 006
  call mpi_isend(fi(1,1,2,6), 1, xyplane_db1, down(2), tag, lbecomm, reqs_down(1), ierr)

  tag = 013
  call mpi_isend(fi(1,1,2,13), 1, xyplane_db1, down(2), tag, lbecomm, reqs_down(2), ierr)

  tag = 014
  call mpi_isend(fi(1,1,2,14), 1, xyplane_db1, down(2), tag, lbecomm, reqs_down(3), ierr)

  tag = 017
  call mpi_isend(fi(1,1,2,17), 1, xyplane_db1, down(2), tag, lbecomm, reqs_down(4), ierr)

  tag = 018
  call mpi_isend(fi(1,1,2,18), 1, xyplane_db1, down(2), tag, lbecomm, reqs_down(5), ierr)



  tag = 005
  call mpi_recv(fi(1,1,1,5), 1, xyplane_db1, down(2), tag, lbecomm, status_up, ierr)

  tag = 011
  call mpi_recv(fi(1,1,1,11), 1, xyplane_db1, down(2), tag, lbecomm, status_up, ierr)

  tag = 012
  call mpi_recv(fi(1,1,1,12), 1, xyplane_db1, down(2), tag, lbecomm, status_up, ierr)

  tag = 015
  call mpi_recv(fi(1,1,1,15), 1, xyplane_db1, down(2), tag, lbecomm, status_up, ierr)

  tag = 016
  call mpi_recv(fi(1,1,1,16), 1, xyplane_db1, down(2), tag, lbecomm, status_up, ierr)


  tag = 006
  call mpi_recv(fi(1,1,mz+2,6), 1, xyplane_db1, up(2), tag, lbecomm, status_down, ierr)

  tag = 013
  call mpi_recv(fi(1,1,mz+2,13), 1, xyplane_db1, up(2), tag, lbecomm, status_down, ierr)

  tag = 014
  call mpi_recv(fi(1,1,mz+2,14), 1, xyplane_db1, up(2), tag, lbecomm, status_down, ierr)

  tag = 017
  call mpi_recv(fi(1,1,mz+2,17), 1, xyplane_db1, up(2), tag, lbecomm, status_down, ierr)

  tag = 018
  call mpi_recv(fi(1,1,mz+2,18), 1, xyplane_db1, up(2), tag, lbecomm, status_down, ierr)


  !syncronize
  call mpi_wait(reqs_up(1), status_up, ierr)
  call mpi_wait(reqs_down(1), status_down, ierr)
  call mpi_wait(reqs_up(2), status_up, ierr)
  call mpi_wait(reqs_down(2), status_down, ierr)
  call mpi_wait(reqs_up(3), status_up, ierr)
  call mpi_wait(reqs_down(3), status_down, ierr)
  call mpi_wait(reqs_up(4), status_up, ierr)
  call mpi_wait(reqs_down(4), status_down, ierr)
  call mpi_wait(reqs_up(5), status_up, ierr)
  call mpi_wait(reqs_down(5), status_down, ierr)

  tag = 503
  call mpi_isend(fi(1,my+1,1,3), 1, xzplane_db1, front(2), tag, lbecomm, reqs_front(1), ierr)

  tag = 507
  call mpi_isend(fi(1,my+1,1,7), 1, xzplane_db1, front(2), tag, lbecomm, reqs_front(2), ierr)

  tag = 508
  call mpi_isend(fi(1,my+1,1,8), 1, xzplane_db1, front(2), tag, lbecomm, reqs_front(3), ierr)

  tag = 515
  call mpi_isend(fi(1,my+1,1,15), 1, xzplane_db1, front(2), tag, lbecomm, reqs_front(4), ierr)

  tag = 517
  call mpi_isend(fi(1,my+1,1,17), 1, xzplane_db1, front(2), tag, lbecomm, reqs_front(5), ierr)


  tag = 504
  call mpi_isend(fi(1,2,1,4), 1, xzplane_db1, rear(2), tag, lbecomm, reqs_rear(1), ierr)

  tag = 509
  call mpi_isend(fi(1,2,1,9), 1, xzplane_db1, rear(2), tag, lbecomm, reqs_rear(2), ierr)

  tag = 510
  call mpi_isend(fi(1,2,1,10), 1, xzplane_db1, rear(2), tag, lbecomm, reqs_rear(3), ierr)

  tag = 516
  call mpi_isend(fi(1,2,1,16), 1, xzplane_db1, rear(2), tag, lbecomm, reqs_rear(4), ierr)

  tag = 518
  call mpi_isend(fi(1,2,1,18), 1, xzplane_db1, rear(2), tag, lbecomm, reqs_rear(5), ierr)


  tag = 503
  call mpi_recv(fi(1,1,1,3), 1, xzplane_db1, rear(2), tag, lbecomm, status_front, ierr)

  tag = 507
  call mpi_recv(fi(1,1,1,7), 1, xzplane_db1, rear(2), tag, lbecomm, status_front, ierr)

  tag = 508
  call mpi_recv(fi(1,1,1,8), 1, xzplane_db1, rear(2), tag, lbecomm, status_front, ierr)

  tag = 515
  call mpi_recv(fi(1,1,1,15), 1, xzplane_db1, rear(2), tag, lbecomm, status_front, ierr)

  tag = 517
  call mpi_recv(fi(1,1,1,17), 1, xzplane_db1, rear(2), tag, lbecomm, status_front, ierr)


  tag = 504
  call mpi_recv(fi(1,my+2,1,4), 1, xzplane_db1, front(2), tag, lbecomm, status_rear, ierr)

  tag = 509
  call mpi_recv(fi(1,my+2,1,9), 1, xzplane_db1, front(2), tag, lbecomm, status_rear, ierr)

  tag = 510
  call mpi_recv(fi(1,my+2,1,10), 1, xzplane_db1, front(2), tag, lbecomm, status_rear, ierr)

  tag = 516
  call mpi_recv(fi(1,my+2,1,16), 1, xzplane_db1, front(2), tag, lbecomm, status_rear, ierr)

  tag = 518
  call mpi_recv(fi(1,my+2,1,18), 1, xzplane_db1, front(2), tag, lbecomm, status_rear, ierr)


  !syncronize
  call mpi_wait(reqs_front(1), status_front, ierr)
  call mpi_wait(reqs_rear(1), status_rear, ierr)
  call mpi_wait(reqs_front(2), status_front, ierr)
  call mpi_wait(reqs_rear(2), status_rear, ierr)
  call mpi_wait(reqs_front(3), status_front, ierr)
  call mpi_wait(reqs_rear(3), status_rear, ierr)
  call mpi_wait(reqs_front(4), status_front, ierr)
  call mpi_wait(reqs_rear(4), status_rear, ierr)
  call mpi_wait(reqs_front(5), status_front, ierr)
  call mpi_wait(reqs_rear(5), status_rear, ierr)


  end subroutine mpi_fi


  subroutine mpi_gi(gi,ncomp)

  use comvar_module
  use mpi
  !mpi parameters
  integer      :: status_up(MPI_STATUS_SIZE)
  integer      :: reqs_up(1)
  integer      :: status_down(MPI_STATUS_SIZE)
  integer      :: reqs_down(1)
  integer      :: status_front(MPI_STATUS_SIZE)
  integer      :: reqs_front(1)
  integer      :: status_rear(MPI_STATUS_SIZE)
  integer      :: reqs_rear(1)
  integer      :: tag,ncomp
  real*8, dimension(1:ncomp,mx,my+2,mz+2,0:6) :: gi   



  tag = 055
  call mpi_isend(gi(1,1,1,mz+1,5), 1, xyplane_db2, up(2), tag, lbecomm, reqs_up(1), ierr)

  tag = 056
  call mpi_isend(gi(1,1,1,2,6), 1, xyplane_db2, down(2), tag, lbecomm, reqs_down(1), ierr)

  tag = 055
  call mpi_recv(gi(1,1,1,1,5), 1, xyplane_db2, down(2), tag, lbecomm, status_up, ierr)

  tag = 056
  call mpi_recv(gi(1,1,1,mz+2,6), 1, xyplane_db2, up(2), tag, lbecomm, status_down, ierr)

  !syncronize
  call mpi_wait(reqs_up(1), status_up, ierr)
  call mpi_wait(reqs_down(1), status_down, ierr)  

  tag = 053
  call mpi_isend(gi(1,1,my+1,1,3), 1, xzplane_db2, front(2), tag, lbecomm, reqs_front(1), ierr)

  tag = 054
  call mpi_isend(gi(1,1,2,1,4), 1, xzplane_db2, rear(2), tag, lbecomm, reqs_rear(1), ierr)

  tag = 053
  call mpi_recv(gi(1,1,1,1,3), 1, xzplane_db2, rear(2), tag, lbecomm, status_front, ierr)

  tag = 054
  call mpi_recv(gi(1,1,my+2,1,4), 1, xzplane_db2, front(2), tag, lbecomm, status_rear, ierr)

  !syncronize
  call mpi_wait(reqs_front(1), status_front, ierr)
  call mpi_wait(reqs_rear(1), status_rear, ierr)

  end subroutine mpi_gi


 


  subroutine mpi_fi_full(fi)

  use comvar_module
  use mpi
  !mpi parameters
  integer      :: status_up(MPI_STATUS_SIZE)
  integer      :: reqs_up(1)
  integer      :: status_down(MPI_STATUS_SIZE)
  integer      :: reqs_down(1)
  integer      :: status_front(MPI_STATUS_SIZE)
  integer      :: reqs_front(1)
  integer      :: status_rear(MPI_STATUS_SIZE)
  integer      :: reqs_rear(1)
  integer      :: tag
  real*8, dimension(mx,my+2,mz+2,0:18) :: fi   



  tag = 530
  call mpi_isend(fi(1,1,mz+1,0), 1, xyplane_db7, up(2), tag, lbecomm, reqs_up(1), ierr)

  tag = 531
  call mpi_isend(fi(1,1,2,0), 1, xyplane_db7, down(2), tag, lbecomm, reqs_down(1), ierr)

  tag = 530
  call mpi_recv(fi(1,1,1,0), 1, xyplane_db7, down(2), tag, lbecomm, status_up, ierr)

  tag = 531
  call mpi_recv(fi(1,1,mz+2,0), 1, xyplane_db7, up(2), tag, lbecomm, status_down, ierr)

  !syncronize
  call mpi_wait(reqs_up(1), status_up, ierr)
  call mpi_wait(reqs_down(1), status_down, ierr)


  tag = 532
  call mpi_isend(fi(1,my+1,1,0), 1, xzplane_db7, front(2), tag, lbecomm, reqs_front(1), ierr)

  tag = 533
  call mpi_isend(fi(1,2,1,0), 1, xzplane_db7, rear(2), tag, lbecomm, reqs_rear(1), ierr)

  tag = 532
  call mpi_recv(fi(1,1,1,0), 1, xzplane_db7, rear(2), tag, lbecomm, status_front, ierr)

  tag = 533
  call mpi_recv(fi(1,my+2,1,0), 1, xzplane_db7, front(2), tag, lbecomm, status_rear, ierr)

  !syncronize
  call mpi_wait(reqs_front(1), status_front, ierr)
  call mpi_wait(reqs_rear(1), status_rear, ierr)

  end subroutine mpi_fi_full


  subroutine mpi_gi_full(gi,ncomp)

  use comvar_module
  use mpi
  !mpi parameters
  integer      :: status_up(MPI_STATUS_SIZE)
  integer      :: reqs_up(1)
  integer      :: status_down(MPI_STATUS_SIZE)
  integer      :: reqs_down(1)
  integer      :: status_front(MPI_STATUS_SIZE)
  integer      :: reqs_front(1)
  integer      :: status_rear(MPI_STATUS_SIZE)
  integer      :: reqs_rear(1)
  integer      :: tag,ncomp
  real*8, dimension(1:ncomp,mx,my+2,mz+2,0:6) :: gi   



  tag = 057
  call mpi_isend(gi(1,1,1,mz+1,0), 1, xyplane_db6, up(2), tag, lbecomm, reqs_up(1), ierr)

  tag = 058
  call mpi_isend(gi(1,1,1,2,0), 1, xyplane_db6, down(2), tag, lbecomm, reqs_down(1), ierr)

  tag = 057
  call mpi_recv(gi(1,1,1,1,0), 1, xyplane_db6, down(2), tag, lbecomm, status_up, ierr)

  tag = 058
  call mpi_recv(gi(1,1,1,mz+2,0), 1, xyplane_db6, up(2), tag, lbecomm, status_down, ierr)

  !syncronize
  call mpi_wait(reqs_up(1), status_up, ierr)
  call mpi_wait(reqs_down(1), status_down, ierr)


  tag = 059
  call mpi_isend(gi(1,1,my+1,1,0), 1, xzplane_db6, front(2), tag, lbecomm, reqs_front(1), ierr)

  tag = 060
  call mpi_isend(gi(1,1,2,1,0), 1, xzplane_db6, rear(2), tag, lbecomm, reqs_rear(1), ierr)

  tag = 059
  call mpi_recv(gi(1,1,1,1,0), 1, xzplane_db6, rear(2), tag, lbecomm, status_front, ierr)

  tag = 060
  call mpi_recv(gi(1,1,my+2,1,0), 1, xzplane_db6, front(2), tag, lbecomm, status_rear, ierr)

  !syncronize
  call mpi_wait(reqs_front(1), status_front, ierr)
  call mpi_wait(reqs_rear(1), status_rear, ierr)

  end subroutine mpi_gi_full


  subroutine mpi_cj(cj,ncomp)

  use comvar_module
  use mpi
  !mpi parameters
  integer      :: status_up(MPI_STATUS_SIZE)
  integer      :: reqs_up(1)
  integer      :: status_down(MPI_STATUS_SIZE)
  integer      :: reqs_down(1)
  integer      :: status_front(MPI_STATUS_SIZE)
  integer      :: reqs_front(1)
  integer      :: status_rear(MPI_STATUS_SIZE)
  integer      :: reqs_rear(1)
  integer      :: tag
  real*8, dimension(1:ncomp,mx,my+2,mz+2) :: cj



  tag = 155
  call mpi_isend(cj(1,1,1,mz+1), 1, xyplane_db2, up(2), tag, lbecomm, reqs_up(1), ierr)

  tag = 156
  call mpi_isend(cj(1,1,1,2), 1, xyplane_db2, down(2), tag, lbecomm, reqs_down(1), ierr)

  tag = 155
  call mpi_recv(cj(1,1,1,1), 1, xyplane_db2, down(2), tag, lbecomm, status_up, ierr)

  tag = 156
  call mpi_recv(cj(1,1,1,mz+2), 1, xyplane_db2, up(2), tag, lbecomm, status_down, ierr)

  !syncronize
  call mpi_wait(reqs_up(1), status_up, ierr)
  call mpi_wait(reqs_down(1), status_down, ierr)

  tag = 153
  call mpi_isend(cj(1,1,my+1,1), 1, xzplane_db2, front(2), tag, lbecomm, reqs_front(1), ierr)

  tag = 154
  call mpi_isend(cj(1,1,2,1), 1, xzplane_db2, rear(2), tag, lbecomm, reqs_rear(1), ierr)

  tag = 153
  call mpi_recv(cj(1,1,1,1), 1, xzplane_db2, rear(2), tag, lbecomm, status_front, ierr)

  tag = 154
  call mpi_recv(cj(1,1,my+2,1), 1, xzplane_db2, front(2), tag, lbecomm, status_rear, ierr)

  !syncronize
  call mpi_wait(reqs_front(1), status_front, ierr)
  call mpi_wait(reqs_rear(1), status_rear, ierr)

  end subroutine mpi_cj


  subroutine mpi_psi(psi_lbm,ncomp)

  use comvar_module
  use mpi
  !mpi parameters
  integer      :: status_up(MPI_STATUS_SIZE)
  integer      :: reqs_up(1)
  integer      :: status_down(MPI_STATUS_SIZE)
  integer      :: reqs_down(1)
  integer      :: status_front(MPI_STATUS_SIZE)
  integer      :: reqs_front(1)
  integer      :: status_rear(MPI_STATUS_SIZE)
  integer      :: reqs_rear(1)
  integer      :: tag,ncomp

  real*8, dimension(1:ncomp,mx,my+2,mz+2) :: psi_lbm



  tag = 159
  call mpi_isend(psi_lbm(1,1,1,mz+1), 1, xyplane_db2x, up(2), tag, lbecomm, reqs_up(1), ierr)

  tag = 160
  call mpi_isend(psi_lbm(1,1,1,2), 1, xyplane_db2x, down(2), tag, lbecomm, reqs_down(1), ierr)

  tag = 159
  call mpi_recv(psi_lbm(1,1,1,1), 1, xyplane_db2x, down(2), tag, lbecomm, status_up, ierr)

  tag = 160
  call mpi_recv(psi_lbm(1,1,1,mz+2), 1, xyplane_db2x, up(2), tag, lbecomm, status_down, ierr)

  !syncronize
  call mpi_wait(reqs_up(1), status_up, ierr)
  call mpi_wait(reqs_down(1), status_down, ierr)

  tag = 157
  call mpi_isend(psi_lbm(1,1,my+1,1), 1, xzplane_db2x, front(2), tag, lbecomm, reqs_front(1), ierr)

  tag = 158
  call mpi_isend(psi_lbm(1,1,2,1), 1, xzplane_db2x, rear(2), tag, lbecomm, reqs_rear(1), ierr)

  tag = 157
  call mpi_recv(psi_lbm(1,1,1,1), 1, xzplane_db2x, rear(2), tag, lbecomm, status_front, ierr)

  tag = 158
  call mpi_recv(psi_lbm(1,1,my+2,1), 1, xzplane_db2x, front(2), tag, lbecomm, status_rear, ierr)

  !syncronize
  call mpi_wait(reqs_front(1), status_front, ierr)
  call mpi_wait(reqs_rear(1), status_rear, ierr)

  end subroutine mpi_psi


  subroutine mpi_3d(tensor_3d)

  use comvar_module
  use mpi
  !mpi parameters
  integer      :: status_up(MPI_STATUS_SIZE)
  integer      :: reqs_up(1)
  integer      :: status_down(MPI_STATUS_SIZE)
  integer      :: reqs_down(1)
  integer      :: status_front(MPI_STATUS_SIZE)
  integer      :: reqs_front(1)
  integer      :: status_rear(MPI_STATUS_SIZE)
  integer      :: reqs_rear(1)
  integer      :: tag
  real*8, dimension(mx,my+2,mz+2) :: tensor_3d



  tag = 065
  call mpi_isend(tensor_3d(1,1,mz+1), 1, xyplane_db3, up(2), tag, lbecomm, reqs_up(1), ierr)

  tag = 066
  call mpi_isend(tensor_3d(1,1,2), 1, xyplane_db3, down(2), tag, lbecomm, reqs_down(1), ierr)

  tag = 065
  call mpi_recv(tensor_3d(1,1,1), 1, xyplane_db3, down(2), tag, lbecomm, status_up, ierr)

  tag = 066
  call mpi_recv(tensor_3d(1,1,mz+2), 1, xyplane_db3, up(2), tag, lbecomm, status_down, ierr)

  !syncronize
  call mpi_wait(reqs_up(1), status_up, ierr)
  call mpi_wait(reqs_down(1), status_down, ierr)

  tag = 063
  call mpi_isend(tensor_3d(1,my+1,1), 1, xzplane_db3, front(2), tag, lbecomm, reqs_front(1), ierr)

  tag = 064
  call mpi_isend(tensor_3d(1,2,1), 1, xzplane_db3, rear(2), tag, lbecomm, reqs_rear(1), ierr)

  tag = 063
  call mpi_recv(tensor_3d(1,1,1), 1, xzplane_db3, rear(2), tag, lbecomm, status_front, ierr)

  tag = 064
  call mpi_recv(tensor_3d(1,my+2,1), 1, xzplane_db3, front(2), tag, lbecomm, status_rear, ierr)

  !syncronize
  call mpi_wait(reqs_front(1), status_front, ierr)
  call mpi_wait(reqs_rear(1), status_rear, ierr)

  end subroutine mpi_3d



!   subroutine mpi_wall(wall)

!   use comvar_module
!   use mpi

!   integer      :: status_up(MPI_STATUS_SIZE)
!   integer      :: status_down(MPI_STATUS_SIZE)
!   integer      :: status_front(MPI_STATUS_SIZE)
!   integer      :: status_rear(MPI_STATUS_SIZE)
!   integer      :: tag
!   logical, dimension(l,mx,0:my+3,0:mz+3) :: wall

!   tag = 105
!   call MPI_Sendrecv(wall(1,1,0,mz), 1, xyplane_logical, up(2), tag, &
!                     wall(1,1,0,0), 1, xyplane_logical, down(2), tag, lbecomm, status_up, ierr)

!   tag = 106
!   call MPI_Sendrecv(wall(1,1,0,2), 1, xyplane_logical, down(2), tag, &
!                     wall(1,1,0,mz+2), 1, xyplane_logical, up(2), tag, lbecomm, status_down, ierr)

!   tag = 107
!   call MPI_Sendrecv(wall(1,1,my,0), 1, xzplane_logical, front(2), tag, &
!                     wall(1,1,0,0), 1, xzplane_logical, rear(2), tag, lbecomm, status_front, ierr)

!   tag = 108
!   call MPI_Sendrecv(wall(1,1,2,0), 1, xzplane_logical, rear(2), tag, &
!                     wall(1,1,my+2,0), 1, xzplane_logical, front(2), tag, lbecomm, status_rear, ierr)

! end subroutine mpi_wall

  subroutine mpi_wall(wall,nkin)

  use comvar_module
  use mpi

  integer      :: status_up(MPI_STATUS_SIZE)
  integer      :: reqs_up(1)
  integer      :: status_down(MPI_STATUS_SIZE)
  integer      :: reqs_down(1)
  integer      :: status_front(MPI_STATUS_SIZE)
  integer      :: reqs_front(1)
  integer      :: status_rear(MPI_STATUS_SIZE)
  integer      :: reqs_rear(1)
  integer      :: tag,nkin
  logical, dimension(nkin,mx,0:my+3,0:mz+3) :: wall



  tag = 105
  call mpi_isend(wall(1,1,0,mz), 1, xyplane_logical, up(2), tag, lbecomm, reqs_up(1), ierr)

  tag = 106
  call mpi_isend(wall(1,1,0,2), 1, xyplane_logical, down(2), tag, lbecomm, reqs_down(1), ierr)

  tag = 105
  call mpi_recv(wall(1,1,0,0), 1, xyplane_logical, down(2), tag, lbecomm, status_up, ierr)

  tag = 106
  call mpi_recv(wall(1,1,0,mz+2), 1, xyplane_logical, up(2), tag, lbecomm, status_down, ierr)

  call mpi_wait(reqs_up(1), status_up, ierr)  
  call mpi_wait(reqs_down(1), status_down, ierr)   


  tag = 107
  call mpi_isend(wall(1,1,my,0), 1, xzplane_logical, front(2), tag, lbecomm, reqs_front(1), ierr)

  tag = 108
  call mpi_isend(wall(1,1,2,0), 1, xzplane_logical, rear(2), tag, lbecomm, reqs_rear(1), ierr)  

  tag = 107
  call mpi_recv(wall(1,1,0,0), 1, xzplane_logical, rear(2), tag, lbecomm, status_front, ierr)

  tag = 108
  call mpi_recv(wall(1,1,my+2,0), 1, xzplane_logical, front(2), tag, lbecomm, status_rear, ierr)

  call mpi_wait(reqs_front(1), status_front, ierr) 
  call mpi_wait(reqs_rear(1), status_rear, ierr)  
 

  end subroutine mpi_wall


subroutine mpi_bnd(bnd,nkin)

  use comvar_module
  use mpi

  integer      :: status_up(MPI_STATUS_SIZE)
  integer      :: reqs_up(1)
  integer      :: status_down(MPI_STATUS_SIZE)
  integer      :: reqs_down(1)
  integer      :: status_front(MPI_STATUS_SIZE)
  integer      :: reqs_front(1)
  integer      :: status_rear(MPI_STATUS_SIZE)
  integer      :: reqs_rear(1)
  integer      :: tag,nkin
  real*8, dimension(nkin,mx,0:my+3,0:mz+3) :: bnd



  tag = 115
  call mpi_isend(bnd(1,1,0,mz+1), 1, xyplane_db5, up(2), tag, lbecomm, reqs_up(1), ierr)

  tag = 116
  call mpi_isend(bnd(1,1,0,2), 1, xyplane_db5, down(2), tag, lbecomm, reqs_down(1), ierr)

  tag = 115
  call mpi_recv(bnd(1,1,0,1), 1, xyplane_db5, down(2), tag, lbecomm, status_up, ierr)

  tag = 116
  call mpi_recv(bnd(1,1,0,mz+2), 1, xyplane_db5, up(2), tag, lbecomm, status_down, ierr)

  call mpi_wait(reqs_up(1), status_up, ierr)  
  call mpi_wait(reqs_down(1), status_down, ierr)   


  tag = 117
  call mpi_isend(bnd(1,1,my+1,0), 1, xzplane_db5, front(2), tag, lbecomm, reqs_front(1), ierr)

  tag = 118
  call mpi_isend(bnd(1,1,2,0), 1, xzplane_db5, rear(2), tag, lbecomm, reqs_rear(1), ierr) 

  tag = 117
  call mpi_recv(bnd(1,1,1,0), 1, xzplane_db5, rear(2), tag, lbecomm, status_front, ierr)

  tag = 118
  call mpi_recv(bnd(1,1,my+2,0), 1, xzplane_db5, front(2), tag, lbecomm, status_rear, ierr)

  call mpi_wait(reqs_front(1), status_front, ierr) 
  call mpi_wait(reqs_rear(1), status_rear, ierr)  
 

  end subroutine mpi_bnd

  

  subroutine mpi_dpst()
  use comvar_module
  use mpi

  integer      :: status_down(MPI_STATUS_SIZE)
  integer      :: reqs_down(1)
  integer      :: status_up(MPI_STATUS_SIZE)
  integer      :: reqs_up(1)
  integer      :: status_front(MPI_STATUS_SIZE)
  integer      :: reqs_front(1)
  integer      :: status_rear(MPI_STATUS_SIZE)
  integer      :: reqs_rear(1)

  integer      :: status_rear_down(MPI_STATUS_SIZE)
  integer      :: reqs_rear_down(1)
  integer      :: status_rear_up(MPI_STATUS_SIZE)
  integer      :: reqs_rear_up(1)
  integer      :: status_front_down(MPI_STATUS_SIZE)
  integer      :: reqs_front_down(1)
  integer      :: status_front_up(MPI_STATUS_SIZE)
  integer      :: reqs_front_up(1)
  integer      :: tag

  tag = 302
  call mpi_isend(dpst, 1, MPI_LOGICAL, up(2), tag, lbecomm, reqs_up(1), ierr)

  tag = 304
  call mpi_isend(dpst, 1, MPI_LOGICAL, down(2), tag, lbecomm, reqs_down(1), ierr)

  tag = 306
  call mpi_isend(dpst, 1, MPI_LOGICAL, front(2), tag, lbecomm, reqs_front(1), ierr)

  tag = 308
  call mpi_isend(dpst, 1, MPI_LOGICAL, rear(2), tag, lbecomm, reqs_rear(1), ierr)


  tag = 312
  call mpi_isend(dpst, 1, MPI_LOGICAL, front_up, tag, lbecomm, reqs_front_up(1), ierr)

  tag = 314
  call mpi_isend(dpst, 1, MPI_LOGICAL, front_down, tag, lbecomm, reqs_front_down(1), ierr)

  tag = 316
  call mpi_isend(dpst, 1, MPI_LOGICAL, rear_up, tag, lbecomm, reqs_rear_up(1), ierr)

  tag = 318
  call mpi_isend(dpst, 1, MPI_LOGICAL, rear_down, tag, lbecomm, reqs_rear_down(1), ierr)


  tag = 302
  call mpi_recv(dpst_up, 1, MPI_LOGICAL, down(2), tag, lbecomm, status_up, ierr)

  tag = 304
  call mpi_recv(dpst_down, 1, MPI_LOGICAL, up(2), tag, lbecomm, status_down, ierr)

  tag = 306
  call mpi_recv(dpst_front, 1, MPI_LOGICAL, rear(2), tag, lbecomm, status_front, ierr)

  tag = 308
  call mpi_recv(dpst_rear, 1, MPI_LOGICAL, front(2), tag, lbecomm, status_rear, ierr)


  tag = 312
  call mpi_recv(dpst_front_up, 1, MPI_LOGICAL, rear_down, tag, lbecomm, status_front_up, ierr)

  tag = 314
  call mpi_recv(dpst_front_down, 1, MPI_LOGICAL, rear_up, tag, lbecomm, status_front_down, ierr)

  tag = 316
  call mpi_recv(dpst_rear_up, 1, MPI_LOGICAL, front_down, tag, lbecomm, status_rear_up, ierr)

  tag = 318
  call mpi_recv(dpst_rear_down, 1, MPI_LOGICAL, front_up, tag, lbecomm, status_rear_down, ierr)


  call mpi_wait(reqs_up(1), status_up, ierr)
  call mpi_wait(reqs_down(1), status_down, ierr)
  call mpi_wait(reqs_front(1), status_front, ierr) 
  call mpi_wait(reqs_rear(1), status_rear, ierr)  

  call mpi_wait(reqs_front_up(1), status_front_up, ierr)
  call mpi_wait(reqs_front_down(1), status_front_down, ierr)
  call mpi_wait(reqs_rear_up(1), status_rear_up, ierr) 
  call mpi_wait(reqs_rear_down(1), status_rear_down, ierr) 

  end subroutine mpi_dpst 

end module mpi_comm_module
