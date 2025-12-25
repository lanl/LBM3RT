module output_module

contains


!******************************************************************
!****************  parallel_output  *****************************
!*************** Author: Ruoyu Li @ LANL 09/05/2024 ***************
!******************************************************************

  subroutine oflow_3d_d(io_x,io_y,io_z,istep,mat_3d,baseFilename,commsize,mpy,mpz) !output a 3d array with double precision
  
    implicit none

    integer :: i,j,k
    integer :: io_x,io_y,io_z
    integer :: istep,commsize,mpy,mpz
    integer :: unit_number,rec_length
    real*8, dimension(io_x,io_y,io_z) :: mat_3d
    real*8, dimension(io_z,io_y,io_x) :: io_reorder
    character(len=:), allocatable :: baseFilename, numStr, finalFilename
    integer, parameter :: MaxIntLength = 20

    character(len=5) :: mpyStr,mpzStr,commstr
    character(len=20) :: filename
    character(len=80) :: header 

    allocate(character(len=MaxIntLength) :: numStr)



    do i=1,io_x
      do j=1,io_y
        do k=1,io_z
          io_reorder(k,j,i) = mat_3d(i,j,k)
        enddo
      enddo
    enddo

    unit_number = 10  ! a free unit number for the file
    rec_length = io_x * io_y * io_z * 2

    WRITE(numStr, '(I9.9)') istep
    WRITE(mpyStr, '(I4.4)') mpy
    WRITE(mpzStr, '(I4.4)') mpz
    WRITE(commstr, '(I4.4)') commsize

    finalFilename = './Results/' // TRIM(baseFilename) // trim(numStr) //'_'// trim(commStr) //'_'// trim(mpyStr) //'_'// trim(mpzStr) // '.raw'

    print *,finalFilename
    
    
    open(unit=unit_number, file=finalFilename, form='unformatted'&
        , access='direct', recl=rec_length, action='write', status='replace')

    write(unit_number,rec=1) io_reorder
    close(unit_number)
  
  end subroutine oflow_3d_d

  subroutine oflow_3d_i(io_x,io_y,io_z,istep,mat_3d,baseFilename,commsize,mpy,mpz) !output a 3d array, integer
  
    implicit none

    integer :: i,j,k
    integer :: io_x,io_y,io_z
    integer :: istep,commsize,mpy,mpz
    integer :: unit_number,rec_length
    integer, dimension(io_x,io_y,io_z) :: mat_3d
    integer, dimension(io_z,io_y,io_x) :: io_reorder
    character(len=:), allocatable :: baseFilename, numStr, finalFilename
    integer, parameter :: MaxIntLength = 20

    character(len=5) :: mpyStr,mpzStr,commstr
    character(len=20) :: filename
    character(len=80) :: header 

    allocate(character(len=MaxIntLength) :: numStr)



    do i=1,io_x
      do j=1,io_y
        do k=1,io_z
          io_reorder(k,j,i) = mat_3d(i,j,k)
        enddo
      enddo
    enddo

    unit_number = 5  ! a free unit number for the file
    rec_length = io_x * io_y * io_z

    WRITE(numStr, '(I9.9)') istep
    WRITE(mpyStr, '(I4.4)') mpy
    WRITE(mpzStr, '(I4.4)') mpz
    WRITE(commstr, '(I4.4)') commsize

    finalFilename = './Results/' // TRIM(baseFilename) // trim(numStr) //'_'// trim(commStr) //'_'// trim(mpyStr) //'_'// trim(mpzStr) // '.raw'
    print *,finalFilename
    
    
    open(unit=unit_number, file=finalFilename, form='unformatted'&
        , access='direct', recl=rec_length, action='write', status='replace')

    write(unit_number,rec=1) io_reorder
    close(unit_number)
  
  end subroutine oflow_3d_i

  subroutine oflow_3d_bo_2(io_x,io_y,io_z,istep,mat_3d,baseFilename,commsize,mpy,mpz) !output a 3d array, logical
  
    implicit none

    integer :: i,j,k
    integer :: io_x,io_y,io_z
    integer :: istep,commsize,mpy,mpz
    integer :: unit_number,rec_length
    logical, dimension(io_x,io_y,0:io_z-1) :: mat_3d
    logical, dimension(io_z,io_y,io_x) :: io_reorder
    character(len=:), allocatable :: baseFilename, numStr, finalFilename
    integer, parameter :: MaxIntLength = 20

    character(len=5) :: mpyStr,mpzStr,commstr
    character(len=20) :: filename
    character(len=80) :: header 

    allocate(character(len=MaxIntLength) :: numStr)



    do i=1,io_x
      do j=1,io_y
        do k=1,io_z
          io_reorder(k,j,i) = mat_3d(i,j,k-1)
        enddo
      enddo
    enddo

    unit_number = 5  ! a free unit number for the file
    rec_length = io_x * io_y * io_z

    WRITE(numStr, '(I9.9)') istep
    WRITE(mpyStr, '(I4.4)') mpy
    WRITE(mpzStr, '(I4.4)') mpz
    WRITE(commstr, '(I4.4)') commsize

    finalFilename = './Results/' // TRIM(baseFilename) // trim(numStr) //'_'// trim(commStr) //'_'// trim(mpyStr) //'_'// trim(mpzStr) // '.raw'
    print *,finalFilename
    
    
    open(unit=unit_number, file=finalFilename, form='unformatted'&
        , access='direct', recl=rec_length, action='write', status='replace')

    write(unit_number,rec=1) io_reorder
    close(unit_number)
  
  end subroutine oflow_3d_bo_2

  subroutine oflow_3d_bo(io_x,io_y,io_z,istep,mat_3d,baseFilename,commsize,mpy,mpz) !output a 3d array, logical
  
    implicit none

    integer :: i,j,k
    integer :: io_x,io_y,io_z
    integer :: istep,commsize,mpy,mpz
    integer :: unit_number,rec_length
    logical, dimension(io_x,io_y,io_z) :: mat_3d
    logical, dimension(io_z,io_y,io_x) :: io_reorder
    character(len=:), allocatable :: baseFilename, numStr, finalFilename
    integer, parameter :: MaxIntLength = 20

    character(len=5) :: mpyStr,mpzStr,commstr
    character(len=20) :: filename
    character(len=80) :: header 

    allocate(character(len=MaxIntLength) :: numStr)



    do i=1,io_x
      do j=1,io_y
        do k=1,io_z
          io_reorder(k,j,i) = mat_3d(i,j,k)
        enddo
      enddo
    enddo

    unit_number = 5  ! a free unit number for the file
    rec_length = io_x * io_y * io_z

    WRITE(numStr, '(I9.9)') istep
    WRITE(mpyStr, '(I4.4)') mpy
    WRITE(mpzStr, '(I4.4)') mpz
    WRITE(commstr, '(I4.4)') commsize

    finalFilename = './Results/' // TRIM(baseFilename) // trim(numStr) //'_'// trim(commStr) //'_'// trim(mpyStr) //'_'// trim(mpzStr) // '.raw'
    print *,finalFilename
    
    
    open(unit=unit_number, file=finalFilename, form='unformatted'&
        , access='direct', recl=rec_length, action='write', status='replace')

    write(unit_number,rec=1) io_reorder
    close(unit_number)
  
  end subroutine oflow_3d_bo

  subroutine oflow_4d_bo(io_x,io_y,io_z,istep,mat_4d,d4,baseFilename,commsize,mpy,mpz) !output a 4d array, integer
  
    implicit none

    integer :: i,j,k,l
    integer :: io_x,io_y,io_z,d4
    integer :: istep,commsize,mpy,mpz
    integer :: unit_number,rec_length
    logical, dimension(d4,io_x,io_y,io_z) :: mat_4d
    logical, dimension(io_z,io_y,io_x,d4) :: io_reorder
    character(len=:), allocatable :: baseFilename, numStr, finalFilename
    integer, parameter :: MaxIntLength = 20

    character(len=5) :: mpyStr,mpzStr,commstr
    character(len=20) :: filename
    character(len=80) :: header 

    allocate(character(len=MaxIntLength) :: numStr)



    do i=1,io_x
      do j=1,io_y
        do k=1,io_z
          do l=1,d4
            io_reorder(k,j,i,l) = mat_4d(l,i,j,k)
          enddo
        enddo
      enddo
    enddo

    unit_number = 5  ! a free unit number for the file
    rec_length = io_x * io_y * io_z * d4

    WRITE(numStr, '(I9.9)') istep
    WRITE(mpyStr, '(I4.4)') mpy
    WRITE(mpzStr, '(I4.4)') mpz
    WRITE(commstr, '(I4.4)') commsize

    finalFilename = './Results/' // TRIM(baseFilename) // trim(numStr) //'_'// trim(commStr) //'_'// trim(mpyStr) //'_'// trim(mpzStr) // '.raw'
    print *,finalFilename
    
    
    open(unit=unit_number, file=finalFilename, form='unformatted'&
        , access='direct', recl=rec_length, action='write', status='replace')

    write(unit_number,rec=1) io_reorder
    close(unit_number)
  
  end subroutine oflow_4d_bo


  subroutine oflow_4d_bo_2(io_x,io_y,io_z,istep,mat_4d,d4,baseFilename,commsize,mpy,mpz) !output a 4d array, integer
  
    implicit none

    integer :: i,j,k,l
    integer :: io_x,io_y,io_z,d4
    integer :: istep,commsize,mpy,mpz
    integer :: unit_number,rec_length
    logical, dimension(d4,io_x,0:io_y-1,0:io_z-1) :: mat_4d
    logical, dimension(io_z,io_y,io_x,d4) :: io_reorder
    character(len=:), allocatable :: baseFilename, numStr, finalFilename
    integer, parameter :: MaxIntLength = 20

    character(len=5) :: mpyStr,mpzStr,commstr
    character(len=20) :: filename
    character(len=80) :: header 

    allocate(character(len=MaxIntLength) :: numStr)



    do i=1,io_x
      do j=1,io_y
        do k=1,io_z
          do l=1,d4
            io_reorder(k,j,i,l) = mat_4d(l,i,j-1,k-1)
          enddo
        enddo
      enddo
    enddo

    unit_number = 5  ! a free unit number for the file
    rec_length = io_x * io_y * io_z * d4

    WRITE(numStr, '(I9.9)') istep
    WRITE(mpyStr, '(I4.4)') mpy
    WRITE(mpzStr, '(I4.4)') mpz
    WRITE(commstr, '(I4.4)') commsize

    finalFilename = './Results/' // TRIM(baseFilename) // trim(numStr) //'_'// trim(commStr) //'_'// trim(mpyStr) //'_'// trim(mpzStr) // '.raw'
    print *,finalFilename
    
    
    open(unit=unit_number, file=finalFilename, form='unformatted'&
        , access='direct', recl=rec_length, action='write', status='replace')

    write(unit_number,rec=1) io_reorder
    close(unit_number)
  
  end subroutine oflow_4d_bo_2



  subroutine oflow_4d_i(io_x,io_y,io_z,istep,mat_4d,d4_start,d4_ed,baseFilename,commsize,mpy,mpz) !output a 4d array, integer
  
    implicit none

    integer :: i,j,k,s
    integer :: io_x,io_y,io_z,d4_start,d4_ed
    integer :: istep,commsize,mpy,mpz
    integer :: unit_number,rec_length
    integer, dimension(d4_start:d4_ed,io_x,io_y,io_z) :: mat_4d
    integer, dimension(io_z,io_y,io_x,d4_start:d4_ed) :: io_reorder
    character(len=:), allocatable :: baseFilename, numStr, finalFilename
    integer, parameter :: MaxIntLength = 20

    character(len=5) :: mpyStr,mpzStr,commstr
    character(len=20) :: filename
    character(len=80) :: header 

    allocate(character(len=MaxIntLength) :: numStr)



    do i=1,io_x
      do j=1,io_y
        do k=1,io_z
          do s=d4_start,d4_ed
            io_reorder(k,j,i,s) = mat_4d(s,i,j,k)
          enddo
        enddo
      enddo
    enddo

    unit_number = 5  ! a free unit number for the file
    rec_length = io_x * io_y * io_z * (d4_ed-d4_start+1)

    WRITE(numStr, '(I9.9)') istep
    WRITE(mpyStr, '(I4.4)') mpy
    WRITE(mpzStr, '(I4.4)') mpz
    WRITE(commstr, '(I4.4)') commsize

    finalFilename = './Results/' // TRIM(baseFilename) // trim(numStr) //'_'// trim(commStr) //'_'// trim(mpyStr) //'_'// trim(mpzStr) // '.raw'
    print *,finalFilename
    
    
    open(unit=unit_number, file=finalFilename, form='unformatted'&
        , access='direct', recl=rec_length, action='write', status='replace')

    write(unit_number,rec=1) io_reorder
    close(unit_number)
  
  end subroutine oflow_4d_i

  subroutine oflow_4d_d(io_x,io_y,io_z,istep,mat_4d,d4,baseFilename,commsize,mpy,mpz) !output a 4d array, integer
  
    implicit none

    integer :: i,j,k,s
    integer :: io_x,io_y,io_z,d4
    integer :: istep,commsize,mpy,mpz
    integer :: unit_number,rec_length
    real*8, dimension(d4,io_x,io_y,io_z) :: mat_4d
    real*8, dimension(io_z,io_y,io_x,d4) :: io_reorder
    character(len=:), allocatable :: baseFilename, numStr, finalFilename
    integer, parameter :: MaxIntLength = 20

    character(len=5) :: mpyStr,mpzStr,commstr
    character(len=20) :: filename
    character(len=80) :: header 

    allocate(character(len=MaxIntLength) :: numStr)



    do i=1,io_x
      do j=1,io_y
        do k=1,io_z
          do s=1,d4
            io_reorder(k,j,i,s) = mat_4d(s,i,j,k)
          enddo
        enddo
      enddo
    enddo

    unit_number = 5  ! a free unit number for the file
    rec_length = io_x * io_y * io_z * d4 * 2

    WRITE(numStr, '(I9.9)') istep
    WRITE(mpyStr, '(I4.4)') mpy
    WRITE(mpzStr, '(I4.4)') mpz
    WRITE(commstr, '(I4.4)') commsize

    finalFilename = './Results/' // TRIM(baseFilename) // trim(numStr) //'_'// trim(commStr) //'_'// trim(mpyStr) //'_'// trim(mpzStr) // '.raw'
    print *,finalFilename
    
    
    open(unit=unit_number, file=finalFilename, form='unformatted'&
        , access='direct', recl=rec_length, action='write', status='replace')

    write(unit_number,rec=1) io_reorder
    close(unit_number)
  
  end subroutine oflow_4d_d


  subroutine oflow_4d_d_2(io_x,io_y,io_z,istep,mat_4d,d4,baseFilename,commsize,mpy,mpz) !output a 4d array, integer
  
    implicit none

    integer :: i,j,k,s
    integer :: io_x,io_y,io_z,d4
    integer :: istep,commsize,mpy,mpz
    integer :: unit_number,rec_length
    real*8, dimension(d4,io_x,0:io_y-1,0:io_z-1) :: mat_4d
    real*8, dimension(io_z,io_y,io_x,d4) :: io_reorder
    character(len=:), allocatable :: baseFilename, numStr, finalFilename
    integer, parameter :: MaxIntLength = 20

    character(len=5) :: mpyStr,mpzStr,commstr
    character(len=20) :: filename
    character(len=80) :: header 

    allocate(character(len=MaxIntLength) :: numStr)



    do i=1,io_x
      do j=1,io_y
        do k=1,io_z
          do s=1,d4
            io_reorder(k,j,i,s) = mat_4d(s,i,j-1,k-1)
          enddo
        enddo
      enddo
    enddo

    unit_number = 5  ! a free unit number for the file
    rec_length = io_x * io_y * io_z * d4 * 2

    WRITE(numStr, '(I9.9)') istep
    WRITE(mpyStr, '(I4.4)') mpy
    WRITE(mpzStr, '(I4.4)') mpz
    WRITE(commstr, '(I4.4)') commsize

    finalFilename = './Results/' // TRIM(baseFilename) // trim(numStr) //'_'// trim(commStr) //'_'// trim(mpyStr) //'_'// trim(mpzStr) // '.raw'
    print *,finalFilename
    
    
    open(unit=unit_number, file=finalFilename, form='unformatted'&
        , access='direct', recl=rec_length, action='write', status='replace')

    write(unit_number,rec=1) io_reorder
    close(unit_number)
  
  end subroutine oflow_4d_d_2


  subroutine oflow_4d_f(io_x,io_y,io_z,istep,mat_4d,d4,baseFilename,commsize,mpy,mpz) !output a 4d array, integer
  
    implicit none

    integer :: i,j,k,s
    integer :: io_x,io_y,io_z,d4
    integer :: istep,commsize,mpy,mpz
    integer :: unit_number,rec_length
    real*8, dimension(io_x,io_y,io_z,0:d4-1) :: mat_4d
    real*8, dimension(d4,io_z,io_y,io_x) :: io_reorder
    character(len=:), allocatable :: baseFilename, numStr, finalFilename
    integer, parameter :: MaxIntLength = 20

    character(len=5) :: mpyStr,mpzStr,commstr
    character(len=20) :: filename
    character(len=80) :: header 

    allocate(character(len=MaxIntLength) :: numStr)



    do i=1,io_x
      do j=1,io_y
        do k=1,io_z
          do s=1,d4
            io_reorder(s,k,j,i) = mat_4d(i,j,k,s-1)
          enddo
        enddo
      enddo
    enddo

    unit_number = 5  ! a free unit number for the file
    rec_length = io_x * io_y * io_z * d4 * 2

    WRITE(numStr, '(I9.9)') istep
    WRITE(mpyStr, '(I4.4)') mpy
    WRITE(mpzStr, '(I4.4)') mpz
    WRITE(commstr, '(I4.4)') commsize

    finalFilename = './Results/' // TRIM(baseFilename) // trim(numStr) //'_'// trim(commStr) //'_'// trim(mpyStr) //'_'// trim(mpzStr) // '.raw'
    print *,finalFilename
    
    
    open(unit=unit_number, file=finalFilename, form='unformatted'&
        , access='direct', recl=rec_length, action='write', status='replace')

    write(unit_number,rec=1) io_reorder
    close(unit_number)
  
  end subroutine oflow_4d_f

  subroutine oflow_5d_g(io_x,io_y,io_z,istep,mat_5d,m,d4,baseFilename,commsize,mpy,mpz) !output a 5d array, integer
  
    implicit none

    integer :: i,j,k,s,mm
    integer :: io_x,io_y,io_z,d4,m
    integer :: istep,commsize,mpy,mpz
    integer :: unit_number,rec_length
    real*8, dimension(m,io_x,io_y,io_z,0:d4-1) :: mat_5d
    real*8, dimension(d4,io_z,io_y,io_x,m) :: io_reorder
    character(len=:), allocatable :: baseFilename, numStr, finalFilename
    integer, parameter :: MaxIntLength = 20

    character(len=5) :: mpyStr,mpzStr,commstr
    character(len=20) :: filename
    character(len=80) :: header 

    allocate(character(len=MaxIntLength) :: numStr)



    do i=1,io_x
      do j=1,io_y
        do k=1,io_z
          do s=1,d4
            do mm=1,m
              io_reorder(s,k,j,i,mm) = mat_5d(mm,i,j,k,s-1)
            enddo
          enddo
        enddo
      enddo
    enddo

    unit_number = 5  ! a free unit number for the file
    rec_length = io_x * io_y * io_z * d4 * m * 2

    WRITE(numStr, '(I9.9)') istep
    WRITE(mpyStr, '(I4.4)') mpy
    WRITE(mpzStr, '(I4.4)') mpz
    WRITE(commstr, '(I4.4)') commsize

    finalFilename = './Results/' // TRIM(baseFilename) // trim(numStr) //'_'// trim(commStr) //'_'// trim(mpyStr) //'_'// trim(mpzStr) // '.raw'
    print *,finalFilename
    
    
    open(unit=unit_number, file=finalFilename, form='unformatted'&
        , access='direct', recl=rec_length, action='write', status='replace')

    write(unit_number,rec=1) io_reorder
    close(unit_number)
  
  end subroutine oflow_5d_g


!******************************************************************
!****************  mpi_gather_output  *****************************
!*************** Author: Ruoyu Li @ LANL 09/05/2024 ***************
!******************************************************************


  subroutine gather_out_3d_dble(mat_3d,save_dir,baseFilename,istep,ix0,ixn,iy0,iyn,iz0,izn,ybuff,zbuff,commsize,npy,npz,myrank)
  
  use comvar_module
  use mpi
  implicit none
  integer :: ix0,ixn,iy0,iyn,iz0,izn,ybuff,zbuff,commsize,myrank,npy,npz
  integer :: block_size, status(MPI_STATUS_SIZE),send_rank
  integer, ALLOCATABLE:: send_request(:),recv_request(:)
  integer :: i,j,k,p,jj,kk,istep,ipy,jpy
  real*8, ALLOCATABLE :: global_mat(:,:,:),io_reorder(:,:,:),recv_buf(:,:,:),send_buf(:,:,:)
  real*8 :: mat_3d(ix0:ixn,iy0:iyn,iz0:izn)

  character(len=:), allocatable :: baseFilename, numStr, finalFilename
  character(len=*), intent(in) :: save_dir
  integer, parameter :: MaxIntLength = 20
  integer :: unit_number,rec_length

  character(len=5) :: commstr
  character(len=20) :: filename
  character(len=80) :: header 
  allocate(character(len=MaxIntLength) :: numStr)
  ALLOCATE(send_request(1:commsize-1))
  ALLOCATE(recv_request(1:commsize-1))

      
  if (myrank == 0) then
     ALLOCATE(recv_buf(1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1))
     ALLOCATE(global_mat(1:lx,1:ly,1:lz))
     ALLOCATE(io_reorder(1:lz,1:ly,1:lx))
  else
     allocate(send_buf(1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1))
    do i=1,lx
      do j=iy0+ybuff,iyn-ybuff
        do k=iz0+zbuff,izn-zbuff
          jj = j-iy0-ybuff+1
          kk = k-iz0-zbuff+1
          send_buf(i,jj,kk) = mat_3d(i,j,k)
        enddo
      enddo
    enddo     
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (myrank == 0) then
    do i=1,mx
      do j=1,my
        do k=1,mz
          jj = iy0+ybuff+j-1
          kk = iz0+zbuff+k-1
          global_mat(i,j,k) = mat_3d(i,jj,kk)
        enddo
      enddo
    enddo
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
!isend&recv
    do ipy = 0, npy-1
        do jpy = 0, npz-1
          send_rank = jpy+ipy*npz         
          if(send_rank/=0) then
              block_size = my_local(ipy+1) * mz_local(jpy+1) * lx
              if (myrank == send_rank) then
                call MPI_Isend(send_buf, block_size, MPI_DOUBLE_PRECISION, 0, 1200+send_rank, lbecomm, send_request(send_rank),ierr)
                call MPI_Wait(send_request(send_rank), status, ierr)   
              endif              
          endif
        enddo
    enddo

    do ipy = 0, npy-1
        do jpy = 0, npz-1
          send_rank = jpy+ipy*npz         
          if(send_rank/=0) then
              block_size = my_local(ipy+1) * mz_local(jpy+1) * lx
              if (myrank == 0) then
                call MPI_Recv(recv_buf, block_size, MPI_DOUBLE_PRECISION, send_rank, 1200+send_rank, lbecomm, status, ierr)

                do i=1,lx
                  do j=1,iyn-iy0-2*ybuff+1
                    do k=1,izn-iz0-2*zbuff+1
                      jj = ny_local(ipy)+j
                      kk = nz_local(jpy)+k
                      global_mat(i,jj,kk) = recv_buf(i,j,k)
                    enddo
                  enddo
                enddo

              endif
          endif
        enddo
    enddo
    
   call MPI_Barrier(MPI_COMM_WORLD, ierr)

  !output file
  if (myrank == 0) then
    do i = 1,lx
      do j = 1,ly
        do k = 1,lz
          io_reorder(k,j,i) = global_mat(i,j,k)
        enddo
      enddo
    enddo


    unit_number = 5  ! a free unit number for the file
    rec_length = lx * ly * lz * 2

    WRITE(numStr, '(I9.9)') istep
    WRITE(commstr, '(I4.4)') commsize

    finalFilename = TRIM(save_dir) // TRIM(baseFilename) // trim(numStr) //'_'// trim(commStr) // '.raw'
    print *,finalFilename
    
    
    open(unit=unit_number, file=finalFilename, form='unformatted'&
        , access='direct', recl=rec_length, action='write', status='replace')

    write(unit_number,rec=1) io_reorder
    close(unit_number)

  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)


  if (myrank == 0) then
     deallocate(global_mat,stat = ierr)
     deallocate(io_reorder,stat = ierr)
     deallocate(recv_buf,stat = ierr)
  else
     deallocate(send_buf,stat = ierr)
  endif

  end subroutine gather_out_3d_dble


  subroutine gather_out_4d_dble(mat_4d,save_dir,baseFilename,istep,mm,ix0,ixn,iy0,iyn,iz0,izn,ybuff,zbuff,commsize,npy,npz,myrank)
  
  use comvar_module
  use mpi
  implicit none
  integer :: mm,ix0,ixn,iy0,iyn,iz0,izn,ybuff,zbuff,commsize,myrank,npy,npz
  integer :: block_size, status(MPI_STATUS_SIZE),send_rank
  integer, ALLOCATABLE:: send_request(:),recv_request(:)
  integer :: i,j,k,p,jj,kk,istep,ipy,jpy
  real*8, ALLOCATABLE :: global_mat(:,:,:,:),io_reorder(:,:,:,:),recv_buf(:,:,:,:),send_buf(:,:,:,:)
  real*8 :: mat_4d(1:mm,ix0:ixn,iy0:iyn,iz0:izn)

  character(len=:), allocatable :: baseFilename, numStr, finalFilename
  character(len=*), intent(in) :: save_dir
  integer, parameter :: MaxIntLength = 20
  integer :: unit_number,rec_length

  character(len=5) :: commstr
  character(len=20) :: filename
  character(len=80) :: header 
  allocate(character(len=MaxIntLength) :: numStr)
  ALLOCATE(send_request(1:commsize-1))
  ALLOCATE(recv_request(1:commsize-1))

      
  if (myrank == 0) then
     ALLOCATE(recv_buf(mm,1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1))
     ALLOCATE(global_mat(mm,1:lx,1:ly,1:lz))
     ALLOCATE(io_reorder(1:lz,1:ly,1:lx,1:mm))
  else
     allocate(send_buf(1:mm,1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1))
     do p=1,mm
      do i=1,lx
        do j=iy0+ybuff,iyn-ybuff
          do k=iz0+zbuff,izn-zbuff
            jj = j-iy0-ybuff+1
            kk = k-iz0-zbuff+1
            send_buf(p,i,jj,kk) = mat_4d(p,i,j,k)
          enddo
        enddo
      enddo
    enddo
     
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (myrank == 0) then
    do p=1,mm
      do i=1,mx
        do j=1,my
          do k=1,mz
            jj = iy0+ybuff+j-1
            kk = iz0+zbuff+k-1
            global_mat(p,i,j,k) = mat_4d(p,i,jj,kk)
          enddo
        enddo
      enddo
    enddo
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
!isend&recv
    do ipy = 0, npy-1
        do jpy = 0, npz-1
          send_rank = jpy+ipy*npz         
          if(send_rank/=0) then
              block_size = my_local(ipy+1) * mz_local(jpy+1) * mm * lx
              if (myrank == send_rank) then
                call MPI_Isend(send_buf, block_size, MPI_DOUBLE_PRECISION, 0, 800+send_rank, lbecomm, send_request(send_rank),ierr)
                call MPI_Wait(send_request(send_rank), status, ierr)   
              endif              
          endif
        enddo
    enddo

    do ipy = 0, npy-1
        do jpy = 0, npz-1
          send_rank = jpy+ipy*npz         
          if(send_rank/=0) then
              block_size = my_local(ipy+1) * mz_local(jpy+1) * mm * lx
              if (myrank == 0) then
                call MPI_Recv(recv_buf, block_size, MPI_DOUBLE_PRECISION, send_rank, 800+send_rank, lbecomm, status, ierr)

                do p=1,mm
                  do i=1,lx
                    do j=1,iyn-iy0-2*ybuff+1
                      do k=1,izn-iz0-2*zbuff+1
                        jj = ny_local(ipy)+j
                        kk = nz_local(jpy)+k
                        global_mat(p,i,jj,kk) = recv_buf(p,i,j,k)
                      enddo
                    enddo
                  enddo
                enddo

              endif
          endif
        enddo
    enddo
    
   call MPI_Barrier(MPI_COMM_WORLD, ierr)

  !output file
  if (myrank == 0) then
    do i = 1,lx
      do j = 1,ly
        do k = 1,lz
          do p = 1,mm
            io_reorder(k,j,i,p) = global_mat(p,i,j,k)
          enddo
        enddo
      enddo
    enddo


    unit_number = 5  ! a free unit number for the file
    rec_length = lx * ly * lz * mm * 2

    WRITE(numStr, '(I9.9)') istep
    WRITE(commstr, '(I4.4)') commsize

    finalFilename = TRIM(save_dir) // TRIM(baseFilename) // trim(numStr) //'_'// trim(commStr) // '.raw'
    print *,finalFilename
    
    
    open(unit=unit_number, file=finalFilename, form='unformatted'&
        , access='direct', recl=rec_length, action='write', status='replace')

    write(unit_number,rec=1) io_reorder
    close(unit_number)

  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)


  if (myrank == 0) then
     deallocate(global_mat,stat = ierr)
     deallocate(io_reorder,stat = ierr)
     deallocate(recv_buf,stat = ierr)
  else
     deallocate(send_buf,stat = ierr)
  endif

  end subroutine gather_out_4d_dble


  subroutine gather_out_fi(mat_4d,save_dir,baseFilename,istep,ix0,ixn,iy0,iyn,iz0,izn,d4i,d4e,ybuff,zbuff,commsize,npy,npz,myrank)
  
  use comvar_module
  use mpi
  implicit none
  integer :: ix0,ixn,iy0,iyn,iz0,izn,d4i,d4e,ybuff,zbuff,commsize,myrank,npy,npz
  integer :: block_size, status(MPI_STATUS_SIZE),send_rank
  integer, ALLOCATABLE:: send_request(:),recv_request(:)
  integer :: i,j,k,p,jj,kk,istep,ipy,jpy,pp
  real*8, ALLOCATABLE :: global_mat(:,:,:,:),io_reorder(:,:,:,:),recv_buf(:,:,:,:),send_buf(:,:,:,:)
  real*8 :: mat_4d(ix0:ixn,iy0:iyn,iz0:izn,d4i:d4e)

  character(len=:), allocatable :: baseFilename, numStr, finalFilename
  character(len=*), intent(in) :: save_dir
  integer, parameter :: MaxIntLength = 20
  integer :: unit_number,rec_length

  character(len=5) :: commstr
  character(len=20) :: filename
  character(len=80) :: header 
  allocate(character(len=MaxIntLength) :: numStr)
  ALLOCATE(send_request(1:commsize-1))
  ALLOCATE(recv_request(1:commsize-1))

      
  if (myrank == 0) then
     ALLOCATE(recv_buf(1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1,1:d4e-d4i+1))
     ALLOCATE(global_mat(1:lx,1:ly,1:lz,1:d4e-d4i+1))
     ALLOCATE(io_reorder(1:d4e-d4i+1,1:lz,1:ly,1:lx))
  else
     allocate(send_buf(1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1,1:d4e-d4i+1))
     do p=d4i,d4e
      do i=1,lx
        do j=iy0+ybuff,iyn-ybuff
          do k=iz0+zbuff,izn-zbuff
            jj = j-iy0-ybuff+1
            kk = k-iz0-zbuff+1
            pp = p-d4i+1
            send_buf(i,jj,kk,pp) = mat_4d(i,j,k,p)
          enddo
        enddo
      enddo
    enddo
     
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (myrank == 0) then
    do p=1,d4e-d4i+1
      do i=1,mx
        do j=1,my
          do k=1,mz
            jj = iy0+ybuff+j-1
            kk = iz0+zbuff+k-1
            pp = p+d4i-1
            global_mat(i,j,k,p) = mat_4d(i,jj,kk,pp)
          enddo
        enddo
      enddo
    enddo
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
!isend&recv
    do ipy = 0, npy-1
        do jpy = 0, npz-1
          send_rank = jpy+ipy*npz         
          if(send_rank/=0) then
              block_size = my_local(ipy+1) * mz_local(jpy+1) * (d4e-d4i+1)  * lx
              if (myrank == send_rank) then
                call MPI_Isend(send_buf, block_size, MPI_DOUBLE_PRECISION, 0, 1100+send_rank, lbecomm, send_request(send_rank),ierr)
                call MPI_Wait(send_request(send_rank), status, ierr)   
              endif              
          endif
        enddo
    enddo

    do ipy = 0, npy-1
        do jpy = 0, npz-1
          send_rank = jpy+ipy*npz         
          if(send_rank/=0) then
              block_size = my_local(ipy+1) * mz_local(jpy+1) * (d4e-d4i+1)  * lx
              if (myrank == 0) then
                call MPI_Recv(recv_buf, block_size, MPI_DOUBLE_PRECISION, send_rank, 1100+send_rank, lbecomm, status, ierr)

                do p=1,d4e-d4i+1
                  do i=1,lx
                    do j=1,iyn-iy0-2*ybuff+1
                      do k=1,izn-iz0-2*zbuff+1
                        jj = ny_local(ipy)+j
                        kk = nz_local(jpy)+k
                        global_mat(i,jj,kk,p) = recv_buf(i,j,k,p)
                      enddo
                    enddo
                  enddo
                enddo

              endif
          endif
        enddo
    enddo
    
   call MPI_Barrier(MPI_COMM_WORLD, ierr)

  !output file
  if (myrank == 0) then
    do i = 1,lx
      do j = 1,ly
        do k = 1,lz
          do p = 1,d4e-d4i+1
            io_reorder(p,k,j,i) = global_mat(i,j,k,p)
          enddo
        enddo
      enddo
    enddo


    unit_number = 5  ! a free unit number for the file
    rec_length = lx * ly * lz * (d4e-d4i+1)  * 2

    WRITE(numStr, '(I9.9)') istep
    WRITE(commstr, '(I4.4)') commsize

    finalFilename = TRIM(save_dir) // TRIM(baseFilename) // trim(numStr) //'_'// trim(commStr) // '.raw'
    print *,finalFilename
    
    
    open(unit=unit_number, file=finalFilename, form='unformatted'&
        , access='direct', recl=rec_length, action='write', status='replace')

    write(unit_number,rec=1) io_reorder
    close(unit_number)

  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)


  if (myrank == 0) then
     deallocate(global_mat,stat = ierr)
     deallocate(io_reorder,stat = ierr)
     deallocate(recv_buf,stat = ierr)
  else
     deallocate(send_buf,stat = ierr)
  endif

  end subroutine gather_out_fi


  subroutine gather_out_gi(mat_5d,save_dir,baseFilename,istep,mm,ix0,ixn,iy0,iyn,iz0,izn,d5i,d5e,ybuff,zbuff,commsize,npy,npz,myrank)
  
  use comvar_module
  use mpi
  implicit none
  integer :: mm,ix0,ixn,iy0,iyn,iz0,izn,d5i,d5e,ybuff,zbuff,commsize,myrank,npy,npz
  integer :: block_size, status(MPI_STATUS_SIZE),send_rank
  integer, ALLOCATABLE:: send_request(:),recv_request(:)
  integer :: i,j,k,p,jj,kk,istep,ipy,jpy,pp,s
  real*8, ALLOCATABLE :: global_mat(:,:,:,:,:),io_reorder(:,:,:,:,:),recv_buf(:,:,:,:,:),send_buf(:,:,:,:,:)
  real*8 :: mat_5d(1:mm,ix0:ixn,iy0:iyn,iz0:izn,d5i:d5e)

  character(len=:), allocatable :: baseFilename, numStr, finalFilename
  character(len=*), intent(in) :: save_dir
  integer, parameter :: MaxIntLength = 20
  integer :: unit_number,rec_length

  character(len=5) :: commstr
  character(len=20) :: filename
  character(len=80) :: header 
  allocate(character(len=MaxIntLength) :: numStr)
  ALLOCATE(send_request(1:commsize-1))
  ALLOCATE(recv_request(1:commsize-1))

      
  if (myrank == 0) then
     ALLOCATE(recv_buf(1:mm,1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1,1:d5e-d5i+1))
     ALLOCATE(global_mat(1:mm,1:lx,1:ly,1:lz,1:d5e-d5i+1))
     ALLOCATE(io_reorder(1:d5e-d5i+1,1:lz,1:ly,1:lx,1:mm))
  else
     allocate(send_buf(1:mm,1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1,1:d5e-d5i+1))
     do s = 1,mm
      do p=d5i,d5e
        do i=1,lx
          do j=iy0+ybuff,iyn-ybuff
            do k=iz0+zbuff,izn-zbuff
              jj = j-iy0-ybuff+1
              kk = k-iz0-zbuff+1
              pp = p-d5i+1
              send_buf(s,i,jj,kk,pp) = mat_5d(s,i,j,k,p)
            enddo
          enddo
        enddo
      enddo
     enddo
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (myrank == 0) then
    do s = 1,mm
      do p=1,d5e-d5i+1
        do i=1,mx
          do j=1,my
            do k=1,mz
              jj = iy0+ybuff+j-1
              kk = iz0+zbuff+k-1
              pp = p+d5i-1
              global_mat(s,i,j,k,p) = mat_5d(s,i,jj,kk,pp)
            enddo
          enddo
        enddo
      enddo
    enddo
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
!isend&recv
    do ipy = 0, npy-1
        do jpy = 0, npz-1
          send_rank = jpy+ipy*npz         
          if(send_rank/=0) then
              block_size = my_local(ipy+1) * mz_local(jpy+1) * (d5e-d5i+1)  * mm * lx
              if (myrank == send_rank) then
                call MPI_Isend(send_buf, block_size, MPI_DOUBLE_PRECISION, 0, 1400+send_rank, lbecomm, send_request(send_rank),ierr)
                call MPI_Wait(send_request(send_rank), status, ierr)   
              endif              
          endif
        enddo
    enddo

    do ipy = 0, npy-1
        do jpy = 0, npz-1
          send_rank = jpy+ipy*npz         
          if(send_rank/=0) then
              block_size = my_local(ipy+1) * mz_local(jpy+1) * (d5e-d5i+1)  * mm * lx
              if (myrank == 0) then
                call MPI_Recv(recv_buf, block_size, MPI_DOUBLE_PRECISION, send_rank, 1400+send_rank, lbecomm, status, ierr)

                do s = 1,mm
                  do p=1,d5e-d5i+1
                    do i=1,lx
                      do j=1,iyn-iy0-2*ybuff+1
                        do k=1,izn-iz0-2*zbuff+1
                          jj = ny_local(ipy)+j
                          kk = nz_local(jpy)+k
                          global_mat(s,i,jj,kk,p) = recv_buf(s,i,j,k,p)
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
              endif
          endif
        enddo
    enddo
    
   call MPI_Barrier(MPI_COMM_WORLD, ierr)

  !output file
  if (myrank == 0) then
    do s = 1,mm
      do i = 1,lx
        do j = 1,ly
          do k = 1,lz
            do p = 1,d5e-d5i+1
              io_reorder(p,k,j,i,s) = global_mat(s,i,j,k,p)
            enddo
          enddo
        enddo
      enddo
    enddo

    unit_number = 5  ! a free unit number for the file
    rec_length = lx * ly * lz * (d5e-d5i+1)  * mm * 2

    WRITE(numStr, '(I9.9)') istep
    WRITE(commstr, '(I4.4)') commsize

    finalFilename = TRIM(save_dir) // TRIM(baseFilename) // trim(numStr) //'_'// trim(commStr) // '.raw'
    print *,finalFilename
    
    
    open(unit=unit_number, file=finalFilename, form='unformatted'&
        , access='direct', recl=rec_length, action='write', status='replace')

    write(unit_number,rec=1) io_reorder
    close(unit_number)

  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)


  if (myrank == 0) then
     deallocate(global_mat,stat = ierr)
     deallocate(io_reorder,stat = ierr)
     deallocate(recv_buf,stat = ierr)
  else
     deallocate(send_buf,stat = ierr)
  endif

  end subroutine gather_out_gi



  subroutine gather_out_4d_int(mat_4d,save_dir,baseFilename,istep,mm,ix0,ixn,iy0,iyn,iz0,izn,ybuff,zbuff,commsize,npy,npz,myrank)
  
  use comvar_module
  use mpi
  implicit none
  integer :: mm,ix0,ixn,iy0,iyn,iz0,izn,ybuff,zbuff,commsize,myrank,npy,npz
  integer :: block_size, status(MPI_STATUS_SIZE),send_rank
  integer, ALLOCATABLE:: send_request(:),recv_request(:)
  integer :: i,j,k,p,jj,kk,istep,ipy,jpy
  integer, ALLOCATABLE :: global_mat(:,:,:,:),io_reorder(:,:,:,:),recv_buf(:,:,:,:),send_buf(:,:,:,:)
  integer :: mat_4d(1:mm,ix0:ixn,iy0:iyn,iz0:izn)

  character(len=:), allocatable :: baseFilename, numStr, finalFilename
  character(len=*), intent(in) :: save_dir
  integer, parameter :: MaxIntLength = 20
  integer :: unit_number,rec_length

  character(len=5) :: commstr
  character(len=20) :: filename
  character(len=80) :: header 
  allocate(character(len=MaxIntLength) :: numStr)
  ALLOCATE(send_request(1:commsize-1))
  ALLOCATE(recv_request(1:commsize-1))

      
  if (myrank == 0) then
     ALLOCATE(recv_buf(mm,1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1))
     ALLOCATE(global_mat(mm,1:lx,1:ly,1:lz))
     ALLOCATE(io_reorder(1:lz,1:ly,1:lx,1:mm))
  else
     allocate(send_buf(1:mm,1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1))
     do p=1,mm
      do i=1,lx
        do j=iy0+ybuff,iyn-ybuff
          do k=iz0+zbuff,izn-zbuff
            jj = j-iy0-ybuff+1
            kk = k-iz0-zbuff+1
            send_buf(p,i,jj,kk) = mat_4d(p,i,j,k)
          enddo
        enddo
      enddo
    enddo
     
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (myrank == 0) then
    do p=1,mm
      do i=1,mx
        do j=1,my
          do k=1,mz
            jj = iy0+ybuff+j-1
            kk = iz0+zbuff+k-1
            global_mat(p,i,j,k) = mat_4d(p,i,jj,kk)
          enddo
        enddo
      enddo
    enddo
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
!isend&recv
    do ipy = 0, npy-1
        do jpy = 0, npz-1
          send_rank = jpy+ipy*npz         
          if(send_rank/=0) then
              block_size = my_local(ipy+1) * mz_local(jpy+1) * mm * lx
              if (myrank == send_rank) then
                call MPI_Isend(send_buf, block_size, MPI_INTEGER, 0, 700+send_rank, lbecomm, send_request(send_rank),ierr)
                call MPI_Wait(send_request(send_rank), status, ierr)   
              endif              
          endif
        enddo
    enddo

    do ipy = 0, npy-1
        do jpy = 0, npz-1
          send_rank = jpy+ipy*npz         
          if(send_rank/=0) then
              block_size = my_local(ipy+1) * mz_local(jpy+1) * mm * lx
              if (myrank == 0) then
                call MPI_Recv(recv_buf, block_size, MPI_INTEGER, send_rank, 700+send_rank, lbecomm, status, ierr)

                do p=1,mm
                  do i=1,lx
                    do j=1,iyn-iy0-2*ybuff+1
                      do k=1,izn-iz0-2*zbuff+1
                        jj = ny_local(ipy)+j
                        kk = nz_local(jpy)+k
                        global_mat(p,i,jj,kk) = recv_buf(p,i,j,k)
                      enddo
                    enddo
                  enddo
                enddo

              endif
          endif
        enddo
    enddo
    
   call MPI_Barrier(MPI_COMM_WORLD, ierr)

  !output file
  if (myrank == 0) then
    do i = 1,lx
      do j = 1,ly
        do k = 1,lz
          do p = 1,mm
            io_reorder(k,j,i,p) = global_mat(p,i,j,k)
          enddo
        enddo
      enddo
    enddo


    unit_number = 5  ! a free unit number for the file
    rec_length = lx * ly * lz * mm

    WRITE(numStr, '(I9.9)') istep
    WRITE(commstr, '(I4.4)') commsize

    finalFilename = TRIM(save_dir) // TRIM(baseFilename) // trim(numStr) //'_'// trim(commStr) // '.raw'
    print *,finalFilename
    
    
    open(unit=unit_number, file=finalFilename, form='unformatted'&
        , access='direct', recl=rec_length, action='write', status='replace')

    write(unit_number,rec=1) io_reorder
    close(unit_number)

  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)


  if (myrank == 0) then
     deallocate(global_mat,stat = ierr)
     deallocate(io_reorder,stat = ierr)
     deallocate(recv_buf,stat = ierr)
  else
     deallocate(send_buf,stat = ierr)
  endif

  end subroutine gather_out_4d_int


  subroutine gather_out_4d_logic(mat_4d,save_dir,baseFilename,istep,mm,ix0,ixn,iy0,iyn,iz0,izn,ybuff,zbuff,commsize,npy,npz,myrank)
  
  use comvar_module
  use mpi
  implicit none
  integer :: mm,ix0,ixn,iy0,iyn,iz0,izn,ybuff,zbuff,commsize,myrank,npy,npz
  integer :: block_size, status(MPI_STATUS_SIZE),send_rank
  integer, ALLOCATABLE:: send_request(:),recv_request(:)
  integer :: i,j,k,p,jj,kk,istep,ipy,jpy
  logical, ALLOCATABLE :: global_mat(:,:,:,:),io_reorder(:,:,:,:),recv_buf(:,:,:,:),send_buf(:,:,:,:)
  logical :: mat_4d(1:mm,ix0:ixn,iy0:iyn,iz0:izn)

  character(len=:), allocatable :: baseFilename, numStr, finalFilename
  character(len=*), intent(in) :: save_dir
  integer, parameter :: MaxIntLength = 20
  integer :: unit_number,rec_length

  character(len=5) :: commstr
  character(len=20) :: filename
  character(len=80) :: header 
  allocate(character(len=MaxIntLength) :: numStr)
  ALLOCATE(send_request(1:commsize-1))
  ALLOCATE(recv_request(1:commsize-1))

      
  if (myrank == 0) then
     ALLOCATE(recv_buf(mm,1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1))
     ALLOCATE(global_mat(mm,1:lx,1:ly,1:lz))
     ALLOCATE(io_reorder(1:lz,1:ly,1:lx,1:mm))
  else
     allocate(send_buf(1:mm,1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1))
     do p=1,mm
      do i=1,lx
        do j=iy0+ybuff,iyn-ybuff
          do k=iz0+zbuff,izn-zbuff
            jj = j-iy0-ybuff+1
            kk = k-iz0-zbuff+1
            send_buf(p,i,jj,kk) = mat_4d(p,i,j,k)
          enddo
        enddo
      enddo
    enddo
     
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (myrank == 0) then
    do p=1,mm
      do i=1,mx
        do j=1,my
          do k=1,mz
            jj = iy0+ybuff+j-1
            kk = iz0+zbuff+k-1
            global_mat(p,i,j,k) = mat_4d(p,i,jj,kk)
          enddo
        enddo
      enddo
    enddo
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
!isend&recv
    do ipy = 0, npy-1
        do jpy = 0, npz-1
          send_rank = jpy+ipy*npz         
          if(send_rank/=0) then
              block_size = my_local(ipy+1) * mz_local(jpy+1) * mm * lx
              if (myrank == send_rank) then
                call MPI_Isend(send_buf, block_size, MPI_LOGICAL, 0, 600+send_rank, lbecomm, send_request(send_rank),ierr)
                call MPI_Wait(send_request(send_rank), status, ierr)   
              endif              
          endif
        enddo
    enddo

    do ipy = 0, npy-1
        do jpy = 0, npz-1
          send_rank = jpy+ipy*npz         
          if(send_rank/=0) then
              block_size = my_local(ipy+1) * mz_local(jpy+1) * mm * lx
              if (myrank == 0) then
                call MPI_Recv(recv_buf, block_size, MPI_LOGICAL, send_rank, 600+send_rank, lbecomm, status, ierr)

                do p=1,mm
                  do i=1,lx
                    do j=1,iyn-iy0-2*ybuff+1
                      do k=1,izn-iz0-2*zbuff+1
                        jj = ny_local(ipy)+j
                        kk = nz_local(jpy)+k
                        global_mat(p,i,jj,kk) = recv_buf(p,i,j,k)
                      enddo
                    enddo
                  enddo
                enddo

              endif
          endif
        enddo
    enddo

   call MPI_Barrier(MPI_COMM_WORLD, ierr)

  !output file
  if (myrank == 0) then
    do i = 1,lx
      do j = 1,ly
        do k = 1,lz
          do p = 1,mm
            io_reorder(k,j,i,p) = global_mat(p,i,j,k)
          enddo
        enddo
      enddo
    enddo


    unit_number = 5  ! a free unit number for the file
    rec_length = lx * ly * lz * mm

    WRITE(numStr, '(I9.9)') istep
    WRITE(commstr, '(I4.4)') commsize

    finalFilename = TRIM(save_dir) // TRIM(baseFilename) // trim(numStr) //'_'// trim(commStr) // '.raw'
    print *,finalFilename
    
    
    open(unit=unit_number, file=finalFilename, form='unformatted'&
        , access='direct', recl=rec_length, action='write', status='replace')

    write(unit_number,rec=1) io_reorder
    close(unit_number)

  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)


  if (myrank == 0) then
     deallocate(global_mat,stat = ierr)
     deallocate(io_reorder,stat = ierr)
     deallocate(recv_buf,stat = ierr)
  else
     deallocate(send_buf,stat = ierr)
  endif

  end subroutine gather_out_4d_logic


!******************************************************************
!****************  mpi_scatter_input  *****************************
!*************** Author: Ruoyu Li @ LANL 09/05/2024 ***************
!******************************************************************



  subroutine scatter_in_4d_dble(mat_4d,save_dir,baseFilename,istep,mm,ix0,ixn,iy0,iyn,iz0,izn,ybuff,zbuff,commsize,npy,npz,myrank)
  
  use comvar_module
  use mpi
  implicit none
  integer :: mm,ix0,ixn,iy0,iyn,iz0,izn,ybuff,zbuff,commsize,myrank,npy,npz
  integer :: block_size, status(MPI_STATUS_SIZE),recv_rank,iounit
  integer, ALLOCATABLE:: send_request(:),recv_request(:)
  integer :: i,j,k,p,jj,kk,istep,ipy,jpy
  real*8, ALLOCATABLE :: global_mat(:,:,:,:),io_reorder(:,:,:,:),recv_buf(:,:,:,:),send_buf(:,:,:,:)
  real*8 :: mat_4d(1:mm,ix0:ixn,iy0:iyn,iz0:izn)

  character(len=:), allocatable :: baseFilename, numStr, finalFilename
  character(len=*), intent(in) :: save_dir
  integer, parameter :: MaxIntLength = 20
  integer :: unit_number,rec_length

  character(len=5) :: commstr
  character(len=20) :: filename
  character(len=80) :: header 
  allocate(character(len=MaxIntLength) :: numStr)
  ALLOCATE(send_request(1:commsize-1))
  ALLOCATE(recv_request(1:commsize-1))

      
  if (myrank == 0) then
     ALLOCATE(send_buf(mm,1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1))
     ALLOCATE(global_mat(mm,1:lx,1:ly,1:lz))
     ALLOCATE(io_reorder(1:lz,1:ly,1:lx,1:mm))
  else
     allocate(recv_buf(1:mm,1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1)) 
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (myrank == 0) then
      rec_length = lx * ly * lz * mm * 2
      WRITE(numStr, '(I9.9)') istep
      WRITE(commstr, '(I4.4)') commsize
      filename = './Results/'//baseFilename//numStr//'.raw'
      finalFilename = TRIM(save_dir) // TRIM(baseFilename) // trim(numStr) //'_'// trim(commStr) // '.raw'
      open(unit=iounit, file=finalFilename, form='unformatted', access='direct', &
             recl=rec_length, status='old', iostat=ierr)
      if (ierr /= 0) then
        print*, 'Error opening file:', finalFilename
        stop
      end if
      read(iounit, rec=1, iostat=ierr) io_reorder
      if (ierr /= 0) then
        print*, 'Error reading file:', finalFilename
        stop
      end if
      close(iounit)

      do i = 1,lx
        do j = 1,ly
          do k = 1,lz
            do p = 1,mm
              global_mat(p,i,j,k) = io_reorder(k,j,i,p)
            enddo
          enddo
        enddo
      enddo
      print *,'Input file loaded: ',finalFilename

  endif


  if (myrank == 0) then
    do p=1,mm
      do i=1,mx
        do j=1,my
          do k=1,mz
            jj = iy0+ybuff+j-1
            kk = iz0+zbuff+k-1
            mat_4d(p,i,jj,kk) = global_mat(p,i,j,k)
          enddo
        enddo
      enddo
    enddo
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
!isend&recv
  if(myrank/=0) then
    do ipy = 0, npy-1
        do jpy = 0, npz-1
          recv_rank = jpy + ipy * npz
          if (myrank == recv_rank) then
              block_size = my_local(ipy+1) * mz_local(jpy+1) * mm * lx
              call MPI_Irecv(recv_buf, block_size, &
                                MPI_DOUBLE_PRECISION, 0, 900+recv_rank, MPI_COMM_WORLD, recv_request(recv_rank),ierr)        
              call MPI_Wait(recv_request(recv_rank), status, ierr) 
              do p=1,mm
                do i=1,lx
                    do j=iy0+ybuff,iyn-ybuff
                      do k=iz0+zbuff,izn-zbuff
                          jj = j-iy0-ybuff+1
                          kk = k-iz0-zbuff+1
                          mat_4d(p,i,j,k) = recv_buf(p,i,jj,kk)
                      enddo
                    enddo
                enddo
              enddo

          endif
        enddo
    enddo
  endif

  if(myrank==0) then
    do ipy = 0, npy-1
        do jpy = 0, npz-1
          recv_rank = jpy + ipy * npz
          if (recv_rank /= 0) then
              block_size = my_local(ipy+1) * mz_local(jpy+1) * mm * lx   

              do p=1,mm
                do i=1,lx
                  do j=1,iyn-iy0-2*ybuff+1
                    do k=1,izn-iz0-2*zbuff+1
                      jj = ny_local(ipy)+j
                      kk = nz_local(jpy)+k
                      send_buf(p,i,j,k) = global_mat(p,i,jj,kk)
                    enddo
                  enddo
                enddo
              enddo
              call MPI_send(send_buf , block_size, &
                                MPI_DOUBLE_PRECISION, recv_rank, 900+recv_rank, MPI_COMM_WORLD, ierr)               
          endif        
        enddo
    enddo
  endif
    
   call MPI_Barrier(MPI_COMM_WORLD, ierr)



  if (myrank == 0) then
     deallocate(global_mat,stat = ierr)
     deallocate(io_reorder,stat = ierr)
     deallocate(send_buf,stat = ierr)
  else
     deallocate(recv_buf,stat = ierr)
  endif

  end subroutine scatter_in_4d_dble


  subroutine scatter_in_fi(mat_4d,save_dir,baseFilename,istep,ix0,ixn,iy0,iyn,iz0,izn,d4i,d4e,ybuff,zbuff,commsize,npy,npz,myrank)
  
  use comvar_module
  use mpi
  implicit none
  integer :: ix0,ixn,iy0,iyn,iz0,izn,d4i,d4e,ybuff,zbuff,commsize,myrank,npy,npz
  integer :: block_size, status(MPI_STATUS_SIZE),recv_rank,iounit
  integer, ALLOCATABLE:: send_request(:),recv_request(:)
  integer :: i,j,k,p,jj,kk,istep,ipy,jpy,pp
  real*8, ALLOCATABLE :: global_mat(:,:,:,:),io_reorder(:,:,:,:),recv_buf(:,:,:,:),send_buf(:,:,:,:)
  real*8 :: mat_4d(ix0:ixn,iy0:iyn,iz0:izn,d4i:d4e)

  character(len=:), allocatable :: baseFilename, numStr, finalFilename
  character(len=*), intent(in) :: save_dir
  integer, parameter :: MaxIntLength = 20
  integer :: unit_number,rec_length

  character(len=5) :: commstr
  character(len=20) :: filename
  character(len=80) :: header 
  allocate(character(len=MaxIntLength) :: numStr)
  ALLOCATE(send_request(1:commsize-1))
  ALLOCATE(recv_request(1:commsize-1))

      
  if (myrank == 0) then
     ALLOCATE(send_buf(1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1,1:d4e-d4i+1))
     ALLOCATE(global_mat(1:lx,1:ly,1:lz,1:d4e-d4i+1))
     ALLOCATE(io_reorder(1:d4e-d4i+1,1:lz,1:ly,1:lx))
  else
     allocate(recv_buf(1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1,1:d4e-d4i+1)) 
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (myrank == 0) then
      rec_length = lx * ly * lz * (d4e-d4i+1) * 2
      WRITE(numStr, '(I9.9)') istep
      WRITE(commstr, '(I4.4)') commsize
      filename = './Results/'//baseFilename//numStr//'.raw'
      finalFilename = TRIM(save_dir) // TRIM(baseFilename) // trim(numStr) //'_'// trim(commStr) // '.raw'
      open(unit=iounit, file=finalFilename, form='unformatted', access='direct', &
             recl=rec_length, status='old', iostat=ierr)
      if (ierr /= 0) then
        print*, 'Error opening file:', finalFilename
        stop
      end if
      read(iounit, rec=1, iostat=ierr) io_reorder
      if (ierr /= 0) then
        print*, 'Error reading file:', finalFilename
        stop
      end if
      close(iounit)

      do i = 1,lx
        do j = 1,ly
          do k = 1,lz
            do p = 1,d4e-d4i+1
              global_mat(i,j,k,p) = io_reorder(p,k,j,i)
            enddo
          enddo
        enddo
      enddo
      print *,'Input file loaded: ',finalFilename

  endif


  if (myrank == 0) then
    do p=1,d4e-d4i+1
      do i=1,mx
        do j=1,my
          do k=1,mz
            jj = iy0+ybuff+j-1
            kk = iz0+zbuff+k-1
            pp = d4i-1+p
            mat_4d(i,jj,kk,pp) = global_mat(i,j,k,p)
          enddo
        enddo
      enddo
    enddo
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
!isend&recv
  if(myrank/=0) then
    do ipy = 0, npy-1
        do jpy = 0, npz-1
          recv_rank = jpy + ipy * npz
          if (myrank == recv_rank) then
              block_size = my_local(ipy+1) * mz_local(jpy+1) * (d4e-d4i+1) * lx
              call MPI_Irecv(recv_buf, block_size, &
                                MPI_DOUBLE_PRECISION, 0, 1000+recv_rank, MPI_COMM_WORLD, recv_request(recv_rank),ierr)        
              call MPI_Wait(recv_request(recv_rank), status, ierr) 
              do p=d4i,d4e
                do i=1,lx
                    do j=iy0+ybuff,iyn-ybuff
                      do k=iz0+zbuff,izn-zbuff
                          jj = j-iy0-ybuff+1
                          kk = k-iz0-zbuff+1
                          pp = p-d4i+1
                          mat_4d(i,j,k,p) = recv_buf(i,jj,kk,pp)
                      enddo
                    enddo
                enddo
              enddo

          endif
        enddo
    enddo
  endif

  if(myrank==0) then
    do ipy = 0, npy-1
        do jpy = 0, npz-1
          recv_rank = jpy + ipy * npz
          if (recv_rank /= 0) then
              block_size = my_local(ipy+1) * mz_local(jpy+1) * (d4e-d4i+1) * lx   

              do p=1,d4e-d4i+1
                do i=1,lx
                  do j=1,iyn-iy0-2*ybuff+1
                    do k=1,izn-iz0-2*zbuff+1
                      jj = ny_local(ipy)+j
                      kk = nz_local(jpy)+k
                      send_buf(i,j,k,p) = global_mat(i,jj,kk,p)
                    enddo
                  enddo
                enddo
              enddo
              call MPI_send(send_buf , block_size, &
                                MPI_DOUBLE_PRECISION, recv_rank, 1000+recv_rank, MPI_COMM_WORLD, ierr)               
          endif        
        enddo
    enddo
  endif
    
   call MPI_Barrier(MPI_COMM_WORLD, ierr)



  if (myrank == 0) then
     deallocate(global_mat,stat = ierr)
     deallocate(io_reorder,stat = ierr)
     deallocate(send_buf,stat = ierr)
  else
     deallocate(recv_buf,stat = ierr)
  endif

  end subroutine scatter_in_fi



  subroutine scatter_in_gi(mat_5d,save_dir,baseFilename,istep,mm,ix0,ixn,iy0,iyn,iz0,izn,d5i,d5e,ybuff,zbuff,commsize,npy,npz,myrank)
  
  use comvar_module
  use mpi
  implicit none
  integer :: mm,ix0,ixn,iy0,iyn,iz0,izn,d5i,d5e,ybuff,zbuff,commsize,myrank,npy,npz
  integer :: block_size, status(MPI_STATUS_SIZE),recv_rank,iounit
  integer, ALLOCATABLE:: send_request(:),recv_request(:)
  integer :: i,j,k,p,jj,kk,istep,ipy,jpy,pp,s
  real*8, ALLOCATABLE :: global_mat(:,:,:,:,:),io_reorder(:,:,:,:,:),recv_buf(:,:,:,:,:),send_buf(:,:,:,:,:)
  real*8 :: mat_5d(1:mm,ix0:ixn,iy0:iyn,iz0:izn,d5i:d5e)

  character(len=:), allocatable :: baseFilename, numStr, finalFilename
  character(len=*), intent(in) :: save_dir
  integer, parameter :: MaxIntLength = 20
  integer :: unit_number,rec_length

  character(len=5) :: commstr
  character(len=20) :: filename
  character(len=80) :: header 
  allocate(character(len=MaxIntLength) :: numStr)
  ALLOCATE(send_request(1:commsize-1))
  ALLOCATE(recv_request(1:commsize-1))

      
  if (myrank == 0) then
     ALLOCATE(send_buf(1:mm,1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1,1:d5e-d5i+1))
     ALLOCATE(global_mat(1:mm,1:lx,1:ly,1:lz,1:d5e-d5i+1))
     ALLOCATE(io_reorder(1:d5e-d5i+1,1:lz,1:ly,1:lx,1:mm))
  else
     allocate(recv_buf(1:mm,1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1,1:d5e-d5i+1)) 
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  if (myrank == 0) then
      rec_length = lx * ly * lz * (d5e-d5i+1) * mm * 2
      WRITE(numStr, '(I9.9)') istep
      WRITE(commstr, '(I4.4)') commsize
      filename = './Results/'//baseFilename//numStr//'.raw'
      finalFilename = TRIM(save_dir) // TRIM(baseFilename) // trim(numStr) //'_'// trim(commStr) // '.raw'
      open(unit=iounit, file=finalFilename, form='unformatted', access='direct', &
             recl=rec_length, status='old', iostat=ierr)
      if (ierr /= 0) then
        print*, 'Error opening file:', finalFilename
        stop
      end if
      read(iounit, rec=1, iostat=ierr) io_reorder
      if (ierr /= 0) then
        print*, 'Error reading file:', finalFilename
        stop
      end if
      close(iounit)
      do s = 1,mm
        do i = 1,lx
          do j = 1,ly
            do k = 1,lz
              do p = 1,d5e-d5i+1
                global_mat(s,i,j,k,p) = io_reorder(p,k,j,i,s)
              enddo
            enddo
          enddo
        enddo
      enddo
      print *,'Input file loaded: ',finalFilename

  endif


  if (myrank == 0) then
    do s = 1,mm
      do p=1,d5e-d5i+1
        do i=1,mx
          do j=1,my
            do k=1,mz
              jj = iy0+ybuff+j-1
              kk = iz0+zbuff+k-1
              pp = d5i-1+p
              mat_5d(s,i,jj,kk,pp) = global_mat(s,i,j,k,p)
            enddo
          enddo
        enddo
      enddo
    enddo
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
!isend&recv
  if(myrank/=0) then
    do ipy = 0, npy-1
        do jpy = 0, npz-1
          recv_rank = jpy + ipy * npz
          if (myrank == recv_rank) then
              block_size = my_local(ipy+1) * mz_local(jpy+1) * (d5e-d5i+1) * mm * lx
              call MPI_Irecv(recv_buf, block_size, &
                                MPI_DOUBLE_PRECISION, 0, 1300+recv_rank, MPI_COMM_WORLD, recv_request(recv_rank),ierr)        
              call MPI_Wait(recv_request(recv_rank), status, ierr) 

              do s = 1,mm
                do p=d5i,d5e
                  do i=1,lx
                      do j=iy0+ybuff,iyn-ybuff
                        do k=iz0+zbuff,izn-zbuff
                            jj = j-iy0-ybuff+1
                            kk = k-iz0-zbuff+1
                            pp = p-d5i+1
                            mat_5d(s,i,j,k,p) = recv_buf(s,i,jj,kk,pp)
                        enddo
                      enddo
                  enddo
                enddo
              enddo

          endif
        enddo
    enddo
  endif

  if(myrank==0) then
    do ipy = 0, npy-1
        do jpy = 0, npz-1
          recv_rank = jpy + ipy * npz
          if (recv_rank /= 0) then
              block_size = my_local(ipy+1) * mz_local(jpy+1) * (d5e-d5i+1) * mm *  lx   

              do s = 1,mm
                do p=1,d5e-d5i+1
                  do i=1,lx
                    do j=1,iyn-iy0-2*ybuff+1
                      do k=1,izn-iz0-2*zbuff+1
                        jj = ny_local(ipy)+j
                        kk = nz_local(jpy)+k
                        send_buf(s,i,j,k,p) = global_mat(s,i,jj,kk,p)
                      enddo
                    enddo
                  enddo
                enddo
              enddo

              call MPI_send(send_buf , block_size, &
                                MPI_DOUBLE_PRECISION, recv_rank, 1300+recv_rank, MPI_COMM_WORLD, ierr)               
          endif        
        enddo
    enddo
  endif
    
   call MPI_Barrier(MPI_COMM_WORLD, ierr)



  if (myrank == 0) then
     deallocate(global_mat,stat = ierr)
     deallocate(io_reorder,stat = ierr)
     deallocate(send_buf,stat = ierr)
  else
     deallocate(recv_buf,stat = ierr)
  endif

  end subroutine scatter_in_gi



  subroutine scatter_in_4d_logic(mat_4d,save_dir,baseFilename,istep,mm,ix0,ixn,iy0,iyn,iz0,izn,ybuff,zbuff,commsize,npy,npz,myrank)
  
  use comvar_module
  use mpi
  implicit none
  integer :: mm,ix0,ixn,iy0,iyn,iz0,izn,ybuff,zbuff,commsize,myrank,npy,npz
  integer :: block_size, status(MPI_STATUS_SIZE),recv_rank,iounit
  integer, ALLOCATABLE:: send_request(:),recv_request(:)
  integer :: i,j,k,p,jj,kk,istep,ipy,jpy
  logical, ALLOCATABLE :: global_mat(:,:,:,:),io_reorder(:,:,:,:),recv_buf(:,:,:,:),send_buf(:,:,:,:)
  logical :: mat_4d(1:mm,ix0:ixn,iy0:iyn,iz0:izn)

  character(len=:), allocatable :: baseFilename, numStr, finalFilename
  character(len=*), intent(in) :: save_dir
  integer, parameter :: MaxIntLength = 20
  integer :: unit_number,rec_length

  character(len=5) :: commstr
  character(len=20) :: filename
  character(len=80) :: header 
  allocate(character(len=MaxIntLength) :: numStr)
  ALLOCATE(send_request(1:commsize-1))
  ALLOCATE(recv_request(1:commsize-1))

      
  if (myrank == 0) then
     ALLOCATE(send_buf(mm,1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1))
     ALLOCATE(global_mat(mm,1:lx,1:ly,1:lz))
     ALLOCATE(io_reorder(1:lz,1:ly,1:lx,1:mm))
  else
     allocate(recv_buf(1:mm,1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1)) 
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (myrank == 0) then
      rec_length = lx * ly * lz * mm
      WRITE(numStr, '(I9.9)') istep
      WRITE(commstr, '(I4.4)') commsize
      filename = './Results/'//baseFilename//numStr//'.raw'
      finalFilename = TRIM(save_dir) // TRIM(baseFilename) // trim(numStr) //'_'// trim(commStr) // '.raw'
      open(unit=iounit, file=finalFilename, form='unformatted', access='direct', &
             recl=rec_length, status='old', iostat=ierr)
      if (ierr /= 0) then
        print*, 'Error opening file:', finalFilename
        stop
      end if
      read(iounit, rec=1, iostat=ierr) io_reorder
      if (ierr /= 0) then
        print*, 'Error reading file:', finalFilename
        stop
      end if
      close(iounit)

      do i = 1,lx
        do j = 1,ly
          do k = 1,lz
            do p = 1,mm
              global_mat(p,i,j,k) = io_reorder(k,j,i,p)
            enddo
          enddo
        enddo
      enddo
      print *,'Input file loaded: ',finalFilename

  endif


  if (myrank == 0) then
    do p=1,mm
      do i=1,mx
        do j=1,my
          do k=1,mz
            jj = iy0+ybuff+j-1
            kk = iz0+zbuff+k-1
            mat_4d(p,i,jj,kk) = global_mat(p,i,j,k)
          enddo
        enddo
      enddo
    enddo
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
!isend&recv
  if(myrank/=0) then
    do ipy = 0, npy-1
        do jpy = 0, npz-1
          recv_rank = jpy + ipy * npz
          if (myrank == recv_rank) then
              block_size = my_local(ipy+1) * mz_local(jpy+1) * mm * lx
              call MPI_Irecv(recv_buf, block_size, &
                                MPI_LOGICAL, 0, 1000+recv_rank, MPI_COMM_WORLD, recv_request(recv_rank),ierr)        
              call MPI_Wait(recv_request(recv_rank), status, ierr) 
              do p=1,mm
                do i=1,lx
                    do j=iy0+ybuff,iyn-ybuff
                      do k=iz0+zbuff,izn-zbuff
                          jj = j-iy0-ybuff+1
                          kk = k-iz0-zbuff+1
                          mat_4d(p,i,j,k) = recv_buf(p,i,jj,kk)
                      enddo
                    enddo
                enddo
              enddo

          endif
        enddo
    enddo
  endif

  if(myrank==0) then
    do ipy = 0, npy-1
        do jpy = 0, npz-1
          recv_rank = jpy + ipy * npz
          if (recv_rank /= 0) then
              block_size = my_local(ipy+1) * mz_local(jpy+1) * mm * lx   

              do p=1,mm
                do i=1,lx
                  do j=1,iyn-iy0-2*ybuff+1
                    do k=1,izn-iz0-2*zbuff+1
                      jj = ny_local(ipy)+j
                      kk = nz_local(jpy)+k
                      send_buf(p,i,j,k) = global_mat(p,i,jj,kk)
                    enddo
                  enddo
                enddo
              enddo
              call MPI_send(send_buf , block_size, &
                                MPI_LOGICAL, recv_rank, 1000+recv_rank, MPI_COMM_WORLD, ierr)               
          endif        
        enddo
    enddo
  endif
    
   call MPI_Barrier(MPI_COMM_WORLD, ierr)



  if (myrank == 0) then
     deallocate(global_mat,stat = ierr)
     deallocate(io_reorder,stat = ierr)
     deallocate(send_buf,stat = ierr)
  else
     deallocate(recv_buf,stat = ierr)
  endif

  end subroutine scatter_in_4d_logic


  subroutine scatter_in_4d_logic2(mat_4d,save_dir,baseFilename,mm,ix0,ixn,iy0,iyn,iz0,izn,ybuff,zbuff,commsize,npy,npz,myrank)
  
  use comvar_module
  use mpi
  implicit none
  integer :: mm,ix0,ixn,iy0,iyn,iz0,izn,ybuff,zbuff,commsize,myrank,npy,npz
  integer :: block_size, status(MPI_STATUS_SIZE),recv_rank,iounit
  integer, ALLOCATABLE:: send_request(:),recv_request(:)
  integer :: i,j,k,p,ii,jj,kk,ipy,jpy
  logical, ALLOCATABLE :: global_mat(:,:,:,:),global_inner(:,:,:),io_reorder(:,:,:,:),recv_buf(:,:,:,:),send_buf(:,:,:,:)
  logical :: mat_4d(1:mm,ix0:ixn,iy0:iyn,iz0:izn)

  character(len=:), allocatable :: baseFilename, numStr, finalFilename
  character(len=*), intent(in) :: save_dir
  integer, parameter :: MaxIntLength = 20
  integer :: unit_number,rec_length
  integer :: pdx_i=3,pdx_o=3,pdy_u=6,pdy_d=5,pdz_u=6,pdz_d=5

  character(len=5) :: commstr
  character(len=80) :: header 
  allocate(character(len=MaxIntLength) :: numStr)
  ALLOCATE(send_request(1:commsize-1))
  ALLOCATE(recv_request(1:commsize-1))

      
  if (myrank == 0) then
     ALLOCATE(send_buf(mm,1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1))
     ALLOCATE(global_mat(mm,1:lx,1:ly,1:lz))
     ALLOCATE(io_reorder(1:lz,1:ly,1:lx,1:mm))
     allocate( global_inner(1:lx-pdx_i-pdx_o,1:ly-pdy_u-pdy_d,1:lz-pdz_u-pdz_d), STAT=istat)
  else
     allocate(recv_buf(1:mm,1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1)) 
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (myrank == 0) then
      rec_length = lx * ly * lz * mm
      WRITE(commstr, '(I4.4)') commsize
      finalFilename = TRIM(save_dir) //'../'// TRIM(baseFilename) // '.raw'
      open(unit=iounit, file=finalFilename, form='unformatted', access='direct', &
             recl=rec_length, status='old', iostat=ierr)
      if (ierr /= 0) then
        print*, 'Error opening file:', finalFilename
        stop
      end if
      read(iounit, rec=1, iostat=ierr) io_reorder
      if (ierr /= 0) then
        print*, 'Error reading file:', finalFilename
        stop
      end if
      close(iounit)

      do i = 1,lx
        do j = 1,ly
          do k = 1,lz
            do p = 1,mm
              global_mat(p,i,j,k) = io_reorder(k,j,i,p)
            enddo
          enddo
        enddo
      enddo
      print *,'Input file loaded: ',finalFilename

  endif


  if (myrank == 0) then
    do p=1,mm
      do i=1,mx
        do j=1,my
          do k=1,mz
            jj = iy0+ybuff+j-1
            kk = iz0+zbuff+k-1
            mat_4d(p,i,jj,kk) = global_mat(p,i,j,k)
          enddo
        enddo
      enddo
    enddo
  endif


  call MPI_Barrier(MPI_COMM_WORLD, ierr)
!isend&recv
  if(myrank/=0) then
    do ipy = 0, npy-1
        do jpy = 0, npz-1
          recv_rank = jpy + ipy * npz
          if (myrank == recv_rank) then
              block_size = my_local(ipy+1) * mz_local(jpy+1) * mm * lx
              call MPI_Irecv(recv_buf, block_size, &
                                MPI_LOGICAL, 0, 1000+recv_rank, MPI_COMM_WORLD, recv_request(recv_rank),ierr)        
              call MPI_Wait(recv_request(recv_rank), status, ierr) 
              do p=1,mm
                do i=1,lx
                    do j=iy0+ybuff,iyn-ybuff
                      do k=iz0+zbuff,izn-zbuff
                          jj = j-iy0-ybuff+1
                          kk = k-iz0-zbuff+1
                          mat_4d(p,i,j,k) = recv_buf(p,i,jj,kk)
                      enddo
                    enddo
                enddo
              enddo

          endif
        enddo
    enddo
  endif

  if(myrank==0) then
    do ipy = 0, npy-1
        do jpy = 0, npz-1
          recv_rank = jpy + ipy * npz
          if (recv_rank /= 0) then
              block_size = my_local(ipy+1) * mz_local(jpy+1) * mm * lx   

              do p=1,mm
                do i=1,lx
                  do j=1,iyn-iy0-2*ybuff+1
                    do k=1,izn-iz0-2*zbuff+1
                      jj = ny_local(ipy)+j
                      kk = nz_local(jpy)+k
                      send_buf(p,i,j,k) = global_mat(p,i,jj,kk)
                    enddo
                  enddo
                enddo
              enddo
              call MPI_send(send_buf , block_size, &
                                MPI_LOGICAL, recv_rank, 1000+recv_rank, MPI_COMM_WORLD, ierr)               
          endif        
        enddo
    enddo
  endif
    

 
   call MPI_Barrier(MPI_COMM_WORLD, ierr)



  if (myrank == 0) then
     deallocate(global_mat,stat = ierr)
     deallocate(io_reorder,stat = ierr)
     deallocate(send_buf,stat = ierr)
  else
     deallocate(recv_buf,stat = ierr)
  endif

  end subroutine scatter_in_4d_logic2


!   subroutine scatter_in_4d_logic3(mat_4d,wall_global,mm,ix0,ixn,iy0,iyn,iz0,izn,ybuff,zbuff,commsize,npy,npz,myrank)
  
!   use comvar_module
!   use mpi
!   implicit none
!   integer :: mm,ix0,ixn,iy0,iyn,iz0,izn,ybuff,zbuff,commsize,myrank,npy,npz
!   integer :: block_size, status(MPI_STATUS_SIZE),recv_rank
!   integer, ALLOCATABLE:: send_request(:),recv_request(:)
!   integer :: i,j,k,p,ii,jj,kk,ipy,jpy
!   logical, ALLOCATABLE :: global_mat(:,:,:,:),global_inner(:,:,:),io_reorder(:,:,:,:),recv_buf(:,:,:,:),send_buf(:,:,:,:)
!   logical :: mat_4d(1:mm,ix0:ixn,iy0:iyn,iz0:izn)

!   character(len=:), allocatable :: baseFilename, numStr, finalFilename
!   character(len=*), intent(in) :: save_dir
!   integer, parameter :: MaxIntLength = 20
!   integer :: unit_number,rec_length
!   integer :: pdx_i=3,pdx_o=3,pdy_u=6,pdy_d=5,pdz_u=6,pdz_d=5

!   character(len=5) :: commstr
!   character(len=80) :: header 
!   allocate(character(len=MaxIntLength) :: numStr)
!   ALLOCATE(send_request(1:commsize-1))
!   ALLOCATE(recv_request(1:commsize-1))

      
!   if (myrank == 0) then
!      ALLOCATE(send_buf(mm,1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1))
!      ALLOCATE(global_mat(mm,1:lx,1:ly,1:lz))
!      ALLOCATE(io_reorder(1:lz,1:ly,1:lx,1:mm))
!      allocate( global_inner(1:lx-pdx_i-pdx_o,1:ly-pdy_u-pdy_d,1:lz-pdz_u-pdz_d), STAT=istat)
!   else
!      allocate(recv_buf(1:mm,1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1)) 
!   endif

!   call MPI_Barrier(MPI_COMM_WORLD, ierr)


! ! if (myrank == 0) then
! !     finalFilename = TRIM(save_dir) // TRIM(baseFilename) // '.dat'
! !     open(unit=iounit, file=finalFilename, form='formatted', status='old', iostat=ierr)
! !     if (ierr /= 0) then
! !         print*, 'Error opening file:', finalFilename
! !         stop
! !     end if

! !     do k = 1, lz
! !         do j = 1, ly
! !             do i = 1, lx
! !                 do p = 1, mm
! !                     read(iounit, '(L1)', iostat=ierr) io_reorder(k, j, i, p)
! !                     if (ierr /= 0) then
! !                         print*, 'Error reading file:', finalFilename
! !                         stop
! !                     end if
! !                 end do
! !             end do
! !         end do
! !     end do
! !     close(iounit)

! !     do i = 1, lx
! !         do j = 1, ly
! !             do k = 1, lz
! !                 do p = 1, mm
! !                     global_mat(p, i, j, k) = io_reorder(k, j, i, p)
! !                 end do
! !             end do
! !         end do
! !     end do
! !     print *, 'Input file loaded: ', finalFilename
! ! endif

!   if (myrank == 0) then
!     global_inner=.false.
!     do i=1,lx-pdx_i-pdx_o
!       do j=1,ly-pdy_u-pdy_d
!         do k=1,lz-pdz_u-pdz_d
!           if(i<28.and.j>=20.and.j<93) global_inner(i,j,k) = .true.!1
!           if(i<73.and.i>=49.and.j<93) global_inner(i,j,k) = .true.!2
!           if(i<248.and.i>=223.and.j<93) global_inner(i,j,k) = .true.!3
!           if(i>=267.and.j>=20.and.j<93) global_inner(i,j,k) = .true.!4
!           if(i>=92.and.i<130.and.j>=66) global_inner(i,j,k) = .true.!5
!           if(i>=135.and.i<181.and.j>=66) global_inner(i,j,k) = .true.!6
!           if(i>=92.and.i<106.and.j>=20) global_inner(i,j,k) = .true.!7
!           if(i>=111.and.i<130.and.j>=20) global_inner(i,j,k) = .true.!8
!           if(i>=135.and.i<157.and.j>=20) global_inner(i,j,k) = .true.!9
!           if(i>=167.and.i<181.and.j>=20) global_inner(i,j,k) = .true.!10
!           if(i>=191.and.i<204.and.j>=20) global_inner(i,j,k) = .true.!11
!         enddo
!       enddo
!     enddo

!     global_mat = .false.
!     global_mat(2,:,:,:)=.true.
!     do i=1,lx-pdx_i-pdx_o
!       do j=1,ly-pdy_u-pdy_d
!         do k=1,lz-pdz_u-pdz_d
!           ii=pdx_i+i
!           jj=pdy_d+j
!           kk=pdz_d+k
!           global_mat(2,ii,jj,kk) = global_inner(i,j,k)
!         enddo
!       enddo
!     enddo

!     global_mat(:,:pdx_i,:,:) = .false.
!     global_mat(:,lx-pdx_o+1:,:,:) = .false.

!   endif

!   print*,'current rank=',myrank
!   if (myrank == 0) then
!     do p=1,mm
!       do i=1,mx
!         do j=1,my
!           do k=1,mz
!             jj = iy0+ybuff+j-1
!             kk = iz0+zbuff+k-1
!             mat_4d(p,i,jj,kk) = global_mat(p,i,j,k)
!           enddo
!         enddo
!       enddo
!     enddo
!   endif


!   call MPI_Barrier(MPI_COMM_WORLD, ierr)
! ! !isend&recv
! !   if(myrank/=0) then
! !     do ipy = 0, npy-1
! !         do jpy = 0, npz-1
! !           recv_rank = jpy + ipy * npz
! !           if (myrank == recv_rank) then
! !               block_size = my_local(ipy+1) * mz_local(jpy+1) * mm * lx
! !               call MPI_Irecv(recv_buf, block_size, &
! !                                 MPI_LOGICAL, 0, 1000+recv_rank, MPI_COMM_WORLD, recv_request(recv_rank),ierr)        
! !               call MPI_Wait(recv_request(recv_rank), status, ierr) 
! !               do p=1,mm
! !                 do i=1,lx
! !                     do j=iy0+ybuff,iyn-ybuff
! !                       do k=iz0+zbuff,izn-zbuff
! !                           jj = j-iy0-ybuff+1
! !                           kk = k-iz0-zbuff+1
! !                           mat_4d(p,i,j,k) = recv_buf(p,i,jj,kk)
! !                       enddo
! !                     enddo
! !                 enddo
! !               enddo

! !           endif
! !         enddo
! !     enddo
! !   endif

! !   if(myrank==0) then
! !     do ipy = 0, npy-1
! !         do jpy = 0, npz-1
! !           recv_rank = jpy + ipy * npz
! !           if (recv_rank /= 0) then
! !               block_size = my_local(ipy+1) * mz_local(jpy+1) * mm * lx   

! !               do p=1,mm
! !                 do i=1,lx
! !                   do j=1,iyn-iy0-2*ybuff+1
! !                     do k=1,izn-iz0-2*zbuff+1
! !                       jj = ny_local(ipy)+j
! !                       kk = nz_local(jpy)+k
! !                       send_buf(p,i,j,k) = global_mat(p,i,jj,kk)
! !                     enddo
! !                   enddo
! !                 enddo
! !               enddo
! !               call MPI_send(send_buf , block_size, &
! !                                 MPI_LOGICAL, recv_rank, 1000+recv_rank, MPI_COMM_WORLD, ierr)               
! !           endif        
! !         enddo
! !     enddo
! !   endif
    

!    call MPI_Barrier(MPI_COMM_WORLD, ierr)



!   if (myrank == 0) then
!      deallocate(global_mat,stat = ierr)
!      deallocate(io_reorder,stat = ierr)
!      deallocate(send_buf,stat = ierr)
!   else
!      deallocate(recv_buf,stat = ierr)
!   endif

!   end subroutine scatter_in_4d_logic3


!   subroutine scatter_in_3d_int_to_logic(mat_3d_logic,baseFilename,ix0,ixn,iy0,iyn,iz0,izn,ybuff,zbuff,commsize,npy,npz,myrank)
  
!   use comvar_module
!   use mpi
!   implicit none
!   integer :: ix0,ixn,iy0,iyn,iz0,izn,ybuff,zbuff,commsize,myrank,npy,npz
!   integer :: block_size, status(MPI_STATUS_SIZE),recv_rank,iounit
!   integer, ALLOCATABLE:: send_request(:),recv_request(:)
!   integer :: i,j,k,p,jj,kk,ipy,jpy
!   integer, ALLOCATABLE :: global_mat(:,:,:),io_reorder(:,:,:),recv_buf(:,:,:),send_buf(:,:,:),mat_3d(:,:,:)
!   logical :: mat_3d_logic(ix0:ixn,iy0:iyn,iz0:izn)

!   character(len=:), allocatable :: baseFilename, finalFilename
!   integer, parameter :: MaxIntLength = 20
!   integer :: unit_number,rec_length

!   character(len=20) :: filename
!   character(len=80) :: header 
!   ALLOCATE(send_request(1:commsize-1))
!   ALLOCATE(recv_request(1:commsize-1))
!   ALLOCATE(mat_3d(ix0:ixn,iy0:iyn,iz0:izn))
      
!   if (myrank == 0) then
!      ALLOCATE(send_buf(1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1))
!      ALLOCATE(global_mat(1:lx,1:ly,1:lz))
!      ALLOCATE(io_reorder(1:lz,1:ly,1:lx))
!   else
!      allocate(recv_buf(1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1)) 
!   endif

!   call MPI_Barrier(MPI_COMM_WORLD, ierr)
!   if (myrank == 0) then
!       rec_length = lx * ly * lz
!       filename = './'//baseFilename//'.raw'
!       finalFilename = TRIM(baseFilename) // '.raw'
!       open(unit=iounit, file=finalFilename, form='unformatted', access='direct', &
!              recl=rec_length, status='old', iostat=ierr)
!       if (ierr /= 0) then
!         print*, 'Error opening file:', finalFilename
!         stop
!       end if
!       read(iounit, rec=1, iostat=ierr) io_reorder
!       if (ierr /= 0) then
!         print*, 'Error reading file:', finalFilename
!         stop
!       end if
!       close(iounit)

!       do i = 1,lx
!         do j = 1,ly
!           do k = 1,lz
!               global_mat(i,j,k) = io_reorder(k,j,i)
!           enddo
!         enddo
!       enddo


!       print *,'Input file loaded: ',finalFilename

!   endif


!   if (myrank == 0) then
!     do i=1,mx
!       do j=1,my
!         do k=1,mz
!           jj = iy0+ybuff+j-1
!           kk = iz0+zbuff+k-1
!           mat_3d(i,jj,kk) = global_mat(i,j,k)
!         enddo
!       enddo
!     enddo
!   endif

!   call MPI_Barrier(MPI_COMM_WORLD, ierr)
! !isend&recv
!   if(myrank/=0) then
!     do ipy = 0, npy-1
!         do jpy = 0, npz-1
!           recv_rank = jpy + ipy * npz
!           if (myrank == recv_rank) then
!               block_size = my_local(ipy+1) * mz_local(jpy+1) * lx
!               call MPI_Irecv(recv_buf, block_size, &
!                                 MPI_INTEGER, 0, 1000+recv_rank, MPI_COMM_WORLD, recv_request(recv_rank),ierr)        
!               call MPI_Wait(recv_request(recv_rank), status, ierr) 
!               do i=1,lx
!                   do j=iy0+ybuff,iyn-ybuff
!                     do k=iz0+zbuff,izn-zbuff
!                         jj = j-iy0-ybuff+1
!                         kk = k-iz0-zbuff+1
!                         mat_3d(i,j,k) = recv_buf(i,jj,kk)
!                     enddo
!                   enddo
!               enddo

!           endif
!         enddo
!     enddo
!   endif

!   if(myrank==0) then
!     do ipy = 0, npy-1
!         do jpy = 0, npz-1
!           recv_rank = jpy + ipy * npz
!           if (recv_rank /= 0) then
!               block_size = my_local(ipy+1) * mz_local(jpy+1) * lx   
!               do i=1,lx
!                 do j=1,iyn-iy0-2*ybuff+1
!                   do k=1,izn-iz0-2*zbuff+1
!                     jj = ny_local(ipy)+j
!                     kk = nz_local(jpy)+k
!                     send_buf(i,j,k) = global_mat(i,jj,kk)
!                   enddo
!                 enddo
!               enddo
!               call MPI_send(send_buf , block_size, &
!                                 MPI_INTEGER, recv_rank, 1000+recv_rank, MPI_COMM_WORLD, ierr)               
!           endif        
!         enddo
!     enddo
!   endif
    
!   call MPI_Barrier(MPI_COMM_WORLD, ierr)


!   do i=1,mx
!     do j=1,my
!       do k=1,mz
!         if(mat_3d(i,j,k)==1) then
!           mat_3d_logic(i,j,k) = .true.
!         else
!           mat_3d_logic(i,j,k) = .false.
!         endif
!       enddo
!     enddo
!   enddo

!   call MPI_Barrier(MPI_COMM_WORLD, ierr)

!   if (myrank == 0) then
!      deallocate(global_mat,stat = ierr)
!      deallocate(io_reorder,stat = ierr)
!      deallocate(send_buf,stat = ierr)
!   else
!      deallocate(recv_buf,stat = ierr)
!   endif
!   deallocate(mat_3d,stat = ierr)

!   end subroutine scatter_in_3d_int_to_logic



!   subroutine scatter_in_3d_int(mat_3d,baseFilename,ix0,ixn,iy0,iyn,iz0,izn,ybuff,zbuff,commsize,npy,npz,myrank)
  
!   use comvar_module
!   use mpi
!   implicit none
!   integer :: ix0,ixn,iy0,iyn,iz0,izn,ybuff,zbuff,commsize,myrank,npy,npz
!   integer :: block_size, status(MPI_STATUS_SIZE),recv_rank,iounit
!   integer, ALLOCATABLE:: send_request(:),recv_request(:)
!   integer :: i,j,k,p,jj,kk,ipy,jpy
!   integer, ALLOCATABLE :: global_mat(:,:,:),io_reorder(:,:,:),recv_buf(:,:,:),send_buf(:,:,:)
!   integer :: mat_3d(ix0:ixn,iy0:iyn,iz0:izn)

!   character(len=:), allocatable :: baseFilename, finalFilename
!   integer, parameter :: MaxIntLength = 20
!   integer :: unit_number,rec_length

!   character(len=20) :: filename
!   character(len=80) :: header 
!   ALLOCATE(send_request(1:commsize-1))
!   ALLOCATE(recv_request(1:commsize-1))

      
!   if (myrank == 0) then
!      ALLOCATE(send_buf(1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1))
!      ALLOCATE(global_mat(1:lx,1:ly,1:lz))
!      ALLOCATE(io_reorder(1:lz,1:ly,1:lx))
!   else
!      allocate(recv_buf(1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1)) 
!   endif

!   call MPI_Barrier(MPI_COMM_WORLD, ierr)
!   if (myrank == 0) then
!       rec_length = lx * ly * lz
!       filename = './'//baseFilename//'.raw'
!       finalFilename = TRIM(baseFilename) // '.raw'
!       open(unit=iounit, file=finalFilename, form='unformatted', access='direct', &
!              recl=rec_length, status='old', iostat=ierr)
!       if (ierr /= 0) then
!         print*, 'Error opening file:', finalFilename
!         stop
!       end if
!       read(iounit, rec=1, iostat=ierr) io_reorder
!       if (ierr /= 0) then
!         print*, 'Error reading file:', finalFilename
!         stop
!       end if
!       close(iounit)

!       do i = 1,lx
!         do j = 1,ly
!           do k = 1,lz
!               global_mat(i,j,k) = io_reorder(k,j,i)
!           enddo
!         enddo
!       enddo


!       print *,'Input file loaded: ',finalFilename

!   endif


!   if (myrank == 0) then
!     do i=1,mx
!       do j=1,my
!         do k=1,mz
!           jj = iy0+ybuff+j-1
!           kk = iz0+zbuff+k-1
!           mat_3d(i,jj,kk) = global_mat(i,j,k)
!         enddo
!       enddo
!     enddo
!   endif

!   call MPI_Barrier(MPI_COMM_WORLD, ierr)
! !isend&recv
!   if(myrank/=0) then
!     do ipy = 0, npy-1
!         do jpy = 0, npz-1
!           recv_rank = jpy + ipy * npz
!           if (myrank == recv_rank) then
!               block_size = my_local(ipy+1) * mz_local(jpy+1) * lx
!               call MPI_Irecv(recv_buf, block_size, &
!                                 MPI_INTEGER, 0, 1000+recv_rank, MPI_COMM_WORLD, recv_request(recv_rank),ierr)        
!               call MPI_Wait(recv_request(recv_rank), status, ierr) 
!               do i=1,lx
!                   do j=iy0+ybuff,iyn-ybuff
!                     do k=iz0+zbuff,izn-zbuff
!                         jj = j-iy0-ybuff+1
!                         kk = k-iz0-zbuff+1
!                         mat_3d(i,j,k) = recv_buf(i,jj,kk)
!                     enddo
!                   enddo
!               enddo

!           endif
!         enddo
!     enddo
!   endif

!   if(myrank==0) then
!     do ipy = 0, npy-1
!         do jpy = 0, npz-1
!           recv_rank = jpy + ipy * npz
!           if (recv_rank /= 0) then
!               block_size = my_local(ipy+1) * mz_local(jpy+1) * lx   
!               do i=1,lx
!                 do j=1,iyn-iy0-2*ybuff+1
!                   do k=1,izn-iz0-2*zbuff+1
!                     jj = ny_local(ipy)+j
!                     kk = nz_local(jpy)+k
!                     send_buf(i,j,k) = global_mat(i,jj,kk)
!                   enddo
!                 enddo
!               enddo
!               call MPI_send(send_buf , block_size, &
!                                 MPI_INTEGER, recv_rank, 1000+recv_rank, MPI_COMM_WORLD, ierr)               
!           endif        
!         enddo
!     enddo
!   endif
    
!   call MPI_Barrier(MPI_COMM_WORLD, ierr)


!   if (myrank == 0) then
!      deallocate(global_mat,stat = ierr)
!      deallocate(io_reorder,stat = ierr)
!      deallocate(send_buf,stat = ierr)
!   else
!      deallocate(recv_buf,stat = ierr)
!   endif

!   end subroutine scatter_in_3d_int


!   subroutine scatter_in_3d_logic(mat_3d,ix0,ixn,iy0,iyn,iz0,izn,ybuff,zbuff,commsize,npy,npz,myrank)
  
!   use comvar_module
!   use mpi
!   implicit none
!   integer :: ix0,ixn,iy0,iyn,iz0,izn,ybuff,zbuff,commsize,myrank,npy,npz
!   integer :: block_size, status(MPI_STATUS_SIZE),recv_rank,iounit
!   integer, ALLOCATABLE:: send_request(:),recv_request(:)
!   integer :: i,j,k,p,jj,kk,ipy,jpy
!   logical, ALLOCATABLE :: global_mat(:,:,:),global_inner(:,:,:),recv_buf(:,:,:),send_buf(:,:,:)
!   logical :: mat_3d(ix0:ixn,iy0:iyn,iz0:izn)
!   integer :: pdx_i=3,pdx_o=3,pdy_u=6,pdy_d=5,pdz_u=6,pdz_d=5
  
!   integer, parameter :: MaxIntLength = 20
!   integer :: unit_number,rec_length

!   character(len=80) :: header 
!   ALLOCATE(send_request(1:commsize-1))
!   ALLOCATE(recv_request(1:commsize-1))

      
!   if (myrank == 0) then
!      ALLOCATE(send_buf(1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1))
!      ALLOCATE(global_mat(1:lx,1:ly,1:lz))
!      ALLOCATE(global_inner(1:lx-pdx_i-pdx_o,1:ly-pdy_u-pdy_d,1:lz-pdz_u-pdz_d))
!   else
!      allocate(recv_buf(1:lx,1:iyn-iy0-2*ybuff+1,1:izn-iz0-2*zbuff+1)) 
!   endif

!   call MPI_Barrier(MPI_COMM_WORLD, ierr)
!   if (myrank == 0) then
!     global_inner=.false.
!     do i=1,lx
!       do j=1,ly
!         do k=1,lz
!           if(i<28.and.j>=20.and.j<93) global_inner(i,j,k) = .true.!1
!           if(i<73.and.i>=49.and.j<93) global_inner(i,j,k) = .true.!2
!           if(i<248.and.i>=223.and.j<93) global_inner(i,j,k) = .true.!3
!           if(i>=267.and.j>=20.and.j<93) global_inner(i,j,k) = .true.!4
!           if(i>=92.and.i<130.and.j>=66) global_inner(i,j,k) = .true.!5
!           if(i>=135.and.i<181.and.j>=66) global_inner(i,j,k) = .true.!6
!           if(i>=92.and.i<106.and.j>=20) global_inner(i,j,k) = .true.!7
!           if(i>=111.and.i<130.and.j>=20) global_inner(i,j,k) = .true.!8
!           if(i>=135.and.i<157.and.j>=20) global_inner(i,j,k) = .true.!9
!           if(i>=167.and.i<181.and.j>=20) global_inner(i,j,k) = .true.!10
!           if(i>=191.and.i<204.and.j>=20) global_inner(i,j,k) = .true.!11
!         enddo
!       enddo
!     enddo
!     global_mat=.true.
!     global_mat(1+pdx_i:lx-pdx_o,1+pdy_d:ly-pdy_u,1+pdz_d:lz-pdz_u) = global_inner
!     global_mat(:pdx_i,:,:) = .false.
!     global_mat(lx-pdx_o+1:,:,:) = .false.


      

!   endif


!   if (myrank == 0) then
!     do i=1,mx
!       do j=1,my
!         do k=1,mz
!           jj = iy0+ybuff+j-1
!           kk = iz0+zbuff+k-1
!           mat_3d(i,jj,kk) = global_mat(i,j,k)
!         enddo
!       enddo
!     enddo
!   endif

!   call MPI_Barrier(MPI_COMM_WORLD, ierr)
! !isend&recv
!   if(myrank/=0) then
!     do ipy = 0, npy-1
!         do jpy = 0, npz-1
!           recv_rank = jpy + ipy * npz
!           if (myrank == recv_rank) then
!               block_size = my_local(ipy+1) * mz_local(jpy+1) * lx
!               call MPI_Irecv(recv_buf, block_size, &
!                                 MPI_LOGICAL, 0, 1000+recv_rank, MPI_COMM_WORLD, recv_request(recv_rank),ierr)        
!               call MPI_Wait(recv_request(recv_rank), status, ierr) 
!               do i=1,lx
!                   do j=iy0+ybuff,iyn-ybuff
!                     do k=iz0+zbuff,izn-zbuff
!                         jj = j-iy0-ybuff+1
!                         kk = k-iz0-zbuff+1
!                         mat_3d(i,j,k) = recv_buf(i,jj,kk)
!                     enddo
!                   enddo
!               enddo

!           endif
!         enddo
!     enddo
!   endif

!   if(myrank==0) then
!     do ipy = 0, npy-1
!         do jpy = 0, npz-1
!           recv_rank = jpy + ipy * npz
!           if (recv_rank /= 0) then
!               block_size = my_local(ipy+1) * mz_local(jpy+1) * lx   
!               do i=1,lx
!                 do j=1,iyn-iy0-2*ybuff+1
!                   do k=1,izn-iz0-2*zbuff+1
!                     jj = ny_local(ipy)+j
!                     kk = nz_local(jpy)+k
!                     send_buf(i,j,k) = global_mat(i,jj,kk)
!                   enddo
!                 enddo
!               enddo
!               call MPI_send(send_buf , block_size, &
!                                 MPI_LOGICAL, recv_rank, 1000+recv_rank, MPI_COMM_WORLD, ierr)               
!           endif        
!         enddo
!     enddo
!   endif
    
!   call MPI_Barrier(MPI_COMM_WORLD, ierr)


!   if (myrank == 0) then
!      deallocate(global_mat,stat = ierr)
!      deallocate(global_inner,stat = ierr)
!      deallocate(send_buf,stat = ierr)
!   else
!      deallocate(recv_buf,stat = ierr)
!   endif

!   end subroutine scatter_in_3d_logic


end module output_module




