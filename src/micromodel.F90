program micromodel

implicit none

  logical, allocatable:: wall(:,:,:,:),global_inner(:,:,:),io_reorder(:,:,:,:)
  integer :: nkin=1,lx,ly,lz,i,j,k,ii,jj,kk,p,lx_trunc
  integer:: istat

  lx=10
  ly=50
  lz=50


  allocate( wall(1:nkin,1:lx,1:ly,1:lz), STAT=istat)
  if (istat /= 0) then
    write(*,*) 'Allocation of {wall} failed.'
  end if

  allocate( io_reorder(1:lz,1:ly,1:lx,1:nkin), STAT=istat)
  if (istat /= 0) then
    write(*,*) 'Allocation of {io_reorder} failed.'
  end if


  wall = .false.
  wall(1,1,:,:)=.true.
  wall(1,lx,:,:)=.true.


    do i = 1,lx
      do j = 1,ly
        do k = 1,lz
          do p = 1,nkin
            io_reorder(k,j,i,p) = wall(p,i,j,k)
          enddo
        enddo
      enddo
    enddo




    
    open(unit=5, file='./wallin_small.raw', form='unformatted'&
        , access='direct', recl=lx * ly * lz * nkin, action='write', status='replace')

    write(5,rec=1) io_reorder
    close(5)

    ! open(unit=5, file='./micromodel_ft91x122x16.raw', form='unformatted'&
    !     , access='direct', recl=lx_trunc * ly * lz * nkin, action='write', status='replace')

    ! write(5,rec=1) io_reorder(:,:,:lx_trunc,:)
    ! close(5)

end program micromodel
