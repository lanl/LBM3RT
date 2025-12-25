module inout_module

contains

      subroutine input(istep,fi,gi,cj,wall,bnd)
      
      use comvar_module
      use mpi
      
      implicit none
      
      integer :: istep

      real*8, dimension(mx,my+2,mz+2,0:18) :: fi
      real*8, dimension(1:m,1:mx,1:my+2,1:mz+2,0:6) :: gi
      real*8, dimension(1:m,1:mx,1:my+2,1:mz+2) :: cj
      real*8, dimension(l,1:mx,0:my+3,0:mz+3) :: bnd 
      logical, dimension(l,1:mx,0:my+3,0:mz+3) :: wall

      real*8, allocatable :: fi_global,gi_global,cj_global,bnd_global,wall_global

      integer :: rec_length
      
      integer :: iounit, ierr
      integer :: i, j, k, s
      character(len=:), allocatable :: numStr
      character(len=:), allocatable :: filename
      integer, parameter :: MaxIntLength = 30
      allocate(character(len=MaxIntLength) :: numStr)
      allocate(character(len=MaxIntLength) :: filename)

      
      
      ! open(11,file='variables',status='unknown',form='unformatted')
      ! read(11) istep,fi,gi,walls,wall,bnd

      
      WRITE(numStr, '(I8.8)') istep
      numStr = adjustl(trim(numStr))

      rec_length = lx * ly * lz * 19 * 2
      filename = './Results/'//'fi'//numStr//'.raw'
      open(unit=iounit, file=filename, form='unformatted', access='direct', &
             recl=rec_length, status='old', iostat=ierr)
      if (ierr /= 0) then
        print*, 'Error opening file:', filename
        stop
      end if
      read(iounit, rec=1, iostat=ierr) io_reorder
      if (ierr /= 0) then
        print*, 'Error reading file:', filename
        stop
      end if
      close(iounit)
    ! Reorder data
      do i = 1, lx
            do j = 1, ly
            do k = 1, lz
                  do s = 1, 19
                        fi(s-1, i, j, k) = io_reorder(k, j, i, s)
                  end do
            end do
            end do
      end do
      


      rec_length = lx * ly * lz * 7 * m * 2
      filename = './Results/'//'gi'//numStr//'.raw'
      open(unit=iounit, file=filename, form='unformatted', access='direct', &
             recl=rec_length, status='old', iostat=ierr)
      if (ierr /= 0) then
        print*, 'Error opening file:', filename
        stop
      end if
      read(iounit, rec=1, iostat=ierr) gi
      if (ierr /= 0) then
        print*, 'Error reading file:', filename
        stop
      end if
      close(iounit)
      


      rec_length = lx * ly * lz * m * 2
      filename = './Results/'//'cj'//numStr//'.raw'
      open(unit=iounit, file=filename, form='unformatted', access='direct', &
             recl=rec_length, status='old', iostat=ierr)
      if (ierr /= 0) then
        print*, 'Error opening file:', filename
        stop
      end if
      read(iounit, rec=1, iostat=ierr) cj
      if (ierr /= 0) then
        print*, 'Error reading file:', filename
        stop
      end if
      close(iounit)


      rec_length = lx * ly * lz
      filename = './Results/'//'walls'//numStr//'.raw'
      open(unit=iounit, file=filename, form='unformatted', access='direct', &
             recl=rec_length, status='old', iostat=ierr)
      if (ierr /= 0) then
        print*, 'Error opening file:', filename
        stop
      end if
      read(iounit, rec=1, iostat=ierr) walls
      if (ierr /= 0) then
        print*, 'Error reading file:', filename
        stop
      end if
      close(iounit)


      rec_length = lx * ly * lz * l * 2
      filename = './Results/'//'bnd'//numStr//'.raw'
      open(unit=iounit, file=filename, form='unformatted', access='direct', &
             recl=rec_length, status='old', iostat=ierr)
      if (ierr /= 0) then
        print*, 'Error opening file:', filename
        stop
      end if
      read(iounit, rec=1, iostat=ierr) bnd
      if (ierr /= 0) then
        print*, 'Error reading file:', filename
        stop
      end if
      close(iounit)

      rec_length = lx * ly * lz * l
      filename = './Results/'//'wall'//numStr//'.raw'
      open(unit=iounit, file=filename, form='unformatted', access='direct', &
             recl=rec_length, status='old', iostat=ierr)
      if (ierr /= 0) then
        print*, 'Error opening file:', filename
        stop
      end if
      read(iounit, rec=1, iostat=ierr) wall
      if (ierr /= 0) then
        print*, 'Error reading file:', filename
        stop
      end if
      close(iounit)


      return
      end subroutine input
      
!-----------------------------------------------------------------

      subroutine output(istep,fi,gi,walls,wall,bnd)
      
      use comvar_module
      
      implicit none
      
      integer :: istep

      real*8, dimension(0:8,lx,ly) :: fi
      real*8, dimension(1:m,0:8,lx,ly) ::gi
      real*8, dimension(l,lx,ly) :: bnd
      logical, dimension(lx,ly) :: walls
      logical, dimension(l,lx,ly) :: wall
      
      open(11,file='variables',status='unknown',form='unformatted')
      write(11) istep,fi,gi,walls,wall,bnd

!     open(11,file='variables',status='unknown')
!     write(11,*) istep,fi,gi,walls,wall,bnd
      
      close(11)
      return
      end subroutine output
      
      subroutine convert(t_phys,u_phys,k_phys,l_phys,nt_lb)
      
      use comvar_module
      use ptran_global_module
      
      implicit none
      
      integer :: k,nt_lb
      real*8 :: t_phys,u_phys,k_phys(*),l_phys
      
      dif_lb = (2.d0*tau0-1.d0)/6.d0
      
      pe = vxt * lx / dif_lb
      
      u_phys = difaq * pe / l_phys
      
      dtlb = dif_lb/difaq*(l_phys/lx)**2
        
      t_phys = dtlb*nt_lb
      
      do k = 1, nkin
        k_phys(k) = l_phys/lx*difaq/dif_lb*kb(k)
      enddo
      
      return
      
      end subroutine convert
      
end module inout_module
