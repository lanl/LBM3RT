module boundaryg_module

contains
 
  subroutine boundaryg(gi,cj,psi,ux,uy,uz,walls,wall,Rm,bnd,bonds,nbonds,istep)

      use equilg_module
      
      use ptran_global_module
      
      use comvar_module
      
      use con_module

      use mpi

      implicit none
      
      integer :: i,j,k,kk,s,p,nn,istep

      real*8, dimension(1:ncomp,1:mx,1:my+2,1:mz+2,0:6) :: gi,gi_tmp
      real*8, dimension(1:ncomp,0:6) :: gi_pre
      real*8, dimension(1:ncomp,1:mx,1:my+2,1:mz+2) :: cj,psi
      real*8, dimension(1:mx,1:my+2,1:mz+2) :: ux,uy,uz
      real*8, dimension(1:ncomp) :: tmp,psi_tmp
      real*8, dimension(1:ncomp,0:6) :: geq
      logical, dimension(1:mx,0:my+3,0:mz+3) :: walls,bonds
      logical, dimension(nkin,1:mx,0:my+3,0:mz+3) :: wall 
      real*8, dimension(nkin,1:mx,0:my+3,0:mz+3) :: bnd
      real*8, dimension(nkin,1:mx,1:my+2,1:mz+2) :: Rm
      ! real*8, dimension(1:nkin,1:mx,1:my,1:mz+2,1:6) :: bndcomp
      integer, dimension(6,1:mx,0:my+3,0:mz+3) :: nbonds
      real*8 :: multip



     !RHS: const concentration bc

    !  i = mx
    !  do j = 1,my+2
    !    do k = 1,mz+2
    !      !if(.not.walls(i,j,k)) then
    !        do p = 1,m
    !          call equilg(gi(p,i,j,k,:),psiin(p),ux(i,j,k),uy(i,j,k),uz(i,j,k))
    !          psi(p,i,j,k) = psiin(p)
    !        enddo 
    !      !endif 
    !    enddo
    !  enddo

      !RHS: extrapolation bc
    !   i = mx
    !   do j = 1,my+2
    !    do k = 1,mz+2
    !      if(.not.walls(i,j,k)) then
    !        do p = 1,ncomp
    !          psi(p,i,j,k) = psi(p,i-1,j,k)
    !          gi(p,i,j,k,:) = gi(p,i-1,j,k,:)
    !        enddo 
    !      endif 
    !    enddo
    !   enddo

    !  !LHS: const concentration bc
    !  i = 1
    !  do j = 1,my+2
    !    do k = 1,mz+2
    !      !if(.not.walls(i,j,k)) then
    !        do p = 1,ncomp
    !          call equilg(gi(p,i,j,k,:),psiin(p),ux(i,j,k),uy(i,j,k),uz(i,j,k))
    !          psi(p,i,j,k) = psiin(p)
    !        enddo 
    !      !endif 
    !    enddo
    !  enddo

    !  ! BOT: const concentration bc
    !  if(mpicoords(3)==npz-1) then
    !  k = mz+1
    !  do i = 1,mx
    !    do j = 2,my+1
    !      !if(.not.walls(i,j,k)) then
    !        do p = 1,m
    !          call equilg(gi(p,i,j,k,:),psiin(p),ux(i,j,k),uy(i,j,k),uz(i,j,k))
    !          psi(p,i,j,k) = psiin(p)
    !        enddo 
    !      !endif 
    !    enddo
    !  enddo
    !  endif

     ! BOT: extrapolation bc
!     if(mpicoords(3)==npz-1) then
!     k = mz+2
!     do i = 1,mx
!       do j = 2,my+1
!         if(.not.walls(i,j,k)) then
!           do p = 1,ncomp
!             psi(p,i,j,k) = psi(p,i,j,k-1)
!             gi(p,i,j,k,:) = gi(p,i,j,k-1,:)
!           enddo 
!         endif 
!       enddo
!     enddo
!     endif

     !TOP: const concentration bc
!     if(mpicoords(3)==0) then
!     k = 1
!     do i = 1,mx
!       do j = 2,my+1
!         !if(.not.walls(i,j,k)) then
!          do p = 1,ncomp
!             call equilg(gi(p,i,j,k,:),psiin(p),ux(i,j,k),uy(i,j,k),uz(i,j,k))
!             psi(p,i,j,k) = psiin(p)
!          enddo 
!         !endif 
!       enddo
!     enddo
!     endif


    !  !BACK: const concentration bc
    !  if(mpicoords(2)==npy-1) then
    !  j = my+1
    !  do i = 1,mx
    !    do k = 2,mz+1
    !      !if(.not.walls(i,j,k)) then
    !        do p = 1,m
    !          call equilg(gi(p,i,j,k,:),psiin(p),ux(i,j,k),uy(i,j,k),uz(i,j,k))
    !          psi(p,i,j,k) = psiin(p)
    !        enddo 
    !      !endif 
    !    enddo
    !  enddo
    !  endif

    !  !FRONT: const concentration bc
    !  if(mpicoords(2)==0) then
    !  j = 2
    !  do i = 1,mx
    !    do k = 2,mz+1
    !      !if(.not.walls(i,j,k)) then
    !        do p = 1,m
    !          call equilg(gi(p,i,j,k,:),psiin(p),ux(i,j,k),uy(i,j,k),uz(i,j,k))
    !          psi(p,i,j,k) = psiin(p)
    !        enddo 
    !      !endif 
    !    enddo
    !  enddo
    !  endif





      do i=1,mx
        do j=1,my+2
          do k=1,mz+2

            if(walls(i,j,k) .and. (.not.bonds(i,j,k))) then
            ! if(walls(i,j,k)) then
              do s = 1,6
                gi_tmp(:,i,j,k,ind_cj_opp(s)) = gi(:,i,j,k,s)
              enddo

              gi(:,i,j,k,1:6) = gi_tmp(:,i,j,k,1:6)
            endif
          
            !if(.false.) then
            if(bonds(i,j,k)) then
              gi_pre(:,:) = gi(:,i,j,k,:)
              
              if (nbonds(1,i,j,k) == 1) then
                tmp=gi_pre(:,2)
                psi(:,i,j,k) = tmp * 8.d0
                call con_cj_ind(psi(:,i,j,k),cj(:,i,j,k),istep,i,j,k)

                call bnd_cj(tmp,psi(:,i,j,k),cj(:,i,j,k),i,j,k,wall(:,i,j,k),Rm,istep)

                gi(:,i,j,k,1) = psi(:,i,j,k)/8.d0

                do p=1,nkin
                  bnd(p,i,j,k)=bnd(p,i,j,k)+vbarkin(p)*Rm(p,i,j,k)
                enddo 


              end if           

              if (nbonds(2,i,j,k) == 1) then
                tmp=gi_pre(:,1)
                psi(:,i,j,k) = tmp * 8.d0
                call con_cj_ind(psi(:,i,j,k),cj(:,i,j,k),istep,i,j,k)

                call bnd_cj(tmp,psi(:,i,j,k),cj(:,i,j,k),i,j,k,wall(:,i,j,k),Rm,istep)

                gi(:,i,j,k,2) = psi(:,i,j,k)/8.d0!tmp!geq(:,1) + geq(:,2) - gi_pre(:,1)

                do p=1,nkin
                  bnd(p,i,j,k)=bnd(p,i,j,k)+vbarkin(p)*Rm(p,i,j,k)
                enddo 

              end if

              if (nbonds(3,i,j,k) == 1) then
                tmp=gi_pre(:,4)
                psi(:,i,j,k) = tmp * 8.d0
                call con_cj_ind(psi(:,i,j,k),cj(:,i,j,k),istep,i,j,k)
                call bnd_cj(tmp,psi(:,i,j,k),cj(:,i,j,k),i,j,k,wall(:,i,j,k),Rm,istep)

                gi(:,i,j,k,3) = psi(:,i,j,k)/8.d0!tmp!geq(:,3) + geq(:,4) - gi_pre(:,4)

                do p=1,nkin
                  bnd(p,i,j,k)=bnd(p,i,j,k)+vbarkin(p)*Rm(p,i,j,k)
                enddo 

              end if

              if (nbonds(4,i,j,k) == 1) then
                tmp=gi_pre(:,3)
                psi(:,i,j,k) = tmp * 8.d0
                call con_cj_ind(psi(:,i,j,k),cj(:,i,j,k),istep,i,j,k)
                call bnd_cj(tmp,psi(:,i,j,k),cj(:,i,j,k),i,j,k,wall(:,i,j,k),Rm,istep)

                gi(:,i,j,k,4) = psi(:,i,j,k)/8.d0!tmp!geq(:,3) + geq(:,4) - gi_pre(:,3)

                do p=1,nkin
                  bnd(p,i,j,k)=bnd(p,i,j,k)+vbarkin(p)*Rm(p,i,j,k)
                enddo 

              end if

              if (nbonds(5,i,j,k) == 1) then
                tmp=gi_pre(:,6)
                psi(:,i,j,k) = tmp * 8.d0
                call con_cj_ind(psi(:,i,j,k),cj(:,i,j,k),istep,i,j,k)
                call bnd_cj(tmp,psi(:,i,j,k),cj(:,i,j,k),i,j,k,wall(:,i,j,k),Rm,istep)

                gi(:,i,j,k,5) = psi(:,i,j,k)/8.d0!tmp!geq(:,5) + geq(:,6) -gi_pre(:,6)

                do p=1,nkin
                  bnd(p,i,j,k)=bnd(p,i,j,k)+vbarkin(p)*Rm(p,i,j,k)
                enddo 

              end if

              if (nbonds(6,i,j,k) == 1) then
                tmp=gi_pre(:,5)
                psi(:,i,j,k) = tmp * 8.d0
                call con_cj_ind(psi(:,i,j,k),cj(:,i,j,k),istep,i,j,k)
                call bnd_cj(tmp,psi(:,i,j,k),cj(:,i,j,k),i,j,k,wall(:,i,j,k),Rm,istep)

                gi(:,i,j,k,6) = psi(:,i,j,k)/8.d0!tmp!geq(:,5) + geq(:,6) - gi_pre(:,5)

                do p=1,nkin
                  bnd(p,i,j,k)=bnd(p,i,j,k)+vbarkin(p)*Rm(p,i,j,k)
                enddo 

              end if


              
            endif
            do p = 1,ncomp
              if (npos(p) == 1 .and. psi(p,i,j,k) < 0.d0) then
                print *,'LBM/boundaryg: warning neg. psi: '
                ! WRITE(*, '(I0, 1X, A, 1X, E12.4)', ADVANCE='NO') istep,TRIM(nam(p)),psi(p,i,j,k)
                ! WRITE(*, '(I0, 1X,I0, 1X,I0, 1X,I0)', ADVANCE='NO') myrank,i,j,k
                WRITE(*, '(I0, 1X, A, 1X, ES12.4, I0)') istep, TRIM(nam(p)), psi(p, i, j, k)
                WRITE(*, '(A, 1X, I0, 1X, I0, 1X, I0, 1X, I0, 1X, I0, 1X, I0)', ADVANCE='NO') 'rank =',&
                     myrank, mpicoords(2), mpicoords(3), i, j, k
                WRITE(*,*)  

              endif
            enddo
          enddo
        enddo
      enddo
      call MPI_Barrier(MPI_COMM_WORLD, ierr)


    end subroutine boundaryg

    subroutine boundaryg_z_inout(gi,psi_lbm,cj,walls,ux,uy,uz,ncomp,npz,myrank,istep)
      use equilg_module

      use comvar_module

      use mpi

      implicit none
      
      integer :: i,j,k,p,pp,ncomp,npz,myrank,istep
      real*8, dimension(0:6) :: geq
      real*8, dimension(1:ncomp,1:mx,1:my+2,1:mz+2,0:6) :: gi
      real*8, dimension(1:ncomp,1:mx,1:my+2,1:mz+2) :: psi_lbm,cj
      logical, dimension(1:mx,0:my+3,0:mz+3) :: walls
      real*8, dimension(1:mx,1:my+2,1:mz+2) :: ux,uy,uz


     ! BOT: extrapolation bc
     if(mpicoords(3)==npz-1) then
     k = mz+1
     do i = 1,mx
       do j = 1,my+2
         if(.not.walls(i,j,k)) then
           do p = 1,ncomp
             psi_lbm(p,i,j,k) = psi_lbm(p,i,j,k-1)
             cj(p,i,j,k) = cj(p,i,j,k-1)
             do pp = 0,6
               gi(p,i,j,k,pp) = gi(p,i,j,k-1,pp)
             enddo
           enddo
         endif 
       enddo
     enddo
     endif
     call MPI_Barrier(MPI_COMM_WORLD, ierr)

     !TOP: const concentration bc
     if(mpicoords(3)==0) then
     k = 2
     do i = 1,mx
       do j = 1,my+2
         !if(.not.walls(i,j,k)) then
          do p = 1,ncomp
             call equilg(geq(:),psiin(p),ux(i,j,k),uy(i,j,k),uz(i,j,k))
             gi(p,i,j,k,:) = geq(:)
             psi_lbm(p,i,j,k) = psiin(p)
          enddo 
         !endif 
       enddo
     enddo
     endif
     call MPI_Barrier(MPI_COMM_WORLD, ierr)

    end subroutine boundaryg_z_inout
end module boundaryg_module
