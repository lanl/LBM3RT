program lbm3d_mpi

  use mpi

  use comvar_module

  use ptran_global_module
  use trdynmem_module
  use ptran_read_module
  use ptran_init_module

  use setupMPI_module
  use mpi_comm_module

  use initialize_flow_module
  use equilf_module
  use denvel_module
  use streaming_flow_module
  use collision_flow_module
  use boundaryf_module

  use bindex_module
  use initialize_solute_module
  use equilg_module
  use output_module
  use streaming_solute_module
  use collision_solute_module
  use boundaryg_module
  use con_module

  use expand_module
  use redistri_bnd_module


  implicit none


  integer :: i,j,k,s,p,ii,il,pp, &
              istep,itermx,plt_ind0,plt_ind1,plt_ind2,plt_ind3,sum_walls,glb_sum_walls
  real*8, allocatable, dimension(:,:,:,:)   :: fi
  real*8, allocatable, dimension(:,:,:)     :: rho_lbm,ux,uy,uz
  logical, allocatable, dimension(:,:,:)   :: walls
  real*8, allocatable, dimension(:,:,:,:,:) :: gi
  real*8, allocatable, dimension(:,:,:,:)   :: psi_lbm,cj,ci,bnd
  logical, allocatable, dimension(:,:,:,:) :: wall
  logical, allocatable, dimension(:,:,:)   :: bonds
  integer, allocatable, dimension(:,:,:,:)   :: nbonds
  integer, allocatable,dimension(:,:,:) :: totbonds
  real*8, allocatable, dimension(:,:,:,:) :: Rm
  character(len=:), allocatable :: baseFilename  
  character(len=:), allocatable :: save_dir
  character(len=255) :: path = './Results/'
  real*8 :: start_time, end_time
  integer :: minutes

  


 
  ! Initialize mpi
  call mpi_init(ierr)
  call MPI_comm_size(MPI_COMM_WORLD, commsize, ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)


  call ptran_read
  call ptran_chem

  call comvar()

  allocate( tau_p(1:ncomp) )
  allocate( cjini(1:ncomp) )
  allocate( psiini(1:ncomp) )
  allocate( cjin(1:ncomp) )
  allocate( psiin(1:ncomp) )
  allocate( c1(1:ncomp) )
  allocate( gama_j(1:ncomp) )
  allocate( ciin(1:ncmplx) )
  allocate( ciini(1:ncmplx) )
  allocate( gama_i(1:ncmplx) )
  allocate( kb(1:nkin) )
  allocate( npos(1:ncomp) )
  allocate( alog(1:ncomp) )
  allocate(ny_local(1:npy))
  allocate(nz_local(1:npz))
  allocate(my_local(1:npy))
  allocate(mz_local(1:npz))

  lx = nx
  ly = ny
  lz = nz  

  call setupMPI(nkin,ncomp,commsize,myrank,npy,npz)

  tau=1.0
  !print*,myrank,tau0,'myrank&tau',dx0(1),lx

  !set tau_p
  do i=1,ncomp
    tau_p(i) = tau0
  enddo

  ! store initial and boundary concentrations
  do j = 1, ncomp
    cjini(j) = c0(j)
    cjin(j) = ccbnd(j,1)
!       write(*,*) nam(j),c0(j),ccbnd(j,1)
  enddo



  ! store mineral logK and stoichiometry
!     l_phys = (nx-2)*dx0(1) ! top and bottom nodes inert boundaries
  l_phys = dble(lx)*dx0(1) ! top and bottom nodes inert boundaries

  dif_lb = (2.d0*tau_p(1)-1.d0)/8.d0

  do i = 1, nkin
    vbarkin(i) = vbarkin(i)*1.0d0
    !convert rate constants to lattice units(M*lt/dt)
    !klb = rkf(i)*(l_phys/dble(lx))*(dif_lb/difaq)
    ! klb = rkf(i)
    !kb(i) = klb
    kb(i)=rkf(i)*dy0(1)*(dif_lb/difaq)

    
    ! Damkoehler number
    !damk = kb(i)*dble(lz)/dif_lb*10.0d0**eqkin(i)   
    damk = rkf(i)*vbarkin(i)*l_phys/difaq
!       damk1 = rkf(i)*l_phys/difaq
  enddo


  if(myrank==0) then !print physical properties
    write(*,'(/,'' Lattice Boltzmann Parameters'')')
    write(*, '(A, I3)') ' Total number of minerals: ', nkin
    write(*, '(A, I3)') ' Total number of primary chemicals: ', ncomp
    write(*, '(A, I3)') ' Total number of secondary chemicals: ', ncmplx
    write(*,'('' dif_lb = '',1pe12.4,'' [l.u.] difaq = '',1pe12.4, &
  &   '' [m^2/s]'','' l = '',1pe12.4,'' [m]'')') &
    dif_lb,difaq,l_phys
    write(*,'('' mineral'',8x,''  log K  '','' rk[mol/dm^3 m/s]'', &
  &   ''   klb      '',''    Da    '',''molar vol.[dm^3/mol]'')')

    do i = 1, nkin
      write(*,'(1x,a12,1p10e12.4)') &
      namk(i),eqkin(i),rkf(i),kb(i),damk,vbarkin(i)
      
      ! set primary mineral: note this only works for one primary phase
      !if (phik_reg(i,1) == 1.d0) iprimary_min = i
    enddo
  endif


  !Memory allocation
  allocate(character(len=len_trim(path)) :: save_dir)
  save_dir = path

  allocate( fi(1:mx,1:my+2,1:mz+2,0:18), STAT=istat)
  if (istat /= 0) then
    write(*,*) 'Allocation of {fi} failed.'
  end if
  allocate( ux(1:mx,1:my+2,1:mz+2), STAT=istat)
  if (istat /= 0) then
    write(*,*) 'Allocation of {ux} failed.'
  end if
  allocate( uy(1:mx,1:my+2,1:mz+2), STAT=istat)
  if (istat /= 0) then
    write(*,*) 'Allocation of {uy} failed.'
  end if
  allocate( uz(1:mx,1:my+2,1:mz+2), STAT=istat)
  if (istat /= 0) then
    write(*,*) 'Allocation of {uz} failed.'
  end if
  allocate( rho_lbm(1:mx,1:my+2,1:mz+2), STAT=istat)
  if (istat /= 0) then
    write(*,*) 'Allocation of {rho_lbm} failed.'
  end if
  allocate( gi(1:ncomp,1:mx,1:my+2,1:mz+2,0:6), STAT=istat)
  if (istat /= 0) then
    write(*,*) 'Allocation of {gi} failed.'
  end if
  allocate( psi_lbm(1:ncomp,1:mx,1:my+2,1:mz+2), STAT=istat )
  if (istat /= 0) then
    write(*,*) 'Allocation of {psi_lbm} failed.'
  end if
  allocate( cj(1:ncomp,1:mx,1:my+2,1:mz+2), STAT=istat )
  if (istat /= 0) then
    write(*,*) 'Allocation of {cj} failed.'
  end if
  allocate( ci(1:ncmplx,1:mx,1:my+2,1:mz+2), STAT=istat )
  if (istat /= 0) then
    write(*,*) 'Allocation of {ci} failed.'
  end if

  allocate( bnd(1:nkin,1:mx,0:my+3,0:mz+3) , STAT=istat)
  if (istat /= 0) then
    write(*,*) 'Allocation of {bnd} failed.'
  end if  
  allocate( Rm(1:nkin,1:mx,1:my+2,1:mz+2) , STAT=istat)
  if (istat /= 0) then
    write(*,*) 'Allocation of {Rm} failed.'
  end if  
  allocate( bonds(1:mx,0:my+3,0:mz+3), STAT=istat )
  if (istat /= 0) then
    write(*,*) 'Allocation of {bonds} failed.'
  end if
  allocate( nbonds(1:6,1:mx,0:my+3,0:mz+3), STAT=istat )
  if (istat /= 0) then
    write(*,*) 'Allocation of {nbonds} failed.'
  end if
  allocate( totbonds(1:mx,0:my+3,0:mz+3), STAT=istat )
  if (istat /= 0) then
    write(*,*) 'Allocation of {totbonds} failed.'
  end if
  allocate( wall(1:nkin,1:mx,0:my+3,0:mz+3), STAT=istat)
  if (istat /= 0) then
    write(*,*) 'Allocation of {walls} failed.'
  end if


  allocate( walls(1:mx,0:my+3,0:mz+3), STAT=istat)
  if (istat /= 0) then
    write(*,*) 'Allocation of {walls} failed.'
  end if

  ! print*,'memory allocated on rank = ', myrank  
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if(myrank==0) write(*,*) 'Memory allocated for LBM.'
  ! write(*,*) 'memory allocated on rank = ', myrank
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  call random_seed()


  plt_ind0=0
  plt_ind1=0
  plt_ind2=0
  plt_ind3=0

 ! set positivity array for total concentration psi
  do j = 1, ncomp
    npos(j) = 1
    do i = 1, ncmplx
      if (shom(j,i) < 0.d0) then
        npos(j) = -1
      endif
    enddo
  enddo

  !initialize activity coefficients
  do i = 1, ncomp
    gama_j(i) = 1.d0
  enddo
  do i = 1, ncmplx
    gama_i(i) = 1.d0
  enddo

!     compute initial total concentration for secondary species     
  if (ncmplx > 0) then
    do j = 1,ncomp
      alog(j) = log(gama_j(j)*cjin(j))
    enddo
    do i=1,ncmplx
      ciin(i) = aln10*eqhom(i)-log(gama_i(i))
      do j=1,ncomp
        ciin(i)=ciin(i)+shom(j,i)*alog(j)
      enddo
      ciin(i) = exp(ciin(i))
    enddo
  endif

  if (ncmplx > 0) then
    do j = 1,ncomp
      alog(j) = log(gama_j(j)*cjini(j))
    enddo
    do i=1,ncmplx
      ciini(i) = aln10*eqhom(i)-log(gama_i(i))
      do j=1,ncomp
        ciini(i)=ciini(i)+shom(j,i)*alog(j)
      enddo
      ciini(i) = exp(ciini(i))
    enddo
  endif

!     compute initial total concentration

  do j=1,ncomp
    psiin(j)=cjin(j)
    do i=1,ncmplx
      psiin(j)=psiin(j)+shom(j,i)*ciin(i)
    enddo
  enddo		

  do j=1,ncomp
    psiini(j)=cjini(j)
    do i=1,ncmplx
      psiini(j)=psiini(j)+shom(j,i)*ciini(i)
    enddo
  enddo		

  if(myrank==0) then
    print*,'ciin'
    do j=1,ncmplx
      WRITE(*, '(ES12.5)') ciin(j)
    enddo
  endif

  if(myrank==0) then
    print*,'cjin'
    do j=1,ncomp
      WRITE(*, '(ES12.5)') cjin(j)
    enddo
  endif

  if(myrank==0) then
    print*,'psiin'
    do j=1,ncomp
      WRITE(*, '(ES12.5)') psiin(j)
    enddo
  endif


  if(new) then
    istep=0
    if(myrank==0) then
      write(*,*) 'initialization start'  
    endif

    call initialize_flow(fi,rho_lbm,ux,uy,uz,walls,wall,nkin,rhol,rhor,myrank,npy,npz)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    baseFilename = 'wallin_small'
 
    call scatter_in_4d_logic2(wall,save_dir,baseFilename,nkin,1,mx,0,my+3,0,mz+3,2,2,commsize,npy,npz,myrank)

    
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call mpi_3d(ux)
    call mpi_3d(uy)
    call mpi_3d(uz)
    call mpi_3d(rho_lbm)
    call mpi_wall(wall,nkin)  
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    walls(:,:,:)=wall(1,:,:,:)	
    do il=2,nkin
      walls(:,:,:)=walls(:,:,:).or.wall(il,:,:,:) 
    enddo
    do i=1,mx
      do j=1,my+2
        do k=1,mz+2
          call equilf(fi(i,j,k,:),rho_lbm(i,j,k),ux(i,j,k),uy(i,j,k),uz(i,j,k))
        enddo
      enddo
    enddo


    call bindex(walls,bonds,nbonds,totbonds) 

    call initialize_solute(cj,walls,ncomp) 


    call MPI_Barrier(MPI_COMM_WORLD, ierr) 
    call mpi_cj(cj,ncomp)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    call con_ci(ci,cj)


    psi_lbm=cj
    do j=1,ncomp
      do i=1,ncmplx
        psi_lbm(j,:,:,:)=psi_lbm(j,:,:,:)+shom(j,i)*ci(i,:,:,:)
      enddo
    enddo
    

    bnd(:,:,:,:)=bsat

    do i=1,mx
      do j=1,my+2
        do k = 1,mz+2
          do p=1,ncomp
            call equilg(gi(p,i,j,k,:),psi_lbm(p,i,j,k),ux(i,j,k),uy(i,j,k),uz(i,j,k))
          enddo
        enddo
      enddo
    enddo

    istep=0
    if ((mod(istep, iwrite) .eq. 0)) then    
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      baseFilename = 'fi'
      call gather_out_fi(fi,save_dir,baseFilename,istep,1,mx,1,my+2,1,mz+2,0,18,1,1,commsize,npy,npz,myrank)      
    endif

    if ((mod(istep, iwrite) .eq. 0)) then    
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      baseFilename = 'gi'
      call gather_out_gi(gi,save_dir,baseFilename,istep,ncomp,1,mx,1,my+2,1,mz+2,0,6,1,1,commsize,npy,npz,myrank)      
    endif

    if ((mod(istep, iwrite) .eq. 0)) then
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      baseFilename = 'cj'
      call gather_out_4d_dble(cj,save_dir,baseFilename,istep,ncomp,1,mx,1,my+2,1,mz+2,1,1,commsize,npy,npz,myrank) 
    endif

    ! if ((mod(istep, iwrite) .eq. 0)) then
    !   call con_ci(ci,cj)
    !   call MPI_Barrier(MPI_COMM_WORLD, ierr)
    !   baseFilename = 'ci'
    !   call gather_out_4d_dble(ci,save_dir,baseFilename,istep,ncmplx,1,mx,1,my+2,1,mz+2,1,1,commsize,npy,npz,myrank)
    ! endif
    
    if ((mod(istep, iwrite) .eq. 0)) then
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      baseFilename = 'bnd'
      call gather_out_4d_dble(bnd,save_dir,baseFilename,istep,nkin,1,mx,0,my+3,0,mz+3,2,2,commsize,npy,npz,myrank) 
    endif

    if ((mod(istep, iwrite) .eq. 0)) then
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      baseFilename = 'wall'
      call gather_out_4d_logic(wall,save_dir,baseFilename,istep,nkin,1,mx,0,my+3,0,mz+3,2,2,commsize,npy,npz,myrank) 
    endif

  elseif(.not.new) then !restart pre-processing

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    baseFilename = 'bnd'
    call scatter_in_4d_dble(bnd,save_dir,baseFilename,iread,nkin,1,mx,0,my+3,0,mz+3,2,2,commsize,npy,npz,myrank)  

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    baseFilename = 'cj'
    call scatter_in_4d_dble(cj,save_dir,baseFilename,iread,ncomp,1,mx,1,my+2,1,mz+2,1,1,commsize,npy,npz,myrank)


    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    baseFilename = 'wall'
    call scatter_in_4d_logic(wall,save_dir,baseFilename,iread,nkin,1,mx,0,my+3,0,mz+3,2,2,commsize,npy,npz,myrank)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    baseFilename = 'fi'
    call scatter_in_fi(fi,save_dir,baseFilename,iread,1,mx,1,my+2,1,mz+2,0,18,1,1,commsize,npy,npz,myrank)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    baseFilename = 'gi'
    call scatter_in_gi(gi,save_dir,baseFilename,iread,ncomp,1,mx,1,my+2,1,mz+2,0,6,1,1,commsize,npy,npz,myrank)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
  
    istep = iread


    call mpi_wall(wall,nkin) 
    call mpi_fi_full(fi)
    call mpi_gi_full(gi,ncomp) 
    call mpi_bnd(bnd,nkin)
    call mpi_cj(cj,ncomp)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    walls(:,:,:)=wall(1,:,:,:)	
    do il=2,nkin	
      walls(:,:,:)=walls(:,:,:).or.wall(il,:,:,:) 
    enddo

    call rho_u(fi,rho_lbm,ux,uy,uz,walls(:,1:my+2,1:mz+2)) 

    call bindex(walls,bonds,nbonds,totbonds) 

    call con_ci(ci,cj)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    do i=1,mx
      do j=1,my+2
        do k=1,mz+2
          do p=1,ncomp
            psi_lbm(p,i,j,k)=sum(gi(p,i,j,k,:))
          enddo
        enddo
      enddo
    enddo

  endif




!-----begin time stepping-----

  if(myrank==0) then
    write(*,*) 'Loop start'  
    call CPU_TIME(start_time)    
  endif


  do ii = 1,ntimes
    istep=istep+1

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call mpi_fi(fi)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call streamingf(fi,walls(:,1:my+2,1:mz+2),myrank)

    call boundaryf_z_incompressible(fi,rho_lbm,ux,uy,uz,rhol,rhor,uin,walls(:,1:my+2,1:mz+2),npz)
    call rho_u(fi,rho_lbm,ux,uy,uz,walls(:,1:my+2,1:mz+2)) 
    call collisionf(fi,rho_lbm,ux,uy,uz,walls(:,1:my+2,1:mz+2))

!    call MPI_Barrier(MPI_COMM_WORLD, ierr)
!    call mpi_fi(fi)
!    call MPI_Barrier(MPI_COMM_WORLD, ierr)

!    call streamingf(fi,walls(:,1:my+2,1:mz+2),myrank)

!    call boundaryf_z(fi,rho_lbm,ux,uy,uz,rhol,rhor,uin,walls(:,1:my+2,1:mz+2),npz)
!    call rho_u_compress(fi,rho_lbm,ux,uy,uz,walls(:,1:my+2,1:mz+2))
!    call collisionf(fi,rho_lbm,ux,uy,uz,walls(:,1:my+2,1:mz+2))

    if(istep>10000) then
    call streamingg(gi,ncomp)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call mpi_gi_full(gi,ncomp)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)    
    call boundaryg_z_inout(gi,psi_lbm,cj,walls,ux,uy,uz,ncomp,npz,myrank,istep)

    do i=1,mx
      do j=1,my+2
        do k=1,mz+2
          if(.not.walls(i,j,k)) then !has to be fluids
            do p=1,ncomp
              psi_lbm(p,i,j,k)=sum(gi(p,i,j,k,:))
            enddo
          endif
        enddo
      enddo
    enddo

    call collisiong(gi,psi_lbm,ux,uy,uz,walls,bonds,ncomp)

    call con_cj(psi_lbm,cj,walls,bonds,ii,itermx,&
                istep)

    call boundaryg(gi,cj,psi_lbm,ux,uy,uz,walls,wall,Rm,&
              bnd,bonds,nbonds,istep)!updated 11.02.2024 for logics

  

    !Random

    dpst=.false.

    do i=1,mx
      do j=1,my+2
        do k=1,mz+2
          do p=1,nkin
            if(bonds(i,j,k).and.bnd(p,i,j,k).ge.nbsat*bsat) then
              call expand(walls,bnd(p,:,:,:),i,j,k,totbonds)

            else if(bonds(i,j,k).and.bnd(p,i,j,k).le.0.d0) then
              call redistri_bnd(bnd,walls,totbonds,i,j,k,bsat,p,nkin)      

            endif 
          enddo
        enddo
      enddo   
    enddo 

    do i=1,mx
      do j=2,my+1
        do k=2,mz+1
          do p=1,nkin
            if((.not.walls(i,j,k)).and.bnd(p,i,j,k).ge.nbsat*bsat) then
              wall(p,i,j,k) = .true.
              dpst=.true.
              bnd(p,i,j,k) = bnd(p,i,j,k) - bsat
            endif
 
            if(bonds(i,j,k).and.bnd(p,i,j,k).le.0.d0) then
              wall(p,i,j,k)=.false.
              dpst=.true.
              bnd(p,i,j,k)=bnd(p,i,j,k) + bsat
              
              call extra_relaxation(gi,totbonds,nbonds,i,j,k,ncomp)!updated 11.02.2024 for logics


            endif
          enddo
        enddo
      enddo   
    enddo 



    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call mpi_wall(wall,nkin)  
    call mpi_bnd(bnd,nkin)   
    call mpi_dpst()
    call mpi_gi_full(gi,ncomp)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)


    if(dpst.or.dpst_up.or.dpst_down.or.dpst_front.or.dpst_rear.or. &
        dpst_front_down.or.dpst_front_up.or.dpst_rear_down.or.dpst_rear_up) then
      walls(:,:,:)=wall(1,:,:,:)	
      do il=2,nkin	
        walls(:,:,:)=walls(:,:,:).or.wall(il,:,:,:) 
      enddo
      call bindex(walls,bonds,nbonds,totbonds)  

    endif  

    endif

    if ((mod(istep, iwrite) .eq. 0)) then
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      baseFilename = 'uz'
      call gather_out_3d_dble(uz,save_dir,baseFilename,istep,1,mx,1,my+2,1,mz+2,1,1,commsize,npy,npz,myrank) 
    endif

    if ((mod(istep, iwrite) .eq. 0)) then
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      baseFilename = 'ux'
      call gather_out_3d_dble(ux,save_dir,baseFilename,istep,1,mx,1,my+2,1,mz+2,1,1,commsize,npy,npz,myrank)
    endif

    if ((mod(istep, iwrite) .eq. 0)) then
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      baseFilename = 'uy'
      call gather_out_3d_dble(uy,save_dir,baseFilename,istep,1,mx,1,my+2,1,mz+2,1,1,commsize,npy,npz,myrank)
    endif

!    if ((mod(istep, iwrite) .eq. 0)) then
!      call MPI_Barrier(MPI_COMM_WORLD, ierr)
!      baseFilename = 'rho'
!      call gather_out_3d_dble(rho_lbm,save_dir,baseFilename,istep,1,mx,1,my+2,1,mz+2,1,1,commsize,npy,npz,myrank)
!    endif

    if ((mod(istep, iwrite) .eq. 0)) then
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      baseFilename = 'psi'
      call gather_out_4d_dble(psi_lbm,save_dir,baseFilename,istep,ncomp,1,mx,1,my+2,1,mz+2,1,1,commsize,npy,npz,myrank) 
    endif

    if ((mod(istep, iwrite) .eq. 0)) then    
     call MPI_Barrier(MPI_COMM_WORLD, ierr)
     baseFilename = 'fi'
     call gather_out_fi(fi,save_dir,baseFilename,istep,1,mx,1,my+2,1,mz+2,0,18,1,1,commsize,npy,npz,myrank)      
    endif

    if ((mod(istep, iwrite) .eq. 0).and.istep>=10000) then    
     call MPI_Barrier(MPI_COMM_WORLD, ierr)
     baseFilename = 'gi'
     call gather_out_gi(gi,save_dir,baseFilename,istep,ncomp,1,mx,1,my+2,1,mz+2,0,6,1,1,commsize,npy,npz,myrank)      
    endif

    if ((mod(istep, iwrite) .eq. 0)) then
     call MPI_Barrier(MPI_COMM_WORLD, ierr)
     baseFilename = 'cj'
     call gather_out_4d_dble(cj,save_dir,baseFilename,istep,ncomp,1,mx,1,my+2,1,mz+2,1,1,commsize,npy,npz,myrank) 
    endif

    !if ((mod(istep, iwrite) .eq. 0)) then
    ! call con_ci(ci,cj)
    ! call MPI_Barrier(MPI_COMM_WORLD, ierr)
    ! baseFilename = 'ci'
    ! call gather_out_4d_dble(ci,save_dir,baseFilename,istep,ncmplx,1,mx,1,my+2,1,mz+2,1,1,commsize,npy,npz,myrank)
    !endif
    
    if ((mod(istep, iwrite) .eq. 0)) then
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      baseFilename = 'bnd'
      call gather_out_4d_dble(bnd,save_dir,baseFilename,istep,nkin,1,mx,0,my+3,0,mz+3,2,2,commsize,npy,npz,myrank) 
    endif

    if ((mod(istep, iwrite) .eq. 0)) then
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      baseFilename = 'wall'
      call gather_out_4d_logic(wall,save_dir,baseFilename,istep,nkin,1,mx,0,my+3,0,mz+3,2,2,commsize,npy,npz,myrank) 
    endif

    ! if ((mod(istep, iwrite) .eq. 0)) then
    !   call MPI_Barrier(MPI_COMM_WORLD, ierr)
    !   baseFilename = 'bnd'
    !   call gather_out_4d_dble(bnd,baseFilename,istep,l,1,mx,0,my+3,0,mz+3,2,2,commsize,npy,npz,myrank) 
    ! endif
!     ! if ((mod(istep, iwrite) .eq. 0)) then
!     !   baseFilename = 'fi'
!     !   call oflow_4d_f(mx,my,mz+2,istep,fi(:,:,:,:),19,baseFilename,myrank,commsize)
!     ! endif 

!     ! if ((mod(istep, iwrite) .eq. 0)) then
!     !   baseFilename = 'ux'
!     !   call oflow_3d_d(mx,my+2,mz+2,istep,ux(:,:,:),baseFilename,commsize,mpicoords(2),mpicoords(3))
!     ! endif  

!     ! if ((mod(istep, iwrite) .eq. 0)) then
!     !   baseFilename = 'rho'
!     !   call oflow_3d_d(mx,my+2,mz+2,istep,rho_lbm(:,:,:),baseFilename,commsize,mpicoords(2),mpicoords(3))
!     ! endif  

!     ! if ((mod(istep, iwrite) .eq. 0)) then
!     !   baseFilename = 'walls'
!     !   call oflow_3d_bo_2(mx,my+4,mz+4,istep,walls,baseFilename,commsize,mpicoords(2),mpicoords(3))
!     ! endif  

!     ! if ((mod(istep, iwrite) .eq. 0)) then
!     !   baseFilename = 'bonds'
!     !   call oflow_3d_bo_2(mx,my+4,mz+4,istep,bonds,baseFilename,commsize,mpicoords(2),mpicoords(3))
!     ! endif  

!     ! if ((mod(istep, iwrite) .eq. 0)) then
!     !   baseFilename = 'gi'
!     !   call oflow_5d_g(mx,my+2,mz+2,istep,gi,m,7,baseFilename,commsize,mpicoords(2),mpicoords(3))
!     ! endif  

    ! call MPI_Barrier(MPI_COMM_WORLD, ierr)
    ! if ((mod(istep, iwrite) .eq. 0).and.istep==349) then
    !   baseFilename = 'wall0'
    !   call oflow_4d_bo_2(mx,my+4,mz+4,istep,wall,l,baseFilename,commsize,mpicoords(2),mpicoords(3))
    ! endif  
    ! call MPI_Barrier(MPI_COMM_WORLD, ierr)

!     ! if ((mod(istep, iwrite) .eq. 0)) then
!     !   baseFilename = 'wall_buff'
!     !   call oflow_4d_bo(mx,my,4,istep,wall_buff,l,baseFilename,myrank,commsize)
!     ! endif 

    ! call MPI_Barrier(MPI_COMM_WORLD, ierr)
    ! if ((mod(istep, iwrite) .eq. 0)) then
    !   baseFilename = 'bnd'
    !   call oflow_4d_d_2(mx,my+4,mz+4,istep,bnd,l,baseFilename,commsize,mpicoords(2),mpicoords(3))
    ! endif 
    ! call MPI_Barrier(MPI_COMM_WORLD, ierr)
!     ! if ((mod(istep, iwrite) .eq. 0)) then
!     !   baseFilename = 'bnd_buff'
!     !   call oflow_4d_d(mx,my,4,istep,bnd_buff,l,baseFilename,myrank,commsize)
!     ! endif  

!     ! if ((mod(istep, iwrite) .eq. 0)) then
!     !   baseFilename = 'totbonds'
!     !   call oflow_3d_i(mx,my,mz+2,istep,totbonds,baseFilename,myrank,commsize)
!     ! endif  

!     ! if ((mod(istep, iwrite) .eq. 0)) then
!     !   baseFilename = 'totbonds_buff'
!     !   call oflow_3d_i(mx,my,4,istep,totbonds_buff,baseFilename,myrank,commsize)
!     ! endif  

!   !   sum_walls = 0
!   !   do i=1,mx
!   !     do j=1,my
!   !       do k=2,mz+1
!   !         if(walls(i,j,k)) then
!   !           sum_walls = sum_walls+1
!   !         endif
!   !       enddo
!   !     enddo
!   !   enddo

!   ! call MPI_Allreduce(sum_walls, glb_sum_walls, 1, MPI_INTEGER, MPI_SUM, lbecomm, ierr)

!   ! if(glb_sum_walls>=29990.and.glb_sum_walls<=30010.and.(plt_ind0.eq.0)) then
!   !   baseFilename = 'walls30000_'
!   !   call oflow_3d_bo_2(mx,my,mz+4,istep,walls,baseFilename,myrank,commsize)
!   !   plt_ind0 = plt_ind0 + 1
!   ! endif

!   ! if(glb_sum_walls>=59990.and.glb_sum_walls<=60010.and.(plt_ind1.eq.0)) then
!   !   baseFilename = 'walls60000_'
!   !   call oflow_3d_bo_2(mx,my,mz+4,istep,walls,baseFilename,myrank,commsize)
!   !   plt_ind1 = plt_ind1 + 1
!   ! endif

!  if(myrank==0.and.istep>=1536589) then
!    print*,istep, 'step is finished!' 
!  endif


  if(myrank==0.and.(mod(istep, icheck) .eq. 0).and.istep>10000) then 
    call CPU_TIME(end_time)  
    minutes = FLOOR((end_time - start_time)/60.0d0)
    WRITE(*, '(A, 1X, I0, A, 1X, I0, 1X, A, 1X, F0.1, 1X, A, 1X)') 'commsize =',commsize,&
      ', elapsed time =',minutes,'mins',end_time - start_time - 60.0d0*dble(minutes),'s.'
  endif


  enddo




  ! Finalize MPI
  ! call MPI_Type_free(xyplane_db1, ierr)
  ! call MPI_Type_free(xzplane_db1, ierr)
  ! call MPI_Type_free(xyplane_db2, ierr)
  ! call MPI_Type_free(xzplane_db2, ierr)
  ! call MPI_Type_free(xyplane_db3, ierr)
  ! call MPI_Type_free(xzplane_db3, ierr)
  ! call MPI_Type_free(xyplane_db5, ierr)
  ! call MPI_Type_free(xzplane_db5, ierr)
  ! call MPI_Type_free(xyplane_logical, ierr)
  ! call MPI_Type_free(xzplane_logical, ierr)

  call mpi_finalize(ierr)

  
end program lbm3d_mpi
