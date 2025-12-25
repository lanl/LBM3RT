module denvel_module

contains

!------------------------- calculate density and local velocities-----
  subroutine rho_u(fi,rho,ux,uy,uz,walls)
    
    use comvar_module
    
    implicit none


    integer :: s
    real*8, dimension(mx,my+2,mz+2,0:18) :: fi      
    real*8, dimension(mx,my+2,mz+2) :: rho,ux,uy,uz
    logical, dimension(mx,my+2,mz+2) ::walls
    

    ux=0.d0
    uy=0.d0
    uz=0.d0

    where(.not.walls) rho=fi(:,:,:,0)
        
    do s=1,18
      where(.not.walls)
        rho=rho+fi(:,:,:,s)
        ux=ux+cix(s)*fi(:,:,:,s)
        uy=uy+ciy(s)*fi(:,:,:,s)
        uz=uz+ciz(s)*fi(:,:,:,s)
      end where
    enddo
        
    where(walls)
      ux=0.d0
      uy=0.d0
      uz=0.d0
    end where
  end subroutine rho_u

  subroutine rho_u_compress(fi,rho,ux,uy,uz,walls)

    use comvar_module

    implicit none


    integer :: s
    real*8, dimension(mx,my+2,mz+2,0:18) :: fi
    real*8, dimension(mx,my+2,mz+2) :: rho,ux,uy,uz
    logical, dimension(mx,my+2,mz+2) ::walls


    ux=0.d0
    uy=0.d0
    uz=0.d0

    where(.not.walls) rho=fi(:,:,:,0)

    do s=1,18
      where(.not.walls)
        rho=rho+fi(:,:,:,s)
        ux=ux+cix(s)*fi(:,:,:,s)
        uy=uy+ciy(s)*fi(:,:,:,s)
        uz=uz+ciz(s)*fi(:,:,:,s)
      end where
    enddo

    where(walls)
      ux=0.d0
      uy=0.d0
      uz=0.d0
    elsewhere
      ux=(ux + 0.5d0*ffx)/rho
      uy=uy/rho
      uz=uz/rho
    end where
  end subroutine rho_u_compress
end module denvel_module

