module collision_solute_module

  contains

  subroutine collisiong(gi,psi,ux,uy,uz,walls,bonds,ncomp)

    use equilg_module
    
    use comvar_module
    
    implicit none
    
    integer :: i,j,k,s,p,ncomp
    
    real*8, dimension(1:ncomp,mx,my+2,mz+2,0:6) :: gi
    real*8, dimension(1:ncomp,mx,my+2,mz+2) :: psi
    real*8, dimension(mx,my+2,mz+2) :: ux,uy,uz
    logical,dimension(mx,0:my+3,0:mz+3) :: walls,bonds
    real*8,dimension(0:6) :: geq


    do i = 1,mx
      do j = 1,my+2
        do k = 1,mz+2
          if(.not.walls(i,j,k)) then!has to be only fluids without bonds to avoid extra relaxation
          
            do p = 1,ncomp

              call equilg(geq,psi(p,i,j,k),ux(i,j,k),uy(i,j,k),uz(i,j,k))
              gi(p,i,j,k,:)=gi(p,i,j,k,:)-1.d0/tau_p(p)*(gi(p,i,j,k,:)-geq(:))

            enddo
          endif
        enddo
      enddo
    enddo

  end subroutine collisiong
end module collision_solute_module