module expand_module

contains

      subroutine expand(wall_ind_ex,bnd_ind_ex,i,j,k,totbonds)

      use comvar_module
      use perind_module

      implicit none

      integer :: icounter,i,j,k,i_p,j_p,k_p,i_n,j_n,k_n
      real*8 :: dchance

      real*8, dimension(mx,0:my+3,0:mz+3) :: bnd_ind_ex
      logical, dimension(mx,0:my+3,0:mz+3) :: wall_ind_ex
      integer, dimension(1:mx,0:my+3,0:mz+3) :: totbonds
      integer,dimension(3) :: pd_vct_p,pd_vct_n

      ! chance = 1.0d0/dble(mx)*dble(i)
      !chance = 0.1d0
      call random_number(chance)

      pd_vct_p = perind_s(i,j,k,1,1,1,mx,my,mz)
      pd_vct_n = perind_s(i,j,k,-1,-1,-1,mx,my,mz)

      i_p = pd_vct_p(1)
      j_p = pd_vct_p(2)
      k_p = pd_vct_p(3)

      i_n = pd_vct_n(1)
      j_n = pd_vct_n(2)
      k_n = pd_vct_n(3)

      icounter = totbonds(i,j,k)

      ! icounter=0
      ! if(.not.wall_ind(i_p,j  ,k  )) icounter=icounter+1
      ! if(.not.wall_ind(i  ,j_p,k  )) icounter=icounter+1
      ! if(.not.wall_ind(i_n,j  ,k  )) icounter=icounter+1
      ! if(.not.wall_ind(i  ,j_n,k  )) icounter=icounter+1
      ! if(.not.wall_ind(i  ,j  ,k_n)) icounter=icounter+1
      ! if(.not.wall_ind(i  ,j  ,k_p)) icounter=icounter+1

      ! if(icounter.ne.totbonds(i,j,k)) print*,i,j,k

      if(icounter.eq.0) return

      dchance=1.0d0/dble(icounter)

      if(.not.wall_ind_ex(i_n,j,k)) then
        if(chance.le.dchance*icounter.and.chance.gt.dchance*(icounter-1)) then
          bnd_ind_ex(i_n,j,k)=bnd_ind_ex(i_n,j,k) + bnd_ind_ex(i,j,k) - 1.0d0
          bnd_ind_ex(i,j,k) = 1.0d0
          if(i-i_n==2) print*,'False!!!!!!'
        endif
        icounter=icounter-1
      endif 

      if(.not.wall_ind_ex(i_p,j,k)) then
        if(chance.le.dchance*icounter.and.chance.gt.dchance*(icounter-1)) then
          bnd_ind_ex(i_p,j,k)=bnd_ind_ex(i_p,j,k) + bnd_ind_ex(i,j,k) - 1.0d0
          bnd_ind_ex(i,j,k) = 1.0d0
        endif
        icounter=icounter-1
      endif 

      if(.not.wall_ind_ex(i,j_n,k)) then
        if(chance.le.dchance*icounter.and.chance.gt.dchance*(icounter-1)) then
          bnd_ind_ex(i,j_n,k)=bnd_ind_ex(i,j_n,k) + bnd_ind_ex(i,j,k) - 1.0d0
          bnd_ind_ex(i,j,k) = 1.0d0
        endif 
        icounter=icounter-1
      endif

      if(.not.wall_ind_ex(i,j_p,k)) then
        if(chance.le.dchance*icounter.and.chance.gt.dchance*(icounter-1)) then
          bnd_ind_ex(i,j_p,k)=bnd_ind_ex(i,j_p,k) + bnd_ind_ex(i,j,k) - 1.0d0
          bnd_ind_ex(i,j,k) = 1.0d0
        endif      
        icounter=icounter-1
      endif 

      if(.not.wall_ind_ex(i,j,k_n)) then
        if(chance.le.dchance*icounter.and.chance.gt.dchance*(icounter-1)) then
          bnd_ind_ex(i,j,k_n)=bnd_ind_ex(i,j,k_n) + bnd_ind_ex(i,j,k) - 1.0d0
          bnd_ind_ex(i,j,k) = 1.0d0
        endif     
        icounter=icounter-1
      endif

      if(.not.wall_ind_ex(i,j,k_p)) then
        if(chance.le.dchance*icounter.and.chance.gt.dchance*(icounter-1)) then
          bnd_ind_ex(i,j,k_p)=bnd_ind_ex(i,j,k_p) + bnd_ind_ex(i,j,k) - 1.0d0
          bnd_ind_ex(i,j,k) = 1.0d0
        endif
        icounter=icounter-1
      endif 


      end subroutine expand




end module expand_module
