module redistri_bnd_module

contains

      subroutine redistri_bnd(bnd,walls,totbonds,i,j,k,bsat,p,nkin)
      
      use comvar_module
      use perind_module
      
      implicit none
      
      integer :: icounter,i,j,k,i_p,j_p,k_p,i_n,j_n,k_n,p,pp,nkin
      real*8 :: bsat
      
      logical, dimension(mx,0:my+3,0:mz+3) :: walls
      real*8, dimension(1:nkin,1:mx,0:my+3,0:mz+3) :: bnd
      integer, dimension(1:mx,0:my+3,0:mz+3) :: totbonds
      real*8, dimension(1:nkin) :: dbnd
      integer,dimension(3) :: pd_vct_p,pd_vct_n

      pd_vct_p = perind_s(i,j,k,1,1,1,mx,my,mz)
      pd_vct_n = perind_s(i,j,k,-1,-1,-1,mx,my,mz)

      i_p = pd_vct_p(1)
      j_p = pd_vct_p(2)
      k_p = pd_vct_p(3)

      i_n = pd_vct_n(1)
      j_n = pd_vct_n(2)
      k_n = pd_vct_n(3)
      
      !icounter=totbonds(i,j,k)
      icounter = 0
      if(walls(i_p,j,k).and.bnd(p,i_p,j,k)>0) icounter=icounter+1
      if(walls(i_n,j,k).and.bnd(p,i_n,j,k)>0) icounter=icounter+1
      if(walls(i,j_p,k).and.bnd(p,i,j_p,k)>0) icounter=icounter+1
      if(walls(i,j_n,k).and.bnd(p,i,j_n,k)>0) icounter=icounter+1
      if(walls(i,j,k_p).and.bnd(p,i,j,k_p)>0) icounter=icounter+1
      if(walls(i,j,k_n).and.bnd(p,i,j,k_n)>0) icounter=icounter+1

      if(icounter.eq.0) return

      do pp=1,nkin
        if(pp.ne.p) then
          dbnd(pp)=(bnd(pp,i,j,k)-bsat)/float(icounter)

          if(walls(i_p,j,k).and.bnd(p,i_p,j,k)>0) bnd(pp,i_p,j,k)=bnd(pp,i_p,j,k)+dbnd(pp)
          if(walls(i_n,j,k).and.bnd(p,i_n,j,k)>0) bnd(pp,i_n,j,k)=bnd(pp,i_n,j,k)+dbnd(pp)
          if(walls(i,j_p,k).and.bnd(p,i,j_p,k)>0) bnd(pp,i,j_p,k)=bnd(pp,i,j_p,k)+dbnd(pp)
          if(walls(i,j_n,k).and.bnd(p,i,j_n,k)>0) bnd(pp,i,j_n,k)=bnd(pp,i,j_n,k)+dbnd(pp)
          if(walls(i,j,k_p).and.bnd(p,i,j,k_p)>0) bnd(pp,i,j,k_p)=bnd(pp,i,j,k_p)+dbnd(pp)
          if(walls(i,j,k_n).and.bnd(p,i,j,k_n)>0) bnd(pp,i,j,k_n)=bnd(pp,i,j,k_n)+dbnd(pp)

          bnd(pp,i,j,k) = bsat
        endif
      enddo        


      end subroutine redistri_bnd


      
      subroutine extra_relaxation(gi,totbonds,nbonds,l1,l2,l3,ncomp)
      use comvar_module

      integer :: ncomp
      real*8, dimension(1:ncomp,1:mx,1:my+2,1:mz+2,0:6) :: gi
      real*8, dimension(1:ncomp) :: gie
      integer, dimension(1:mx,0:my+3,0:mz+3) :: totbonds
      integer, dimension(6,1:mx,0:my+3,0:mz+3) :: nbonds
      integer :: l1,l2,l3,pp

      if(totbonds(l1,l2,l3)==1) then
        if(nbonds(1,l1,l2,l3)==1) then
          gie = gi(:,l1,l2,l3,2)
        endif
        if(nbonds(2,l1,l2,l3)==1) then
          gie = gi(:,l1,l2,l3,1)
        endif
        if(nbonds(3,l1,l2,l3)==1) then
          gie = gi(:,l1,l2,l3,4)
        endif
        if(nbonds(4,l1,l2,l3)==1) then
          gie = gi(:,l1,l2,l3,3)
        endif
        if(nbonds(5,l1,l2,l3)==1) then
          gie = gi(:,l1,l2,l3,6)
        endif
        if(nbonds(6,l1,l2,l3)==1) then
          gie = gi(:,l1,l2,l3,5)
        endif

        do pp = 1,6
          gi(:,l1,l2,l3,pp) = gie(:)
        enddo
        gi(:,l1,l2,l3,0) = gie * 2.d0
      endif


      if(totbonds(l1,l2,l3)==2) then
        gie = 0.d0
        if(nbonds(1,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,2)
        endif
        if(nbonds(2,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,1)
        endif
        if(nbonds(3,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,4)
        endif
        if(nbonds(4,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,3)
        endif
        if(nbonds(5,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,6)
        endif
        if(nbonds(6,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,5)
        endif
        

        do pp = 1,6
          gi(:,l1,l2,l3,pp) = gie(:)/2.d0
        enddo
        gi(:,l1,l2,l3,0) = gie
      endif

      if(totbonds(l1,l2,l3)==3) then
        gie = 0.d0
        if(nbonds(1,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,2)
        endif
        if(nbonds(2,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,1)
        endif
        if(nbonds(3,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,4)
        endif
        if(nbonds(4,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,3)
        endif
        if(nbonds(5,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,6)
        endif
        if(nbonds(6,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,5)
        endif
        

        do pp = 1,6
          gi(:,l1,l2,l3,pp) = gie(:)/3.d0
        enddo
        gi(:,l1,l2,l3,0) = gie/3.d0*2.d0
      endif

      if(totbonds(l1,l2,l3)==4) then
        gie = 0.d0
        if(nbonds(1,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,2)
        endif
        if(nbonds(2,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,1)
        endif
        if(nbonds(3,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,4)
        endif
        if(nbonds(4,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,3)
        endif
        if(nbonds(5,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,6)
        endif
        if(nbonds(6,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,5)
        endif
        

        do pp = 1,6
          gi(:,l1,l2,l3,pp) = gie(:)/4.d0
        enddo
        gi(:,l1,l2,l3,0) = gie/4.d0*2.d0
      endif


      if(totbonds(l1,l2,l3)==5) then
        gie = 0.d0
        if(nbonds(1,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,2)
        endif
        if(nbonds(2,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,1)
        endif
        if(nbonds(3,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,4)
        endif
        if(nbonds(4,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,3)
        endif
        if(nbonds(5,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,6)
        endif
        if(nbonds(6,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,5)
        endif
        

        do pp = 1,6
          gi(:,l1,l2,l3,pp) = gie(:)/5.d0
        enddo
        gi(:,l1,l2,l3,0) = gie/5.d0*2.d0
      endif

      if(totbonds(l1,l2,l3)==6) then
        gie = 0.d0
        if(nbonds(1,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,2)
        endif
        if(nbonds(2,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,1)
        endif
        if(nbonds(3,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,4)
        endif
        if(nbonds(4,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,3)
        endif
        if(nbonds(5,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,6)
        endif
        if(nbonds(6,l1,l2,l3)==1) then
          gie = gie + gi(:,l1,l2,l3,5)
        endif
        

        do pp = 1,6
          gi(:,l1,l2,l3,pp) = gie(:)/6.d0
        enddo
        gi(:,l1,l2,l3,0) = gie/6.d0*2.d0
      endif

      end subroutine extra_relaxation



      ! subroutine mass_supplement(cj,psi,gi,i,j,k,walls,uxe,uye,uze)
      ! use comvar_module
      ! use perind_module
      ! use equilg_module

      ! logical, dimension(mx,1:my+2,1:mz+2) :: walls
      ! real*8, dimension(1:m,1:mx,1:my+2,1:mz+2) :: cj,psi
      ! real*8, dimension(1:m,1:mx,1:my+2,1:mz+2,0:6) :: gi
      ! integer,dimension(3) :: pd_vct_p,pd_vct_n
      ! integer :: icounter,i,j,k,i_p,j_p,k_p,i_n,j_n,k_n,pp
      ! real*8 :: uxe,uye,uze


      ! pd_vct_p = perind(i,j,k,1,1,1,mx,my,mz)
      ! pd_vct_n = perind(i,j,k,-1,-1,-1,mx,my,mz)

      ! icounter=0
      ! if(walls(i_p,j,k)) icounter=icounter+1
      ! if(walls(i_n,j,k)) icounter=icounter+1
      ! if(walls(i,j_p,k)) icounter=icounter+1
      ! if(walls(i,j_n,k)) icounter=icounter+1
      ! if(walls(i,j,k_p)) icounter=icounter+1
      ! if(walls(i,j,k_n)) icounter=icounter+1

      ! if(icounter.eq.0) return

      ! cj(:,i,j,k) = 0.0d0
      ! psi(:,i,j,k) = 0.0d0


      ! if(.not.walls(i_p,j,k)) then
      !   do pp=1,m
      !     cj(pp,i,j,k) = cj(pp,i,j,k) + cj(pp,i_p,j,k)/dble(icounter)
      !     psi(pp,i,j,k) = psi(pp,i,j,k) + psi(pp,i_p,j,k)/dble(icounter)
      !   enddo
      ! endif

      ! if(.not.walls(i_n,j,k)) then
      !   do pp=1,m
      !     cj(pp,i,j,k) = cj(pp,i,j,k) + cj(pp,i_n,j,k)/dble(icounter)
      !     psi(pp,i,j,k) = psi(pp,i,j,k) + psi(pp,i_n,j,k)/dble(icounter)
      !   enddo
      ! endif

      ! if(.not.walls(i,j_p,k)) then
      !   do pp=1,m
      !     cj(pp,i,j,k) = cj(pp,i,j,k) + cj(pp,i,j_p,k)/dble(icounter)
      !     psi(pp,i,j,k) = psi(pp,i,j,k) + psi(pp,i,j_p,k)/dble(icounter)
      !   enddo
      ! endif

      ! if(.not.walls(i,j_n,k)) then
      !   do pp=1,m
      !     cj(pp,i,j,k) = cj(pp,i,j,k) + cj(pp,i,j_n,k)/dble(icounter)
      !     psi(pp,i,j,k) = psi(pp,i,j,k) + psi(pp,i,j_n,k)/dble(icounter)
      !   enddo
      ! endif

      ! if(.not.walls(i,j,k_p)) then
      !   do pp=1,m
      !     cj(pp,i,j,k) = cj(pp,i,j,k) + cj(pp,i,j,k_p)/dble(icounter)
      !     psi(pp,i,j,k) = psi(pp,i,j,k) + psi(pp,i,j,k_p)/dble(icounter)
      !   enddo
      ! endif

      ! if(.not.walls(i,j,k_n)) then
      !   do pp=1,m
      !     cj(pp,i,j,k) = cj(pp,i,j,k) + cj(pp,i,j,k_n)/dble(icounter)
      !     psi(pp,i,j,k) = psi(pp,i,j,k) + psi(pp,i,j,k_n)/dble(icounter)
      !   enddo
      ! endif

      ! call equilg(gi(p,i,j,k,:),psi(p,i,j,k),uxe,uye,uze)



      ! end subroutine mass_supplement




      ! subroutine interpolate_bonds(cj,psi,totbonds,bonds,walls)

      ! use comvar_module
      ! use perind_module

      ! real*8, dimension(1:m,1:mx,1:my+2,1:mz+2) :: cj,psi
      ! integer, dimension(1:mx,0:my+3,0:mz+3) :: totbonds
      ! logical, dimension(1:mx,0:my+3,0:mz+3) :: bonds,walls


      ! do i=1,mx
      !   do j=2,my+1
      !     do k=2,mz+1

      !       if(bonds(i,j,k)) then

      !         pd_vct_p = perind_s(i,j,k,1,1,1,mx,my,mz)
      !         pd_vct_n = perind_s(i,j,k,-1,-1,-1,mx,my,mz)

      !         i_p = pd_vct_p(1)
      !         j_p = pd_vct_p(2)
      !         k_p = pd_vct_p(3)

      !         i_n = pd_vct_n(1)
      !         j_n = pd_vct_n(2)
      !         k_n = pd_vct_n(3)
              
      !         icounter=totbonds(i,j,k)


      !         cj(:,i,j,k) = 0.0d0
      !         psi(:,i,j,k) = 0.0d0


      !         if(.not.walls(i_p,j,k)) then
      !           do pp=1,m
      !             cj(pp,i,j,k) = cj(pp,i,j,k) + cj(pp,i_p,j,k)/dble(icounter)
      !             psi(pp,i,j,k) = psi(pp,i,j,k) + psi(pp,i_p,j,k)/dble(icounter)
      !           enddo
      !         endif

      !         if(.not.walls(i_n,j,k)) then
      !           do pp=1,m
      !             cj(pp,i,j,k) = cj(pp,i,j,k) + cj(pp,i_n,j,k)/dble(icounter)
      !             psi(pp,i,j,k) = psi(pp,i,j,k) + psi(pp,i_n,j,k)/dble(icounter)
      !           enddo
      !         endif

      !         if(.not.walls(i,j_p,k)) then
      !           do pp=1,m
      !             cj(pp,i,j,k) = cj(pp,i,j,k) + cj(pp,i,j_p,k)/dble(icounter)
      !             psi(pp,i,j,k) = psi(pp,i,j,k) + psi(pp,i,j_p,k)/dble(icounter)
      !           enddo
      !         endif

      !         if(.not.walls(i,j_n,k)) then
      !           do pp=1,m
      !             cj(pp,i,j,k) = cj(pp,i,j,k) + cj(pp,i,j_n,k)/dble(icounter)
      !             psi(pp,i,j,k) = psi(pp,i,j,k) + psi(pp,i,j_n,k)/dble(icounter)
      !           enddo
      !         endif

      !         if(.not.walls(i,j,k_p)) then
      !           do pp=1,m
      !             cj(pp,i,j,k) = cj(pp,i,j,k) + cj(pp,i,j,k_p)/dble(icounter)
      !             psi(pp,i,j,k) = psi(pp,i,j,k) + psi(pp,i,j,k_p)/dble(icounter)
      !           enddo
      !         endif

      !         if(.not.walls(i,j,k_n)) then
      !           do pp=1,m
      !             cj(pp,i,j,k) = cj(pp,i,j,k) + cj(pp,i,j,k_n)/dble(icounter)
      !             psi(pp,i,j,k) = psi(pp,i,j,k) + psi(pp,i,j,k_n)/dble(icounter)
      !           enddo
      !         endif

      !       endif
      !     enddo
      !   enddo
      ! enddo

      ! end subroutine interpolate_bonds

      ! subroutine equig_bond(gi,psi,walls,i,j,k)
      ! use comvar_module
      ! use perind_module
      ! use equilg_module

      ! logical, dimension(1:mx,0:my+3,0:mz+3) :: walls
      ! real*8, dimension(1:m,1:mx,1:my+2,1:mz+2,0:6) :: gi
      ! real*8, dimension(1:m,1:mx,1:my+2,1:mz+2) :: psi
      ! integer,dimension(3) :: pd_vct_p,pd_vct_n
      ! integer :: i,j,k,i_p,j_p,k_p,i_n,j_n,k_n,pp



      ! pd_vct_p = perind(i,j,k,1,1,1,mx,my,mz)
      ! pd_vct_n = perind(i,j,k,-1,-1,-1,mx,my,mz)

      ! ! if(walls(i_p,j,k)) then
      ! !   do pp=1,m
      ! !     call equilg(gi(i_p,i,j,k,:),psi(p,i_p,j,k),0.d0,0.d0,0.d0)
      ! !   enddo
      ! ! endif

      ! ! if(walls(i_n,j,k)) then
      ! !   do pp=1,m
      ! !     call equilg(gi(p,i_n,j,k,:),psi(p,i_n,j,k),0.d0,0.d0,0.d0)
      ! !   enddo
      ! ! endif

      ! ! if(walls(i,j_p,k)) then
      ! !   do pp=1,m
      ! !     call equilg(gi(p,i,j_p,k,:),psi(p,i,j_p,k),0.d0,0.d0,0.d0)
      ! !   enddo
      ! ! endif

      ! ! if(walls(i,j_n,k)) then
      ! !   do pp=1,m
      ! !     call equilg(gi(p,i,j_n,k,:),psi(p,i,j_n,k),0.d0,0.d0,0.d0)
      ! !   enddo
      ! ! endif

      ! ! if(walls(i,j,k_p)) then
      ! !   do pp=1,m
      ! !     call equilg(gi(p,i,j,k_p,:),psi(p,i,j,k_p),0.d0,0.d0,0.d0)
      ! !   enddo
      ! ! endif


      ! ! if(walls(i,j,k_n)) then
      ! !   do pp=1,m
      ! !     call equilg(gi(p,i,j,k_n,:),psi(p,i,j,k_n),0.d0,0.d0,0.d0)
      ! !   enddo
      ! ! endif

      ! end subroutine equig_bond




end module redistri_bnd_module
