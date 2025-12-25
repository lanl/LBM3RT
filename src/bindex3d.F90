module bindex_module

contains
    subroutine bindex(walls,bonds,nbonds,totbonds)

    use comvar_module
    use perind_module

    implicit none

    integer :: i,j,k,s,liq_num,i_p,j_p,k_p,i_n,j_n,k_n
    integer,dimension(3) :: pd_vct,pd_vct_p,pd_vct_n

    logical, dimension(1:mx,0:my+3,0:mz+3) :: walls,bonds
    integer, dimension(1:mx,0:my+3,0:mz+3) :: totbonds !number of liq neighbors
    integer, dimension(6,1:mx,0:my+3,0:mz+3) :: nbonds !orientation index for fluids: dim(1(x+), 2(x-), 3(y+), 4(y-), 5(z+), 6(z-)),1-true,0-false

    bonds = .false.
    totbonds = 0

    do i = 1,mx
      do j = 1,my+2
        do k = 1,mz+2
          if(walls(i,j,k)) then
            liq_num = 0
            do s = 1,6
              pd_vct = perind_s(i,j,k,cjxi(s),cjyi(s),cjzi(s),mx,my,mz)
              if(.not.walls(pd_vct(1),pd_vct(2),pd_vct(3))) then
                liq_num = liq_num + 1
              endif
            enddo
            totbonds(i,j,k) = liq_num 

            if (liq_num > 0) then
              bonds(i,j,k) = .true.
            endif

            pd_vct_p = perind_s(i,j,k,1,1,1,mx,my,mz)
            pd_vct_n = perind_s(i,j,k,-1,-1,-1,mx,my,mz)

            i_p = pd_vct_p(1)
            j_p = pd_vct_p(2)
            k_p = pd_vct_p(3)

            i_n = pd_vct_n(1)
            j_n = pd_vct_n(2)
            k_n = pd_vct_n(3)

            if (bonds(i,j,k) .eqv. .true.) then
              if(.not.walls(i_p,j,k)) then
                nbonds(1,i,j,k) = 1
              else
                nbonds(1,i,j,k) = 0
              endif

              if(.not.walls(i_n,j,k)) then
                nbonds(2,i,j,k) = 1
              else
                nbonds(2,i,j,k) = 0
              endif

              if(.not.walls(i,j_p,k)) then
                nbonds(3,i,j,k) = 1
              else
                nbonds(3,i,j,k) = 0
              endif

              if(.not.walls(i,j_n,k)) then
                nbonds(4,i,j,k) = 1
              else
                nbonds(4,i,j,k) = 0
              endif

              if(.not.walls(i,j,k_p)) then
                nbonds(5,i,j,k) = 1
              else
                nbonds(5,i,j,k) = 0
              endif

              if(.not.walls(i,j,k_n)) then
                nbonds(6,i,j,k) = 1
              else
                nbonds(6,i,j,k) = 0
              endif
            endif

          endif
        enddo
      enddo
    enddo




    end subroutine bindex


    subroutine bindex2(walls,lbuffers)

    use comvar_module
    use perind_module

    implicit none

    integer :: i,j,k,s,liq_num,i_p,j_p,k_p,i_n,j_n,k_n
    integer,dimension(3) :: pd_vct,pd_vct_p,pd_vct_n

    logical, dimension(1:mx,0:my+3,0:mz+3) :: walls,lbuffers!lbuffers is the fluid-side bonds matrix


    lbuffers=.false.

    do i = 1,mx
      do j = 1,my+2
        do k = 1,mz+2
          if(.not.walls(i,j,k)) then

            pd_vct_p = perind_s(i,j,k,1,1,1,mx,my,mz)
            pd_vct_n = perind_s(i,j,k,-1,-1,-1,mx,my,mz)

            i_p = pd_vct_p(1)
            j_p = pd_vct_p(2)
            k_p = pd_vct_p(3)

            i_n = pd_vct_n(1)
            j_n = pd_vct_n(2)
            k_n = pd_vct_n(3)

            if(walls(i_p,j,k)) then
              lbuffers(i,j,k)=.true.
            endif

            if(walls(i_n,j,k)) then
              lbuffers(i,j,k)=.true.
            endif

            if(walls(i,j_p,k)) then
              lbuffers(i,j,k)=.true.
            endif   

            if(walls(i,j_n,k)) then
              lbuffers(i,j,k)=.true.
            endif   

            if(walls(i,j,k_p)) then
              lbuffers(i,j,k)=.true.
            endif  

            if(walls(i,j,k_n)) then
              lbuffers(i,j,k)=.true.
            endif  
          endif
        enddo
      enddo
    enddo

    end subroutine bindex2


end module bindex_module