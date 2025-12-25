module solve_module

      integer :: nmax=1000

      contains

      subroutine solve (nc,r,amat)
      
      implicit none

      integer :: nc,m,i,i1,imax,j,k,l
      real*8 :: r,amat,v,zmax

      dimension r(*),amat(nc,*)

!-----elimination
      do m = 1, nc-1
        zmax = 0.d0 
        imax = 0

!-------find max of column    
        do i  = m, nc 
          if(abs(amat(i,m)).gt.zmax) then
            imax = i
            zmax = abs(amat(i,m))
          endif
        enddo

        if(imax .eq. 0) goto 20

!-------row interchange       
        if(imax.eq.m) goto 10
        v = r(m)      
        r(m) = r(imax)
        r(imax) = v   

        do j = m, nc
          v = amat(m,j)    
          amat(m,j) = amat(imax,j)      
          amat(imax,j) = v 
        enddo

   10   continue

!-------diagonalize
        do i = m+1, nc
          v = amat(i,m)/amat(m,m)
          r(i) = r(i)-v*r(m)
          do j = m, nc
            amat(i,j) = amat(i,j)-v*amat(m,j)
          enddo
        enddo
      enddo

!-----back substitute
      r(nc)=r(nc)/amat(nc,nc)
      do k = 1, nc-1 
        i = nc-k
        i1 = i+1
        do j = i1, nc
          r(i) = r(i)-r(j)*amat(i,j)
        enddo
        r(i) = r(i)/amat(i,i)
      enddo

      return

   20 continue
      write(*,'("singular jacobian matrix in solve.f: stop")') 
      do j = 1, nc
        write(*,'(a12,1pe12.4/(1p10e12.4))') j,r(j),(amat(j,l),l=1,nc)
      enddo
      stop ' in solve.f'

      end subroutine solve

      subroutine ludcmp(a,n,np,indx,d)
      
      implicit none

!     parameter (nmax=1000,tiny=1.d-20)

      integer :: i,imax,indx,j,k,n,np
      real*8 :: a,aamax,d,dum,sum,tiny=1.d-20,vv
      
      dimension a(np,np),indx(n),vv(nmax)

      d=1.d0
      do i=1,n
        aamax=0.d0
        do j=1,n
          if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
        enddo
        if (aamax.eq.0.d0) pause 'singular matrix. [ludcmp]'
        vv(i)=1.d0/aamax
      enddo
      do j=1,n
        if (j.gt.1) then
          do i=1,j-1
            sum=a(i,j)
            if (i.gt.1)then
              do k=1,i-1
                sum=sum-a(i,k)*a(k,j)
              enddo
              a(i,j)=sum
            endif
          enddo
        endif
        aamax=0.d0
        do i=j,n
          sum=a(i,j)
          if (j.gt.1)then
            do k=1,j-1
              sum=sum-a(i,k)*a(k,j)
            enddo
            a(i,j)=sum
          endif
          dum=vv(i)*dabs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
        enddo
        if (j.ne.imax)then
          do k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          enddo
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(j.ne.n)then
          if(a(j,j).eq.0.d0)a(j,j)=tiny
          dum=1.d0/a(j,j)
          do i=j+1,n
            a(i,j)=a(i,j)*dum
          enddo
        endif
      enddo
      if(a(n,n).eq.0.d0)a(n,n)=tiny
      return
      end subroutine ludcmp

      subroutine lubksb(a,n,np,indx,b)

      implicit none
      
      integer :: ii,i,indx,j,ll,n,np
      real*8 :: sum,a,b

      dimension a(np,np),indx(n),b(n)

      ii=0
      do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do j=ii,i-1
            sum=sum-a(i,j)*b(j)
          enddo
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
      enddo
      do i=n,1,-1
        sum=b(i)
        if(i.lt.n)then
          do j=i+1,n
            sum=sum-a(i,j)*b(j)
          enddo
        endif
        b(i)=sum/a(i,i)
      enddo
      return
      end subroutine lubksb

      subroutine mprove (a,alud,n,np,indx,b,x)

      implicit none
      
!     parameter (nmax=1000)

      integer :: i,indx,j,n,np
      real*8 :: a,alud,b,x,r,sdp,dble
      
      dimension a(np,np),alud(np,np),indx(n),b(n),x(n),r(nmax)

      do i=1,n
        sdp=-b(i)
        do j=1,n
          sdp=sdp+dble(a(i,j))*dble(x(j))
        enddo
        r(i)=sdp
      enddo
      call lubksb(alud,n,np,indx,r)
      do i=1,n
        x(i)=x(i)-r(i)
      enddo
      return
      end subroutine mprove
      
end module solve_module
