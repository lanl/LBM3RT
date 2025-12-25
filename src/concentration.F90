  module con_module

      contains

      subroutine con_cj(psi,cj,walls,bonds,mcyc,itermx,istep)!,lbuffers)

      use ptran_global_module
      use solve_module
      use perind_module

      use comvar_module
      use mpi
      
      implicit none

      real*8, dimension(1:ncomp,mx,my+2,mz+2) :: cj,psi
      logical, dimension(mx,0:my+3,0:mz+3) :: walls,bonds!,lbuffers
      
      real*8 :: alogjn(ncomp),r(ncomp),dpsi(ncomp,ncomp),cx(ncmplx)
      
      integer :: l1,l2,l3,i,iter,itmax,j,k,mcyc,indx(ncomp),itermx,istep
      real*8 :: delcmx,prod1,sumc,udelc,ucc,dd,dcc
      integer,dimension(3) :: pd_vct_p,pd_vct_n
      integer :: i_p,j_p,k_p,i_n,j_n,k_n

      
      if (ncmplx == 0) then
        cj=psi
        return
      endif
      
      tol  = 1.d-10
      itmax = 100
      icut = 0
!     loglin = 0 ! solve for log C
      loglin = 1 ! solve for  C
      iwarn = 2
      tolexp = 10.d0
      itermx = 0

!     begin loop over nodes

      do l1 = 1, mx
        do l2 = 1, my+2
          do l3 = 1, mz+2
            if(.not.walls(l1,l2,l3)) then
            !if(.not.walls(l1,l2,l3).and.(.not.lbuffers(l1,l2,l3))) then
              ! if(istep>29790.and.myrank==104.and.l1==11.and.l2==4.and.l3==4) print*,'illegal!',lbuffers(l1,l2,l3),walls(l1,l2,l3)
              iter = 0
            
    10       continue
              iter = iter + 1

              do j = 1, ncomp
  !             alogjn(j) = log(cj(j,l1,l2)*gam(j,n))
                alogjn(j) = log(cj(j,l1,l2,l3))
              enddo

              do i = 1, ncmplx
  !             prod1 = eqhom(i)*aln10-log(gamx(i,n))
                prod1 = eqhom(i)*aln10
                do k = 1, ncomp
                  prod1 = prod1 + shom(k,i)*alogjn(k)
                enddo
                cx(i) = exp(prod1)
              enddo

  !-----------compute psi for primary species
              do j = 1, ncomp 
                sumc = cj(j,l1,l2,l3)
                do i = 1, ncmplx
                  sumc = sumc + shom(j,i)*cx(i) 
                enddo
  !             sumc = sumc*rho(n)
                r(j) = sumc - psi(j,l1,l2,l3)
  !             print *,'sum-r: ',iter,j,r(j),sumc,psi(j,l1,l2)
          
  !             diagonal terms
                sumc = cj(j,l1,l2,l3)
  !             print *,'sum: ',j,sum,r(j)
                do i = 1, ncmplx
                  sumc = sumc + shom(j,i)*shom(j,i)*cx(i)
  !               print *,'sum: ',i,cx(i),shom(j,i),sumc
                enddo
  !             dpsi(j,j) = sumc*rho(n)
                dpsi(j,j) = sumc
              
  !             print *,'diag: ',j,dpsi(j,j),cj(j,l1,l2)

  !             off-diagonal terms
                do k = j+1, ncomp
                  sumc = 0.d0
                  do i = 1, ncmplx
                    sumc = sumc + shom(j,i)*shom(k,i)*cx(i)
                  enddo
  !               dpsi(j,k) = sumc*rho(n)
                  dpsi(j,k) = sumc
                  dpsi(k,j) = sumc
                enddo
              enddo
        
              if (loglin == 1) then
                do j = 1, ncomp
                  do k = 1, ncomp
                    dpsi(j,k) = dpsi(j,k)/cj(k,l1,l2,l3)
                  enddo
                enddo
              endif
            
  !           check convergence
              if (iter > 1) then
                do j = 1, ncomp
                  if (abs(r(j)) > tol) goto 20
                enddo
        
  !             stop
          
                goto 30
              endif
            
    20       continue
              
              if (mcyc==-2) then
                print *,'res0: ',mcyc,iter,i1,i2,(r(j),j=1,ncomp)
                print *,'cj: ',mcyc,iter,i1,i2,(cj(j,l1,l2,l3),j=1,ncomp)
                print *,'psi: ',mcyc,iter,(psi(j,l1,l2,l3),j=1,ncomp)
                do j = 1, ncomp
                  print *,'diag: ',(dpsi(j,k),k=1,ncomp)
                enddo
                ! WRITE(*, '(A, 1X, I0)') 'checkpoint1: step =',istep
                ! WRITE(*, '(A, 1X, I0)') 'checkpoint1: myrank =',myrank
              endif
              
              call solve (ncomp,r,dpsi)
  !           call ludcmp(dpsi,ncomp,ncomp,indx,dd)
  !           call lubksb(dpsi,ncomp,ncomp,indx,r)
              
  !           print *,'res: ',mcyc,iter,i1,i2,(r(j),j=1,ncomp)

  !-----------update field variables corresponding to new iteration
              delcmx = 0.d0
              if (loglin .eq. 0) then           ! solution is log(delc) 
                do j = 1, ncomp
                  ucc = cj(j,l1,l2,l3)
  !               udelc = -r(j)
                  udelc = -r(j) + log(ucc)
                  if (abs(udelc) .gt. abs(tolexp)) then
                    if (iwarn.ge.2) write (*,100) mcyc,l1,l2,l3,iter,icut, &
                    nam(j),cj(j,l1,l2,l3),udelc
    100             format ('OSRXN: delc > tolexp, mcyc =',i5,' n = ',2i4, &
      &            ' iter =',i3,' icuts =',i2,1x,a8, &
      &            '=',1pe12.4,' delc =',1pe12.4)
                    ! WRITE(*, '(A, 1X, I0)') 'checkpoint2: step =',istep
                    ! WRITE(*, '(A, 1X, I0)') 'checkpoint2: myrank =',myrank
                    if (tolexp.lt.zero) then
  !                   goto 70  ! cut time step
                    else
                      udelc = tolexp*udelc/abs(udelc)
                    endif
                  endif
  !               ucc = cj(j,l1,l2)
  !               udelc = udelc + log(ucc)
  !               cj(j,l1,l2) = ucc*exp(udelc)
                  cj(j,l1,l2,l3) = exp(udelc)
                  udelc = cj(j,l1,l2,l3)-ucc
                  if(abs(udelc).gt.abs(delcmx)) delcmx = udelc
                enddo

              else                           ! solution is delc

                do j = 1, ncomp 
                  udelc = -r(j)
                  if(abs(udelc).gt.abs(delcmx)) delcmx = udelc !误差取最大值
                  ucc = cj(j,l1,l2,l3)
  !               cj(j,l1,l2) = ucc+udelc
                  dcc = max (one, -two*udelc/ucc)
                  cj(j,l1,l2,l3) = ucc+udelc/dcc
                  if(cj(j,l1,l2,l3).le.0.d0) then ! check for negative conc
                    cj(j,l1,l2,l3) = ucc
                    if (iwarn.ge.2) &
                    write(*,'("mcyc= ",i4," iter= ",i3, &
      &            " OSRXN:--negative or zero concentration for species: ", &
      &            a8," at node: ",2i3,1p10e12.4)') &
                    mcyc,iter,nam(j),l1,l2,l3,udelc,ucc
                    ! WRITE(*, '(A, 1X, I0)') 'checkpoint3: step =',istep
                    ! WRITE(*, '(A, 1X, I0)') 'checkpoint3: myrank =',myrank
  !                 goto 70
                  endif
                enddo
              endif

  !           print *,'cj1:   ',mcyc,iter,(cj(j,l1,l2),j=1,ncomp)
  !           print *,'psi1: ',mcyc,iter,(psi(j,l1,l2),j=1,ncomp)

              if (iter < itmax) goto 10
              print *,'Max iters at node: (',myrank,', ',l1,', ',l2,', ',l3,')',' its = ', &
              iter,' mcyc = ',istep


              print *,'Res = ',(r(j),j=1,ncomp)
              print *,'cj = ',(cj(j,l1,l2,l3),j=1,ncomp)
              print *,'psi = ',(psi(j,l1,l2,l3),j=1,ncomp)
              do j = 1, ncomp
                print *,'dpsi = ',(dpsi(j,k),k=1,ncomp)
              enddo
              ! WRITE(*, '(A, 1X, I0)') 'checkpoint4: step =',istep
              ! WRITE(*, '(A, 1X, I0)') 'checkpoint4: myrank =',myrank
              CALL MPI_ABORT(MPI_COMM_WORLD, 0, ierr)
    30       continue
              if (iter > itermx) itermx = iter
  !           write(*,'(2i3,(a8,1pe12.4))') l1,l2,(nam(j),r(j),j=1,ncomp)
            endif
          enddo
        enddo
      enddo

      end subroutine con_cj

!--------------------------

!-----------------------------------------------------

      subroutine con_cj_ind(psi_ind,cj_ind,istep,l1,l2,l3)!,lbuffers)

      use ptran_global_module
      use solve_module
      use perind_module

      use comvar_module
      use mpi
      
      implicit none

      real*8, dimension(1:ncomp) :: cj_ind,psi_ind

      
      real*8 :: alogjn(ncomp),r(ncomp),dpsi(ncomp,ncomp),cx(ncmplx)
      
      integer :: i,iter,itmax,j,k,indx(ncomp),itermx,istep,mcyc,l1,l2,l3
      real*8 :: delcmx,prod1,sumc,udelc,ucc,dd,dcc
      integer,dimension(3) :: pd_vct_p,pd_vct_n
      integer :: i_p,j_p,k_p,i_n,j_n,k_n

      
      mcyc = istep
      if (ncmplx == 0) then
        cj_ind=psi_ind
        return
      endif
      
      tol  = 1.d-10
      itmax = 105
      icut = 0
!     loglin = 0 ! solve for log C
      loglin = 1 ! solve for  C
      iwarn = 2
      tolexp = 10.d0
      itermx = 0

!     begin loop over nodes


              iter = 0
            
    10       continue
              iter = iter + 1

              do j = 1, ncomp
                alogjn(j) = log(cj_ind(j))
              enddo

              do i = 1, ncmplx
  !             prod1 = eqhom(i)*aln10-log(gamx(i,n))
                prod1 = eqhom(i)*aln10
                do k = 1, ncomp
                  prod1 = prod1 + shom(k,i)*alogjn(k)
                enddo
                cx(i) = exp(prod1)
              enddo

  !-----------compute psi for primary species
              do j = 1, ncomp 
                sumc = cj_ind(j)
                do i = 1, ncmplx
                  sumc = sumc + shom(j,i)*cx(i) 
                enddo
  !             sumc = sumc*rho(n)
                r(j) = sumc - psi_ind(j)
          
  !             diagonal terms
                sumc = cj_ind(j)
  !             print *,'sum: ',j,sum,r(j)
                do i = 1, ncmplx
                  sumc = sumc + shom(j,i)*shom(j,i)*cx(i)
  !               print *,'sum: ',i,cx(i),shom(j,i),sumc
                enddo
  !             dpsi(j,j) = sumc*rho(n)
                dpsi(j,j) = sumc


  !             off-diagonal terms
                do k = j+1, ncomp
                  sumc = 0.d0
                  do i = 1, ncmplx
                    sumc = sumc + shom(j,i)*shom(k,i)*cx(i)
                  enddo
  !               dpsi(j,k) = sumc*rho(n)
                  dpsi(j,k) = sumc
                  dpsi(k,j) = sumc
                enddo
              enddo
        
              if (loglin == 1) then
                do j = 1, ncomp
                  do k = 1, ncomp
                    dpsi(j,k) = dpsi(j,k)/cj_ind(k)
                  enddo
                enddo
              endif
            
  !           check convergence
              if (iter > 1) then
                do j = 1, ncomp
                  if (abs(r(j)) > tol) goto 20
                enddo
        
  !             stop
          
                goto 30
              endif
            
    20       continue
              
              if (mcyc==-2) then
                print *,'res0: ',mcyc,iter,i1,i2,(r(j),j=1,ncomp)
                print *,'cj: ',mcyc,iter,i1,i2,(cj_ind(j),j=1,ncomp)
                print *,'psi: ',mcyc,iter,(psi_ind(j),j=1,ncomp)
                do j = 1, ncomp
                  print *,'diag: ',(dpsi(j,k),k=1,ncomp)
                enddo
                ! WRITE(*, '(A, 1X, I0)') 'checkpoint1: step =',istep
                ! WRITE(*, '(A, 1X, I0)') 'checkpoint1: myrank =',myrank
              endif
              
              call solve (ncomp,r,dpsi)
  !           call ludcmp(dpsi,ncomp,ncomp,indx,dd)
  !           call lubksb(dpsi,ncomp,ncomp,indx,r)
              
  !           print *,'res: ',mcyc,iter,i1,i2,(r(j),j=1,ncomp)

  !-----------update field variables corresponding to new iteration
              delcmx = 0.d0
              if (loglin .eq. 0) then           ! solution is log(delc) 
                do j = 1, ncomp
                  ucc = cj_ind(j)
  !               udelc = -r(j)
                  udelc = -r(j) + log(ucc)
                  if (abs(udelc) .gt. abs(tolexp)) then
                    if (iwarn.ge.2) write (*,100) mcyc,l1,l2,l3,iter,icut, &
                    nam(j),cj_ind(j),udelc
    100             format ('OSRXN: delc > tolexp, mcyc =',i5,' n = ',2i4, &
      &            ' iter =',i3,' icuts =',i2,1x,a8, &
      &            '=',1pe12.4,' delc =',1pe12.4)
                    ! WRITE(*, '(A, 1X, I0)') 'checkpoint2: step =',istep
                    ! WRITE(*, '(A, 1X, I0)') 'checkpoint2: myrank =',myrank
                    if (tolexp.lt.zero) then
  !                   goto 70  ! cut time step
                    else
                      udelc = tolexp*udelc/abs(udelc)
                    endif
                  endif
                  cj_ind(j) = exp(udelc)
                  udelc = cj_ind(j)-ucc
                  if(abs(udelc).gt.abs(delcmx)) delcmx = udelc
                enddo

              else                           ! solution is delc

                do j = 1, ncomp 
                  udelc = -r(j)
                  if(abs(udelc).gt.abs(delcmx)) delcmx = udelc !误差取最大值
                  ucc = cj_ind(j)
                  dcc = max (one, -two*udelc/ucc)
                  cj_ind(j) = ucc+udelc/dcc
                  if(cj_ind(j).le.0.d0) then ! check for negative conc
                    cj_ind(j) = ucc
                    if (iwarn.ge.2) &
                    write(*,'("mcyc= ",i4," iter= ",i3, &
      &            " OSRXN:--negative or zero concentration for species: ", &
      &            a8," at node: ",2i3,1p10e12.4)') &
                    mcyc,iter,nam(j),l1,l2,l3,udelc,ucc
                    ! WRITE(*, '(A, 1X, I0)') 'checkpoint3: step =',istep
                    ! WRITE(*, '(A, 1X, I0)') 'checkpoint3: myrank =',myrank
  !                 goto 70
                  endif
                enddo
              endif


              if (iter < itmax) goto 10
              print *,'Max iters at node: (',myrank,', ',l1,', ',l2,', ',l3,')',' its = ', &
              iter,' mcyc = ',istep



              print *,'Res = ',(r(j),j=1,ncomp)
              print *,'cj = ',(cj_ind(j),j=1,ncomp)
              print *,'psi = ',(psi_ind(j),j=1,ncomp)
              do j = 1, ncomp
                print *,'dpsi = ',(dpsi(j,k),k=1,ncomp)
              enddo
              ! WRITE(*, '(A, 1X, I0)') 'checkpoint4: step =',istep
              ! WRITE(*, '(A, 1X, I0)') 'checkpoint4: myrank =',myrank
              CALL MPI_ABORT(MPI_COMM_WORLD, 0, ierr)
              
    30       continue
              if (iter > itermx) itermx = iter

\


      end subroutine con_cj_ind


!       subroutine con_cj_ind(psi,cj,istep,l1,l2,l3)!,lbuffers)

!       use ptran_global_module
!       use solve_module
!       use perind_module

!       use comvar_module
!       use mpi

!       implicit none

!       real*8, dimension(1:m,mx,my+2,mz+2) :: cj,psi


!       real*8 :: alogjn(ncomp),r(ncomp),dpsi(ncomp,ncomp),cx(ncmplx)

!       integer :: l1,l2,l3,i,iter,itmax,j,k,indx(ncomp),itermx,istep,mcyc
!       real*8 :: delcmx,prod1,sumc,udelc,ucc,dd,dcc
!       integer,dimension(3) :: pd_vct_p,pd_vct_n
!       integer :: i_p,j_p,k_p,i_n,j_n,k_n


!       mcyc = istep
!       if (ncmplx == 0) then
!         cj=psi
!         return
!       endif

!       tol  = 1.d-10
!       itmax = 105
!       icut = 0
! !     loglin = 0 ! solve for log C
!       loglin = 1 ! solve for  C
!       iwarn = 2
!       tolexp = 10.d0
!       itermx = 0

! !     begin loop over nodes


!               iter = 0

!     10       continue
!               iter = iter + 1

!               do j = 1, ncomp
!   !             alogjn(j) = log(cj(j,l1,l2)*gam(j,n))
!                 alogjn(j) = log(cj(j,l1,l2,l3))
!               enddo

!               do i = 1, ncmplx
!   !             prod1 = eqhom(i)*aln10-log(gamx(i,n))
!                 prod1 = eqhom(i)*aln10
!                 do k = 1, ncomp
!                   prod1 = prod1 + shom(k,i)*alogjn(k)
!                 enddo
!                 cx(i) = exp(prod1)
!               enddo

!   !-----------compute psi for primary species
!               do j = 1, ncomp
!                 sumc = cj(j,l1,l2,l3)
!                 do i = 1, ncmplx
!                   sumc = sumc + shom(j,i)*cx(i)
!                 enddo
!   !             sumc = sumc*rho(n)
!                 r(j) = sumc - psi(j,l1,l2,l3)
!   !             print *,'sum-r: ',iter,j,r(j),sumc,psi(j,l1,l2)

!   !             diagonal terms
!                 sumc = cj(j,l1,l2,l3)
!   !             print *,'sum: ',j,sum,r(j)
!                 do i = 1, ncmplx
!                   sumc = sumc + shom(j,i)*shom(j,i)*cx(i)
!   !               print *,'sum: ',i,cx(i),shom(j,i),sumc
!                 enddo
!   !             dpsi(j,j) = sumc*rho(n)
!                 dpsi(j,j) = sumc

!   !             print *,'diag: ',j,dpsi(j,j),cj(j,l1,l2)

!   !             off-diagonal terms
!                 do k = j+1, ncomp
!                   sumc = 0.d0
!                   do i = 1, ncmplx
!                     sumc = sumc + shom(j,i)*shom(k,i)*cx(i)
!                   enddo
!   !               dpsi(j,k) = sumc*rho(n)
!                   dpsi(j,k) = sumc
!                   dpsi(k,j) = sumc
!                 enddo
!               enddo

!               if (loglin == 1) then
!                 do j = 1, ncomp
!                   do k = 1, ncomp
!                     dpsi(j,k) = dpsi(j,k)/cj(k,l1,l2,l3)
!                   enddo
!                 enddo
!               endif

!   !           check convergence
!               if (iter > 1) then
!                 do j = 1, ncomp
!                   if (abs(r(j)) > tol) goto 20
!                 enddo

!   !             stop

!                 goto 30
!               endif

!     20       continue

!               if (mcyc==-2) then
!                 print *,'res0: ',mcyc,iter,i1,i2,(r(j),j=1,ncomp)
!                 print *,'cj: ',mcyc,iter,i1,i2,(cj(j,l1,l2,l3),j=1,ncomp)
!                 print *,'psi: ',mcyc,iter,(psi(j,l1,l2,l3),j=1,ncomp)
!                 do j = 1, ncomp
!                   print *,'diag: ',(dpsi(j,k),k=1,ncomp)
!                 enddo
!                 ! WRITE(*, '(A, 1X, I0)') 'checkpoint1: step =',istep
!                 ! WRITE(*, '(A, 1X, I0)') 'checkpoint1: myrank =',myrank
!               endif

!               call solve (ncomp,r,dpsi)
!   !           call ludcmp(dpsi,ncomp,ncomp,indx,dd)
!   !           call lubksb(dpsi,ncomp,ncomp,indx,r)

!   !           print *,'res: ',mcyc,iter,i1,i2,(r(j),j=1,ncomp)

!   !-----------update field variables corresponding to new iteration
!               delcmx = 0.d0
!               if (loglin .eq. 0) then           ! solution is log(delc) 
!                 do j = 1, ncomp
!                   ucc = cj(j,l1,l2,l3)
!   !               udelc = -r(j)
!                   udelc = -r(j) + log(ucc)
!                   if (abs(udelc) .gt. abs(tolexp)) then
!                     if (iwarn.ge.2) write (*,100) mcyc,l1,l2,l3,iter,icut, &
!                     nam(j),cj(j,l1,l2,l3),udelc
!     100             format ('OSRXN: delc > tolexp, mcyc =',i5,' n = ',2i4, &
!       &            ' iter =',i3,' icuts =',i2,1x,a8, &
!       &            '=',1pe12.4,' delc =',1pe12.4)
!                     ! WRITE(*, '(A, 1X, I0)') 'checkpoint2: step =',istep
!                     ! WRITE(*, '(A, 1X, I0)') 'checkpoint2: myrank =',myrank
!                     if (tolexp.lt.zero) then
!   !                   goto 70  ! cut time step
!                     else
!                       udelc = tolexp*udelc/abs(udelc)
!                     endif
!                   endif
!   !               ucc = cj(j,l1,l2)
!   !               udelc = udelc + log(ucc)
!   !               cj(j,l1,l2) = ucc*exp(udelc)
!                   cj(j,l1,l2,l3) = exp(udelc)
!                   udelc = cj(j,l1,l2,l3)-ucc
!                   if(abs(udelc).gt.abs(delcmx)) delcmx = udelc
!                 enddo

!               else                           ! solution is delc

!                 do j = 1, ncomp
!                   udelc = -r(j)
!                   if(abs(udelc).gt.abs(delcmx)) delcmx = udelc !误差取最大值
!                   ucc = cj(j,l1,l2,l3)
!   !               cj(j,l1,l2) = ucc+udelc
!                   dcc = max (one, -two*udelc/ucc)
!                   cj(j,l1,l2,l3) = ucc+udelc/dcc
!                   if(cj(j,l1,l2,l3).le.0.d0) then ! check for negative conc
!                     cj(j,l1,l2,l3) = ucc
!                     if (iwarn.ge.2) &
!                     write(*,'("mcyc= ",i4," iter= ",i3, &
!       &            " OSRXN:--negative or zero concentration for species: ", &
!       &            a8," at node: ",2i3,1p10e12.4)') &
!                     mcyc,iter,nam(j),l1,l2,l3,udelc,ucc
!                     ! WRITE(*, '(A, 1X, I0)') 'checkpoint3: step =',istep
!                     ! WRITE(*, '(A, 1X, I0)') 'checkpoint3: myrank =',myrank
!   !                 goto 70
!                   endif
!                 enddo
!               endif

!   !           print *,'cj1:   ',mcyc,iter,(cj(j,l1,l2),j=1,ncomp)
!   !           print *,'psi1: ',mcyc,iter,(psi(j,l1,l2),j=1,ncomp)

!               if (iter < itmax) goto 10
!               print *,'Max iters at node: (',myrank,', ',l1,', ',l2,', ',l3,')',' its = ', &
!               iter,' mcyc = ',istep



!               print *,'Res = ',(r(j),j=1,ncomp)
!               print *,'cj = ',(cj(j,l1,l2,l3),j=1,ncomp)
!               print *,'psi = ',(psi(j,l1,l2,l3),j=1,ncomp)
!               do j = 1, ncomp
!                 print *,'dpsi = ',(dpsi(j,k),k=1,ncomp)
!               enddo
!               ! WRITE(*, '(A, 1X, I0)') 'checkpoint4: step =',istep
!               ! WRITE(*, '(A, 1X, I0)') 'checkpoint4: myrank =',myrank
!               CALL MPI_ABORT(MPI_COMM_WORLD, 0, ierr)
!     30       continue
!               if (iter > itermx) itermx = iter
!   !           write(*,'(2i3,(a8,1pe12.4))') l1,l2,(nam(j),r(j),j=1,ncomp)
! \


!       end subroutine con_cj_ind
!--------------------------

      subroutine con_ci(ci,cj)

      use ptran_global_module
      
      use comvar_module
      
      implicit none
      
      integer :: i,j

      real*8,dimension(1:ncomp,mx,my+2,mz+2):: cj
      real*8,dimension(1:ncmplx,mx,my+2,mz+2):: ci  

      do i=1,ncmplx
        ci(i,:,:,:)=aln10*eqhom(i) !-log(gama_i(i))
        do j=1,ncomp
      	  ci(i,:,:,:)=ci(i,:,:,:)+log(gama_j(j)*cj(j,:,:,:))*shom(j,i)
        enddo
        ci(i,:,:,:)=exp(ci(i,:,:,:))
      enddo	      
      
      end subroutine con_ci
!-----------------------------------------------------
     subroutine bnd_cj(gi_tmp,psi_tmp,cj_tmp,l1,l2,l3,wall_ind,Rm,istep)

      use ptran_global_module
      use solve_module

      use comvar_module
      use mpi
      
      implicit none

      real*8, dimension(1:ncomp) :: gi_tmp,psi_tmp,cj_tmp
      
      real*8 :: alogjn(ncomp),r(ncomp),dpsi(ncomp,ncomp),cx(ncmplx), &
                drate(ncomp,ncomp)
      real*8, dimension(1:ncomp) :: Rmm
      real*8, dimension(nkin,1:mx,1:my+2,1:mz+2) :: Rm

      
      real*8 :: qk(nkin)

      integer :: npri,nsec,lp,kk,jj,istep
      real*8 :: prefac,cxin
      
      integer :: l1,l2,l3,i,iter,itmax,j,k,indx(ncomp),itermx
      real*8 :: delcmx,prod1,sumc,sumi,udelc,ucc,dd,dcc
      
      logical, dimension(nkin) :: wall_ind 
      
      tol  = 1.d-8
      itmax = 500
      icut = 0
!     loglin = 0 ! solve for log C
      loglin = 1 ! solve for  C
      iwarn = 2
      tolexp = 10.d0
      itermx = 0

      iter = 0
          
   10 continue

      iter = iter + 1
      

            do j = 1, ncomp
!             alogjn(j) = log(cj(j,l1,l2)*gam(j,n))
              alogjn(j) = log(cj_tmp(j))
              ! print*,alogjn(j),cj(j,l1,l2,l3),'0'
            enddo

            do i = 1, ncmplx
!             prod1 = eqhom(i)*aln10-log(gamx(i,n))
              prod1 = eqhom(i)*aln10

              do k = 1, ncomp
                prod1 = prod1 + shom(k,i)*alogjn(k)
              enddo
              cx(i) = exp(prod1)
            enddo
            
!           compute mineral reaction rates
            ! print*, iter,l1,l2,l3,qk(i),'1',cj(j,l1,l2,l3)
            do i=1,nkin
              qk(i) = aln10*eqkin(i)
              
              ! print*, iter,l1,l2,l3,eqkin(i),'2'
              do j=1,ncomp
                qk(i)=qk(i)+alogjn(j)*skin(j,i)
                ! if (l1==1) then
                !   print*, iter,l1,l2,l3,qk(i),'3'
                ! endif
                ! print*,alogjn(j),skin(j,i)
              enddo
              qk(i) = exp(qk(i))


!*** prefactor start

!               kb(i) = 0.d0

! !-----------compute rate for each parallel reaction and add contribution
! !           to jacobian
!               do lp = npar1(i), npar2(i)

!   !-------------compute prefactor
!                 npri = nkinpri(lp)
!                 nsec = nkinsec(lp)
!                 prefac = one
!                 if (npri.gt.0 .or. nsec.gt.0) then
!                   prefac = zero
!                   if (npri .gt. 0) then
!                     do j = 1, npri
!                       jj = jpri(j,lp)
!                       prefac = prefac + skinpri(j,lp)*alogjn(jj)
!   !                   write(*,*) 'trkinmin: ',nam(j),prefac,jpri(j,lp),
!   !    .              skinpri(j,lp),ajnlog(jpri(j,lp))
!                     enddo
!                   endif
!                   if (nsec .gt. 0) then
!                     do k = 1, nsec
!                       kk = isec(k,lp)
!                       cxin = log(cx(kk))
!                       prefac = prefac + skinsec(k,lp)*cxin
!                     enddo
!                   endif
!                   prefac = exp(prefac)
!                 else
!                   prefac = 1.d0
!                 endif
!                 kb(i) = kb(i) + rkf(lp)*prefac
!               enddo
              
              
              ! kb(i)=kb(i)*dy0(1)*(dif_lb/difaq)

!*** prefactor end
              
              Rm(i,l1,l2,l3)=-kb(i)*(1.d0-qk(i))
              
              ! print*,Rm(i,l1,l2,l3),kb(i),qk(i),myrank

              if(.not.wall_ind(i).and.Rm(i,l1,l2,l3).lt.0.d0) Rm(i,l1,l2,l3)=0.d0

              ! if (l1==1) then
              !   print*,Rm(i,l1,l2,l3)
              ! endif
            enddo    
            ! if((l1==10).and.(l2==1).and.(l3==1)) then
            !   print*,kb(1),eqkin(1),skin(1,1),alogjn(1),aln10,qk(1)
            ! endif
          
            Rmm=0.d0
            do j=1,ncomp
              do i=1,nkin
              	Rmm(j)=Rmm(j)+skin(j,i)*Rm(i,l1,l2,l3)
              enddo
            enddo  

!-----------compute psi for primary species
            do j = 1, ncomp 
              sumc = cj_tmp(j)
              do i = 1, ncmplx
                sumc = sumc + shom(j,i)*cx(i)
              enddo
!             sumc = sumc*rho(n)
              
              psi_tmp(j) = sumc
              
              r(j) = sumc/8.0 + Rmm(j) - gi_tmp(j) 
              !print*,sumc,Rmm(j),gi_tmp(j) 
        
!             compute derivatives: dpsi_j/dlnc_k
!             diagonal terms
              sumc = cj_tmp(j)
              do i = 1, ncmplx
                sumc = sumc + shom(j,i)*shom(j,i)*cx(i)
              enddo
!             dpsi(j,j) = sumc*rho(n)
              dpsi(j,j) = sumc

!             off-diagonal terms
              do k = j+1, ncomp
                sumc = 0.d0
                do i = 1, ncmplx
                  sumc = sumc + shom(j,i)*shom(k,i)*cx(i)
                enddo
!               dpsi(j,k) = sumc*rho(n)
                dpsi(j,k) = sumc
                dpsi(k,j) = sumc
              enddo

!             rate term
!             diagonal
              sumi = 0.d0
              do i = 1, nkin
                sumi = sumi + skin(j,i)*skin(j,i)*kb(i)*qk(i)
              enddo
              drate(j,j) = sumi
              
              !off diagonal
              do k = j+1, ncomp
                sumi = 0.d0
                do i = 1, nkin
                  sumi = sumi + skin(j,i)*skin(k,i)*kb(i)*qk(i)
                enddo
                drate(j,k) = sumi
                drate(k,j) = sumi
              enddo
            enddo
            
            do j = 1, ncomp
              do k = 1, ncomp
!               print *,'conc: ',iter,j,k,omega,c1(j),dpsi(j,k),drate(j,k)
                
                dpsi(j,k) = dpsi(j,k)/8.0 + drate(j,k)
!               print *,'dpsi: ',iter,j,k,dpsi(j,k)
              enddo
			      enddo
            if (loglin == 1) then
              do j = 1, ncomp
                do k = 1, ncomp
                  dpsi(j,k) = dpsi(j,k)/cj_tmp(k)
                enddo
              enddo
            endif
          
!           check convergence
            if (iter > 1) then
              do j = 1, ncomp
                if (abs(r(j)) > tol) goto 20
              enddo			  
              goto 30
            endif
          
   20       continue
            
            if (iter==-1) then
              print *,'res0: ',iter,i1,i2,(r(j),j=1,ncomp)
              print *,'cj: ',iter,i1,i2,(cj_tmp(j),j=1,ncomp)
              print *,'psi: ',iter,(psi_tmp(j),j=1,ncomp)
              do j = 1, ncomp
                print *,'diag: ',(dpsi(j,k),k=1,ncomp)
              enddo
            endif
            
            call solve (ncomp,r,dpsi)
!           call ludcmp(dpsi,ncomp,ncomp,indx,dd)
!           call lubksb(dpsi,ncomp,ncomp,indx,r)
            
!           print *,'res: ',iter,i1,i2,(r(j),j=1,ncomp)

!-----------update field variables corresponding to new iteration
            delcmx = 0.d0
            if (loglin .eq. 0) then           ! solution is log(delc) 
              do j = 1, ncomp
                ucc = cj_tmp(j)
!               udelc = -r(j)
                udelc = -r(j) + log(ucc)
                if (abs(udelc) .gt. abs(tolexp)) then
                  if (iwarn.ge.2) write (*,100) l1,l2,l3,iter,icut, &
                  nam(j),cj_tmp(j),udelc
  100             format ('OSRXN: delc > tolexp, n = ',2i4, &
     &            ' iter =',i3,' icuts =',i2,1x,a8, &
     &            '=',1pe12.4,' delc =',1pe12.4)
                  if (tolexp.lt.zero) then
!                   goto 70  ! cut time step
                  else
                    udelc = tolexp*udelc/abs(udelc)
                  endif
                endif
!               ucc = cj(j,l1,l2,l3)
!               udelc = udelc + log(ucc)
!               cj(j,l1,l2,l3) = ucc*exp(udelc)
                cj_tmp(j) = exp(udelc)
                udelc = cj_tmp(j)-ucc
                if(abs(udelc).gt.abs(delcmx)) delcmx = udelc
              enddo

            else                           ! solution is delc

              do j = 1, ncomp 
                udelc = -r(j)
                if(abs(udelc).gt.abs(delcmx)) delcmx = udelc
                ucc = cj_tmp(j)
!               cj(j,l1,l2) = ucc+udelc
                dcc = max (one, -two*udelc/ucc)
                cj_tmp(j) = ucc+udelc/dcc
                if(cj_tmp(j).le.0.d0) then ! check for negative conc
                  cj_tmp(j) = ucc
                  if (iwarn.ge.2) &
                  write(*,'(" iter= ",i3, &
     &            " OSRXN:--negative or zero concentration for species: ", &
     &            a8," at node: ",2i3,1p10e12.4)') &
                  iter,nam(j),l1,l2,l3,udelc,ucc
!                 goto 70
                endif
              enddo
            endif

!           print *,'cj1:   ',iter,(cj(j,l1,l2,l3),j=1,ncomp)
!           print *,'psi1: ',iter,(psi(j,l1,l2,l3),j=1,ncomp)

            if (iter < itmax) goto 10

            print *,'Max iters at node: (',l1,', ',l2,', ',l3,')',' its = ', &
            iter
            print *,'Res = ',tol,(r(j),j=1,ncomp)
            print *,'cj = ',(cj_tmp(j),j=1,ncomp)
            print *,'psi = ',(psi_tmp(j),j=1,ncomp)
            do j = 1, ncomp
              print *,'dpsi = ',(dpsi(j,k),k=1,ncomp)
            enddo
            ! WRITE(*, '(A, 1X, I0)') 'checkpoint5: step =',istep
            ! WRITE(*, '(A, 1X, I0)') 'checkpoint5: myrank =',myrank
            CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)

   30       continue

            if (iter > itermx) itermx = iter
!           write(*,'(2i3,(a8,1pe12.4))') l1,l2,l3,(nam(j),r(j),j=1,ncomp)
      
      end subroutine bnd_cj


end module con_module


