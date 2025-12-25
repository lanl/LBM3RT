module comvar_module

  implicit none

      integer :: lx,ly,lz
      real*8, dimension(0:18) :: cix,ciy,ciz,wts
      integer, dimension(0:18) :: cixi,ciyi,cizi,ind_ci_opp
      real*8, dimension(0:6) :: cjx,cjy,cjz,wts_sol
      integer, dimension(0:6) :: cjxi,cjyi,cjzi,ind_cj_opp
      logical :: dpst,dpst_up,dpst_down,dpst_front,dpst_rear
      logical :: dpst_front_down,dpst_front_up,dpst_rear_down,dpst_rear_up
      real*8 :: dif_lb,l_phys,damk,chance,klb,tau,ffx=0.0d0
      
      real*8, allocatable, dimension(:)     :: tau_p,cjini,psiini,cjin,psiin, &
                                               c1,gama_j,gama_i,npos,alog
      real*8, allocatable, dimension(:)     :: ciin,ciini
      real*8, allocatable, dimension(:)     :: kb

      ! MPI Parameters
      integer:: ierr
      integer:: xyplane_db1,xzplane_db1,xyplane_db2,xzplane_db2,xyplane_db2x,xzplane_db2x
      integer:: xyplane_db3,xzplane_db3,xyplane_db5,xzplane_db5
      integer:: xyplane_db6,xzplane_db6,xyplane_db7,xzplane_db7
      integer:: xyplane_int
      integer:: xyplane_logical,xzplane_logical
      integer, parameter::  mpid=3 ! proc grid dimensions
      integer:: front(2),rear(2),right(2)
      integer:: up(2),down(2),left(2)
      integer ::front_down,front_up,rear_down,rear_up
      logical remdims(mpid)
      logical Tperiodic(mpid) !periodicity of procs in each dimension
      logical rreorder ! whether reorder processor rank to optimize
      integer:: lbecomm,localcomm
      integer:: mydev 
      integer:: prgrid(mpid) !number of processes on each dimension
      integer:: mpicoords(mpid)
      integer:: gsizes(3),lsizes(3),start_idx(3)
      integer:: offset(3)
      integer:: blocklength, stride
      integer:: mx,my,mz,zrem,yrem
      integer:: istat

      integer, ALLOCATABLE :: ny_local(:), nz_local(:) !end-point position
      integer, ALLOCATABLE :: my_local(:), mz_local(:) !mpi grid width
      
  contains

  subroutine comvar()


  implicit none

      cix(0) = 0.d0
      cix(1) = 1.d0
      cix(2) = -1.d0
      cix(3) = 0.d0
      cix(4) = 0.d0
      cix(5) = 0.d0
      cix(6) = 0.d0
      cix(7) = 1.d0
      cix(8) = -1.d0
      cix(9) = 1.d0
      cix(10) = -1.d0
      cix(11) = 1.d0
      cix(12) = -1.d0
      cix(13) = 1.d0
      cix(14) = -1.d0
      cix(15) = 0.d0
      cix(16) = 0.d0
      cix(17) = 0.d0
      cix(18) = 0.d0

      ciy(0) = 0.d0
      ciy(1) = 0.d0
      ciy(2) = 0.d0
      ciy(3) = 1.d0
      ciy(4) = -1.d0
      ciy(5) = 0.d0
      ciy(6) = 0.d0
      ciy(7) = 1.d0
      ciy(8) = 1.d0
      ciy(9) = -1.d0
      ciy(10) = -1.d0
      ciy(11) = 0.d0
      ciy(12) = 0.d0
      ciy(13) = 0.d0
      ciy(14) = 0.d0
      ciy(15) = 1.d0
      ciy(16) = -1.d0
      ciy(17) = 1.d0
      ciy(18) = -1.d0

      ciz(0) = 0.d0
      ciz(1) = 0.d0
      ciz(2) = 0.d0
      ciz(3) = 0.d0
      ciz(4) = 0.d0
      ciz(5) = 1.d0
      ciz(6) = -1.d0
      ciz(7) = 0.d0
      ciz(8) = 0.d0
      ciz(9) = 0.d0
      ciz(10) = 0.d0
      ciz(11) = 1.d0
      ciz(12) = 1.d0
      ciz(13) = -1.d0
      ciz(14) = -1.d0
      ciz(15) = 1.d0
      ciz(16) = 1.d0
      ciz(17) = -1.d0
      ciz(18) = -1.d0

      wts(0) = 1.0d0/3.0d0
      wts(1) = 1.0d0/18.0d0
      wts(2) = 1.0d0/18.0d0
      wts(3) = 1.0d0/18.0d0
      wts(4) = 1.0d0/18.0d0
      wts(5) = 1.0d0/18.0d0
      wts(6) = 1.0d0/18.0d0
      wts(7) = 1.0d0/36.0d0
      wts(8) = 1.0d0/36.0d0
      wts(9) = 1.0d0/36.0d0
      wts(10) = 1.0d0/36.0d0
      wts(11) = 1.0d0/36.0d0
      wts(12) = 1.0d0/36.0d0
      wts(13) = 1.0d0/36.0d0
      wts(14) = 1.0d0/36.0d0
      wts(15) = 1.0d0/36.0d0
      wts(16) = 1.0d0/36.0d0
      wts(17) = 1.0d0/36.0d0
      wts(18) = 1.0d0/36.0d0

      cixi(0) = 0
      cixi(1) = 1
      cixi(2) = -1
      cixi(3) = 0
      cixi(4) = 0
      cixi(5) = 0
      cixi(6) = 0
      cixi(7) = 1
      cixi(8) = -1
      cixi(9) = 1
      cixi(10) = -1
      cixi(11) = 1
      cixi(12) = -1
      cixi(13) = 1
      cixi(14) = -1
      cixi(15) = 0
      cixi(16) = 0
      cixi(17) = 0
      cixi(18) = 0

      ciyi(0) = 0
      ciyi(1) = 0
      ciyi(2) = 0
      ciyi(3) = 1
      ciyi(4) = -1
      ciyi(5) = 0
      ciyi(6) = 0
      ciyi(7) = 1
      ciyi(8) = 1
      ciyi(9) = -1
      ciyi(10) = -1
      ciyi(11) = 0
      ciyi(12) = 0
      ciyi(13) = 0
      ciyi(14) = 0
      ciyi(15) = 1
      ciyi(16) = -1
      ciyi(17) = 1
      ciyi(18) = -1

      cizi(0) = 0
      cizi(1) = 0
      cizi(2) = 0
      cizi(3) = 0
      cizi(4) = 0
      cizi(5) = 1
      cizi(6) = -1
      cizi(7) = 0
      cizi(8) = 0
      cizi(9) = 0
      cizi(10) = 0
      cizi(11) = 1
      cizi(12) = 1
      cizi(13) = -1
      cizi(14) = -1
      cizi(15) = 1
      cizi(16) = 1
      cizi(17) = -1
      cizi(18) = -1    

      ind_ci_opp(0) = 0
      ind_ci_opp(1) = 2
      ind_ci_opp(2) = 1
      ind_ci_opp(3) = 4
      ind_ci_opp(4) = 3
      ind_ci_opp(5) = 6
      ind_ci_opp(6) = 5
      ind_ci_opp(7) = 10
      ind_ci_opp(8) = 9
      ind_ci_opp(9) = 8
      ind_ci_opp(10) = 7
      ind_ci_opp(11) = 14
      ind_ci_opp(12) = 13
      ind_ci_opp(13) = 12
      ind_ci_opp(14) = 11
      ind_ci_opp(15) = 18
      ind_ci_opp(16) = 17
      ind_ci_opp(17) = 16
      ind_ci_opp(18) = 15

      cjx(0) = 0.d0
      cjx(1) = 1.d0
      cjx(2) = -1.d0
      cjx(3) = 0.d0
      cjx(4) = 0.d0
      cjx(5) = 0.d0
      cjx(6) = 0.d0

      cjy(0) = 0.d0
      cjy(1) = 0.d0
      cjy(2) = 0.d0
      cjy(3) = 1.d0
      cjy(4) = -1.d0
      cjy(5) = 0.d0
      cjy(6) = 0.d0

      cjz(0) = 0.d0
      cjz(1) = 0.d0
      cjz(2) = 0.d0
      cjz(3) = 0.d0
      cjz(4) = 0.d0
      cjz(5) = 1.d0
      cjz(6) = -1.d0

      cjxi(0) = 0
      cjxi(1) = 1
      cjxi(2) = -1
      cjxi(3) = 0
      cjxi(4) = 0
      cjxi(5) = 0
      cjxi(6) = 0

      cjyi(0) = 0
      cjyi(1) = 0
      cjyi(2) = 0
      cjyi(3) = 1
      cjyi(4) = -1
      cjyi(5) = 0
      cjyi(6) = 0

      cjzi(0) = 0
      cjzi(1) = 0
      cjzi(2) = 0
      cjzi(3) = 0
      cjzi(4) = 0
      cjzi(5) = 1
      cjzi(6) = -1  

      ind_cj_opp(0) = 0
      ind_cj_opp(1) = 2
      ind_cj_opp(2) = 1
      ind_cj_opp(3) = 4
      ind_cj_opp(4) = 3
      ind_cj_opp(5) = 6
      ind_cj_opp(6) = 5

      wts_sol(0) = 0.25d0
      wts_sol(1) = 0.125d0
      wts_sol(2) = 0.125d0
      wts_sol(3) = 0.125d0
      wts_sol(4) = 0.125d0
      wts_sol(5) = 0.125d0
      wts_sol(6) = 0.125d0

  end subroutine comvar

end module comvar_module
