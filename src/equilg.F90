module equilg_module

contains

    subroutine equilg(geq,psie,uxe,uye,uze)

      use comvar_module

      implicit none

      integer :: s
      real*8 :: uxe,uye,uze,psie
      real*8 :: tmpg
      real*8, dimension(0:6) :: geq

      geq(0)=wts_sol(0)*psie

      do s = 1,6
        tmpg = cjx(s)*uxe + cjy(s)*uye + cjz(s)*uze
        geq(s) = wts_sol(s)*psie*(1.0d0 + 4.0d0*tmpg)

      enddo
        
    end subroutine equilg

end module equilg_module