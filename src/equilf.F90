module equilf_module

contains

    subroutine equilf(feq,rhoe,uxe,uye,uze)

      use comvar_module

      implicit none

      integer :: s
      real*8 :: uxe,uye,uze,rhoe
      real*8 :: usqr,pres
      real*8, dimension(0:18) :: feq,tmp

      usqr = uxe*uxe + uye*uye + uze*uze
      

      feq(0)=wts(0)*(rhoe-1.5d0*usqr)

      do s = 1,18
        tmp(s) = cix(s)*uxe + ciy(s)*uye + ciz(s)*uze
        feq(s) = wts(s)*(rhoe + 3.0d0*tmp(s) + 4.5d0*tmp(s)*tmp(s) - 1.5d0*usqr)

      enddo
        
    end subroutine equilf

end module equilf_module
