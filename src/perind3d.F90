module perind_module

  contains

  function perind(i,j,k,dx,dy,dz,xe,ye,ze) result(pd_vct) !periodic index for fluids

    implicit none

    integer :: i,j,k,dx,dy,dz,xe,ye,ze
    integer, dimension(3) :: pd_vct 

    pd_vct(1) = i + dx
    pd_vct(2) = j + dy
    pd_vct(3) = k + dz

    if (pd_vct(1)<1) then
      pd_vct(1) = xe
    endif

    if (pd_vct(1)>xe) then
      pd_vct(1) = 1
    endif

    if (pd_vct(2)<1) then
      pd_vct(2) = ye+2
    endif
    
    if (pd_vct(2)>ye+2) then
      pd_vct(2) = 1
    endif

    if (pd_vct(3)<1) then
      pd_vct(3) = ze+2
    endif
    
    if (pd_vct(3)>ze+2) then
      pd_vct(3) = 1
    endif

  end function perind

  function perind_s(i,j,k,dx,dy,dz,xe,ye,ze) result(pd_vct) !periodic index for solids

    implicit none

    integer :: i,j,k,dx,dy,dz,xe,ye,ze
    integer, dimension(3) :: pd_vct 

    pd_vct(1) = i + dx
    pd_vct(2) = j + dy
    pd_vct(3) = k + dz

    if (pd_vct(1)<1) then
      pd_vct(1) = xe
    endif

    if (pd_vct(1)>xe) then
      pd_vct(1) = 1
    endif

    if (pd_vct(2)<0) then
      pd_vct(2) = ye+3
    endif
    
    if (pd_vct(2)>ye+3) then
      pd_vct(2) = 0
    endif

    if (pd_vct(3)<0) then
      pd_vct(3) = ze+3
    endif
    
    if (pd_vct(3)>ze+3) then
      pd_vct(3) = 0
    endif

  end function perind_s

end module perind_module
