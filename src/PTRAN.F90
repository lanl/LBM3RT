!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! VERSION/REVISION HISTORY
 
! $Id: ptran.F90,v 1.1.1.1 2004/07/30 21:49:42 lichtner Exp $
! $Log: ptran.F90,v $
! Revision 1.1.1.1  2004/07/30 21:49:42  lichtner
! initial import
!
! Revision 1.3  2004/04/06 17:39:55  lichtner
! Added source term time stepping check.
!
! Revision 1.2  2004/01/10 18:32:06  lichtner
! Began work on 2 phase capability.
!
! Revision 1.1.1.1  2003/11/23 20:12:46  lichtner
! initial entry
!

!  pFLOTRAN Version 1.0 LANL
!-----------------------------------------------------------------------
!  Date             Author(s)                Comments/Modifications
!-----------------------------------------------------------------------
!  May  2003        Peter C. Lichtner        Initial Implementation
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  program ptran

  use ptran_global_module
  use trdynmem_module
  use ptran_read_module
  use ptran_init_module
   
  implicit none 

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 Beginning of Program: PFLOTRAN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  call ptran_read

  call ptran_chem

  write(*,*) '--> initialization complete'
  
  end program ptran
