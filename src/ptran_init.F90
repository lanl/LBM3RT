!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! VERSION/REVISION HISTORY
 
! $Id: ptran_init.F90,v 1.1.1.1 2004/07/30 21:49:42 lichtner Exp $
! $Log: ptran_init.F90,v $
! Revision 1.1.1.1  2004/07/30 21:49:42  lichtner
! initial import
!
! Revision 1.3  2004/04/06 17:33:01  lichtner
! Revised storage of phik and surf. Set rho=1 temporarily.
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

module ptran_init_module

  public

contains
  
  subroutine ptran_chem

  use ptran_global_module
  use trdynmem_module
  use ptran_speciation_module
  use trgamdh_module
  use ptran_setbnd_module

  implicit none 
  integer :: ierr,i,j,l,m
! integer :: k,ir,ii1,ii2,jj1,jj2,kk1,kk2,n,ng,nr
  integer :: isumjl
  real*8, pointer :: ccloc_p(:),temploc_p(:),pressloc_p(:),ssat_loc_p(:)

  if (myrank == 0) write(*,*) '--> initialize field variables'
  call trinit

  if (myrank == 0) write(*,*) '--> speciate initial fluid composition'
  call trstartup

  if (myrank==0) write(*,*) '--> speciate boundary/source fluid composition'
  call trsetbnd

  end subroutine ptran_chem

!===================================================================

  subroutine trinit
  
  use ptran_global_module
  use trdynmem_module
  use water_eos_module
  
  implicit none
  
  integer :: i,ii,j,iireg,jj,jaq,k,kk,kcount,kcountp,l,lp,m,mm,ml, &
             ns,ntemp,ncsq,ncsqdp,ncsqdpg,nr
  
  real*8 :: den,atm,btm,bdtm,frc,dadebye,dbdebye,dbext,bextend0,fjl

!-----set activity coefficient algorithm
      ntemp = ntmpmx
      i = 1
      do ii=1,ntemp
      
!       write(*,*) 'trinit: ',ii,tempini,temptab(ii)
        
        if(tempini.le.temptab(ii)) then
          i = ii
          goto 10
        endif
      enddo
      if (myrank == 0) write(iunit2,1001) tempini
 1001 format(' illegal temperature (',1pg12.4,') deg.c')
      if(tempini.le.300.) then
        if (myrank == 0) &
        write(*,*) 'error in temperature interpolation! STOP'
        stop
      endif

   10 continue

      if (i .eq. 1) then
        adebye  = atab(i)
        bdebye  = btab(i)
        bextend0 = bdotab(i)
      else
        atm = atab(i-1)
        btm = btab(i-1)
        bdtm = bdotab(i-1)
        frc   = (tempini-temptab(i-1))/(temptab(i)-temptab(i-1))
        dadebye = atab(i)   - atm
        dbdebye = btab(i)   - btm
        dbext   = bdotab(i) - bdtm
        adebye  = frc*dadebye + atm
        bdebye  = frc*dbdebye + btm
        bextend0 = frc*dbext  + bdtm
      endif
      if (iact.eq.1) then
        do j = 1, ncomp
          bextend(j) = bextend0
        enddo
        do i = 1, ncmplx
          bextendx(i) = bextend0
        enddo
      else if (iact .eq. 2) then ! Davies algorithm
        adebye = half
        bdebye = one
        do j = 1, ncomp
          a0(j) = one
          bextend(j) = 0.3d0*z(j)*z(j)*half
        enddo
        do i = 1, ncmplx
          ax0(i) = one
          bextendx(i) = 0.3d0*zx(i)*zx(i)*half
        enddo
      else if (iact .eq. 5) then
        iact = 1
        bextend0 = zero
        do j = 1, ncomp
          bextend(j) = bextend0
        enddo
        do i = 1, ncmplx
          bextendx(i) = bextend0
        enddo
      endif

!-----set initial temperature and density fields
      if (mode .eq. 2) then
!       do n = 1,nmax
!         temp(n) = tempini
!         if (pref0 .gt. zero) then
!           press(n) = pref0
!           call density(temp(n),pref0,den)
!           rho(n) = den*1.d-3
!         else
!           press(n) = -pref0
!           rho(n) = one
!         endif
!       enddo

        do ibc = 1,nblkbc
!         ibc = iregn(mmbc)
          if (pref0 .gt. zero) then
            call density (tempbc(ibc),pref0,den)
            dwbc(ibc) = den*1.d-3
          else
            dwbc(ibc) = 1.d0
          endif
        enddo
      endif
      
      ireg = initreg
      if (ibcreg(2).gt.0) ireg = ireg + ibcreg(2)-ibcreg(1) + 1
      if (ibcreg(4).gt.0) ireg = ireg + ibcreg(4)-ibcreg(3) + 1

      do iireg = 1, ireg
        do j = 1, ncomp
          if (itype(j,iireg) .eq. 0) then
            ncon(j,iireg) = 'log'
          else if (itype(j,iireg) .eq. 1) then
            ncon(j,iireg) = 'Total'
          else if (itype(j,iireg) .eq. 2) then
            ncon(j,iireg) = 'Aq+Sorp'
          else if (itype(j,iireg) .eq. 21) then
            ncon(j,iireg) = 'Aq+CEC'
          else if (itype(j,iireg) .eq. 5) then
            ncon(j,iireg) = 'Wt%'
          else if (itype(j,iireg) .eq. 7) then
            ncon(j,iireg) = 'Conc'
          else if (itype(j,iireg) .eq. 8) then
            ncon(j,iireg) = 'pH'
          else if (itype(j,iireg) .eq. 10) then
            ncon(j,iireg) = 'Constr. Qty'
          else if (itype(j,iireg) .eq. -1) then
            ncon(j,iireg) = 'Charge Bal.'
          endif
        enddo
      enddo

!-----set mineral kinetic rate law translation
      itypkini(20) = 1
      itypkini(21) = 2
      itypkini(22) = 3
      itypkini(25) = 4
      itypkini(30) = 5
      itypkini(90) = 6
      itypkini(91) = 7

!-----set aqueous kinetic rate law translation
      itypkiniaq(19) =  1
      itypkiniaq(20) =  2
      itypkiniaq(21) =  3
      itypkiniaq(22) =  4
      itypkiniaq(23) =  5
      itypkiniaq(25) =  6
      itypkiniaq(30) =  7
      itypkiniaq(60) =  8
      itypkiniaq(90) =  9
      itypkiniaq(91) = 10
      itypkiniaq(95) = 11
      
!-----load irreversible mineral properties
      do nr = 1, nkin
        do m = 1, mnrl
          if (namk(nr) .eq. namrl(m)) then
            ndxkin(nr) = m
            eqkin(nr) = alnk(m)
!           ze(nr) = ze0(m)
            do j = 1, ncpri
              skin(j,nr) = smnrl(j,m)
            enddo

!-----------store molar volume - liter/mole
            vbarkin(nr) = vbar(m)
            wtkin(nr)   = wtmin(m)
            goto 20
          endif
        enddo
        if (myrank == 0) write(*,*) 'mineral name not found: ',namk(nr)
        stop
   20   continue
   
!-----rate constant - mole/cm**2/sec

!-----convert units of kinetic reaction rate constant
!     aeqkin=ten**eqkin(nr)

        if (rlim0(nr) .lt. zero) rlim0(nr) = ten**rlim0(nr)
!       rlim(nr) = 1.d4*rlim0(nr)
        rlim(nr) = 1.d2*rlim0(nr) ! mol/cm^2/s -> mol/dm^2/s
        do l = npar1(nr), npar2(nr)
          if (rkf00(l) .lt. zero) rkf00(l) = ten**rkf00(l)
!         rkf0(l)  = 1.d4*rkf00(l) ! mol/cm^2/s -> mol/m^2/s
!         rkf0(l)  = 1.d2*rkf00(l) ! mol/cm^2/s -> mol/dm^2/s
          rkf0(l)  = 1.d1*rkf00(l) ! mol/cm^2/s -> mol/L m/s
          rkf(l)   = rkf0(l)
          !print *,'ptran_init: ',nr,l,npar1(nr),npar2(nr),rkf(l)
        enddo
        if (rkfa00(nr) .lt. zero) rkfa00(nr) = ten**rkfa00(nr)
        if (rkfb00(nr) .lt. zero) rkfb00(nr) = ten**rkfb00(nr)
!       rkfa0(nr) = rkfa00(nr)*1.d4
!       rkfb0(nr) = rkfb00(nr)*1.d4
        rkfa0(nr) = rkfa00(nr)*1.d2
        rkfb0(nr) = rkfb00(nr)*1.d2
        rkfa(nr)  = rkfa0(nr)
        rkfb(nr)  = rkfb0(nr)
      enddo

!-----compress stoichiometric matrix storage
!-----homogeneous reactions

!-----store product of stoichiometric coefficients
      jj = 0
      do l = 1, ncomp
        do j = l, ncomp
          jj = jj + 1
          do i = 1, ncmplx
            sshom(i,jj) = shom(j,i)*shom(l,i)
          enddo
        enddo
      enddo

!-----compute nonzero dpsi(j,l) matrix indices: store in jpsi
      ncsqdp = 0
      jj = 0
      do l = 1, ncomp
        do j = 1, ncomp
          jj = jj + 1
          fjl = zero
          if (j.ne.l) then ! always include diagonal terms
            do i = 1, ncmplx
              fjl = shom(j,i)*shom(l,i)
              if (fjl .ne. zero) goto 555
            enddo
            goto 556
          endif
  555     continue
          ncsqdp = ncsqdp + 1
          jpsi(ncsqdp) = jj
  556     continue
        enddo
      enddo

!-----compute nonzero dpsig(j,l) matrix indices: store in jpsig
      if (ngas > 0 .and. iphase == 2) then
      ncsqdpg = 0
      jj = 0
      do l = 1, ncomp
        do j = 1, ncomp
          jj = jj + 1
          fjl = zero
          if (j.ne.l) then ! always include diagonal terms
            do i = 1, ngas
              fjl = sgas(j,i)*sgas(l,i)
              if (fjl .ne. zero) goto 557
            enddo
            goto 558
          endif
  557     continue
          ncsqdpg = ncsqdpg + 1
          jpsig(ncsqdpg) = jj
  558     continue
        enddo
      enddo
      endif
      
      ncsq = ncomp*ncomp
!     write(iunit2,*) 'trdatint-storage: ',ncsq,ncsqdp,ncsqdpg
!     do jj = 1, ncsqdp
!       jjj = jpsi(jj)
!       write(*,*) 'trdatint: ',jj,jjj,ncsqdp,ncomp*ncomp
!     enddo

      do j = 1, ncomp
        k = 0
        do i = 1, ncmplx
          if (shom(j,i).ne.zero) then
            k = k + 1
            fshom(j,k) = shom(j,i)
            ki(j,k) = i
          endif
        enddo
        lc(j) = k

        do l = j+1, ncomp
          kk = 0
          do i = 1, ncmplx
            if (shom(j,i).ne.zero .and. shom(l,i).ne.zero) then
              kk = kk + 1
              kki(j,l,kk) = i
              ffshom(j,l,kk) = shom(j,i)*shom(l,i)
            endif
          enddo
          llc(j,l) = kk
        enddo
      enddo

!-----gaseous reactions

!-----store product of stoichiometric coefficients
      if (ngas > 0 .and. iphase == 2) then
      k = 0
      do l = 1, ncomp
        do j = l, ncomp
          k = k + 1
          do i = 1, ngas
            ssgas(i,k) = sgas(j,i)*sgas(l,i)
          enddo
        enddo
      enddo

      do i = 1, ngas
        njg(i) = 0
        do j = 1, ncomp
          if (sgas(j,i).ne.zero) then
            njg(i) = njg(i) + 1
            jg(njg(i),i) = j
          endif
        enddo
      enddo

      do j = 1, ncomp
        k = 0
        do i = 1, ngas
          if (sgas(j,i).ne.zero) then
            k = k + 1
            fsgas(j,k) = sgas(j,i)
            kg(j,k) = i
          endif
        enddo
        lg(j) = k

        do l = j+1, ncomp
          kk = 0
          do i = 1, ngas
            if (sgas(j,i).ne.zero .and. sgas(l,i).ne.zero) then
              kk = kk + 1
              kkg(j,l,kk) = i
              ffsgas(j,l,kk) = sgas(j,i)*sgas(l,i)
            endif
          enddo
          llg(j,l) = kk
        enddo
      enddo
      endif
      
!-----mineral reactions
      do nr = 1, nkin
        k = 0
        do j = 1, ncomp
          if (skin(j,nr) .ne. zero) then
            k = k + 1
            kinj(nr,k) = j
            fskin(nr,k) = skin(j,nr)*skin(j,nr)
          endif
        enddo
        kmind(nr) = k
        k = 0
        do j = 1, ncomp
          do l = j+1, ncomp
            if (skin(j,nr).ne.zero .and. skin(l,nr).ne.zero) then
              k = k + 1
              kminj(nr,k) = j
              kminl(nr,k) = l
              ffskin(nr,k) = skin(j,nr)*skin(l,nr)
            endif
          enddo
        enddo
        kkmin(nr) = k
      enddo

      idebug = 0
      if (idebug > 0 .and. myrank == 0) then
        write(iunit2,*)
        write(iunit2,*) ' stoichiometrix matrix compression'
        kcount  = 0
        kcountp = 0
        do j = 1, ncomp
          write(iunit2,'(1x,a12,1x,"no. primary species: ",i4)') &
          nam(j),lc(j)
          write(iunit2,*) ' complex  stoichiometric coef.'
          do k = 1, lc(j)
            kcount = kcount + 1
            write(iunit2,'(a12,1x,1pe12.4)') namcx(ki(j,k)),fshom(j,k)
          enddo
          do l = j+1, ncomp
            write(iunit2,'("no. pairs: ",2x,2a12,1x,i4)') nam(j), &
            nam(l),llc(j,l)
            write(iunit2,*) ' complex  stoichiometric coef.'
            do k = 1, llc(j,l)
              kcountp = kcountp + 1
              write(iunit2,'(a12,1x,1pe12.4)') namcx(kki(j,l,k)), &
              ffshom(j,l,k)
            enddo
          enddo
        enddo

        write(iunit2,*) ' nonzero terms: ',kcount,' pairs: ',kcountp
        write(iunit2,*) ' nc x ncx = ',ncomp*ncmplx, &
        ', nc x nc x ncx = ',ncomp*ncomp*ncmplx
        write(iunit2,'()')
      endif

      do i = 1, ncmplx
        k = 0
        do j = 1, ncomp
          if (shom(j,i).ne.zero) then
            k = k + 1
            cshom(k,i) = shom(j,i)
            jcmpr(k,i) = j
          endif
        enddo
        ncmpr(i) = k
      enddo

!-----water stoichiometric coefficients
      if (myrank == 0) then
        write(iunit2,*)
        write(iunit2,*) 'stoichiometric coefficients for H2O'
        write(iunit2,*) 'species                   nH2O'
        if (jh2o.eq.0) then
          jaq = ncomp+1
        else
          jaq = jh2o
        endif
        do i = 1, ncmplx
          write(iunit2,'(1x,a20,1x,1pe12.4)') namcx(i),shom(jaq,i)
        enddo
        do i = 1, ngas
          write(iunit2,'(1x,a20,1x,1pe12.4)') namg(i),sgas(jaq,i)
        enddo
        do i = 1, nkin
          write(iunit2,'(1x,a20,1x,1pe12.4)') namk(i),skin(jaq,i)
        enddo
        do i = 1, nsrfmx
          write(iunit2,'(1x,a20,1x,1pe12.4)') namscx(i),ssorp(jaq,i)
        enddo
        write(iunit2,*)
      endif

!-----store mineral kinetic prefactor stoichiometry in proper order
      do m = 1, nkin
        do ml = npar1(m), npar2(m)

!---------set primary species prefactor coefficients
          do j = 1, nkinpri(ml)
            do l = 1, ncomp
              if (namprik(j,ml) .eq. nam(l)) then
                jpri(j,ml) = l
                goto 222
              endif
            enddo
            write(*,*) 'error finding prefactor primary species: stop', &
            myrank
            stop
  222       continue
          enddo

!---------set secondary species prefactor coefficients
          do i = 1, nkinsec(ml)
            do l = 1, ncmplx
              if (namseck(i,ml) .eq. namcx(l)) then
                isec(i,ml) = l
                goto 223
              endif
            enddo
            write(*,*) 'error finding prefactor secondary species: stop', &
            myrank
            stop
  223       continue
          enddo
        enddo
      enddo

   if (myrank == 0) then
      if (nkin > 0) then
        write(iunit2,'(/,"parallel reaction stoichiometry")')
        write(iunit2,'("mineral",10x,"      npri  nsec")')
      endif
      do nr = 1, nkin
        do lp = npar1(nr), npar2(nr)
          write(iunit2,'(a20,2i6)') namk(nr),nkinpri(lp),nkinsec(lp)
          do l = 1, nkinpri(lp)
            write(iunit2,'(10x,a20,1pe12.4)') nam(jpri(l,lp)), &
            skinpri(l,lp)
          enddo
          do l = 1, nkinsec(lp)
            write(iunit2,'(10x,a20,1pe12.4)') namcx(isec(l,lp)), &
            skinsec(l,lp)
          enddo
        enddo
      enddo
      write(iunit2,*)
   endif

!-----setup ion-exchange ordering of minerals, colloids, and cations
      ncollex = 0
      do m = 1, nexsolid
!-------set species indices for minerals
        do ns = 1, nkin
          if (namex(m) .eq. namk(ns)) then
            mex(m) = ns
            goto 50
          endif
        enddo

!-------set species indices for ion exchange colloids
        do j = ncomp-ncoll+1, ncomp
          if (namex(m) .eq. nam(j)) then
            mex(m) = j
            ncollex = ncollex + 1
            goto 50
          endif
        enddo

        write(*,*) 'error in trdatint setting ion exchange colloid ', &
        'and mineral indices: stop'
        stop

   50   continue

!-------set species indices for cations
        do jj = 1, nexmax
          do l = 1, ncomp
            if (namcat(jj) .eq. nam(l)) then
              jex(jj) = l
              goto 60
            endif
          enddo
          write(*,*) 'error in trdatint setting cation exchange ', &
          'indices: stop'
          stop
   60     continue
        enddo
      enddo

!-----determine nr. of colloid surface complexation and exchange solids
      ncolsrf = 0
      do m = 1, nsrfmin
        do mm = 1, ncoll
          if (namcoll(mm).eq.namsrf(m)) then
            ncolsrf = ncolsrf + 1
          endif
        enddo
      enddo
      
      do i = 1, nsrfmx
!       print *,i,nsrfmx,namscx(i),eqsorp0(i),eqsorp(i)
        if (eqsorp0(i).ne.zero) eqsorp(i) = eqsorp0(i)
      enddo

!-----setup mineral ordering for surface complexation reactions
      do mm = 1, nsrfmin
        do m = 1, nkin
          if (namsrf(mm) .eq. namk(m)) then
            msorp(mm) = m
            if (wtkin(m).le.zero .or. vbarkin(m).le.zero) then
              write(*,'(''error: zero formula weight or molar volume '', &
    &         ''for mineral '',a20,'': '',1p2e12.4)') namk(m),wtkin(m), &
              vbarkin(m)
              stop
            endif
            goto 70
          endif
        enddo

!-------set species indices for surface complexation colloids
        do j = ncomp-ncoll+1, ncomp
          if (namsrf(mm) .eq. nam(j)) then
            msorp(mm) = j
            if (wt(j).le.zero) then
              write(*,'(''error: zero formula weight '', &
    &         ''for colloid '',a20,'': '',1p2e12.4)') nam(j),wt(j)
              stop
            endif
            goto 70
          endif
        enddo

        write(*,*) 'error in trdatint: surface complex mineral not ', &
        'found! STOP',namsrf(mm),nam(j)
        stop
   70   continue
      enddo

  end subroutine trinit

end module ptran_init_module
