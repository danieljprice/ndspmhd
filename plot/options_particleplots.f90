!!
!! sub-menu with utilities relating to particle plots
!!
subroutine options_particleplots
  use exact
  use labels
  use settings
  use particle_data
  use prompting
  implicit none
  integer :: iaction,n,itype
  character(LEN=1) :: ans

  iaction = 0      
  print 10, iplotline,iplotlinein,iplotav,ilabelpart,plotcirc, &
        iplotpartoftype,imarktype,icoordsnew,iexact
10  format(' 0) exit ',/, 		&
         ' 1) toggle plot line                   ( ',L1,',',1x,L1,' ) ',/, &
         ' 2) toggle plot average line           ( ',L1,' ) ',/,           &
         ' 3) toggle label particles             ( ',L1,' ) ',/,           &
         ' 4) toggle circles of interaction      ( ',L1,' ) ',/,           &
         ' 5) toggle plot particles by type      ( ',6(L1,',',1x)' )',/,  &
         ' 6) change graph markers for each type ( ',6(i2,',',1x)' )',/,  &
         ' 7) change coordinate systems          ( ',i2,' ) ',/,           &
	 ' 8) toggle exact solution              ( ',i2,' ) ')
    call prompt('enter option',iaction,0,8)
!
  select case(iaction)

!------------------------------------------------------------------------
  case(1)
     print*,' Plot initial only(i), all(a), both(b) or not (n)?'
     read*,ans
     iplotline = .false.
     iplotlinein = .false.
     if (ans.eq.'i'.or.ans.eq.'b') iplotlinein = .true.
     if (ans.eq.'a'.or.ans.eq.'b') iplotline = .true.
     if (iplotlinein) then
        call prompt('Enter PGPLOT line style',linestylein,0,5)
     endif
     print*,' Plot line = ',iplotline,iplotlinein
     return 
!------------------------------------------------------------------------
  case(2)
     iplotav=.not.iplotav
     if (iplotav) then
        call prompt('Enter no. of bins for averaging ',nbins,1,1000)
     endif
     print*,' Plot average, nbins = ',iplotav,nbins
     return 		  	    	  	  
!-----------------------------------------------------------------------
  case(3)
     !	  label particles with particle numbers
     ilabelpart=.not.ilabelpart
     print*,' label particles = ',ilabelpart
     return 	  
!------------------------------------------------------------------------
  case(4)
     plotcirc=.not.plotcirc
     print*,' Plot circles of interaction = ',plotcirc
     if (plotcirc) then	     
        call prompt('Plot all circles?',plotcircall)	     
        if (.not.plotcircall) then
           call prompt('Enter number of circles to draw',ncircpart, &
                1,size(icircpart))
           do n=1,ncircpart
              call prompt('Enter particle number to plot circle around', &
                   icircpart(n),1,maxval(ntot))
           enddo
        endif
     endif
     return 	  
!------------------------------------------------------------------------
  case(5)
     !	  plot particles by type?
     do itype=1,ntypes
        call prompt('Plot '//trim(labeltype(itype))//' particles?',iplotpartoftype(itype))
     enddo
     return 	  
!------------------------------------------------------------------------
  case(6)
     print*,'(0: Square 1: . 2: + 3: * 4: o 5: x 17: bold circle)'
     do itype=1,ntypes
        call prompt(' Enter PGPLOT marker for '//trim(labeltype(itype)) &
	     //' particles:',imarktype(itype))
     enddo
     return
!------------------------------------------------------------------------
  case(7)
     print 20,icoords
20   format(' 0) reset (=',i2,')',/, &
	    ' 1) cartesian ',/,            &
            ' 2) cylindrical polars ',/,   &
	    ' 3) spherical polars ')
     call prompt(' Enter coordinate system to plot in:',icoordsnew,0,3)
     if (icoordsnew.eq.0) icoordsnew = icoords
     return
!------------------------------------------------------------------------
  case(8)
     call submenu_exact(iexact)
     return
!------------------------------------------------------------------------
  case default
     return

  end select

  return      
end subroutine options_particleplots
