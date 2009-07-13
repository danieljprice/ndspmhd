!--------------------
!  SPLASH MAIN MENU
!--------------------
module mainmenu
 implicit none
 public :: menu,allowrendering,set_coordlabels,set_extracols
 integer, private :: icoordsprev = -1
 
 private

contains

subroutine menu
  use filenames, only:defaultsfile,limitsfile,animfile,fileprefix,set_filenames
  use labels, only:label,labelvec,iamvec,isurfdens,itoomre,ipdf,icolpixmap
  use limits, only:write_limits,lim2,lim,reset_lim2,lim2set
  use options_data, only:submenu_data
  use settings_data, only:ndim,numplot,ndataplots,nextra,ncalc,ivegotdata, &
                     buffer_data,ncolumns
  use settings_limits, only:submenu_limits,iadapt
  use settings_part, only:submenu_particleplots
  use settings_page, only:submenu_page,submenu_legend,interactive
  use settings_render, only:submenu_render,iplotcont_nomulti
  use settings_vecplot, only:submenu_vecplot,iplotpartvec
  use settings_xsecrot, only:submenu_xsecrotate,write_animfile
  use multiplot
  use prompting, only:prompt,print_logical
  use transforms, only:transform_label
  use defaults, only:defaults_write
  use geometry, only:labelcoord
  use getdata, only:get_data
  use timestepping
  implicit none
  integer :: i,icol,ihalf,iadjust,indexi,ierr
  integer :: ipicky,ipickx,irender,ivecplot,icontourplot
  integer :: iamvecprev, ivecplottemp,ichoose
  character(len=2) :: ioption
  character(len=100) :: vecprompt
  logical :: iAllowRendering

  irender = 0
  icontourplot = 0
  ivecplot = 0
  ipickx = 1
  ipicky = 1

  menuloop: do
!---------------------------------------------------------------------------
!  preliminaries
!---------------------------------------------------------------------------
!
!--make sure the number of columns is set appropriately
!  (nextra can change depending on what options are set)
!
  !
  !--numplot is the total number of data columns (read + calculated)
  !   not including the particle co-ordinates
  !  nextra are extra graphs to plot (e.g. convergence plots, power spectrum)
  !
  !  note that numplot and ndataplots should *only* be set here
  !  this means that even if ncolumns changes during data reads while plotting
  !  we don't start plotting new quantities
  !
  call set_extracols(ncolumns,ncalc,nextra,numplot,ndataplots)
!
!--set the coordinate and vector labels 
!  if working in a different coordinate system
!
  call set_coordlabels(numplot)

!--set contents of the vector plotting prompt
  vecprompt(1:6) = '0=none'
  indexi = 7
  iamvecprev = 0
  do icol=1,numplot
     if (iamvec(icol).ne.0 .and. iamvec(icol).ne.iamvecprev) then
        iamvecprev = iamvec(icol)
        if (iamvec(icol).ge.10) then
           write(vecprompt(indexi:),"(',',1x,i2,'=',a)") &
                 iamvec(icol),trim(labelvec(icol))
        else
           write(vecprompt(indexi:),"(',',1x,i1,'=',a)") &        
                 iamvec(icol),trim(labelvec(icol))
        endif
        indexi = len_trim(vecprompt) + 1
     endif
  enddo 

  ichoose = 0

!---------------------------------------------------------------------------
!  print menu
!---------------------------------------------------------------------------

  if (numplot.gt.0) then
!
!--data columns
!
     print "(/a)",' You may choose from a delectable sample of plots '
     print 12
     ihalf = numplot/2                ! print in two columns
     iadjust = mod(numplot,2)
     print 11, (i,transform_label(label(i),itrans(i)), &
          ihalf + i + iadjust, transform_label(label(ihalf + i + iadjust), &
          itrans(ihalf+i+iadjust)),i=1,ihalf)
     if (iadjust.ne.0) then
        print 13, ihalf + iadjust,transform_label(label(ihalf + iadjust), &
             itrans(ihalf+iadjust))
     endif 
!
!--multiplot
!  
     print 12
     print 18,numplot+1,'multiplot ',nyplotmulti,'m','set multiplot '
  else
!
!--if no data
!
     print "(/a)",' No data: You may choose from the options below '
  endif
  
11 format(1x,i2,')',1x,a20,1x,i2,')',1x,a20)
12 format(55('-'))
13 format(1x,i2,')',1x,a20)
18 format(1x,i2,')',1x,a,'[ ',i2,' ]',5x,a2,') ',a)

!
!--options 
! 
  print 12
  if (ndim.le.1) then
     print "(a)",' d(ata) p(age) o(pts) l(imits) le(g)end s,S(ave) q(uit)'  
  else
     print "(a)",' d(ata) p(age) o(pts) l(imits) le(g)end h(elp)'
     print "(a)",' r(ender) v(ector) x(sec/rotate) s,S(ave) q(uit)'
  endif
  print 12

!
!--prompt user for selection
!
  write(*,"(a)",ADVANCE='NO') 'Please enter your selection now (y axis or option):'
  read(*,*,iostat=ierr) ioption
  if (ierr < 0) stop 'reached end of input' ! end of input (e.g. in script)
  if (ierr > 0) stop !'error reading input' 

!------------------------------------------------------------
!  if input is an integer and within range, plot data
!------------------------------------------------------------
  read(ioption,*,iostat=ierr) ipicky
  if (ierr /= 0) ipicky = -1

  if (ipicky.gt.0 .and. ipicky.le.numplot+1) then
     
     if (.not.ivegotdata) then
        !
        !--do not allow plotting if no data - instead try to read data
        !
        print*,' no data '
        if (buffer_data) then
           call get_data(-1,.false.)
        else
           call get_data(1,.false.,firsttime=.true.)
        endif
     else        
        !
        !--if needed prompt for x axis selection
        !
        if (ipicky.le.(numplot-nextra)) then
           if (ipickx.eq.0) ipickx = 1 ! do not allow zero as default
           call prompt(' (x axis) ',ipickx)
           !--go back to y prompt if out of range
           if (ipickx.gt.numplot .or. ipickx.le.0) cycle menuloop
           !
           !--work out whether rendering is allowed
           !
           iAllowRendering = allowrendering(ipickx,ipicky)
           !
           !--prompt for render and vector plots 
           ! -> only allow if in "natural" coord system, otherwise h's would be wrong)
           ! (a future feature might be to interpolate in icoord then translate the pixels
           !  to icoordsnew, or alternatively plot non-cartesian pixel shapes)
           ! -> also do not allow if transformations are applied
           !
           if (ipicky.le.ndim .and. ipickx.le.ndim .and. iAllowRendering) then
              call prompt('(render) (0=none)',irender,0,numplot)
              if (irender.gt.0 .and. iplotcont_nomulti) then
                 call prompt('(contours) (0=none)',icontourplot,0,numplot)
                 if (icontourplot.eq.irender) then
                    if (iadapt) then
                       print "(a)",' contour limits are adaptive '
                    else
                       if (.not.lim2set(icontourplot)) lim2(icontourplot,:) = lim(icontourplot,:)
                       call prompt(' enter min for contours:',lim2(icontourplot,1))
                       call prompt(' enter max for contours:',lim2(icontourplot,2))
                       if (all(lim2(icontourplot,:).eq.lim(icontourplot,:))) call reset_lim2(icontourplot)
                    endif
                 endif
              endif
              if (any(iamvec(1:numplot).ne.0)) then
                 ivecplottemp = -1
                 ierr = 1
                 do while(ierr.ne.0 .and. ivecplottemp.ne.0)
                    ivecplottemp = ivecplot
                    ierr = 0
                    call prompt('(vector plot) ('//trim(vecprompt)//')',ivecplottemp,0,maxval(iamvec))
                    if (.not.any(iamvec(1:numplot).eq.ivecplottemp)) then
                       print "(a)",'Error, value not in list'
                       ierr = 1
                    endif
                 enddo
                 ivecplot = ivecplottemp
              else
                 ivecplot = 0
              endif
              if (ivecplot.gt.0 .and. irender.eq.0) then
                 call prompt('plot particles?',iplotpartvec)
              endif
           else
              irender = 0
              ivecplot = 0
           endif
        elseif (ipicky.gt.0 .and. ipicky.eq.itoomre .or. ipicky.eq.isurfdens) then
            if (ipicky.eq.isurfdens) print "(a)",' setting x axis to r for surface density plot'
            if (ipicky.eq.itoomre) print "(a)",' setting x axis to r for Toomre Q plot'
            ipickx = 1
            irender = 0
            ivecplot = 0
        elseif (ipicky.gt.0 .and. ipicky.eq.ipdf) then
            call prompt(' enter x axis for PDF calculation ',ipickx,1,ndataplots)
            irender = 0
            ivecplot = 0
        elseif (ipicky.gt.0 .and. ipicky.eq.icolpixmap) then
            call prompt(' enter corresponding SPH column for particle data ',irender,0,ndataplots)
            ipickx = 0
            ivecplot = 0
        endif
        !
        !--call main plotting routine
        !
        call timestep_loop(ipicky,ipickx,irender,icontourplot,ivecplot)
     endif
!------------------------------------------------------------------------
!  if input is an integer > numplot+1, quit
!------------------------------------------------------------------------
  elseif (ipicky.gt.numplot+1) then
     return
  else
!------------------------------------------------------------------------
!  if input is a string, use menu options
!------------------------------------------------------------------------
!--  Menu shortcuts; so you can type e.g. o2 and get the o)ptions menu, item 2
     read(ioption(2:2),*,iostat=ierr) ichoose
     if (ierr /= 0) ichoose = 0

     select case(ioption(1:1))
!------------------------------------------------------------------------
!+ Sets up plotting of (m)ultiple quantities per timestep
     case('m','M')
        call options_multiplot
!------------------------------------------------------------------------
!+ This submenu sets options relating to the (d)ata read
     case('d','D')
        call submenu_data(ichoose)
!------------------------------------------------------------------------
!+ This option turns (i)nteractive mode on/off
     case('i','I')
        interactive = .not.interactive
        print "(a)",' Interactive mode is '//print_logical(interactive)
!------------------------------------------------------------------------
!+ This submenu sets (p)age setup options
     case('p','P')
        call submenu_page(ichoose)
!------------------------------------------------------------------------
!+ This submenu sets particle plot (o)ptions
     case('o','O')
        call submenu_particleplots(ichoose)
!------------------------------------------------------------------------
!+ This submenu sets le(g)end and title options
     case('g','G')
        call submenu_legend(ichoose)
!------------------------------------------------------------------------
!+ This submenu sets (r)endering options
     case('r','R')
        if (ndim.le.1) print "(a)",'WARNING: these options have no effect in < 2D'
        call submenu_render(ichoose)
!------------------------------------------------------------------------
!+ This submenu sets (v)ector plotting options
     case('v','V')
        if (ndim.le.1) print "(a)",'WARNING: these options have no effect in < 2D'
        call submenu_vecplot(ichoose)
!------------------------------------------------------------------------
!+ This submenu sets cross section and rotation options
     case('x','X')
        if (ndim.le.1) print "(a)",'WARNING: these options have no effect in < 2D'
        call submenu_xsecrotate(ichoose)
!------------------------------------------------------------------------
!+ This submenu sets options relating to the plot limits
     case('l','L')
        call submenu_limits(ichoose)
!------------------------------------------------------------------------
!+ The (s)ave option saves the default options to a
!+ file called `splash.defaults'' in the current directory which
!+ is read automatically upon the next invocation of splash.
!+ This file uses namelist formatting and may be edited
!+ manually prior to startup if so desired. This is quite
!+ useful for setting multiplots with many plots per page
!+ The (S)ave option writes both the defaults file and
!+ also saves the current plot limits to a file called
!+ 'splash.limits' which is also read automatically
!+ at startup.
     case('s')
        if (ioption(2:2).eq.'a') then
           call prompt('enter prefix for defaults file: ',fileprefix,noblank=.true.)
           if (index(fileprefix,'.defaults').eq.0) then
              defaultsfile = trim(fileprefix)//'.defaults'
           else
              defaultsfile = trim(fileprefix)
           endif
        endif
        call defaults_write(defaultsfile)
     case('S')
        if (ioption(2:2).eq.'a' .or. ioption(2:2).eq.'A') then
           call prompt('enter prefix for filenames: ',fileprefix,noblank=.true.)
           call set_filenames(trim(fileprefix))
        endif
        call defaults_write(defaultsfile)
        call write_limits(limitsfile)
        call write_animfile(animfile)
!------------------------------------------------------------------------
!+ Slightly obsolete: prints whatever help may be helpful
     case('h','H')
        print "(10(/a))",' Hint: menu items can be shortcut by typing, e.g. o2 for ',&
                 ' the o)ptions menu, item 2.',' ', &
                 ' for detailed help, consult the user guide',' ',&
                 '  (splash/docs/splash.pdf ',&
                 '   or http://users.monash.edu.au/~dprice/splash/userguide/)', &
                 ' ',' and/or the online FAQ. If you''re really stuck, email me! '
        read*
!------------------------------------------------------------------------
!+ (q)uit, unsurprisingly, quits. Typing a number greater
!+ than the number of data columns also exits the program
!+ (e.g. I often simply type 99 to exit).
     case('q','Q')
        return
!------------------------------------------------------------------------
     case DEFAULT
        print "(a)",'unknown option '//trim(ioption) 
     end select
     
  endif

  enddo menuloop

  return
  
 contains

!----------------------------------------------------
! multiplot setup
!----------------------------------------------------
  subroutine options_multiplot
   use settings_page, only: nacross, ndown
   use settings_render, only: iplotcont_nomulti
   use limits, only:lim,lim2,lim2set,reset_lim2
   implicit none
   integer :: ifac,ierr
   logical :: isamex, isamey, icoordplot, anycoordplot, imultisamepanel
   
   call prompt('Enter number of plots per timestep:',nyplotmulti,1,numplot)

   isamey = all(multiploty(1:nyplotmulti).eq.multiploty(1))
   if (ndim.ge.2) call prompt('Same y axis for all?',isamey)
   if (isamey) then
      call prompt('Enter y axis for all plots',multiploty(1),1,numplot)
      multiploty(2:nyplotmulti) = multiploty(1)
   endif
   
   isamex = all(multiplotx(1:nyplotmulti).eq.multiplotx(1))
   call prompt('Same x axis for all?',isamex)
   if (isamex) then
      call prompt('Enter x axis for all plots',multiplotx(1),1,numplot)
      multiplotx(2:nyplotmulti) = multiplotx(1)
   endif
   
   anycoordplot = .false.
   do i=1,nyplotmulti
      print*,'-------------- Plot number ',i,' --------------'
      if (.not.isamey) then
         call prompt(' y axis ',multiploty(i),1,numplot)
      endif
      if (multiploty(i).le.ndataplots .and. .not.isamex) then
         call prompt(' x axis ',multiplotx(i),1,ndataplots)
      else
         if (multiploty(i).eq.isurfdens) then
            print "(a)",' setting x axis to r for surface density plot'
            multiplotx(i) = 1
         elseif (multiploty(i).eq.itoomre) then
            print "(a)",' setting x axis to r for Toomre Q plot'
            multiplotx(i) = 1
         elseif (multiploty(i).eq.ipdf) then
            call prompt(' enter x axis for PDF calculation ',multiplotx(i),1,ndataplots)         
         elseif (multiploty(i).eq.icolpixmap) then
            call prompt(' enter corresponding SPH column for particle data ',irendermulti(i),0,ndataplots)
            multiplotx(i) = 0
         elseif(.not.isamex) then
            multiplotx(i) = multiploty(i)
         endif
      endif
      !
      !--work out whether rendering is allowed
      !
      iAllowRendering = allowrendering(multiplotx(i),multiploty(i))
      
      icoordplot = (multiplotx(i).le.ndim .and. multiploty(i).le.ndim)
      if (icoordplot) anycoordplot = icoordplot
      
      if (icoordplot .and.iAllowRendering) then
         call prompt('(render) (0=none)',irendermulti(i),0,numplot)
         if (irendermulti(i).gt.0 .and. iplotcont_nomulti) then
            call prompt('(contours) (0=none)',icontourmulti(i),0,numplot)
            if (icontourmulti(i).eq.irendermulti(i)) then
               if (iadapt) then
                  print "(a)",' contour limits are adaptive '
               else
                  if (.not.lim2set(icontourmulti(i))) lim2(icontourmulti(i),:) = lim(icontourmulti(i),:)
                  call prompt(' enter min for contours:',lim2(icontourmulti(i),1))
                  call prompt(' enter max for contours:',lim2(icontourmulti(i),2))
                  if (all(lim2(icontourmulti(i),:).eq.lim(icontourmulti(i),:))) call reset_lim2(icontourmulti(i))
               endif
            endif
         else
            icontourmulti(i) = 0
         endif
         !iplotcontmulti(i) = iplotcont_nomulti

         if (any(iamvec(1:numplot).gt.0)) then
            ivecplottemp = -1
            ierr = 1
            do while(ierr.ne.0 .and. ivecplottemp.ne.0)
               ivecplottemp = ivecplotmulti(i)
               ierr = 0
               call prompt('(vector plot) ('//trim(vecprompt)//')',ivecplottemp,0,maxval(iamvec))
               if (.not.any(iamvec(1:numplot).eq.ivecplottemp)) then
                  print "(a)",'Error, value not in list'
                  ierr = 1
               endif
            enddo
            ivecplotmulti(i) = ivecplottemp
         else
            ivecplotmulti(i) = 0
         endif
         if (ivecplotmulti(i).gt.0 .and. irendermulti(i).eq.0) then
            call prompt('plot particles?',iplotpartvec)
         endif
      else
         !
         !--set irender, icontour and ivecplot to zero if no rendering allowed
         !
         if (multiploty(i).ne.icolpixmap) irendermulti(i) = 0
         icontourmulti(i) = 0
         ivecplotmulti(i) = 0
      endif

      if (icoordplot .and. ndim.ge.2) then
         call prompt(' is this a cross section (no=projection)? ',x_secmulti(i))
         if (x_secmulti(i)) then
            call prompt('enter co-ordinate location of cross section slice',xsecposmulti(i))
         endif
      endif
      
   enddo
   
   if (isamex .and. .not.anycoordplot) then
      imultisamepanel = .false.
      !call prompt('plot all plots in same panel? (default is different panels)',imultisamepanel)
   else
      imultisamepanel = .false.
   endif

   if (nyplotmulti.eq.1 .or. imultisamepanel) then
      nacross = 1
      ndown = 1
      print*,'setting nacross,ndown = ',nacross,ndown 
   elseif (nyplotmulti.ne.nacross*ndown) then
      !--guess nacross,ndown based on largest factor
      ifac = nyplotmulti/2
      do while (mod(nyplotmulti,ifac).ne.0 .and. ifac.gt.1)
         ifac = ifac - 1
      end do
      if (ifac.le.1) then
         nacross = nyplotmulti/2
      else
         nacross = ifac
      endif
      if (nacross.le.0) nacross = 1
      ndown = nyplotmulti/nacross
      print*,'setting nacross,ndown = ',nacross,ndown 
   else
      print*,'nacross = ',nacross,' ndown = ',ndown
   endif
   
   return
   end subroutine options_multiplot
end subroutine menu

!----------------------------------------------
! utility function which determines whether
! or not rendering is allowed or not
!----------------------------------------------
logical function allowrendering(iplotx,iploty)
 use labels, only:ih,irho !,ipmass
 use multiplot, only:itrans
 use settings_data, only:icoords,icoordsnew,ndataplots
 use settings_render, only:icolour_particles
 implicit none
 integer, intent(in) :: iplotx,iploty
!
!--work out whether rendering is allowed based on presence of rho, h & m in data read
!  also must be in base coordinate system and no transformations applied
!
 if ((ih.gt.0 .and. ih.le.ndataplots) &
    .and.(irho.gt.0 .and. irho.le.ndataplots) &
    .and.(icoords.eq.icoordsnew .or. icolour_particles) &
    .and.((itrans(iplotx).eq.0 .and. itrans(iploty).eq.0).or.icolour_particles)) then
 
    allowrendering = .true.
 else
    allowrendering = .false.
 endif

end function allowrendering

!----------------------------------------------
! utility function which sets up the "extra"
! plot columns and returns the total number
! of allowed columns for plotting
!----------------------------------------------
subroutine set_extracols(ncolumns,ncalc,nextra,numplot,ndataplots)
 use params, only:maxplot
 use labels, only:ipowerspec,iacplane,isurfdens,itoomre,iutherm,ipdf,label,icolpixmap
 use settings_data, only:ndim,icoordsnew,ivegotdata
 use settings_part, only:iexact
 use write_pixmap, only:ireadpixmap
 implicit none
 integer, intent(in) :: ncolumns
 integer, intent(inout) :: ncalc
 integer, intent(out) :: nextra,numplot,ndataplots

 nextra = 0
 ipowerspec = 0
 iacplane = 0
 isurfdens = 0
 itoomre = 0
 if (ndim.eq.3 .and. icoordsnew.eq.2 .or. icoordsnew.eq.3) then
    nextra = nextra + 1
    isurfdens = ncolumns + ncalc + nextra
    label(isurfdens) = 'Surface density'
    if (iutherm.gt.0 .and. iutherm.le.ncolumns) then
       nextra = nextra + 1
       itoomre = ncolumns + ncalc + nextra
       label(itoomre) = 'Toomre Q parameter'
    endif
 endif
 if (ndim.eq.3) then  !--Probability Density Function
    nextra = nextra + 1
    ipdf = ncolumns + ncalc + nextra
    label(ipdf) = 'PDF'
 endif

 if (ndim.le.1) then !! .or. ndim.eq.3) then ! if 1D or no coord data (then prompts for which x)
    nextra = nextra + 1      ! one extra plot = power spectrum
    ipowerspec = ncolumns + ncalc + nextra
    label(ipowerspec) = '1D power spectrum'
 else
    ipowerspec = 0
 endif
 if (iexact.eq.4) then       ! toy star plot a-c plane
    nextra = nextra + 1
    iacplane = ncolumns + ncalc + nextra
    label(iacplane) = 'a-c plane'
 else
    iacplane = 0
 endif
 !nextra = nextra + 1
 !label(ncolumns+ncalc+nextra) = 'gwaves'
 if (ndim.ge.2) then
    if (ireadpixmap) then
       nextra = nextra + 1
       icolpixmap = ncolumns + ncalc + nextra
       label(icolpixmap) = '2D pixel map'
    endif
 endif

!
!--now that we know nextra, set the total number of allowed plots (numplot).
!
 if (ivegotdata) then
    numplot = ncolumns + ncalc + nextra
    if (numplot.gt.maxplot) then
       print "(a,i3,a)",' ERROR: total number of columns = ',numplot,' is greater '
       print "(a,i3,a)",'        than the current allowed maximum (',maxplot,').'
       print "(a)",'        This is set by the parameter "maxplot" in the params module'
       print "(a)",'        in the file globaldata.f90 -- edit this and recompile splash'
       print "(a)",'        (or email me if this happens to increase the default limit)'
       stop
    endif
    ndataplots = ncolumns + ncalc
 else
    numplot = 0
    ndataplots = 0
    ncalc = 0
 endif

 return
end subroutine set_extracols

subroutine set_coordlabels(numplot)
 use geometry, only:labelcoord
 use labels, only:label,iamvec,labelvec,ix,labeldefault
 use settings_data, only:icoords,icoordsnew,ndim,iRescale
 use settings_units, only:unitslabel
 use getdata, only:get_labels
 implicit none
 integer, intent(in) :: numplot
 integer :: i
 
 if (icoordsprev.lt.0) icoordsprev = icoordsnew

!--set coordinate and vector labels (depends on coordinate system)
! if (icoordsnew.gt.0 .or. icoords.gt.0) then
 if (icoordsnew.ne.icoords) then
    do i=1,ndim
       label(ix(i)) = labelcoord(i,icoordsnew)
       if (iRescale .and. icoords.eq.icoordsnew) then
          label(ix(i)) = trim(label(ix(i)))//trim(unitslabel(ix(i)))
       endif
    enddo
 elseif (icoordsnew.ne.icoordsprev) then
    call get_labels
 endif
!
!--set vector labels if iamvec is set and the labels are the default
! 
 if (icoordsnew.gt.0) then
    do i=1,numplot
       if (iamvec(i).ne.0 .and. &
          (icoordsnew.ne.icoords .or. icoordsnew.ne.icoordsprev &
           .or. index(label(i),trim(labeldefault)).ne.0)) then
          label(i) = trim(labelvec(iamvec(i)))//'\d'//trim(labelcoord(i-iamvec(i)+1,icoordsnew))
          if (iRescale) then
             label(i) = trim(label(i))//'\u'//trim(unitslabel(i))
          endif
       endif
    enddo
 endif
 icoordsprev = icoordsnew
 
 return
end subroutine set_coordlabels

end module mainmenu
