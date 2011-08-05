module timestep_plotting
  implicit none

  integer, private :: ninterp
  integer, private :: iplotx,iploty,iplotz,irenderplot,ivectorplot,ivecx,ivecy
  integer, private :: nyplots,npartdim      
  integer, private :: ngrid
  integer, private :: just, ntitles
  integer, private :: iplots,ipanel

  real, dimension(:), allocatable, private :: datpix1D, xgrid
  real, dimension(:,:), allocatable, private :: datpix
  real, dimension(:,:,:), allocatable, private :: datpix3D
  real, private :: xmin,xmax,ymin,ymax,zmin
  real, private :: rendermin,rendermax,vecmax
  real, private :: dz,zslicepos,dobserver,dscreenfromobserver
  real, private :: charheight
  real, private :: dxgrid,xmingrid,xmaxgrid
  real, private :: angletempx, angletempy, angletempz

  logical, private :: iplotpart,iplotcont,x_sec,isamexaxis,isameyaxis
  logical, private :: inewpage, tile_plots, lastplot
  logical, private :: imulti,iChangeRenderLimits,irerender,iAllowspaceforcolourbar
  
  public :: initialise_plotting, plotstep
  private

contains

!
! initialise plotting options
! called once for all steps
!
subroutine initialise_plotting(ipicky,ipickx,irender_nomulti)
  use params
  use colours, only:colour_set
  use labels, only:label,ipowerspec,ih,ipmass
  use limits, only:lim
  use multiplot, only:multiplotx,multiploty,irendermulti,nyplotmulti,x_secmulti
  use prompting
  use titles, only:read_titles,read_steptitles
  use settings_data, only:ndim,numplot
  use settings_page, only:nacross,ndown,ipapersize,tile,papersizex,aspectratio,&
                     colour_fore,colour_back,iadapt,iadaptcoords,linewidth
  use settings_part, only:linecolourthisstep,linecolour,linestylethisstep,linestyle
  use settings_render, only:icolours,iplotcont_nomulti,iPlotColourBar
  use settings_xsecrot, only:xsec_nomulti,xsecpos_nomulti,flythru,nxsec, &
                        use3Dperspective,use3Dopacityrendering,zobserver,dzscreenfromobserver,taupartdepth,rkappa
  use settings_powerspec, only:options_powerspec
  use particle_data, only:npartoftype
  use projections3D, only:coltable
  implicit none
  real, parameter :: pi=3.1415926536
  integer, intent(in) :: ipicky,ipickx,irender_nomulti
  integer :: i,j,ierr,ifirst,iplotzprev
  logical :: iadapting,iamrendering,icoordplot,iallrendered
  real :: hav,pmassav
  
  !------------------------------------------------------------------------
  ! initialisations
  !------------------------------------------------------------------------

  isamexaxis = .true.  ! same x axis on all plots? (only relevant for >1 plots per page)
  isameyaxis = .true.  ! same y axis on all plots?
  tile_plots = .false.
  iplots = 0 ! counter for how many plots have been plotted in total
  ipanel = 0  ! counter for which panel we are in on plotting page
  irenderplot = 0
  ivectorplot = 0
  x_sec = xsec_nomulti
  iplotcont = iplotcont_nomulti
  lastplot = .false.
  iplotpart = .true.
  iChangeRenderLimits = .false.
  irerender = .false.
  
  xmin = 0.
  xmax = 0.
  ymin = 0.
  ymax = 0.

  if (ndim.eq.1) x_sec = .false. ! can't have xsec in 1D
  nxsec = 1

  iamrendering = .false.
  if (ipicky.eq.numplot+1) then   ! multiplot
     imulti = .true.
     nyplots = nyplotmulti
     !
     !--if doing multiplot can only take a single cross section slice      
     !
     flythru = .false.
     nxsec = 1
     !
     !--work out whether to tile plots and make labelling decisions
     !
     if (any(multiplotx(1:nyplotmulti).ne.multiplotx(1))) isamexaxis = .false.
     if (any(multiploty(1:nyplotmulti).ne.multiploty(1))) isameyaxis = .false.
     if (any(irendermulti(1:nyplotmulti).gt.ndim)) iamrendering = .true.
     if (any(x_secmulti(1:nyplotmulti))) x_sec = .true.
  else
     !
     !--or else set number of plots = 1 and use ipicky and ipickx
     !
     imulti = .false.
     nyplots = 1 
     iploty = ipicky
     iplotx = ipickx
     if (irender_nomulti.gt.ndim) iamrendering = .true.
  endif

  !------------------------------------------------------------------------
  ! initialise options to be set before plotting

  icoordplot = iploty.le.ndim .and. iplotx.le.ndim
  iallrendered = iamrendering
  iplotz = 0
  if (imulti) then
     do i=1,nyplotmulti
        if (multiplotx(i).le.ndim .and. multiploty(i).le.ndim) then
           icoordplot = .true.
           !--this check is to see if any co-ordinate plots involve just particles
           !  (if so need to initialise the cross section slice width)
           if (irendermulti(i).le.ndim) iallrendered = .false.
           iplotzprev = iplotz
           !!--work out coordinate that is not being plotted on cross-section/ 3D plots
           iplotz = 0
           if (ndim.ge.3 .and. (x_sec .or. use3Dperspective)) then
              do j=1,ndim
                 if ((multiplotx(i).ne.multiploty(i)).and. &
                     (j.ne.multiplotx(i)).and.(j.ne.multiploty(i))) iplotz = j
              enddo
              !--use only first iplotz in the case of multiple slices
              !  (only effect is on default values for slice thickness etc below)
              if (iplotzprev.gt.0) iplotz = iplotzprev
           endif
        endif
     enddo
  elseif (icoordplot) then
     !!--work out coordinate that is not being plotted 
     if (ndim.ge.3 .and. (x_sec .or. use3Dperspective)) then
        do j=1,ndim
           if ((iplotx.ne.iploty).and. &
               (j.ne.iplotx).and.(j.ne.iploty)) iplotz = j
        enddo
     endif
  endif

  !
  !--work out whether or not to tile plots on the page
  !  if plots are coord plots, make tiling decisions based on iadaptcoords
  !  otherwise use iadapt
  !
  if (icoordplot) then
     iadapting = iadaptcoords
  else
     iadapting = iadapt
  endif
  tile_plots = tile .and. (isamexaxis.and.isameyaxis .or. isameyaxis.and.ndown.eq.1  &
                                                     .or. isamexaxis.and.nacross.eq.1)
  !--do not tile if limits are adaptive
  if (tile_plots .and. (iadapting .or. (iamrendering .and. iadapt))) then
     print "(a)",'WARNING: cannot tile plots because limits are set to adaptive'
     tile_plots = .false.
  endif
  
  !--( a further constraint on plot tiling is required in the case of 
  !    multiple renderings which would involve different colour bars )
  if (iamrendering .and. icolours.ne.0 .and. iPlotColourbar) then
     !--this option means that a margin is set aside for a colour bar on tiled plots
     iAllowspaceforcolourbar = .true.
     !--do not allow tiled plots if multiple (different) colour bars are plotted
     if (tile_plots) then
        ifirst = 0
        do i=1,nyplots
           if (irendermulti(i).gt.ndim .and. ifirst.eq.0) ifirst = i
           if (ifirst.gt.0) then
              if (irendermulti(i).gt.ndim .and. irendermulti(i).ne.irendermulti(ifirst)) then
                 if (tile_plots) print "(a)",'WARNING: cannot tile plots because of multiple colour bars'
                 tile_plots = .false.
              endif
           endif
        enddo
     endif
  else
     iAllowspaceforcolourbar = .false.
  endif  
    
  if (icoordplot) then
     if (x_sec .and. iplotz.gt.0) then
!
!--if series of cross sections (flythru), set position of first one
!
        if (flythru) then
           print 32,label(iplotz)
32         format('enter number of ',a1,' cross-section slices')
           read*,nxsec
           !!--dz is the distance between slices            
           dz = (lim(iplotz,2)-lim(iplotz,1))/float(nxsec)
           zslicepos = lim(iplotz,1) - 0.5*dz
           xsecpos_nomulti = zslicepos
        else
!
!--if single cross-section, read position of cross-section slice
!
           if (.not.imulti) call prompt(' enter '//trim(label(iplotz))// &
                       ' position for cross-section slice:', &
                       xsecpos_nomulti,lim(iplotz,1),lim(iplotz,2))
!
!--set thickness if plotting particles
!  (default thickness is half of the average particle spacing)
!
           if (.not.iallrendered) then
              npartdim = int(maxval(npartoftype(:,1))**(1./real(ndim)))
              print*,'average # of particles in each dimension = ',npartdim
              if (npartdim.gt.0) then
                 dz = (lim(iplotz,2)-lim(iplotz,1))/float(npartdim)
              else
                 dz = 0.
              endif
              if (imulti) then
                 call prompt(' enter thickness for cross section slice(s):', &
                           dz,0.0,lim(iplotz,2)-lim(iplotz,1))
              else
                 call prompt(' enter thickness of cross section slice:', &
                           dz,0.0,lim(iplotz,2)-lim(iplotz,1))           
              endif
           elseif (ndim.eq.3) then
!
!--for rendered cross sections in 3D, set thickness to 10%
!  this is the distance slices are moved up and down in interactive mode
!           
              dz = 0.1*(lim(iplotz,2)-lim(iplotz,1))
           endif
        endif ! flythru or single
     endif
     
     if (iplotz.gt.0 .and. use3Dperspective) then
!
!--initialise 3D perspective
!
       !--set default values if none set
       if (abs(zobserver).lt.tiny(zobserver)) zobserver = 10.*lim(iplotz,2)
       if (abs(dzscreenfromobserver).lt.tiny(dzscreenfromobserver)) dzscreenfromobserver = lim(iplotz,2)
       call prompt('enter z coordinate of observer ',zobserver)
       call prompt('enter distance between observer and projection screen ',dzscreenfromobserver,0.)
!
!--initialise opacity for 3D opacity rendering
!       
       if (use3Dopacityrendering .and. iamrendering) then
          hav = 0.5*(lim(ih,2) + lim(ih,1))
          pmassav = 0.5*(lim(ipmass,2) + lim(ipmass,1))
          print*,'using current h and pmass limits to calculate kappa (cross section/unit mass)'
          print*,'taking average h = ',hav,' average particle mass = ',pmassav
          print*,'[ kappa = pi*h_mean**2/(particle_mass*n_smoothing_lengths) ]'
          call prompt('enter approximate surface depth (number of smoothing lengths):',taupartdepth)          
          rkappa = pi*hav*hav/(pmassav*coltable(0)*taupartdepth)
          print*,'kappa (particle cross section per unit mass) = ',rkappa
       endif
    endif

  endif

  !!--prompt for options if plotting power spectrum      
  if (iploty.eq.ipowerspec .and. .not. imulti &
     .or. (imulti.and.any(multiploty(1:nyplotmulti).eq.ipowerspec))) then
     call options_powerspec
  endif

  !!--read step titles (don't need to store ntitles for this)
  ntitles = 0
  call read_steptitles(ntitles)
  !!--read plot titles
  ntitles = 0
  call read_titles(ntitles)

  !!------------------------------------------------------------------------
  ! initialise PGPLOT

  !!--start PGPLOT
  if (tile_plots) then
     call pgbegin(0,'?',1,1)
  else
     call pgbegin(0,'?',nacross,ndown)
  endif

  !!--set paper size if necessary
  if (ipapersize.gt.0 .and. papersizex.gt.0.0 .and. aspectratio.gt.0.0 ) then
     call pgpaper(papersizex,aspectratio)
  endif
  !!--turn off page prompting
  call pgask(.false.)
  !!if (.not. interactive) call pgbbuf !! start buffering output
  
  !!--set background/foreground colours
  ierr = 0
  if (len_trim(colour_back).gt.0) call pgscrn(0,colour_back,ierr)
  if (ierr /= 0) print 10,'background',trim(colour_back)
  if (len_trim(colour_fore).gt.0) call pgscrn(1,colour_fore,ierr)
  if (ierr /= 0) print 10,'foreground',trim(colour_fore)
10 format(' error: ',a,' colour "',a,'" not found in table')
  
  !!--set colour table
  if (iamrendering .and.(icolours.ne.0)) call colour_set(icolours)
    
  !!--set line width to something visible
  call pgslw(linewidth)
  
  linecolourthisstep = linecolour
  linestylethisstep = linestyle

end subroutine initialise_plotting

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Internal subroutines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine plotstep(ipos,istep,istepsonpage,irender_nomulti,ivecplot, &
                    npartoftype,dat,timei,gammai,ipagechange,iadvance)
  use params
  use filenames, only:nsteps
  use exact, only:exact_solution,atstar,ctstar,sigma
  use toystar1D, only:exact_toystar_ACplane
  use toystar2D, only:exact_toystar_ACplane2D
  use labels, only:label,labeltype,labelvec,iamvec, &
              ih,irho,ipmass,ix,iacplane,ipowerspec
  use limits, only:lim
  use multiplot,only:multiplotx,multiploty,irendermulti,ivecplotmulti,itrans, &
                iplotcontmulti,x_secmulti,xsecposmulti
  use particle_data, only:maxpart,icolourme
  use rotation
  use settings_data, only:numplot,ndataplots,icoords,ndim,ndimV,nfreq,iRescale,units,&
                     unitslabel,iendatstep
  use settings_limits
  use settings_part, only:icoordsnew,iexact,iplotpartoftype,imarktype,PlotOnRenderings, &
                     iplotline,linecolourthisstep,linestylethisstep
  use settings_page, only:nacross,ndown,iadapt,interactive,iaxis,iPlotLegend,iPlotStepLegend, &
                     charheightmm,iPlotTitles,vpostitle,hpostitle,fjusttitle,nstepsperpage, &
                     hposlegend,vposlegend,fjustlegend,legendtext
  use settings_render, only:npix,ncontours,icolours,iplotcont_nomulti, &
      iPlotColourBar,icolour_particles,inormalise_interpolations,ifastrender
  use settings_vecplot, only:npixvec, iplotpartvec
  use settings_xsecrot
  use settings_powerspec
!
!--subroutines called from this routine
!
  use colourparts
  use transforms
  use interactive_routines
  use geometry
  use legends, only:legend,legend_markers
  use particleplots
  use powerspectrums, only:powerspectrum
  use interpolations1D, only:interpolate1D
  use interpolations2D, only:interpolate2D, interpolate2D_xsec
  use projections3D, only:interpolate3D_projection
  use opacityrendering3D, only:interpolate3D_proj_opacity
  use xsections3D, only:interpolate3D, interpolate3D_fastxsec, &
                        interpolate3D_xsec_vec
  use titles, only:pagetitles,steptitles
  use render, only:render_pix,colourbar,colourbarfull
  use pagesetup, only:redraw_axes

  implicit none
  integer, intent(in) :: ipos, istep, istepsonpage, irender_nomulti, ivecplot
  integer, dimension(maxparttypes), intent(in) :: npartoftype
  real, dimension(:,:), intent(in) :: dat
  real, intent(in) :: timei,gammai
  logical, intent(in) :: ipagechange
  integer, intent(inout) :: iadvance
  
  integer :: ntoti,iz
  integer :: i,j,k,icolumn,irow
  integer :: nyplot
  integer :: irender,irenderpart
  integer :: npixx,npixy,npixz,ipixxsec
  integer :: npixyvec,nfreqpts,itranstemp
  integer :: index1,index2,itype,icolourprev,linestyleprev

  real, parameter :: pi = 3.1415926536
  real, parameter :: tol = 1.e-10 ! used to compare real numbers
  real, dimension(ndim) :: xcoords,vecnew
  real, dimension(max(maxpart,2000)) :: xplot,yplot,zplot
  real, dimension(maxpart) :: renderplot,hh,pmass,rho,weight
  real :: angleradx, anglerady, angleradz
  real :: rendermintemp,rendermaxtemp
  real :: xsecmin,xsecmax,dummy
  real :: pixwidth,dxfreq

  character(len=len(label(1))+20) :: labelx,labely,labelz,labelrender,labelvecplot
  character(len=120) :: title
  character(len=20) :: string,labeltimeunits
  
  logical :: iColourBar, rendering, inormalise
  
34   format (25(' -'))

  !--set labels to blank (just in case)
  labelx = ' '
  labely = ' '
  labelz = ' '
  labelrender = ' '
  labelvecplot = ' '
  xplot = 0.
  yplot = 0.
  dummy = 0.
  hh = 0.
  rho = 0.
  pmass = 0.
  labeltimeunits = ' '
  if (iReScale) labeltimeunits = unitslabel(0)
    
  !--set the arrays needed for rendering if they are present
  if (ih.gt.0 .and. ih.le.ndataplots) hh(:) = dat(:,ih)
  if (irho.gt.0 .and. irho.le.ndataplots) rho(:) = dat(:,irho)
  if (ipmass.gt.0 .and. ipmass.le.ndataplots) pmass(:) = dat(:,ipmass)
  !
  !--set number of particles to use in the interpolation routines
  !  (ie. including only gas particles and ghosts)
  !--if plotting ghost particles, set ntotplot = ntot, else ntot=npart
  !
  ntoti = sum(npartoftype)
  ninterp = npartoftype(1)
  if (labeltype(2)(1:5).eq.'ghost') then
     ninterp = ninterp + npartoftype(2)
  endif
  !
  !--set weight factor for interpolation routines
  !
  if (ipmass.gt.0) then
  !  make sure this is done in code units (ie. a consistent set)
     if (iRescale) then
        weight(1:ninterp) = (dat(1:ninterp,ipmass)/units(ipmass))/ &
                            (dat(1:ninterp,irho)/units(irho)*(dat(1:ninterp,ih)/units(ih))**ndim)
     else
        weight(1:ninterp) = (dat(1:ninterp,ipmass))/(dat(1:ninterp,irho)*dat(1:ninterp,ih)**ndim)
     endif
     inormalise = inormalise_interpolations
  else
     weight(1:ninterp) = 1.0
     inormalise = .true.
  endif
  
  !-------------------------------------
  ! loop over plots per timestep
  !-------------------------------------
  over_plots: do nyplot=1,nyplots

     if (nyplot.gt.1) print 34 
     !--make sure character height is set correctly
     call danpgsch(charheightmm,2) ! set in mm
     call pgqch(charheight) ! in PGPLOT scaled units
     iColourBar = .false.   ! should be false by default until set to true

     !--set current x, y, render and vector plot from multiplot array
     if (imulti) then
        iploty = multiploty(nyplot)
        iplotx = multiplotx(nyplot)
        irender = irendermulti(nyplot)
        ivectorplot = ivecplotmulti(nyplot)
        iplotcont = iplotcontmulti(nyplot)
        x_sec = x_secmulti(nyplot)
        zslicepos = xsecposmulti(nyplot)
     else
        irender = irender_nomulti
        ivectorplot = ivecplot
        iplotcont = iplotcont_nomulti
        x_sec = xsec_nomulti
        if (iadvance.ne.0 .and. x_sec) zslicepos = xsecpos_nomulti
     endif

     if (icolour_particles) then
        irenderpart = irender
        irenderplot = 0
     else
        irenderpart = 0
        irenderplot = irender
     endif

     if (ivectorplot.gt.0) iplotpart = iplotpartvec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! initialisation for plots of particle data
     ! copy from main dat array into xplot, yplot 
     ! also set labels and plot limits
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     initdataplots: if (iploty.le.ndataplots .and. iplotx.le.ndataplots) then
        xplot(1:ntoti) = dat(1:ntoti,iplotx)
        yplot(1:ntoti) = dat(1:ntoti,iploty)
        zplot = 0.   !--set later if x-sec
        xsecmin = 0. !-- " " 
        xsecmax = 0.
        labelx = label(iplotx)
        labely = label(iploty)
        if (iadvance.ne.0) then
           xmin = lim(iplotx,1)
           xmax = lim(iplotx,2)
           ymin = lim(iploty,1)
           ymax = lim(iploty,2)
           angletempx = anglex
           angletempy = angley
           angletempz = anglez
        endif
        !
        !--flag for whether or not we have raw particle plot or not
        !  (not allowed to use transformations on coords otherwise)
        !
        rendering = (iplotx.le.ndim .and. iploty.le.ndim .and. &
                     (irenderplot.gt.ndim .or. ivectorplot.gt.0))

        !
        !--change coordinate system if relevant
        !        
        if (icoordsnew.ne.icoords) then
           !--do this if one is a coord but not if rendering
           if ((iplotx.le.ndim .or. iploty.le.ndim).and..not.rendering) then
              print*,'changing coords from ',trim(labelcoordsys(icoords)), &
                     ' to ',trim(labelcoordsys(icoordsnew))
              do j=1,ntoti
                 call coord_transform(dat(j,ix(1:ndim)),ndim,icoords, &
                                      xcoords(1:ndim),ndim,icoordsnew)
                 if (iplotx.le.ndim) xplot(j) = xcoords(iplotx)
                 if (iploty.le.ndim) yplot(j) = xcoords(iploty)
              enddo
           endif
           if (iamvec(iplotx).gt.0) then
              if (iplotx-iamvec(iplotx)+1 .le. ndim) then
                 print*,'changing vector component from ', &
                  trim(labelcoordsys(icoords)),' to ',trim(labelcoordsys(icoordsnew))
                 do j=1,ntoti
                    call vector_transform(dat(j,ix(1:ndim)), &
                         dat(j,iamvec(iplotx):iamvec(iplotx)+ndim-1), &
                         ndim,icoords,vecnew(1:ndim),ndim,icoordsnew)
                    xplot(j) = vecnew(iplotx-iamvec(iplotx)+1)
                 enddo
              else
                 print*,'error: can''t convert vector components with ndimV > ndim'
              endif
           endif
           if (iamvec(iploty).gt.0) then
              if (iploty-iamvec(iploty)+1 .le.ndim) then
                 print*,'changing vector component from ', &
                  trim(labelcoordsys(icoords)),' to ',trim(labelcoordsys(icoordsnew))
                 do j=1,ntoti
                    call vector_transform(dat(j,ix(1:ndim)), &
                         dat(j,iamvec(iploty):iamvec(iploty)+ndim-1), &
                         ndim,icoords,vecnew(1:ndim),ndim,icoordsnew)
                    yplot(j) = vecnew(iploty-iamvec(iploty)+1)
                 enddo
              else
                 print*,'error: can''t convert vector components with ndimV > ndim'
              endif
           endif
        endif
        !--apply transformations (log, 1/x etc) if appropriate
        !  also change labels and limits appropriately
        if (.not.(rendering)) then
           if (itrans(iplotx).ne.0) then
              call transform(xplot(1:ntoti),itrans(iplotx))
              if (iaxis.eq.10 .or. iaxis.eq.30) then ! logarithmic axes
                 write(string,*) itrans(iplotx)
                 string = adjustl(string)
                 itranstemp = 0
                 if (string(len_trim(string):len_trim(string)).eq.'1') then  
                    if (len_trim(string).gt.1) read(string(1:len_trim(string)-1),*) itranstemp
                    labelx = transform_label(labelx,itranstemp)
                 else
                    labelx = transform_label(labelx,itrans(iplotx))
                 endif
              else
                 labelx = transform_label(labelx,itrans(iplotx))
              endif
              if (iadvance.ne.0) call transform_limits(xmin,xmax,itrans(iplotx))
           endif
           if (itrans(iploty).ne.0) then
              call transform(yplot(1:ntoti),itrans(iploty))
              if (iaxis.eq.20 .or. iaxis.eq.30) then
                 write(string,*) itrans(iploty)
                 string = adjustl(string)
                 itranstemp = 0
                 if (string(len_trim(string):len_trim(string)).eq.'1') then
                    if (len_trim(string).gt.1) read(string(1:len_trim(string)-1),*) itranstemp
                    labely = transform_label(labely,itranstemp)
                 else
                    labely = transform_label(labely,itrans(iploty))
                 endif              
              else
                 labely = transform_label(labely,itrans(iploty))
              endif
              if (iadvance.ne.0) call transform_limits(ymin,ymax,itrans(iploty))
           endif
        endif

        !--write username, date on plot
        !         if (nacross.le.2.and.ndown.le.2) call pgiden

        !
        !--adjust plot limits if adaptive plot limits set
        !  (find minimum/maximum only on particle types actually plotted)
        !
        if (ipagechange .and. iadvance.ne.0) then
           !--x axis
           if ((iplotx.le.ndim .and. iadaptcoords) &
           .or.(iplotx.gt.ndim .and. iadapt)) then
              xmin = huge(xmin)
              xmax = -huge(xmax)
              index1 = 1
              print*,'adapting x limits'
              do itype=1,maxparttypes
                 index2 = index1 + npartoftype(itype) - 1
                 if (iplotpartoftype(itype).and.npartoftype(itype).gt.0 &
                    .or. (iplotline.and.itype.eq.1)) then
                    xmin = min(xmin,minval(xplot(index1:index2)))
                    xmax = max(xmax,maxval(xplot(index1:index2))*scalemax)
                 endif
                 index1 = index2 + 1
              enddo
           endif
           !--y axis
           if ((iploty.le.ndim .and. iadaptcoords) &
           .or.(iploty.gt.ndim .and. iadapt)) then
              ymin = huge(ymin)
              ymax = -huge(ymax)
              index1 = 1
              print*,'adapting y limits'
              do itype=1,maxparttypes
                 index2 = index1 + npartoftype(itype) - 1
                 if (iplotpartoftype(itype).and.npartoftype(itype).gt.0 &
                    .or. (iplotline.and.itype.eq.1)) then
                    ymin = min(ymin,minval(yplot(index1:index2)))
                    ymax = max(ymax,maxval(yplot(index1:index2))*scalemax)
                 endif
                 index1 = index2 + 1
              enddo
           endif           
        endif

        !!-reset co-ordinate plot limits if particle tracking           
        if (itrackpart.gt.0 .and. iadvance.ne.0) then
           if (iplotx.le.ndim) then
              xmin = xplot(itrackpart) - xminoffset_track(iplotx)
              xmax = xplot(itrackpart) + xmaxoffset_track(iplotx)
           endif
           if (iploty.le.ndim) then
              ymin = yplot(itrackpart) - xminoffset_track(iploty)
              ymax = yplot(itrackpart) + xmaxoffset_track(iploty)           
           endif
        endif

     endif initdataplots

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! plots with co-ordinates as x and y axis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if ((iploty.le.ndim).and.(iplotx.le.ndim)) then

        npixx = npix

        !!--work out coordinate that is not being plotted         
        iz = 0
        if (ndim.ge.3) then
           do j=1,ndim
              if ((iplotx.ne.iploty).and. &
                  (j.ne.iplotx).and.(j.ne.iploty)) iz = j
           enddo
        endif
        
        !!iplotz = 0
        iplotz = iz ! this is used as cross sectioned quantity
        if (iplotz.gt.0 .and. iplotz.le.ndataplots) then
           zplot(1:ntoti) = dat(1:ntoti,iplotz)
           labelz = label(iplotz)
        endif

        if (iadvance.ne.0) then
           !--set 3D perspective values to off (zero) by default
           dscreenfromobserver = 0.
           dobserver = 0.
           irerender = .false.
           !
           !--this is a flag to say whether or not rendered limits have been changed
           !  interactively. False by default, but must retain value whilst
           !  iadvance = 0
           !
           iChangeRenderLimits = .false.
        endif
        !
        !--rotate the particles about the z (and y) axes
        !  only applies to particle plots at the moment
        !                
        if (ndim.ge.2 .and. (irotate .or. (ndim.eq.3 .and.use3Dperspective))) then
           !
           !--convert angles to radians
           !
           angleradz = angletempz*pi/180.
           anglerady = angletempy*pi/180.
           angleradx = angletempx*pi/180.
           print "(1x,a,f6.2)",'rotating particles about z by ',angletempz
           if (ndim.eq.3) then
              print "(1x,a,f6.2)",'rotating particles about y by ',angletempy
              print "(1x,a,f6.2)",'rotating particles about x by ',angletempx
           endif
           if (ndim.eq.3 .and. use3Dperspective) then
              if (iadvance.ne.0) then
                 dscreenfromobserver = dzscreenfromobserver
                 dobserver = zobserver
              endif
              if (.not.x_sec .or. x_sec.and.use3Dopacityrendering) then
                 print*,' observer height = ',dobserver,', screen at ',dobserver-dscreenfromobserver
              endif
           endif
           do j=1,ntoti
              xcoords(1:ndim) = dat(j,ix(1:ndim)) - xorigin(1:ndim)
               if (ndim.eq.2) then
                  call rotate2D(xcoords(:),angleradz)
               elseif (ndim.eq.3) then
                  call rotate3D(xcoords(:),angleradx,anglerady,angleradz,dobserver,dscreenfromobserver)
               endif
              xplot(j) = xcoords(iplotx) + xorigin(iplotx)
              yplot(j) = xcoords(iploty) + xorigin(iploty)
              if (iplotz.gt.0) then
                 zplot(j) = xcoords(iplotz) + xorigin(iplotz)
              endif
           enddo
        endif 

        !------------------------------------------------------------------
        !  rendering setup and interpolation (this is the rendering done
        !  *before* the cross sections are taken, e.g. to 3D grid)
        !------------------------------------------------------------------
        if ((irenderplot.gt.ndim).and. &
             ((ndim.eq.3).or.(ndim.eq.2.and..not.x_sec))) then
           
           !!--determine number of pixels in rendered image (npix = pixels in x direction)
           pixwidth = (xmax-xmin)/real(npix)
           npixx = int((xmax-xmin)/pixwidth) + 1
           npixy = int((ymax-ymin)/pixwidth) + 1

           !!--only need z pixels if working with interpolation to 3D grid
           !  (then number of z pixels is equal to number of cross sections)
           if ((ndim.ge.3).and.(x_sec.and.nxsec.gt.2)) then
              zmin = lim(iplotz,1)
              npixz = nxsec
           endif

           if (iadvance.ne.0 .or. irerender) then
              if (allocated(datpix)) then
                 if (npixx.ne.size(datpix(:,1)) .or. npixy.ne.size(datpix(1,:))) then
                    deallocate(datpix)           
                    allocate (datpix(npixx,npixy))
                 endif
              else
                 allocate (datpix(npixx,npixy))
              endif

              select case(ndim)
              case(2)
                 !!--interpolate to 2D grid
                 !!  allocate memory for rendering array
                 if (.not. x_sec) then
                    call interpolate2D(xplot(1:ninterp),yplot(1:ninterp), &
                         hh(1:ninterp),weight(1:ninterp),dat(1:ninterp,irenderplot), &
                         ninterp,xmin,ymin,datpix,npixx,npixy,pixwidth,inormalise)
                 endif
              case(3)
                 !!--interpolation to 3D grid - then take multiple cross sections/projections
                 !!  do this if taking more than 2 cross sections, otherwise use fast xsec
                 if (x_sec.and.nxsec.gt.2) then
                    !!--allocate memory for 3D rendering array
                    if (allocated(datpix3D)) deallocate(datpix3D)
                    allocate ( datpix3D(npixx,npixy,npixz) )
                    !!--interpolate from particles to 3D grid
                    call interpolate3D(xplot(1:ninterp),yplot(1:ninterp), &
                         zplot(1:ninterp),hh(1:ninterp),weight(1:ninterp), &
                         dat(1:ninterp,irenderplot), &
                         ninterp,xmin,ymin,zmin,datpix3D,npixx,npixy,npixz,pixwidth,dz, &
                         inormalise)
                 endif
              end select
           endif

        endif
        !
        !--if vector plot determine whether or not to plot the particles as well
        !
        iplotpart = .true.
        if (ivectorplot.gt.0) iplotpart = iplotpartvec
        if (irenderplot.gt.0) iplotpart = .false.

        !
        !%%%%%%%%%%%%%%% loop over cross-section slices %%%%%%%%%%%%%%%%%%%%%%%
        !
        over_cross_sections: do k=1,nxsec

           if (k.gt.1) print 34 

           if (x_sec) then
              !!--for multislice cross section (flythru)
              !!  increment the position of the current cross section slice          
              if (flythru) zslicepos = zslicepos + dz
              !!--for cross sections of particle plots, need range of co-ordinates in which
              !!  particles may lie
              if (iplotpart) then
                 xsecmin = zslicepos-0.5*dz
                 xsecmax = zslicepos+0.5*dz
              endif
           endif

           !------------take projections/cross sections through 3D data-----------------!
           if (irenderplot.gt.ndim .and. ndim.eq.3) then

              !!--allocate memory for 2D rendered array
              if (iadvance.ne.0) then
                 if (allocated(datpix)) then
                    if (npixx.ne.size(datpix(:,1)) .or. npixy.ne.size(datpix(1,:))) then
                       deallocate(datpix)
                       print*,'reallocating...'
                       allocate ( datpix(npixx,npixy) )
                    endif
                 else
                    print*,'allocating...'
                    allocate ( datpix(npixx,npixy) )
                 endif
              endif

              !------------------------------------------------------------------------
              ! if we have rendered to a 3D grid, take cross sections from this array
              !------------------------------------------------------------------------
              if (x_sec .and. nxsec.gt.2) then
                 ipixxsec = int((zslicepos-zmin)/dz) + 1
                 if (ipixxsec.gt.npixz) ipixxsec = npixz
                 print*,TRIM(label(iplotz)),' = ',zslicepos, &
                      ' cross section, pixel ',ipixxsec
                 datpix = datpix3D(:,:,ipixxsec)    ! slices are in 3rd dimension

              else
                 !-------------------------------------------------------------------
                 !  or do a fast projection/cross section of 3D data to 2D array
                 !-------------------------------------------------------------------
                 
                 !--set limits for opacity rendering 
                 !  (these must be known before the interpolate call)
                 if (use3Dperspective .and. use3Dopacityrendering) then
                    if (iadvance.ne.0 .or. .not.iChangeRenderLimits) then
                       if (iadapt) then
                          !!--if adaptive limits, find limits of rendered array
                          print*,'adapting render limits for opacity rendering'
                          rendermin = minval(dat(1:npartoftype(1),irenderplot))
                          rendermax = maxval(dat(1:npartoftype(1),irenderplot))
                       else
                          !!--or use fixed limits
                          rendermin = lim(irenderplot,1)
                          rendermax = lim(irenderplot,2)
                       endif
                       !!--apply transformations to limits
                       call transform_limits(rendermin,rendermax,itrans(irenderplot))
                    endif
                 endif

                 !--only rerender if absolutely necessary
                 if (iadvance.ne.0 .or. irerender) then
                    if (x_sec) then
                       if (use3Dperspective .and. use3Dopacityrendering) then
                          !!--do surface-rendered cross-section with opacity
                          print*,trim(label(ix(iplotz))),' = ',zslicepos,  &
                               ' : opacity-rendered cross section', xmin,ymin                    
                          call interpolate3D_proj_opacity( &
                               xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                               pmass(1:ninterp),hh(1:ninterp),dat(1:ninterp,irenderplot), &
                               dat(1:ninterp,iz), &
                               ninterp,xmin,ymin,datpix,npixx,npixy,pixwidth,dobserver, &
                               dscreenfromobserver,rkappa,zslicepos, &
                               rendermin,rendermax,itrans(irenderplot),istep)                    
                       elseif (use3Dperspective) then
                          print*,'ERROR: X_SEC WITH 3D PERSPECTIVE NOT IMPLEMENTED'
                          datpix = 0.
                       else
                          !!--do fast cross-section
                          print*,trim(label(ix(iplotz))),' = ',zslicepos,  &
                               ' : fast cross section', xmin,ymin
                          call interpolate3D_fastxsec( &
                               xplot(1:ninterp),yplot(1:ninterp), &
                               zplot(1:ninterp),hh(1:ninterp), &
                               weight(1:ninterp),dat(1:ninterp,irenderplot), &
                               ninterp,xmin,ymin,zslicepos,datpix,npixx,npixy,pixwidth, &
                               inormalise)
                       endif
                    else                 
                       if (use3Dperspective .and. use3Dopacityrendering) then
                          !!--do fast projection with opacity
                          call interpolate3D_proj_opacity( &
                               xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                               pmass(1:ninterp),hh(1:ninterp),dat(1:ninterp,irenderplot), &
                               dat(1:ninterp,iz), &
                               ninterp,xmin,ymin,datpix,npixx,npixy,pixwidth,dobserver, &
                               dscreenfromobserver,rkappa,huge(zslicepos), &
                               rendermin,rendermax,itrans(irenderplot),istep)                    
                       else
                          !!--do fast projection of z integrated data (e.g. column density)
                          call interpolate3D_projection( &
                               xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                               hh(1:ninterp),weight(1:ninterp),dat(1:ninterp,irenderplot), &
                               ninterp,xmin,ymin,datpix,npixx,npixy,pixwidth, &
                               dobserver,dscreenfromobserver,ifastrender)
                       endif
                    endif
                 endif

              endif ! whether 3D grid or fast renderings

              !-------------take cross sections through 2D data------------------!
           elseif (irenderplot.gt.ndim .and. ndim.eq.2 .and. x_sec) then
              !-------------------------------------------------------------------
              !  or do a fast cross section through 2D data to 1D array
              !-------------------------------------------------------------------
              !!--interpolate from 2D data to 1D line
              !!  line is specified by giving two points, (x1,y1) and (x2,y2)
              !--set up 1D grid and allocate memory for datpix1D
              if (iadvance.ne.0) then
                 xmin = 0.   ! distance (r) along cross section
                 xmax = SQRT((xseclineY2-xseclineY1)**2 + (xseclineX2-xseclineX1)**2)
              endif
              dxgrid = (xmax-xmin)/REAL(npixx)
              call set_grid1D(xmin,dxgrid,npixx)

              call interpolate2D_xsec( &
                   dat(1:ninterp,iplotx),dat(1:ninterp,iploty),&
                   hh(1:ninterp),weight(1:ninterp),dat(1:ninterp,irenderplot), &
                   ninterp,xseclineX1,xseclineY1,xseclineX2,xseclineY2, &
                   datpix1D,npixx,inormalise)
              !
              !--find limits of datpix1D for plotting
              !  do transformations on rendered array where appropriate
              !  set these as ymin,ymax and set labels of plot
              !
              call transform(datpix1D,itrans(irenderplot))
              labely = transform_label(label(irenderplot),itrans(irenderplot))
              labelx = 'cross section'
              !!--if adaptive limits, find limits of datpix
              if (iadvance.ne.0) then               
                 if (iadapt) then
                    ymin = minval(datpix1D)
                    ymax = maxval(datpix1D)
                 else
                    !!--or apply transformations to fixed limits
                    ymin = lim(irenderplot,1)
                    ymax = lim(irenderplot,2)
                    call transform_limits(ymin,ymax,itrans(irenderplot))
                 endif
              endif

           endif ! 2 or 3D and rendering
           !-----end of preliminary muff for 2D/3D cross sections/renderings ------------------

           !---------------------------------
           ! setup page
           !---------------------------------

           just = 1  ! x and y axis have same scale
           ! unless 1D xsec through 2D data or non-cartesian
           if ((irender.gt.ndim .and. ndim.eq.2 .and. x_sec) &
               .or.(icoordsnew.gt.1)) then
              just = 0 
           endif
           title = ' '
           !--work out if colour bar is going to be plotted 
           !  (leave space in page setup if so)
           iColourBar = .false.
           if (irender.gt.ndim) iColourBar = iPlotColourBar

           call page_setup
           !--only plot colour bar at the end of first row on tiled plots
           if (tile_plots .and..not.(icolumn.eq.nacross .and. irow.eq.1)) iColourBar = .false.

           !--add to log
           if (x_sec.and.iplotpart.and.iplotz.gt.0) print 35,label(iplotz),xsecmin,label(iplotz),xsecmax
35            format('cross section: ',a1,' = ',f7.3,' to ',a1,' = ',f7.3)

           !------------------------------
           ! now actually plot the data
           !------------------------------
           if (irenderplot.gt.ndim) then
              if ((ndim.eq.3).or.(ndim.eq.2.and. .not.x_sec)) then
                 !---------------------------------------------------------------
                 ! scalar quantity which has been rendered to a 2D pixel array (datpix)
                 !---------------------------------------------------------------
                 !!--do transformations on rendered array (but only the first time!)
                 if (iadvance.ne.0 .or. irerender) then
                    call transform2(datpix,itrans(irenderplot))
                 endif
                 
                 !!--set label for rendered quantity
                 labelrender = label(irenderplot)
                 !!--set label for column density (projection) plots (2268 or 2412 for integral sign)
                 if (ndim.eq.3 .and..not. x_sec .and..not.(use3Dperspective.and.use3Dopacityrendering)) then
                    labelrender = '\(2268) '//trim(labelrender)//' d'//trim(label(ix(iz)))
                    if (irenderplot.eq.irho) labelrender = 'column density'
                 endif
                 !!--apply transformations to the label for the rendered quantity 
                 labelrender = transform_label(labelrender,itrans(irenderplot))
                 
                 !!--limits for rendered quantity
                 write(string,"(i8)") itrans(irenderplot) ! used to determine whether logged or not
                 if (iadvance.ne.0 .or. .not.iChangeRenderLimits) then
                    if (iadapt .and. .not.(use3Dperspective.and.use3Dopacityrendering.and.ndim.eq.3)) then
                       !!--if adaptive limits, find limits of rendered array
                       if (index(string,'1').ne.0) then
                          rendermin = minval(datpix,mask=datpix.gt.-666.) ! see below
                       else
                          rendermin = minval(datpix)
                       endif
                       rendermax = maxval(datpix)
                       print*,'adapting render limits'
                    elseif (.not.iadapt) then
                       !!--or apply transformations to fixed limits
                       rendermin = lim(irenderplot,1)
                       rendermax = lim(irenderplot,2)
                       call transform_limits(rendermin,rendermax,itrans(irenderplot))
                    endif
                 endif
                 
                 !!--if log, then set zero values to minimum
                 !!  (must be done after minimum is known)
                 !!
                 if (index(string,'1').ne.0) then
                    !!print*,'setting zero values to ',rendermin,' on log array'
                    where (abs(datpix).lt.tiny(datpix))
                       datpix = -666. ! just some very negative number !rendermin
                    end where
                    !!  also do not let max=0 as this is suspiciously wrong
                    if (iadapt .and. abs(rendermax).lt.tiny(datpix)) then
                       !!print*,'max=0 on log plot, fixing'
                       rendermax = maxval(datpix)
                    endif
                 endif

                 !!--print plot limits to screen
                 print*,trim(labelrender),' min, max = ',rendermin,rendermax
                                       
                 !!--call subroutine to actually render the image
                 call render_pix(datpix,rendermin,rendermax,trim(labelrender), &
                   npixx,npixy,xmin,ymin,pixwidth,    &
                   icolours,iplotcont,.false.,ncontours,.false.)
                 
                 !!--plot colour bar
                 if (iColourBar) then
                    !--for tiled plots only on last plot in first row,
                    !  and use full viewport size in the y direction
                    if (tile_plots) then
                       if (icolumn.eq.nacross .and. irow.eq.1) &
                          call colourbarfull(icolours,rendermin,rendermax, &
                          trim(labelrender),.false.,nacross,ndown,iaxis)
                    else
                       !!--plot colour bar, but only if last in row
                       call colourbar(icolours,rendermin,rendermax, &
                                      trim(labelrender),.false.)
                    endif
                 endif
                 
                 !!--plot other particle types (e.g. sink particles) on top
                 call particleplot(xplot(1:ntoti),yplot(1:ntoti), &
                   zplot(1:ntoti),hh(1:ntoti),ntoti,iplotx,iploty, &
                   icolourme(1:ntoti),npartoftype(:),PlotOnRenderings(:), &
                   x_sec,xsecmin,xsecmax,labelz)

              elseif (ndim.eq.2 .and. x_sec) then
                 !---------------------------------------------------------------
                 ! plot 1D cross section through 2D data (contents of datpix) 
                 !---------------------------------------------------------------  
                 call pgline(npixx,xgrid,datpix1D)
              endif
           else
              !-----------------------
              ! particle plots
              !-----------------------
              if (iplotpart) then
                 !
                 !--sort out particle colouring
                 !
                 if (irenderpart.gt.0 .and. irenderpart.le.numplot) then
                    renderplot(1:ntoti) = dat(1:ntoti,irenderpart)
                    call transform(renderplot(1:ntoti),itrans(irenderpart))
                    labelrender = label(irenderpart)
                    labelrender = transform_label(labelrender,itrans(irenderpart))
                    
                    !!--limits for rendered quantity
                    if (iadvance.ne.0 .or. .not.iChangeRenderLimits) then
                       if (iadapt) then
                          !!--if adaptive limits, find limits of rendered array
                          rendermin = minval(renderplot(1:ntoti))
                          rendermax = maxval(renderplot(1:ntoti))
                          print*,'adapting render limits'
                       else
                          !!--or apply transformations to fixed limits
                          rendermin = lim(irenderpart,1)
                          rendermax = lim(irenderpart,2)
                          call transform_limits(rendermin,rendermax,itrans(irenderpart))
                       endif
                    endif
                    !!--print plot limits to screen
                    print*,trim(labelrender),' min, max = ',rendermin,rendermax       

                    call colour_particles(renderplot(1:ntoti), &
                         rendermin,rendermax, &
                         icolourme(1:ntoti),ntoti)

                    !!--plot colour bar
                    if (iColourBar) then
                       !--for tiled plots only on last plot in first row,
                       !  and use full viewport size in the y direction
                       if (tile_plots) then
                          if (icolumn.eq.nacross .and. irow.eq.1) &
                             call colourbarfull(icolours,rendermin,rendermax, &
                             trim(labelrender),.false.,nacross,ndown,iaxis)
                       else
                          !!--plot colour bar, but only if last in row
                          call colourbar(icolours,rendermin,rendermax, &
                                         trim(labelrender),.false.)
                       endif
                    endif
                 endif
                 !
                 !--do particle plot
                 !
                 call particleplot(xplot(1:ntoti),yplot(1:ntoti), &
                   zplot(1:ntoti),hh(1:ntoti),ntoti,iplotx,iploty, &
                   icolourme(1:ntoti),npartoftype(:),iplotpartoftype(:), &
                   x_sec,xsecmin,xsecmax,labelz)
              else
                 !!--plot other particle types on top of vector plots (e.g. sinks)
                 call particleplot(xplot(1:ntoti),yplot(1:ntoti), &
                   zplot(1:ntoti),hh(1:ntoti),ntoti,iplotx,iploty, &
                   icolourme(1:ntoti),npartoftype(:),PlotOnRenderings(:), &
                   x_sec,xsecmin,xsecmax,labelz)
                   
              endif
           endif

           !--------------------------------------------------------------
           ! vector maps (can be on top of particle plots and renderings)
           !--------------------------------------------------------------
           if (ivectorplot.ne.0 .and. ndim.ge.2) then
             if (iamvec(ivectorplot).ne.0) then
                !!--choose quantity to be plotted
                ivecx = iamvec(ivectorplot) + iplotx - 1
                ivecy = iamvec(ivectorplot) + iploty - 1

                labelvecplot = trim(labelvec(ivectorplot))//trim(unitslabel(ivectorplot))
                !!--set label for projection plots (2268 or 2412 for integral sign)
                if (ndim.eq.3 .and..not. x_sec) then
                   labelvecplot = '\(2268) '//trim(labelvecplot)//' d'//trim(label(ix(iz)))
                endif
                pixwidth = (xmax-xmin)/real(npixvec)
                npixyvec = int((ymax-ymin)/pixwidth) + 1
                if (iadvance.ne.0) then
                   if (iadapt) then
                      vecmax = -1.0  ! plot limits then set in vectorplot
                   else
                      vecmax = max(lim(ivecx,2),lim(ivecy,2))
                   endif
                endif

                call vector_plot(ivecx,ivecy,npixvec,npixyvec,pixwidth,vecmax,labelvecplot)
             endif
           endif
           !---------------------------------
           ! plot rotated axes
           !---------------------------------
           if (irotate .and. irotateaxes.gt.0) then
              if (ndim.eq.3) then
                 call rotate_axes3D(irotateaxes,iplotx,iploty, &
                      xminrotaxes(1:ndim),xmaxrotaxes(1:ndim),xorigin(1:ndim), &
                      angleradx,anglerady,angleradz,dobserver,dscreenfromobserver)
              elseif (ndim.eq.2) then
                 call rotate_axes2D(irotateaxes,xminrotaxes(1:ndim), &
                                   xmaxrotaxes(1:ndim),xorigin(1:ndim),angleradz)
              endif
           endif          
           !
           !--plot exact solution if relevant (before going interactive)
           !
           if (iexact.ne.0) then
               call exact_solution(iexact,iplotx,iploty, &
                    itrans(iplotx),itrans(iploty),icoordsnew, &
                    ndim,ndimV,timei,xmin,xmax,gammai, &
                    xplot(1:npartoftype(1)),yplot(1:npartoftype(1)), &
                    pmass(1:npartoftype(1)),npartoftype(1),imarktype(1), &
                    units(iplotx),units(iploty),irescale)
           endif
           !
           !--redraw axes over what has been plotted
           !
           call redraw_axes(iaxis)
           !
           !--annotate with time / marker legend and title
           !
           call legends_and_title
           !
           !--enter interactive mode
           !
           lastplot = ((ipos.eq.iendatstep .or. istep.eq.nsteps) &
                       .and. nyplot.eq.nyplots .and. k.eq.nxsec)

           if (interactive) then
              if (nacross*ndown.eq.1) then
                 iadvance = nfreq
                 rendermintemp = rendermin
                 rendermaxtemp = rendermax
                 iChangeRenderLimits = .false.
                 call interactive_part(ninterp,iplotx,iploty,iplotz,irender,ivecx,ivecy, &
                      xplot(1:ninterp),yplot(1:ninterp),zplot(1:ninterp), &
                      hh(1:ninterp),icolourme(1:ninterp), &
                      xmin,xmax,ymin,ymax,rendermin,rendermax,vecmax, &
                      angletempx,angletempy,angletempz,ndim,x_sec,zslicepos,dz, &
                      dobserver,dscreenfromobserver,irerender, &
                      itrackpart,icolours,iadvance,ipos,iendatstep)
                 !--turn rotation on if necessary
                 if (abs(angletempx-anglex).gt.tol) irotate = .true.
                 if (abs(angletempy-angley).gt.tol) irotate = .true.
                 if (abs(angletempz-anglez).gt.tol) irotate = .true.
                 if (abs(rendermintemp-rendermin).gt.tol) iChangeRenderLimits = .true.
                 if (abs(rendermaxtemp-rendermax).gt.tol) iChangeRenderLimits = .true.
                 if (iadvance.eq.-666) exit over_plots
              elseif ((ipanel.eq.nacross*ndown .and. istepsonpage.eq.nstepsperpage) .or. lastplot) then
                 !
                 !--timestep control only if multiple plots on page
                 !
                 iadvance = nfreq
                 call interactive_step(iadvance,ipos,iendatstep,xmin,xmax,ymin,ymax)
                 if (iadvance.eq.-666) exit over_plots
              endif
           endif

           !
           !--%%%%%%%%%%%%% end loop over cross-section slices %%%%%%%%%%%%%%%%%%%%%%%
           !
        enddo over_cross_sections

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! not both coordinates - these are just particle plots
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     elseif ((iploty.gt.ndim .or. iplotx.gt.ndim)  &
          .and.(iploty.le.ndataplots .and. iplotx.le.ndataplots)) then

        !--------------------------------
        ! setup page
        !--------------------------------
        just = 0
        title = ' '
        call page_setup

        !--------------------------------
        ! now plot particles
        !--------------------------------
        !
        !--sort out particle colouring
        !
        if (irenderpart.gt.0 .and. irenderpart.le.numplot) then
           renderplot(1:ntoti) = dat(1:ntoti,irenderpart)
           call transform(renderplot(1:ntoti),itrans(irenderpart))

           !!--limits for rendered quantity
           if (iadapt .and. iadvance.ne.0) then
              !!--if adaptive limits, find limits of rendered array
              rendermin = minval(renderplot(1:ntoti))
              rendermax = maxval(renderplot(1:ntoti))
           elseif (iadvance.ne.0) then                   
              !!--or apply transformations to fixed limits
              rendermin = lim(irenderpart,1)
              rendermax = lim(irenderpart,2)
              call transform_limits(rendermin,rendermax,itrans(irenderpart))
           endif
           !!--print plot limits to screen
           print*,trim(labelrender),' min, max = ',rendermin,rendermax       

           call colour_particles(renderplot(1:ntoti), &
                rendermin,rendermax, &
                icolourme(1:ntoti),ntoti)
        endif
        !
        !--do the particle plot
        !

        call particleplot(xplot(1:ntoti),yplot(1:ntoti), &
             zplot(1:ntoti),hh(1:ntoti),ntoti,iplotx,iploty, &
             icolourme(1:ntoti),npartoftype(:),iplotpartoftype,.false.,0.0,0.0,' ')

        if (iexact.ne.0) then
           call exact_solution(iexact,iplotx,iploty,itrans(iplotx),itrans(iploty), &
                icoordsnew,ndim,ndimV,timei,xmin,xmax,gammai, &
                xplot(1:npartoftype(1)),yplot(1:npartoftype(1)), &
                pmass(1:npartoftype(1)),npartoftype(1),imarktype(1), &
                units(iplotx),units(iploty),irescale)
        endif
        !
        !--redraw axes over what has been plotted
        !
        call redraw_axes(iaxis)
        !
        !--annotate with time / marker legend and title
        !
        call legends_and_title
        !
        !--enter interactive mode
        !
        lastplot = ((ipos.eq.iendatstep .or. istep.eq.nsteps) .and. nyplot.eq.nyplots)

        if (interactive) then
           if (nacross*ndown.eq.1) then
              iadvance = nfreq
              call interactive_part(ntoti,iplotx,iploty,0,irenderpart,0,0, &
                   xplot(1:ntoti),yplot(1:ntoti),zplot(1:ntoti), &
                   hh(1:ntoti),icolourme(1:ntoti), &
                   xmin,xmax,ymin,ymax,rendermin,rendermax,vecmax, &
                   angletempx,angletempy,angletempz,ndim, &
                   .false.,dummy,dummy,dummy,dummy,irerender, &
                   itrackpart,icolours,iadvance,ipos,iendatstep)
              if (iadvance.eq.-666) exit over_plots
           elseif ((ipanel.eq.nacross*ndown .and. istepsonpage.eq.nstepsperpage) .or. lastplot) then
              !
              !--timestep control only if multiple plots on page
              !
              iadvance = nfreq
              call interactive_step(iadvance,ipos,iendatstep,xmin,xmax,ymin,ymax)
              if (iadvance.eq.-666) exit over_plots
           endif
        endif

     elseif (iploty.le.numplot) then! ie iploty = extra
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! additional plots (not plots of particle data - e.g where some additional 
! information is read from a file and plotted on the same page as the 
! particle plots, or where some additional plot is calculated
! from the particle data, such as errors etc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !--------------------------------------------------------------
        !  then call subroutine to plot the additional plot
        ! e.g. call routine to do convergence plot here
        !--------------------------------------------------------------

        if (iexact.eq.4 .and. iploty.eq.iacplane) then
        !
        !--A vs C for exact toystar solution
        !
           if (ndim.eq.1) then
              call exact_toystar_acplane(atstar,ctstar,sigma,gammai)
           elseif (ndim.eq.2) then
              call exact_toystar_acplane2D(atstar,ctstar,sigma,gammai)           
           endif
           !--increment page counter as setpage is not called
           iplots = iplots + 1
           ipanel = ipanel + 1
           if (ipanel.gt.nacross*ndown) ipanel = 1
        elseif (iploty.eq.ipowerspec) then 
        !
        !--power spectrum plots (uses x and data as yet unspecified)
        !
           labelx = 'frequency'
           labely = 'power'
           if (iadvance.ne.0) then
              xmin = 1./wavelengthmax  ! freq min
              xmax = 1./wavelengthmin  ! freq max
           endif
           if (iadvance.ne.0 .and. itrans(iploty).gt.0) then
              call transform_limits(xmin,xmax,itrans(iploty))
           endif       
           !
           !--setup frequency grid (evenly spaced in transformed grid)
           !
           nfreqpts = nfreqspec
           if (nfreqpts.ge.size(xplot)) then
              nfreqpts = size(xplot)
              print*,' WARNING: nfreqpts > array size, restricting to ',nfreqpts
           else
              print "(a,i6)",' number of frequency points = ',nfreqpts
           endif
           dxfreq = (xmax - xmin)/real(nfreqpts)
           do i=1,nfreqpts
              xplot(i) = xmin + (i-1)*dxfreq
           enddo
           !
           !--transform back to frequency space
           !
           if (itrans(iploty).gt.0) &
              call transform_inverse(xplot(1:nfreqpts),itrans(iploty))

           if (.not.idisordered) then! interpolate first
              !!--allocate memory for 1D grid (size = 2*npart)
              ngrid = 2*npartoftype(1)
              !!--set up 1D grid
              xmingrid = lim(ipowerspecx,1)
              xmaxgrid = lim(ipowerspecx,2)
              dxgrid = (xmaxgrid-xmingrid)/ngrid
              call set_grid1D(xmingrid,dxgrid,ngrid)

              ninterp = ntoti
              !!--interpolate to 1D grid  
              call interpolate1D(dat(1:ninterp,ipowerspecx), & 
                   pmass(1:ninterp),rho(1:ninterp), &
                   hh(1:ninterp),dat(1:ninterp,ipowerspecy), & 
                   ninterp,xmingrid,datpix1D,ngrid,dxgrid)
              !!--plot interpolated 1D data to check it
              !!print*,minval(datpix1D),maxval(datpix1D)
              !call pgswin(xmin,xmax,minval(datpix1D),maxval(datpix1D),0,1)
              !call pgbox('BCNST',0.0,0,'1BVCNST',0.0,0)      
              !call pglabel('x',label(ipowerspecy),'1D interpolation')
              !call pgline(ngrid,xgrid,datpix1D)
              !read*
              !call pgpage! change page

              !!--call power spectrum calculation on the even grid
              call powerspectrum(ngrid,xgrid,datpix1D,nfreqpts,xplot(1:nfreqpts), &
                                 yplot(1:nfreqpts),idisordered)
              if (allocated(datpix1D)) deallocate(datpix1D)
              if (allocated(xgrid)) deallocate(xgrid)             
           else
              !!--or else call power spectrum calculation on the particles themselves    
              call powerspectrum(ntoti,dat(1:ntoti,ipowerspecx), &
                   dat(1:ntoti,ipowerspecy),nfreqpts, &
                   xplot(1:nfreqpts),yplot(1:nfreqpts),idisordered)
           endif

           if (iadvance.ne.0) then
              ymin = minval(yplot(1:nfreqspec))
              ymax = maxval(yplot(1:nfreqspec))
           endif

           !!--uncomment next few lines to plot wavelengths instead
           !labelx = 'wavelength'
           !zplot(1:nfreqspec) = 1./xplot(1:nfreqspec)
           !xplot(1:nfreqspec) = zplot(1:nfreqspec)
           !if (iadvance.ne.0) then
           !   xmin = minval(xplot(1:nfreqspec))
           !   xmax = maxval(xplot(1:nfreqspec))
           !endif

           if (itrans(iploty).ne.0) then
              call transform(xplot(1:nfreqpts),itrans(iploty))
              labelx = transform_label(labelx,itrans(iploty))

              call transform(yplot(1:nfreqpts),itrans(iploty))
              labely = transform_label(labely,itrans(iploty))
              if (iadvance.ne.0) then
                 call transform_limits(xmin,xmax,itrans(iploty))
                 call transform_limits(ymin,ymax,itrans(iploty))
              endif
           endif
           
           just = 0
           title = ' '
           call page_setup

           call pgqci(icolourprev)    ! query line style and colour
           call pgqls(linestyleprev)
           if (nstepsperpage.gt.1) then
              call pgsci(linecolourthisstep) ! set appropriate colour and style if multiple steps per page
              call pgsls(linestylethisstep)
           endif
           
           call pgline(nfreqpts,xplot(1:nfreqpts),yplot(1:nfreqpts))
           print*,' maximum power at '//trim(labelx)//' = ',xplot(maxloc(yplot(1:nfreqpts)))

           call pgsci(icolourprev)
           call pgsls(linestyleprev)
           
           !
           !--redraw axes over what has been plotted
           !
           call redraw_axes(iaxis)           
           !
           !--annotate with time / marker legend and title
           !
           call legends_and_title

        endif

        lastplot = (ipos.eq.iendatstep .or. istep.eq.nsteps)

        if (interactive .and.((ipanel.eq.nacross*ndown .and. istepsonpage.eq.nstepsperpage) &
           .or. lastplot)) then
           iadvance = nfreq
           call interactive_step(iadvance,ipos,iendatstep,xmin,xmax,ymin,ymax)
           if (iadvance.eq.-666) exit over_plots
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! if plot not in correct range
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     else
        print*,' error in plotting : iplotx = ',iplotx,' iploty =',iploty, 'numplot =',numplot
        call pgpage! just skip to next plot

     endif ! ploty = whatever


  enddo over_plots ! over plots per timestep (nyplot)
  
  if (iadvance.ne.0) then
     if (allocated(datpix)) deallocate(datpix)
     if (allocated(datpix3D)) deallocate(datpix3D)
  endif
  
  return
    
contains

!----------------------------------------------
! interfaces to the page setup routines
! this is called just before a plot is
! actually plotted
!----------------------------------------------
  subroutine page_setup
    use pagesetup
    use settings_render, only:ColourBarWidth
    implicit none
    real :: barwidth, TitleOffset,xch,ych
    logical :: ipanelchange
        
    !--------------------------------------------------------------
    ! output some muff to the screen
    !--------------------------------------------------------------

    if (interactive) then
       print*,trim(labelx),' min, max = ',xmin,xmax
       print*,trim(labely),' min, max = ',ymin,ymax
    endif

    !---------------------
    ! increment counters
    !---------------------
    iplots = iplots + 1
    
    ipanelchange = .true.
    if (iplots.gt.1 .and. nyplots.eq.1 .and. nacross*ndown.gt.1.and..not.ipagechange) ipanelchange = .false.
    if (ipanelchange) ipanel = ipanel + 1
    if (ipanel.gt.nacross*ndown) ipanel = 1
    !--set counter for where we are in row, col
    icolumn = ipanel - ((ipanel-1)/nacross)*nacross
    irow = (ipanel-1)/nacross + 1

    !--------------------------------------------------------------
    ! set up pgplot page
    !--------------------------------------------------------------



    if (tile_plots) then
       !--leave space for colour bar if necessary (at end of row only)
       if (iAllowspaceforcolourbar) then
          call pgqcs(0,xch,ych)
          !--be CAREFUL changing this - there needs to be an identical version in
          !  colourbarfull
          barwidth = (ColourBarWidth + 0.25)*ych
       else
          barwidth = 0.
       endif

       inewpage = ipanel.eq.1 .and. ipanelchange .and. ipagechange
       if (inewpage) call pgpage
       call danpgtile(ipanel,nacross,ndown,xmin,xmax,ymin,ymax, &
                      trim(labelx),trim(labely),trim(title),just,iaxis,0.0,barwidth,0.0,0.0)
    else
       !--work out whether to leave space for colour bar (for each plot)
       if (iColourBar) then
          barwidth = ColourBarWidth
       else
          barwidth = 0.
       endif
        !--change the page if pagechange set
       !  or, if turned off, between plots on first page only
       inewpage = ipagechange .or. (iplots.le.nacross*ndown .and. ipanelchange)
       
       !--work out whether or not to leave space above plots for titles
       TitleOffset = 0.
       if (iPlotTitles .and. nstepsperpage.eq.1 .and. vpostitle.gt.0.) TitleOffset = vpostitle + 1.5
       
       call setpage(ipanel,nacross,ndown,xmin,xmax,ymin,ymax, &
         trim(labelx),trim(labely),trim(title), &
         just,iaxis,barwidth,TitleOffset,isamexaxis,inewpage)
    endif

    return
  end subroutine page_setup
  
!------------------------------------------------------
! draws legend(s), titles etc
! (must be called after rendering otherwise rendering
!  will overwrite plot area)
!------------------------------------------------------
  subroutine legends_and_title
    implicit none

    !--plot time on plot
    if (iPlotLegend .and. nyplot.eq.1) call legend(legendtext,timei,labeltimeunits, &
                                                   hposlegend,vposlegend,fjustlegend)
    !--line/marker style/colour legend for multiple timesteps on same page
    if (iPlotStepLegend .and. nyplot.eq.1 .and. istepsonpage.gt.0) then
       if (iploty.gt.ndataplots) then
          call legend_markers(istepsonpage,linecolourthisstep,imarktype(1),linestylethisstep, &
            .false.,.true.,trim(steptitles(istepsonpage)),hposlegend,vposlegend)           
       else
          call legend_markers(istepsonpage,linecolourthisstep,imarktype(1),linestylethisstep, &
            iplotpartoftype(1),iplotline,trim(steptitles(istepsonpage)),hposlegend,vposlegend)
       endif
    endif
    !--print title if appropriate
    if (iPlotTitles .and. istepsonpage.eq.1 .and. ipanel.le.ntitles) then
       if (len_trim(pagetitles(ipanel)).gt.0) then
          call pgmtxt('T',vpostitle,hpostitle,fjusttitle,trim(pagetitles(ipanel)))
       endif
    endif
    
    return
  end subroutine legends_and_title

!--------------------------------------------
! sets up a one dimensional grid of pixels
! and allocates memory for datpix1D
!--------------------------------------------
  subroutine set_grid1D(xmin1D,dxgrid1D,ngridpts)
    implicit none
    integer, intent(in) :: ngridpts
    real, intent(in) :: xmin1D, dxgrid1D
    integer :: igrid

    if (allocated(datpix1D)) deallocate(datpix1D)
    if (allocated(xgrid)) deallocate(xgrid)
    allocate (datpix1D(ngridpts))
    allocate (xgrid(ngridpts))

    do igrid = 1,ngridpts
       xgrid(igrid) = xmin1D + igrid*dxgrid1D - 0.5*dxgrid1D
    enddo

  end subroutine set_grid1D

!-------------------------------------------------------------------
! interface to vector plotting routines
! so that pixel arrays are allocated appropriately
!-------------------------------------------------------------------
  subroutine vector_plot(ivecx,ivecy,numpixx,numpixy,pixwidth,vmax,label)
   use fieldlines
   use settings_vecplot, only:UseBackgndColorVecplot
   use interpolations2D, only:interpolate2D_vec
   use projections3D, only:interpolate3D_proj_vec
   use render, only:render_vec
   implicit none
   integer, intent(in) :: ivecx,ivecy,numpixx,numpixy
   real, intent(in) :: pixwidth
   real, intent(inout) :: vmax
   character(len=*), intent(in) :: label
   real, dimension(numpixx,numpixy) :: vecpixx, vecpixy

   !print*,'plotting vector field ',trim(label)
   if ((ivecx.le.ndim).or.(ivecx.gt.ndataplots) &
        .or.(ivecy.le.ndim).or.(ivecy.gt.ndataplots)) then
      print*,'error finding location of vector plot in array'
   else
      !!--plot arrows in either background or foreground colour
      if (UseBackgndColorVecplot) call pgsci(0)

      !
      !--interpolate using appropriate routine for number of dimensions
      !
      select case(ndim)
      case(3)
         if (x_sec) then ! take vector plot in cross section
            call interpolate3D_xsec_vec(xplot(1:ninterp), &
              yplot(1:ninterp),zplot(1:ninterp), &
              hh(1:ninterp),weight(1:ninterp), &
              dat(1:ninterp,ivecx),dat(1:ninterp,ivecy), &
              ninterp,xmin,ymin,zslicepos, &
              vecpixx,vecpixy,numpixx,numpixy,pixwidth,inormalise)
         else
            call interpolate3D_proj_vec(xplot(1:ninterp), &
              yplot(1:ninterp),zplot(1:ninterp),hh(1:ninterp), &
              weight(1:ninterp),dat(1:ninterp,ivecx),dat(1:ninterp,ivecy), &
              ninterp,xmin,ymin, &
              vecpixx,vecpixy,numpixx,numpixy,pixwidth,dobserver,dscreenfromobserver)
         endif
      case(2)
         !
         !--or interpolate (via averaging) to coarser grid
         !
         !call fieldlines2D(ninterp,xplot(1:ninterp),yplot(1:ninterp), &
         !     dat(1:ninterp,ivecx),dat(1:ninterp,ivecy), &
         !     hh(1:ninterp),pmass(1:ninterp), &
         !     rho(1:ninterp),xmin,xmax,ymin,ymax)
         !call interpolate_vec(xplot(1:ninterp),yplot(1:ninterp), &
         !  dat(1:ninterp,ivecx),dat(1:ninterp,ivecy), &
         !  xmin,ymin,pixwidth,vecpixx,vecpixy, &
         !  ninterp,numpixx,numpixy)
         
         call interpolate2D_vec(xplot(1:ninterp),yplot(1:ninterp), &
              hh(1:ninterp),weight(1:ninterp),dat(1:ninterp,ivecx), &
              dat(1:ninterp,ivecy),ninterp,xmin,ymin, &
              vecpixx,vecpixy,numpixx,numpixy,pixwidth,inormalise)
      
      case default
         print "(a,i1,a)",'ERROR: Cannot do vector plotting in ',ndim,' dimensions'
         return
      end select
      !
      !--plot it
      !
      call render_vec(vecpixx,vecpixy,vmax, &
           numpixx,numpixy,xmin,ymin,pixwidth,trim(label),' ')

      if (UseBackgndColorVecplot) call pgsci(1)

   endif
  
  end subroutine vector_plot

end subroutine plotstep

end module timestep_plotting
