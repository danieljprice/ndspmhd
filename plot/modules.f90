!----------------------------------------------------------------------------
!
!  modules containing global variables
!
!----------------------------------------------------------------------------
!
!--global parameters  (should really allocate memory appropriately)
!
module params
 implicit none
 integer, parameter :: maxplot=40   ! maximum number of plots (for multiplot arrays)
end module params
!
!--particle data
!
module particle_data
 use params
 implicit none
 integer :: maxpart,maxstep,maxcol  ! dimensions of dat array
 integer, allocatable, dimension(:) :: npart,ntot,nghost,ntotplot
 integer, allocatable, dimension(:,:) :: iam
 real, allocatable, dimension(:) :: time, gamma
 real, allocatable, dimension(:,:,:) :: dat
 real :: hfact 
end module particle_data
!
!--filename
!
module filenames
 implicit none
 integer, parameter :: maxfile = 10
 integer :: ifile,nfilesteps,nfiles
 character(len=20), dimension(maxfile) :: rootname
end module filenames
!
!--labels for all plots and the locations of certain useful variables
!
module labels
 use params
 implicit none
 character(len=20), dimension(maxplot+2) :: label
 character(len=7), dimension(3,3) :: labelcoord
 integer, dimension(3) :: ix
 integer :: ivx,ivlast,irho,iutherm,ipr,ih,irad,ibfirst,iblast
 integer :: ipmass
 integer :: ientrop,ipmag,ibeta,itotpr,ike,idivb,idivberr,iJfirst
 integer :: iacplane,itimestep,ipowerspec
 integer :: irad2,ivpar,ivperp,iBpar,iBperp
end module labels
!
!--plot limits
!
module limits
 use params
 implicit none
 real, dimension(maxplot,2) :: lim
 real :: bmin,bmax 
end module limits
!
!--module containing plot settings
!
module settings
 use params
 implicit none
 integer :: numplot,ncalc,ncolumns,nextra,nacross,ndown
 integer :: ndataplots
 integer :: ndim, ndimv 
 integer :: imark, imarkg, imarksink
 integer :: ixsec,nxsec, nbins,nc
 integer :: linestylein, iexact
 integer :: irenderplot
 integer :: ncolours,nstart,n_end,nfreq
 integer :: icoords,icoordsnew,iaxis
 integer :: ncircpart, itrackpart
 integer, dimension(10) :: icircpart
 
 integer :: ncontours_nomulti,npix_nomulti,npixvec_nomulti
 integer :: ivecplot_nomulti,icolours
 integer :: ipapersize
 
 real :: scalemax,zoom,hposlegend,vposlegend
 real, dimension(3) :: xminoffset_track, xmaxoffset_track
 real :: xsecpos_nomulti
 real :: papersizex,aspectratio
!--plot options
 logical :: animate, interactive
 logical :: iadapt,ihavereadfilename
 logical :: plotcirc,plotcircall,flythru,imulti,ipagechange
 logical :: iplotline,iplotlinein,iplotav,ilabelpart
 logical :: iplotpart,iplotghost,iplotsink
 logical :: ishowopts, ivegotdata, tile
 
 logical :: backgnd_vec_nomulti
 logical :: iplotcont_nomulti,xsec_nomulti,iplotpartvec_nomulti
!
!--power spectrum options
!
 integer :: ipowerspecy, nfreqspec
 real :: wavelengthmax
 logical :: idisordered

!
!--sort these into a namelist for input/output
!
 namelist /plotopts/ iaxis,tile, &
   animate,iadapt,xsec_nomulti,flythru, &
   plotcirc,iplotline,iplotlinein,linestylein,          &
   imark, imarkg, imarksink,                            &
   nacross,ndown,                                       &
   iexact,iplotav,nbins,                                &
   ivecplot_nomulti,iplotpartvec_nomulti,       &
   npix_nomulti,npixvec_nomulti,                        &
   iplotcont_nomulti,ncontours_nomulti,                 &
   icolours,iplotghost,iplotsink,                       &
   ipapersize,papersizex,aspectratio,                   &
   ipowerspecy,idisordered,wavelengthmax,nfreqspec,icoordsnew, &
   ncircpart,icircpart,hposlegend,vposlegend
     
end module settings
!
!--multiplot settings
!
module multiplot
 use params
 implicit none
 integer :: nyplotmulti 
 integer, dimension(maxplot) :: multiplotx,multiploty
 integer, dimension(maxplot) :: irendermulti,ivecplotmulti
 integer, dimension(maxplot) :: npixmulti,npixvecmulti 
 integer, dimension(maxplot) :: ncontoursmulti,itrans
 logical, dimension(maxplot) :: iplotcontmulti, iplotpartvecmulti     
 logical, dimension(maxplot) :: x_secmulti,backgnd_vec_multi
 real, dimension(maxplot) :: xsecposmulti
!
!--sort these into a namelist for input/output
!
 namelist /multi/ nyplotmulti,                                  &
    itrans,multiplotx,multiploty,irendermulti,                  &
    iplotcontmulti,ncontoursmulti,ivecplotmulti,npixmulti,      &
    npixvecmulti,iplotpartvecmulti,x_secmulti,xsecposmulti
 
end module multiplot
!
!--exact solution parameters
!
module exact_params
 implicit none
!--toy star
 integer :: norder ! for toy star
 real :: htstar,atstar,ctstar,totmass,sigma,sigma0
!--sound wave
 integer :: iwaveplot ! linear wave
 real :: ampl,lambda,period
!--sedov blast wave
 real :: rhosedov,esedov
!--polytrope
 real :: polyk
!--mhd shock solutions
 integer :: ishk
!--from file
 integer, parameter :: maxexactpts = 1001
 integer :: iexactpts, iexactplotx, iexactploty
 real, dimension(maxexactpts) :: xexact,yexact
!--shock tube
 real :: rho_L, rho_R, pr_L, pr_R, v_L, v_R
!
!--sort these into a namelist for input/output
!
 namelist /exactparams/ ampl,lambda,period,iwaveplot, &
          htstar,atstar,ctstar,sigma0,norder,rhosedov,esedov, &
	  rho_L, rho_R, pr_L, pr_R, v_L, v_R, &
	  iexactplotx,iexactploty 

end module exact_params
!
!--tabulated column density through the kernel 
!  (used in interpolate3d_projection)
!
module column
 implicit none
 integer, parameter :: maxcoltable = 1000
 real, dimension(maxcoltable) :: coltable
end module column
