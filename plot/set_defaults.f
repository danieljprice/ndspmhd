!!
!!     set initial default options
!!      these are used if no defaults file is found
!!
      SUBROUTINE set_defaults
      USE exact_params
      USE labels
      USE multiplot
      USE settings
      IMPLICIT NONE
!
!--set default options
!
      numplot=maxplot 	! reset if read from file
      ncalc = 0		! number of columns to calculate(e.g. radius)
      nextra = 0	! extra plots aside from particle data
      iextra = 0	! label position of extra plot
      ncolumns=maxplot-ncalc	! number of columns in data file
      ndim = ndimmax		! number of coordinate dimensions
      ndimV = ndim	! default velocity same dim as coords
      nstart = 1	! timestep to start from
      n_end = maxstep	! timestep to finish on
      nfreq = 1		! frequency of timesteps to read
      animate = .true.	! turns off/on prompt between page changes
      magfield = .true.	! historical - can be used to set whether MHD or not
      iadapt = .true.	! adaptive plot limits
      plotcirc = .false.	! plot circle of radius 2h around particles
      plotcircall = .false.	!  " " around all particle
      icircpart = 1		!  " " around a specific particle
      xsec_nomulti = .false.		! take cross section of data / particles
      flythru = .false.		! take series of cross sections through data
      ipagechange = .true.	! if false plots graphs on top of each other
      scalemax = 1.0	! for rescaling adaptive limits
      zoom = 1.0	! for rescaling fixed limits
      imark = 1	! PGPLOT marker for particles
      imarkg = 4	! PGPLOT marker for ghost particles
      imarksink = 17	! PGPLOT marker for sink particles 
      nacross = 1	! number of plots across page
      ndown = 1		! number of plots down page
      ipapersize = 0	! paper size option
      papersizex = 0.0	! size of x paper (no call to PGPAP if zero)
      aspectratio = 0.0	! aspect ratio of paper (no call to PGPAP if zero)
      iplotline = .false.	! plot line joining the particles
      iplotlinein = .false.	! " " but on first step only
      linestylein = 4		! PGPLOT line style for above
      iexact = 0		! exact solution to plot
      iplotav = .false.		! plot average line through particles
      nbins = 24		! number of bins for this
      ilabelpart = .false.	! plot particle numbers
      iplotpart = .true.	! flag whether or not to plot actual SPH particles
      iplotpartvec_nomulti = .true.	! whether to plot particles on vector plot
      iplotghost = .true.	! plot ghost particles
      iplotsink = .true.	! plot sink particles
      npix_nomulti = 20		! pixels in x direction for rendering
      npixvec_nomulti = 20	! pixels in x direction on vector plots
      ivecplot_nomulti = 0	! choice of vector plot
      irender = 0	! choice of rendering plot
      iplotcont_nomulti = .true.	! plot contours
      ncontours_nomulti = 30		! number of contours to plot
      icolours = 0		! colour scheme to use
      ncolours=10		! number of colours in colour table
      itrans(:) = 0		! no transformations (log10 etc)
      backgnd_vec_nomulti = .false. ! plot vector plot using black/white
      
      lambda = 1.0	! sound wave exact solution : wavelenght
      delta = 0.005	! sound wave exact solution : amplitude

      nyplotmulti = 1		! number of plots in multiplot
      multiploty(1) = ndim+1 	! first plot : y axis
      multiplotx(1) = 1		! first plot : x axis
      irendermulti(1) = 0	! first plot : rendering
      ivecplotmulti(1) = 0	! first plot : vector plot
      x_secmulti(1) = .true.	! first plot : take cross section
      xsecposmulti(1) = 0.0	! first plot : position of cross section
      
      RETURN    
      END SUBROUTINE set_defaults
