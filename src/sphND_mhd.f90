!!---------------------------------------------------------------------!!
!!                                                    _         _      !!
!!      ___ _   _ _ __   ___ _ __ ___ _ __  _ __ ___ | |__   __| |     !!
!!     / __| | | | '_ \ / _ \ '__/ __| '_ \| '_ ` _ \| '_ \ / _` |     !!
!!     \__ \ |_| | |_) |  __/ |  \__ \ |_) | | | | | | | | | (_| |     !!
!!     |___/\__,_| .__/ \___|_|  |___/ .__/|_| |_| |_|_| |_|\__,_|     !!
!!               |_|                 |_|                               !!
!!                                                                     !!
!!       _   _     _   _   _   _   _   _     _   _   _   _   _         !!
!!      / \ / \   / \ / \ / \ / \ / \ / \   / \ / \ / \ / \ / \        !!
!!     ( B | y ) ( D | a | n | i | e | l ) ( P | r | i | c | e )       !!
!!     \_/ \_/   \_/ \_/ \_/ \_/ \_/ \_/   \_/ \_/ \_/ \_/ \_/         !!
!!								       !!
!!---------------------------------------------------------------------!!
!! An N-D SPH code to handle compressible gas dynamics with MHD        !!
!!                                                                     !!
!! Written in Fortran 90                                               !!
!! By Daniel Price, Institute of Astronomy, Cambridge, UK, 2002-2003   !!
!! Email: dprice@ast.cam.ac.uk				               !!
!!								       !!
!! This version is designed to be as modular (and thus as adaptable)   !!
!! as possible, as a testbed for SPH algorithms		               !!
!!                                                                     !!
!! Specific features include:                                          !!
!!                                                                     !!
!!  * options to choose continuity equation or density by              !!
!!    direct summation                                                 !!
!!                                                                     !!
!!  * options to choose different energy eqn implementations           !!
!!                                                                     !!
!!  * different artificial viscosity prescriptions                     !!
!!                                                                     !!
!!  * choice of whether to use average h, average gradient             !!
!!    of kernel or Springel/Hernquist (2001) type correction           !!
!!    terms for varying h		                               !!
!!                                                                     !!
!!  * Morris and Monaghan (1997) artificial viscosity switch           !!
!!    (turns off artificial viscosity away from shocks)                !!
!!                                                                     !!
!!  * ghost particle boundaries (reflective/periodic)                  !!
!!                                                                     !!
!! Although this code is original, parts of it are inspired            !!
!! by other SPH codes written by Joe Monaghan and Matthew Bate. These  !!
!! bits are (hopefully) acknowledged.                                  !!
!!                                                                     !!
!!---------------------------------------------------------------------!!

PROGRAM SUPERSPMHD_ND
 USE dimen_mhd
 USE debug
 USE loguns
 USE options
 USE timestep
 USE versn
!
!--define local variables
!
 IMPLICIT NONE
 REAL :: tprint
 INTEGER :: noutput,nevwrite
!
!--version number
!
   version = 'SPMHD-ND-v3.7'
!
!   version = 'SPMHD-ND-v3.6_22_10_2003'
!   * bug fix in thermal conductivity in thermal energy equation
!   * bug fix in source term for artificial stress in total energy equation
!   * something goes wrong with mshk2 when hfact = 1.3 (art. stress?)
!   * use inflow/outflow boundaries in 1D (boundaryND_1D).
!   * bug fix in this for outflow boundaries (linkND)
!   * another bug fix in inflow/outflow (copies forces etc when reshuffling)
!   * setkern not called with ikernel (set inside setkern itself)
!   * new header
!   * uses dhdt in step/rates
!   * uniform setups work in 1D
!   * setup_wave2D, setup_wave_x_ND, initialise calls ghosts properly
!   *** this is the version used for 1D results in revised papers ***
!   version = 'SPMHD-ND-v3.5_10_10_2003'
!   * timestep criterion agrees with paper I (and version SPMHD-f90-v2.5)
!   * found the bug in iener=3!!! (rates)
!   * fixed offset bug in set_uniform_cartesian
!   * bug with -ve densities in unused ghosts
!   * changes to sametime/multirun
!   version = 'SPMHD-ND-v3.4_07_08_2003'
!   * interpolate_kernel function - setkern,density,rates,density_partial
!   * minor changes to step
!   * no explicit references to idim - bugs with this fixed (nlistdim)
!   * smoothing in blastwave, random particle setups
!   * smooth.f90 -> SPH smoothing of initial conditions
!   * eos can do either individual elements or whole arrays
!   * experimenting with namelist input
!   * gravity by direct summation
!   * multirun compiles with SPH source files
!   * set_uniform_spherical
!   * setkern cleaned up
!   * polytrope setup
!   * BUG IN LINK FOR NO BOUNDARIES
!   * beginnings of maketree
!   * trace and idebug not itrace
!   version = 'SPMHD-ND-v3.3_22_07_2003'
!   * density iteration - only iterates particles which are not converged
!     (density_partial.f90, iterate_density.f90) something wrong however
!   * toy star: setup using different particle masses - force bug fixed
!   * initialise: no longer reads file runname, opens .dat,.ev earlier
!   * hconst removed, different particle masses didn't work
!   * toy star solution in supersphplot
!   * toy star results obtained using this code
!   * grad h terms bug fixed in ND
!   * setup for 2D MHD blast wave
!   * log plots in supersphplot
!   version = 'SPMHD-ND-v3.2-25-06-03'
!   * rendering for supersphplot
!   * rates cleaned up and split into separate subroutines
!   * bug fix in MHD dissipation terms, also in boundaryND for periodic bc's
!   * initialise can read rootname off command line or from file
!   * optimization of density/rates - divisions done only once
!   * catches vsig error
!   * write header in two parts
!   * bug in density iteration fixed (+ Newton Raphson)
!   * toy star setup 1-x^2
!   * new Makefiles - make .tar.gz's, saves versions, compiles in 1 or 2D
!   * processplot.pl to script movie generation
!   * supersphplot split into subroutines, new Makefile
!   * NB: there is *something* wrong with the iner=3 results for 1D shocks
!   version = 'SPMHD-ND-v3.1-16-05-03'
!   * uniform density setups in ND (setup_unifdis)
!   * 2D shock setup
!   * 1D shock setup for ND code - bug fix in get_neighbour_list for 1D
!   * ghost particle setup rehashed in ND + bug fixed (dx > 0 only now)
!   * supersphplot changes
!   * some changes to linkND
!   * rates doesn't calculate interactions between ghost particles
!   * orszag-tang vortex seems to work OK
!   * bug fix in Joe's anticlumping term (2D) - uses 1/hfact, define hconst
!  version = 'SPMHD-ND-v3.0-29-04-03'
!   * ND (not fully tested)
!   * setup for Orszag-Tang vortex in 2D  
!   * ghost particles setup is simpler - must call set_ghost_particles
!     before the call to link (link finds ghosts rather than ghosts creating
!     new link list cells)
!   * allocate reallocates arrays if particle # > array size
!   * ghost particles work in 2D, need to sort out corners in 3D
!   * bug in prvaniso with anticlumping term
!   * need to fix periodic BC's, 1D problems
! version = 'SPMHD-f90-v2.5-lots-of-options'
!   * isothermal MHD
! version = 'SPMHD-f90-v2.4-lots-of-options'
!   * allocation of arrays is done in allocate.f90
!   * link lists can handle empty cells
!   * bug in link fixed (ibound=1)
!   * bug in density/rates when only one particle in a link list cell fixed
!   * inflow/outflow boundaries bug fixed (Srojeen test + MHD shocks)
!   * pressure on fixed particles near edges fixed
!   ** this version used to obtain 1D results with inflow/outflow
! version = 'SPMHD-f90-v2.3-lots-of-options'
!   * inflow/outflow boundaries work in 1D
!   * wrinsph writes an infile. 
!   * batch running - domulti.pl, multirun.f90, sametime.f90
!   ** this is the version used to obtain all the 
!     results in the 1D case
! version = 'SPMHD-f90-v2.2-lots-of-options'
!   * time step condition uses correct signal velocity (previously too slow)
!   * energy equation includes grad h terms and Joe's anti-clumping term
!   * initialise: ghost particles set to agree initially (t=0 graphs look better)
!   * equation_of_state no longer calculates alfven speed (unnecessary)
!   * limits/plots improved in supersphplot
!   * bug in periodic boundaries fixed (causing spurious currents)
!   * fast/slow wave setups properly done
!   * inflow/outflow boundary conditions (boundary1D) - not tested
! version = 'SPMHD-f90-v2.1-lots-of-options'
!   * MHD works!! Uses AV in direction of vunit for 1.5D problems
!   * Joes anticlumping term fixed - now works as it should
!   * all setups in f90 
! version = 'SPMHD-f90-v2.0-lots-of-options'
!   * fortran 90 (uses modules, allocatable arrays)
! version = 'SPMHD-v1.6-lots-of-options'
!   * follows Joe's notes - viscosity term in induction equation
!   * bug fix in ghost1D - reflects all velocities, not just vx.
! version = 'SPMHD-v1.5-lots-of-options'
!   * setup for advection test(ok) and ryu/jones shock tubes (still crap)
!   * spits out current density
!   * evolves h in evolution equations (if ihvar = 2)
!   * Joes clumping correction only in x-direction (line between particles)
!   * fixed bug in periodic boundary conditions (link1D/step1D)
! version = 'SPMHD-v1.4-lots-of-options'
!   * reverted to old h update (blast wave didn't work) 
!   * can evolve using either B/rho or B 
!    (setup is now always just the basic field B)
!   * bug check in read_infile
!   * neighbour array size checked in get_neigh...(set in COMMONS/linklist)
!   * now uses correct signal speed!
!   * anisotropic forces with vectors
! version = 'SPMHD-v1.3-lots-of-options'
!   * ensures h proportional to 1/rho
!   * bug in iav = 2 fixed
!   * uses COMMONS/dimen_mhd separate to SUPERSPH1D
!   * alternative induction equations (incl. v*div B term)
!   * more mag force options (6=Morris' hybrid force)
!   * alternative kernel is used only when idivBzero=2
!   * iener=3 seems to work with MHD
! version = 'SPMHD-v1.2-lots-of-options' 
!   * kernel table changed to try the thomas/couchman anti-clumping kernel
!    (stores grkern instead of grkern/(r/h) ). 
!   * changed magnetic field variable to be of any dimension e.g. Brho(3,idim)
!   * fixed bug in adjust_bound and in sum_density. 
!   * Velocity variable also changed to be of arbitrary dimension.
!   * fixed bug in read_infile (didn't read hfact)
!   * however there is still a bug in iav = 2
! version = 'SPMHD-v1.1-lots-of-options' 
!   * polytropic equation of state
!   * prints out loads of MHD parameters (use evsupersph.f)
!   * setup for Brio/Wu shock tubes
!   * MHD for iener = 3 in step
!   * Joe's clumping correction term (doesn't fix things though)
! version = 'SPMHD-v1.0-lots-of-options' 
!   * MHD
! version = 'SPH-v1.1-lots-of-options'
!   * individual common blocks 
 
 trace = .false.		! set tracing flow (prints entry into subroutine)
! trace = .true.
! idebug = 'density'
! idebug = 'none'
! idebug = 'link'
 
! itemp = 40000		! debug one particular particle

 CALL logun     	! set logical unit numbers to use for input/output
 CALL initialise	! read files and parameters and setup particles
!
!--Set initial timestep
!
 dt = 0.
 tprint = 0.0
 time = 0.0
 nsteps = 0
 nevwrite = 1	! frequency of writing to .ev file (could be read as parameter)
!
!--write the initial conditions to the output and evolution files
!
 CALL output(time,nsteps)
 CALL evwrite(time)
 noutput = 1
 tprint = tout
! CALL quit
!
! --------------------- Main loop ----------------------------------------
!
 dostep: DO WHILE ((time.LT.tmax).AND.(nsteps.LT.nmax))

    hdt = 0.5*dt
    time = time + dt
    nsteps = nsteps + 1

    CALL step 	 		!  Evolve data for one timestep

    dt = min(dtforce,dtcourant) ! new timestep from force/Courant condition
!
!--write log every step in 2D/3D
!
    IF (ndim.GE.2) THEN
       WRITE(iprint,10) time,dtforce,dtcourant
10     FORMAT(' t = ',f9.4,' dtforce = ',1pe10.3,' dtcourant = ',1pe10.3)
    ENDIF
    
    IF (dt.LT.1e-8) THEN
       WRITE(iprint,*) 'main loop: timestep too small, dt = ',dt
       CALL quit
    ENDIF
!
!--Write to data file if time is right
! 
    IF (   (time.GE.tprint) 				&
       .OR.(time.GE.tmax)				&
       .OR.((MOD(nsteps,nout).EQ.0).AND.(nout.GT.0))	&
       .OR.(nsteps.GE.nmax)  ) THEN
 
!--if making movies and need ghosts to look right uncomment the line below, 
       IF (idumpghost.EQ.1) CALL set_ghost_particles
       CALL output(time,nsteps)
       noutput = noutput + 1
       tprint = noutput*tout
    ENDIF
!    
!--calculate total energy etc and write to ev file    
!
    IF (MOD(nsteps,nevwrite).EQ.0) CALL evwrite(time)

    IF (dt.GE.(tprint-time)) dt = tprint-time	! reach tprint exactly
	 	 
 ENDDO dostep

!------------------------------------------------------------------------

!
!--close all open files and exit
!
 IF (iprint.NE.6) THEN
    CLOSE(UNIT=iprint)
 ENDIF
 CLOSE(UNIT=ievfile)
 CLOSE(UNIT=idatfile)
 STOP 
 
END PROGRAM

!!---------------------------------------------------------------------
!! This subroutine performs a graceful exit if the program has crashed
!!---------------------------------------------------------------------

SUBROUTINE quit
 USE loguns
 USE timestep
 IMPLICIT NONE

 WRITE(iprint,*) 'performing graceful exit...'
 CALL output(time,nsteps) ! dump particles before crashing
!
!--Close all open files and exit
!
 IF (iprint.NE.6) THEN
    CLOSE(UNIT=iprint)
 ENDIF
 CLOSE(UNIT=ievfile)
 CLOSE(UNIT=idatfile)
 STOP
 
END
