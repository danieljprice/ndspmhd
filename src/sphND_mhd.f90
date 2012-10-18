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
!!                                                                     !!
!!---------------------------------------------------------------------!!
!! An N-D SPH code to handle compressible gas dynamics with MHD        !!
!!                                                                     !!
!! Written in Fortran 90                                               !!
!! By Daniel Price, Institute of Astronomy, Cambridge, UK, 2002-2004   !!
!! Email: dprice@ast.cam.ac.uk                                         !!
!!                                                                     !!
!! This version is designed to be as modular (and thus as adaptable)   !!
!! as possible, as a testbed for SPH algorithms                        !!
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
!!    terms for varying h                                              !!
!!                                                                     !!
!!  * Morris and Monaghan (1997) artificial viscosity switch           !!
!!    (turns off artificial viscosity away from shocks)                !!
!!                                                                     !!
!!  * ghost particle boundaries (reflective/periodic)                  !!
!!                                                                     !!
!! Although this code is original, I learnt my SPH from                !!
!! other SPH codes written by Joe Monaghan and Matthew Bate, and so    !!
!! some parts bear similarities to these codes.                        !!
!!                                                                     !!
!!---------------------------------------------------------------------!!

PROGRAM SUPERSPMHD_ND
 USE debug
 USE loguns
 USE versn
 IMPLICIT NONE
 INTEGER, PARAMETER :: maxruns = 200
 INTEGER :: i,iprev,irun, nruns
 CHARACTER(LEN=120), DIMENSION(maxruns) :: runname
 runname(1) = 'a'
!
!--version number
!
    version = 'NDSPMHD-3D-SR-v5-4'
!   * special relativity
!   * a LOT of other changes over the years
!    version = 'NDSPMHD-3D-v5-3'
!   * multiple runnames off command line (does them in order)
!   * dissipation switches for resistivity and conductivity
!   * anticlumping term implemented as a modified kernel gradient
!   * gradh terms calculated for anticlumping kernel (**removed)
!   * particle info explicitly passed to density, density_partial
!   * external_forces subroutine cleaned up + potential calculation
!    version = 'NDSPMHD-3D-v5-2_18_05_2004' (saved 2:50pm)
!    *** this version used to obtain swave and hydro shocks for thesis 10/5/04 ***
!    *** this version used to obtain 1D MHD shock tube results for thesis 13/5/04
!    *** also used for Bx peak test in thesis/paper III 17/5/04
!   * 1D turbulence setup revised/set_powerspec
!   * fixed particle boundaries work in >1D, initialise cleaned up
!   * hyperbolic/parabolic divergence cleaning
!   * div B = 0 by projection method in 2D - poisson eq. by direct sum in 2D,3D
!   * accretion disc setups
!   * step (predictor-corrector) could be bad for disks - look out!
!   * lots of crap to do with MHD AV/switches in rates
!   * lots of things for GR code (grutils, conservative2primitive_gr etc)
!   * tested on shear flows
!   * preliminary riemann solver routine
!   * compiler warnings fixed in lots of subroutines (mostly unused variables)
!     -> link renamed so doesn't conflict with internal function
!   * kernel plotting utility, more kernels added
!   * quite a few changes to supersphplot
!   * major revamp of rates - split into separate subroutines (faster)
!    version = 'NDSPMHD-3D-v5-1'
!   * smoothing length iteration on single particles works 
!   * conservative2primitive and primitive2conservative
!   * eos rehashed
!   * setup is on primitive variables  
!   * psep sent to set_uniform_cartesian
!   * external forces subroutine and option instead of itoystar
!   * main loop in subroutine evolve
!   * potential energy from external forces
!   * symplectic integrator (not yet for MHD, not as good as leapfrog)
!   * ihvar = 3 predictor step only in h update
!   * setup_unifsph
!    version = 'NDSPMHD-3D-v5-0'
!   *** saved as working 3D version ***
!   * equations use general alternative formulation
!   * compiles in 3D
!   * set_ghosts totally rewritten -> works in up to 3D
!   * 3D neighbour finding fixed - density calculated OK in 3D
!   * fixed particle boundaries in > 1D (itype)
!   * handles zero pressure
!   * fix for negative thermal energies
!   * my symmetrisation of vsig - but not for MHD yet
!   * can have different boundary options in different dimensions
!   * set_uniform cartesian can be called multiple times
!   * shock setups in 3D
!   version = 'NDSPMHD-v4-0'        ! make sure there are no .'s in version name
!   *** versioning now done with CVS ***
!   use 'make tag' to tag a working copy of the code in CVS
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
!   * rates cleaned up and tried to split into separate subroutines
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
 
 trace = .false.                ! set tracing flow (prints entry into subroutine)
! trace = .true.
 idebug = 'none' 
! idebug = 'fixed'
! idebug = 'density'
! idebug = 'neighb'
! idebug = 'link'
! idebug = 'divB'
  
! itemp = 40000                ! debug one particular particle

 CALL logun     ! set logical unit numbers to use for input/output
!
!--get runname(s) off command line
!  
 iprev = 1
 i = 1
 do while (runname(iprev)(1:1).ne.' ' .and. i.le.maxruns)
    call getarg(i,runname(i))
    !!print*,i,runname(i)
    iprev = i
    i = i + 1
    if (i.gt.maxruns .and. runname(iprev)(1:1).ne.' ') then
       print*,'WARNING: number of runs >= array size: setting nruns = ',maxruns
    endif
 enddo
 if (i.gt.maxruns .and. runname(maxruns)(1:1).ne.' ') then
    nruns = maxruns
 else
    nruns = iprev - 1
 endif
 print*,'number of runs = ',nruns
! 
!--If nothing on command exit and print usage
!
 IF (runname(1)(1:1).EQ.' ') THEN
    nruns = 1
10  WRITE(6,*) 'Enter name of run:'
    READ(*,*,ERR=10) runname(1)
!    STOP 'Usage: spmhd runname '
 ENDIF

!
!--for each runname, perform a simulation
!
 DO irun = 1,nruns
    rootname = runname(irun)
    print*,'run = ',irun,' runname = ',rootname
 
    CALL initialise  ! read files and parameters and setup particles
    !
    !--now call the main timestepping loop
    !
    CALL evolve
    !
    !--close all open files and exit
    !
    IF (iprint.NE.6) THEN
       CLOSE(UNIT=iprint)
    ENDIF
    CLOSE(UNIT=ievfile)
    CLOSE(UNIT=idatfile)
 ENDDO

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
 CALL output(time,-nsteps) ! dump particles before crashing
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
