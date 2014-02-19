!------------------------------------------------------------------------------!
! NDSPMHD: A Smoothed Particle (Magneto)Hydrodynamics code for (astrophysical) !
! fluid dynamics simulations in 1, 2 and 3 spatial dimensions.                 !
!                                                                              !
! (c) 2002-2014 Daniel Price                                                   !
!                                                                              !
! http://users.monash.edu.au/~dprice/ndspmhd                                   !
! daniel.price@monash.edu -or- dprice@cantab.net (forwards to current address) !
!                                                                              !
!  NDSPMHD comes with ABSOLUTELY NO WARRANTY.                                  !
!  This is free software; and you are welcome to redistribute                  !
!  it under the terms of the GNU General Public License                        !
!  (see LICENSE file for details) and the provision that                       !
!  this notice remains intact. If you modify this file, please                 !
!  note section 2a) of the GPLv2 states that:                                  !
!                                                                              !
!  a) You must cause the modified files to carry prominent notices             !
!     stating that you changed the files and the date of any change.           !
!                                                                              !
!  ChangeLog: see below                                                        !
!------------------------------------------------------------------------------!

!!---------------------------------------------------------------------!!
!!                     _                     _         _               !!
!!           _ __   __| |___ _ __  _ __ ___ | |__   __| |              !!
!!          | '_ \ / _` / __| '_ \| '_ ` _ \| '_ \ / _` |              !!
!!          | | | | (_| \__ \ |_) | | | | | | | | | (_| |              !!
!!          |_| |_|\__,_|___/ .__/|_| |_| |_|_| |_|\__,_|              !!
!!                          |_|                                        !!
!!        _   _     _   _   _   _   _   _     _   _   _   _   _        !!
!!       / \ / \   / \ / \ / \ / \ / \ / \   / \ / \ / \ / \ / \       !!
!!      ( B | y ) ( D | a | n | i | e | l ) ( P | r | i | c | e )      !!
!!       \_/ \_/   \_/ \_/ \_/ \_/ \_/ \_/   \_/ \_/ \_/ \_/ \_/       !!
!!                                                                     !!
!!---------------------------------------------------------------------!!
!! An N-D SPH code to handle compressible gas dynamics with MHD        !!
!!                                                                     !!
!! Written in Fortran 90                                               !!
!!                                                                     !!
!! By Daniel Price, Institute of Astronomy, Cambridge, UK, 2002-2004   !!
!!                  University of Exeter, UK 2004-2008                 !!
!!                  Monash University, Melbourne, Australia 2008-      !!
!!                                                                     !!
!! Email: daniel.price@monash.edu                                      !!
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

PROGRAM ndspmhd
 USE debug
 USE loguns
 USE versn
 IMPLICIT NONE
 INTEGER, PARAMETER :: maxruns = 200
 INTEGER :: i,iprev,irun, nruns
 CHARACTER(LEN=120), DIMENSION(maxruns) :: runname
!
!--version number
!
    version = 'v2.0 [21st Feb 2014]'
!   * major new release containing dust algorithms
!    version = 'v1.0.1 [21st Dec 2010]'
!   * minor bug fix with build
!    version = 'v1.0 [25th Oct 2010]'
!   * first public version
!    version = 'NDSPMHD-basic-Torun2010-v1.0'
!   * basic version released for ASTROSIM summer school in Torun, Poland
 
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
 
END PROGRAM ndspmhd

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
