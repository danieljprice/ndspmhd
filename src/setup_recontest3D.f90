!----------------------------------------------------------------
!     Set up a uniform density cartesian grid of particles in ND
!----------------------------------------------------------------

subroutine setup
!
!--include relevant global variables
!
 use dimen_mhd
 use debug
 use loguns
 use bound
 use options
 use part
 use setup_params
 
 use uniform_distributions
!
!--define local variables
!            
 implicit none
 integer :: i
 real :: massp,volume,totmass
 real :: denszero
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(unifdis)'
!
!--set boundaries
! 	    
 ibound = 3	! boundaries
 nbpts = 0	! use ghosts not fixed
 xmin(:) = -0.5	! set position of boundaries
 xmax(:) = 0.5
!
!--set up the uniform density grid
!
 call set_uniform_cartesian(1,psep,xmin,xmax,.false.)
 npart = ntotal
 print*,'npart =',npart
!
!--determine particle mass
!
 denszero = 1.0
 volume = product(xmax(:)-xmin(:))
 totmass = denszero*volume
 massp = totmass/float(ntotal) ! average particle mass
!
!--now assign particle properties
! 
 do i=1,ntotal
    vel(:,i) = 0.
    if (x(1,i).gt.0.45 .and. x(1,i).lt.0.55) then
       vel(1,i) = 150.*SIN(pi*(x(3,i)-0.5))**36
    endif
    !!!vel(1,i) = x(1,i)
    dens(i) = denszero
    pmass(i) = massp
    uu(i) = 1.5	! isothermal
    Bfield(:,i) = 0.
    Bfield(3,i) = 5.
 enddo 
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
  
 return
end
