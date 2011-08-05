!-------------------------------------------------------------------------
! Set up a uniform density grid of particles in cylindrical co-ordinates
!-------------------------------------------------------------------------

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
 if (trace) write(iprint,*) ' Entering subroutine setup(unifdis)'
!
!--set co-ordinate system
!
 igeom = 2
!
!--set boundaries
! 	    
 nbpts = 0	! use ghosts not fixed
 ibound(1) = 2	! reflective in r
 xmin(1) = 0.0	! rmin
 xmax(1) = 0.5*pi  ! rmax

 if (ndim.ge.2)then
    ibound(2) = 3  ! periodic in phi
    xmin(2) = 0.0  ! phi_min
    xmax(2) = 2.*pi  ! phi_max
 endif
 if (ndim.ge.3) then
    ibound(3) = 0  ! nothing in z
 endif
!
!--set up the uniform density grid
!
 call set_uniform_cartesian(1,psep,xmin,xmax,.false.)
 ntotal = npart
!
!--determine particle mass (WORK THIS OUT)
!
 denszero = 1.0
 volume = product(xmax(:)-xmin(:))
 totmass = denszero*volume
 massp = totmass/FLOAT(ntotal) ! average particle mass
!
!--now assign particle properties
!  (note *do not* setup smoothing length as it depends on the conservative
!   density rho which has not yet been calculated) 
!
 do i=1,ntotal
    vel(:,i) = 0.
    dens(i) = denszero
    pmass(i) = massp
    uu(i) = 1.0	! isothermal
    Bfield(:,i) = 0.
 enddo 
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  Exiting subroutine setup'
  
 return
end
