!-------------------------------------------------------------------------
! Set up a Keplerian disc with a radial density profile in 2D
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
 real :: denszero,rr
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' Entering subroutine setup(unifdis)'
!
!--set co-ordinate system
!
 igeom = 2
 iexternal_force = 2
!
!--set boundaries
! 	    
 nbpts = 0	! use ghosts not fixed
 ibound(1) = 2	! reflective in r
 xmin(1) = 0.5	! rmin
 xmax(1) = 1.5  ! rmax

 if (ndim.ge.2)then
    ibound(2) = 3  ! periodic in phi
    xmin(2) = 0.0  ! phi_min
    xmax(2) = 2.*pi  ! phi_max
 endif
 if (ndim.ge.3) then
    ibound(3) = 0  ! nothing in z
 endif
!
!--set up the uniform density grid (uniform in r,phi)
!
 call set_uniform_cartesian(1,psep,xmin,xmax,.false.)
 ntotal = npart
!
!--if not cylindrical coordinates, translate to cartesians
!
 do i=1,npart
    CALL coord_transform(x(:,i),ndim,igeom,x(:,i),ndim,1)
 enddo
 igeom = 1
 ibound = 0    ! no boundaries

!
!--determine particle mass (WORK THIS OUT)
!
 denszero = 1.0
 volume = pi*(xmax(1)**2 - xmin(1)**2)
 totmass = denszero*volume
 massp = totmass/FLOAT(ntotal) ! average particle mass
!
!--now assign particle properties
!  (note *do not* setup smoothing length as it depends on the conservative
!   density dens which has not yet been calculated) 
!
 do i=1,ntotal
!
!--vel is keplerian
!
    rr = SQRT(DOT_PRODUCT(x(:,i),x(:,i)))
    vel(:,i) = 0.
    dens(i) = denszero*(1./rr)
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
