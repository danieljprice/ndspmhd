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
 use eos
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
 real :: denszero,uuzero,cs0,polyk0,Rcentre
 real :: przero,Bzeroz,asize,betamhd
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(unifdis)'
 write(iprint,*) '2D magneto-rotational instability'
 if (ndim.ne.2) stop 'error: this is a 2D problem and ndim.ne.2'
 if (ndimV.ne.3) stop 'error: we need ndimV=3 for this problem'
!
!--geometry is cylindrical r-z
!
 geom = 'cylrzp'
 geomsetup = 'cylrzp'
!
!--set position of disc patch for coriolis & centrifugal forces
!
 Rcentre = 100.
 Omega = (1./Rcentre**1.5)
 Omega2 = (1./Rcentre**3)
 iexternal_force = 2 ! 1/r^2 force
!
!--set boundaries
!
 ibound(1) = 2  ! reflecting in x (==R)  
 ibound(2) = 3	! periodic in y   (==z)
 nbpts = 0	! use ghosts not fixed
 asize = 1.0    ! box size (corresponds to a in Hawley/Balbus 1992)
 xmin(1) = Rcentre-asize  ! set position of boundaries
 xmax(1) = Rcentre+asize
 xmin(2) = -0.5*asize
 xmax(2) = 0.5*asize
!
!--set up the uniform density grid
!
 call set_uniform_cartesian(11,psep,xmin,xmax,.false.,perturb=0.1)
 ntotal = npart
!
!--determine particle mass
!
 denszero = 1.0
 volume = product(xmax(:)-xmin(:))
 totmass = denszero*volume
 massp = totmass/float(ntotal) ! average particle mass
!
!--sound speed
!
 przero = 1.e-5
 cs0 = gamma*przero/denszero
 uuzero = przero/(denszero*(gamma-1.))
 polyk0 = (gamma-1.)*denszero**(1.-gamma)*uuzero
 write(iprint,*) ' cs0 = ',cs0, ' u = ',uuzero
 write(iprint,*) ' polyk = ',polyk,' should be = ',polyk0
 polyk = polyk0

 write(iprint,*) ' Omega_c = ',Omega, ' R_c = ',Rcentre

!
!--field strength
!
 betamhd = 1000.
 Bzeroz = sqrt(2.*przero/betamhd)
 write(iprint,*) ' beta(MHD) = ',betamhd,' Bz = ',Bzeroz,' P = ',przero
!
!--now assign particle properties
! 
 do i=1,ntotal
    vel(:,i) = 0.
    dens(i) = denszero
    pmass(i) = massp
    uu(i) = uuzero ! isothermal
    Bfield(:,i) = 0.
    !
    !--set constant z field in region between -a/5 -> a/5
    !
    if (x(1,i).gt.(Rcentre-0.2*asize) .and. x(1,i).lt.(Rcentre+0.2*asize)) then
       Bfield(2,i) = Bzeroz
    endif
    vel(3,i) = -1.5*Omega*x(1,i)
 enddo 
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
  
 return
end
