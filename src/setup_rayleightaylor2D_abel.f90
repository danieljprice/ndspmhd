!----------------------------------------------------------------
!     Set up a Rayleigh-Taylor instability test as
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
 use setup_params, only:psep,pi
 use eos, only:gamma
 use mem_allocation, only:alloc
 
 use uniform_distributions
 use cons2prim, only:primitive2conservative
!
!--define local variables
!            
 implicit none
 integer :: i,ipart
 real :: massp,volume,totmass,massmedium
 real :: denszero,densmedium,przero,psepmedium,pri
 logical, parameter :: equalmass = .false.
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(RT-Abel)'
 if (ndim.ne.2) stop 'error: need ndim=2 for this problem'
!
!--set boundaries
! 	    
 ibound(1) = 2 ! periodic in x
 ibound(2) = 1 ! fixed in y
 nbpts = 0     ! use ghosts not fixed
 xmin(1) = 0. ! set position of boundaries
 xmax(1) = 0.5
 xmin(2) = 0. !- 6.*psep	! set position of boundaries
 xmax(2) = 1. !+ 6.*psep
 if (iexternal_force.ne.9) stop 'need iexternal force = 9 for this problem'
!
!--set up the uniform density grid
!
 denszero = 1.0
 densmedium = 2.0
 przero = denszero/gamma
 psepmedium = psep*(denszero/densmedium)**(1./ndim)
 !massmedium = massp*(densmedium/denszero)
 volume = product(xmax(:)-xmin(:))
 totmass = denszero*volume

!
!--unequal masses setup whole grid
!
 call set_uniform_cartesian(1,psep,xmin,xmax,fill=.true.)
 !call set_uniform_cartesian(2,psep,xmin,xmax,fill=.true.)
 massp = totmass/float(npart)
 massmedium = massp*(densmedium/denszero)
 npart = ntotal
 print*,'npart =',npart
!
!--now assign particle properties
!
 do i=1,ntotal
    vel(:,i) = 0.
    if (x(2,i).gt.0.3 .and. x(2,i).lt.0.7) then
       vel(2,i) = 0. !0.01*(1. + cos(8.*pi*(x(1,i)+0.25)))*&
                     !  (1. + cos(2./0.4*pi*(x(2,i)-0.5)))/4.
    else
       vel(2,i) = 0.
    endif
    dens(i) = densmedium + (denszero-densmedium)/ &
                           (1. + exp(-2.*(x(2,i)-0.5)/0.05))
    pmass(i) = massmedium + (massp-massmedium)/ &
                           (1. + exp(-2.*(x(2,i)-0.5)/0.05))
!    if (x(2,i).gt.0.5) then
!       dens(i) = densmedium
!       pmass(i) = massmedium
!    else
!       dens(i) = denszero
!       pmass(i) = massp
!    endif
    pri = przero - 0.5*dens(i)*(x(2,i) - 0.5)
    uu(i) = pri/((gamma-1.)*dens(i))
    !--set fixed particles in y
    if (x(2,i).gt.0.9 .or. x(2,i).lt.0.1) then
       nbpts = nbpts + 1
       itype(i) = 1
    endif
!    Bfield(:,i) = 0.
!    Bfield(1,i) = 0.5
 enddo
!
!--get rho from a sum and then set u to give a
!  smooth pressure
!
 write(iprint,*) 'calling density to make smooth pressure jump...'
 call primitive2conservative
 do i=1,npart
    pri = przero - 0.5*rho(i)*(x(2,i) - 0.5)
    uu(i) = pri/((gamma-1.)*dens(i))
 enddo
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
  
 return
end

subroutine modify_dump
 use loguns, only:iprint
 use part
 use options, only:imhd
 use timestep, only:time
 use setup_params, only:pi
 implicit none
 integer :: i
!
!--now assign particle properties
!
 write(iprint,*) 'modifying dump with velocities for Rayleigh-Taylor run'
 do i=1,ntotal
    vel(1,i) = 0.
    if (abs(x(2,i)).lt.1./3.) then
       vel(2,i) = 0.5*(1. + cos(4.*pi*x(1,i)))*(1. + cos(3.*pi*x(2,i)))/4.
    else
       vel(2,i) = 0.
    endif
 enddo
 time = 0.
 
end subroutine modify_dump
