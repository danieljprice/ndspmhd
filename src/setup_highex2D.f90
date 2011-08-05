!!-------------------------------------------------------------------------!!
!!                                                                         !!
!!  setup for 2D high explosives test                                      !!
!!                                                                         !!
!!  density is set to unity all over, whilst the pressure (or equivalently !!
!!  the thermal energy) is set to some large quantity in a small circle    !!
!!  around the origin                                                      !!
!!                                                                         !!
!!  magnetic field of strength 10g in the x-direction                      !!
!!                                                                         !!                                                                        !!
!!-------------------------------------------------------------------------!!

subroutine setup
!
!--include relevant global variables
!
 use dimen_mhd, only:ndim,ndimV
 use debug, only:trace
 use loguns, only:iprint
 use bound
 use eos
 use kernels, only:interpolate_kernel
 use options, only:damp,ibound
 use part
 use setup_params, only:psep,pi
 
 use uniform_distributions
!
!--define local variables
!      
 implicit none
 integer :: ipart
 real :: denszero,przero,densc,uuc,uuzero,rc,prc
 real :: totmass,gam1,massp,psepc,radius
 real, dimension(ndim) :: xmaxtemp

!
!--set boundaries
!            	    
 ibound = 2     ! fixed ghosts
 nbpts = 0
 xmin(:) = 0.     ! same xmin in all dimensions
 xmax(:) = 2.
 denszero = 1.293
 densc = 1630.
 rc = 0.0527
 przero = 10.1325   ! 1 atm in pascal
 prc = 27603.94
 uuc = prc/densc  !4.29e6 J ! should 
 

 write(iprint,10) ndim
 write(iprint,20) uuc,rc,denszero,przero
10 format(/,1x,i1,'-dimensional high explosives test')
20 format(/,' central energy  = ',f10.3,', blast radius = ',f6.3,/, &
            ' ext. density = ',f6.3,', ext. pressure = ',f6.3,/)
 
 gam1 = gamma - 1.
 if (abs(gam1).lt.1.e-3) stop 'eos cannot be isothermal for this setup'
 uuzero = przero/(gam1*denszero)
!
!--setup uniform density grid of particles
!  (determines particle number and allocates memory)
!
 psepc = psep*(denszero/densc)**(1./ndim)
 xmaxtemp = rc
 call set_uniform_cartesian(1,psepc,xmin,xmaxtemp,rmax=rc) ! 2 = close packed arrangement
 call set_uniform_cartesian(1,psep,xmin,xmax,rmin=rc)
 ntotal = npart
!
!--determine particle mass
!
 totmass = denszero*(product(xmax(:)-xmin(:)) - pi*rc**2)   ! assumes cartesian boundaries
 massp = totmass/float(ntotal)                  ! average particle mass
!
!--now assign particle properties
! 
 do ipart=1,ntotal
    vel(:,ipart) = 0.
!--uniform density and smoothing length
    dens(ipart) = denszero

    radius = sqrt(dot_product(x(:,ipart),x(:,ipart)))
    if (radius.le.rc) then
       dens(ipart) = densc
       uu(ipart) = uuc
    else
       dens(ipart) = denszero
       uu(ipart) = uuzero
    endif
    pmass(ipart) = massp
 enddo

!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
            
 return
end
