!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Setup for the blob evaporation problem                                !!
!!                                                                        !!
!!  dense disc of fluid at origin                                         !!
!!                                                                        !!
!!------------------------------------------------------------------------!!

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
 use mem_allocation, only:alloc
 
 use uniform_distributions
 use cons2prim, only:primitive2conservative
!
!--define local variables
!      
 implicit none
 integer :: i,ipart,nphi,nr,ir,iphi
 real :: rmin,rmax,rdisk,phi,deltaphi,ri
 real :: denszero,przero,vphi,pri
 real :: massp,gam1,omega,spsoundi
 real, dimension(ndim) :: xorigin, dx
 real, dimension(ndimv) :: Bzero
 logical, parameter :: equalmass = .true.
!
!--set boundaries
!                        
 ibound = 1       ! fixed
 nbpts = 0        ! no fixed particles initially
 xmin(:) = -5.0        ! unit square
 xmax(:) = 5.0
!
!--setup parameters for the problem
! 
 xorigin(:) = 0.0        ! co-ordinates of the centre of the initial blast
 rdisk = 1.0             ! radius of the initial disk
 rmax  = 5.0             ! maximum radius
 omega = omegafixed
 Bzero(:) = 0.
 massp  = 1.25e-4
 spsoundi = 0.045
 gam1 = gamma - 1.
 !pext = przero

 write(iprint,*) 'Couette-Taylor Instability Problem '
 write(iprint,10) spsoundi,rdisk,omega,gamma
10 format(/,' cs  = ',f10.3,', disk radius = ',f6.3,/, &
            ' rmax   = ',f6.3,', Omega = ',f6.3,', gamma = ',f6.3,/)
!
!--setup uniform density grid of particles (2d) 
!  (determines particle number and allocates memory)
!
 massp = 1.25e-4
 rmin = psep
 nr = int((rmax - rmin)/psep) + 5
 deltaphi = 2.*asin(psep/2.*rdisk)
 nphi = nint((2.*pi)/deltaphi)
 deltaphi = 2.*pi/nphi
 
 call alloc(nr*nphi)
 
 ipart = 0
 do ir = 1,nr
    ri = rmin + (ir-1)*psep
    !deltaphi = 2.*asin(psep/2.*ri)
    !nphi = nint((2.*pi)/deltaphi)
    !deltaphi = 2.*pi/nphi
    !print*,'ri = ',ri,'deltaphi = ',deltaphi,' nphi = ',nphi

    do iphi = 1,nphi
       phi = (iphi-1)*deltaphi
       ipart = ipart + 1
       if (ipart.gt.size(dens)) call alloc(10*size(dens))
       x(1,ipart) = ri*cos(phi) + xorigin(1)
       x(2,ipart) = ri*sin(phi) + xorigin(2)
    enddo
 enddo
 npart = ipart
 ntotal = npart
 
 denszero = (npart*massp)/(pi*rmax**2)
 przero = spsoundi**2*denszero/gamma
 print*,' dens = ',denszero
!
!--now assign particle properties
!
 do ipart=1,ntotal
    dx(:) = x(:,ipart)  !-xorigin(:)
    ri = sqrt(dot_product(dx,dx))
    phi = atan2(x(2,ipart),x(1,ipart))
    vel(:,ipart) = 0.
    if (ri.le.rdisk) then
       itype(ipart) = 2
       nbpts = nbpts + 1
       vphi = omega*ri
       vel(1,ipart) = -vphi*sin(phi)
       vel(2,ipart) = vphi*cos(phi)
    elseif (ri.ge.rmax) then
       itype(ipart) = 1
       nbpts = nbpts + 1
       vel(:,ipart) = 0.
    else
       itype(ipart) = 0
       vel(:,ipart) = 0.
    endif
    hh(ipart) = 1.2*psep
    pmass(ipart) = massp
    dens(ipart) = denszero
    pri = przero
    Bfield(:,ipart) = Bzero(:)
    uu(ipart) = przero/(gam1*dens(ipart))
 enddo
!
!--make sure it is *really* in pressure equilibrium
!
 call set_fixedbound
 call primitive2conservative
 do ipart=1,ntotal
    przero = spsoundi**2*rho(ipart)/gamma
    uu(ipart) = przero/(gam1*rho(ipart))
!!    print*,ipart
    if (sqrt(dot_product(x(:,ipart),x(:,ipart))).lt.1.e-5) then
       print*,'particle very close to zero!!!',i,x(:,i)
    endif
 enddo
 
 Bconst(:) = Bzero(:)
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
            
 return
end
