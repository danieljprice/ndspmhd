!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Setup for Zisis test problem (SPHERIC 2013)  !!
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
 
 use uniform_distributions
!
!--define local variables
!      
 implicit none
 integer :: ipart
 real :: dens1,dens2
 real :: pri,pr1,pr2,xi,vol2
 real :: totmass,gam1,massp,psep2,totvol
 real, dimension(ndimv) :: Bzero
 logical, parameter :: equalmass=.true.
!
!--set boundaries
!                        
 ibound(:) = 3 ! periodic
 ibound(1) = 1        ! fixed parts
 nbpts = 0        ! no fixed particles
 xmin(:) = -0.3        ! unit square
 xmax(:) = 1.4
!
!--setup parameters for the problem
! 
 dens1 = 1.0           ! density of region 1
 dens2 = 0.5           ! density of region 2
 pr1 = 1.*dens1
 pr2 = 1.*dens2
 Bzero = 0.
 gam1 = gamma - 1.

 write(iprint,*) 'Zisis test problem'
 write(iprint,10) dens2,pr1,pr2
10 format(/,' dens2 = ',f10.3,/, &
            ' pressure(1)  = ',1pe10.3,', pressure(2)=',1pe10.3,/)
!
!--setup uniform density grid of particles (2D)
!  (determines particle number and allocates memory)
!
 psep2 = psep*(dens1/dens2)**(1./ndim)
 write(iprint,*) 'psep in region2 = ',psep2
 if (equalmass) then
    call set_uniform_cartesian(2,psep,xmin,xmax,mask=-6)
    call set_uniform_cartesian(2,psep2,xmin,xmax,mask=6)
 else
    call set_uniform_cartesian(1,psep,xmin,xmax)  ! 2 = close packed arrangement
 endif
 ntotal = npart
!
!--determine particle mass in ambient medium
!
 if (equalmass) then
    vol2 = (1.1 - 0.9)
    if (ndim.ge.2) vol2 = vol2*product(xmax(2:)-xmin(2:))
    totvol = product(xmax(:)-xmin(:)) - vol2
    totmass = dens1*totvol + dens2*vol2
    massp = totmass/float(ntotal)
 else
    totmass = dens1*product(xmax(:)-xmin(:))
    massp = totmass/float(ntotal) ! particle mass
 endif
!
!--now assign particle properties
!
 nbpts = 0
 do ipart=1,ntotal
    xi = x(1,ipart)
    if (xi > 0.9) then
       dens(ipart) = dens2
       if (equalmass) then
          pmass(ipart) = massp
       else
          pmass(ipart) = massp*dens2/dens1
       endif
       pri = pr2
       itype(ipart) = itypegas2
    else
       pmass(ipart) = massp
       dens(ipart) = dens1
       pri = pr1
       itype(ipart) = itypegas1
    endif  
    vel(:,ipart) = 0.
    if (xi.lt.0.5) vel(1,ipart) = 1.
    Bfield(:,ipart) = Bzero(:)
    uu(ipart) = 1. !pri/(gam1*dens(ipart))
    if (xi > 1.2 .or. xi < -0.1) then
       itype(ipart) = itypebnd
       nbpts = nbpts + 1
    endif
 enddo
 Bconst(:) = Bzero(:)
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
            
 return
end

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
