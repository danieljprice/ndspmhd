!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Setup for cylindrical interacting blast wave as in Borve et al. 2005  !!
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
 real :: denszero,densdisc,przero
 real :: pri,rdisc,radius,radius2,prdisc1,prdisc2
 real :: totmass,gam1,massp,psepdisc,totvol
 real, dimension(ndim) :: xorigin1, xorigin2, dx
 real, dimension(ndimv) :: Bzero
 logical, parameter :: equalmass=.true.
!
!--check number of dimensions is right
!
 if (ndim.ne.2) stop ' ndim must be = 2 for this problem'
 if (ndimv.lt.2) stop ' ndimv>=2 for this problem'
!
!--set boundaries
!                        
 ibound = 2        ! reflective ghosts (boundaries not important in this problem)
 nbpts = 0        ! no fixed particles
 xmin(:) = 0.0        ! unit square
 xmax(:) = 1.0
!
!--setup parameters for the problem
! 
 xorigin1(:) = 0.0        ! co-ordinates of the centre of the initial blast
 xorigin2(:) = 1.0        ! co-ordinates of the centre of the initial blast
 rdisc = 0.4             ! radius of the initial disc
 przero = 0.01             ! initial pressure
 denszero = 1.0            ! ambient density
 densdisc = 10.0           ! density of rotating disc
 prdisc1 = 1.e4
 prdisc2 = 1.e3
 Bzero = 0.
 gam1 = gamma - 1.

 write(iprint,*) 'two dimensional cylindrical blast wave problem '
 write(iprint,10) densdisc,rdisc,prdisc1,prdisc2,przero
10 format(/,' central density  = ',f10.3,', disc radius = ',f6.3,/, &
            ' disc pressure(1)  = ',1pe10.3,', disc pressure(2)=',1pe10.3,', p_0 = ',0pf6.3,/)
!
!--setup uniform density grid of particles (2D)
!  (determines particle number and allocates memory)
!
 psepdisc = psep*(denszero/densdisc)**(1./ndim)
 write(iprint,*) 'psep in disc = ',psepdisc
 if (equalmass) then
    call set_uniform_cartesian(2,psep,xmin,xmax,mask=-1)
    call set_uniform_cartesian(2,psepdisc,xmin,xmax,mask=1)
 else
    call set_uniform_cartesian(1,psep,xmin,xmax)  ! 2 = close packed arrangement
 endif
 ntotal = npart
!
!--determine particle mass in ambient medium
!
 if (equalmass) then
    totvol = product(xmax(:)-xmin(:)) - 0.5*pi*rdisc**2
    totmass = denszero*totvol + densdisc*0.5*pi*rdisc**2
    massp = totmass/float(ntotal)
 else
    totmass = denszero*product(xmax(:)-xmin(:))
    massp = totmass/float(ntotal) ! particle mass
 endif
!
!--now assign particle properties
! 
 do ipart=1,ntotal
    dx(:) = x(:,ipart)-xorigin1(:) 
    radius = sqrt(dot_product(dx,dx))
    dx(:) = x(:,ipart)-xorigin2(:) 
    radius2 = sqrt(dot_product(dx,dx))
    if (radius.le.rdisc) then
       dens(ipart) = densdisc
       if (equalmass) then
          pmass(ipart) = massp
       else
          pmass(ipart) = massp*densdisc/denszero
       endif
       pri = prdisc1 
    elseif (radius2.le.rdisc) then 
       dens(ipart) = densdisc
       if (equalmass) then
          pmass(ipart) = massp
       else
          pmass(ipart) = massp*densdisc/denszero
       endif
       pri = prdisc2
    else
       pmass(ipart) = massp
       dens(ipart) = denszero
       pri = przero 
    endif  
    vel(:,ipart) = 0.
    Bfield(:,ipart) = Bzero(:)
    uu(ipart) = pri/(gam1*dens(ipart))
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
