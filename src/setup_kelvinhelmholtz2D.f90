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
!  ChangeLog:                                                                  !
!------------------------------------------------------------------------------!

!----------------------------------------------------------------
!     Set up a Kelvin-Helmholtz instability test as
!     in Jim Stone et al. (Athena test suite)
!     Also does Wang et al. (2010)
!----------------------------------------------------------------

!--stuff we need in other routines
module khsetup
 implicit none
 real, parameter :: denszero = 1.
 real, parameter :: densmedium = 2.
 real, parameter :: vzero = 0.5
 real, parameter :: vmedium = -0.5
 real, parameter :: smoothl = 0.025
 real, parameter :: przero = 2.5
 logical, parameter :: wang = .false.     ! use EITHER wang OR robertson
 logical, parameter :: robertson = .false.

contains
!
!--y profile of density and velocity for Wang setup
!
real function yprofile(yi,azero,amedium,smoothl)
 implicit none
 real, intent(in) :: yi,azero,amedium,smoothl
 real :: amid,expterm

 amid = 0.5*(azero - amedium)
 if (yi.gt.0.75) then
    expterm = exp(-(yi-0.75)/smoothl)
    yprofile = azero - amid*expterm
 elseif (yi.gt.0.5) then
    expterm = exp(-(0.75-yi)/smoothl)
    yprofile = amedium + amid*expterm
 elseif (yi.gt.0.25) then
    expterm = exp((-yi + 0.25)/smoothl)
    yprofile = amedium + amid*expterm
 else
    expterm = exp((yi - 0.25)/smoothl)
    yprofile = azero - amid*expterm
 endif

end function yprofile

end module khsetup

!
!--main setup routine
!
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
 use khsetup, only:denszero,densmedium,vzero,vmedium,smoothl,yprofile,&
                   przero,wang,robertson
!
!--define local variables
!            
 implicit none
 integer :: i,iseed,ipart
 real :: massp,volume,totmass !,ran1
 real :: psepmedium
 real :: densmid,vmid,expterm,yi
 real, dimension(ndim) :: xminregion,xmaxregion
 logical, parameter :: equalmass = .true.
 logical, parameter :: symmetric = .false.
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(kh)'
 if (ndim.ne.2) stop 'error: need ndim=2 for this problem'
 if (robertson .and. wang) stop 'error, need robertson OR wang, not both'
!
!--set boundaries
! 	    
 nbpts = 0      ! use ghosts not fixed
 if (robertson) then
    xmin(:) = 0.
    xmax(:) = 1.
 else
    xmin(:) = -0.5 ! set position of boundaries
    xmax(:) = 0.5
 endif
!
!--set up the uniform density grid
!
 if (.not.equalmass) then
    if (wang) stop 'wang needs equal mass'
!
!--unequal masses setup whole grid
! 
    call set_uniform_cartesian(2,psep,xmin,xmax,fill=.true.)
    volume = product(xmax(:)-xmin(:))
    totmass = denszero*volume
    massp = totmass/float(npart)
    if (robertson) write(iprint,*) 'using robertson setup'
 else
    if (robertson) stop 'robertson setup with equal mass not implemented'
!
!--for equal mass particles setup each region separately
!
    xminregion(1) = xmin(1)
    xmaxregion(1) = xmax(1)
    xminregion(2) = -0.5
    xmaxregion(2) = -0.25
    call set_uniform_cartesian(2,psep,xminregion,xmaxregion,fill=.true.)
!
!--determine particle mass from npart in first region
!
    volume = product(xmaxregion(:)-xminregion(:))
    totmass = denszero*volume
    massp = totmass/float(npart) ! average particle mass

    if (.not.symmetric) then
       xminregion(2) = 0.25
       xmaxregion(2) = 0.5
       call set_uniform_cartesian(2,psep,xminregion,xmaxregion,fill=.true.)
    endif
!
!--setup -0.25 < y < 0.25
! 
    xminregion(2) = -0.25
    if (symmetric) then
       xmaxregion(2) = 0.
    else
       xmaxregion(2) = 0.25
    endif
    psepmedium = psep*(denszero/densmedium)**(1./ndim)
    call set_uniform_cartesian(2,psepmedium,xminregion,xmaxregion,fill=.true.)

    if (symmetric) then
    !
    !--reallocate memory to new size of list
    !
       call alloc(2*npart)
    !
    !--reflect particles above and below the y axis
    !
       ipart = npart
       do i=1,npart
          ipart = ipart + 1
          x(1,ipart) = x(1,i)
          x(2,ipart) = -x(2,i)
       enddo
       npart = ipart
       ntotal = npart
    endif
 endif
 npart = ntotal
 print*,'npart =',npart
!
!--now declare periodic boundaries
!
 ibound = 3     ! boundaries
 iseed = -23864

 vmid = 0.5*(vzero - vmedium)
 densmid = 0.5*(denszero - densmedium)
!
!--now assign particle properties
! 
 do i=1,ntotal
    vel(:,i) = 0.
    if (robertson) then
       vel(1,i) = vzero + Rfunc(x(2,i))*(vmedium-vzero)
       vel(2,i) = 0.01*sin(2.*pi*x(1,i))
       dens(i) = denszero + Rfunc(x(2,i))*(densmedium-denszero)
       pmass(i) = massp*(dens(i)/denszero)
    elseif (wang) then
       yi = x(2,i) - xmin(2)
       vel(1,i) = yprofile(yi,vzero,vmedium,smoothl)
       dens(i)  = yprofile(yi,denszero,densmedium,smoothl)       
       vel(2,i) = 0.01*sin(2.*2.*pi*(x(1,i)-xmin(1)))
       pmass(i) = 1.5/real(npart)
    else
       if (abs(x(2,i)).lt.0.25) then
          vel(1,i) = vmedium
          dens(i) = densmedium
          if (equalmass) then
             pmass(i) = massp
          else
             pmass(i) = massp*(densmedium/denszero)
          endif
       else
          vel(1,i) = vzero
          dens(i) = denszero
          pmass(i) = massp
       endif
       !
       !--add velocity perturbation
       !
       !!vel(1,i) = vel(1,i)*(1.0 + 0.01*(ran1(iseed)-0.5))
       !!vel(2,i) = 0.01*(ran1(iseed)-0.5)
       if (abs(x(2,i)-0.25).lt.0.025) then
          vel(2,i) = 0.025*sin(-2.*pi*(x(1,i)+0.5)*6.)
       elseif (abs(x(2,i)+0.25).lt.0.025) then
          vel(2,i) = 0.025*sin(2.*pi*(x(1,i)+0.5)*6.)
       endif
    endif
    if (iener.gt.0) then
       uu(i) = przero/((gamma-1.)*dens(i))
    else
       uu(i) = 1.5
    endif
    Bfield(:,i) = 0.
    Bfield(1,i) = 0.5
 enddo

 if (wang) call stretch_to_make_density_profile(x,npart,xmin(2))
 
 !!omegakh = sqrt(densmedium/denszero)*1.0/(denszero + densmedium)
 print*,' tau_kh = ',sqrt(densmedium/denszero)*(1./6.)
!
!--get rho from a sum and then set u to give a
!  smooth pressure
!
 if (iener.gt.0) then
    write(iprint,*) 'calling density to make smooth pressure jump...'
    call primitive2conservative
    do i=1,npart
       uu(i) = przero/((gamma - 1.)*rho(i))
    enddo
 endif
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
  
 return
 
contains

real function Rfunc(y)
 implicit none
 real, parameter  :: delta = 0.05
 real, intent(in) :: y
 
 Rfunc = 1./(1. + exp(2.*(y-0.25)/delta)) !* 1./(1. + exp(2.*(0.75-y)/delta))
 print*,' Rfunc = ',Rfunc
 
end function Rfunc

subroutine stretch_to_make_density_profile(x,np,ymin)
 implicit none
 real, dimension(:,:), intent(inout) :: x
 integer, intent(in) :: np
 real, intent(in) :: ymin
 integer :: i,its
 real :: yi,ynew
 
 do i=1,np
    yi = x(2,i) - ymin
    !--solve M(y) = 0.
    ynew = yi
    do its=1,12 ! overkill, but not costly
       ynew = ynew - (Massprofile(ynew) - Massprofile_orig(yi))/yprofile(ynew,denszero,densmedium,smoothl)
    enddo
    x(2,i) = ymin + ynew
 enddo

end subroutine stretch_to_make_density_profile
!
!--desired mass profile
!
real function Massprofile(yi)
 implicit none
 real, intent(in) :: yi
 real :: sum,sum1,sum2,sum3,sum4,densmid

 densmid = 0.5*(denszero - densmedium)
 sum1 = 0.25*denszero + (exp(-0.25/smoothl) - 1.)*smoothl*densmid
 sum2 = 0.25*densmedium + (1. - exp(-0.25/smoothl))*smoothl*densmid
 sum3 = 0.25*densmedium + (1. - exp(-0.25/smoothl))*smoothl*densmid

 if (yi.gt.0.75) then
    expterm = exp(-(yi-0.75)/smoothl)
    sum4 = denszero*(yi-0.75) + (expterm - 1.)*smoothl*densmid
    sum  = sum1 + sum2 + sum3 + sum4
    !print*,'sum = ',sum1,sum2,sum3,sum4
 elseif (yi.gt.0.5) then
    expterm = exp(-(0.75-yi)/smoothl)
    sum3 = densmedium*(yi-0.5) + (expterm - exp(-0.25/smoothl))*smoothl*densmid
    sum  = sum1 + sum2 + sum3
 elseif (yi.gt.0.25) then
    expterm = exp((-yi + 0.25)/smoothl)
    sum2 = densmedium*(yi-0.25) + (1. - expterm)*smoothl*densmid
    sum  = sum1 + sum2
 else
    expterm = exp((yi - 0.25)/smoothl)
    sum = denszero*yi - exp(-0.25/smoothl)*(exp(yi/smoothl) - 1.)*smoothl*densmid
 endif
 Massprofile = sum

end function Massprofile

!
!--mass profile that we started from
!
real function Massprofile_orig(yi)
 implicit none
 real, intent(in) :: yi
 real :: sum

 if (yi.gt.0.75) then
    sum = 0.25*denszero + 0.5*densmedium + (yi-0.75)*denszero
 elseif (yi.gt.0.25) then
    sum = 0.25*denszero + (yi-0.25)*densmedium
 else
    sum = yi*denszero
 endif
 Massprofile_orig = sum

end function Massprofile_orig
 
end subroutine setup

subroutine modify_dump
 use loguns,       only:iprint
 use part
 use timestep,     only:time
 use setup_params, only:pi
 use bound,        only:xmin
 use eos,          only:gamma
 use khsetup,      only:smoothl,vzero,vmedium,yprofile,wang,przero
 implicit none
 integer :: i
 real :: yi
!
!--now assign particle properties
!
 write(iprint,*) 'modifying dump with velocities for Kelvin-Helmholtz run'
 !print*,vzero,vmedium,wang,przero
 !read*
 do i=1,ntotal
    vel(:,i) = 0.
    if (wang) then
       yi = x(2,i) - xmin(2)
       vel(1,i) = yprofile(yi,vzero,vmedium,smoothl)
       vel(2,i) = 0.01*sin(2.*2.*pi*(x(1,i)-xmin(1)))
       uu(i) = przero/((gamma - 1.)*dens(i))
       !print*,uu(i),dens(i),przero
    else
       if (abs(x(2,i)).lt.0.25) then
          vel(1,i) = 0.5 
       else
          vel(1,i) = -0.5
       endif
       !
       !--add random velocity perturbation
       !
       !!vel(1,i) = vel(1,i)*(1.0 + 0.01*(ran1(iseed)-0.5))
       !!vel(2,i) = 0.01*(ran1(iseed)-0.5)
       if (abs(x(2,i)-0.25).lt.0.025 .or. abs(x(2,i)+0.25).lt.0.025) then
          vel(2,i) = 0.025*sin(2.*pi*(x(1,i)+0.5)*6.)
       endif
    endif
 enddo

 time = 0.
 print*,'time = ',time
 
end subroutine modify_dump
