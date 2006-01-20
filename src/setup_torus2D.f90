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
 use eos
 use rates, only:force
 
 use uniform_distributions
!
!--define local variables
!            
 implicit none
 integer :: i,nrings,iring,ipart,npartphi
 real :: Rtorus,dfac,Mstar,Mtorus
 real :: massp,r_in,r_out,deltar,polyn,sumA
 real :: ri,zi,rhofac,deltaphi,densi,phii,pri
 real :: deltartemp,denstemp,rtemp,deltar0,dens0
 real :: rhat(3), omegai,rpart,frad,v2onr,deltarprev
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(torus)'
!
!--set boundaries
! 	    
 ibound = 0	! boundaries
 nbpts = 0	! use ghosts not fixed
 xmin(:) = 0.	! set position of boundaries
 xmax(:) = 0.
 iexternal_force = 6
!
!--set up the uniform density grid
!
 Rtorus = 1.0
 dfac = 1.1
 Mstar = 1.0
 Mtorus = 4.e-4
  
 write(iprint,*) 'Setup for beer-related torus thing (Owen, Laure)'
 write(iprint,*) 'enter particle number '
 read*,npart
 ntotal = npart

 call alloc(2*ntotal)

 massp = Mtorus/real(ntotal)
!
!--integrate to get total mass (in order to set A in P = A*rho^gamma)
!  (Owen was particularly fussy about using A not K)
!
 r_in = 0.1
 r_out = 1.5
 nrings = 50
 deltar = (r_out - r_in)/real(nrings-1)
 polyn = 1./(gamma-1.)
 sumA = 0.
 do iring=1,nrings
    ri = r_in + (iring-1)*deltar
    zi = 0.
    sumA = sumA + 2.*pi*ri*deltar*rhofunc(ri,zi,polyn,dfac,Mstar,Rtorus)
 enddo
!
!--rearrange to get A (which I call polyk)
!
 polyk = ((sumA/Mtorus)**(1./polyn))/(polyn + 1.)
 print*,' polyk = ',polyk
!--rhofac is the factor which multiplies rhofunc to give rho
!  this is 1/B in Owen lingo
 rhofac = 1./(polyk*(polyn+1.))**polyn

!
!--setup first ring at r=Rtorus
!
 ipart = 0
 ri = Rtorus
 zi = 0.
 densi = rhofac*rhofunc(ri,zi,polyn,dfac,Mstar,Rtorus)
!--calculate delta r and delta phi from rho
 deltar = sqrt(massp/densi)
 deltaphi = deltar/ri
 npartphi = int(2.*pi/deltaphi)
!--correct deltaphi to integer number of particles
 deltaphi = 2.*pi/npartphi
!--correct deltar to match new deltaphi
 deltar = ri*deltaphi
!--now setup ring of particles
 do i = 1,npartphi
    ipart = ipart + 1
    phii = (i-1)*deltaphi
    x(1,ipart) = ri*COS(phii)
    x(2,ipart) = ri*SIN(phii)
    dens(ipart) = densi
    pmass(ipart) = massp
    pri = polyk*densi**gamma
    uu(ipart) = pri/(densi*(gamma-1.))
 enddo
 deltar0 = deltar
 dens0 = densi

!
!--setup rings from Rtorus outwards until dens < 0
!
 do while (densi > tiny(dens))
    !--take half step in r using current deltar
    rtemp = ri + 0.5*deltar
    !--get density
    denstemp = rhofac*rhofunc(rtemp,zi,polyn,dfac,Mstar,Rtorus)
    !--calculate delta r and delta phi from rho
    deltartemp = sqrt(massp/denstemp)
    !--corrector step on r using midpoint density
    ri = ri + deltartemp
    !--get density
    densi = rhofac*rhofunc(ri,zi,polyn,dfac,Mstar,Rtorus)
    if (densi.gt.tiny(densi)) then
       !--get new deltar at new position
       deltar = sqrt(massp/densi)    
       deltaphi = deltar/ri
       npartphi = int(2.*pi/deltaphi)
       !--correct deltaphi to integer number of particles
       deltaphi = 2.*pi/npartphi

       if (ri*deltaphi.lt.2.0*deltartemp) then
          !--correct deltar to match new deltaphi
          deltar = ri*deltaphi
          !--now setup ring of particles
          do i = 1,npartphi
             ipart = ipart + 1
             phii = (i-1)*deltaphi
             x(1,ipart) = ri*COS(phii)
             x(2,ipart) = ri*SIN(phii)
             dens(ipart) = densi
             pmass(ipart) = massp
             pri = polyk*densi**gamma
             uu(ipart) = pri/(densi*(gamma-1.))
          enddo
       else
          deltar = ri*deltaphi
          print*,'skipping ring, too few particles : ',npartphi
       endif
    endif
    
 enddo
! 
!--setup rings from Rtorus inwards until dens < 0
!
 ri = Rtorus
 deltar = deltar0
 densi = dens0
 do while (densi > tiny(dens))
    !--take half step in r using current deltar
    rtemp = ri - 0.5*deltar
    !--get density
    denstemp = rhofac*rhofunc(rtemp,zi,polyn,dfac,Mstar,Rtorus)
    !--calculate delta r and delta phi from rho
    deltartemp = sqrt(massp/denstemp)
    !--corrector step on r using midpoint density
    ri = ri - deltartemp
    !--get density
    densi = rhofac*rhofunc(ri,zi,polyn,dfac,Mstar,Rtorus)
    if (densi.gt.tiny(densi)) then
       !--get new deltar at new position
       deltar = sqrt(massp/densi)    
       deltaphi = deltar/ri
       npartphi = int(2.*pi/deltaphi)
       !--correct deltaphi to integer number of particles
       deltaphi = 2.*pi/npartphi
       if (ri*deltaphi.lt.2.0*deltartemp) then
          !--correct deltar to match new deltaphi
          deltar = ri*deltaphi
          !--now setup ring of particles
          do i = 1,npartphi
             ipart = ipart + 1
             phii = (i-1)*deltaphi
             x(1,ipart) = ri*COS(phii)
             x(2,ipart) = ri*SIN(phii)
             dens(ipart) = densi
             pmass(ipart) = massp
             pri = polyk*densi**gamma
             uu(ipart) = pri/(densi*(gamma-1.))
          enddo
       else
          deltar = ri*deltaphi
          print*,'skipping ring, too few particles : ',npartphi
       endif
    endif
 enddo
 npart = ipart
 ntotal = ipart

 call primitive2conservative
 !
 !--balance pressure forces with centripedal acceleration
 !
 vel = 0.
! do i=1,npart
!    rpart = sqrt(dot_product(x(1:2,i),x(1:2,i)))
!    rhat(1:2) = x(1:2,i)/rpart
!    rhat(3) = 0.
!    frad = abs(dot_product(force(:,i),rhat(:)))
!    omegai = sqrt(frad)/rpart
!    vel(1,i) = -omegai*x(2,i)
!    vel(2,i) = omegai*x(1,i)
!    if (ndimV.ge.3) vel(3,i) = 0.
    !!print*,'forcez = ',force(3,i),force(2,i)
    !v2onr = dot_product(vel(1:2,i),vel(1:2,i))/rpart
    !print*,'v2onr = ',v2onr*rhat(2),force(2,i),v2onr*rhat(1),force(1,i)
    
! enddo

!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
  
 return

contains

real function rhofunc(rr,zz,polyn,dd,Mstar,R0)
 implicit none
 real, intent(in) :: rr,zz,polyn,dd,Mstar,R0
 real :: term
!
!--functional form of rho/(A(n+1))^n
!
 term = Mstar/R0*(R0/sqrt(rr**2 + zz**2) - 0.5*(R0/rr)**2 - 1./(2.*dd))
 if (term.gt.tiny(term)) then
    rhofunc = term**polyn
 else
    rhofunc = 0.
 endif

end function rhofunc

end subroutine setup

