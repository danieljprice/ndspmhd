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
 use geometry
 
 use uniform_distributions
 use mem_allocation, only:alloc
!
!--define local variables
!            
 implicit none
 integer :: i,nrings,nlayers,iring,iz,ipart,npartphi
 real :: Rtorus,dfac,Mstar,Mtorus,zmax,deltaz
 real :: massp,r_in,r_out,deltar,polyn,sumA
 real :: ri,zi,rhofac,deltaphi,densi,pri
 real :: deltartemp,denstemp,rtemp,deltar0,dens0
 real :: rhat(3), omegai,rpart,frad,v2onr,rcyl2,rcyl,rsph
 real :: beta,Bzi,dbeta,densmax
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
 iexternal_force = 2
!
!--set up the uniform density grid
!
 Rtorus = 1.0
 dfac = 1.1
 Mstar = 1.0
 Mtorus = 4.e-4
 densmax = 1.0 ! set maximum density in torus
  
 write(iprint,*) 'Setup for beer-related torus thing (Owen, Laure)'
 write(iprint,*) 'enter approximate particle number '
 read*,npart
 ntotal = npart

 call alloc(2*ntotal)

!
!--integrate to get total mass (in order to set A in P = A*rho^gamma)
!  (Owen was particularly fussy about using A not K)
!
 r_in = 0.1
 r_out = 1.5
 nrings = 50
 nlayers = 20
 deltar = (r_out - r_in)/real(nrings-1)
 zmax = 0.5
 deltaz = 2.*zmax/real(nlayers-1)
 polyn = 1./(gamma-1.)
! sumA = 0.
! do iring=1,nrings
!    ri = r_in + (iring-1)*deltar
!    zi = 0.
!    sumA = sumA + 2.*pi*ri*deltar*rhofunc(ri,zi,polyn,dfac,Mstar,Rtorus)
! enddo
!--rearrange to get A (which I call polyk)
! polyk = ((sumA/Mtorus)**(1./polyn))/(polyn + 1.)
!
!--alternatively determine polyk from maximum density
!
 polyk = 1./((polyn + 1.)*densmax**(gamma-1.))*(dfac - 1.)/(2.*dfac)
 print*,' polyk = ',polyk
!
!--rhofac is the factor which multiplies rhofunc to give rho
!  this is 1/B in Owen lingo
!
 rhofac = 1./(polyk*(polyn+1.))**polyn
!
!--work out total mass in torus and set particle mass
!
 Mtorus = 0.
 do iring=1,nrings
    ri = r_in + (iring-1)*deltar
    do iz = 1,nlayers
       zi = (iz - 1)*deltaz - zmax
       Mtorus = Mtorus + 2.*pi*ri*deltar*deltaz*rhofunc(ri,zi,polyn,dfac,Mstar,Rtorus)
!       print*,'ri = ',ri,'zi = ',zi, ' rho = ',rhofac*rhofunc(ri,zi,polyn,dfac,Mstar,Rtorus)
    enddo
 enddo
 Mtorus = Mtorus*rhofac
 print*,'torus mass = ',Mtorus
 if (Mtorus.gt.Mstar) stop 'error Mtorus > Mstar'
 
 massp = Mtorus/real(ntotal)

!
!--setup first ring at r=Rtorus
!
 ipart = 0
 ri = Rtorus
 zi = 0.
 densi = rhofac*rhofunc(ri,zi,polyn,dfac,Mstar,Rtorus)
!--calculate delta r and delta phi from rho
 deltar = (massp/densi)**dndim
 deltaphi = deltar/ri
 npartphi = int(2.*pi/deltaphi)
!--correct deltaphi to integer number of particles
 deltaphi = 2.*pi/npartphi
!--correct deltar to match new deltaphi
 deltar = ri*deltaphi
!--setup particles in z and phi at this r
 call setring(npartphi,ipart,ri,deltar,deltaphi,densi)
 deltar0 = deltar
 dens0 = densi

!
!--setup rings from Rtorus outwards until dens < 0
!
 ri = Rtorus
 do while (densi > tiny(dens))
    zi = 0.
    !--take half step in r using current deltar
    rtemp = ri + 0.5*deltar
    !--get density
    denstemp = rhofac*rhofunc(rtemp,zi,polyn,dfac,Mstar,Rtorus)
    !--calculate delta r and delta phi from rho
    if (denstemp.gt.tiny(denstemp)) then
       deltartemp = (massp/denstemp)**dndim
       !--corrector step on r using midpoint density
       ri = ri + deltartemp
       !--get density
       densi = rhofac*rhofunc(ri,zi,polyn,dfac,Mstar,Rtorus)
    else
       densi = 0.
    endif
    if (densi.gt.tiny(densi)) then
       zi = 0.
       !--get new deltar at new position
       deltar = (massp/densi)**dndim 
       deltaphi = deltar/ri
       npartphi = int(2.*pi/deltaphi)
       !--correct deltaphi to integer number of particles
       deltaphi = 2.*pi/npartphi

       if (ri*deltaphi.lt.2.0*deltartemp) then
          !--correct deltar to match new deltaphi
          deltar = ri*deltaphi
          call setring(npartphi,ipart,ri,deltar,deltaphi,densi)
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
    zi = 0.
    !--take half step in r using current deltar
    rtemp = ri - 0.5*deltar
    !--get density
    denstemp = rhofac*rhofunc(rtemp,zi,polyn,dfac,Mstar,Rtorus)
    if (denstemp.gt.tiny(denstemp)) then
       !--calculate delta r and delta phi from rho
       deltartemp = (massp/denstemp)**dndim
       !--corrector step on r using midpoint density
       ri = ri - deltartemp
       !--get density
       densi = rhofac*rhofunc(ri,zi,polyn,dfac,Mstar,Rtorus)
    else
       densi = 0.
    endif
    if (densi.gt.tiny(densi)) then
       !--get new deltar at new position
       deltar = (massp/densi)**dndim   
       deltaphi = deltar/ri
       npartphi = int(2.*pi/deltaphi)
       !--correct deltaphi to integer number of particles
       deltaphi = 2.*pi/npartphi
       if (ri*deltaphi.lt.2.0*deltartemp) then
          !--correct deltar to match new deltaphi
          deltar = ri*deltaphi
          call setring(npartphi,ipart,ri,deltar,deltaphi,densi)
       else
          deltar = ri*deltaphi
          print*,'skipping ring, too few particles : ',npartphi
       endif
    endif
 enddo
 npart = ipart
 ntotal = ipart
 print*,'Mtorus = ',Mtorus,sum(pmass(1:npart))

! call primitive2conservative
!call set_linklist
! call iterate_density
! call conservative2primitive
 !
 !--balance pressure forces with centripedal acceleration
 !
 vel = 0.
 !--set magnetic field using plasma beta
 beta = 1.e2
 Bzi = sqrt(2.*pr(1)/beta)
 print*,' using beta = ',beta
 if (imhd.ne.0) then
    print*,'initial Bz field = ',Bzi
 else
    Bzi = 0.
 endif
 dbeta = 1./beta
 
 if (iexternal_force.eq.2) then
    print*,'setting v to balance pressure gradients'
!    do i=1,npart
!       rpart = sqrt(dot_product(x(1:2,i),x(1:2,i)))
!       rhat(1:2) = x(1:2,i)/rpart
!       rhat(3) = 0.
!       frad = abs(dot_product(force(1:2,i),rhat(1:2)))
!       omegai = sqrt(frad/rpart)
!       vel(1,i) = -omegai*x(2,i)
!       vel(2,i) = omegai*x(1,i)
!       if (ndimV.ge.3) vel(3,i) = 0.
       !!print*,'forcez = ',force(3,i),force(2,i)
       !v2onr = dot_product(vel(1:2,i),vel(1:2,i))/rpart
       !print*,'v2onr = ',v2onr*rhat(2),force(2,i),v2onr*rhat(1),force(1,i)

!    enddo
    !
    !--analytic velocities
    !
    do i=1,npart
       rcyl2 = DOT_PRODUCT(x(1:2,i),x(1:2,i))
       rcyl = SQRT(rcyl2)
       if (ndim.eq.3) then
          rsph = sqrt(rcyl2 + x(3,i)*x(3,i))
       else
          rsph = rcyl
       endif
       v2onr = 1./(Rtorus)*(-Rtorus*rcyl/rsph**3 + Rtorus**2/(rcyl2*rcyl)) + rcyl/rsph**3
       !--compare to frad from SPH forces
       frad = abs(dot_product(force(1:2,i),x(1:2,i)/rcyl))
       
       omegai = sqrt(v2onr/rcyl)
       vel(1,i) = -omegai*x(2,i)
       vel(2,i) = omegai*x(1,i)
!       call vector_transform(x(1:2,i),vel(1:2,i),2,1,rhat(1:2),2,2)
!       print*,'omegai = ',omegai,rhat(2)
       if (ndimV.ge.3 .and. imhd.ne.0) then
          vel(3,i) = 0.
          densi = dens(i)
          Bfield(3,i) = dbeta*(densi**2/rcyl &
                      + 2.*Mstar/(polyk*gamma)*densi**(3.-gamma) &
                      *(Rtorus/rcyl**3 - rcyl/rsph**3))
          Bfield(1,i) = dbeta*((2.*Mstar/(polyk*gamma))*densi**(3.-gamma) &
                       *(x(3,i)/rsph**3))*x(1,i)/rcyl
          Bfield(2,i) = dbeta*((2.*Mstar/(polyk*gamma))*densi**(3.-gamma) &
                       *(x(3,i)/rsph**3))*x(2,i)/rcyl
!          print*,'B ',i,' = ',Bfield(:,i)
       endif
    enddo
 endif

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

!
!--sets several ring of particles above and below the midplane
!  and around in phi
!
subroutine setring(npartphi,ipart,ri,deltar,deltaphi,densi)
 implicit none
 integer, intent(in) :: npartphi
 integer, intent(inout) :: ipart
 real, intent(in) :: ri,deltar,deltaphi,densi
 real :: deltaz,phii,zi,ztemp,denszi,denstemp,deltaztemp
 integer :: i
 
 deltaz = deltar
 denszi = densi
 zi = 0.
!--upwards from and including the midplane
 do while (denszi > tiny(denszi))
   !--now setup ring of particles
    do i = 1,npartphi
       ipart = ipart + 1
       phii = (i-1)*deltaphi
       x(1,ipart) = ri*COS(phii)
       x(2,ipart) = ri*SIN(phii)
       if (ndim.ge.3) x(3,ipart) = zi
       dens(ipart) = denszi
       pmass(ipart) = massp
       pri = polyk*denszi**gamma
       uu(ipart) = pri/(denszi*(gamma-1.))
    enddo
    !--take half step in r using current deltaz
    ztemp = zi + 0.5*deltaz
    denstemp = rhofac*rhofunc(ri,ztemp,polyn,dfac,Mstar,Rtorus)
    if (denstemp.gt.tiny(denstemp)) then
       deltaztemp = massp/(denstemp*deltar**2)
       zi = zi + deltaztemp
       !--get density
       denszi = rhofac*rhofunc(ri,zi,polyn,dfac,Mstar,Rtorus)
       deltaz = massp/(denstemp*deltar**2)
    else
       denszi = 0.
    endif
    !print*,'zi = ',ri,zi,denszi,ipart
 enddo
 deltaz = deltar
 denszi = densi
 zi = 0.
!--downwards from midplane
 do while (denszi > tiny(denszi))
   !--take half step in r using current deltaz
    ztemp = zi - 0.5*deltaz
    denstemp = rhofac*rhofunc(ri,ztemp,polyn,dfac,Mstar,Rtorus)
    if (denstemp.gt.tiny(denstemp)) then
       deltaztemp = massp/(denstemp*deltar**2)
       zi = zi - deltaztemp
       !--get density
       denszi = rhofac*rhofunc(ri,zi,polyn,dfac,Mstar,Rtorus)
       deltaz = massp/(denstemp*deltar**2)
    else
       denszi = 0.
    endif
    if (denszi > tiny(dens)) then
       !--now setup ring of particles
        do i = 1,npartphi
           ipart = ipart + 1
           phii = (i-1)*deltaphi
           x(1,ipart) = ri*COS(phii)
           x(2,ipart) = ri*SIN(phii)
           if (ndim.ge.3) x(3,ipart) = zi
           dens(ipart) = denszi
           pmass(ipart) = massp
           pri = polyk*denszi**gamma
           uu(ipart) = pri/(denszi*(gamma-1.))
        enddo
    endif
     !print*,'zi = ',ri,zi,denszi,ipart
 enddo

end subroutine setring

end subroutine setup

!----------------------------------------------------
! this subroutine modifies the static configuration
! (from the dumpfile) and gives it the appropriate
! velocity perturbation for the rotating torus
!----------------------------------------------------
subroutine modify_dump
 use dimen_mhd
 use debug
 use loguns
 use options, only:iexternal_force
 use part, only:x,vel,npart
 use rates, only:force
 use timestep, only:time
 implicit none
 integer :: i
 real :: omegai,rcyl2,rcyl,rsph,v2onr
 real, parameter :: Rtorus = 1.0

 
 write(iprint,*) 'MODIFYING INITIAL SETUP with torus rotation '
 if (iexternal_force.ne.2) then
    write(iprint,*) '*** setting external force to 1/r^2 ***'
 endif
 iexternal_force = 2
 time = 0.

 call primitive2conservative
 !
 !--balance pressure forces with centripedal acceleration
 !
! do i=1,npart
!    rpart = sqrt(dot_product(x(1:2,i),x(1:2,i)))
!    rhat(1:2) = x(1:2,i)/rpart
!    rhat(3) = 0.
!    frad = abs(dot_product(force(:,i),rhat(:)))
!    omegai = sqrt(frad/rpart)
!    vel(1,i) = -omegai*x(2,i)  !!/sqrt(rpart)
!    vel(2,i) = omegai*x(1,i) !!!/sqrt(rpart)
!    if (ndimV.ge.3) vel(3,i) = 0.
    !!print*,'forcez = ',force(3,i),force(2,i)
    !!v2onr = dot_product(vel(1:2,i),vel(1:2,i))/rpart
    !!print*,'v2onr = ',v2onr*rhat(2),force(2,i)
! enddo
 !
 !--analytic velocities
 !
 do i=1,npart
    rcyl2 = DOT_PRODUCT(x(1:2,i),x(1:2,i))
    rcyl = SQRT(rcyl2)
    if (ndim.eq.3) then
       rsph = sqrt(rcyl2 + x(3,i)*x(3,i))
    else
       rsph = rcyl
    endif
    v2onr = 1./(Rtorus)*(-Rtorus*rcyl/rsph**3 + Rtorus**2/(rcyl2*rcyl)) + 1./rsph**2
    omegai = sqrt(v2onr/rcyl)
    vel(1,i) = -omegai*x(2,i)
    vel(2,i) = omegai*x(1,i)
    if (ndimV.ge.3) vel(3,i) = 0.
 enddo


end subroutine modify_dump
