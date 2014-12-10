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
!  Set up the sedimenting dust problem from Monaghan (1997b)
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
 use externf,   only:external_forces
 use cons2prim, only:primitive2conservative 
 use uniform_distributions
 use eos, only:polyk
!
!--define local variables
!            
 implicit none
 integer :: i,ngas,ndust
 real :: massp,volume,totmass,oldtolh,rhodust,masspdust,voldust
 real :: denszero,fext(3),cs,gx,xmindust(ndim),xmaxdust(ndim)
 real :: vdust,vgas,epsi,ts
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(sediment)'
!
!--set boundaries
!
 print*,'Setup for sedimenting dust layer '
 if (iexternal_force.ne.9) stop 'need iexternal_force = 9 for this setup'
 ibound(1) = 3
      ! boundaries
 ibound(2) = 0     ! boundaries
 nbpts = 0      ! use ghosts not fixed
 xmin(:) = 0.   ! set position of boundaries
 xmax(:) = 1.
 xmin(2) = xmin(2) - 8.*psep
 xmax(2) = xmax(2) + 8.*psep
 
 call set_uniform_cartesian(1,psep,xmin,xmax,adjustbound=.true.)
 npart = ntotal
 !xmin(2) = 0.
 !xmax(2) = 1.
 print*,'npart =',npart

!
!--determine particle mass
!
 denszero = 1.0
 volume = product(xmax(:)-xmin(:))
 totmass = denszero*volume
 massp = totmass/float(ntotal) ! average particle mass
!
!--now assign particle properties
!
 cs = 10.
 polyk = cs*cs
 !
 !--get the actual value of the gravity force from the external force routine
 !
 call external_forces(9,x(:,1),fext,ndim,ndimV,vel(:,1),hh(1),cs,itypegas)
 gx = -fext(2)
 itype = itypegas
 
 do i=1,ntotal
    vel(:,i) = 0.
    dens(i) = denszero*exp(-gx*x(2,i)/cs**2)
    pmass(i) = massp*exp(-gx*x(2,i)/cs**2)
    uu(i) = 1.5*polyk ! isothermal
    Bfield(:,i) = 0.
 enddo 

 !
 !--calculate density assuming free boundaries
 !  use very small htol to get it right
 !
 oldtolh = tolh
 tolh = 1.e-12
 call primitive2conservative
 do i=1,ntotal
    dens(i) = pmass(i)/(hh(i)/hfact)**ndim
 enddo
 tolh = oldtolh

 !
 !--now fix boundary particles
 !
 ibound(2) = 1
 do i=1,ntotal
    if (x(2,i) < 0. .or. x(2,i) > 1.) then
       itype(i) = itypebnd
       nbpts = nbpts + 1
    endif
 enddo
 
 !
 !--add dust particles
 !
 if (idust > 0) then
    xmindust = xmin
    xmaxdust = xmax
    xmindust(2) = 0.6
    xmaxdust(2) = 0.8
    rhodust = 0.1 ! in Monaghan 97 this is rhod*thetad = 1000*0.001 = 1.0
    print*,' rhodust = ',rhodust,' rhogas = ',denszero
    voldust = product(xmaxdust-xmindust)
    print*,' assuming Kdrag = ',Kdrag, ' gives vdust = ',-rhodust*gx/Kdrag

    select case(idust)
    case(2)
       ngas = npart
       call set_uniform_cartesian(1,psep,xmindust,xmaxdust,fill=.true.)
       npart = ntotal
       ndust = npart - ngas
       masspdust = rhodust*voldust/real(ndust)
       do i=ngas+1,npart
          itype(i) = itypedust
          pmass(i) = masspdust
          dens(i)  = rhodust
          uu(i)    = 0.
          Bfield(:,i) = 0.
          vel(:,i) = 0.
          vel(2,i) = -rhodust*gx/Kdrag
       enddo
    case default
       ndust = 0
       do i=1,npart
          if (x(2,i) >= xmindust(2) .and. x(2,i) <= xmaxdust(2)) then
             ndust = ndust + 1
          endif
       enddo

       masspdust = rhodust*voldust/real(ndust)
       do i=1,npart
          if (x(2,i) >= xmindust(2) .and. x(2,i) <= xmaxdust(2)) then
             epsi = rhodust/(rhodust + dens(i))
             vdust = -rhodust*gx/Kdrag
             vgas  = 0.
             dustfrac(i) = epsi
             deltav(:,i) = 0.
             deltav(2,i) = vdust - vgas
             vel(2,i) = epsi*vdust + (1. - epsi)*vgas
             !print*,' vy = ',vel(2,i)
             !dens(i)  = dens(i) !+ rhodust
             pmass(i) = pmass(i) + masspdust
          else
             dustfrac(i) = 0.
             deltav(:,i) = 0.
          endif
       enddo
    end select
    ts = (rhodust*denszero)/(Kdrag*(rhodust + denszero))
    print*,' stopping time ts = ',ts
    print*,' resolution criterion is h < ',cs*ts
    if (idust.eq.2 .and. psep > cs*ts) then
       print*,' warning: this is not being met '
       read*
    endif
 endif
 
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
  
 return
end
