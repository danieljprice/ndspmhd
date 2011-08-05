!----------------------------------------------------------------
!     Set up a Kelvin-Helmholtz instability test as
!     in Jim Stone et al. (Athena test suite)
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
 use setup_params, only:psep
 use eos, only:gamma
 
 use uniform_distributions
!
!--define local variables
!            
 implicit none
 integer :: i,iseed
 real :: massp,volume,totmass,ran1
 real :: denszero,densmedium,przero,psepmedium
 real, dimension(ndim) :: xminregion,xmaxregion
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(unifdis)'
 if (ndim.ne.2) stop 'error: need ndim=2 for this problem'
!
!--set boundaries
! 	    
 nbpts = 0	! use ghosts not fixed
 xmin(:) = -0.5	! set position of boundaries
 xmax(:) = 0.5
!
!--set up the uniform density grid
!
 denszero = 1.0
 densmedium = 2.0
 przero = 2.5
!
!--setup -0.5 < y < -0.25 and 0.25 < y < 0.5
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

 xminregion(2) = 0.25
 xmaxregion(2) = 0.5
 call set_uniform_cartesian(2,psep,xminregion,xmaxregion,fill=.true.)
!
!--setup -0.25 < y < 0.25
! 
 xminregion(2) = -0.25
 xmaxregion(2) = 0.25
 psepmedium = psep*(denszero/densmedium)**(1./ndim)
 call set_uniform_cartesian(2,psepmedium,xminregion,xmaxregion,fill=.true.)
 
 npart = ntotal
 print*,'npart =',npart
!
!--now declare periodic boundaries
!
 ibound = 3	! boundaries
 iseed = -23864
!
!--now assign particle properties
! 
 do i=1,ntotal
    vel(:,i) = 0.
    if (abs(x(2,i)).lt.0.25) then
       vel(1,i) = 0.5 
       dens(i) = densmedium
    else
       vel(1,i) = -0.5
       dens(i) = denszero
    endif
    pmass(i) = massp
    uu(i) = przero/((gamma-1.)*dens(i))
    Bfield(:,i) = 0.
    Bfield(1,i) = 0.5
!
!--add random velocity perturbation
!
    vel(1,i) = vel(1,i)*(1.0 + 0.01*(ran1(iseed)-0.5))
    vel(2,i) = 0.01*(ran1(iseed)-0.5)
 enddo
!
!--get rho from a sum and then set u to give a
!  smooth pressure
!
 write(iprint,*) 'calling density to make smooth pressure jump...'
 call primitive2conservative
 do i=1,npart
    uu(i) = przero/((gamma - 1.)*rho(i))
 enddo
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
  
 return
end
