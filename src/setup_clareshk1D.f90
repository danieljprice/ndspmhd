!----------------------------------------------------------------
!     Set up a uniform density cartesian grid of particles in ND
!----------------------------------------------------------------

subroutine setup
!
!--include relevant global variables
!
 use dimen_mhd, only:ndim
 use debug, only:trace
 use loguns, only:iprint
 use bound,only:xmin,xmax,pext
 use options, only:ibound
 use part
 use setup_params
 
 use uniform_distributions
!
!--define local variables
!            
 implicit none
 integer :: i
 real :: massp,volume,totmass
 real :: denszero,spsound,betamhd,theta,Bzero
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(unifdis)'
!
!--set boundaries
! 	    
 ibound = 0	! boundaries
 nbpts = 0	! use ghosts not fixed
 xmin(:) = -1.5	! set position of boundaries
 xmax(:) = 1.5
 spsound = 0.3
 betamhd = 1.0
 theta = pi/4.  ! angle of magnetic field to shock
!
!--set up the uniform density grid
!
 call set_uniform_cartesian(1,psep,xmin,xmax,.false.)
 npart = ntotal
 print*,'npart =',npart
!
!--determine particle mass
!
 !!!denszero = 0.1
 volume = product(xmax(:)-xmin(:))
 totmass = 1.0 !!!denszero*volume
 denszero = totmass/volume
 massp = totmass/float(ntotal) ! average particle mass

! betamhd = pr/0.5*Bzero**2
 Bzero = sqrt(2.*spsound**2*denszero/betamhd)
 Bconst(1) = Bzero*cos(theta)
 Bconst(2) = Bzero*sin(theta)
!
!--now assign particle properties
! 
 do i=1,ntotal
    vel(1,i) = 50.*spsound
    dens(i) = denszero
    pmass(i) = massp
    uu(i) = 1.5*spsound**2	! isothermal
    Bfield(:,i) = Bconst(:)
 enddo
 pext = spsound**2*denszero
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
