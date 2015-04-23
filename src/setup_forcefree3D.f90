!----------------------------------------------------------------
!     Set up for force-free magnetic cylinder in 3D
!     (with Chiara Toniolo)
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
 use eos, only:polyk
 
 use uniform_distributions
!
!--define local variables
!            
 implicit none
 integer :: i
 real :: massp,volume,totmass
 real :: denszero,rmin,rmax,rcyl,Bphi
 real, parameter :: rmu = 2.5
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(unifdis)'
!
!--set boundaries
! 	    
 ibound = 1	! boundaries
! ibound(3) = 3
 nbpts = 0	! use ghosts not fixed
 xmin(:) = -1.0 - 2.*hfact*psep	! set position of boundaries
 xmax(:) = 1.0 + 2.*hfact*psep
 xmin(3) = -1.5- 2.*hfact*psep
 xmax(3) = 1.5+ 2.*hfact*psep
!
!--set up the uniform density grid
!
 rmin = 0.
 rmax = 1.0

 call set_uniform_cartesian(2,psep,xmin,xmax,adjustbound=.true.)
 npart = ntotal
 print*,'npart =',npart
!
!--determine particle mass
!
 denszero = 1.0e-3
 polyk = 1.0
 volume = product(xmax(:)-xmin(:))
 totmass = denszero*volume
 massp = totmass/float(ntotal) ! average particle mass
!
!--now assign particle properties
! 
 do i=1,ntotal
    vel(:,i) = 0.
    dens(i) = denszero
    pmass(i) = massp
    uu(i) = 1.5	! isothermal
    rcyl = sqrt(dot_product(x(1:2,i),x(1:2,i)))
    Bphi = besj1(rmu*rcyl)
    Bfield(1,i) = Bphi*(-x(2,i)/rcyl)
    Bfield(2,i) = Bphi*(x(1,i)/rcyl)
    Bfield(3,i) = besj0(rmu*rcyl)
    if (rcyl.gt.rmax  .or. &
       (x(3,i).gt.(xmax(3)-2.*hfact*psep) .or. x(3,i).lt.(xmin(3)+2.*hfact*psep))) then
       itype(i) = 1
       nbpts = nbpts + 1
    endif
 enddo 
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
