!----------------------------------------------------------------
!  Setup for the Gresho Vortex problem                         !!
!                                                              !!
!  Gresho & Chan (1990), Liska & Wendroff (2003)               !!
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
 use eos, only:gamma
 use mem_allocation
 
 use uniform_distributions
!
!--define local variables
!            
 implicit none
 integer :: i,ir,iphi,nr,ipart,nphi
 real :: massp,volume,totmass,deltaphi
 real :: denszero,r,ri,pri,vphi,phi,rmin,rmax
 logical, parameter :: use_rings = .false.
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(Gresho Vortex)'

 print "(a)",' Setup for Gresho Vortex problem'
!
!--set boundaries
! 	    
 ibound = 3     ! boundaries
 nbpts = 0      ! use ghosts not fixed
 xmin(:) = -0.5   ! set position of boundaries
 xmax(:) = 0.5
 if (use_rings) then
    if (ndim.eq.3) stop 'gresho vortex with rings not implemented in 3D'
    rmax = 0.4
 else
    rmax = 0.
 endif
!
!--set up the uniform density grid
!
 call set_uniform_cartesian(2,psep,xmin,xmax,fill=.true.,rmin=rmax)
 npart = ntotal
 print*,'npart =',npart

 if (use_rings) then
!
!--set up particles on rings
!
    rmin = 0.
    nr = int((rmax - rmin)/psep)

    call alloc(5*npart)

    ipart = npart
    do ir = 1,nr
       ri = rmin + (ir-1)*psep
       deltaphi = psep
       nphi = int((2.*pi*ri)/deltaphi) + 1
       deltaphi = 2.*pi/real(nphi)
       print*,'nphi = ',nphi,npart+ipart

       do iphi = 1,nphi
          phi = (iphi-1)*deltaphi
          ipart = ipart + 1
          if (ipart.gt.size(dens)) call alloc(10*size(dens))
          x(1,ipart) = ri*cos(phi)
          x(2,ipart) = ri*sin(phi)
       enddo
    enddo
    npart = ipart
    ntotal = npart
 endif
!
!--determine particle mass
!
 denszero = 1.0
 print*,' dens = ',denszero
 volume = product(xmax(:)-xmin(:))
 totmass = denszero*volume
 massp = totmass/float(ntotal) ! average particle mass
!
!--now assign particle properties
! 
 do i=1,ntotal
    phi = atan2(x(2,i),x(1,i))
    vel(:,i) = 0.
    r = sqrt(dot_product(x(1:2,i),x(1:2,i)))
    if (r.lt.0.2) then
       vphi = 5.*r
       pri  = 5. + 25./(2.)*r*r
    elseif (r.lt.0.4) then
       vphi = 2. - 5.*r
       pri  = 9. + 25./(2.)*r*r - 20.*r + 4.*log(5.*r)
    else
       vphi = 0.
       pri = 3. + 4.*log(2.)
    endif
    vel(1,i) = -vphi*sin(phi)
    vel(2,i) = vphi*cos(phi)
    !vel(1,i) = 0.01*sin(2.*pi*(x(1,i)-xmin(1)))
    dens(i) = denszero
    pmass(i) = massp
    uu(i) = pri/((gamma - 1.)*denszero)
    Bfield(:,i) = 0.
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
