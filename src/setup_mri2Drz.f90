!----------------------------------------------------------------
!     Set up for the 2D (r-z) MRI test problem described in 
!     Hawley & Balbus, ApJ 376, 223-233 (1991)
!----------------------------------------------------------------

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
 use mem_allocation, only:alloc
!
!--define local variables
!            
 implicit none
 integer :: i,j,ipart,npartr,npartz
 real :: massp,volume,totmass,deltar,deltaz,deltarav,rpos,zpos
 real :: denszero,uuzero,cs0,polyk0,Rcentre
 real :: przero,Bzeroz,asize,betamhd
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(unifdis)'
 write(iprint,*) '2D magneto-rotational instability'
 if (ndim.ne.2) stop 'error: this is a 2D problem and ndim.ne.2'
 if (ndimV.ne.3) stop 'error: we need ndimV=3 for this problem'
!
!--geometry is cylindrical r-z
!
 geom = 'cylrzp'
 geomsetup = 'cylrzp'
!
!--set position of disc patch for coriolis & centrifugal forces
!
 Rcentre = 100.
 Omega = (1./Rcentre**1.5)
 Omega2 = (1./Rcentre**3)
 iexternal_force = 2 ! 1/r^2 force
!
!--set boundaries
!
 ibound(1) = 2  ! reflecting in x (==R)  
 ibound(2) = 3	! periodic in y   (==z)
 nbpts = 0	! use ghosts not fixed
 asize = 1.0    ! box size (corresponds to a in Hawley/Balbus 1992)
 xmin(1) = Rcentre-asize  ! set position of boundaries
 xmax(1) = Rcentre+asize
 xmin(2) = -0.5*asize
 xmax(2) = 0.5*asize
!
!--set up the uniform density grid
!  NB: uniform density in r-z means r spacing decreases linearly
!
 npartr = int((xmax(1)-xmin(1))/psep)    
 npartz = int((xmax(2)-xmin(2))/psep)
 deltarav = (xmax(1)-xmin(1))/npartr
 deltaz = (xmax(2)-xmin(2))/npartz
 print*,'deltarav,deltaz = ',deltarav,deltaz
 print*,'npartr, npartz = ',npartr,npartz
 npart = npartr*npartz
 call alloc(npart)
     
 ipart = 0
 do j=1,npartz
    zpos = xmin(2) + (j-0.5)*deltaz
    rpos = xmin(1)
    do i=1,npartr
       ipart = ipart + 1
       deltar = deltarav*Rcentre/rpos
       !print*,'deltar = ',deltar, ' at r = ',rpos
       rpos = rpos + 0.5*deltar
       x(1,ipart) = rpos
       x(2,ipart) = zpos
       rpos = rpos + 0.5*deltar
    enddo
 enddo
 !--adjust outer r boundary to fall halfway between lattice points
 xmax(1) = rpos
 print*,ipart,npart
! call set_uniform_cartesian(11,psep,xmin,xmax,.false.)
 ntotal = npart
!
!--determine particle mass
!
 denszero = 1.0
!--volume is difference between 2 circles times z height
!  utterly no idea where the factor of 2 comes from
 volume = 0.5*(xmax(1)**2 - xmin(1)**2)*(xmax(2)-xmin(2)) 
 totmass = denszero*volume
 massp = totmass/float(ntotal) ! average particle mass
!
!--sound speed
!
 przero = 1.e-5
 cs0 = gamma*przero/denszero
 uuzero = przero/(denszero*(gamma-1.))
 polyk0 = przero/denszero**gamma
 write(iprint,*) ' cs0 = ',cs0, ' u = ',uuzero
 write(iprint,*) ' polyk = ',polyk,' should be = ',polyk0
 polyk = polyk0

 write(iprint,*) ' Omega_c = ',Omega, ' R_c = ',Rcentre

!
!--field strength (beta = pr/(0.5*B^2))
!
 betamhd = 1000.
 Bzeroz = sqrt(2.*przero/betamhd)
 write(iprint,*) ' beta(MHD) = ',betamhd,' Bz = ',Bzeroz,' P = ',przero
!
!--now assign particle properties
! 
 do i=1,ntotal
    vel(:,i) = 0.
    dens(i) = denszero
    pmass(i) = massp
    uu(i) = uuzero ! isothermal
    Bfield(:,i) = 0.
    !
    !--set constant z field in region between -a/5 -> a/5
    !
    if (x(1,i).gt.(Rcentre-0.2*asize) .and. x(1,i).lt.(Rcentre+0.2*asize)) then
       Bfield(2,i) = Bzeroz
    endif
    vel(3,i) = 1./x(1,i)**1.5 !!Omega
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
