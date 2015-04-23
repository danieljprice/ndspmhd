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
!     Set up a wave in 2 or 3 dimensions
!     propagating in some arbitrary direction relative to the 
!     x, y and z axes
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
!
!--define local variables
!            
 implicit none
 integer :: i
 integer, parameter :: itsmax = 100
! real, parameter :: pi = 3.1415926536
 real, parameter :: tol = 1.e-8
 integer :: its
 real, dimension(ndim) :: runit,dxi,dxmax
 real, dimension(ndimv) :: velzero
 real, dimension(ndimb) :: bzero
 real :: massp,totmass,denszero,gam1,uuzero,przero
 real :: anglexy,ampl,wk,xlambda,rmax,denom
 real :: ri,rprev,xmassfrac,func,fderiv
 real :: pri,spsoundi,valfven2i
 real :: vampl_par,vampl_perp,vamplz,vamply,vamplx
 real :: vparallel,vperp,vz
 real :: bparallel,bperp,bz
 real :: bamplx,bamply,bamplz,bampl_perp
 real :: vfast,vslow,vcrap,vwave,term,dens1
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup'
 write(iprint,*) '------------ wave setup ----------------'
!
!--set direction of wave propagation (runit is unit vector in this direction)
!
 anglexy = 30	! angle in degrees x,y plane
! anglez = 45	! angle in degrees z plane
 anglexy = anglexy*pi/180.	! convert to radians

 if (ndim.eq.2) then
    runit(1) = cos(anglexy)
    runit(2) = sin(anglexy)
 else 
    stop 'this wave setup only for 2d'       
!      runit(3) = 0.
 endif
 write(iprint,*) ' runit = ',runit
!
!--read/set wave parameters
! 
 ampl = 0.001
! write (*,*) 'enter amplitude of disturbance'
! read (*,*) ampl
 
 xlambda = 1.0
! write (*,*) 'enter wavelength lambda'
! read (*,*) xlambda
    
 wk = 2.0*pi/xlambda	! 	wave number
 
!
!--set boundaries
! 	    
 ibound = 3	! periodic boundaries
 nbpts = 0	! no fixed particles
 xmin(:) = 0.0	! set position of boundaries
 xmax(:) = 1.0/runit(:)
! print*,'xmin,xmax = ',xmin,xmax
!
!--setup parameters
!
 vperp = 0.1
 vparallel = 0.1
 vz = 0.1
 denszero = 1.0
 przero = 0.1
 bparallel = 1.0
 bperp = 0.1
 bz = 0.1
 uuzero = 0.3
!
!--work out dependent parameters
!
 gam1 = gamma - 1.
! uuzero = przero/(gam1*denszero)
 velzero(1) = vparallel*runit(1) - vperp*runit(2)
 velzero(2) = vparallel*runit(2) + vperp*runit(1)
 velzero(3) = vz
 bzero(1) = bparallel*runit(1) - bperp*runit(2)
 bzero(2) = bparallel*runit(2) + bperp*runit(1)
 bzero(3) = bz
!
!--initially set up a uniform density grid (also determines npart)
!
 call set_uniform_cartesian(2,psep,xmin,xmax,.false.)	! 2 = close packed
!
!--determine particle mass
!
 dxmax(:) = xmax(:) - xmin(:)
 totmass = denszero*product(dxmax)
 massp = totmass/float(npart) ! average particle mass
 print*,'npart,massp = ',npart,massp
 
 do i=1,npart
    vel(:,i) = velzero
    dens(i) = denszero
    pmass(i) = massp
    uu(i) = uuzero
    if (imhd.ge.1) then 
       bfield(:,i) = bzero
    endif 
 enddo

 ntotal = npart
  
 rmax = dot_product(dxmax,runit)
! print*,'rmax = ',rmax
 denom = rmax - ampl/wk*(cos(wk*rmax)-1.0)
!
!--get sound speed from equation of state (want average sound speed, so
!  before the density is perturbed)
!
 call equation_of_state(przero,spsoundi,uuzero,denszero)
!
!--work out mhd wave speeds
!
 dens1 = 1./denszero
 
 valfven2i = bparallel**2*dens1
 vfast = sqrt(0.5*(spsoundi**2 + valfven2i			&
                 + sqrt((spsoundi**2 + valfven2i)**2	&
                 - 4.*(spsoundi*bparallel)**2*dens1)))
 vslow = sqrt(0.5*(spsoundi**2 + valfven2i			&
                 - sqrt((spsoundi**2 + valfven2i)**2	&
                 - 4.*(spsoundi*bparallel)**2*dens1)))
 vcrap = sqrt(spsoundi**2 + valfven2i)

!----------------------------
! set the wave speed to use
!----------------------------
		    
 vwave = vfast

!----------------------------

!
!--now perturb the particles appropriately
!    
 do i=1,npart
 
    dxi(:) = x(:,i) - xmin(:)
    ri = dot_product(dxi,runit)

    rprev = 2.*rmax
    xmassfrac = ri/rmax	! current mass fraction(for uniform density)
				! determines where particle should be
!
!--use rootfinder on the integrated density perturbation
!  to find the new position of the particle
!    
    its = 0
	 
    do while ((abs(ri-rprev).gt.tol).and.(its.lt.itsmax))
       rprev = ri
       func = xmassfrac*denom - (ri - ampl/wk*(cos(wk*ri)-1.0))
       fderiv = -1.0 - ampl*sin(wk*ri)
       ri = ri - func/fderiv	! newton-raphson iteration
       its = its + 1 
!      print*,'iteration',its,'ri =',ri 
    enddo
	 	 	 
    if (its.ge.itsmax) then
       write(iprint,*) 'error: soundwave - too many iterations'
       call quit 
    endif
	  
!    if (idebug(1:5).eq.'sound') then
!       write(*,99002) i,its,ri-rprev,x(1,i),xmin(1)+ri*runit(1)
99002  format('particle',i5,' converged in ',i4,	&
       ' iterations, error in x =',1(1pe10.2,1x),/,	&
       'previous r = ',1(0pf8.5,1x),			&
       'moved to r = ',1(0pf8.5,1x)) 
!    endif
    x(:,i) = xmin(:) + ri*runit(:)
!
!--multiply by the appropriate amplitudes
!
    dens1 = 1./dens(i)
    term = 1./(vwave**2 - bparallel**2*dens1)
    vampl_par = vwave*ampl
    vampl_perp = -vampl_par*bparallel*bperp*dens1*term
    vamplz = -vampl_par*bparallel*bfield(3,i)*dens1*term
    
    vamplx = vampl_par*runit(1) - vampl_perp*runit(2)
    vamply = vampl_par*runit(2) + vampl_perp*runit(1)

    vel(1,i) = vamplx*sin(wk*ri)
    vel(2,i) = vamply*sin(wk*ri)
    vel(3,i) = vamplz*sin(wk*ri)
!
!--perturb internal energy if not using a polytropic equation of state 
!  (do this before density is perturbed)
!
   uu(i) = uu(i) + pri/dens(i)*ampl*sin(wk*ri)	! if not polytropic
!    
!--perturb density if not using summation
!
    bampl_perp = vwave*bperp*vampl_par*term*sin(wk*ri)
    bamplz = vwave*bz*vampl_par*term*sin(wk*ri)
    
    bamplx = -bampl_perp*runit(2)
    bamply = bampl_perp*runit(1)
    
    bfield(1,i) = bzero(1) + bamplx*sin(wk*ri)
    bfield(2,i) = bzero(2) + bamply*sin(wk*ri)
    bfield(3,i) = bzero(3) + bamplz*sin(wk*ri)
    dens(i) = dens(i)*(1.+ampl*sin(wk*ri))

 enddo

 write(iprint,*) ' wave set: amplitude = ',ampl,' wavelength = ',xlambda,' k = ',wk
 write(iprint,*) ' sound speed  = ',spsoundi
 write(iprint,*) ' alfven speed = ',sqrt(valfven2i)
 write(iprint,*) ' fast speed   = ',vfast
 write(iprint,*) ' slow speed   = ',vslow
 write(iprint,*) ' wave speed   = ',vwave
 
 return
end

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
