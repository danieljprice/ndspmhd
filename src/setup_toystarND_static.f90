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

!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  setup for the toy star static solution in 1, 2 or 3 dimensions        !!
!!                                                                        !!
!!  sets up a uniform spherical distribution in 1, 2 or 3 dimensions,     !!
!!  which should then be damped to the static toy star solution           !!
!!                                                                        !!
!!------------------------------------------------------------------------!!

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
 use part, only:npart,x,vel,dens,uu,Bfield,pmass
 use setup_params
 
 use geometry
 use uniform_distributions
!
!--define local variables
!      
 implicit none
 integer :: i
 real :: rmax,totmass,totvol,gamm1,rr2
 real :: denszero,uuzero,massp,denscentre,volpart !!,aa,Jmag
 real, dimension(ndim) :: xnew
 logical :: iuserings,iequalmass

 write(iprint,*) 'uniform spherical distribution (for toy star)'
 iuserings = .false.
 iequalmass = .true.
!
!--set bounds of initial setup
!                   
 denscentre = 1.0                ! toy star central density 
!
!--reset polyk to give r = 1
!
 gamm1 = gamma - 1.
 polyk = gamm1/(2.*gamma*denscentre**gamm1)
 write(iprint,*) 'resetting polyk = ',polyk
 rmax = 1.0
 totmass = pi*rmax**2*gamm1/gamma
 ibound = 0        ! no boundaries 
 iexternal_force = 1        ! use toy star force
!
!--setup a uniform sphere of particles
! 
 if (iuserings .and. iequalmass) then
    xmin(1) = 2.*psep
    xmax(1) = 1.0
    if (ndim.ge.2) then 
       ibound(2) = 3 ! periodic in phi
       xmin(2) = -pi ! phi min
       xmax(2) = pi ! phi max
    endif
    call set_uniform_cartesian(1,psep,xmin,xmax,perturb=0.5)
    do i=1,npart
       call coord_transform(x(:,i),ndim,2,xnew(:),ndim,1)
       x(:,i) = xnew(:)
    enddo
 elseif (iequalmass) then
    call set_uniform_spherical(1,rmax,perturb=0.5)        ! 4 = random
 else
    call set_uniform_spherical(2,rmax,centred=.true.,trim=0.5*psep) 
 endif
!
!--set particle properties
! 
 select case (ndim)
  case(1)
    totvol = 2.*rmax
  case(2)
    if (iequalmass) then
       totvol = pi*rmax**2
    else
       totvol = pi*(rmax)**2
    endif
  case(3)
    totvol = 4./3.*pi*rmax**3
 end select

 denszero = totmass/totvol        ! initial density
 massp = totmass/real(npart)
 volpart  = totvol/real(npart)
 print*,' volpart = ',volpart
! volpart = psep**2
! print*,' new one = ',volpart
 uuzero = 0.1
 write(iprint,10) denscentre,totmass
10 format(/,' Toy star static solution ',/, &
            '     central density: ',f7.3, ' total mass = ',f7.3,/)
!
!--set these for all particles
! 
 vel(:,:) = 0.
 dens(:) = denszero
 uu(:) = uuzero
 if (iequalmass) then
    pmass(:) = massp
 else
    !--assign variable masses (assumes uniform distribution)
    do i=1,npart
       rr2 = dot_product(x(:,i),x(:,i))
       dens(i) = denscentre*(1.-rr2)**(1./gamm1)
       pmass(i) = dens(i)*volpart
    enddo
 endif
 Bfield(:,:) = 0.
!
!--reset centre of mass to zero
!
 if (iequalmass) call reset_centre_of_mass(x(:,1:npart),pmass(1:npart))
 
! if (imhd.ne.0) then
!    xnew(:) = 0.
!    xnew(2) = 1. 
!    aa = 0.5
!    Jmag = 1.0
!    call set_dipole(aa,Jmag,xnew,x(:,1:npart),Bfield(1:ndim,1:npart),ndim,npart)
! endif
 
 return
end subroutine setup

!----------------------------------------------------
! this subroutine modifies the static configuration
! (from the dumpfile) and gives it the appropriate
! velocity perturbation for the toy star
!
! NB: modification is done in cartesian coords
!----------------------------------------------------
subroutine modify_dump
 use dimen_mhd
 use debug
 use loguns
 use eos, only:gamma,polyk
 use part, only:npart,x,vel,pmass
 use geometry
 use timestep, only:time
 use setup_params, only:pi
 implicit none
 integer :: i,ierr,jmode,smode
 real :: Ctstar,Atstar,scalefac,sigma2,sigma,rstar,denscentre,gamm1
 real :: omegasq,cs2centre,ekin,ekin_norm,amplitude,alpha,betatstar
 real :: ctstar1,ctstar2
 real, dimension(ndim) :: xcyl,velcyl,dvel
 character(len=len(rootname)+6) :: tstarfile
 character(len=30) :: dummy
 logical :: oscills,symmetric

 time = 0.
 amplitude = 0.05
 alpha = 1.0
 betatstar = 2.*pi

 jmode = 2
 smode = 0
 
 write(iprint,*) 'MODIFYING INITIAL SETUP with toystar oscillations'

!
!--read parameters from file
!  
 tstarfile = rootname(1:len_trim(rootname))//'.tstar2D'
 open(unit=ireadf,err=11,file=tstarfile,status='old',form='formatted')
    read(ireadf,*,err=12) dummy
    read(ireadf,*,err=12) alpha, betatstar, ctstar1,ctstar2
    read(ireadf,*,err=12) jmode,smode
 close(unit=ireadf)
 oscills = .true.
   goto 13
11 continue
   write(iprint,*) tstarfile,' not found, using default options '
   goto 13
12 continue
   write(iprint,*) ' error reading ',tstarfile
13 continue   
   symmetric = .true.
   if (jmode.lt.0) oscills = .false.
   if (smode.lt.0) symmetric = .false.

 if (oscills) then
    write(iprint,*) 'radial mode = ',jmode,' theta mode = ',smode
 else
    write(iprint,*) 'NONLINEAR MODES: alpha = ',alpha,' beta = ',betatstar
    if (.not.symmetric) write(iprint,*) '               : c = ',ctstar1,' d = ',ctstar2   
 endif
 
 gamm1 = gamma - 1.
 if (gamm1.lt.1.e-5) then
    stop 'error: gamma - 1 <= 0'
 endif
!
!--work out frequency of oscillation
!
 if (oscills) then
    omegasq = 1.0
    sigma2 = 0.5*omegasq*(gamm1)*((jmode+smode)*(jmode+smode + 2./gamm1) - smode**2)
    if (sigma2.lt.1.e-5) then
       print*,'ERROR sigma2 < 0 in perturbation'
    else
       sigma = sqrt(sigma2)
    endif
 endif

 denscentre = 1.0
 Ctstar = 1.0
 scalefac = polyk*gamma/(sigma*gamm1)
 rstar = sqrt((2.*polyk*gamma*denscentre**gamm1)/gamm1)
 
 write(iprint,*) 'polyk = ',polyk,' rstar = ',rstar,' period = ',2.*pi/sigma
 cs2centre = gamma*polyk*denscentre**gamm1
 write(iprint,*) 'denscentre = ',denscentre,' cs_0 = ',sqrt(cs2centre)

!
!--set velocity perturbation
!
 if (oscills) then
    ekin = 0.
    ekin_norm = 0.
    do i=1,npart
       !--get r,theta
       call coord_transform(x(:,i),ndim,1,xcyl(:),ndim,2)

       !--set v_r
       velcyl(1) = scalefac*detadr(jmode,smode,xcyl(1)/rstar,gamma)*COS(smode*xcyl(2))
       !--set theta_dot
       velcyl(2) = -scalefac*etar(jmode,smode,xcyl(1)/rstar,gamma)*smode*SIN(smode*xcyl(2))/xcyl(1)**2
       !!print*,'v_phi = ',velcyl(2),xcyl(2),etar(jmode,smode,xcyl(1)/rstar,gamma)
       !--now transform back to get vx, vy
       call vector_transform(xcyl(1:ndim),velcyl(1:ndim),ndim,2,dvel(1:ndim),ndim,1)
       if (xcyl(1).lt.1.e-5) then
          print*,' r = 0 on particle ',i,' xcyl(1) = ',xcyl(1), &
                 ' v_cyl = ',velcyl,' v_cart = ',dvel
       endif
       !--now perturb v with appropriate amplitude
       vel(1:ndim,i) = dvel(1:ndim)
       ekin = ekin + 0.5*pmass(i)*dot_product(vel(1:ndim,i),vel(1:ndim,i))
       ekin_norm = ekin_norm + 0.5*pmass(i)
    enddo
   !
   !--normalise the amplitude
   !
    ekin_norm = (amplitude)**2*cs2centre*ekin_norm
    write(iprint,*) ' ekin = ',ekin, ' ekin_norm = ',ekin_norm
    vel = vel*sqrt(ekin_norm/ekin)

    Atstar = scalefac*sqrt(ekin_norm/ekin)
    write(iprint,*) ' v = ',Atstar,'*detadr(r)'
 else
    if (symmetric) then
       do i=1,npart
          !--get r, theta
          call coord_transform(x(:,i),ndim,1,xcyl(:),ndim,2)
          !--set v_r
          velcyl(1) = alpha*xcyl(1)
          !--set theta_dot
          velcyl(2) = betatstar
          !--now transform back to get vx, vy
          call vector_transform(xcyl(1:ndim),velcyl(1:ndim),ndim,2,dvel(1:ndim),ndim,1)
          !--now perturb v with appropriate amplitude
          vel(1:ndim,i) = dvel(1:ndim)
       enddo
    else
       do i=1,npart
          vel(1,i) = alpha*x(1,i) + ctstar1*x(2,i)
          vel(2,i) = ctstar2*x(1,i) + betatstar*x(2,i)
       enddo
    endif
 endif
 
!
!--rewrite the tstar2D file giving the amplitude
!  
 tstarfile = rootname(1:len_trim(rootname))//'.tstar2D'
 write(iprint,*) ' writing to file ',trim(tstarfile)
 open(unit=ireadf,iostat=ierr,file=tstarfile,status='replace',form='formatted')
 if (ierr.eq.0) then
    write(ireadf,*,iostat=ierr) denscentre,Ctstar,Atstar
    write(ireadf,*,iostat=ierr) alpha,betatstar,ctstar1,ctstar2
    write(ireadf,*,iostat=ierr) jmode,smode
    if (ierr /= 0) write(iprint,*) 'ERROR WRITING TO ',trim(tstarfile)
    close(unit=ireadf)
 else
    write(iprint,*) 'ERROR OPENING ',trim(tstarfile)
 endif

 return

contains

!
!--function that evaluates the polynomial for rho(r/re) for a given radial mode
!  (from the power series solution to the 2nd order ODE)
!
!  rad = r/r_star
!  j = radial (axisymmetric) mode
!  m = theta mode 
!
!  solution is for delta(rho**(gamma-1))
!  ie. rho**(gamma-1) = rho_0**(gamma-1) + etar
!
!  and takes the form
!
!  etar = rad**m sum_k a_k rad**k
!
real function etar(j,m,rad,gamma)
  implicit none 
  integer :: j,m,k,kprev   ! j is the radial mode, m is the theta mode
  real :: rad,gamma,denom
  real :: ak,akprev,gamm1,freqsq
!
!--this solution is for arbitrary gamma
!
  gamm1 = gamma - 1.
  if (gamm1.lt.1.e-3) then
     print*,'error gamma -1 <= 0'
     etar = 0.
     return
  endif
!
!--the solution is of the form
!  drhor = a_0 + a_2 (r/re)**2 + a_4 (r/re)**4 + ...
!  where for j = k, coefficients >= a_k+2 are zero
!  
  freqsq = (j+m)*(j+m + 2./gamm1) - m**2

  akprev = 1.0  ! this is a_0 which is the amplitude
  etar = akprev
  !!print*,'mode = ',j,m,' nu^2 = ',freqsq,' a_0 = ',akprev
!
!--the co-efficients for the terms above a_0 are calculated using
!  the recurrence relation between the a_k's
!
  do k = 2,j,2
     kprev = k-2
     denom = real((kprev + 2 + m)**2 - m**2)
     ak = akprev*(kprev**2 + 2.*kprev*m + 2.*(kprev+m)/gamm1 - freqsq)/denom
     !!print*,'coeff ',k,' = ',ak,k**2,2.*k/gamm1
     etar = etar + ak*rad**k
     akprev = ak
  enddo
  
  etar = etar * rad**m

end function etar

!
!--function that evaluates the polynomial for v(r/re) for a given radial mode
!  (from the power series solution to the 2nd order ODE)
!
real function detadr(j,m,rad,gamma)
  implicit none
  integer :: j,m,k,kprev   ! j is the radial mode, m is the theta mode
  real :: rad,gamma,denom,term1,term2
  real :: ak,akprev,gamm1,freqsq
!
!--this solution is for arbitrary gamma
!
  gamm1 = gamma - 1.
  if (gamm1.lt.1.e-3) then
     print*,'error gamma -1 <= 0'
     detadr = 0.
     return
  endif
!
!--the solution is of the form
!  drhor = a_0 + a_2 (r/re)**2 + a_4 (r/re)**4 + ...
!  where for j = k, coefficients >= a_k+2 are zero
!  
  freqsq = (j+m)*(j+m + 2./gamm1) - m**2

  detadr = 0.
  akprev = 1.0  ! this is a_0 which is the amplitude
  term1 = akprev
  term2 = 0.
!  print*,'mode = ',j,m,' nu^2 = ',freqsq,' a_0 = ',akprev
!
!--the co-efficients for the terms above a_0 are calculated using
!  the recurrence relation between the a_k's
!
  do k = 2,j,2
     kprev = k-2
     denom = real((kprev + 2 + m)**2 - m**2)
     ak = akprev*(kprev**2 + 2.*kprev*m + 2.*(kprev+m)/gamm1 - freqsq)/denom
     !!print*,'coeff ',k,' = ',ak,k*ak,rad,(k-1)
     term1 = term1 + ak*rad**k
     term2 = term2 + k*ak*rad**(k-1)
     akprev = ak
  enddo
  
  if (m.eq.0) then
     detadr = term2
  else
     detadr = m*rad**(m-1)*term1 + rad**m*term2
  endif
  
end function detadr

end subroutine modify_dump
