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
 use part
 use setup_params
 
 use geometry
 use uniform_distributions
!
!--define local variables
!      
 implicit none
 integer :: i,j
 real :: rmax,totmass,totvol,gamm1,rr2
 real :: denszero,uuzero,massp,denscentre,volpart
 real, dimension(ndim) :: xnew
 logical :: iuserings,iequalmass

 write(iprint,*) 'uniform spherical distribution (for toy star)'
 iuserings = .false.
 iequalmass = .false.
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
    call set_uniform_spherical(2,rmax,centred=.true.) 
 endif
!
!--set particle properties
! 
 select case (ndim)
  case(1)
    totvol = 2.*rmax
  case(2)
    totvol = pi*rmax**2
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
 call reset_centre_of_mass(x(:,1:npart),pmass(1:npart))
 
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
 use part
 use toystar2D_utils
 use geometry
 use timestep, only:time
 use setup_params, only:pi
 implicit none
 integer :: i,jmode,smode
 real :: rr,H,C,A,scalefac,sigma2,sigma,rstar,denscentre,gamm1
 real :: omegasq,cs2centre,ekin,ekin_norm
 real, dimension(ndim) :: xcyl,velcyl,dvel
 character(len=len(rootname)+6) :: tstarfile
 logical :: oscills

 time = 0.

 jmode = 2
 smode = 0
 
 write(iprint,*) 'MODIFYING INITIAL SETUP with toystar oscillations'

!
!--read parameters from file
!  
 tstarfile = rootname(1:len_trim(rootname))//'.tstar'
 open(unit=ireadf,err=11,file=tstarfile,status='old',form='formatted')
!!    read(ireadf,*,err=12) H,C,A
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
   if (jmode.lt.0) oscills = .false.

 write(iprint,*) 'radial mode = ',jmode,' theta mode = ',smode
 
 gamm1 = gamma - 1.
 if (gamm1.lt.1.e-5) then
    stop 'error: gamma - 1 <= 0'
 endif

 omegasq = 1.0
 sigma2 = 0.5*omegasq*(gamm1)*((jmode+smode)*(jmode+smode + 2./gamm1) - smode**2)
 if (sigma2.lt.1.e-5) then
    print*,'ERROR sigma2 < 0 in perturbation'
 else
    sigma = sqrt(sigma2)
 endif

 denscentre = 1.0
 scalefac = polyk*gamma/(sigma*gamm1)
 rstar = sqrt((2.*polyk*gamma*denscentre**gamm1)/gamm1)
 
 write(iprint,*) 'polyk = ',polyk,' rstar = ',rstar,' period = ',2.*pi/sigma
 cs2centre = gamma*polyk*denscentre**gamm1
 write(iprint,*) 'denscentre = ',denscentre,' cs_0 = ',sqrt(cs2centre)

!
!--work out frequency of oscillation
!
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
 ekin_norm = (0.05)**2*cs2centre*ekin_norm
 print*,' ekin = ',ekin, ' ekin_norm = ',ekin_norm
 vel = vel*sqrt(ekin_norm/ekin)

 return
end subroutine modify_dump
