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
 real :: rmax,totmass,totvol
 real :: denszero,uuzero,massp,denscentre
 real, dimension(ndim) :: xnew
 logical :: iuserings

 write(iprint,*) 'uniform spherical distribution (for toy star)'
 iuserings = .false.
!
!--set bounds of initial setup
!                   
 denscentre = 1.0                ! toy star central density 
!
!--reset polyk to give r = 1
!
 polyk = (gamma - 1.)/(2.*gamma*denscentre**(gamma-1.))
 write(iprint,*) 'resetting polyk = ',polyk
 rmax = 1.0
 totmass = pi*rmax**2*(gamma-1.)/gamma
 ibound = 0        ! no boundaries 
 iexternal_force = 1        ! use toy star force
!
!--setup a uniform sphere of particles
! 
 if (iuserings) then
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
 else
    call set_uniform_spherical(1,rmax,perturb=0.5)        ! 4 = random
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
 pmass(:) = massp
 Bfield(:,:) = 0.
!
!--reset centre of mass to zero
!
 call reset_centre_of_mass(x(:,1:npart),pmass(1:npart))
 
 return
end subroutine setup

subroutine modify_dump
 implicit none

 return
end subroutine modify_dump
