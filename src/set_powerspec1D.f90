!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Velocity field for periodic MHD turbulence simulations in 1D          !!
!!                                                                        !!
!!  Power spectrum of Gaussian random perturbations                       !!
!!                                                                        !!
!!  Power spectrum is normalised according to the total specific          !!
!!  kinetic energy of the particles.                                      !!
!!                                                                        !!
!!  Input  : x(npart)   : particle positions
!!           npart      : number of particles
!!           xmin, xmax : box dimensions
!!           pindex     : index of power spectrum
!!           nfreq      : number of frequencies to set
!!
!!  Output : vel(npart) 1D velocity field (note this is NOT normalised)
!!------------------------------------------------------------------------!!

SUBROUTINE set_powerspec1D(x,vel,npart,xmin,xmax,pindex,nfreq,iseedin,iseed2in)
!
!--global variables
!
 USE debug
 USE loguns
 use random, only:ran1,rayleigh_deviate
!
!--local variables
! 
 IMPLICIT NONE
 REAL, PARAMETER :: pi=3.1415926536
 INTEGER, INTENT(IN) :: npart, nfreq, iseedin, iseed2in
 REAL, INTENT(IN), DIMENSION(npart) :: x
 REAL, INTENT(OUT), DIMENSION(npart) :: vel
 REAL, INTENT(IN) :: xmin,xmax,pindex
 INTEGER :: i,j,iseed,iseed2,igrid,ifreq
 REAL, DIMENSION(nfreq) :: amplk,phase	! amplitude, phase for each frequency
 REAL :: dxx,dvel 
 REAL :: amplin			! normalisation factor
 REAL :: pindex			! index of power spectrum
 REAL :: sigma			! dispersion of Gaussian
 REAL :: wk,wk_min,wkmod,wkdotx	! wave number
 REAL :: Lbox,twopi
 REAL :: ekin
!
!--allow for tracing flow
! 
 IF (trace) WRITE(iprint,*) 'Entering subroutine set_powerspec1D'
!
!--set parameters of the turbulent velocity field
!
 Lbox = xmax-xmin		! assumes cube in > 1D
 amplin = 1.0			! normalisation factor of power spectrum
!
!--set random seed for random number generator
!  (same seed gives same random number sequence each time)
!
 iseed = iseedin		! seed for phases
 iseed2 = iseed2in 	! seed for amplitudes
!
!--work out the amplitude of the perturbation at each wavenumber (frequency)
!  do this for wavenumbers from kmin = 2*pi/L to kmax = nfreq*(2*pi/L)
!  The amplitude is a random number drawn from a Gaussian distribution,
!  where sigma is the dispersion of the power spectrum and is proportional to the wavenumber to 
!  some power.
!
 twopi = 2.*pi
 wk_min = twopi/Lbox	! so division is only done once

 WRITE(iprint,*) ' Setting amplitudes and phases'

 DO ifreq = 1,nfreq		! looping from max wavelength to min wavelength
    wk = REAL(ifreq)*wk_min	! this is the wavenumber
    wkmod = wk
    sigma = SQRT(amplin*(wkmod**pindex))
!
!--choose a random phase for each frequency (betw. -pi and pi)
!    
    phase(ifreq) = -pi + twopi*ran1(iseed)
!
!--then choose a random amplitude from a gaussian with dispersion sigma
!  (since we are drawing the modulus of the amplitude, 
!   this comes from a Rayleigh distribution)
!
    amplk(ifreq) = sigma*rayleigh_deviate(iseed2)
    PRINT*,ifreq,amplk(ifreq),sigma,phase(ifreq)

 ENDDO

 WRITE(iprint,*) ' Constructing velocity field'
!
!--now interpolate from the fixed grid to the SPH particles
!
 ekin = 0.

 DO i=1,npart
!
!--having got the amplitude of each frequency, we need to convert this back
!  from k-space to real space, to work out the actual velocity perturbation
!  at each point in space.
!
!  So, for each point in space, we sum over all the frequencies.
!
!--to do this fast we would need to insert an inverse fourier transform here
!  and interpolate from a grid
!
    DO ifreq=1,nfreq
       wk = REAL(ifreq)*wk_min	! this is the wavenumber (kx,ky,kz)
       wkdotx = wk*(x(i)-xmin)
       vel(i) = vel(i) + amplk(ifreq)*SIN(wkdotx + phase(ifreq))
    ENDDO

    ekin = ekin + 0.5*vel(i)**2   
 
 ENDDO
 
 WRITE(iprint,*) ' Kinetic energy = ',ekin
!
!--lastly, normalise the power spectrum using the kinetic energy
! 
! vel(:) = vel(:)*SQRT(ekin_in/ekin)
 
 IF (trace) WRITE(iprint,*) 'Exiting subroutine set_powerspec1D'
 RETURN
END SUBROUTINE set_powerspec1D
