!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Velocity field for periodic MHD turbulence simulations in 1D          !!
!!                                                                        !!
!!  Power spectrum of Gaussian random perturbations                       !!
!!  Based partly on a similar code used by M. Bate & Volker Bromm         !!
!!
!!  Power spectrum is normalised according to the total specific
!!  kinetic energy of the particles.
!!                                                                        !!
!!------------------------------------------------------------------------!!

SUBROUTINE set_vperp(xmin,xmax,ekin_in)
!
!--global variables
!
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE part
 USE part_in
 USE random	! random number generators
!
!--local variables
! 
 IMPLICIT NONE
 REAL, INTENT(IN), DIMENSION(ndim) :: xmin,xmax
 REAL, INTENT(IN) :: ekin_in
 REAL, PARAMETER :: pi=3.1415926536
 INTEGER, PARAMETER :: ngrid = 32	! size of grid
 INTEGER, PARAMETER :: nfreq = ngrid 	! number of frequencies
 INTEGER :: i,j,iseed,iseed2,igrid,ifreq
 REAL, DIMENSION(ngrid) :: xgrid	! fixed grid to calculate vperp on
 REAL, DIMENSION(nfreq) :: amplk,phase	! amplitude, phase for each frequency
 REAL, DIMENSION(ngrid) :: vperp 	! velocity field on grid
 REAL :: dxx,dvel 
 REAL :: amplin			! normalisation factor
 REAL :: pindex			! index of power spectrum
 REAL :: sigma			! dispersion of Gaussian
 REAL :: wk,wk_min,wkmod,wkdotx	! wave number
 REAL :: dxgrid,Lbox,twopi
 REAL :: ekin
!
!--allow for tracing flow
! 
 IF (trace) WRITE(iprint,*) 'Entering subroutine set_vperp'
!
!--check number of dimensions and other settings are OK
! 
 IF (ndim.GT.1) STOP 'set_vperp not yet implemented for > 1D'
!
!--set parameters of the turbulent velocity field
!
 pindex = 17./3.		! index of power spectrum
 Lbox = xmax(1)-xmin(1)		! assumes cube in > 1D
 amplin = 1.0			! normalisation factor of power spectrum
 dxgrid = Lbox/REAL(ngrid)	! separation of grid points
!
!--set random seed for random number generator
!  (same seed gives same random number sequence each time)
!
 iseed = -3507		! seed for phases
 iseed2 = -2394 	! seed for amplitudes
!
!--work out the amplitude of the perturbation at each wavenumber (frequency)
!  do this for wavenumbers from kmin = 2*pi/L to kmax = 32*(2*pi/L)
!  The amplitude is a random number drawn from a Gaussian distribution,
!  so that we have e^(-ampl^2/sigma^2) = random, where sigma is the 
!  dispersion of the power spectrum and is proportional to the wavenumber to 
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
!--then choose a random amplitude
!
    amplk(ifreq) = sigma*SQRT(-log(ran1(iseed2)))
    PRINT*,ifreq,amplk(ifreq),sigma

 ENDDO
!
!--having got the amplitude of each frequency, we need to convert this back
!  from k-space to real space, to work out the actual velocity perturbation
!  at each point in space.
!
!  So, for each grid point, we sum over all the frequencies.
!
!
!--to do this fast we would need to insert an inverse fourier transform here
!
 WRITE(iprint,*) ' Constructing velocity field'
 
 DO i=1,ngrid
    xgrid(i) = xmin(1) + (i-1)*dxgrid + 0.5*dxgrid
    DO ifreq=1,nfreq
       wk = REAL(ifreq)*wk_min	! this is the wavenumber (kx,ky,kz)
       wkdotx = wk*(xgrid(i)-xmin(i))
       vperp(i) = vperp(i) + amplk(ifreq)*SIN(wkdotx + phase(ifreq))
    ENDDO
 ENDDO

 WRITE(iprint,*) ' Interpolating to particles...'
!
!--now interpolate from the fixed grid to the SPH particles
!
 ekin = 0.

 DO i=1,npart
!
!--work out nearest grid point
!    
    igrid = INT((xin(1,i) - xmin(1))/dxgrid) + 1 	! closest grid point
    PRINT*,'particle ',i,' grid = ',igrid,xin(1,i)
    dxx = xin(1,i) - xgrid(igrid) 		! distance of particle from pt 
    dvel = (vperp(igrid+1) - vperp(igrid))/dxgrid
!
!--interpolate contribution from neighbouring grid points
!   
    velin(2,i) = velin(2,i) + (vperp(igrid) + dxx*dvel)
    ekin = ekin + 0.5*velin(2,i)**2   
 
 ENDDO
 
 WRITE(iprint,*) ' Kinetic energy = ',ekin,ekin_in
!
!--lastly, normalise the power spectrum using the kinetic energy
! 
 velin(2,:) = velin(2,:)*SQRT(ekin_in/ekin)
 
 IF (trace) WRITE(iprint,*) 'Exiting subroutine set_vperp'
 RETURN
END SUBROUTINE set_vperp
