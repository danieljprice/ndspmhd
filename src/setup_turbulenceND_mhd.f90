!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Setup for periodic MHD turbulence simulations in 1, 2 and 3D          !!
!!                                                                        !!
!!  Density initially uniform, then call routine to put on a power        !!
!!  spectrum of Gaussian random perturbations                             !!
!!                                                                        !!
!!  Note for all MHD setups, only the magnetic field should be setup      !!
!!  Similarly the thermal energy is setup even if using total energy.     !!
!!                                                                        !!
!!------------------------------------------------------------------------!!

SUBROUTINE setup
!
!--include relevant global variables
!
 USE dimen_mhd
 USE debug
 USE loguns
 USE bound
 USE eos
 USE options
 USE part
 USE setup_params
!
!--define local variables
!      
 IMPLICIT NONE
 INTEGER :: i,j,ntot,ipart,nfreq,iseed1,iseed2
 REAL :: denszero,przero,spsound2
 REAL, DIMENSION(ndimV) :: velzero,Bzero
 REAL :: totmass,volume,massp
 REAL :: rbuffer, exx,Ekin,ekinpart,pindex
 REAL :: Ewave,Emag,beta,beta_gammie,length
 REAL :: vrms,Brms
!
!--check number of dimensions is right
!
 IF (ndim.EQ.3) STOP ' code not yet working for ndim=3'
 IF (ndimV.NE.3) STOP ' Need ndimV=3 for this problem'
!! IF (iener.NE.0) STOP ' eos should be isothermal for this problem (iener=0)'
!
!--set boundaries
!            	    
 ibound = 3	! periodic ghosts
 nbpts = 0	! must use fixed particles if inflow/outflow at boundaries
 xmin(:) = -0.5
 xmax(:) = 0.5
 length = 1.0	! length of box in each direction
!
!--setup parameters for the problem (could read these from file)
! 
 beta_gammie = 0.01	! Gammie and Ostriker use a funny beta
 beta = 0.5*beta_gammie	! this is the usual beta (ratio of pressures)

 denszero = 1.0		! initial density
 velzero(:) = 0.
 velzero(1) = 0.	! initial velocities
 przero = 1.0		! initial pressure 
 spsound2 = przero/denszero
 Bzero(:) = 0.
 Bzero(1) = SQRT(przero/beta)
 Bconst(:) = Bzero(:)
 
 Ewave = 100*denszero*length*spsound2		! initial wave energy
 Ekin = 0.5*Ewave				! initial kinetic energy
 
 WRITE(iprint,10) ndim
10 FORMAT(/,1x,i1,'D compressible MHD turbulence ')
 WRITE(iprint,20) beta_gammie, Ewave, Bzero(1), denszero, przero
20 FORMAT(/,' Beta (gammie)  = ',f10.3,', wave energy = ',f10.3,/, &
            ' Initial Bx     = ',f6.3,', density = ',f6.3,', pressure = ',f6.3,/)
!
!--setup uniform density grid of particles (2D) 
!  (determines particle number and allocates memory)
!
 CALL set_uniform_cartesian(1,psep,xmin,xmax,.false.)	! 2 = close packed arrangement

 ntotal = npart
!
!--determine particle mass
!
 volume = PRODUCT(xmax(:)-xmin(:))
 totmass = denszero*volume
 massp = totmass/FLOAT(ntotal) ! average particle mass
!
!--now assign particle properties
! 
 DO ipart=1,ntotal
    vel(:,ipart) = velzero(:)
    dens(ipart) = denszero
    pmass(ipart) = massp
    uu(ipart) = 1.0	! isothermal
    Bfield(:,ipart) = Bzero(:)
 ENDDO
!
!--now call routine to set up the velocity field
!
 pindex = -2  ! P(k) = k**pindex
 nfreq = 32   ! wavelengths down to lambda/32
 iseed1 = -3507
 iseed2 = -2394
 CALL set_powerspec1D(x(1,:),vel(2,:),npart,xmin(1),xmax(1), &
                      pindex,nfreq,iseed1,iseed2)
!
!--normalise according to kinetic energy
!
 ekinpart = 0.
 DO i = 1,npart
    ekinpart = ekinpart + 0.5*pmass(i)*DOT_PRODUCT(vel(:,i),vel(:,i))
 ENDDO
 vel(2,1:npart) = vel(2,1:npart)*sqrt(Ekin/ekinpart)
!
!
!--set up magnetic field
!
 iseed1 = -732
 iseed2 = -8973
 CALL set_powerspec1D(x(1,:),Bfield(2,:),npart,xmin(1),xmax(1), &
                      pindex,nfreq,iseed1,iseed2)
!
!--normalise this using magnetic energy
!
 ekinpart = 0.
 DO i = 1,npart
    ekinpart = ekinpart + 0.5*pmass(i)*(Bfield(2,i)**2)
 ENDDO
 Bfield(2,1:npart) = Bfield(2,1:npart)*sqrt(Ekin/ekinpart)

!
!--work out rms velocity and mag field perturbations
!
 Brms = 0.
 vrms = 0.
 Ewave = 0.
 Ekin = 0.
 Emag = 0.
 DO i=1,npart
    Ekin = Ekin + 0.5*pmass(i)*DOT_PRODUCT(vel(:,i),vel(:,i))
    Emag = Emag + 0.5*pmass(i)*Bfield(2,i)**2/denszero
    Brms = Brms + SQRT((Bfield(2,i)/Bzero(1))**2)
    vrms = vrms + SQRT(vel(2,i)**2/spsound2)
 ENDDO
 Ewave = Ekin + Emag
 vrms = vrms/REAL(npart)
 Brms = Brms/REAL(npart)
 WRITE(iprint,*) ' rms velocity perturbation = ',vrms,' mag field = ',Brms
 WRITE(iprint,*) ' Ekin = ',Ekin,' Emag = ',Emag,' Ewave = ',Ewave
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) '  Exiting subroutine setup'
            
 RETURN
END
