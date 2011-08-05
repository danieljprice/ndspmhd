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
 USE part_in
 USE setup_params
!
!--define local variables
!      
 IMPLICIT NONE
 INTEGER :: i,j,ntot,ipart
 REAL :: rhozero,przero,spsound2
 REAL, DIMENSION(ndimV) :: velzero,Bzero
 REAL :: totmass,volume,massp,const
 REAL :: rbuffer, exx,Ekin
 REAL :: Ewave,beta,beta_gammie,length
!
!--check number of dimensions is right
!
 IF (ndim.EQ.3) STOP ' code not yet working for ndim=3'
 IF (ndimV.NE.3) STOP ' Need ndimV=3 for this problem'
 IF (iener.NE.0) STOP ' eos should be isothermal for this problem (iener=0)'
!
!--set boundaries
!            	    
 ibound = 3	! periodic ghosts
 nbpts = 0	! must use fixed particles if inflow/outflow at boundaries
 xmin(:) = -0.5
 xmax(:) = 0.5
 length = 1.0	! length of box in each direction
 const = 4.*pi 
!
!--setup parameters for the problem (could read these from file)
! 
 beta_gammie = 0.01	! Gammie and Ostriker use a funny beta
 beta = 2.*beta_gammie	! this is the usual beta (ratio of pressures)

 rhozero = 1.0		! initial density
 velzero(1) = 0.	! initial velocities
 przero = 1.0		! initial pressure 
 spsound2 = przero/rhozero
 Bzero(1) = SQRT(przero/beta)
 
 Ewave = 100*rhozero*length*spsound2		! initial wave energy
 Ekin = 0.5*Ewave				! initial kinetic energy
 
 WRITE(iprint,10) ndim
10 FORMAT(/,1x,i1,'D compressible MHD turbulence ')
 WRITE(iprint,20) beta_gammie, Ewave, Bzero(1), rhozero, przero
20 FORMAT(/,' Beta (gammie)  = ',f10.3,', wave energy = ',f10.3,/, &
            ' Initial Bx     = ',f6.3,', density = ',f6.3,', pressure = ',f6.3,/)
!
!--setup uniform density grid of particles (2D) 
!  (determines particle number and allocates memory)
!
 CALL set_uniform_cartesian(1,xmin,xmax,.false.)	! 2 = close packed arrangement

 ntotal = npart
!
!--determine particle mass
!
 volume = PRODUCT(xmax(:)-xmin(:))
 totmass = rhozero*volume
 massp = totmass/FLOAT(ntotal) ! average particle mass
!
!--now assign particle properties
! 
 DO ipart=1,ntotal
    velin(:,ipart) = velzero(:)
    rhoin(ipart) = rhozero
    pmass(ipart) = massp
    uuin(ipart) = 1.0	! isothermal
    hhin(ipart) = hfact*(massp/rhoin(ipart))**hpower	 ! ie constant everywhere
    Bin(:,ipart) = Bzero(:)
 ENDDO
!
!--now call routine to do the Gaussian random perturbation
!
 CALL set_vperp(xmin,xmax,Ekin)
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) '  Exiting subroutine setup'
            
 RETURN
END
