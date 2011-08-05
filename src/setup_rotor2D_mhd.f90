!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Setup for the MHD Rotor problems in Toth (2000)                       !!
!!                                                                        !!
!!  dense, rotating disk of fluid at origin                               !!
!!  magnetic field is initially uniform in x-direction                    !!
!!  -> simplest way to set this up is to vary the particle mass           !!
!!                                                                        !!                                                                        !!
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
 INTEGER :: i,j,ntot,npartx,nparty,ipart
 REAL :: rhozero,rhodisk,przero,vzero,ftaper
 REAL :: pri,rdisk,rbuffer,radius
 REAL :: totmass,gam1,massp,const
 REAL, DIMENSION(ndim) :: xorigin, dx
 REAL, DIMENSION(ndimV) :: Bzero
!
!--check number of dimensions is right
!
 IF (ndim.NE.2) STOP ' ndim must be = 2 for MHD rotor problem'
 IF (ndimV.NE.2) WRITE(iprint,*) ' WARNING: best if ndimV=2 for this problem'
!
!--set boundaries
!            	    
 ibound = 2	! reflective ghosts (boundaries not important in this problem)
 nbpts = 0	! no fixed particles
 xmin(:) = -0.5	! unit square
 xmax(:) = 0.5
 const = SQRT(4.*pi) 
!
!--setup parameters for the problem
! 
 xorigin(:) = 0.0	! co-ordinates of the centre of the initial blast
 rdisk = 0.1		! radius of the initial disk
 rbuffer = rdisk	!0.115	! radius of the smoothed front
 vzero = 2.0		! rotation speed of initial disk
 Bzero(:) = 0.
 IF (imhd.NE.0) Bzero(1) = 5.0/const	! uniform field in Bx direction
 przero = 1.0		! initial pressure
 rhozero = 1.0		! ambient density
 rhodisk = 10.0		! density of rotating disk
 
 gam1 = gamma - 1.

 WRITE(iprint,*) 'Two dimensional MHD rotor problem '
 WRITE(iprint,10) rhodisk,rdisk,Bzero(1),vzero,przero
10 FORMAT(/,' Central density  = ',f10.3,', disk radius = ',f6.3,/, &
            ' Initial Bx   = ',f6.3,', rotation = ',f6.3,', pressure = ',f6.3,/)
!
!--setup uniform density grid of particles (2D) 
!  (determines particle number and allocates memory)
!
 CALL set_uniform_cartesian(1,xmin,xmax,.false.)	! 2 = close packed arrangement

 ntotal = npart
!
!--determine particle mass in ambient medium
!
 totmass = rhozero*PRODUCT(xmax(:)-xmin(:))
 massp = totmass/FLOAT(ntotal) ! average particle mass
!
!--now assign particle properties
! 
 DO ipart=1,ntotal
    dx(:) = xin(:,ipart)-xorigin(:) 
    radius = SQRT(DOT_PRODUCT(dx,dx))
    IF (radius.LE.rdisk) THEN
       rhoin(ipart) = rhodisk
       pmass(ipart) = massp*rhodisk/rhozero
       velin(1,ipart) = -vzero*(xin(2,ipart)-xorigin(2))/rdisk
       velin(2,ipart) = vzero*(xin(1,ipart)-xorigin(1))/rdisk
    ELSEIF (radius.LE.rbuffer) THEN	! smooth edge with taper function (Toth)
       ftaper = (rbuffer-radius)/(rbuffer - rdisk)
       rhoin(ipart) = rhozero + (rhodisk-rhozero)*ftaper
       pmass(ipart) = massp*rhoin(ipart)/rhozero
       velin(1,ipart) = -ftaper*vzero*(xin(2,ipart)-xorigin(2))/radius
       velin(2,ipart) = ftaper*vzero*(xin(1,ipart)-xorigin(1))/radius
    ELSE
       pmass(ipart) = massp
       rhoin(ipart) = rhozero
       velin(:,ipart) = 0.
    ENDIF  
    pri = przero 
    uuin(ipart) = pri/(gam1*rhozero)
    hhin(ipart) = hfact*(pmass(ipart)/rhoin(ipart))**hpower	 ! ie constant everywhere
    Bin(:,ipart) = Bzero(:)
 ENDDO
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) '  Exiting subroutine setup'
            
 RETURN
END
