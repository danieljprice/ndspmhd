!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Generic setup for the 2D Orszag-Tang vortex test in MHD               !!
!!                                                                        !!
!!  Gives particles a variable separation so the mass per SPH particle    !!
!!  is constant. Shock is smoothed slightly (simple smoothing).           !!
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
 INTEGER :: i,j,ntot,npartx,nparty,ipart
 REAL :: betazero,rhozero,przero,vzero,Bzero,uuzero,machzero
 REAL :: totmass,gam1,massp,const
!
!--check number of dimensions is right
!
 IF (ndim.NE.2) STOP ' ndim must be = 2 for Orszag-Tang vortex'
!
!--set boundaries
!            	    
 ibound = 3	! periodic
 nbpts = 0	! must use fixed particles if inflow/outflow at boundaries
 xmin(1) = 0.0		! x
 xmax(1) = 1.0
 xmin(2) = 0.0		! y
 xmax(2) = 1.0
!
!--setup parameters
!
 const = 4.*pi
 betazero = 10./3.
 machzero = 1.0
 vzero = 1.0
 Bzero = 1.0/SQRT(const)
 przero = 0.5*Bzero**2*betazero
 rhozero = gamma*przero*machzero
! rhozero = 1.0
 
 gam1 = gamma - 1.
! npartx = INT((xmax(1)-xmin(1))/psep)
! nparty = INT((xmax(2)-xmin(2))/psep)
 uuzero = przero/(gam1*rhozero)

 WRITE(iprint,*) 'Two dimensional Orszag-Tang vortex problem '
 WRITE(iprint,10) betazero,machzero,Bzero,rhozero,przero
10 FORMAT(/,' beta        = ',f6.3,', mach number = ',f6.3,/, &
            ' Initial B   = ',f6.3,', density = ',f6.3,', pressure = ',f6.3,/)
!
!--allocate memory here
!
! ntot = npartx*nparty
! CALL alloc(ntot,1)
! npart = ntot
! ntotal = ntot
!
!--setup uniform density grid of particles (2D) with sinusoidal field/velocity
!  determines particle number and allocates memory
 CALL set_uniform_cartesian(1,xmin,xmax,.false.)	! 2 = close packed arrangement

 ntotal = npart
!
!--determine particle mass
!
 totmass = rhozero*(xmax(2)-xmin(2))*(xmax(1)-xmin(1))
 massp = totmass/FLOAT(ntotal) ! average particle mass
!
!--now assign particle properties
! 
 DO ipart=1,ntotal
    velin(1,ipart) = -vzero*SIN(2.*pi*xin(2,ipart))
    velin(2,ipart) = vzero*SIN(2.*pi*xin(1,ipart))
    IF (ndimV.EQ.3) velin(3,ipart) = 0.
    rhoin(ipart) = rhozero
    pmass(ipart) = massp
    uuin(ipart) = uuzero
    hhin(ipart) = hfact*(massp/rhoin(ipart))**hpower	 ! ie constant everywhere
    IF (imhd.GE.1) THEN 
       Bin(1,ipart) = -Bzero*SIN(2.*pi*xin(2,ipart))
       Bin(2,ipart) = Bzero*SIN(4.*pi*xin(1,ipart))
       IF (ndimV.EQ.3) Bin(3,ipart) = 0.0	
    ELSE
       Bin(:,ipart) = 0.
    ENDIF 
!       print*,ipart,xin(:,ipart),rhoin(ipart),uuin(ipart),pmass(ipart)
!       IF (MOD(i,1000).EQ.0) read*
 ENDDO

!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) '  Exiting subroutine setup'
            
 RETURN
END
