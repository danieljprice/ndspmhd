!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Generic setup for 2D shock tests in MHD                               !!
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
 USE setup_params
!
!--define local variables
!      
 IMPLICIT NONE
 INTEGER :: i,j,ntot,npartx,nparty,ipart,iregion,nregions
 REAL :: betazero,denszero,przero,vzero,Bzero,uuzero,machzero
 REAL :: totmass,gam1,massp,const
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) '  Entering subroutine setup'

!
!--check number of dimensions is right
!
 IF (ndim.NE.2) STOP ' ndim must be = 2 for this setup'
!
!--set boundaries
!            	    
 ibound = 3	! periodic
 nbpts = 0	! must use fixed particles if inflow/outflow at boundaries
 xmin(1) = 0.0		! x
 xmax(1) = 1.0
 xmin(2) = 0.0		! y
 xmax(2) = 1.0
 nregions = 1		! number of distinct physical regions
!
!--setup parameters
!
 const = 4.*pi
 betazero = 10./3.
 machzero = 1.0
 vzero = 1.0
 Bzero = 1.0/SQRT(const)
 przero = 0.5*Bzero**2*betazero
 denszero = gamma*przero*machzero
! denszero = 1.0
 
 gam1 = gamma - 1.
 npartx = INT((xmax(1)-xmin(1))/psep)
 nparty = INT((xmax(2)-xmin(2))/psep)
 totmass = denszero*(xmax(2)-xmin(2))*(xmax(1)-xmin(1))
 massp = totmass/FLOAT(npartx*nparty) ! average particle mass
 uuzero = przero/(gam1*denszero)
 
 PRINT*,' totmass, massp, hfact = ',totmass,massp,hfact
 PRINT*,' dens,pr,uu,v,B,m = ',denszero,przero,uuzero,vzero,Bzero,massp
 PRINT*,' npartx, nparty, npart = ',npartx,nparty,npartx*nparty
!
!--allocate memory here
!
 ntot = npartx*nparty
 CALL alloc(ntot)
 npart = ntot
 ntotal = ntot

 DO iregion = 1,nregions
!
!--setup uniform density grid (close packed arrangement)
!
    CALL set_uniform_cartesian(2,psep,xmin,xmax,.false.)
!
!--set particle properties
!
    DO ipart=1,ntotal
       vel(1,ipart) = -vzero*SIN(2.*pi*x(2,ipart))
       vel(2,ipart) = vzero*SIN(2.*pi*x(1,ipart))
       vel(3,ipart) = 0.
       dens(ipart) = denszero
       pmass(ipart) = massp
       uu(ipart) = uuzero
       IF (imhd.GE.1) THEN 
          Bfield(1,ipart) = -Bzero*SIN(2.*pi*x(2,ipart))
          Bfield(2,ipart) = Bzero*SIN(4.*pi*x(1,ipart))
          Bfield(3,ipart) = 0.0	
       ELSE
          Bfield(:,ipart) = 0.
       ENDIF 
    ENDDO
 
 ENDDO
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) '  Exiting subroutine setup'
            
 RETURN
END SUBROUTINE
