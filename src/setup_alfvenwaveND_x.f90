!----------------------------------------------------------------
!     Set up an Alfven wave in 2 or 3 dimensions
!     density is constant so this is easy
!     >> should be compiled with ndimB = 3
!----------------------------------------------------------------

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
 USE polyconst
 USE setup_params
!
!--define local variables
!            
 IMPLICIT NONE
 INTEGER :: i
! REAL, PARAMETER :: pi = 3.1415926536
 REAL, DIMENSION(3) :: rvec
 REAL, DIMENSION(ndim) :: runit
 REAL :: massp,totmass,denszero,gam1,uuzero,przero
 REAL :: anglexy,ampl,wk,xlambda,rmax
 REAL :: valfven
 REAL :: vampl_par,vampl_perp,vamplz,vamply,vamplx
 REAL :: vparallel,vperp,vz,vperp0,vz0
 REAL :: Bparallel,Bperp,Bz,Bperp0,Bz0
 REAL :: perturb_sin, perturb_cos
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) ' Entering subroutine setup'
!
!--set direction of wave propagation (runit is unit vector in this direction)
!
! anglexy = 30.	! angle in degrees x,y plane
! anglez = 45.	! angle in degrees z plane
! anglexy = anglexy*pi/180.	! convert to radians
 rvec(1) = 1.0	        !0.5*SQRT(3.) !COS(anglexy)
 rvec(2) = 0.
 rvec(3) = 0.
 IF (ndim.EQ.1) THEN ! in 1D must be in x direction
    rvec(1) = 1.0
    rvec(2) = 0.	!0.5	!SIN(anglexy)
 ENDIF   

 runit(1:ndim) = rvec(1:ndim) 

 WRITE(iprint,*) ' runit = ',runit
!
!--set boundaries
! 	    
 ibound = 3	! periodic boundaries
 nbpts = 0	! no fixed particles
 xmin(:) = 0.0	! set position of boundaries
 xmax(:) = 1.0
 if (ndim.ge.2) xmax(2:ndim) = 0.5		!/runit(:)
 PRINT*,'xmin,xmax = ',xmin,xmax
!
!--read/set wave parameters
! 
 rmax = DOT_PRODUCT((xmax(:)-xmin(:)),runit)
 ampl = 0.001
! WRITE (*,*) 'Enter amplitude of disturbance'
! READ (*,*) ampl
 
 xlambda = 1.0	!/COS(anglexy)	!*rmax
! WRITE (*,*) 'Enter wavelength lambda'
! READ (*,*) xlambda
    
 wk = 2.0*pi/xlambda	! 	wave number


!
!--setup parameters
!
 vperp0 = 0.1
 vparallel = 0.0
 vz0 = 0.1
 denszero = 1.0
 przero = 0.1
 Bparallel = 1.0
 Bperp0 = 0.1 
 Bz0 = 0.1
!
!--work out dependent parameters
!
 gam1 = gamma - 1.
 uuzero = przero/(gam1*denszero)
 polyk = przero/(denszero**gamma)	! override setting in input
!
!--initially set up a uniform density grid (also determines npart)
!
 PRINT*,' setting up uniform density grid'
 CALL set_uniform_cartesian(2,psep,xmin,xmax,.false.)	! 2 = close packed
!
!--determine particle mass
!
 totmass = denszero*PRODUCT(xmax(:)-xmin(:))
 massp = totmass/FLOAT(npart) ! average particle mass
 PRINT*,'npart,massp = ',npart,massp
 
 DO i=1,npart        
    perturb_sin = SIN(wk*DOT_PRODUCT(x(:,i),runit))
    perturb_cos = COS(wk*DOT_PRODUCT(x(:,i),runit))
    vperp = vperp0*perturb_sin
    vz = vz0*perturb_cos
    Bperp = Bperp0*perturb_sin
    Bz = Bz0*perturb_cos
    
    
    vel(1,i) = vparallel*rvec(1) - vperp*rvec(2)
    vel(2,i) = vparallel*rvec(2) + vperp*rvec(1)
    vel(3,i) = vz
    
    dens(i) = denszero
    pmass(i) = massp
!
!--perturb internal energy if not using a polytropic equation of state 
!  (do this before density is perturbed)
!
    uu(i) = uuzero !+ pri/dens(i)*ampl*SIN(wk*ri)	! if not polytropic

    hh(i) = hfact*(pmass(i)/dens(i))**hpower	 ! ie constant everywhere
    IF (imhd.GE.1) THEN 
       Bfield(1,i) = Bparallel*rvec(1) - Bperp*rvec(2)
       Bfield(2,i) = Bparallel*rvec(2) + Bperp*rvec(1)
       Bfield(3,i) = Bz
    ENDIF 
 ENDDO

 ntotal = npart
 
 valfven = SQRT(Bparallel**2/denszero)

!----------------------------

 WRITE(iprint,*) ' Wave set: Amplitude = ',ampl,' Wavelength = ',xlambda,' k = ',wk
 WRITE(iprint,*) ' alfven speed = ',valfven
 
 RETURN
END
