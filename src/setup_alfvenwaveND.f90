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
 USE part_in
 USE setup_params
!
!--define local variables
!            
 IMPLICIT NONE
 INTEGER :: i
! REAL, PARAMETER :: pi = 3.1415926536
 REAL, DIMENSION(ndim) :: runit
 REAL :: massp,totmass,rhozero,gam1,uuzero,przero
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
 anglexy = 30.	! angle in degrees x,y plane
! anglez = 45.	! angle in degrees z plane
 anglexy = anglexy*pi/180.	! convert to radians
 runit(1) = COS(anglexy)
 runit(2) = SIN(anglexy)
! runit(3) = 0.
 WRITE(iprint,*) ' runit = ',runit
!
!--set boundaries
! 	    
 ibound = 3	! periodic boundaries
 nbpts = 0	! no fixed particles
 xmin(:) = 0.0	! set position of boundaries
 xmax(:) = 1.0/runit(:)
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
 rhozero = 1.0
 przero = 0.1
 Bparallel = 1.0
 Bperp0 = 0.1 
 Bz0 = 0.1
!
!--work out dependent parameters
!
 gam1 = gamma - 1.
 uuzero = przero/(gam1*rhozero)
!
!--initially set up a uniform density grid (also determines npart)
!
 PRINT*,' setting up uniform density grid'
 CALL set_uniform_cartesian(2,xmin,xmax,.false.)	! 2 = close packed
!
!--determine particle mass
!
 totmass = rhozero*(xmax(2)-xmin(2))*(xmax(1)-xmin(1))
 massp = totmass/FLOAT(npart) ! average particle mass
 PRINT*,'npart,massp = ',npart,massp
 
 DO i=1,npart        
    perturb_sin = SIN(wk*DOT_PRODUCT(xin(:,i),runit))
    perturb_cos = COS(wk*DOT_PRODUCT(xin(:,i),runit))
    vperp = vperp0*perturb_sin
    vz = vz0*perturb_cos
    Bperp = Bperp0*perturb_sin
    Bz = Bz0*perturb_cos
    
    
    velin(1,i) = vparallel*runit(1) - vperp*runit(2)
    velin(2,i) = vparallel*runit(2) + vperp*runit(1)
    velin(3,i) = vz
    
    rhoin(i) = rhozero
    pmass(i) = massp
!
!--perturb internal energy if not using a polytropic equation of state 
!  (do this before density is perturbed)
!
    uuin(i) = uuzero !+ pri/rhoin(i)*ampl*SIN(wk*ri)	! if not polytropic

    hhin(i) = hfact*(pmass(i)/rhoin(i))**hpower	 ! ie constant everywhere
    IF (imhd.GE.1) THEN 
       Bin(1,i) = Bparallel*runit(1) - Bperp*runit(2)
       Bin(2,i) = Bparallel*runit(2) + Bperp*runit(1)
       Bin(3,i) = Bz
    ENDIF 
 ENDDO

 ntotal = npart
 
 valfven = SQRT(Bparallel**2/rhozero)

!----------------------------

 WRITE(iprint,*) ' Wave set: Amplitude = ',ampl,' Wavelength = ',xlambda,' k = ',wk
 WRITE(iprint,*) ' alfven speed = ',valfven
 
 RETURN
END
