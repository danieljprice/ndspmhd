!------------------------------------------------------------------------------!
! NDSPMHD: A Smoothed Particle (Magneto)Hydrodynamics code for (astrophysical) !
! fluid dynamics simulations in 1, 2 and 3 spatial dimensions.                 !
!                                                                              !
! (c) 2002-2014 Daniel Price                                                   !
!                                                                              !
! http://users.monash.edu.au/~dprice/ndspmhd                                   !
! daniel.price@monash.edu -or- dprice@cantab.net (forwards to current address) !
!                                                                              !
!  NDSPMHD comes with ABSOLUTELY NO WARRANTY.                                  !
!  This is free software; and you are welcome to redistribute                  !
!  it under the terms of the GNU General Public License                        !
!  (see LICENSE file for details) and the provision that                       !
!  this notice remains intact. If you modify this file, please                 !
!  note section 2a) of the GPLv2 states that:                                  !
!                                                                              !
!  a) You must cause the modified files to carry prominent notices             !
!     stating that you changed the files and the date of any change.           !
!                                                                              !
!  ChangeLog:                                                                  !
!------------------------------------------------------------------------------!

!!----------------------------------------------------------------
!!     Set up a toy star (with oscillations) in 1D
!!     Starts with particles on a uniform density grid
!!----------------------------------------------------------------

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
 INTEGER :: i,j
 INTEGER, PARAMETER :: itsmax = 100
 REAL, PARAMETER :: tol = 1.e-7
 INTEGER :: imax,its,ipart,norder
 REAL :: massp,totmass,radstar,sigma,denszero
 REAL :: ampl,wk,xlambda,dxmax,denom,gam1,dgam1
 REAL :: H,C,A		! toy star parameters
 REAL :: Gn
 LOGICAL :: oscills
 CHARACTER(LEN=40) :: tstarfile
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) ' Entering subroutine setup (toystar)'

 IF (ndim.GT.1) STOP 'Toy star not implemented for ndim > 1'
!
!--initially set up a uniform density grid
! 	    
 ibound = 0	! no boundaries
 nbpts = 0	! no fixed particles
 gam1 = gamma - 1.
 dgam1 = 1./gam1
!
!--toy star parameters
! 
 H = 1.0	! choose H=1 so central density is unity 
 C = 1.0	! C > 0
 A = 0.05/SQRT(2.)	! 
 sigma = 1./SQRT(2.)
 oscills = .true.	! with oscillations ?
 norder = 3 		! mode of oscillation
 
!
!--read parameters from file
! 
 tstarfile = rootname(1:LEN_TRIM(rootname))//'.tstar'
 OPEN(UNIT=ireadf,ERR=11,FILE=tstarfile,STATUS='old',FORM='formatted')
    READ(ireadf,*,ERR=12) H,C,A
    READ(ireadf,*,ERR=12) sigma
    READ(ireadf,*,ERR=12) norder
 CLOSE(unit=ireadf)
 oscills = .true.
   GOTO 13
11 CONTINUE
   WRITE(iprint,*) tstarfile,' not found, using default options '
   GOTO 13
12 CONTINUE
   WRITE(iprint,*) ' error reading ',tstarfile
13 CONTINUE   
   IF (norder.LT.0) oscills = .false.
!
!--set dependent parameters
! 
 totmass = 4./3.*SQRT(H**3/C) 
 radstar = SQRT(H/C)
 xmin(1) = -radstar	! set position of boundaries
 xmax(1) = radstar
 imax = INT((xmax(1)-xmin(1))/psep)
 massp = psep
 gam1 = gamma - 1.
 iexternal_force = 1
!
!--allocate memory here
!
 CALL alloc(imax)
 
 WRITE(iprint,*) ' Toy star : total mass      = ',totmass
 WRITE(iprint,*) '            toy star radius = ',radstar
 WRITE(iprint,*) '            particle mass   = ',massp
 WRITE(iprint,*) '            A = ',A,' C = ',C,' H = ',H
 IF (oscills) WRITE(iprint,*) ' -> with oscills, order = ',norder 
 WRITE(iprint,*)
!
!--setup uniform density distribution
!
 DO i=1,imax
    x(1,i) = xmin(1) + (i-1)*psep + 0.5*psep
    dens(i) = (H - C*x(1,i)**2)**dgam1    
    IF (oscills) THEN
       vel(1,i) = A*Gn(x(1,i),norder)
    ELSE
       vel(1,i) = A*x(1,i)
    ENDIF 
    IF (ndimV.GT.1) vel(2:3,i) = 0.   
    pmass(i) = psep*dens(i)
    uu(i) = polyk*dens(i)**(gam1)/gam1
    IF (imhd.GE.1) THEN 
       Bfield(1,i) = 0.0	!SQRT(1.5)
       Bfield(2,i) = sigma*dens(i)
       Bfield(3,i) = 0.0	
    ELSE
       Bfield(:,i) = 0.
    ENDIF 
!    PRINT*,x(1,i),dens(i),pmass(i),hh(i)
 ENDDO
 
 npart = imax  
 WRITE(iprint,*) ' finished toy star setup, npart = ',npart
 ntotal = npart
 
 RETURN
END
!
!--function to evaluate the Gegenbauer polynomial of index n given x
!  (the parameter lambda = 3/2 in this case)
!
FUNCTION Gn(x,n)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: n
 REAL, INTENT(IN) :: x
 INTEGER :: i
 REAL :: Gn,Gnminus1,Gnminus2
 REAL :: fnorm
 
 fnorm = 2.*(n+1)*(n+2)/REAL(2.*n + 3.)
! PRINT*,' fnorm = ',fnorm
!
!--specify first two Gegenbauer polynomials
!
 Gnminus2 = 1.
 Gnminus1 = 3.*x
 
 SELECT CASE (n)
    CASE (0) 
       Gn = Gnminus2
    CASE (1)
       Gn = Gnminus1
    CASE (2:)
       DO i=2,n
          Gn = ((2*i+1)*x*Gnminus1 - (i+1)*Gnminus2)/REAL(i)
	  Gnminus2 = Gnminus1
	  Gnminus1 = Gn
       ENDDO
 END SELECT
 
 Gn = Gn/fnorm
END FUNCTION Gn

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
