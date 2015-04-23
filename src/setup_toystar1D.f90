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
 REAL :: Gn, Pm		! polynomial functions
 LOGICAL :: oscills
 CHARACTER(LEN=40) :: tstarfile
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) ' Entering subroutine setup (toystar)'

 IF (ndim.GT.1) STOP 'Toy star not implemented for ndim > 1'

 ibound = 0	! no boundaries
 nbpts = 0	! no fixed particles
 gam1 = gamma - 1.
 dgam1 = 1./gam1
!
!--toy star parameters (could read from file)
! 
 H = 1.0	! choose H=1 so central density is unity 
 C = 1.0	! C > 0
 A = 0.05/SQRT(2.)
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
 massp = psep*H		! this means that psep is the minimum separation
 imax = 1.1*totmass/massp   ! max number of particles
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
!--start at centre
!
 ipart = 1
 x(1,ipart) = 0.
 dens(ipart) = (H - C*(x(1,ipart))**2)**dgam1
!
!--setup particles x > 0 with 1-x^2 density profile
!
 i = ipart
 DO WHILE (x(1,i).LT.radstar)
    i = i + 1
    x(1,i) = x(1,i-1) + massp/dens(i-1)
    denszero = H - C*x(1,i)**2
    IF (denszero.GE.0) dens(i) = denszero**dgam1
!   print*,i,x(1,i),dens(i),uu(i),pmass(i),x(1,i)-x(1,i-1)
 ENDDO
 ipart = i-1
 IF (dens(ipart).LT.0.01) ipart = ipart-1
 npart = 2*(ipart - 1) + 1
!
!--copy these values to second half of domain (x > 0)
!
 DO i = ipart+1,npart
    j = i - ipart + 1   ! particle to copy from
    x(1,i) = x(1,j)
    dens(i) = dens(j)
 ENDDO 
!
!--copy these values to first half of domain (x > 0)
!
 DO i = 1,ipart-1
    j = npart - i + 1	! particle to copy from
    x(1,i) = -x(1,j)
    dens(i) = dens(j)
 ENDDO 
!
!--write over central particle
!
 x(1,ipart) = 0.
 dens(ipart) = (H - C*(x(1,ipart))**2)**dgam1 
!
!--set particle properties
!
 DO i=1,npart
    IF (oscills) THEN
       vel(1,i) = A*Gn(x(1,i),norder)
    ELSE
       vel(1,i) = A*x(1,i)
    ENDIF 
    IF (ndimV.GT.1) vel(2:3,i) = 0.
    pmass(i) = massp
    uu(i) = polyk*dens(i)**(gam1)/gam1
    IF (imhd.GE.1) THEN 
       Bfield(1,i) = 0.0	!SQRT(1.5)
       Bfield(2,i) = dens(i)*sigma
       Bfield(3,i) = 0.0	
    ELSE
       Bfield(:,i) = 0.
    ENDIF 
 ENDDO
  
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

FUNCTION Pm(x,m)
!
!--calculates a Legendre Polynomial of order m
!
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: m
 REAL, INTENT(IN) :: x
 INTEGER :: i
 REAL :: Pmminus1,Pmminus2,Pm
!
!--specify first two Legendre polynomials
!
 Pmminus2 = 1.
 Pmminus1 = x
 
 SELECT CASE(m)
    CASE (0)
       Pm = 1.
    CASE (1)
       Pm = x
    CASE (2:)	! use recurrence relation to calculate the rest
       DO i=2,m
          Pm = ((2.*(i-1.)+1.)*x*Pmminus1 - (i-1.)*Pmminus2)/REAL(i)
	  Pmminus2 = Pmminus1
	  Pmminus1 = Pm
       ENDDO
 END SELECT
 
END FUNCTION Pm

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
