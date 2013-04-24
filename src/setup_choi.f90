!----------------------------------------------------------------
!     Set up a travelling soundwave in the x direction
!     perturbs particles from a uniform density grid
!     should work in 1, 2 and 3 dimensions
!----------------------------------------------------------------

subroutine setup
!
!--include relevant global variables
!
 use dimen_mhd
 use debug
 use loguns
 use bound
 use eos
 use options
 use part
 use setup_params
 use uniform_distributions, only:set_uniform_cartesian
!
!--define local variables
!            
 implicit none
 integer :: i
 integer, parameter :: itsmax = 100
! real, parameter :: pi = 3.1415926536
 real, parameter :: tol = 1.e-8
 integer :: its,iwave
 real, dimension(ndimV) :: Bzero
 real :: massp
 real :: ampl,wk,xlambda,dxmax,denom
 real :: dxi,dxprev,xmassfrac,func,fderiv
 real :: spsoundi,valfven2i,vamplx,vamply,vamplz
 real :: vfast,vslow,vwave,term,dens1
 real :: denszero,uuzero,przero,Rzero,rhoi

 real :: L, kx, ky, k, k2, xi, yi, D
 real :: Rgamma

 real :: a1, b1, c1
 complex*8 :: aa1, bb1, cc1
 complex*8 :: omega, rootable, root1, root2
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' Entering subroutine setup'
10 format(/,'-------------- ',a,' ----------------')

 iwave = 0                 ! preset wave parameters
 !--isothermal equation of state P=c_s^2*rho
! gamma = 1.

 !--(1-D) is the ionisation fraction
 D = 0.90

 !--Ion-Neutral coupling parameter
 Rgamma = 100.

 !--Crude switch for ambipolar diffusion (on/off)
 iambipolar = .false.!true.


!--set boundaries
!
 ibound(:) = 3        ! periodic boundaries
! ibound(3) = 0        ! outflow boundary (nothing).
 nbpts = 0                ! use ghosts not fixed
 xmin(:) = 0.   ! set position of boundaries
 xmax(1) = 1.0 
 if (ndim.GE.2) then
    xmax(2:ndim) = 1.0!6.*psep ! would need to adjust this depending on grid setup
 endif

!--wave number for perturbations
 L = xmax(1) - xmin(1)
 kx = 2.*pi/L
 ky = 2.*pi/L
 k2 = kx**2 + ky**2
 k = sqrt(k2)
!
!--setup parameters for model

 denszero = 4.
 Bzero(:) = 0.
 Bzero(1) = 1.
 print *, 'Bzero = ', Bzero
 ampl = 0.1
 uuzero = 0.9
!
!--initially set up a uniform density grid (also determines npart)
!  (the call to set_uniform_cartesian means this works in 1,2 and 3D)
!
 call set_uniform_cartesian(2,psep,xmin,xmax,.true.)
 
 if (ndim.eq.1) then
    massp = (xmax(1)-xmin(1))*denszero/FLOAT(npart)
 elseif (ndim.eq.2) then
    massp = (xmax(1)-xmin(1))*(xmax(2)-xmin(2))*denszero/FLOAT(npart)
 elseif (ndim.eq.3) then
    massp = (xmax(1)-xmin(1))*(xmax(2)-xmin(2))*(xmax(3)-xmin(3))* &
         denszero/FLOAT(npart)
 endif

 print *, 'Mass p ', massp
!
!--setup uniform density grid of particles
! 
 do i=1,npart
    vel(:,i) = 0.
    dens(i) = denszero
    pmass(i) = massp
    uu(i) = uuzero
    if (imhd.GT.0) then 
       Bfield(:,i) = Bzero
    else
       Bfield(:,i) = 0.
    endif
 ENDDO

 ntotal = npart
!
!--get sound speed from equation of state (want average sound speed, so
!  before the density is perturbed)
!

 call equation_of_state1(przero,spsoundi,uuzero,denszero)
 print*,' gamma = ',gamma
 print*,' pr = ',przero,' cs = ',spsoundi,' u = ',uu(1),' dens = ',dens(1)
!
!--work out MHD wave speeds
!
 dens1 = 1./denszero
 
 valfven2i = dot_product(Bzero,Bzero)*dens1

 vfast = sqrt(0.5*(spsoundi**2 + valfven2i             &
                 + sqrt((spsoundi**2 + valfven2i)**2   &
                 - 4.*(spsoundi*Bzero(1))**2*dens1)))
 vslow = sqrt(0.5*(spsoundi**2 + valfven2i             &
                 - sqrt((spsoundi**2 + valfven2i)**2   &
                 - 4.*(spsoundi*Bzero(1))**2*dens1)))

 eta_a = 1./(Rgamma*(1.-D))

 print *, 'eta_a = ', eta_a

!----------------------------
! set the wave speed to use
!----------------------------

! vwave = sqrt(valfven2i)
 vwave = 1.0
 if (vwave.le.0.) then
    write(iprint,*) 'Error in setup: vwave = ',vwave
    stop
 endif
  
!
!--multiply by appropriate wave speed
!
 vamplx = vwave*ampl    
 print *, 'V amplitude ', vamplx

 do i = 1, npart
    xi = x(1,i)
    yi = x(2,i)
    vel(3,i) = vamplx*SIN(kx*xi + ky*yi)
 enddo

 !
 !--Write out parameters for analytic form of wave decay
 !
 rhoi = denszero*(1.-D)

 print *, 'sqrt ', valfven2i*k2/(Rgamma**2*rhoi**2)

 !
 !--Calculating directly WITHOUT y-substitution method (they match).
 !
 aa1 = cmplx(1., 0., kind=8)
 bb1 = cmplx(0., valfven2i*k2/(Rgamma*rhoi), kind=8)
 cc1 = cmplx(-valfven2i*k2, 0., kind=8)

 rootable = bb1**2 - 4.*aa1*cc1

 root1 = (-bb1 + sqrt(rootable))/(2.*aa1)
 root2 = (-bb1 - sqrt(rootable))/(2.*aa1)

 print *, 'root 1', root1
 print *, 'root 2', root2
 write(76,*) 'root 1', root1
 write(76,*) 'root 2', root2

 print *, 'square ', root1**2, real(root1)**2 - aimag(root1)**2, 2.*real(root1)*aimag(root1)
 print *, 'linear ', -valfven2i*k2/(Rgamma*rhoi)*aimag(root1), valfven2i*k2/(Rgamma*rhoi)*real(root1)
 print *, 'c ', valfven2i*k2

! root1 = cmplx(5.45, -0.02)

 print *, 'real total ', real(root1)**2 - aimag(root1)**2 -valfven2i*k2/(Rgamma*rhoi)*aimag(root1) - valfven2i*k2
 print *, 'Imag total ', 2.*real(root1)*aimag(root1) + valfven2i*k2/(Rgamma*rhoi)*real(root1)

! stop

! omegaI = 0.5*(valfcen2i*k2/(Rgamma**2*rhoi) + sqrt(

 write(iprint,*) ' Wave set: Amplitude = ',ampl,' Wavelength = ',L,' k = ',k
 write(iprint,*) ' sound speed  = ',spsoundi
 write(iprint,*) ' alfven speed = ',sqrt(valfven2i)
 write(iprint,*) ' fast speed   = ',vfast
 write(iprint,*) ' slow speed   = ',vslow
 write(iprint,*) ' wave speed   = ',vwave,iwave
 
 return
end subroutine setup
