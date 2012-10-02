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
 use mem_allocation, only:alloc
!
!--define local variables
!            
 implicit none
 integer :: i
 integer, parameter :: itsmax = 100
 integer, parameter :: ntypes = 1
! real, parameter :: pi = 3.1415926536
 real, parameter :: tol = 1.e-8
 integer :: its,iwave,ngas,ndust,jtype
 real, dimension(ndimV) :: Bzero
 real :: massp,masspdust
 real :: ampl,wk,xlambda,dxmax,denom
 real :: dxi,dxprev,xmassfrac,func,fderiv
 real :: spsoundi,valfven2i,vamplx,vamply,vamplz
 real :: vfast,vslow,vcrap,vwave,term,dens1
 real :: denszero,uuzero,przero,Rzero,denszerodust
 real :: dust_to_gas_ratio, ffac

 real :: L,kx,kz,Bigkx,Bigkz,tausmode
 real :: Rux,Iux,Ruy,Iuy,Ruz,Iuz,Rwx,Iwx,Rwy,Iwy,Rwz,Iwz
 real :: Rrhog,Irhog,Rrhod,Irhod
 real :: onepluseps,denomNSH,vgxNSH,vgyNSH,vdxNSH,vdyNSH
 real :: xi,zi,cokx0,sikx0,cokz0,fx0
 real :: cokx,sikx,cokz,sikz, eta, totmass
 logical :: iperturb = .true.

 eta = 0.005
! eta = 1.0

!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' Entering subroutine setup'
10 format(/,'-------------- ',a,' ----------------')

 iwave = 0                 ! preset wave parameters
!
!--setup parameters (could read in from a file)
!
 Bzero = 0.
 if (iwave.EQ.1) then        ! MHD slow wave
    ampl = 0.006
    denszero = 1.0
    uuzero = 4.5
    if (imhd.ne.0) Bzero(1:3) = sqrt(2.)
    write(iprint,10) 'MHD slow wave'
 elseif (iwave.EQ.2) then    ! MHD fast wave
    ampl = 0.0055
    denszero = 1.0
    uuzero = 0.3
    if (imhd.ne.0) Bzero(1:3) = 0.5  
    write(iprint,10) 'MHD fast wave'
 else                        ! Sound wave
    ampl = 100!--0.01
    denszero = 1.0
    if (abs(gamma-1.).gt.1e-3) then
       uuzero = 1.0/((gamma-1.)*gamma)
    else
       uuzero = 1.0
    endif
    Rzero = -1.0  !  negative stress parameter
    open(unit=20,err=30,file=trim(rootname)//'.rstress',status='old')
      read(20,*) Rzero
    close(unit=20)
    print*,'Read stress file: R parameter = ',Rzero
30  continue
    Bzero(1) = sqrt(2.*(1.-Rzero))
    write(iprint,10) 'sound wave'
 endif
 xlambda = 1.0    
 wk = 2.0*pi/xlambda        !         wave number
!
!--set boundaries
!
 ibound = 3        ! periodic boundaries
 nbpts = 0                ! use ghosts not fixed
 xmin(:) = -0.5   ! set position of boundaries
 xmax(1) = 0.5 
 if (ndim.GE.2) then
!    xmax(2:ndim) = 11.*psep ! would need to adjust this depending on grid setup
    xmax(2:ndim) = xmax(1)
 endif
!
!--initially set up a uniform density grid (also determines npart)
!  (the call to set_uniform_cartesian means this works in 1,2 and 3D)
!
 ngas = 0
 ndust = 0
 do jtype=1,ntypes
    call set_uniform_cartesian(1,psep,xmin,xmax,adjustbound=.true.)
    if (jtype.eq.1) then
       ngas = npart
       itype(1:ngas) = itypegas
    elseif (jtype.eq.2) then
       ndust = npart - ngas
       itype(ngas+1:ngas+ndust) = itypedust
    endif
 enddo
 
! npart = INT((xmax(1)-xmin(1))/psep)
! massp = 1.0/FLOAT(ngas)        ! average particle mass
 L           = (xmax(1) - xmin(1))!*eta/Bigkx
 totmass     =     denszero*L*L
 massp       =     totmass/FLOAT(ngas)
 kx          = 2.*pi/L
 kz          = 2.*pi/L
 
 masspdust = 0.
 dust_to_gas_ratio = 1.
 if (ndust.gt.0) masspdust = dust_to_gas_ratio*1.0/FLOAT(ndust) ! average particle mass
 denszerodust = dust_to_gas_ratio*denszero
 if (ntypes.gt.1) print*,' ngas = ',ngas,' ndust = ',ndust
!
!--allocate memory here
!
 call alloc(npart)
!
!--setup uniform density grid of particles
! 
 do i=1,npart
!    x(1,i) = xmin(1) + (i-1)*psep  + 0.5*psep 
    vel(:,i) = 0.
    if (itype(i).eq.itypedust) then
       dens(i) = denszerodust
       pmass(i) = masspdust
       uu(i) = 0. 
    else
       dens(i) = denszero
       pmass(i) = massp
       uu(i) = uuzero
    endif
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
 vcrap = sqrt(spsoundi**2 + valfven2i)

!----------------------------
! set the wave speed to use
!----------------------------
                    
 if (iwave.EQ.1) then
    vwave = vslow
    if (imhd.le.0) write(iprint,*) 'Error: can''t have slow wave if no mhd'
 elseif (iwave.EQ.2) then
    vwave = vfast
 else
    vwave = spsoundi
 endif
 if (vwave.le.0.) then
    write(iprint,*) 'Error in setup: vwave = ',vwave
    stop
 endif
  
 dxmax = xmax(1) - xmin(1)
 denom = dxmax - ampl/wk*(COS(wk*dxmax)-1.0)
 write (iprint,*) 'Wave number = ',wk
 write (iprint,*) 'Amplitude = ',ampl
!
!--now perturb the particles appropriately
!

 tausmode          = 0.1
 Bigkx             = 30.
 Bigkz             = 30.
 dust_to_gas_ratio = 3.d0
! Rux   = -0.1691398
 Rux   =  0.0000224*1000.
! Iux   =  0.0361553
 Iux   =  0.0000212*1000.
 Ruy   =  0.1336704
 Iuy   =  0.0591695
! Ruz   =  0.1691389
 Ruz   = -0.0000224*1000.
! Iuz   = -0.0361555
 Iuz   =  0.0000212*1000.
 Rrhog =  0.0000224
 Irhog =  0.0!000212
 Rwx   = -0.1398623
 Iwx   =  0.0372951
 Rwy   =  0.1305628
 Iwy   =  0.0640574
 Rwz   =  0.1639549
 Iwz   = -0.0233277

 ffac = 1000.0

 Rrhog = Rrhog*denszero*ampl
 Irhog = Irhog*denszero*ampl
 Rrhod = ampl*denszero
 Irhod = 0. 
!--convert the data YJ in code units (v(code) = eta v(YJ)) 
 Rux   = ampl*Rux*eta
 Iux   = ampl*Iux*eta
 Ruy   = ampl*Ruy*eta
 Iuy   = ampl*Iuy*eta
 Ruz   = ampl*Ruz*eta
 Iuz   = ampl*Iuz*eta
 Rwx   = ampl*Rwx*eta
 Iwx   = ampl*Iwx*eta
 Rwy   = ampl*Rwy*eta
 Iwy   = ampl*Iwy*eta
 Rwz   = ampl*Rwz*eta
 Iwz   = ampl*Iwz*eta

 do i=1,ntotal
       xi      = x(1,i)
       zi      = x(2,i)
       cokx0   = cos(kx*xi)
       sikx0   = sin(kx*xi)
       cokz0   = cos(kz*zi)
!-------------------------------------------------------------
! start with the gas particles
!-------------------------------------------------------------
    if (itype(i).eq.itypegas) then
       if (iperturb.eqv. .true.) then
   !setup a density perturbation (Rrhog*cokx - Irhog*sikx)*cokz       
   !--setup the perturbation on x
          call deltarhox(xi,cokz0,kx,Rrhog,Irhog)
          x(1,i) = xi
   !--setup the perturbation on z
          fx0    = Rrhog*cokx0 - Irhog*sikx0   
          call deltarhoz(zi,fx0,kz)
          x(2,i) = zi
       endif
       !--setup the kinematic quantities (v = vk + vNSH86 + vpert)
       cokx = cos(kx*xi)
       sikx = sin(kx*xi)
       cokz = cos(kz*zi)
       sikz = sin(kz*zi)
       vel(:,i) = 0.

       if (iperturb.eqv. .true.) then
!          vel(1,i) = vel(1,i) + (Rux*cokx - Iux*sikx)*cokz
!          vel(2,i) = vel(2,i) - (Ruz*sikx + Iuz*cokx)*sikz
       endif
       
       dens(i)  = denszero
       pmass(i) = massp
       uu(i)    = uuzero
    endif


!    
!--perturb density if not using summation
!
!    Bfield(2,i) = Bfield(2,i) + vwave*Bfield(2,i)*vamplx/term*SIN(wk*dxi)
!    Bfield(3,i) = Bfield(3,i) + vwave*Bfield(3,i)*vamplx/term*SIN(wk*dxi)
!    dens(i) = dens(i)*(1.+ampl*SIN(wk*dxi))

 enddo

 write(iprint,*) ' Wave set: Amplitude = ',ampl,' Wavelength = ',xlambda,' k = ',wk
 write(iprint,*) ' sound speed  = ',spsoundi
 write(iprint,*) ' alfven speed = ',sqrt(valfven2i)
 write(iprint,*) ' fast speed   = ',vfast
 write(iprint,*) ' slow speed   = ',vslow
 write(iprint,*) ' wave speed   = ',vwave,iwave
 
 return
end subroutine setup


!--------------------------------------------------------------------------------
!-setup the density perturbation on x 
!--------------------------------------------------------------------------------
 subroutine deltarhox(xi,cokz0,kx,Rrho,Irho)
  implicit none
  
  real, intent(inout) :: xi
  real, intent(in)    :: cokz0,kx,Rrho,Irho
  
  integer :: its
  real    :: xi0,xprev,func,fderiv
  integer, parameter :: itsmax = 20
  real, parameter    :: tol = 1.e-12
 
   xi0 = xi
   its = 0
   xprev = 1.0E6
   do while ((abs(xi-xprev).gt.tol).and.(its.lt.itsmax))
      xprev  = xi
      func   = xi0 - xi -cokz0/kx*(Rrho*sin(kx*xi) + Irho*(cos(kx*xi) - 1.))
      fderiv = -1. - cokz0*(Rrho*cos(kx*xi) - Irho*sin(kx*xi))
      xi     = xi - func/fderiv        ! Newton-Raphson iteration
      its    = its + 1
   enddo
   if (its.GE.itsmax) then
      print*,'Error: SI on x - too many iterations'
      stop
   endif 
 end subroutine deltarhox


!--------------------------------------------------------------------------------
!-setup the density perturbation on z 
!--------------------------------------------------------------------------------
 subroutine deltarhoz(zi,fx0,kz)
  implicit none

  real, intent(inout) :: zi
  real, intent(in)    :: fx0,kz

  integer :: its
  real    :: zi0,zprev,func,fderiv
  integer, parameter :: itsmax = 20
  real, parameter    :: tol = 1.e-12
  
  zi0 = zi
  its = 0
  zprev = 1.0E6
  do while ((abs(zi-zprev).gt.tol).and.(its.lt.itsmax))
     zprev  = zi
     func   = zi0 - zi - fx0/kz*sin(kz*zi)
     fderiv = -1. - fx0*cos(kz*zi)
     zi     = zi - func/fderiv        ! Newton-Raphson iteration
     its    = its + 1
  enddo
  if (its.GE.itsmax) then
     print*,'Error: SI on z(gas) - too many iterations'
     stop
  endif
  
 end subroutine deltarhoz
