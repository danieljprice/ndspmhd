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
 integer :: ntypes
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
 real :: dust_to_gas_ratio
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
    ampl = 0.0001
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
 xmin(:) = 0.   ! set position of boundaries
 xmax(1) = 1.0 
 if (ndim.GE.2) then
    xmax(2:ndim) = 11.*psep ! would need to adjust this depending on grid setup
 endif
!
!--initially set up a uniform density grid (also determines npart)
!  (the call to set_uniform_cartesian means this works in 1,2 and 3D)
!
 if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
    ntypes = 1
 else
    ntypes = 2
 endif

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
 massp = 1.0/FLOAT(ngas)        ! average particle mass
 
 masspdust = 0.
 dust_to_gas_ratio = 1.
 if (ndust.gt.0) masspdust = dust_to_gas_ratio*1.0/FLOAT(ndust) ! average particle mass
 denszerodust = dust_to_gas_ratio*denszero
 if (ntypes.gt.1) print*,' ngas = ',ngas,' ndust = ',ndust
!
!--allocate memory here
!
! call alloc(npart)
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
    if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
       dustfrac(i) = dust_to_gas_ratio/(1. + dust_to_gas_ratio)
       deltav(:,i) = 0.
       pmass(i) = pmass(i)/(1. - dustfrac(i))
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
 do i=1,npart
    
    dxi = x(1,i)-xmin(1)
    dxprev = dxmax*2.
    xmassfrac = dxi/dxmax        ! current mass fraction(for uniform density)
                                ! determines where particle should be
!
!--Use rootfinder on the integrated density perturbation
!  to find the new position of the particle
!    
    its = 0
         
    do while ((abs(dxi-dxprev).gt.tol).and.(its.lt.itsmax))
       dxprev = dxi
       func = xmassfrac*denom - (dxi - ampl/wk*(COS(wk*dxi)-1.0))
       fderiv = -1.0 - ampl*SIN(wk*dxi)
       dxi = dxi - func/fderiv        ! Newton-Raphson iteration
       its = its + 1 
!      PRINT*,'iteration',its,'dxi =',dxi 
    enddo
                           
    if (its.GE.itsmax) then
       write(iprint,*) 'Error: soundwave - too many iterations'
       call quit 
    endif
          
!    if (idebug(1:5).EQ.'sound') then
!       write(*,99002) i,its,dxi-dxprev,x(1,i),xmin(1)+dxi
!99002  FORMAT('Particle',i5,' converged in ',i4,        &
!       ' iterations, error in x =',1(1pe10.2,1x),/,        &
!       'previous x = ',1(0pf8.5,1x),                        &
!       'moved to x = ',1(0pf8.5,1x)) 
!    endif
    x(1,i) = xmin(1)+dxi

!
!--multiply by appropriate wave speed
!
    term = vwave**2 - Bfield(1,i)**2/dens(i)
    vamplx = vwave*ampl
    vamply = -vamplx*Bfield(1,i)*Bfield(2,i)/(dens(i)*term)
    vamplz = -vamplx*Bfield(1,i)*Bfield(3,i)/(dens(i)*term)
    
    vel(1,i) = vamplx*SIN(wk*dxi)
!    vel(2,i) = vamply*SIN(wk*dxi)
!    vel(3,i) = vamplz*SIN(wk*dxi)
!
!--perturb internal energy if not using a polytropic equation of state 
!  (do this before density is perturbed)
!
    if (itype(i).eq.itypegas) then
       uu(i) = uu(i) + przero/dens(i)*ampl*SIN(wk*dxi)        ! if not polytropic
    endif
!    
!--perturb density if not using summation
!
!    Bfield(2,i) = Bfield(2,i) + vwave*Bfield(2,i)*vamplx/term*SIN(wk*dxi)
!    Bfield(3,i) = Bfield(3,i) + vwave*Bfield(3,i)*vamplx/term*SIN(wk*dxi)
    dens(i) = dens(i)*(1.+ampl*SIN(wk*dxi))

 enddo

 write(iprint,*) ' Wave set: Amplitude = ',ampl,' Wavelength = ',xlambda,' k = ',wk
 write(iprint,*) ' sound speed  = ',spsoundi
 write(iprint,*) ' alfven speed = ',sqrt(valfven2i)
 write(iprint,*) ' fast speed   = ',vfast
 write(iprint,*) ' slow speed   = ',vslow
 write(iprint,*) ' wave speed   = ',vwave,iwave
 
 return
end subroutine setup

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
