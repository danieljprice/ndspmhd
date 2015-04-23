!----------------------------------------------------------------
!  Set up for the 2D epicyclic oscillations in a 
! "(r-z)" shearing box. 
! All the quantities are zero, except the radial and the azimuthal velocites
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
 use streaming
 
 use uniform_distributions
 use mem_allocation, only:alloc
!
!--define local variables
!            
 implicit none
 integer :: i,ngas,ndust,jtype
 real    :: massp,masspdust,totmass,totmassdust
 real    :: denszero,uuzero,cs0,polyk0
 real    :: przero,asize
 real    :: denszerodust,dust_to_gas_ratio
 real    :: L,kz
 real    :: xi,zi

 integer, parameter          :: ntypes   = 2
 real, parameter             :: ampl     = 1.d-7
 real, parameter             :: scaling  = 0.001 !--0.025
 
!-------------------------------------------------------------
! set the external forces and the boundary conditions
!-------------------------------------------------------------

!--allow for tracing flow
 if (trace) write(iprint,*) ' entering subroutine setup (SI)'
 write(iprint,*) '2D streaming instability'
 if (ndim.ne.2) stop 'error: this is a 2D problem and ndim.ne.2'
 if (ndimV.ne.3) stop 'error: we need ndimV=3 for this problem'

!--set position of disc patch for coriolis & centrifugal forces
 iexternal_force = 11 !--5

!-setup the initial quantities
 dust_to_gas_ratio = 1.

!--set boundaries
 ibound(:) = 3     ! periodic in y,z
 ibound(1) = 5     ! shearing box in x (==R)
 nbpts = 0         ! use ghosts not fixed
 asize = 1.0
 xmin(:) = -asize*scaling
 xmax(:) = asize*scaling
 psep    = psep*scaling

 if (abs(asize-1.).lt.tiny(0.)) print*,'WARNING: if asize .ne. 1, change cs'

!-------------------------------------------------------------
! set the particles and their masses
!-------------------------------------------------------------

!--setup uniform density grid of gas and dust particles
 ngas  = 0
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
 ntotal = ngas + ndust

 denszero     = 1.0
 denszerodust = denszero*dust_to_gas_ratio
 Kdrag = 0.
 print*,'enforce Kdrag =',Kdrag

!--determine particle mass
 L           = (xmax(1) - xmin(1))
 totmass     =     denszero*L*L
 totmassdust = denszerodust*L*L
 massp       =     totmass/FLOAT(ngas)
 masspdust   = totmassdust/FLOAT(ndust)
 !--kx in code units (one wavelenght per box)
 kz          = 2.*pi/L

!-------------------------------------------------------------
! setup the gas quantities
!-------------------------------------------------------------
 cs0    = 0.1
 przero = denszero*cs0*cs0/gamma
 uuzero = przero/(denszero*(gamma-1.))
 polyk0 = przero/denszero**gamma

 write(iprint,*) ' gamma = ', gamma
 write(iprint,*) ' cs0 = ',cs0,' cs2 = ',cs0**2,' u = ',uuzero
 write(iprint,*) ' Omega_c = ',Omega0, ' R_c = ',Rcentre
 write(iprint,*) ' polyk = ',polyk,' setting to = ',polyk0
 polyk = polyk0
 write(iprint,*) ' orbital time = ',2.*pi/Omega0
 write(iprint,*) ' c / (R omega) = ',cs0/(Rcentre*Omega0)

 do i=1,ntotal
       xi      = x(1,i)
       zi      = x(2,i)
!-------------------------------------------------------------
! start with the gas particles
!-------------------------------------------------------------
    if (itype(i).eq.itypegas) then
    
       !--setup the kinematic quantities
       vel(:,i) = 0.
       vel(3,i) = -domegadr*Omega0*xi 

       vel(1,i) = vel(1,i) + ampl*cos(kz*zi)
       vel(3,i) = vel(3,i) + ampl*cos(kz*zi)        
      
       dens(i)  = denszero
       pmass(i) = massp        
       uu(i)    = uuzero
          
!-------------------------------------------------------------
! now setup the dust particles
!-------------------------------------------------------------
    elseif(itype(i).eq.itypedust) then

    !--setup the kinematic quantities
       vel(:,i) = 0.
       vel(3,i) = -domegadr*Omega0*xi
        
       vel(1,i) = vel(1,i) + ampl*cos(kz*zi)
       vel(3,i) = vel(3,i) + ampl*cos(kz*zi)        

       dens(i)  = denszerodust
       pmass(i) = masspdust     
       uu(i)    = 0.

    endif
 enddo

!--allow for tracing flow
 if (trace) write(iprint,*) '  exiting subroutine setup'

 return
end subroutine setup

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
