!----------------------------------------------------------------
!  Set up for the 2D cartesian Streaming Imstability MRI
!  test problem described in Youdin and Johansen (2007) +
!  the additional problems of Bai and Stone (2010)
!  It is possible to add magnetic field
!
! Be careful to the dimensionless quantities:
! l(phys)   = l(code)R           = l(YJ)eta R    
! t(phys)   = t(code)/Omega      = t(YJ)/Omega   
! k(pkys)   = k(code)/R          = k(YJ)/(eta R) 
! v(phys)   = v(code) vk         = v(YJ)(eta vk) 
! rho(phys) = rho(code) Msol/R^3 = rho(code) Msol/(eta R)^3
!----------------------------------------------------------------
!--la soluce NSH86 est calculkee  une precision plus faible que la perturbation.
!--tester la perturbation sans drag (linAwg): est ce qu'on a un increase de la slution ??


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
 real :: massp,masspdust,totmass,totmassdust
 real :: denszero,uuzero,cs0,polyk0
 real :: przero,Bzeroz,asize,betamhd,wavekmin,valfvenz
 real :: denszerodust,dust_to_gas_ratio
 real :: L,kx,kz,Bigkx,Bigkz,tausmode
 real :: Rux,Iux,Ruy,Iuy,Ruz,Iuz,Rwx,Iwx,Rwy,Iwy,Rwz,Iwz
 real :: Rrhog,Irhog,Rrhod,Irhod
 real :: onepluseps,denomNSH,vgxNSH,vgyNSH,vdxNSH,vdyNSH
 real :: xi,zi,cokx0,sikx0,cokz0,sikz0
 real :: hL
! real :: x1,x2,x3,x4,z1,z2,z3,z4,pert,xstart,zstart
 real :: kxfac, kzfac
 logical :: iperturb

 integer, parameter          :: ntypes = 2
 integer, parameter          :: itsmax = 20
 real, parameter             :: ampl   = 1.d-3
 real, parameter             :: tol    = 1.e-8
! real, parameter             :: corr   = 0.99975747 !--correction to handle the real SPH density
 real, parameter             :: corr   = 0.99994416703041011
 real, parameter             :: vcorrg   = 0.99930182
 real, parameter             :: vcorrd   = 0.99930212
 character(len=10),parameter :: lin    = 'appBsimple'
 
!-------------------------------------------------------------
! set the external forces and the boundary conditions
!-------------------------------------------------------------

!--allow for tracing flow
 if (trace) write(iprint,*) ' entering subroutine setup (SI)'
 write(iprint,*) '2D streaming instability'
 if (ndim.ne.2) stop 'error: this is a 2D problem and ndim.ne.2'
 if (ndimV.ne.3) stop 'error: we need ndimV=3 for this problem'

!--set position of disc patch for coriolis & centrifugal forces
 iexternal_force = 11 ! additional correction due to the background pressure gradient
! iexternal_force = 5

!-get the drag quantities and the SI perturbations
 call getperturb(lin,Bigkx,Bigkz,tausmode,dust_to_gas_ratio,Rux,Iux, &
                 Ruy,Iuy,Ruz,Iuz,Rrhog,Irhog,Rwx,Iwx,Rwy,Iwy,Rwz,Iwz,iperturb, &
                 cs0,kxfac,kzfac,Rrhod)

!--set boundaries
 ibound(:) = 3     ! periodic in y,z
 ibound(1) = 5     ! shearing box in x (==R)
 nbpts = 0         ! use ghosts not fixed
 asize = 1.0!pi       ! box size (corresponds to Bai/Stone 2010)
! xmin(:) = -asize*eta/Bigkx  ! set position of boundaries
! xmax(:) = asize*eta/Bigkx
! psep    = psep*eta/Bigkx
 xmin(:) = -asize!_-*0.05
 xmax(:) = asize!--*0.05
 psep    = psep!--*0.05

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
 if (trim(lin).eq.'nodrag' .or. trim(lin).eq.'linAwg' &
      .or. trim(lin).eq.'nodragpert') then
    Kdrag = 0.
 else
    Kdrag = (denszerodust*corr)*Omega0/tausmode
 endif
 print*,'WARNING: Kdrag has been modified to correspond to the SI test. Kdrag=',Kdrag
!--determine particle mass
 L           = (xmax(1) - xmin(1))!*eta/Bigkx
 hL          = 0.5*L
 totmass     =     denszero*L*L
 totmassdust = denszerodust*L*L
 massp       =     totmass/FLOAT(ngas)
 masspdust   = totmassdust/FLOAT(ndust)
 !--kx,kz in code units (one wavelenght per box)
 kx          = 0. !--kxfac*2.*pi/L
 kz          = kzfac*2.*pi/L

!-------------------------------------------------------------
! setup the gas quantities
!-------------------------------------------------------------
 cs0    = 0.1   !--cfJY07
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

!-------------------------------------------------------------
! setup the parameters for the SI test
!-------------------------------------------------------------
 onepluseps   = 0.
 denomNSH     = 0.
 select case(trim(lin))
 case('nodrag','nodragpert')
    vgxNSH = 0.        
    vgyNSH = -eta      
    vdxNSH = 0.        
    vdyNSH = 0.     
 case('NSH','NSHpert','linA','linAwg','linB') !--enforce the SPH solution for the background
 !quantities calculated for eta=0.01
    vgxNSH = 1.873813705237E-04 !3.75094634d-4*0.5*vcorrg ! 1.87381344594547123e-04
    vgyNSH = -2.50469218d-3*0.5 ! -1.25235010733846901d-3
    vdxNSH = -6.246053793239E-05 !-1.25031547d-4*0.5*vcorrd !-6.24605585029068132e-05
    vdyNSH = -2.49843602d-3*0.5 ! -1.24921927251258388

 case('appB','appBsimple')
    vgxNSH = 0.
    vgyNSH = 0.
    vdxNSH = 0.
    vdyNSH = 0.

 case default
    onepluseps   = 1. + dust_to_gas_ratio
    denomNSH     = onepluseps**2 + tausmode**2
    vgxNSH       = 2.*dust_to_gas_ratio*tausmode/denomNSH*eta
    vgyNSH       = -(1. + dust_to_gas_ratio*tausmode**2/denomNSH)*eta/onepluseps
    vdxNSH       = -2.*tausmode/denomNSH*eta
    vdyNSH       = -(1. - tausmode**2/denomNSH)*eta/onepluseps
 end select
 write(iprint,*) ' vgxNSH = ',vgxNSH
 write(iprint,*) ' vgyNSH = ',vgyNSH
 write(iprint,*) ' vdxNSH = ',vdxNSH
 write(iprint,*) ' vdyNSH = ',vdyNSH
 
!-------------------------------------------------------------
! setup the magnetic field (beta = pr/(0.5*B^2))
!-------------------------------------------------------------

 if (imhd.ne.0) then
    betamhd = 4000.
    Bzeroz = sqrt(2.*przero/betamhd)
    write(iprint,*) ' beta(MHD) = ',betamhd,' Bz = ',Bzeroz,' P = ',przero

    valfvenz = sqrt(Bzeroz**2/denszero)
    wavekmin = sqrt(15.)/4.*Omega0/valfvenz
    write(iprint,*) 'most unstable wavelength = ',2.*pi/wavekmin
 endif

!-------------------------------------------------------------
! assign particles properties
!-------------------------------------------------------------
 Rrhog = Rrhog*denszero*ampl
 Irhog = Irhog*denszero*ampl
 Rrhod = Rrhod*ampl*denszero
 
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
       sikz0   = sin(kz*zi)
!-------------------------------------------------------------
! start with the gas particles
!-------------------------------------------------------------
    if (itype(i).eq.itypegas) then
       if (iperturb.eqv. .true. .AND. trim(lin).ne.'appB' &
           .and. trim(lin).ne.'appBsimple') then
       !setup a density perturbation Rrhog*cos(kx*x + kz*z) - Irhog*sin(kx*x + kz*z) 
!          xstart  = xi 
!          zstart  = zi
!          !CAREFUL: keep 0.5 factor only for perturbations in 2D !
!          !--term 1 : Rrhog*cos(kx*x)*cos(kz*z)
!           pert = 0.5*Rrhog*cokz0/denszero
!           call deltarho_cos(xi,pert,kx,-hL,hL) !--pert in x
!           x1   = xi - xstart
!           xi   = xstart          
!           pert = 0.5*Rrhog*cokx0/denszero  
!           call deltarho_cos(zi,pert,kz,-hL,hL) !--pert in z
!           z1   = zi - zstart
!           zi   = zstart
!          
!          !--term 2 : -Rrhog*sin(kx*x)*sin(kz*z)
!            pert = -0.5*Rrhog*sikz0/denszero
!            call deltarho_sin(xi,pert,kx,-hL,hL) !--pert in x
!            x2   = xi - xstart
!            xi   = xstart          
!            pert = -0.5*Rrhog*sikx0/denszero 
!            call deltarho_sin(zi,pert,kz,-hL,hL) !--pert in z
!            z2   = zi - zstart
!            zi   = zstart
!                    
!          !--term 3 : -Irhog*sin(kx*x)*cos(kz*z)
!            call deltarho_sin(xi,pert,kx,-hL,hL) !--pert in x
!            x3   = xi - xstart
!            xi   = xstart          
!            pert = -0.5*Irhog*sikx0/denszero  
!            call deltarho_cos(zi,pert,kz,-hL,hL) !--pert in z
!            z3   = zi - zstart
!            zi   = zstart         
!          
!           !--term 4 : -Irhog*cos(kx*x)*sin(kz*z)
!            pert = -0.5*Irhog*sikz0/denszero
!            call deltarho_cos(xi,pert,kx,-hL,hL) !--pert in x
!            x4   = xi - xstart
!            xi   = xstart          
!            pert = -0.5*Irhog*cokx0/denszero  
!            call deltarho_sin(zi,pert,kz,-hL,hL) !--pert in z
!            z4   = zi - zstart
!            zi   = zstart          
          
!          !--sum all the contributions
!          x(1,i) = xstart + x1 + x2 + x3 + x4
!          x(2,i) = zstart + z1 + z2 + z3 + z4
!          xi     = x(1,i)
!          zi     = x(2,i)
  
!----for adding a pert in kz only--------------------------- 
         !pert = 1.d-7
         !call deltarho_cos(zi,pert,kz,-hL,hL) !--pert in z
         !x(2,i) = zi    
!-----------------------------------------------------------------          
       endif
       !--setup the kinematic quantities (v = vk + vNSH86 + vpert)
       vel(:,i) = 0.
       vel(1,i) = vgxNSH
       vel(3,i) = -domegadr*Omega0*xi + vgyNSH
                         
!----for adding a pert in kz only-------------------------------- 
         vel(1,i) = vel(1,i) + 1.d-7*cos(kz*zi)
         vel(2,i) = vel(2,i) + 0.*cos(kz*zi)
         vel(3,i) = vel(3,i) + 1.d-7*cos(kz*zi)        
!----------------------------------------------------------------- 
       
!       if (iperturb.eqv. .true.) then 
!          vel(1,i) = vel(1,i) + Rux*cos(kx*xi + kz*zi) - Iux*sin(kx*xi + kz*zi)
!          vel(2,i) = vel(2,i) + Ruz*cos(kx*xi + kz*zi) - Iuz*sin(kx*xi + kz*zi)
!          vel(3,i) = vel(3,i) + Ruy*cos(kx*xi + kz*zi) - Iuy*sin(kx*xi + kz*zi)          
!       endif
       
       dens(i)  = denszero
       pmass(i) = massp
       uu(i)    = uuzero

    !--if present, setup the magnetic field
       if (imhd.gt.0) then
          Bfield(:,i) = 0.
          Bfield(2,i) = Bzeroz
          Bfield(3,i) = 0.
       elseif (imhd.lt.0) then
          Bevol(:,i) = 0.
          Bevol(3,i) = 0.
       endif
          
!-------------------------------------------------------------
! now setup the dust particles
!-------------------------------------------------------------
    elseif(itype(i).eq.itypedust) then
       if (iperturb .eqv. .true. .AND. trim(lin).ne.'appB' &
           .and. trim(lin).ne.'appBsimple') then
       !setup a density perturbation Rrhod*cos(kx*x + kz*z) - Irhod*sin(kx*x + kz*z) 
!          xstart  = xi 
!          zstart  = zi
!          !--term 1 : Rrhod*cos(kx*x)*cos(kz*z)
!          pert = 0.5*Rrhod*cokz0/denszerodust
!          call deltarho_cos(xi,pert,kx,-hL,hL) !--pert in x
!          x1   = xi - xstart
!          xi   = xstart          
!          pert = 0.5*Rrhod*cokx0/denszerodust  
!          call deltarho_cos(zi,pert,kz,-hL,hL) !--pert in z
!          z1   = zi - zstart
!          zi   = zstart
!          
!         !--term 2 : -Rrhod*sin(kx*x)*sin(kz*z)
!          pert = -0.5*Rrhod*sikz0/denszerodust
!          call deltarho_sin(xi,pert,kx,-hL,hL) !--pert in x
!          x2   = xi - xstart
!          xi   = xstart          
!          pert = -0.5*Rrhod*sikx0/denszerodust  
!          call deltarho_sin(zi,pert,kz,-hL,hL) !--pert in z
!          z2   = zi - zstart
!          zi   = zstart
!
!        !--term 3 : -Irhod*sin(kx*x)*cos(kz*z)
!          pert = -0.5*Irhod*cokz0/denszerodust
!          call deltarho_sin(xi,pert,kx,-hL,hL) !--pert in x
!          x3   = xi - xstart
!          xi   = xstart          
!          pert = -0.5*Irhod*sikx0/denszerodust 
!          call deltarho_cos(zi,pert,kz,-hL,hL) !--pert in z
!          z3   = zi - zstart
!          zi   = zstart         
!
!        !--term 4 : -Irhod*cos(kx*x)*sin(kz*z)
!          pert = -0.5*Irhod*sikz0/denszerodust
!          call deltarho_cos(xi,pert,kx,-hL,hL) !--pert in x
!          x4   = xi - xstart
!          xi   = xstart          
!          pert = -0.5*Irhod*cokx0/denszerodust 
!          call deltarho_sin(zi,pert,kz,-hL,hL) !--pert in z
!          z4   = zi - zstart
!          zi   = zstart          
          
          !--sum all the contributions
!          x(1,i) = xstart + x1 + x2 + x3 + x4
!          x(2,i) = zstart + z1 + z2 + z3 + z4
!          xi     = x(1,i)
!          zi     = x(2,i)


!----for adding a pert in kz only--------------------------- 
        ! pert = 1.d-7/denszerodust
        ! call deltarho_cos(zi,pert,kz,-hL,hL) !--pert in z
        ! x(2,i) = zi         
!-----------------------------------------------------------------        

       endif
    !--setup the kinematic quantities (v = vk + vNSH86 + vpert)
       vel(:,i) = 0.
       vel(1,i) = vdxNSH
       vel(3,i) = -domegadr*Omega0*xi + vdyNSH 
       
!----for adding a pert in kz only-------------------------------- 
!         vel(1,i) = vel(1,i) + 1.d-7*cos(kz*zi)
!         vel(2,i) = vel(2,i) + 0.*cos(kz*zi)
!         vel(3,i) = vel(3,i) + 1.d-7*cos(kz*zi)        
!-----------------------------------------------------------------      
       
       
!       if (iperturb.eqv. .true.) then
!          vel(1,i) = vel(1,i) + Rwx*cos(kx*xi + kz*zi) - Iwx*sin(kx*xi + kz*zi)
!          vel(2,i) = vel(2,i) + Rwz*cos(kx*xi + kz*zi) - Iwz*sin(kx*xi + kz*zi)
!          vel(3,i) = vel(3,i) + Rwy*cos(kx*xi + kz*zi) - Iwy*sin(kx*xi + kz*zi)
!       endif
       dens(i)  = denszerodust
       pmass(i) = masspdust
       uu(i)    = 0.

    endif
 enddo

!--allow for tracing flow
 if (trace) write(iprint,*) '  exiting subroutine setup'

 return
end subroutine setup

!--------------------------------------------------------------------------------
!-initialise the different SI tests
!--------------------------------------------------------------------------------

subroutine getperturb(lin,Bigkx,Bigkz,tausmode,dust_to_gas_ratio,Rux,Iux, &
                      Ruy,Iuy,Ruz,Iuz,Rrhog,Irhog,Rwx,Iwx,Rwy,Iwy,Rwz,Iwz, &
                      iperturb,cs0,kxfac,kzfac,Rrhod)
 use setup_params, only:pi
 implicit none
 character(len=10), intent(in) :: lin
 real, intent(out)    :: Bigkx,Bigkz,tausmode,dust_to_gas_ratio, cs0
 real, intent(out)    :: Rrhod,kxfac,kzfac
 real, intent(out)    :: Rux,Iux,Ruy,Iuy,Ruz,Iuz,Rrhog,Irhog,Rwx,Iwx,Rwy,Iwy,Rwz,Iwz
 logical, intent(out) :: iperturb

 cs0    = 0.1   !--cfJY07
 kxfac = 1.0
 kzfac = 1.0
 Rrhod = 1.0

 tausmode          = 1.0
 Bigkx             = 0.
 Bigkz             = 0.
 dust_to_gas_ratio = 1.
 Rux   = 0.
 Iux   = 0.
 Ruy   = 0.
 Iuy   = 0.
 Ruz   = 0.
 Iuz   = 0.
 Rrhog = 0.
 Irhog = 0.
 Rwx   = 0.
 Iwx   = 0.
 Rwy   = 0.
 Iwy   = 0.
 Rwz   = 0.
 Iwz   = 0.
 iperturb = .true.

 select case(trim(lin))
  
 case('nodrag')
    iperturb = .false.
    Bigkx             = 30.
    Bigkz             = 30.
    dust_to_gas_ratio = 3.d0
    
 case('nodragpert')
    Bigkx             = 30.
    Bigkz             = 30.
    dust_to_gas_ratio = 3.d0
    Rux   = -0.1691398
    Iux   =  0.0361553
    Ruy   =  0.1336704
    Iuy   =  0.0591695
    Ruz   =  0.1691389
    Iuz   = -0.0361555
    Rrhog =  0.0000224
    Irhog =  0.0000212
    Rwx   = -0.1398623
    Iwx   =  0.0372951
    Rwy   =  0.1305628
    Iwy   =  0.0640574
    Rwz   =  0.1639549
    Iwz   = -0.0233277 
  
 case('NSH')
    iperturb = .false.
    tausmode          = 0.1
    Bigkx             = 30.
    Bigkz             = 30.
    dust_to_gas_ratio = 3.d0
    
 case('NSHpert')
    iperturb = .true.
    tausmode          = 0.1
    Bigkx             = 30.
    Bigkz             = 30.
    dust_to_gas_ratio = 3.d0 

 case('appB')
    tausmode = 1.0
    cs0 = 1.0
    kxfac = -1.0
    kzfac = 1.0

    Rux = 1.0d-3
!    Ruy = 1.0d-3

    Rrhog = 0.0
    IRhog = 0.0
    Rrhod = 0.0
    Bigkx = 1.0
    Bigkz = 1.0
    
 case('appBsimple')
    iperturb = .true.
    tausmode          = 3.d0
    Bigkx             = 30.
    Bigkz             = 30.
    dust_to_gas_ratio = 3.d0
     
 case('linAwg')
    Bigkx             = 30.
    Bigkz             = 30.
    dust_to_gas_ratio = 3.d0
    Rux   = -0.1691398
    Iux   =  0.0361553
    Ruy   =  0.1336704
    Iuy   =  0.0591695
    Ruz   =  0.1691389
    Iuz   = -0.0361555
    Rrhog =  0.0000224
    Irhog =  0.0000212
    Rwx   = -0.1398623
    Iwx   =  0.0372951
    Rwy   =  0.1305628
    Iwy   =  0.0640574
    Rwz   =  0.1639549
    Iwz   = -0.0233277

 case('linA')
    tausmode          = 0.1
    Bigkx             = 30.
    Bigkz             = 30.
    dust_to_gas_ratio = 3.d0
    Rux   = -0.1691398
    Iux   =  0.0361553
    Ruy   =  0.1336704
    Iuy   =  0.0591695
    Ruz   =  0.1691389
    Iuz   = -0.0361555
    Rrhog =  0.0000224
    Irhog =  0.0000212
    Rwx   = -0.1398623
    Iwx   =  0.0372951
    Rwy   =  0.1305628
    Iwy   =  0.0640574
    Rwz   =  0.1639549
    Iwz   = -0.0233277


 case('linB')
    tausmode          = 0.1
    Bigkx             = 6.
    Bigkz             = 6.
    dust_to_gas_ratio = 2.d-1
    Rux   = -0.0174121
    Iux   = -0.2770347
    Ruy   =  0.2767976
    Iuy   = -0.0187568
    Ruz   =  0.0174130
    Iuz   =  0.2770423
    Rrhog = -0.0000067
    Irhog = -0.0000691
    Rwx   =  0.0462916
    Iwx   = -0.2743072
    Rwy   =  0.2739304
    Iwy   =  0.0039293
    Rwz   =  0.0083263
    Iwz   =  0.2768866

 case('linC')
    tausmode          = 0.01
    Bigkx             = 1.5d3
    Bigkz             = 1.5d3
    dust_to_gas_ratio = 2.d0
    Rux   = -0.1598751
    Iux   =  0.0079669
    Ruy   =  0.1164423
    Iuy   =  0.0122377
    Ruz   =  0.1598751
    Iuz   = -0.0079669
    Rrhog =  8.684872d-8
    Irhog =  5.350037d-7
    Rwx   = -0.1567174
    Iwx   =  0.0028837
    Rwy   =  0.1159782
    Iwy   =  0.0161145
    Rwz   =  0.1590095
    Iwz   = -0.0024850

 case('linD')
    tausmode          = 0.001
    Bigkx             = 2.d3
    Bigkz             = 2.d3
    dust_to_gas_ratio = 2.
    Rux   = -0.1719650
    Iux   =  0.0740712
    Ruy   =  0.1918893
    Iuy   =  0.0786519
    Ruz   =  0.1719650
    Iuz   = -0.0740712
    Rrhog =  2.954631d-7
    Irhog =  1.141385d-7
    Rwx   = -0.1715840
    Iwx   =  0.0740738
    Rwy   =  0.1918542
    Iwy   =  0.0787371
    Rwz   =  0.1719675
    Iwz   = -0.0739160

 case default
    print*,'this setup does not exist for the SI',trim(lin)
    stop
 end select
end subroutine getperturb


!--------------------------------------------------------------------------------
!-setup to add a 'Ampl*cos(ki*xi)' type density perturbation
!--------------------------------------------------------------------------------
 subroutine deltarho_cos(xi,pert,ki,ximin,ximax)
  implicit none
  
  real, intent(inout) :: xi
  real, intent(in)    :: pert,ki,ximin,ximax
  
  integer :: its
  real    :: deltaxi0,xprev,func,fderiv,denom,deltamax
  integer, parameter :: itsmax = 20
  real, parameter    :: tol = 1.e-12

!--Sanity check
   deltamax   = ximax - ximin
   if (deltamax.le.0) then
      print*,'setup: deltamax<0 !'
      stop
   endif 

   if (abs(ki).gt.tiny(0.)) then 
      deltaxi0   = xi - ximin  
      denom      = deltamax + pert/ki*(sin(ki*ximax) - sin(ki*ximin) )  
      its   = 0
      xprev = 1.E6
      do while ((abs(xi-xprev).gt.tol).and.(its.lt.itsmax))
         xprev  = xi
         func   = deltaxi0*denom/deltamax - (xi-ximin) -pert/ki*(sin(ki*xi) - sin(ki*ximin))
         fderiv = -1. - pert*cos(ki*xi)
         xi     = xi - func/fderiv        ! Newton-Raphson iteration
         its    = its + 1
      enddo
      if (its.GE.itsmax) then
         print*,'Error: SI on x - too many iterations'
         stop
      endif 
   endif

 end subroutine deltarho_cos


!--------------------------------------------------------------------------------
!-setup to add a 'Ampl*sin(ki*xi)' type density perturbation 
!--------------------------------------------------------------------------------
 subroutine deltarho_sin(xi,pert,ki,ximin,ximax)
  implicit none
  
  real, intent(inout) :: xi
  real, intent(in)    :: pert,ki,ximin,ximax
  
  integer :: its
  real    :: deltaxi0,xprev,func,fderiv,denom,deltamax
  integer, parameter :: itsmax = 20
  real, parameter    :: tol = 1.e-12
  
!--Sanity check
   deltamax   = ximax - ximin
   if (deltamax.le.0) then
      print*,'setup: deltamax<0 !'
      stop
   endif  
 
   if (abs(ki).gt.tiny(0.)) then
      deltaxi0   = xi - ximin
      denom      = deltamax - pert/ki*(cos(ki*ximax)-cos(ki*ximin))
      its   = 0
      xprev = 1.E6
      do while ((abs(xi-xprev).gt.tol).and.(its.lt.itsmax))
         xprev  = xi
         func   = deltaxi0*denom/deltamax - (xi-ximin) +pert/ki*(cos(ki*xi) - cos(ki*ximin))
         fderiv = -1. - pert*sin(ki*xi)
         xi     = xi - func/fderiv        ! Newton-Raphson iteration
         its    = its + 1
      enddo
      if (its.GE.itsmax) then
         print*,'Error: SI on x - too many iterations'
         stop
      endif
   endif
  
 end subroutine deltarho_sin

