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
 integer :: i,ngas,npdust,jtype
 real :: massp,masspdust,totmass,totmassdust
 real :: denszero,uuzero,cs0,polyk0,denszeroc
 real :: przero,Bzeroz,asize,betamhd,wavekmin,valfvenz
 real :: denszerodust,dust_to_gas_ratio
 real :: L,kx,kz,Bigkx,Bigkz,tausmode
 real :: Rux,Iux,Ruy,Iuy,Ruz,Iuz,Rwx,Iwx,Rwy,Iwy,Rwz,Iwz
 real :: Rrhog,Irhog,Rrhod,Irhod
 real :: onepluseps,denomNSH,vgxNSH,vgyNSH,vdxNSH,vdyNSH
 real :: xi,zi,cokx0,sikx0,cokz0,sikz0
 real :: hL,eta_local
 real :: x1,x2,x3,x4,z1,z2,z3,z4,pert,xstart,zstart
 real :: kxfac, kzfac, combined_amp
 real :: Rrhogc, Irhogc, ddens, gdens
 real :: rhodonrho
 real :: vgx,vgy,vgz,vdx,vdy,vdz
 real :: cokx, sikx, cokz, sikz, tau_drag
 logical :: iperturb, vparity

 integer                     :: ntypes
 integer, parameter          :: itsmax = 20
 real, parameter             :: ampl   = 1.d-3
 real, parameter             :: tol    = 1.e-8
! real, parameter             :: corr   = 0.99975747 !--correction to handle the real SPH density
 real, parameter             :: corr   = 0.99994416703041011
 real, parameter             :: vcorrg   = 0.99930182
 real, parameter             :: vcorrd   = 0.99930212
 !character(len=20),parameter :: lin    = 'kxkzgasIdRvz'
 character(len=20),parameter :: lin    = 'linA'
 
!--one fluid approach
 if (idust.eq.1) then
    ntypes = 1
 else
    ntypes = 2
 endif

 vparity = .false.
! eta_local = eta
 print *, 'eta local ', eta_local, eta
 eta_local = 0.005
 if (trim(lin).eq.'shalloweta') eta_local = 0.001

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
 if (trim(lin).eq.'eigenNoShearNoPNoD') iexternal_force = 0

!-get the drag quantities and the SI perturbations
 call getperturb(lin,Bigkx,Bigkz,tausmode,dust_to_gas_ratio,Rux,Iux, &
                 Ruy,Iuy,Ruz,Iuz,Rrhog,Irhog,Rwx,Iwx,Rwy,Iwy,Rwz,Iwz,iperturb, &
                 cs0,kxfac,kzfac,Rrhod,Irhod,ampl,eta_local)

!--set boundaries
 ibound(:) = 3     ! periodic in y,z
 ibound(1) = 5     ! shearing box in x (==R)

 if (trim(lin).eq.'eigenNoShearNoPNoD') ibound(1) = 3

 nbpts = 0         ! use ghosts not fixed 
 asize = 0.005/30. !0.3 !0.3 !1./3. !0.005/30. !0.01 !1.0!pi       ! box size (corresponds to Bai/Stone 2010)
! xmin(:) = -asize*eta/Bigkx  ! set position of boundaries
! xmax(:) = asize*eta/Bigkx
! psep    = psep*eta/Bigkx
 xmin(:) = -asize!_-*0.05
 xmax(:) = asize!--*0.05
 psep    = psep*asize!--*0.05

 if (abs(asize-1.).lt.tiny(0.)) print*,'WARNING: if asize .ne. 1, change cs'

!-------------------------------------------------------------
! set the particles and their masses
!-------------------------------------------------------------

!--setup uniform density grid of gas and dust particles
 ngas  = 0
 npdust = 0
 do jtype=1,ntypes
  call set_uniform_cartesian(2,psep,xmin,xmax,adjustbound=.true.)
  if (jtype.eq.1) then
    ngas = npart
    itype(1:ngas) = itypegas
    print *, 'gas ', ngas
    print *, x(1,1:10)
  elseif (jtype.eq.2) then
    npdust = npart - ngas
    itype(ngas+1:ngas+npdust) = itypedust
    print *, 'dust ', ngas, npdust
    print *, x(1,ngas+1:ngas+10)
  endif
 enddo
 ntotal = ngas + npdust

 denszero     = 1.0
 denszerodust = denszero*dust_to_gas_ratio
 select case(trim(lin))
 case('nodrag','nodragpert','linAwg','kxepi','kxkzsimple','kxvd','kzvd', &
      'kxkzdensgas','kxkzdens','eigenNoPNoDrag','BBeigenNoPNoDrag', &
      'eigenNoShearNoPNoD')
    Kdrag = 0.
    tausmode = 0.

 case('dpkz','simpledrag','kxsimpledrag','dpkx','kxkzdrag','dpkxkz', &
      'kzvddrag','kxvddrag','kzvddragp','kxkzdensgasdrag','kxkzdensgasdragpr', &
      'kxkzdenswdrag','prestest')
    Kdrag = 1.
    tausmode = denszerodust/Kdrag

 case('strongdrag','kxstrongdrag','dpkxstrong','kxkzdragstrong', &
      'kxkzdenssdrag','kxkzdensgasI','kxkzgasIdRvz','linAsimple', &
      'kxkzdensgasIP','k10asym','damping','kxkzdustdens','kxkzvrgas', &
      'kxkzdustdensp')
    Kdrag = 10.
    tausmode = denszerodust/Kdrag

 case('dpkxkzstrong','kxkzdenssdragp','kxkzdustdensp30','kxkzdustdensp30_3')
    Kdrag = 30.
    tausmode = denszerodust/Kdrag

 case('smallcontrast','eigen1','eigen2','eigen4','eigennopres','eigensettle', &
      'mixedeigen2')
    Kdrag = 0.01
    tausmode = denszerodust/Kdrag

 case default
    Kdrag = (denszerodust*corr)*Omega0/tausmode
    tausmode = denszerodust/Kdrag
    tau_drag = denszerodust/Kdrag

 end select

 print *,'WARNING: Kdrag has been modified to correspond to the SI test. Kdrag=',Kdrag
 print *, 'Tau_drag = ', tau_drag
 print *, 'Pressure eta = ', eta, eta_local
!--determine particle mass
 L           = (xmax(1) - xmin(1))!*eta/Bigkx
 hL          = 0.5*L
 totmass     =     denszero*L*L
 totmassdust = denszerodust*L*L
 massp       =     totmass/FLOAT(ngas)
 if (idust.eq.1) then ! one fluid dust
    massp = massp*(1. + dust_to_gas_ratio)
 else
    masspdust   = totmassdust/FLOAT(npdust)
 endif

 print *, 'Particle mass ', massp, totmass, ngas, ntotal

 !--kx,kz in code units (one wavelenght per box)
 select case (trim(lin))
 case ('kxkzdrag','kxkzsimple','kxkzdragstrong','dpkxkz','dpkxkzstrong',&
      'kxkzdensgas','kxkzdensgasdrag','kxkzdensgasdragpr','kxkzdens', &
      'kxkzdenswdrag','kxkzdenssdrag','kxkzdenssdragp','kxkzdensgasI', &
      'kxkzgasIdRvz','linAsimple','linA','kxkzdensgasIP','damping', &
      'kxkzdustdens','kxkzvrgas','kxkzdustdensp','kxkzdustdensp30', &
      'kxkzdustdensp30_3','linAnew','shalloweta','smallcontrast','eigen1', &
      'eigen2','eigen4','eigennopres','contrasttest','mixedeigen', &
      'eigenNoPNoDrag','BBeigenNoPNoDrag','eigenNoShearNoPNoD','mixedeigen2')
    kx          = kxfac*2.*pi/L
 case default
    kx          = 0.0 !kxfac*2.*pi/L
 end select

 kz          = kzfac*2.*pi/L

 if (trim(lin).eq.'k10asym') kz = 0.

 select case (trim(lin))
 case ('kxepi','kxsimpledrag','kxstrongdrag','dpkx','dpkxstrong', &
      'kxvd','kxvddrag')
    kx          = kxfac*2.*pi/L
 case default
    
 end select

!-------------------------------------------------------------
! setup the gas quantities
!-------------------------------------------------------------
 cs0    = 0.1   !--cfJY07
 if (gamma > 1.) then
    przero = denszero*cs0*cs0/gamma
    uuzero = przero/(denszero*(gamma-1.))
    polyk0 = przero/denszero**gamma
 else
    przero = denszero*cs0*cs0
    uuzero = 1.5*cs0*cs0
    polyk0 = cs0*cs0 
 endif

 write(iprint,*) ' kx, kz = ', kx, kz
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
    vgyNSH = 10.*eta      
    vdxNSH = 0.        
    vdyNSH = 0.     
 case('NSH','NSHpert','linA','linAwg','linB','dpkz','dpkx','dpkxstrong', &
      'dpkxkz','dpkxkzstrong','kxkzdenssdragp','linAnew','mixedeigen') !--enforce the SPH solution for the background
 !quantities calculated for eta=0.01
!--Older values (from original asymptotic approach)
!    vgxNSH = 1.873813705237E-04
!    vgyNSH = -2.50469218d-3*0.5
!    vdxNSH = -6.246053793239E-05
!    vdyNSH = -2.49843602d-3*0.5

!    vgxNSH =  1.873813505437570E-04 ! 1.873813705237E-04 !3.75094634d-4*0.5*vcorrg ! 1.87381344594547123e-04
!    vgyNSH = -1.252342266461199E-03 !-2.50469218d-3*0.5  !-1.25235010733846901d-3
!    vdxNSH = -6.246055479643788E-05 !-6.246053793239E-05 !-1.25031547d-4*0.5*vcorrd !-6.24605585029068132e-05
!    vdyNSH = -1.249219266205729E-03 !-2.49843602d-3*0.5  !-1.24921927251258388

    vgxNSH =  1.87409978489950E-04 !1.873813505437570E-04 ! 1.873813705237E-04 !3.75094634d-4*0.5*vcorrg ! 1.87381344594547123e-04
    vgyNSH = -1.25233666827186E-03 !-1.252342266461199E-03 !-2.50469218d-3*0.5  !-1.25235010733846901d-3
    vdxNSH = -6.24699926597789E-05 !-6.246055479643788E-05 !-6.246053793239E-05 !-1.25031547d-4*0.5*vcorrd !-6.24605585029068132e-05
    vdyNSH = -1.24922109681910E-03 !-1.249219266205729E-03 !-2.49843602d-3*0.5  !-1.24921927251258388


 case('eigensettle','eigen2')
    vgxNSH =  7.936174510944069E-05
    vgyNSH = -4.979761466375418E-03
    vdxNSH = -3.968087255472082E-03
    vdyNSH = -1.011926681229109E-03

 case('appB','appBsimple','simpledrag','strongdrag','kxkzdrag','kxepi', &
      'kxsimpledrag','kxstrongdrag','kxkzsimple','kxkzdragstrong','kxvd', &
      'kzvd','kzvddrag','kxvddrag','kxkzdensgas','kxkzdensgasdrag','kxkzdens',&
      'kxkzdenswdrag','kxkzdenssdrag','kxkzdensgasI','kxkzgasIdRvz', &
      'kxkzdustdens','kxkzvrgas','eigenNoPNoDrag','BBeigenNoPNoDrag',&
      'eigenNoShearNoPNoD')
    vgxNSH = 0.
    vgyNSH = 0.
    vdxNSH = 0.
    vdyNSH = 0.

 case default
    onepluseps   = 1. + dust_to_gas_ratio
    denomNSH     = onepluseps**2 + tausmode**2
    vgxNSH       = 2.*dust_to_gas_ratio*tausmode/denomNSH*eta_local
    vgyNSH       = -(1. + dust_to_gas_ratio*tausmode**2/denomNSH)*eta_local/onepluseps
    vdxNSH       = -2.*tausmode/denomNSH*eta_local
    vdyNSH       = -(1. - tausmode**2/denomNSH)*eta_local/onepluseps

!    vgxNSH = 0.
!    vgyNSH = 0.
!    vdxNSH = 0.
!    vdyNSH = 0.

!    vgxNSH = 0. !2.49360859689E-04 
!    vgyNSH = 0. !-2.50623239911E-03
!    vdxNSH = 0. !-2.49360859689E-04
!    vdyNSH = 0. !-2.49376760089E-03

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
 combined_amp = ampl!*eta_local

 Rrhog = Rrhog*denszero*ampl!*eta_local
 Irhog = Irhog*denszero*ampl!*eta_local
 Rrhod = Rrhod*denszero*combined_amp !ampl*eta_local
 
 Irhod = 0. 
!--convert the data YJ in code units (v(code) = eta v(YJ)) 
 Rux   = Rux*combined_amp
 Iux   = Iux*combined_amp
 Ruy   = Ruy*combined_amp
 Iuy   = Iuy*combined_amp
 Ruz   = Ruz*combined_amp
 Iuz   = Iuz*combined_amp
 Rwx   = Rwx*combined_amp
 Iwx   = Iwx*combined_amp
 Rwy   = Rwy*combined_amp
 Iwy   = Iwy*combined_amp
 Rwz   = Rwz*combined_amp
 Iwz   = Iwz*combined_amp

 print *, 'V Coeffs ', Rux, Iux, Ruy, Iuy, Ruz, Iuz, Rwx, Iwx, Rwy, Iwy, Rwz, Iwz
 print *, 'Dens coeffs ', Rrhog, Irhog, Rrhod, Irhod
 print *, 'Dens amps ', sqrt(Rrhog**2 + Irhog**2), Rrhod

 print *, 'ntotal ', ntotal

 do i=1,ntotal
    xi      = x(1,i)
    zi      = x(2,i)
    cokx0   = cos(kx*xi)
    sikx0   = sin(kx*xi)
    cokz0   = cos(kz*zi)
    sikz0   = sin(kz*zi)
    if (idust.eq.1) dustfrac(1,i) = dust_to_gas_ratio
!-------------------------------------------------------------
! start with the gas particles
!-------------------------------------------------------------
    gdens = denszero
    ddens = denszerodust
    denszeroc = denszero
    
    if (itype(i).eq.itypegas) then
       if (iperturb.eqv. .true.) then

          if (idust.eq.1) then
             Rrhogc = Rrhog + Rrhod
             Irhogc = Irhog + Irhod
             denszeroc = denszero + denszerodust
          else
             Rrhogc = Rrhog
             Irhogc = Irhog
          endif
          
          select case(trim(lin))
          case('kxkzdensgas','kxkzdensgasdrag','kxkzdensgasdragpr','kxkzdens', &
               'kxkzdenswdrag','kxkzdenssdrag','kxkzdenssdragp','kxkzdensgasI', &
               'kxkzgasIdRvz','linAsimple','kxkzdensgasIP','prestest', &
               'linAnew','shalloweta','linA','smallcontrast','eigen1', &
               'eigen2','eigen4','eigennopres','contrasttest','mixedeigen', &
               'eigenNoPNoDrag','BBeigenNoPNoDrag','eigenNoShearNoPNoD', &
               'mixedeigen2')

       !setup a density perturbation Rrhog*cos(kx*x + kz*z) - Irhog*sin(kx*x + kz*z) 
             xstart  = xi 
             zstart  = zi
             !CAREFUL: keep 0.5 factor only for perturbations in 2D !
             !--term 1 : Rrhog*cos(kx*x)*cos(kz*z)
             pert = 0.5*Rrhogc*cokz0/denszeroc
             call deltarho_cos(xi,pert,kx,-hL,hL) !--pert in x
             x1   = xi - xstart
             xi   = xstart          
             pert = 0.5*Rrhogc*cokx0/denszeroc  
             call deltarho_cos(zi,pert,kz,-hL,hL) !--pert in z
             z1   = zi - zstart
             zi   = zstart

             if (.NOT.vparity) then
                !--term 2 : -Rrhogc*sin(kx*x)*sin(kz*z)
                pert = -0.5*Rrhogc*sikz0/denszeroc
                call deltarho_sin(xi,pert,kx,-hL,hL) !--pert in x
                x2   = xi - xstart
                xi   = xstart          
                pert = -0.5*Rrhogc*sikx0/denszeroc 
                call deltarho_sin(zi,pert,kz,-hL,hL) !--pert in z
                z2   = zi - zstart
                zi   = zstart
             else
                x2 = 0.
                z2 = 0.
             endif
             
             !--term 3 : -Irhog*sin(kx*x)*cos(kz*z)
             pert = -0.5*Irhogc*cokz0/denszeroc  
             if (vparity) pert = -1.*pert
             call deltarho_sin(xi,pert,kx,-hL,hL) !--pert in x
             x3   = xi - xstart
             xi   = xstart          
             pert = -0.5*Irhogc*sikx0/denszeroc
             if (vparity) pert = -1.*pert
             call deltarho_cos(zi,pert,kz,-hL,hL) !--pert in z
             z3   = zi - zstart
             zi   = zstart         
             
             if (.NOT.vparity) then
                !--term 4 : -Irhog*cos(kx*x)*sin(kz*z)
                pert = -0.5*Irhogc*sikz0/denszeroc
                call deltarho_cos(xi,pert,kx,-hL,hL) !--pert in x
                x4   = xi - xstart
                xi   = xstart          
                pert = -0.5*Irhogc*cokx0/denszeroc  
                call deltarho_sin(zi,pert,kz,-hL,hL) !--pert in z
                z4   = zi - zstart
                zi   = zstart          
             else
                x4 = 0.
                z4 = 0.
             endif

             !--sum all the contributions
             x(1,i) = xstart + x1 + x2 + x3 + x4
             x(2,i) = zstart + z1 + z2 + z3 + z4
             xi     = x(1,i)
             zi     = x(2,i)

          case default

          end select

          if (idust.eq.1) then
             cokx = cos(kx*xi + kz*zi)
             sikx = sin(kx*xi + kz*zi)
             
             gdens = denszero + (Rrhog*cokx - Irhog*sikx)
             ddens = denszerodust + (Rrhod*cokx - Irhod*sikx)
             
             dustfrac(1,i) = ddens/gdens
          endif

!----for adding a pert in kz only--------------------------- 
         !pert = 1.d-7
         !call deltarho_cos(zi,pert,kz,-hL,hL) !--pert in z
         !x(2,i) = zi    
!-----------------------------------------------------------------          
       endif
       
       if (trim(lin).eq.'kxvd' .or. trim(lin).eq.'kxvddrag') then
!----for adding a pert in kz only--------------------------- 
          pert = 1.d-7
          xi = x(1,i)
          call deltarho_cos(xi,pert,kx,-hL,hL) !--pert in x
          x(1,i) = xi
!-----------------------------------------------------------------          
       endif

       if (trim(lin).eq.'kzvd' .or. trim(lin).eq.'kzvddrag' .or. &
            trim(lin).eq.'kzvddragp') then
!----for adding a pert in kz only--------------------------- 
          pert = 1.d-7
          zi = x(2,i)
          call deltarho_cos(zi,pert,kz,-hL,hL) !--pert in z
          x(2,i) = zi
!-----------------------------------------------------------------          
       endif

       !--setup the kinematic quantities (v = vk + vNSH86 + vpert)
       vel(:,i) = 0.
       if (idust.eq.1) then
          vgx = vgxNSH
          vgy = 0.
          vgz = -domegadr*Omega0*xi + vgyNSH
          
          vdx = vdxNSH
          vdy = 0.
          vdz = -domegadr*Omega0*xi + vdyNSH 
          
          deltav(1,i) = vdx - vgx ! deltav = vdust - vgas
          deltav(2,i) = vdy - vgy
          deltav(3,i) = vdz - vgz
          
          rhodonrho = dustfrac(1,i)/(1.+dustfrac(1,i))
          
          vel(1,i) = vgx + rhodonrho*deltav(1,i) ! v = vg + rhod/rho*deltav
          vel(2,i) = vgy + rhodonrho*deltav(2,i)
          vel(3,i) = vgz + rhodonrho*deltav(3,i)

!          print *, 'vels ', vel(1,i), (gdens*vgx + ddens*vdx)/(gdens+ddens)
!--equivalent to above
!          vel(1,i) = (gdens*vgx + ddens*vdx)/(gdens+ddens)
!          vel(2,i) = (gdens*vgy + ddens*vdy)/(gdens+ddens)
!          vel(3,i) = (gdens*vgz + ddens*vdz)/(gdens+ddens)
       else
          vel(1,i) = vgxNSH
          if (trim(lin).eq.'eigenNoShearNoPNoD') then
             vel(3,i) = 0.
          else
             vel(3,i) = -domegadr*Omega0*xi + vgyNSH
          endif
       endif

       select case(trim(lin))
       case('kxkzdensgas','kxkzdensgasdrag','kxkzdensgasdragpr','kxkzdens', &
            'kxkzdenswdrag','kxkzdenssdrag','kxkzdenssdragp','kxkzdensgasIP', &
            'kxkzdensgasI')
          if(i.eq.1) print *, 'No velocity perturbation'
          
       case('kxepi','kxsimpledrag','kxstrongdrag','dpkx','dpkxstrong', &
            'kxvd','kxvddrag')
!----for adding a pert in kx only-------------------------------- 
          vel(1,i) = vel(1,i) + 1.d-7*cos(kx*xi)
          vel(2,i) = vel(2,i) + 0.*cos(kx*xi)
          vel(3,i) = vel(3,i) + 1.d-7*cos(kx*xi)
!----------------------------------------------------------------- 

       case('kxkzdrag','kxkzsimple','kxkzdragstrong','dpkxkz','dpkxkzstrong', &
            'kxkzgasIdRvz','linAsimple','prestest','linAnew','shalloweta', &
            'linA','smallcontrast','eigen1', &
               'eigen2','eigen4','eigennopres','contrasttest','mixedeigen', &
               'eigenNoPNoDrag','BBeigenNoPNoDrag','eigenNoShearNoPNoD', &
               'mixedeigen2')

          select case(trim(lin))
          case('linAsimple','linA','linAnew','shalloweta','prestest', &
               'smallcontrast','eigen1', &
               'eigen2','eigen4','eigennopres','contrasttest','mixedeigen', &
               'eigenNoPNoDrag','BBeigenNoPNoDrag','eigenNoShearNoPNoD', &
               'mixedeigen2')

          case default
             Iux = 0.
             Iuy = 0.
             Iuz = 0.
             Ruz = 0.
             Rux = 1.d-7
             Ruy = 1.d-7
             print *, 'NOT HERE'
          end select

          if (trim(lin).eq.'kxkzgasIdRvz') then 
             Ruz = 1.0E-7
             Iuz = 0.
!             vel(2,i) = vel(2,i) + Ruz*cos(kx*xi + kz*zi) - Iuz*sin(kx*xi + kz*zi)
          endif
!          else
             if (idust.eq.1) then
                vel(:,i) = 0.
                !--here x,y,z are just that.
                vgx = vgxNSH
                vgy = 0.
                vgz = -domegadr*Omega0*xi + vgyNSH
       
                vdx = vdxNSH
                vdy = 0.
                vdz = -domegadr*Omega0*xi + vdyNSH 
                
                if (iperturb.eqv. .true.) then
                   vgx = vgx + Rux*cos(kx*xi + kz*zi) - Iux*sin(kx*xi + kz*zi)
                   vdx = vdx + Rwx*cos(kx*xi + kz*zi) - Iwx*sin(kx*xi + kz*zi)
                   if (.NOT.vparity) then
                      vgy = vgy + Ruz*cos(kx*xi + kz*zi) - Iuz*sin(kx*xi + kz*zi)
                      vdy = vdy + Rwz*cos(kx*xi + kz*zi) - Iwz*sin(kx*xi + kz*zi)
                   else
                      vgy = vgy - (Ruz*cos(kx*xi + kz*zi) + Iuz*sin(kx*xi + kz*zi))
                      vdy = vdy - (Rwz*cos(kx*xi + kz*zi) + Iwz*sin(kx*xi + kz*zi))
                   endif
                   vgz = vgz + Ruy*cos(kx*xi + kz*zi) - Iuy*sin(kx*xi + kz*zi)
                   vdz = vdz + Rwy*cos(kx*xi + kz*zi) - Iwy*sin(kx*xi + kz*zi)
                   
                endif
                deltav(1,i) = vdx - vgx ! deltav = vdust - vgas
                deltav(2,i) = vdy - vgy
                deltav(3,i) = vdz - vgz
                
                rhodonrho = dustfrac(1,i)/(1.+dustfrac(1,i))

                vel(1,i) = vgx + rhodonrho*deltav(1,i) ! v = vg + rhod/rho*deltav
                vel(2,i) = vgy + rhodonrho*deltav(2,i)
                vel(3,i) = vgz + rhodonrho*deltav(3,i)
                print *, vel(3,i), deltav(3,i)

!                print *, 'vels ', vel(1,i), (gdens*vgx + ddens*vdx)/(gdens+ddens)
                

             else !--not one fluid
                vel(1,i) = vel(1,i) + Rux*cos(kx*xi + kz*zi) - Iux*sin(kx*xi + kz*zi)
                if (.NOT.vparity) then
                   vel(2,i) = vel(2,i) + Ruz*cos(kx*xi + kz*zi) - Iuz*sin(kx*xi + kz*zi)
                else
                   vel(2,i) = vel(2,i) - (Ruz*cos(kx*xi + kz*zi) + Iuz*sin(kx*xi + kz*zi))
                endif
                vel(3,i) = vel(3,i) + Ruy*cos(kx*xi + kz*zi) - Iuy*sin(kx*xi + kz*zi)
             endif
!          endif

       case('kxkzdustdens','kxkzdustdensp','kxkzdustdensp30_3', &
            'kxkzdustdensp30')
          !Do nothing

       case default
          if (iperturb .eqv. .true.) then
!----for adding a pert in kz only-------------------------------- 
             vel(1,i) = vel(1,i) + 1.d-7*cos(kz*zi)
             vel(2,i) = vel(2,i) + 0.*cos(kz*zi)
             vel(3,i) = vel(3,i) + 1.d-7*cos(kz*zi)        
!----------------------------------------------------------------- 
          endif
       end select

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
       if (iperturb .eqv. .true.) then

          select case(trim(lin))
          case('kxkzdens','kxkzdenswdrag','kxkzdenssdrag','kxkzdenssdragp', &
               'kxkzdensgasI','kxkzdensgasIP','linAsimple','prestest', &
               'kxkzdustdens','kxkzdustdensp','kxkzdustdensp30', &
               'kxkzdustdensp30_3','linAnew','shalloweta','linA', &
               'smallcontrast','eigen1', &
               'eigen2','eigen4','eigennopres','contrasttest','mixedeigen', &
               'mixedeigen2')
             
       !setup a density perturbation Rrhod*cos(kx*x + kz*z) - Irhod*sin(kx*x + kz*z) 
             xstart  = xi 
             zstart  = zi

             !--term 1 : Rrhod*cos(kx*x)*cos(kz*z)
             pert = 0.5*Rrhod*cokz0/denszerodust
             call deltarho_cos(xi,pert,kx,-hL,hL) !--pert in x
             x1   = xi - xstart
             xi   = xstart          
             pert = 0.5*Rrhod*cokx0/denszerodust  
             call deltarho_cos(zi,pert,kz,-hL,hL) !--pert in z
             z1   = zi - zstart
             zi   = zstart
             
             if (.NOT.vparity) then
                !--term 2 : -Rrhod*sin(kx*x)*sin(kz*z)
                pert = -0.5*Rrhod*sikz0/denszerodust
                call deltarho_sin(xi,pert,kx,-hL,hL) !--pert in x
                x2   = xi - xstart
                xi   = xstart          
                pert = -0.5*Rrhod*sikx0/denszerodust  
                call deltarho_sin(zi,pert,kz,-hL,hL) !--pert in z
                z2   = zi - zstart
                zi   = zstart
             else
                x2 = 0.
                z2 = 0.
             endif

             !--term 3 : -Irhod*sin(kx*x)*cos(kz*z)
             pert = -0.5*Irhod*cokz0/denszerodust
             if (vparity) pert = -1.*pert
             call deltarho_sin(xi,pert,kx,-hL,hL) !--pert in x
             x3   = xi - xstart
             xi   = xstart          
             pert = -0.5*Irhod*sikx0/denszerodust
             if (vparity) pert = -1.*pert
             call deltarho_cos(zi,pert,kz,-hL,hL) !--pert in z
             z3   = zi - zstart
             zi   = zstart         
            
             if (.NOT.vparity) then
                !--term 4 : -Irhod*cos(kx*x)*sin(kz*z)
                pert = -0.5*Irhod*sikz0/denszerodust
                call deltarho_cos(xi,pert,kx,-hL,hL) !--pert in x
                x4   = xi - xstart
                xi   = xstart          
                pert = -0.5*Irhod*cokx0/denszerodust 
                call deltarho_sin(zi,pert,kz,-hL,hL) !--pert in z
                z4   = zi - zstart
                zi   = zstart        
             else
                x4 = 0.
                z4 = 0.
             endif
             
             !--sum all the contributions
             x(1,i) = xstart + x1 + x2 + x3 + x4
             x(2,i) = zstart + z1 + z2 + z3 + z4
             xi     = x(1,i)
             zi     = x(2,i)

          end select
       endif

!----for adding a pert in kz only--------------------------- 
        ! pert = 1.d-7/denszerodust
        ! call deltarho_cos(zi,pert,kz,-hL,hL) !--pert in z
        ! x(2,i) = zi         
!-----------------------------------------------------------------        

    !--setup the kinematic quantities (v = vk + vNSH86 + vpert)
       vel(:,i) = 0.
       vel(1,i) = vdxNSH
       if (trim(lin).eq.'eigenNoShearNoPNoD') then
          vel(3,i) = 0.
       else
          vel(3,i) = -domegadr*Omega0*xi + vdyNSH 
       endif
       
!----for adding a pert in kz only-------------------------------- 
!         vel(1,i) = vel(1,i) + 1.d-7*cos(kz*zi)
!         vel(2,i) = vel(2,i) + 0.*cos(kz*zi)
!         vel(3,i) = vel(3,i) + 1.d-7*cos(kz*zi)        
!-----------------------------------------------------------------      
       
       
       if (iperturb.eqv. .true.) then
          select case(trim(lin))
          case('linAsimple','linA','linAnew','shalloweta','prestest', &
               'smallcontrast','eigen1', &
               'eigen2','eigen4','eigennopres','contrasttest','mixedeigen', &
               'mixedeigen2')
             vel(1,i) = vel(1,i) + Rwx*cos(kx*xi + kz*zi) - Iwx*sin(kx*xi + kz*zi)
             if (.NOT.vparity) then
                vel(2,i) = vel(2,i) + Rwz*cos(kx*xi + kz*zi) - Iwz*sin(kx*xi + kz*zi)
             else
                vel(2,i) = vel(2,i) - (Rwz*cos(kx*xi + kz*zi) + Iwz*sin(kx*xi + kz*zi))
             endif
             vel(3,i) = vel(3,i) + Rwy*cos(kx*xi + kz*zi) - Iwy*sin(kx*xi + kz*zi)
          end select
       endif
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
                      iperturb,cs0,kxfac,kzfac,Rrhod,Irhod,ampl,eta_local)
 use setup_params, only:pi
 implicit none
 character(len=20), intent(in) :: lin
 real, intent(in)     :: ampl,eta_local
 real, intent(out)    :: Bigkx,Bigkz,tausmode,dust_to_gas_ratio, cs0
 real, intent(out)    :: Rrhod,Irhod,kxfac,kzfac
 real, intent(out)    :: Rux,Iux,Ruy,Iuy,Ruz,Iuz,Rrhog,Irhog,Rwx,Iwx,Rwy,Iwy,Rwz,Iwz
 logical, intent(out) :: iperturb
 
 real :: Rurhog, Iurhog, Ruvrg, Iuvrg, Ruazig, Iuazig, Ruvzg, Iuvzg
 real :: Rurhod, Ruvrd, Iuvrd, Ruazid, Iuazid, Ruvzd, Iuvzd
 real :: Rsrhog, Isrhog, Rsvrg, Isvrg, Rsazig, Isazig, Rsvzg, Isvzg
 real :: Rsrhod, Rsvrd, Isvrd, Rsazid, Isazid, Rsvzd, Isvzd
 real :: sfac, ufac
 
 cs0    = 0.1   !--cfJY07
 kxfac = 1.0
 kzfac = 1.0
 Rrhod = 1.0

 tausmode          = 1.0
 Bigkx             = 1.
 Bigkz             = 1.
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

 case('damping')
    iperturb = .false.
    tausmode          = 1.0
    Bigkx             = 30.
    Bigkz             = 30.
    dust_to_gas_ratio = 1.d0
    
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
    Irhog = 0.0
    Rrhod = 0.0
    Bigkx = 1.0
    Bigkz = 1.0
    
 case('appBsimple')
    iperturb = .true.
    tausmode          = 3.d0
    Bigkx             = 30.
    Bigkz             = 30.
    dust_to_gas_ratio = 3.d0

 case('simpledrag','kxsimpledrag','k10asym')
    iperturb = .true.
    tausmode          = 1.d0
    Bigkx             = 1.
    Bigkz             = 1.
    dust_to_gas_ratio = 1.d0

 case('strongdrag','kxstrongdrag')
    iperturb = .true.
    tausmode          = 3.d-1
    Bigkx             = 1.
    Bigkz             = 1.
    dust_to_gas_ratio = 3.d0
     
 case('kxkzdrag','kxkzsimple','kxkzdragstrong','kxkzdensgas','kxkzdensgasdrag', &
      'kxkzdensgasdragpr','kxkzdens','kxkzdenswdrag','kxkzdenssdrag')
    iperturb = .true.
    tausmode          = 1.d0
    Bigkx             = 1.
    Bigkz             = 1.
    dust_to_gas_ratio = 1.d0
    Irhog = 0.
    Irhod = 0.
    Rrhog = 1.0E-7/ampl
    Rrhod = 1.0E-7/ampl

 case('kxkzdenssdragp')
    iperturb = .true.
    tausmode          = 1.d0
    Bigkx             = 1.
    Bigkz             = 1.
    dust_to_gas_ratio = 3.d0
    Irhog = 0.
    Irhod = 0.
    Rrhog = 1.0E-7/ampl
    Rrhod = 1.0E-7/ampl

 case('prestest')
    iperturb = .true.
    tausmode          = 1.d0
    Bigkx             = 30.
    Bigkz             = 30.
    dust_to_gas_ratio = 1.d0

    Rux   =  1.0E-7/(ampl*eta_local)
    Iux   =  0.
    Ruy   =  1.0E-7/(ampl*eta_local)
    Iuy   =  0.
    Ruz   =  1.0E-7/(ampl*eta_local)
    Iuz   =  0.
    Rwx   =  1.0E-7/(ampl*eta_local)
    Iwx   =  0.
    Rwy   =  0. !1.0E-7/(ampl*eta_local)
    Iwy   =  0.
    Rwz   =  1.0E-7/(ampl*eta_local)
    Iwz   =  0.

    Rrhog =  1.0E-7/ampl
    Rrhod =  1.0E-7/ampl
    Irhog =  0.
    Irhod =  0.

 case('kxkzdensgasI','kxkzgasIdRvz','kxkzdensgasIP')
    iperturb = .true.
    tausmode          = 1.d0
    Bigkx             = 1.
    Bigkz             = 1.
    dust_to_gas_ratio = 1.d0
    Irhog = 1.0E-7/ampl
    Irhod = 0.
    Rrhog = 0.
    Rrhod = 0.

 case('kxkzdustdens','kxkzdustdensp','kxkzdustdensp30')
    iperturb = .true.
    tausmode          = 1.d0
    Bigkx             = 1.
    Bigkz             = 1.
    dust_to_gas_ratio = 1.d0
    Irhog = 0.0
    Irhod = 0.
    Rrhog = 0.
    Rrhod = 1.0E-7/ampl

 case('kxkzdustdensp30_3')
    iperturb = .true.
    tausmode          = 1.d0
    Bigkx             = 1.
    Bigkz             = 1.
    dust_to_gas_ratio = 3.d0
    Irhog = 0.0
    Irhod = 0.
    Rrhog = 0.
    Rrhod = 1.0E-7/ampl

 case('kxkzvrgas')
    iperturb = .true.
    tausmode          = 1.d0
    Bigkx             = 1.
    Bigkz             = 1.
    dust_to_gas_ratio = 1.d0
    Irhog = 0.
    Irhod = 0.
    Rrhog = 0.
    Rrhod = 0.
    Rux = 1.0E-7/ampl

 case('dpkx','dpkxstrong','dpkxkz','kxvd','dpkz','kzvd','kzvddrag','kxvddrag', &
      'kzvddragp')
    iperturb = .true.
    tausmode          = 1.d0
    Bigkx             = 1.
    Bigkz             = 1.
    dust_to_gas_ratio = 1.d0

 case('dpkxkzstrong')
    iperturb = .true.
    tausmode          = 1.d0
    Bigkx             = 1.
    Bigkz             = 1.
    dust_to_gas_ratio = 3.d0

 case('kxepi')
    Bigkx             = 1.
    Bigkz             = 1.
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

 case('contrasttest')
    tausmode          = 0.1
    Bigkx             = 30.
    Bigkz             = 30.
    dust_to_gas_ratio = 3.d0
    Rux   =  1.0
    Iux   =  0.0
    Ruy   =  1.0
    Iuy   =  0.0
    Ruz   =  1.0
    Iuz   =  0.0
    Rrhod =  1.0E3
    Rrhog =  1.0E-4
    Irhog =  0.0
    Rwx   =  1.0
    Iwx   =  0.0
    Rwy   =  1.0
    Iwy   =  0.0
    Rwz   =  1.0
    Iwz   =  0.0

 case('mixedeigen')
    tausmode          = 0.1
    Bigkx             = 30.
    Bigkz             = 30.
    dust_to_gas_ratio = 3.d0

    !--unstable mode
    Rurhog = 1.5223E-9/0.001
    Iurhog = -5.54806E-9/0.001
    Ruvrg = 4.8981E-8/0.001
    Iuvrg =9.11494E-8/0.001
    Ruazig =8.18265E-9/0.001
    Iuazig = -2.14408E-8/0.001
    Ruvzg = -4.89808E-8/0.001
    Iuvzg =-9.11491E-8/0.001
    Rurhod = 0.001/0.001
    Ruvrd = 8.2769E-8/0.001
    Iuvrd = 1.0698E-7/0.001
    Ruazid = 9.88631E-9/0.001
    Iuazid = -2.0535E-8/0.001
    Ruvzd = -1.35957E-8/0.001
    Iuvzd =-8.64193E-8/0.001
    
    !--stable mode
    
    Rsrhog = 1.11889E-8/0.001
    Isrhog = -3.36682E-8/0.001
    Rsvrg = -2.49019E-6/0.001
    Isvrg =7.50562E-6/0.001
    Rsazig =5.30903E-6/0.001
    Isazig = 1.80086E-6/0.001
    Rsvzg = 2.49012E-6/0.001
    Isvzg =-7.50564E-6/0.001
    Rsrhod = 0.001/0.001
    Rsvrd = 5.62796E-7/0.001
    Isvrd = -2.5645E-6/0.001
    Rsazid = -1.8209E-6/0.001
    Isazid = -3.8673E-7/0.001
    Rsvzd = -4.93623E-7/0.001
    Isvzd =2.58506E-6/0.001

    !--combined mode
    ufac = 1.
    sfac = 1.

    Rrhog = ufac*Rurhog + sfac*Rsrhog 
    Irhog = ufac*Iurhog + sfac*Isrhog 
    Rux   = ufac*Ruvrg + sfac*Rsvrg 
    Iux   = ufac*Iuvrg + sfac*Isvrg 
    Ruy   = ufac*Ruazig + sfac*Rsazig 
    Iuy   = ufac*Iuazig + sfac*Isazig 
    Ruz   = ufac*Ruvzg + sfac*Rsvzg 
    Iuz   = ufac*Iuvzg + sfac*Isvzg 
    Rrhod = ufac*Rurhod + sfac*Rsrhod
    Rwx   = ufac*Ruvrd + sfac*Rsvrd
    Iwx   = ufac*Iuvrd + sfac*Isvrd
    Rwy   = ufac*Ruazid + sfac*Rsazid
    Iwy   = ufac*Iuazid + sfac*Isazid
    Rwz   = ufac*Ruvzd + sfac*Rsvzd
    Iwz   = ufac*Iuvzd + sfac*Isvzd

 case('mixedeigen2')
    tausmode          = 0.1
    Bigkx             = 30.
    Bigkz             = 30.
    dust_to_gas_ratio = 0.02

    Rrhog = 1.53514E-6/0.002
    Irhog = 5.10243E-7/0.002
    Rux   = 2.59029E-6/0.002
    Iux   = 1.55755E-6/0.002
    Ruy   = -1.46236E-7/0.002
    Iuy   = 2.41065E-6/0.002
    Ruz   = -2.59864E-6/0.002
    Iuz   = -1.56202E-6/0.002
    Rrhod = 0.002/0.002
    Rwx   = -0.000278517/0.002
    Iwx   = 0.0000494724/0.002
    Rwy   = -0.000136106/0.002
    Iwy   = 3.95301E-6/0.002
    Rwz   = 0.000115908/0.002
    Iwz   = 0.000046371/0.002

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

 case('linAnew')
    tausmode          = 0.1
    Bigkx             = 30.
    Bigkz             = 30.
    dust_to_gas_ratio = 3.d0

    Rux   =  4.898102782E-5
    Iux   =  9.114941884E-5
    Ruy   =  8.18265117E-6
    Iuy   = -2.144079186E-5
    Ruz   = -4.898075002E-5
    Iuz   = -9.114909012E-5

    Rrhog =  1.522304315E-6
    Irhog = -5.548056730E-6

    Rwx   =  8.276901438E-5
    Iwx   =  1.0698028972E-4
    Rwy   =  9.88630951E-6
    Iwy   = -2.053502206E-5
    Rwz   = -1.359568792E-5
    Iwz   = -8.641932839E-5

 case('smallcontrast')
    tausmode          = 0.1
    Bigkx             = 0.3
    Bigkz             = 0.3
    dust_to_gas_ratio = 0.02

    Rux   =  1.092719495E-3
    Iux   =  1.290464335E-3
    Ruy   = -7.765051011E-4
    Iuy   =  1.1729729655E-3
    Ruz   = -1.097581853E-3
    Iuz   = -1.293308458E-3

    Rrhog =  7.509276027E-4
    Irhog =  6.281679118E-4

    Rwx   = -8.027827660E-2
    Iwx   =  4.840970077E-2
    Rwy   = -8.560455755E-2
    Iwy   =  3.84962367E-3
    Rwz   = -1.0260193043E-3
    Iwz   = -4.880099103E-4

 case('eigennopres')
    tausmode          = 0.1
    Bigkx             = 0.3
    Bigkz             = 0.3
    dust_to_gas_ratio = 0.02

    Rrhog = 2.53800666E-4
    Irhog = 1.749600695E-3
    Rux   = 1.695283560E-3
    Iux   = 2.445408991E-3
    Ruy   = -1.971833775E-3
    Iuy   = 1.152944107E-3
    Ruz   = -1.691933816E-3
    Iuz   = -2.455341243E-3
    Rwx   = -2.631661687E-1
    Iwx   = -1.349376995E-1
    Rwy   = 6.75458140E-2
    Iwy   = -1.313957073E-1
    Rwz   = -1.226985812E-3
    Iwz   = 8.55478379E-4

 case('eigensettle')
    iperturb = .false.
    tausmode          = 0.1
    Bigkx             = 0.3
    Bigkz             = 0.3
    dust_to_gas_ratio = 0.02

 case('eigen1')
    tausmode          = 0.1
    Bigkx             = 0.3
    Bigkz             = 0.3
    dust_to_gas_ratio = 0.02

    Rrhog =  4.9520207E1
    Irhog =  2.512672306E3

    Rux   = -3.4690207
    Iux   = -1.778602126E2
    Ruy   =  3.331116002
    Iuy   = -6.4364639E-2
    Ruz   = -3.4680509
    Iuz   = -1.776103081E2

    Rwx   = -3.442851743
    Iwx   = -2.812640E-3
    Rwy   =  1.23615976E-3
    Iwy   = -1.2689989421E-1
    Rwz   = -3.428311581
    Iwz   =  1.447956E-3

 case('eigen2')
    tausmode          = 0.1
    Bigkx             = 0.3
    Bigkz             = 0.3
    dust_to_gas_ratio = 0.02

    Rrhog = -7.250499112E-4 
    Irhog =  1.1045893567E-3

    Rux   = -1.237748952E-4
    Iux   =  1.3060430501E-3 
    Ruy   = -8.769542444E-4
    Iuy   = -1.865336123E-4
    Ruz   =  1.353629433E-4
    Iuz   = -1.3124504129E-3

    Rwx   = -2.403232352E-1
    Iwx   = -2.337346807E-1
    Rwy   =  1.0131648605E-1
    Iwy   = -6.502347641E-2
    Rwz   = -6.411439358E-4
    Iwz   =  1.939945290E-4

 case('eigen4')
    tausmode          = 0.1
    Bigkx             = 0.3
    Bigkz             = 0.3
    dust_to_gas_ratio = 0.02

    Rrhog = 3.702201927E-4
    Irhog = -6.725950827E-4

    Rux   = 2.316714120E-3
    Iux   = -5.69471770E-4
    Ruy   = 1.265653046E-3
    Iuy   = 1.365519310E-3
    Ruz   = -2.320009864E-3
    Iuz   = 5.71196538E-4

    Rwx   = -1.976058901E-1
    Iwx   = 1.1917882E-3
    Rwy   = -5.070805943E-2
    Iwy   = 3.1167572E-4
    Rwz   = 1.980679035E-1
    Iwz   = -1.350593799E-1

 case('shalloweta')
    tausmode          = 0.1
    Bigkx             = 30.
    Bigkz             = 30.
    dust_to_gas_ratio = 3.d0

    Rux   =  1.482860776E-5
    Iux   =  8.838252077E-5
    Ruy   =  3.865733677E-5
    Iuy   = -1.405514535E-5
    Ruz   = -1.482860384E-5
    Iuz   = -8.838256919E-5

    Rrhog = -2.0128000E-8
    Irhog = -2.132422097E-6

    Rwx   =  2.692901894E-5
    Iwx   =  8.909573689E-5
    Rwy   =  3.875395153E-5
    Iwy   = -1.315135449E-5
    Rwz   = -2.70936412E-6
    Iwz   = -8.841160509E-5

 case('linAsimple')
    tausmode          = 0.1
    Bigkx             = 30.
    Bigkz             = 30.
    dust_to_gas_ratio = 1.d0
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

 case('eigenNoPNoDrag')
    Rux   = 1.0E-6
    Iux   = 0.0
    Ruy   = 0.0
    Iuy   = -7.071068060632488E-7  !-7.071068061E-7
    Ruz   = -1.000000140723876E-6 !-1.000000141E-6
    Iuz   = 0.0
    Rrhog = -3.751318379912935E-9 !-3.751318380E-9
    Irhog =  0.
    Rwx   =  0.
    Iwx   =  0.
    Rwy   =  0.
    Iwy   =  0.
    Rwz   =  0.
    Iwz   =  0.
    Rrhod =  0.

 case('eigenNoShearNoPNoD')
    Rrhog =  0.00001414213562
    Irhog =  0.
    Rux   =  4.0E-6
    Iux   =  0.
    Ruy   =  0.
    Iuy   =  0.
    Ruz   =  1.0E-6
    Iuz   =  0.
    Rwx   =  0.
    Iwx   =  0.
    Rwy   =  0.
    Iwy   =  0.
    Rwz   =  0.
    Iwz   =  0.
    Rrhod =  0.

 case('BBeigenNoPNoDrag')
    dust_to_gas_ratio = 3.0
!Big box (1./3.) case of above
    Rrhog = -0.00001102341800
    Irhog =  0.
    Rux   =  1.0E-6
    Iux   =  0.0
    Ruy   =  0.0
    Iuy   = -8.231702257E-7
    Ruz   = -1.710436882E-6
    Iuz   =  0.0
    Rwx   =  0.
    Iwx   =  0.
    Rwy   =  0.
    Iwy   =  0.
    Rwz   =  0.
    Iwz   =  0.
    Rrhod =  0.

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

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
