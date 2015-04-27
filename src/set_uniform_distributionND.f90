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

module uniform_distributions
 implicit none

contains
!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  generic setup for a uniform density distribution of particles         !!
!!   in cartesian geometry for 1, 2 and 3 dimensions                      !!
!!                                                                        !!
!!  in 3D, lattice can be close packed, body centred or cubic or random   !!
!!     2D, lattice can be close packed or cubic or random                 !!
!!     1D, uniform or random                                              !!
!!                                                                        !!
!!------------------------------------------------------------------------!!

subroutine set_uniform_cartesian(idistin,psep,xmin,xmax, &
           fill,adjustbound,offset,perturb,psepx,psepy,psepz,rmin,rmax,mask,&
           stretchfunc,stretchdim)
!
!--include relevant global variables
!
 use dimen_mhd,      only:ndim,ndimV
 use debug,          only:trace
 use loguns,         only:iprint
 use linklist,       only:iamincell
 use bound,          only:hhmax
 use get_neighbour_lists, only:get_neighbour_list_partial
 use mem_allocation, only:alloc

 use options
 use part
!
!--define local variables
!  (note we read boundaries of region as input, so that more than one region
!  can be defined in the particle setup)
! 
 implicit none
 integer, intent(in) :: idistin
 real, dimension(ndim), intent(inout) :: xmin, xmax
 real, intent(in) :: psep
 real, intent(in), optional :: perturb,psepx,psepy,psepz,rmin,rmax
 logical, intent(in), optional :: offset,fill,adjustbound
 integer, intent(in), optional :: mask
 real, external, optional :: stretchfunc
 integer,        optional :: stretchdim
 
 integer :: i,j,k,ntot,npartin,npartx,nparty,npartz,ipart,iseed,imask
 integer :: idist,ineigh,nneigh,jpart,nskipx,nskipy,nskipz,istartx,istarty
 integer, dimension(1000) :: listneigh
 real :: xstart,ystart,zstart,deltax,deltay,deltaz
 real :: ran1,ampl,radmin,radmax,rr,rr2,phi,dangle,hpart
 real :: xi,xprev,dxmax,func,fracmassold,totmass,xminbisect,xmaxbisect,xminstretch !,dfunc
 real, dimension(ndim) :: xran,xcentre,dx
 logical :: adjustdeltas,offsetx
 integer :: idims,its
 integer, parameter :: maxits = 50
 real, parameter    :: tol = 1.e-8
 character(len=1)   :: chdim(3)
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine uniform_cartesian ',idistin
 write(iprint,*) 'uniform cartesian distribution '

 if (ndim.eq.1 .and. idistin.ne.4) then 
    idist = 1   ! in 1D use default version
 else
    idist = idistin
 endif
!
!--trim to radius
! 
 if (present(rmin)) then
    radmin = rmin
    write(iprint,*) 'rmin = ',radmin
 else
    radmin = -huge(radmin)
 endif
 if (present(rmax)) then
    radmax = rmax
    write(iprint,*) 'rmax = ',radmax
 else
    radmax = huge(radmax)
 endif
 if (present(mask)) then
    imask = mask
 else
    imask = 0
 endif

 npartin = npart    ! this is how many particles have already been set up
 
 select case(idist)
!
!--hexagonal close packed arrangement 
!  (22 is with symmetry-breaking)
!
 case(2,22)
!
!--determine number of particles
! 
    deltax = psep
    deltay = 0.5*sqrt(3.)*psep
    deltaz = sqrt(6.)/3.*psep
    if (present(psepx)) deltax = psepx
    if (present(psepy)) deltay = 0.5*sqrt(3.)*psepy
    if (present(psepz)) deltaz = sqrt(6.)/3.*psepz
    
    npartx = int(0.999*(xmax(1)-xmin(1))/deltax) + 1
    nparty = int(0.999*(xmax(2)-xmin(2))/deltay) + 1
    if (ndim.ge.3) then
       npartz = int(0.999*(xmax(3)-xmin(3))/deltaz) + 1
    else
       npartz = 1
    endif
    
!--for periodic boundaries, ymax needs to be divisible by 2
    if (ibound(2).eq.3 .and..not.present(rmax)) then
       nparty = 2*int(nparty/2)
       print*,' periodic boundaries: adjusting nparty = ',nparty
    endif
!--for periodic boundaries, ymax needs to be divisible by 2
    if (ndim.eq.3) then
       if (ibound(3).eq.3 .and..not.present(rmax)) then
          npartz = 3*int(npartz/3)
          print*,' periodic boundaries: adjusting npartz = ',npartz
       endif
    endif
!
!--adjust psep so that particles fill the volume
!
    print*,' npartx,y = ',npartx,nparty,npartz  !!,deltax,deltay
    print*,' partx = ',(xmax(1)-xmin(1))/deltax,(xmax(2)-xmin(2))/deltay  !!,deltax,deltay
    print*,' delta x,y initial  = ',deltax,deltay,deltaz
    adjustdeltas = .true.
    if (present(fill)) then
       adjustdeltas = fill
    endif
!
!--or adjust the boundaries appropriately
!
    if (present(adjustbound)) then
       if (adjustbound) then
          xmax(2) = xmin(2) + nparty*deltay
          print*,' adjusted y boundary : ymax  = ',xmax(2)
          if (ndim.ge.3) then
             xmax(3) = xmin(3) + npartz*deltaz
             print*,' adjusted z boundary : zmax  = ',xmax(3)
          endif
          print*,'xmin,max = ',xmin(:),xmax(:)
       endif
    endif
    
    if (present(offset)) then
       if (offset) then
          xcentre(:) = 0.5*(xmin + xmax)
          print*,'xcentre = ',xcentre
          npartx = int((xmax(1)-xcentre(1))/deltax) &
                 - int((xmin(1)-xcentre(1))/deltax) + 1
          nparty = int((xmax(2)-xcentre(2))/deltay) &
                 - int((xmin(2)-xcentre(2))/deltay) + 1
          if (ndim.ge.3) then
             npartz = int((xmax(3)-xcentre(3))/deltaz) &
                    - int((xmin(3)-xcentre(3))/deltaz) + 1
          endif
          write(iprint,*) &
           ' centred lattice : adjusting npart  = ',npartx,'x',nparty,'x',npartz
          !deltax = (xmax(1) - xmin(1))/npartx
          !deltay = (xmax(2) - xmin(2))/nparty
          write(iprint,*) '                 : adjusting deltax = ',deltax
          write(iprint,*) '                 : adjusting deltay = ',deltay
          
          xmin(1) = xcentre(1) - (npartx - 1)/2*deltax
          xmax(1) = xcentre(1) + (npartx - 1)/2*deltax
          xmin(1) = xmin(1) - 0.25*deltax
          xmax(1) = xmax(1) - 0.25*deltax

          xmin(2) = xcentre(2) - (nparty - 1)/2*deltay
          xmax(2) = xcentre(2) + (nparty - 1)/2*deltay
          xmin(2) = xmin(2) - 0.5*deltay
          xmax(2) = xmax(2) - 0.5*deltay

          if (ndim.ge.3) then
             xmin(3) = xcentre(3) - (npartz - 1)/2*deltaz
             xmax(3) = xcentre(3) + (npartz - 1)/2*deltaz
          endif
       endif
    elseif (adjustdeltas) then
       deltax = (xmax(1)-xmin(1))/(float(npartx))
       deltay = (xmax(2)-xmin(2))/(float(nparty))
       if (ndim.ge.3) deltaz = (xmax(3)-xmin(3))/(float(npartz))
       print*,' delta x,y,z adjusted = ',deltax,deltay,deltaz
    endif
!
!--allocate memory here
!
    ntot = npartx*nparty*npartz + npartin
    
    write(iprint,*) ' hexagonal close packed distribution, npart = ',ntot
    
    call alloc(ntot)

    ipart = npartin
    do k=1,npartz
       do j=1,nparty
          ystart = deltay/6.
          zstart = 0.5*deltaz
          xstart = 0.25*deltax
          if (mod(k,3).eq.0) then  ! 3rd layer
             ystart = ystart + 2./3.*deltay
             if (mod(j,2).eq.0) xstart = xstart + 0.5*deltax
          elseif (mod(k,3).eq.2) then ! 2nd layer
             ystart = ystart + 1./3.*deltay
             if (mod(j,2).eq.1) xstart = xstart + 0.5*deltax
          elseif (mod(j,2).eq.0) then
             xstart = xstart + 0.5*deltax
          endif
          !!xstart = 0.25*deltax
          !!if (mod(j,2).eq.0) xstart = xstart + 0.5*deltax
          do i = 1,npartx
             ipart = ipart + 1      !(k-1)*nparty + (j-1)*npartx + i
             x(1,ipart) = xmin(1) + (i-1)*deltax + xstart
             if (ndim.ge.2) x(2,ipart) = xmin(2) + (j-1)*deltay + ystart
             if (ndim.ge.3) x(3,ipart) = xmin(3) + (k-1)*deltaz + zstart
!--do not use if not within radial cuts
             rr = sqrt(dot_product(x(1:ndim,ipart),x(1:ndim,ipart))) 
             if (rr .lt. radmin .or. &
                 rr .gt. (radmax - 0.01*psep)) ipart = ipart - 1
             if (imask.ne.0) call applymask(imask,x(:,ipart),ipart)
          enddo
       enddo
    enddo
    
!
!--transform to asymmetric "glass-like" distribution 
!  by rotating lattice crystals by random amounts
!
    if (idist.eq.22 .and. ndim.ge.2) then
       hpart = 1.2*psep
       if (ndim.eq.3) then
          nskipz = int(4.5*hpart/deltaz)
          print*,'rotating every ',nskipz,'nd/rd z layer'
       else
          nskipz = 1
       endif
       nskipy = int(4.5*hpart/deltay)
       print*,'rotating every ',nskipy,'nd/rd y layer'
       nskipx = int(4.5*hpart/deltax)
       nskipy = nskipx
       print*,'rotating every ',nskipx,'nd/rd x layer'

   !--set "smoothing length" (size of lattice crystals)
       hh(1:ntot) = hpart
       hhmax = hpart
       npart = ipart
       ntotal = npart
       call set_linklist
   !--initialise random number generator for angle
       iseed = -630
       dangle = ran1(iseed)

       do k=nskipz,npartz,nskipz
          istarty = nskipy
          offsetx = .false.
          do j=istarty,nparty-nskipy,nskipy
             offsetx = .not.offsetx
             if (offsetx) then
                istartx = nskipx
             else
                istartx = nskipx + nskipx/2
             endif
             do i=istartx,npartx-nskipx,nskipx
                ipart = (j-1)*npartx + i - 1
                if (ndim.ge.3) ipart = ipart + (k-1)*nparty*npartx
                call get_neighbour_list_partial(iamincell(ipart),listneigh,nneigh)
                !print*,'ipart= ',ipart,' in cell ',iamincell(ipart),' nneighbours = ',nneigh
                !--get random angle for rotation of lattice points
                dangle = 2.*3.14159*(ran1(iseed) - 0.5)
                !print*,'rotating neighbours by ',dangle
                !--loop over neighbours             
                do ineigh = 1,nneigh
                   jpart = listneigh(ineigh)
                   dx(1:2) = x(1:2,jpart)-x(1:2,ipart)
                   rr2 = dot_product(dx,dx)
                   if (rr2.lt.4.*hpart**2 .and. ipart.ne.jpart) then
                      rr = sqrt(rr2)
                      phi = ATAN2(dx(2),dx(1))
                      phi = phi + dangle
                      dx(1) = rr*COS(phi)
                      dx(2) = rr*SIN(phi)
                      x(1:2,jpart) = x(1:2,ipart) + dx(1:2)
                   endif
                enddo
             enddo
          enddo
       enddo
      !--reset ipart (oops) 
       ipart = ntotal
    
    endif
!
!--body centred lattice
!
 case(3)
    if (ndim.eq.3) stop 'body centred not implemented in 3D'
!
!--determine number of particles
! 
    npartx = int((xmax(1)-xmin(1))/(sqrt(2.)*psep))
    nparty = int((xmax(2)-xmin(2))/(sqrt(2.)*psep))
    npartz = 1
!
!--adjust psep so that particles fill the volume
!
    deltax = sqrt(2.)*(xmax(1)-xmin(1))/(float(npartx)*sqrt(2.))
    deltay = 0.5*sqrt(2.)*(xmax(2)-xmin(2))/(float(nparty)*sqrt(2.))
!    print*,'psep = ',psepx,psepy,psep
!
!--allocate memory here
!
    ntot = 2*npartx*nparty
    write(iprint,*) 'body centred distribution, npart = ',ntot
    
    call alloc(ntot)
    npart = ntot

    ystart = 0.5*deltay      !psepy
    ipart = npartin
    do k=1,npartz
       do j=1,2*nparty
          xstart = 0.25*deltax
          if (mod(j,2).eq.0) xstart = xstart + 0.5*deltax
          do i = 1,npartx
             ipart = ipart + 1      !(k-1)*nparty + (j-1)*npartx + i
             x(1,ipart) = xmin(1) + (i-1)*deltax + xstart
             if (ndim.ge.2) x(2,ipart) = xmin(2) + (j-1)*deltay + ystart
!--do not use if not within radial cuts
             rr = sqrt(dot_product(x(1:ndim,ipart),x(1:ndim,ipart)))
             if (rr .lt. radmin .or. rr .gt. radmax) ipart = ipart - 1
             if (imask.ne.0) call applymask(imask,x(:,ipart),ipart)
          enddo
       enddo
    enddo 
!
!--random particle distribution
!  (uses random number generator ran1 in module random)
! 
 case(4)
    ntot = int(product((xmax(:)-xmin(:))/psep))
    npart = ntot
    write(iprint,*) 'random particle distribution, npart = ',ntot 
    call alloc(ntot)
!
!--initialise random number generator
!    
    iseed = -87682
    write(iprint,*) ' random seed = ',iseed
    xran(1) = ran1(iseed)
    
    ipart = npartin
    do i=1,ntot
       ipart = ipart + 1
       do j=1,ndim
          xran(j) = ran1(iseed)
       enddo
       x(:,ipart) = xmin(:) + xran(:)*(xmax(:)-xmin(:))
!--do not use if not within radial cuts
       rr = sqrt(dot_product(x(1:ndim,ipart),x(1:ndim,ipart)))
       if (rr .lt. radmin .or. rr .gt. radmax) ipart = ipart - 1
    enddo
!
!--random particle distribution
!  (uses random number generator ran1 in module random)
! 
 case(5)
     ipart = npartin
     ntot = int(product((xmax(:)-xmin(:))/psep))
     ipart = ipart + ntot
     write(iprint,*) 'quasi-random (sobol) particle distribution, npart = ',ntot 
     call alloc(ipart)
!
!--initialise quasi-random sequence
!
     !call sobolsequence(-ndim,xran(:))
     call sobseq(-ndim,xran(:))
     
     do i=1,ntot
        call sobseq(ndim,xran(:))
        x(:,i) = xmin(:) + xran(:)*(xmax(:)-xmin(:))
!        print*,i,xran(:),' x = ',x(:,i)
!	read*
     enddo

 case default
!----------------------
!  cubic lattice
!----------------------

    if (present(offset)) then
       if (offset) then
          write(iprint,*) 'offset lattice'
          xmin(:) = xmin(:) + 0.5*psep
          xmax(:) = xmax(:) + 0.5*psep
       endif
    endif
!
!--determine number of particles
! 
    deltax = psep
    deltay = psep
    deltaz = psep
    if (present(psepx)) deltax = psepx
    if (present(psepy)) deltay = psepy
    if (present(psepz)) deltaz = psepz
    
    npartx = int((xmax(1)-xmin(1))/deltax)    
    if (ndim.ge.2) then
       nparty = int((xmax(2)-xmin(2))/deltay)    
    else
       nparty = 1
    endif
    if (ndim.ge.3) then
       npartz = int((xmax(3)-xmin(3))/deltaz)    
    else
       npartz = 1    
    endif
!
!--adjust psep so that particles fill the volume
!
    print*,' npartx,y,z = ',npartx,nparty,npartz !!,deltax,deltay
    print*,' delta x,y,z initial  = ',deltax,deltay,deltaz
    if (present(fill)) then
       if (fill) then
          deltax = (xmax(1)-xmin(1))/(float(npartx))
          if (ndim.ge.2) deltay = (xmax(2)-xmin(2))/(float(nparty))
          if (ndim.ge.3) deltaz = (xmax(3)-xmin(3))/(float(npartz))
          print*,' delta x,y,z adjusted = ',deltax,deltay,deltaz
       endif
    endif
!
!--allocate memory here
!
    ntot = npartx*nparty*npartz
    ipart = npartin
    npart = npartin + ntot   ! add to particles already setup
    write(iprint,*) 'cubic lattice, adding ',ntot,', npart = ',npart,npartx,nparty,npartz
    write(iprint,*) 'xmin = ',xmin, ' xmax = ',xmax
    write(iprint,*) 'psep = ',psep
    call alloc(npart)
    
    do k=1,npartz
       do i=1,npartx
          do j = 1,nparty
             ipart = ipart + 1   !(k-1)*nparty + (j-1)*npartx + i
             x(1,ipart) = xmin(1) + (i-1)*deltax + 0.5*deltax
             if (ndim.ge.2) then
                x(2,ipart) = xmin(2) + (j-1)*deltay + 0.5*deltay
                if (ndim.ge.3) then
                   x(3,ipart) = xmin(3) + (k-1)*deltaz + 0.5*deltaz
                endif
             endif
!--do not use if not within radial cuts
             rr = sqrt(dot_product(x(1:ndim,ipart),x(1:ndim,ipart)))
             if (rr .lt. radmin .or. rr .gt. radmax) ipart = ipart - 1
             if (imask.ne.0) call applymask(imask,x(:,ipart),ipart)
         enddo
       enddo
    enddo
 end select 
!
!--if perturb is set add a small random perturbation to the particle positions
! 
 if (present(perturb)) then
    write(iprint,"(a,f6.1,a)") &
       ' random perturbation of ',perturb*100,'% of particle spacing'
!
!--initialise random number generator
!    
    iseed = -268
    xran(1) = ran1(iseed)
    ampl = perturb*psep   ! perturbation amplitude
    
    do i=npartin+1,ipart  ! apply to new particles only
       do j=1,ndim
          xran(j) = ran1(iseed)
       enddo
       x(:,i) = x(:,i) + ampl*(xran(:)-0.5)
    enddo
 endif

!------------------------------------------------------------------------
! STRETCH distribution in one direction to map to a 1D density function
!------------------------------------------------------------------------
 if (present(stretchfunc)) then
    if (present(stretchdim)) then
       idims = stretchdim
    else
       idims = 1
    endif
    chdim = (/'x','y','z'/)
    print*,' stretching in '//chdim(idims)//' direction to match density function'
    
    !
    ! pivot point for stretch transformation
    ! use the origin if limits straddle it, otherwise pivot about xmin
    !
    if (xmax(idims) > 0. .and. xmin(idims) < 0.) then
       xminstretch = 0.
    else
       xminstretch = xmin(idims)
    endif
    
    !
    ! preliminaries 
    !
    totmass  = get_mass(stretchfunc,xmax(idims),xminstretch)
    dxmax    = xmax(idims) - xminstretch
    !
    ! solve for each particle position by iteration
    ! use bisection method as guarantees convergence
    !
    do i=npartin+1,ipart
       xi   = x(idims,i)
       fracmassold = (xi - xminstretch)/dxmax
       if (abs(fracmassold) > 0.) then
          xprev = 0.
          its = 0
          func = 1.e2
          xminbisect = xmin(idims)
          xmaxbisect = xmax(idims)
          do while(abs(func) > tol .and. its < maxits)
             xi = 0.5*(xminbisect + xmaxbisect)
             !xprev  = xi
             its    = its + 1
             func   = get_mass(stretchfunc,xi,xminstretch)/totmass - fracmassold
             if (func > 0.) then
                xmaxbisect = xi
             else
                xminbisect = xi
             endif
             !dfunc  = stretchfunc(xi)/totmass
             !xi     = xprev - func/dfunc ! Newton-Raphson
          enddo
          if (its >= maxits) then
             write(iprint,*) 'ERROR: iterations not converged during stretch'
          endif
          x(idims,i) = xi
       endif
    enddo
 endif
!
!--set total number of particles = npart
!
 write(iprint,*) ' Number of particles added = ',ipart - npartin
 ntotal = ipart
 npart = ntotal
!--only allocate memory for what is needed
 if (any(ibound.ge.1)) then
    call alloc(int(1.1*ntotal))
 else
    call alloc(ntotal)
 endif
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine uniform_cartesian'
            
 return
end subroutine

!----------------------------------------------------------------
!     Set up a uniform density sphere of particles
!     centred on the origin. This subroutine sets up
!     a uniform cube of particles and trims it to be spherical.
!----------------------------------------------------------------

subroutine set_uniform_spherical(idist,rmax,rmin,perturb,centred,trim)
!
!--include relevant global variables
!
 use dimen_mhd
 use debug
 use loguns
 use part
 use setup_params, only:psep
 use mem_allocation, only:alloc
!
!--define local variables
!            
 implicit none
 integer, intent(in) :: idist
 real, intent(in) :: rmax
 real, intent(in), optional :: rmin
 real, intent(in), optional :: perturb
 real, intent(in), optional :: trim
 logical, intent(in), optional :: centred

 integer :: i,ierr,ntemp 
 integer, dimension(:), allocatable :: partlist
 real :: rad,radmin,radmax
 real, dimension(ndim) :: xmin,xmax
 real, dimension(:,:), allocatable :: xtemp
 logical :: offset
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(sphdis)'
!
!--check for errors
!
 if (rmax.le.0.) stop 'error: set_uniform_spherical: rmax <= 0'
!
!--initially set up a uniform density cartesian grid
!
 xmin(:) = -rmax
 xmax(:) = rmax
 if (present(centred)) then
    offset = centred
 else
    offset = .false.
 endif

  
 if (present(perturb)) then
    call set_uniform_cartesian(idist,psep,xmin,xmax,perturb=perturb,offset=offset,rmax=rmax)
 else
    call set_uniform_cartesian(idist,psep,xmin,xmax,offset=offset,rmax=rmax) 
 endif

 if (present(rmin)) then
    if (present(trim)) then
       radmin = rmin + abs(trim)
    else
       radmin = rmin
    endif
 else
    radmin = 0.
 endif
 if (present(trim)) then
    radmax = rmax - abs(trim)
 else
    radmax = rmax
 endif
!
!--construct list of particles with r < rmax
!  
 allocate ( partlist(npart), stat=ierr )
 if (ierr.ne.0) stop 'error allocating memory in uniform_spherical'
 
 ntemp = 0 ! actual number of particles to use
 do i=1,npart
    rad = sqrt(dot_product(x(:,i),x(:,i)))
    if (rad.lt.radmax .and. rad.ge.radmin) then
       ntemp = ntemp + 1
       partlist(ntemp) = i
    endif   
 enddo 
!
!--reorder particles
!
 write(iprint,*) 'npart (cropped) = ',ntemp, ' cropping factor = ',ntemp/real(npart)
 npart = ntemp
 ntotal = npart
 if (allocated(xtemp)) deallocate(xtemp)
 allocate(xtemp(ndim,npart))
 
 do i=1,npart
    xtemp(:,i) = x(:,partlist(i))
 enddo

 x(:,1:npart) = xtemp(:,1:npart)

 deallocate(xtemp)
 
!
!--reallocate memory to new size of list
!
 call alloc(npart)

 return
end subroutine

subroutine reset_centre_of_mass(x,pmass)
 use dimen_mhd
 use debug
 use loguns
 implicit none
 real, dimension(:,:), intent(inout) :: x
 real, dimension(:), intent(in) :: pmass
 integer :: i,npart
 real, dimension(ndim) :: xcentre
 
 if (trace) write(iprint,*) 'entering subroutine reset_centre_of_mass'
 
 npart = size(pmass)
 
 xcentre = 0.
 do i=1,npart
    xcentre(:) = xcentre(:) + pmass(i)*x(:,i)
 enddo
 
 xcentre = xcentre/SUM(pmass)
 
 write(iprint,*) 'centre of mass is at ',xcentre(:)
 write(iprint,*) 'resetting to zero'
 
 do i=1,npart
    x(:,i) = x(:,i) - xcentre(:)
 enddo
 
end subroutine reset_centre_of_mass

!
!--this subroutine applies a variety of masks to the uniform setup
!  by defining a new mask, it gives a way of setting up density distributions
!  for equal mass particles
!
!  if imask is negative, this should supply the inverse mask
!
subroutine applymask(imask,xpart,ipart)
 use setup_params, only:pi
 implicit none
 integer, intent(in) :: imask
 real, dimension(:), intent(in) :: xpart
 integer, intent(inout) :: ipart
 real, dimension(size(xpart)) :: dx,xorigin
 real :: radius1,radius2,phi,z,Rfunc,Rfunc2
 
 Rfunc(z) = -50. + 10.*cos(2.*pi*z/100.)
 Rfunc2(z) = -0.01*cos(6.*pi*z)
 
 select case(imask)
 case(1,-1)
!
!   TWO CYLINDRICAL BLAST WAVES
!
    xorigin(:) = 0.0
    dx(:) = xpart(:) - xorigin(:)
    radius1 = sqrt(dot_product(dx,dx))
    xorigin(:) = 1.0
    dx(:) = xpart(:) - xorigin(:)
    radius2 = sqrt(dot_product(dx,dx))
    if (imask.gt.0) then
!
!--cut around two radial points
!
       if (radius1.gt.0.4 .and. radius2.gt.0.4) then
          ipart = ipart - 1
       endif
!
!--inverse of above
!
    else
       if (radius1.le.0.4 .or. radius2.le.0.4) then
          ipart = ipart - 1
       endif    
    endif
 case(2,-2)
!
!   DENSE LAYER
!
    if (imask.gt.0 .and. abs(xpart(2)).lt.0.25) then
       ipart = ipart - 1
    elseif (abs(xpart(2)).gt.0.25) then
       ipart = ipart - 1
    endif
 case(3,-3)
!
!   FUNNY SHAPE BLOB
!
    xorigin(:) = 0.0
    dx(:) = xpart(:) - xorigin(:)
    radius1 = sqrt(dot_product(dx,dx))
    phi = atan2(xpart(2),xpart(1))

    if (imask.gt.0) then
!
!--cut around funny blob shape
!
       if (radius1.gt.0.4*(1.+0.2*sin(0.2*phi/(2.*pi)))) then
          ipart = ipart - 1
       endif
!
!--inverse of above
!
    else
       if (radius1.le.0.4*(1.+0.2*sin(0.2*phi/(2.*pi)))) then
          ipart = ipart - 1
       endif
    endif
 case(4,-4)
 
!
!--NOVA test
!
    if (imask.gt.0) then
       if (xpart(1).lt.Rfunc(xpart(2))) then
          ipart = ipart - 1
       endif
    else
       if (xpart(1).gt.Rfunc(xpart(2))) then
          ipart = ipart - 1
       endif
    endif
 case(5,-5)
!
!--Rayleigh Taylor interface
!
    if (imask.gt.0) then
       if (xpart(2).le.Rfunc2(xpart(1))) then
          ipart = ipart - 1
       endif
    else
       if (xpart(2).gt.Rfunc2(xpart(1))) then
          ipart = ipart - 1
       endif
    endif

 end select
 
end subroutine applymask

!----------------------------------------------------
!+
!  Integrate to get total mass as function of x
!  used when stretching distribution to match
!  1D density profile
!+
!----------------------------------------------------
real function get_mass(rhofunc,x,xmin)
 integer, parameter :: ngrid = 1024
 real, intent(in) :: x,xmin
 real, external   :: rhofunc
 real :: dx,xi,dmi,dmprev
 integer :: i

 dx = (x - xmin)/real(ngrid)
 dmprev   = 0.
 get_mass = 0.
 do i=1,ngrid
    xi     = xmin + i*dx
    dmi    = rhofunc(xi)*dx
    get_mass = get_mass + 0.5*(dmi + dmprev) ! trapezoidal rule
    dmprev = dmi
 enddo

end function get_mass

end module uniform_distributions
