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

!!-----------------------------------------------------------------
!!  Sets ghost particles if needed (ND)
!!  This version just loops over all the particles, and copies
!!  those which lie near the boundaries.
!!
!!  (previous version tried to use the link list to find particles
!!   near the boundary then created new link list cells)
!!
!!-----------------------------------------------------------------

subroutine set_ghost_particles
! use dimen_mhd
  use debug
  use loguns
  
  use bound
  use kernels,only:radkern
  use derivB
  ! use linklist
  use options
  use part, only:x,vel,hh,npart,ntotal,rho,uu,itype
  use setup_params, only:Omega0,domegadr
  use timestep, only:time
!
!--define local variables
!
  implicit none
  integer, parameter :: maxbound = 2 ! maximum no. of boundaries in each dim.
  integer, dimension(ndim) :: nbound ! number of boundaries in each dim.
  integer :: i
  integer :: imaxmin,imaxminprev,imaxminprevprev
  integer :: idimen,idimenprev,idimenprevprev
  integer :: jpart!,jtemp
  real, dimension(ndim) :: dxbound,xpart
  real, dimension(ndim,maxbound) :: xnew
  real, dimension(ndimV) :: vpart
  real :: dx,dxshift,xbound,xperbound,xtemp,dy,xmodbound,yi
  logical, dimension(ndim,maxbound) :: imakeghost
  real, dimension(maxbound) :: yoffset
  logical :: ireflect
!
!--allow for tracing flow
!      
  if (trace) write(iprint,*) ' Entering subroutine set_ghost_particles'
!
!--reset to zero ghost particles
!      
  ntotal = npart
  nbound(:) = 2   ! 2 boundaries in each dimension (ie. cartesian)
   ! this could be changed elsewhere
   ! must not be > maxbound
!  jtemp = 1861
! PRINT*,'jtemp = ',jtemp
!
!--use maximum value of h - check 2h is not bigger than box size
! 
  hhmax = maxval(hh(1:npart))
  dxbound(:) = radkern*hhmax ! this is the maximum distance away from boundary
! PRINT*,' hhmax, hhmin = ',hhmax,MAXLOC(hh(1:npart)),MINVAL(hh(1:npart)),MINLOC(hh(1:npart))
  if (any ((xmax(:)-xmin(:)) < dxbound)) then
     write(iprint,*) 'WARNING: ghost: hhmax too large = ',hhmax,maxloc(hh(1:npart))
  endif
!
!--shearing box ghosts
!
  if (ibound(1).eq.5) then
     if (ndim.gt.2) stop 'ghosts not implemented for 3D shearing box'
     do i=1,npart
        xpart(:) = x(:,i)
        vpart(:) = vel(:,i)
        !
        !--xmin boundary
        !
        dx = x(1,i) - xmin(1)
        if ((dx.lt.dxbound(1)).and.(dx.gt.0)) then
           xpart(1) = xmax(1) + dx
           xpart(2) = x(2,i) !- domegadr*Omega0*(xmax(1)-xmin(1))*time
           xpart(2) = xmodbound(xpart(2),xmin(2),xmax(2),0.)
           vpart(3) = vel(3,i) - domegadr*Omega0*(xmax(1)-xmin(1))
           call makeghost(i,xpart,vpart)
           yi = xpart(2)
           !--check ymin and ymax
           dy = yi - xmin(2)
           if ((dy.lt.dxbound(2)).and.(dy.gt.0)) then
              xpart(2) = xmax(2) + dy
              call makeghost(i,xpart,vpart)
           endif
           dy = xmax(2) - yi
           if ((dy.lt.dxbound(2)).and.(dy.gt.0)) then
              xpart(2) = xmin(2) - dy
              call makeghost(i,xpart,vpart)
           endif
        endif        
        !
        !--xmax boundary
        !        
        vpart(:) = vel(:,i)
        dx = xmax(1) - x(1,i)
        if ((dx.lt.dxbound(1)).and.(dx.gt.0)) then
           xpart(1) = xmin(1) - dx
           xpart(2) = x(2,i) !+ domegadr*Omega0*(xmax(1)-xmin(1))*time
           xpart(2) = xmodbound(xpart(2),xmin(2),xmax(2),0.)
           vpart(3) = vel(3,i) + domegadr*Omega0*(xmax(1)-xmin(1))
           call makeghost(i,xpart,vpart)
           yi = xpart(2)
           !--check ymin and ymax
           dy = yi - xmin(2)
           if ((dy.lt.dxbound(2)).and.(dy.gt.0)) then
              xpart(2) = xmax(2) + dy
              call makeghost(i,xpart,vpart)
           endif
           dy = xmax(2) - yi
           if ((dy.lt.dxbound(2)).and.(dy.gt.0)) then
              xpart(2) = xmin(2) - dy
              call makeghost(i,xpart,vpart)
           endif
        endif
        
        vpart(:) = vel(:,i)
        !--ymin boundary (just periodic)
        dy = x(2,i) - xmin(2)
        if ((dy.lt.dxbound(2)).and.(dy.gt.0)) then
           xpart(1) = x(1,i)
           xpart(2) = xmax(2) + dy
           call makeghost(i,xpart,vpart)
        endif        
        !--ymax boundary (just periodic)
        dy = xmax(2) - x(2,i)
        if ((dy.lt.dxbound(2)).and.(dy.gt.0)) then
           xpart(1) = x(1,i)
           xpart(2) = xmin(2) - dy
           call makeghost(i,xpart,vpart)
        endif
     enddo
     return
  endif
!
!--loop over all the particles (with link list could supply a list of 
!  particles within 2h of the boundary to search)
!
  over_part: do jpart=1,npart
     yoffset = 0.
!
!--if using reflecting ghosts, reflect if dx < 2*hi, as particle sees own ghost
!  for periodic must use hmax as particle sees different particles as ghosts
!  (for efficiency in this case should use max h of all particles near boundary)
!
     where (ibound.eq.2 .or. ibound.eq.4 .or. ibound.eq.6) dxbound(:) = radkern*hh(jpart)
!    IF (jpart.EQ.jtemp) PRINT*,jpart,' x = ',x(:,jpart)
    
     over_dimen: do idimen = 1, ndim  ! over spatial dimensions

        if (ibound(idimen).le.1) then
           imakeghost(idimen,:) = .false.
        else

           ireflect = .false.
      
           over_maxmin: do imaxmin = 1, nbound(idimen) ! over xmax, xmin  
!
!--set ghost position initially equal to particle position
!
              xpart(:) = x(:,jpart)
              vpart(:) = vel(:,jpart)
!
!--compute distance from boundary (dx)
!
              if (imaxmin.eq.1) then ! max
                 xbound = xmax(idimen) ! xbound is current boundary
                 xperbound = xmin(idimen) ! boundary where periodic ghosts go
                 dx = xmax(idimen) - x(idimen,jpart)
              else    ! min
                 xbound = xmin(idimen)
                 xperbound = xmax(idimen)
                 dx = x(idimen,jpart) - xmin(idimen)
              endif
!-----------------------------------------------------------------------------------
! make first ghost copy if dx < 2h, don't copy if on or over the boundary (dx <= 0)
!-----------------------------------------------------------------------------------
              imakeghost(idimen,imaxmin) = ((dx.lt.dxbound(idimen)).and.(dx.gt.0))
          
              if (imakeghost(idimen,imaxmin)) then
!
!--dxshift is now the amount to shift the particle from xbound/xperbound
!  (zero shift in other dimensions)
!      
                 dxshift = x(idimen,jpart) - xbound
!
!--xnew is the shifted position of the ghost particle
!  save this for each boundary to use for edges/corners
! 
                 if (ibound(idimen).eq.5) then ! shearing box
                    if (idimen.ne.1) stop 'error in shearing box'
                    xnew(1,imaxmin) = xperbound + dxshift
                    yoffset(imaxmin) = -SIGN(1.0,dxshift)*domegadr*Omega0*(xmax(1)-xmin(1))*time
                    vpart(2) = vel(2,jpart) - SIGN(1.0,dxshift)*domegadr*Omega0*(xmax(1)-xmin(1))
                    !print*,'HERE',vel(2,jpart),vpart(2),x(:,jpart),'xnew = ',xnew(:,imaxmin),SIGN(1.0,dxshift)
                    !read*
                 elseif (ibound(idimen).eq.3) then ! periodic
                    xnew(idimen,imaxmin) = xperbound + dxshift
                 elseif (ibound(idimen).eq.2 .or. ibound(idimen).eq.4 .or. ibound(idimen).eq.6) then ! reflective
                    xnew(idimen,imaxmin) = xbound - dxshift
                    vpart(idimen) = -vel(idimen,jpart)
                 endif
!
!--set ghost position in current dimension equal to shifted position
!      
                 xpart(idimen) = xnew(idimen,imaxmin)
                 if (ibound(idimen).eq.5) then
                    xpart(2) = x(2,jpart) + yoffset(imaxmin)
                    if (xpart(2).gt.(xmax(2)+dxbound(2))) then
                       xpart(2) = xpart(2) - (xmax(2)-xmin(2))
                    elseif (xpart(2).lt.(xmin(2)-dxbound(2))) then
                       xpart(2) = xpart(2) + (xmax(2)-xmin(2))
                    endif
                 endif
                 !if (jpart.eq.jtemp) then
                 !   print*,jpart,' dim=',idimen,'bnd=',imaxmin,' ghost = ',xpart
                 !endif
!
!--make the ghost particle at this position
!
                 call makeghost(jpart,xpart,vpart)

!---------------------------------------------------------------------
!  make additional ghost(s) if near an edge (up to 4 in 2D, 12 in 3D)     
!---------------------------------------------------------------------
! loop over previous dimensions
! (so if currently doing y, check y-x; if doing z, check z-y, z-x)
!
                 if (idimen.gt.1) then ! only in 2D or 3D
                    do idimenprev = 1,idimen-1 ! over previous dimension(s)
                       do imaxminprev=1,nbound(idimenprev) ! over boundaries
!
!--make an edge ghost if ghosts set in this previous dimension
!        
                          if (imakeghost(idimenprev,imaxminprev)) then
!
!--set ghost particle position in previous dimension to shifted position
!        
                             xpart(idimenprev) = xnew(idimenprev,imaxminprev)
                             xtemp = xpart(2)
                             if (ibound(idimenprev).eq.5) then
                                xpart(2) = xpart(2) + yoffset(imaxminprev)
                                if (xpart(2).gt.(xmax(2)+dxbound(2))) then
                                   xpart(2) = xpart(2) - (xmax(2)-xmin(2))
                                elseif (xpart(2).lt.(xmin(2)-dxbound(2))) then
                                   xpart(2) = xpart(2) + (xmax(2)-xmin(2))
                                endif
                             endif
                             !if (jpart.eq.jtemp) then
                             !   print*,'edge', idimen,idimenprev,'bnd=',imaxmin,imaxminprev,xpart
                             !endif
!
!--make the ghost particle at this position
!
                             !ireflect = .false.
                             if (ibound(idimenprev).eq.2 .or. ibound(idimenprev).eq.4 .or. ibound(idimenprev).eq.6) &
                                vpart(:) = vel(:,jpart)
                             if (ibound(idimenprev).eq.5) &
                                vpart(2) = vel(2,jpart) - SIGN(1.0,xpart(idimenprev)-x(1,jpart))*domegadr*Omega0*(xmax(1)-xmin(1))
                             call makeghost(jpart,xpart,vpart)

!----------------------------------------------------------------------
!  make a corner ghost if ghosts set in both previous dimensions
!  (could be up to 8 corner reflections of this particle in 3D
!----------------------------------------------------------------------
                             if (idimenprev.ge.2) then
                                idimenprevprev = idimenprev-1 ! in 3D
                                do imaxminprevprev = 1,nbound(idimenprevprev)
                                   if (imakeghost(idimenprevprev,imaxminprevprev)) then
! 
!--set ghost particle position in previous dimension to shifted position
!            
                                      xpart(idimenprevprev) = xnew(idimenprevprev,imaxminprevprev)
!
!--make the ghost particle at this position
!
                                      !if (jpart.eq.jtemp) then
                                      !   print*,'corner',idimen,idimenprev,idimenprevprev             
                                      !endif
                                      if (ibound(idimenprevprev).eq.2 .or. ibound(idimenprevprev).eq.4 &
                                          .or. ibound(idimenprevprev).eq.6) &
                                         vpart(:) = -vel(:,jpart)
                                      if (ibound(idimenprevprev).eq.5) then
                                         xpart(2) = xnew(2,imaxminprevprev)
                                         vpart(2) = vel(2,jpart) &
                                              - SIGN(1.0,xpart(idimenprevprev)-x(1,jpart))*domegadr*Omega0*(xmax(1)-xmin(1))
                                      endif

                                      call makeghost(jpart,xpart,vpart)
!--reset boundary
                                      xpart(idimenprevprev) = x(idimenprevprev,jpart)
                                      if (ibound(idimenprevprev).eq.5) xpart(2) = xnew(2,imaxminprev)

                                   endif ! made ghost in previous dimension
                                enddo ! over boundaries
                             endif ! if > 2D
!
!--reset position in idimenprev
!         
                             xpart(idimenprev) = x(idimenprev,jpart)      
                             if (ibound(idimenprev).eq.5) xpart(2) = xtemp

        
                          endif   ! made ghost in previous dimension
                       enddo ! over boundaries
                       
                    enddo ! over previous dimensions      
                 endif ! if > 1D
       
      
              endif ! make ghost for this boundary
           enddo over_maxmin
        endif
     enddo over_dimen
   
  enddo over_part
!
!--set type of ghost particles to be same as gas
  itype(npart+1:ntotal) = itype(ireal(npart+1:ntotal))
!
!--set unused elements of the array to zero (can cause errors in eos)
! 
  if (size(rho).gt.ntotal) then
     do i = ntotal+1,size(rho)
        rho(i) = 0.
        uu(i) = 0.
     enddo
  endif

  return
end subroutine set_ghost_particles

!
! subroutine to make a ghost particle
!
subroutine makeghost(jpart,xghost,vghost)
  use dimen_mhd
  use bound
  use loguns
  use part
  use eos, only:gamma
  use options, only:geom,ibound,imhd
  use mem_allocation, only:alloc
  implicit none
  integer, intent(in) :: jpart ! index of particle to be ghosted
  real, intent(in), dimension(ndim) :: xghost ! position of ghost particle
  real, intent(in), dimension(ndimV) :: vghost ! velocity of ghost particle
!  logical, intent(in) :: ireflect
  integer :: ipart
  real, parameter :: pi = 3.141592653589 
  real :: pri
!
!--create new particle and reallocate memory if needed
! 
  ipart = ntotal + 1
  if (ipart.gt.size(rho)) then
     write(iprint,*) 'ghost: ntotal > array size, re-allocating... '
     call alloc(max(size(rho)+1,int(0.9*(size(rho)+ntotal-jpart))))
  endif
  ntotal = ipart
!
!--set position of new ghost particle
!
  x(:,ipart) = xghost(:)
 !!print*,ipart,' ghost of ',jpart,' x = ',x(:,ipart),x(:,jpart)
!
!--copy particle properties
!
  call copy_particle(ipart,jpart)
  if (any(ibound.eq.6)) then
     pri = 2.5 - 0.1*dens(ipart)*x(2,ipart)
     uu(ipart) = pri/((gamma-1.)*dens(ipart))
     en(ipart) = uu(ipart)
  endif
!
!--for ghosts, copy everything including Bevol
!
  if (imhd.lt.0) then
     Bevol(:,ipart) = Bevol(:,jpart)
  endif
!
!--overwrite velocities if reflecting
!
!  if (ireflect) vel(:,ipart) = -vel(:,jpart) ! should reflect all vels
  vel(:,ipart) = vghost(:)
!
!--in cylindrical coords and ibound=4, shift angle by 180 degrees
!  
  if (geom(1:6).eq.'cylrpz' .and. ndim.ge.2 .and. ibound(1).eq.4) then
     x(2,ipart) = x(2,ipart) + pi
     if (x(2,ipart).gt.xmax(2)) x(2,ipart) = x(2,ipart) - 2.*pi
  endif
!
!--ireal for ghosts refers to the real particle of which they are ghosts
!  
  ireal(ipart) = jpart  

! PRINT*,'copying  old particle x(',jpart,'),vel,rho =', &
!        x(:,jpart),vel(:,jpart),rho(jpart)
! PRINT*,'creating new particle x(',ipart,'),vel,rho =', &
!        x(:,ipart),vel(:,ipart),rho(ipart)
 
  return
end subroutine makeghost

real function xmodbound(x,xmin,xmax,dxbound)
 implicit none
 real, intent(in) :: x,xmin,xmax,dxbound
 real :: deltax
 
 deltax = xmax - xmin + 2.*dxbound
 xmodbound = x
 if (x.lt.(xmin-dxbound)) then
    xmodbound = x + deltax - int((x - (xmin-dxbound))/deltax)*deltax
!    print*,'crossed min ',xmodbound,x
 elseif (x.gt.(xmax+dxbound)) then
    xmodbound = x - deltax - int((x - (xmax-dxbound))/deltax)*deltax
!    print*,'crossed max ',xmodbound,x
 endif

end function xmodbound
    
