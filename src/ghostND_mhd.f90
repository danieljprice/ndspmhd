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
  use part, only:x,hh,npart,ntotal,rho,uu,itype
!
!--define local variables
!
  implicit none
  integer, parameter :: maxbound = 2 ! maximum no. of boundaries in each dim.
  integer, dimension(ndim) :: nbound ! number of boundaries in each dim.
  integer :: i
  integer :: imaxmin,imaxminprev,imaxminprevprev
  integer :: idimen,idimenprev,idimenprevprev
  integer :: jpart,jtemp
  real, dimension(ndim) :: dxbound,xpart
  real, dimension(ndim,maxbound) :: xnew
  real :: dx,dxshift,xbound,xperbound
  logical, dimension(ndim,maxbound) :: imakeghost
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
!--loop over all the particles (with link list could supply a list of 
!  particles within 2h of the boundary to search)
!
  over_part: do jpart=1,npart
!
!--if using reflecting ghosts, reflect if dx < 2*hi, as particle sees own ghost
!  for periodic must use hmax as particle sees different particles as ghosts
!  (for efficiency in this case should use max h of all particles near boundary)
!
     where (ibound.eq.2 .or. ibound.eq.4) dxbound(:) = radkern*hh(jpart)
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
                 if (ibound(idimen).eq.3) then ! periodic
                    xnew(idimen,imaxmin) = xperbound + dxshift
                 elseif (ibound(idimen).eq.2 .or. ibound(idimen).eq.4) then ! reflective
                    ireflect = .true.
                    xnew(idimen,imaxmin) = xbound - dxshift
                 endif
!
!--set ghost position in current dimension equal to shifted position
!      
                 xpart(idimen) = xnew(idimen,imaxmin)
                 !if (jpart.eq.jtemp) then
                 !   print*,jpart,' dim=',idimen,'bnd=',imaxmin,' ghost = ',xpart
                 !endif
!
!--make the ghost particle at this position
!
                 call makeghost(jpart,xpart,ireflect)

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
                             !if (jpart.eq.jtemp) then
                             !   print*,'edge', idimen,idimenprev,'bnd=',imaxmin,imaxminprev,xpart
                             !endif
!
!--make the ghost particle at this position
!
                             if (ibound(idimenprev).eq.2) ireflect = .true.
                             call makeghost(jpart,xpart,ireflect)

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
                                      if (ibound(idimenprevprev).eq.2) ireflect = .true.
                                      call makeghost(jpart,xpart,ireflect)
!--reset boundary
                                      xpart(idimenprevprev) = x(idimenprevprev,jpart)
                                   endif ! made ghost in previous dimension
                                enddo ! over boundaries
                             endif ! if > 2D
!
!--reset position in idimenprev
!         
                             xpart(idimenprev) = x(idimenprev,jpart)      

        
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
!--set type of ghost particles to zero
!
  itype(npart+1:ntotal) = 0
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
subroutine makeghost(jpart,xghost,ireflect)
  use dimen_mhd
  use bound
  use loguns
  use part
  use options, only:geom,ibound
  use mem_allocation, only:alloc
  implicit none
  integer, intent(in) :: jpart ! index of particle to be ghosted
  real, intent(in), dimension(ndim) :: xghost ! position of ghost particle
  logical, intent(in) :: ireflect
  integer :: ipart
  real, parameter :: pi = 3.141592653589 
!
!--create new particle and reallocate memory if needed
! 
  ipart = ntotal + 1
  if (ipart.gt.size(rho)) then
     write(iprint,*) 'ghost: ntotal > array size, re-allocating... '
     call alloc(size(rho)+1)
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
!
!--overwrite velocities if reflecting
!
  if (ireflect) vel(:,ipart) = -vel(:,jpart) ! should reflect all vels
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

      
