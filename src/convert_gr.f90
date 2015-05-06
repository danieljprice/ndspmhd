!------------------------------------------------------------------------------!
! NDSPMHD: A Smoothed Particle (Magneto)Hydrodynamics code for (astrophysical) !
! fluid dynamics simulations in 1, 2 and 3 spatial dimensions.                 !
!                                                                              !
! (c) 2002-2015 Daniel Price                                                   !
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

!-----------------------------------------------------------------------------
!
!  converts the setup from one co-ordinate system to another
!  obviously this only works for co-ordinate systems which map the same space.
!
!-----------------------------------------------------------------------------
module convert
 implicit none
 
contains
 
subroutine convert_setup(geomold,geom)
 use dimen_mhd
 use debug
 use loguns
 
 use bound
 use options, only:ibound
 use part
 use setup_params, only:pi

 use geometry
 use grutils, only:metric_diag
 implicit none
 character(len=*), intent(in) :: geomold
 character(len=*), intent(inout) :: geom
 integer :: i,igeomold,igeom,icart
 real, dimension(ndim) :: xnew,vecnew
 real, dimension(ndimV) :: gdiag

 icart = 0
!
!--conversion must be between two recognisable co-ordinate systems
!
 print*,'old geometry = ',trim(geomold),' new geometry = ',trim(geom)
 select case(geomold(1:6))
   case('sphrpt')
      igeomold = 3
   case('cylrpz')
      igeomold = 2
   case('cartes')
      igeomold = 1
   case default
      write(iprint,*) 'ERROR: cannot convert co-ordinate systems (incompatible)'
      write(iprint,*) '=> Retaining original co-ordinate system'
      geom = geomold
      return
 end select
 select case(geom(1:6))
   case('sphrpt')
      igeom = 3
   case('cylrpz')
      igeom = 2
   case('cartes')
      igeom = 1
   case default
      write(iprint,*) 'ERROR: cannot convert co-ordinate systems (incompatible)'
      stop
 end select
!
!--if not to/from cartesians need to do an intermediate conversion
!  stop as this is not yet implemented (easy to do so if you want to...)
!
 if (igeomold.gt.1 .and. igeom.gt.1 .and. (igeom.ne.igeomold)) then
    write(iprint,*) 'ERROR in coord system conversion: intermediate step not implemented'
    stop
 endif

 write(iprint,*) 'CONVERTING setup from coord system ',igeomold,' to ',igeom
 do i=1,npart
    call coord_transform(x(:,i),ndim,igeomold,xnew(:),ndim,igeom)
    call vector_transform(x(:,i),vel(1:ndim,i),ndim,igeomold, &
                                 vecnew(1:ndim),ndim,igeom)                                    
    vel(1:ndim,i) = vecnew(1:ndim)
    
    call vector_transform(x(:,i),Bfield(1:ndim,i),ndim,igeomold, &
                                 vecnew(1:ndim),ndim,igeom)                                    
    Bfield(1:ndim,i) = vecnew(1:ndim)
    x(:,i) = xnew(:)
 enddo
 if (ndimV.gt.ndim) write(iprint,*)'WARNING: DOES NOT DO 2.5D YET'

 call coord_transform_limits(xmin,xmax,igeomold,igeom,ndim)
 select case(trim(geom))
 case('cylrpz')
    !!ibound(1) = 2 ! reflective in r
    if (ndim.ge.2) ibound(2) = 3 ! periodic in phi
 case('sphrpt')
    if (xmin(1) < tiny(xmin)) ibound(1) = 0 ! reflective in r
    if (ndim.ge.2) ibound(2) = 3 ! periodic in phi
    !--map densities from the cartesian space to spherical space
    
    xmin(1) = 0.1 !print*,' rmin = ',xmin(1)
    do i=1,npart
       call metric_diag(x(:,i),gdiag,sqrtg(i),ndim,ndimV,geom)
       x(1,i) = x(1,i)*sqrtg(i) !xmin(1) + (x(1,i) - xmin(1))**2/(xmax(1) - xmin(1))
       print*,i,'x = ',x(1,i)
    enddo
 end select

end subroutine convert_setup

end module convert
