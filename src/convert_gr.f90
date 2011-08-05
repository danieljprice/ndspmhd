!
!--converts the setup from one co-ordinate system to another
!  obviously this only works for co-ordinate systems which map the same space.
!
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
 implicit none
 character(len=*), intent(in) :: geomold
 character(len=*), intent(inout) :: geom
 integer :: i,igeomold,igeom,icart
 real, dimension(ndim) :: xnew,vecnew

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
 if (geom(1:6) == 'cylrpz') then
    !!ibound(1) = 2 ! reflective in r
    xmin(1) = 0.0
    xmax(1) = 1.e10 ! a long way away
    if (ndim.ge.2) then 
       ibound(2) = 3 ! periodic in phi
       xmin(2) = -pi ! phi min
       xmax(2) = pi ! phi max
    endif
 endif 

end subroutine convert_setup

end module convert
