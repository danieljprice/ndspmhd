!
!--converts the setup from one co-ordinate system to another
!
module convert
 implicit none
 
contains
 
subroutine convert_setup(igeomold,igeom)
 use dimen_mhd
 use debug
 use loguns
 
 use bound
 use options, only:ibound
 use part
 use setup_params, only:pi

 use geometry
 implicit none
 integer, intent(in) :: igeomold,igeom
 integer :: i
 real, dimension(ndim) :: xnew,vecnew

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
 if (igeom.eq.2) then
    !!ibound(1) = 2 ! reflective in r
    xmin(1) = 0.0
    xmax(1) = 10000.0 ! a long way away
    if (ndim.ge.2) then 
       ibound(2) = 3 ! periodic in phi
       xmin(2) = -pi ! phi min
       xmax(2) = pi ! phi max
    endif
 endif 

end subroutine convert_setup

end module convert
