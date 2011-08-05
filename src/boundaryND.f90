!!-------------------------------------------------------------------------
!! update of particles at end of step for periodic (particles crossing
!! domain) and inflow/outflow (adds / subtracts particles from the domain)
!!
!! for ND case only periodic boundaries implemented
!!-------------------------------------------------------------------------
	 
subroutine boundary
 use dimen_mhd, only:ndim
 use debug, only:trace
 use loguns, only:iprint
 
 use bound, only:xmin,xmax
 use options, only:ibound
 use part, only:x,vel,npart
! use rates, only:force
! use part_in
! use setup_params
!
!--define local variables
!
 implicit none
 integer :: i,jdim,ncorrect
 logical :: debugging
!
!--allow for tracing flow
!      
 if (trace) then
    write(iprint,*) ' entering subroutine boundary'
    debugging = .true.
 endif   
 
!---------------------------------------------------------------------------
!  inflow/outflow when fixed particles are used.
!---------------------------------------------------------------------------     
 if (any(ibound.eq.1)) then
    if (ndim.gt.1) write(iprint,*) ' inflow/outflow not implemented in nd'
 endif
!----------------------------------------------------------------------
! periodic boundary conditions - allow particles to cross the domain
!----------------------------------------------------------------------
 if (any(ibound.eq.3)) then   
    do i=1,npart
       do jdim=1,ndim
          if (ibound(jdim).eq.3) then
	     if (x(jdim,i).gt.xmax(jdim)) then
!	        print*,'ss xold,xmax,xnew = ',jdim,x(jdim,i),xmax(jdim),xmin(jdim) + x(jdim,i) - xmax(jdim)
	        x(jdim,i) = xmin(jdim) + x(jdim,i) - xmax(jdim)
	     elseif(x(jdim,i).lt.xmin(jdim)) then
!	        print*,'ss xold,xmin,xnew = ',jdim,x(jdim,i),xmin(jdim),xmax(jdim) + x(jdim,i) - xmin(jdim)	     
!	        read*
                x(jdim,i) = xmax(jdim) - (xmin(jdim) - x(jdim,i))
             endif	  
	  endif
       enddo    
    enddo   
 elseif (any(ibound.gt.1)) then
!-------------------------------------------------------------------------
! other ghost boundaries - correct particles which have crossed boundary
!-------------------------------------------------------------------------
    ncorrect = 0
    do i=1,npart
       do jdim=1,ndim
          if (ibound(jdim).gt.1) then
             if (x(jdim,i).gt.xmax(jdim)) then
                !--move particle back inside boundary
                x(jdim,i) = xmax(jdim) - (x(jdim,i) - xmax(jdim))
                !--flip normal velocity component
                vel(jdim,i) = -vel(jdim,i)
!                force(jdim,i) = -force(jdim,i)
                ncorrect = ncorrect + 1
!                print*,'xmax',jdim,i,'xnew = ',x(jdim,i),vel(jdim,i)
	     elseif(x(jdim,i).lt.xmin(jdim)) then
!                print*,'xminold',jdim,i,'x,v = ',x(jdim,i),vel(jdim,i)
                !--move particle back inside boundary
                x(jdim,i) = xmin(jdim) + (xmin(jdim) - x(jdim,i))
                !--flip normal velocity component
                vel(jdim,i) = -vel(jdim,i)
!                force(jdim,i) = -force(jdim,i)
                ncorrect = ncorrect + 1
!                print*,'xminnew',jdim,i,'x,v = ',x(jdim,i),vel(jdim,i)
             endif
          endif
       enddo
    enddo
    if (ncorrect.gt.0) write(iprint,*) ' warning: ',ncorrect,' particles '// &
                                       ' adjusted from boundary'
 endif

 return 
end subroutine boundary
