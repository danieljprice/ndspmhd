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
 use part_in, only:velin,xin
 use setup_params, only:Omega0,domegadr
 use rates, only:force
! use part_in
! use setup_params
!
!--define local variables
!
 implicit none
 integer :: i,jdim,ncorrect,ncross
 logical :: debugging
 real :: xmodbound
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
    !if (ndim.gt.1) write(iprint,*) ' inflow/outflow not implemented in nd'
 endif
!----------------------------------------------------------------------
! periodic boundary conditions - allow particles to cross the domain
!----------------------------------------------------------------------
 if (any(ibound.eq.3)) then
    ncross = 0 
    do i=1,npart
       do jdim=1,ndim
          if (ibound(jdim).eq.3) then
             if (x(jdim,i).gt.xmax(jdim)) then
!                print*,'ss xold,xmax,xnew = ',jdim,x(jdim,i),xmax(jdim),xmin(jdim) + x(jdim,i) - xmax(jdim)
                x(jdim,i) = xmin(jdim) + x(jdim,i) - xmax(jdim)
                ncross = ncross + 1
             elseif(x(jdim,i).lt.xmin(jdim)) then
!                print*,'ss xold,xmin,xnew = ',jdim,x(jdim,i),xmin(jdim),xmax(jdim) + x(jdim,i) - xmin(jdim)             
!                read*
                x(jdim,i) = xmax(jdim) - (xmin(jdim) - x(jdim,i))
                ncross = ncross + 1
             endif          
          endif
       enddo
    enddo
    if (ncross.gt.0 .and. ibound(1).eq.5) &
       write(iprint,*) ncross,' particles crossing y (shearing box)'
 endif
!----------------------------------------------------------------------------
! shearing box boundary in the x (=r) direction
!----------------------------------------------------------------------------
 if (ibound(1).eq.5 .and. ndim.ge.2) then
    ncross = 0
    do i=1,npart
       if (x(1,i).gt.xmax(1)) then
          x(1,i) = x(1,i) - (xmax(1)-xmin(1))
          x(2,i) = x(2,i) !+ domegadr*Omega0*(xmax(1)-xmin(1))*time
          x(2,i) = xmodbound(x(2,i),xmin(2),xmax(2),0.)
!          if (x(2,i).gt.xmax(2)) then
!             x(2,i) = xmin(2) + (x(2,i)-xmax(2))
!          endif
          vel(3,i) = vel(3,i) + domegadr*Omega0*(xmax(1)-xmin(1))
          velin(3,i) = velin(3,i) + domegadr*Omega0*(xmax(1)-xmin(1))
          ncross = ncross + 1
       elseif(x(1,i).lt.xmin(1)) then
          x(1,i) = x(1,i) + (xmax(1) - xmin(1))
          x(2,i) = x(2,i) !- domegadr*Omega0*(xmax(1)-xmin(1))*time
          x(2,i) = xmodbound(x(2,i),xmin(2),xmax(2),0.)
!          if (x(2,i).lt.xmin(2)) then
!             x(2,i) = xmax(2) - (xmin(2) - x(2,i))
!          endif
          vel(3,i) = vel(3,i) - domegadr*Omega0*(xmax(1)-xmin(1))
          velin(3,i) = velin(3,i) - domegadr*Omega0*(xmax(1)-xmin(1))
          ncross = ncross + 1
       endif
    enddo
    if (ncross.gt.0) print*,ncross,' particles crossing x (shearing box)'
    return
 endif

!-------------------------------------------------------------------------
! other ghost boundaries - correct particles which have crossed boundary
!------------------------------------------------------------------------- 
 if (any(ibound.gt.1)) then
    ncorrect = 0
    do i=1,npart
       do jdim=1,ndim
          if (ibound(jdim).gt.1) then
             if (x(jdim,i).gt.xmax(jdim)) then
                !--move particle back inside boundary
                x(jdim,i) = xmax(jdim) - (x(jdim,i) - xmax(jdim))
                xin(jdim,i) = x(jdim,i)
                !--flip normal velocity component
                vel(jdim,i) = -vel(jdim,i)
                velin(jdim,i) = vel(jdim,i)
                force(jdim,i) = 0. !-force(jdim,i)
                ncorrect = ncorrect + 1
!                print*,'xmax',jdim,i,'xnew = ',x(jdim,i),vel(jdim,i)
             elseif(x(jdim,i).lt.xmin(jdim)) then
                print*,'xminold',jdim,i,'x,v = ', x(jdim,i),vel(jdim,i)
                !--move particle back inside boundary
                x(jdim,i) = xmin(jdim) + (xmin(jdim) - x(jdim,i))
                xin(jdim,i) = x(jdim,i)
                !--flip normal velocity component
                vel(jdim,i) = -vel(jdim,i)
                velin(jdim,i) = vel(jdim,i)
                force(jdim,i) = 0. !-force(jdim,i)
                ncorrect = ncorrect + 1
                print*,'xminnew',jdim,i,'x,v = ',x(jdim,i),vel(jdim,i)
             endif
          endif
       enddo
    enddo
    if (ncorrect.gt.0) write(iprint,*) ' warning: ',ncorrect,' particles '// &
                                       ' adjusted from boundary'
 endif

 return 
end subroutine boundary
