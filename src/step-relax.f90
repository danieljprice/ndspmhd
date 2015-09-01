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

!----------------------
! Relax the particles
!----------------------
subroutine step
 use dimen_mhd
 use debug
 use loguns 
 use bound
 use options
 use part
 use part_in
 use rates
 use eos
 use timestep, only:C_cour,dtcourant,dt
!
!--define local variables
!
 implicit none
 integer :: i
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' Entering subroutine relax'

 do i=1,npart
    xin(:,i) = x(:,i)
 enddo
 polyk = 1.
 gamma = 2.
 iener = 0
!
!--shift positions
!
 do i=1,npart
    if (itype(i).eq.itypebnd .or. itype(i).eq.itypebnd2) then ! fixed particles
       x(:,i) = xin(:,i)
    else
       x(:,i) = xin(:,i) + 0.5*dt**2*force(1:ndim,i)
    endif
    vel(:,i) = 0.
 enddo 
!
!--calculate all derivatives
!
 call derivs

 dt = C_cour*dtcourant

 if (any(ibound.ne.0)) call boundary ! inflow/outflow/periodic boundary conditions

 if (trace) write (iprint,*) ' Exiting subroutine relax'
      
 return
end subroutine step
