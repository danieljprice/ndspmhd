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
!! this subroutine returns the prsure and sound speed
!! from rho and/or u_therm via the equation of state
!!
!! size of array is specified on input, so can send in either whole
!! arrays or individual elements (just be consistent!)
!!-----------------------------------------------------------------
!-------------------------------------------------------------------
!  equation of state related quantities
!-------------------------------------------------------------------

module eos
 implicit none
 real :: gamma, polyk

contains

subroutine equation_of_state(pr,vsound,uu,rho,gammai)
 use options, only:iener
 use loguns
 use part, only:itype,itypegas,itypegas1,itypegas2,itypebnd
!
!--define local variables
!
 implicit none
 integer :: i,isize
 real, intent(in), dimension(:) :: rho
 real, intent(out), dimension(size(rho)) :: pr
 real, intent(inout), dimension(size(rho)) :: uu,vsound
 real, intent(in), dimension(size(rho)), optional :: gammai
 real :: gamma1
 
 isize = size(rho)
 gamma1 = gamma - 1.
!
!--exit gracefully if rho is negative
!
 if (any(rho.lt.0.)) then
    write(iprint,*) 'eos: rho -ve, exiting'
    do i=1,isize
       if (rho(i).lt.0.) write(iprint,*) i,rho(i),uu(i)
    enddo
    !call quit
 elseif ((iener.ne.0).and.any(uu.lt.0.)) then
    write(iprint,*) 'eos: u_therm -ve, exiting',isize
    do i=1,isize
       if (uu(i).lt.0.) write(iprint,*) i,rho(i),uu(i)
    enddo    
    stop
    !call quit
 endif

 if (iener.eq.0) then   ! polytropic (isothermal when gamma=1)
    where (rho > 0. .AND. itype(1:isize).EQ.itypegas)
       pr = polyk*rho**gamma
       vsound = sqrt(gamma*pr/rho)
    elsewhere (itype(1:isize).eq.itypegas1)
       pr = polyk*(rho - 1.)
    elsewhere (itype(1:isize).eq.itypegas2)
       pr =  polyk*(rho - 1.) !4.*polyk*((rho/0.5) - 1.)
    elsewhere(itype(1:isize).ne.itypebnd)
      pr = 0.
    end where
    if (abs(gamma1).gt.1.e-3) then       
       where (rho > 0.) 
       uu = pr/(gamma1*rho)
       end where
    endif   
 else      ! adiabatic
    if (present(gammai)) then
       where (rho > 0.)
         pr = (gammai-1)*uu*rho
         vsound = sqrt(gammai*pr/rho)
       end where    
    else
       where (rho > 0.)
         pr = gamma1*uu*rho
         vsound = sqrt(gamma*pr/rho)
       end where
    endif
 endif
      
 return
end subroutine equation_of_state

subroutine equation_of_state1(pr,vsound,uu,rho,gammai)
 use options, only:iener
 use loguns
!
!--define local variables
!
 implicit none
 real, intent(in) :: rho
 real, intent(out) :: pr
 real, intent(inout) :: uu,vsound
 real, intent(in), optional :: gammai
 real :: gamma1
 
 gamma1 = gamma - 1.
!
!--exit gracefully if rho is negative
!
 if (rho.lt.0.) then
    write(iprint,*) 'eos1: rho -ve, exiting'
    !call quit
 elseif ((iener.ne.0).and.uu.lt.0.) then
    write(iprint,*) 'eos1: u_therm -ve, exiting',1    
    !call quit
 endif

 if (iener.eq.0) then   ! polytropic (isothermal when gamma=1)
    if (rho > 0.) then
      pr = polyk*rho**gamma
      vsound = sqrt(gamma*pr/rho)
    endif
    if (abs(gamma1).gt.1.e-3) then       
       if (rho > 0.) then
          uu = pr/(gamma1*rho)
       endif
    endif   
 else      ! adiabatic
    if (present(gammai)) then
       if (rho > 0.) then
         pr = (gammai-1.)*uu*rho
         vsound = sqrt(gammai*pr/rho)
       endif
    else
       if (rho > 0.) then
         pr = gamma1*uu*rho
         vsound = sqrt(gamma*pr/rho)
       endif
    endif
    !print *,'here ',uu,rho
 endif
      
 return
end subroutine equation_of_state1

end module eos
