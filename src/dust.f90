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

module dust
 implicit none
 real, private :: epstein_fac
 real, public :: grain_size, grain_mass, grain_dens
 logical :: done_init_drag = .false.

 real :: grain_size_cm = 0.1    ! grain size in cm (default=mm)
 real, parameter :: grain_dens_cgs = 3. ! g/cm^3
! constants to translate to physical units
 real, parameter :: mm = 0.1       ! 1 mm in cm
 real, parameter :: micron = 1.e-4 ! 1 micron in cm
 real, parameter :: au = 1.496d13
 real, parameter :: solarm = 1.9891d33
 real, parameter :: gg = 6.672d-8
 real, parameter :: years = 3.1556926d7
 real, parameter :: umass = solarm
 real, parameter :: udist = 10.*au
 real, parameter :: utime = sqrt(udist**3/(gg*umass))

contains

!------------------------------------------------------------------------
! Subroutine to pre-calculate various quantities required to get the
! stopping time
!------------------------------------------------------------------------
subroutine init_drag(ierr,gamma)
 use options,      only:Kdrag
 use setup_params, only:pi
 integer, intent(out) :: ierr
 real, intent(in)     :: gamma

 ierr = 0
 print*,' initialising drag terms'
 if (gamma < 1.) then
    print*,' error: gamma < 1 in drag initialisation'
    ierr = 1
    return
 endif
 grain_size_cm = Kdrag

 grain_dens = grain_dens_cgs/(umass/udist**3)
 grain_size = grain_size_cm/udist
 print*,' grain size    =',grain_size,' in cm     =',grain_size*udist,' in m =',grain_size*udist*0.01
 print*,' grain density =',grain_dens,' in g/cm^3 =',grain_dens*(umass/udist**3)

 epstein_fac = 4./3.*sqrt(8.*pi/gamma)
 grain_mass  = 4./3.*pi*grain_dens*grain_size**3
 done_init_drag = .true.

end subroutine init_drag

!------------------------------------------------------------------------
! Function to return the stopping time for various dust prescriptions
!------------------------------------------------------------------------
real function get_tstop(idrag_nature,rhogas,rhodust,cs,Kdrag)
 integer, intent(in) :: idrag_nature
 real, intent(in)    :: rhogas,rhodust,cs,Kdrag
 real :: rho,ts,ts1
 
 if (.not.done_init_drag) stop 'drag not initialised before call to get_tstop'
 rho = rhogas + rhodust
 select case(idrag_nature)
  case(1) !--constant drag
     ts = rhodust*rhogas/(Kdrag*rho)
  case(2,4) !--constant ts
     ts = Kdrag
  case(3) !--Epstein regime
     if (grain_mass <= 0.) stop 'grain mass < 0; Epstein drag not initialised'
     if (epstein_fac <= 0.) stop 'fac < 0; Epstein drag not initialised'
     ts1 = epstein_fac*grain_size**2*cs*rho/grain_mass
     ts  = 1./ts1
     !print*,' DEBUG: equivalent Kdrag = ',rhogas*rhodust/(ts*rho)
  case default
     ts = huge(ts)
 end select
 get_tstop = ts

 return
end function get_tstop
      
end module dust
