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
!! Sets the logical unit numbers to use for output
!!-----------------------------------------------------------------

SUBROUTINE logun
 USE dimen_mhd
 USE loguns
 IMPLICIT NONE
      
! IF (ndim.EQ.1) THEN
    iprint = 6
! ELSE
!    iprint = 8  ! log file / screen
! ENDIF
 ievfile = 13
 idatfile = 12
 iread = 21  ! for reading from input file
 ireadf = 15 ! for reading filenames
      
 RETURN
END SUBROUTINE logun
