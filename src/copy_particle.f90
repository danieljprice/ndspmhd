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

!
!--copies all the particle properties from particle j to particle i
!  just so you don't have to remember everything carried by the particles
!  *does not copy particle position or velocity*
!
subroutine copy_particle(i,j)
  use debug
  use loguns
  use eos
  use derivB
  use fmagarray
  use hterms
  use part
  use rates
  use xsph
  use options, only:imhd,idust
  !use matrixcorr
  implicit none
  integer :: i,j

  if (idebug(1:4).eq.'copy') write(iprint,*) 'copying particle ',j,' to particle ',i
!  vel(:,i) = vel(:,j)
  pmass(i) = pmass(j)
  rho(i) = rho(j)
  hh(i) = hh(j)
  uu(i) = uu(j)
  en(i) = en(j)
  if (.not.(imhd.lt.0 .and. itype(i).gt.0)) Bevol(:,i) = Bevol(:,j)
  Bfield(:,i) = Bfield(:,j)
  alpha(:,i) = alpha(:,j)
  psi(i) = psi(j)
  gradh(i) = gradh(j)
  gradhn(i) = gradhn(j)
  gradsoft(i) = gradsoft(j)
  gradgradh(i) = gradgradh(j)
  if (allocated(zeta)) zeta(i) = zeta(j)
  sqrtg(i) = sqrtg(j)
  spsound(i) = spsound(j)
  pr(i) = pr(j)
  dens(i) = dens(j)
  if (itype(i).ne.itypebnd .and. itype(i).ne.itypebnd2 .and. itype(j).eq.itypegas) itype(i) = itype(j)
  if (allocated(pmom)) pmom(:,i) = pmom(:,j)

  force(:,i) = force(:,j)
  drhodt(i) = drhodt(j)
  dudt(i) = dudt(j)
  dendt(i) = dendt(j)  
  dBevoldt(:,i) = dBevoldt(:,j)
  gradpsi(:,i) = gradpsi(:,j)
  dhdt(i) = dhdt(j)
  daldt(:,i) = daldt(:,j)
  dpsidt(i) = dpsidt(j)
  xsphterm(:,i) = xsphterm(:,j) ! after here not crucial
  fmag(:,i) = fmag(:,j)
  divB(i) = divB(j)
  curlB(:,i) = curlB(:,j)
!
!  dust
!
  if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
     dustfrac(i)    = dustfrac(j)
     dustevol(i)    = dustevol(j)
     ddustevoldt(i) = ddustevoldt(j)
     deltav(:,i)     = deltav(:,j)
     ddeltavdt(:,i)  = ddeltavdt(:,j)
  endif

end subroutine copy_particle
