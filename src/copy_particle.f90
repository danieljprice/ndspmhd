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
  use options, only:imhd,onef_dust
  !use matrixcorr
  implicit none
  integer :: i,j

  if (idebug(1:4).eq.'copy') write(iprint,*) 'copying particle ',j,' to particle ',i
!  vel(:,i) = vel(:,j)
  pmass(i) = pmass(j)
  rho(i)    = rho(j)
  rhoalt(i) = rhoalt(j)
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
  if (onef_dust) then
     dustfrac(i)    = dustfrac(j)
     dustevol(i)    = dustevol(j)
     ddustevoldt(i) = ddustevoldt(j)
     deltav(:,i)     = deltav(:,j)
     ddeltavdt(:,i)  = ddeltavdt(:,j)
     rhodust(i)     = rhodust(j)
     rhogas(i)      = rhogas(j)
  endif
!
! quantum
!
  if (allocated(P_Q)) then
     P_Q(:,:,i) = P_Q(:,:,j)
  endif
!
! viscosity
!
  if (allocated(del2v)) del2v(:,i) = del2v(:,j)
  if (allocated(graddivv)) graddivv(:,i) = graddivv(:,j)

  if (allocated(rho0)) rho0(i) = rho0(j)
  !gradmatrix(:,:,i) = gradmatrix(:,:,j)
!
!  Note to self: don't forget, that if you *really*
!  need to reassign particles, then you may also need
!  to copy the part_in arrays used in the timestepping
!  as well (I've spent a few days rediscovering this
!  several times)
! 

end subroutine copy_particle
