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
  use options, only:imhd
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
  
  !gradmatrix(:,:,i) = gradmatrix(:,:,j)
 
end subroutine copy_particle
