!
!--copies all the particle properties from particle j to particle i
!  just so you don't have to remember everything carried by the particles
!  *does not copy particle position or velocity*
!
subroutine copy_particle(i,j)
  use debug
  use loguns
  use derivB
  use fmagarray
  use hterms
  use part
  use rates
  use unityfunc
  use xsph
  implicit none
  integer :: i,j

  if (idebug(1:4).eq.'copy') write(iprint,*) 'copying particle ',j,' to particle ',i
  vel(:,i) = vel(:,j)
  pmass(i) = pmass(j)
  rho(i) = rho(j)
  hh(i) = hh(j)
  uu(i) = uu(j)
  en(i) = en(j)
  Bcons(:,i) = Bcons(:,j)
  Bfield(:,i) = Bfield(:,j)
  alpha(:,i) = alpha(:,j)
  psi(i) = psi(j)
  gradh(i) = gradh(j)
  gradhaniso(i) = gradhaniso(j)
  sqrtg(i) = sqrtg(j)

!!  unity(i) = unity(j)
!!  gradunity(:,i) = gradunity(:,j)

  force(:,i) = force(:,j)
  drhodt(i) = drhodt(j)
  dudt(i) = dudt(j)
  dendt(i) = dendt(j)  
  dBconsdt(:,i) = dBconsdt(:,j)
  dhdt(i) = dhdt(j)
  daldt(:,i) = daldt(:,j)
  dpsidt(i) = dpsidt(j)
  xsphterm(:,i) = xsphterm(:,j) ! after here not crucial
  fmag(:,i) = fmag(:,j)
  divB(i) = divB(j)
  curlB(:,i) = curlB(:,j) 
 
end subroutine copy_particle
