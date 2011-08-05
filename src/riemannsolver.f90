!-------------------------------------------------------------------
! Implementation of the iterative Riemann solver of Van Leer (1997)
! Follows Cha & Whitworth (2003), MNRAS 340, 73-90
!-------------------------------------------------------------------
subroutine riemannsolver(gamma,pra,prb,va,vb,ca,cb,pr,vstar)
  implicit none
  real, parameter :: tol = 1.5e-2
  real :: gamma,pra,prb,va,vb,ca,cb,pr,vstar
  real :: wa,wb,za,zb,vstara,vstarb,prnew
!
!--get an initial starting estimate of intermediate pressure
!
  prnew = (ca*prb + cb*pra - ca*cb*(va-vb))/(ca+cb)
  pr = pra

  do while (abs(prnew-pr).lt.tol)

     pr = prnew
!
!--find Lagrangian shock speed and tangential slope
!
     wa = wspeed(gamma,pr,pra,ca)
     wb = wspeed(gamma,pr,prb,cb)
     za = zslope(wa,gamma,pr,pra,ca)
     zb = zslope(wb,gamma,pr,prb,cb)
!
!--update vstar
!
     vstara = va + (pr-pra)/wa
     vstarb = vb - (pr-prb)/wb
!
!--then get new estimate of pr
!
     prnew = pr - zb*za*(vstara - vstarb)/(za+zb)

  enddo

  print*,'its = ',its,' pr = ',prnew,' pra,prb = ',pra,prb
  pr = prnew
  vstar = (zb*vstarb + za*vstara)/(za+zb)

end subroutine riemannsolver

!
!--Lagrangian shock speed
!
real function wspeed(gamma,pr,pra,ca)
  implicit none
  real :: gamma, pr, pra, ca, gamfac
  
  if (pr.ge.pra) then
     wspeed = ca*sqrt(1.+(gamma+1.)/(2.*gamma)*(pr-pra)/pra)
  else
     gamfac = gamma-1./(2.*gamma)
     wspeed = ca*gamfac*(1.-pr/pra)/(1.-(pr/pra)**gamfac)
  endif
  
end function wspeed
!
!--tangential slope
!  
real function zslope(wa,gamma,pr,pra,ca)
  implicit none
  real :: gamma,pr,pra,ca,wa,w2
  
  if (pr.ge.pra) then
     w2 = wa*wa
     zslope = 2.*w2/(w2 + ca**2)*wa
  else
     zslope = ca*(pr/pra)**(1-(gamma-1)/(2.*gamma))
  endif
end function zslope
