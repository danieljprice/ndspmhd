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

!-------------------------------------------------------------------
! Implementation of the exact Riemann solver given in Toro (1992)
!
! Solves for the post-shock pressure (pr) and velocity (vstar)
! given the initial left and right states
!
! Does not matter if high P / high rho is on left or right
!
! Daniel Price, Institute of Astronomy, Cambridge, UK, 2004
! dprice@ast.cam.ac.uk
!-------------------------------------------------------------------
subroutine riemannsolver(gamma,p_L,p_R,v_L,v_R,rho_L,rho_R,pr,vstar)
  implicit none
  real, parameter :: tol = 1.5e-2
  real, intent(in) :: gamma,p_L,p_R,v_L,v_R,rho_L,rho_R
  real, intent(out) :: pr,vstar
  integer, parameter :: maxits = 30
  integer :: its
  real :: c_L, c_R
  real :: prnew, f_L, f_R, dfdp_L, dfdp_R, f, df, dp
  real :: power, denom, cs2
!
!--use isothermal solver if appropriate
!
  c_L = sqrt(gamma*p_L/rho_L)
  c_R = sqrt(gamma*p_R/rho_R)
  if (gamma.lt.1.0001) then
     cs2 = p_L/rho_L
     call get_pstar_isothermal(cs2,v_L,v_R,rho_L,rho_R,pr,vstar)
     return
  endif
!
!--get an initial starting estimate of intermediate pressure
!  this one is from Toro(1992) - gives basically the right answer
!  for pressure jumps below about 4
!
  power = (gamma-1.)/(2.*gamma)
  denom = c_L/p_L**power + c_R/p_R**power
  prnew = ((c_L + c_R + (v_L - v_R)*0.5*(gamma-1.))/denom)**(1./power)
  pr = p_L
  its = 0

  !!print*,'initial guess = ',prnew

  do while (abs(prnew-pr).gt.tol .and. its.lt.maxits)

     its = its + 1
     pr = prnew
!
!--evaluate the function and its derivatives
!
     call f_and_df(pr,p_L,c_L,gamma,f_L,dfdp_L)
     call f_and_df(pr,p_R,c_R,gamma,f_R,dfdp_R)
!
!--then get new estimate of pr
!
     f = f_L + f_R + (v_R - v_L)
     df = dfdp_L + dfdp_R
!
!--Newton-Raphson iterations
!     
     dp = -f/df
     prnew = pr + dp
     
  enddo

  if (its.eq.maxits) print*,'WARNING: its not converged in riemann solver'
  if (prnew.le.0.) then
     print*,'ERROR: pr < 0 in riemann solver'
     print*,'its = ',its,'p_L, p_R = ',p_L,p_R,' v_L, v_R = ',v_L,v_R,' p* = ',prnew,'v = ',vstar,v_R + f_R
  endif
  pr = prnew
  vstar = v_L - f_L
  
!  if (its.gt.0) then
!     print*,'its = ',its,'p_L, p_R = ',p_L,p_R,' v_L, v_R = ',v_L,v_R,' p* = ',prnew,'v = ',vstar,v_R + f_R
!  endif
  
end subroutine riemannsolver

!
!--pressure function
!  H is pstar/p_L or pstar/p_R
!
subroutine f_and_df(prstar,pr,cs,gam,fp,dfdp)
  implicit none
  real, intent(in) :: prstar, pr, gam, cs
  real, intent(out) :: fp,dfdp
  real :: H,term, power, gamm1, denom

  H = prstar/pr
  gamm1 = gam - 1.
  
  if (H.gt.1.) then  ! shock
     denom = gam*((gam+1.)*H + gamm1)
     term = sqrt(2./denom)
     fp = (H - 1.)*cs*term
         
     dfdp = cs*term/pr + (H - 1.)*cs/term*(-1./denom**2)*gam*(gam+1.)/pr
  else               ! rarefaction
     power = gamm1/(2.*gam)
     fp = (H**power - 1.)*(2.*cs/gamm1)
  
     dfdp = 2.*cs/gamm1*power*H**(power-1.)/pr
  endif
  
end subroutine f_and_df

!-------------------------------------------------------------
! Non-iterative isothermal Riemann solver 
! from Balsara (1994), ApJ 420, 197-212
!
! See also Cha & Whitworth (2003), MNRAS 340, 73-90
!-------------------------------------------------------------
subroutine get_pstar_isothermal(cs2,v_L,v_R,rho_L,rho_R,pstar,vstar)
  implicit none
  real, intent(in) :: cs2,v_L,v_R,rho_L,rho_R
  real, intent(out) :: pstar,vstar
  real :: sqrtrho_L, sqrtrho_R, X, vdiff, determinant, vstar2

  sqrtrho_L = sqrt(rho_L)
  sqrtrho_R = sqrt(rho_R)
  
  X = sqrtrho_L*sqrtrho_R/(sqrtrho_L + sqrtrho_R)
  vdiff = v_L - v_R
  determinant = (X*vdiff)**2 + 4.*cs2*X*(sqrtrho_L + sqrtrho_R)
  
  pstar = 0.25*(X*vdiff + sqrt(determinant))**2  
  vstar = v_L - (pstar - cs2*rho_L)/(sqrt(pstar*rho_L))
  vstar2 = v_R + (pstar - cs2*rho_R)/(sqrt(pstar*rho_R))
  if (abs(vstar2-vstar).gt.1.e-5) print*,'error: vstar = ',vstar,vstar2
  !print*,' pstar = ',pstar,' v_R,v_L = ',v_L,v_R,cs2,' vstar = ',vstar
  
end subroutine get_pstar_isothermal
