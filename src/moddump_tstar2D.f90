!----------------------------------------------------
! this subroutine modifies the static configuration
! (from the dumpfile) and gives it the appropriate
! velocity perturbation for the toy star
!
! NB: modification is done in cartesian coords
!----------------------------------------------------
subroutine modify_dump
 use dimen_mhd
 use debug
 use loguns
 use eos, only:gamma,polyk
 use part
 use geometry
 implicit none
 integer :: i,ierr,jmode,smode
 real, parameter :: pi=3.1415926536
 real :: Ctstar,Atstar,scalefac,sigma2,sigma,rstar,denscentre,gamm1
 real :: omegasq,cs2centre,ekin,ekin_norm
 real, dimension(ndim) :: xcyl,velcyl,dvel
 character(len=len(rootname)+8) :: tstarfile
 character(len=30) :: dummy

 jmode = 2
 smode = 0
 
 write(iprint,*) 'MODIFYING INITIAL SETUP with toystar oscillations'

!
!--read parameters from file
!  
 tstarfile = rootname(1:len_trim(rootname))//'.tstar2D'
 open(unit=ireadf,err=11,file=tstarfile,status='old',form='formatted')
    read(ireadf,*,err=12) dummy
    read(ireadf,*,err=12) dummy
    read(ireadf,*,err=12) jmode,smode
 close(unit=ireadf)
   goto 13
11 continue
   write(iprint,*) tstarfile,' not found, using default options '
   goto 13
12 continue
   write(iprint,*) ' error reading ',tstarfile
13 continue   

 write(iprint,*) 'radial mode = ',jmode,' theta mode = ',smode
 
 gamm1 = gamma - 1.
 if (gamm1.lt.1.e-5) then
    stop 'error: gamma - 1 <= 0'
 endif

 omegasq = 1.0
 sigma2 = 0.5*omegasq*(gamm1)*((jmode+smode)*(jmode+smode + 2./gamm1) - smode**2)
 if (sigma2.lt.1.e-5) then
    print*,'ERROR sigma2 < 0 in perturbation'
 else
    sigma = sqrt(sigma2)
 endif

 denscentre = 1.0
 Ctstar = 1.0
 scalefac = polyk*gamma/(sigma*gamm1)
 rstar = sqrt((2.*polyk*gamma*denscentre**gamm1)/gamm1)
 
 write(iprint,*) 'polyk = ',polyk,' rstar = ',rstar,' period = ',2.*pi/sigma
 cs2centre = gamma*polyk*denscentre**gamm1
 write(iprint,*) 'denscentre = ',denscentre,' cs_0 = ',sqrt(cs2centre)

!
!--work out frequency of oscillation
!
 ekin = 0.
 ekin_norm = 0.
 do i=1,npart
    !--get r,theta
    call coord_transform(x(:,i),ndim,1,xcyl(:),ndim,2)
    
    !--set v_r
    velcyl(1) = scalefac*detadr(jmode,smode,xcyl(1)/rstar,gamma)*COS(smode*xcyl(2))
    !--set theta_dot
    velcyl(2) = -scalefac*etar(jmode,smode,xcyl(1)/rstar,gamma)*smode*SIN(smode*xcyl(2))/xcyl(1)**2
    !!print*,'v_phi = ',velcyl(2),xcyl(2),etar(jmode,smode,xcyl(1)/rstar,gamma)
    !--now transform back to get vx, vy
    call vector_transform(xcyl(1:ndim),velcyl(1:ndim),ndim,2,dvel(1:ndim),ndim,1)
    if (xcyl(1).lt.1.e-5) then
       print*,' r = 0 on particle ',i,' xcyl(1) = ',xcyl(1), &
              ' v_cyl = ',velcyl,' v_cart = ',dvel
    endif
    !--now perturb v with appropriate amplitude
    vel(1:ndim,i) = dvel(1:ndim)
    ekin = ekin + 0.5*pmass(i)*dot_product(vel(1:ndim,i),vel(1:ndim,i))
    ekin_norm = ekin_norm + 0.5*pmass(i)
 enddo
!
!--normalise the amplitude
!
 ekin_norm = (0.05)**2*cs2centre*ekin_norm
 write(iprint,*) ' ekin = ',ekin, ' ekin_norm = ',ekin_norm
 vel = vel*sqrt(ekin_norm/ekin)
 
 Atstar = scalefac*sqrt(ekin_norm/ekin)
 write(iprint,*) ' v = ',Atstar,'*detadr(r)'
!
!--rewrite the tstar2D file giving the amplitude
!  
 tstarfile = rootname(1:len_trim(rootname))//'.tstar2D'
 write(iprint,*) ' writing to file ',trim(tstarfile)
 open(unit=ireadf,iostat=ierr,file=tstarfile,status='replace',form='formatted')
 if (ierr.eq.0) then
    write(ireadf,*,iostat=ierr) denscentre,Ctstar,Atstar
    write(ireadf,*,iostat=ierr) 0.
    write(ireadf,*,iostat=ierr) jmode,smode
    if (ierr /= 0) write(iprint,*) 'ERROR WRITING TO ',trim(tstarfile)
    close(unit=ireadf)
 else
    write(iprint,*) 'ERROR OPENING ',trim(tstarfile)
 endif

 return

contains

!
!--function that evaluates the polynomial for rho(r/re) for a given radial mode
!  (from the power series solution to the 2nd order ODE)
!
!  rad = r/r_star
!  j = radial (axisymmetric) mode
!  m = theta mode 
!
!  solution is for delta(rho**(gamma-1))
!  ie. rho**(gamma-1) = rho_0**(gamma-1) + etar
!
!  and takes the form
!
!  etar = rad**m sum_k a_k rad**k
!
real function etar(j,m,rad,gamma)
  implicit none 
  integer :: j,m,k,kprev   ! j is the radial mode, m is the theta mode
  real :: rad,gamma,denom
  real :: ak,akprev,gamm1,freqsq
!
!--this solution is for arbitrary gamma
!
  gamm1 = gamma - 1.
  if (gamm1.lt.1.e-3) then
     print*,'error gamma -1 <= 0'
     etar = 0.
     return
  endif
!
!--the solution is of the form
!  drhor = a_0 + a_2 (r/re)**2 + a_4 (r/re)**4 + ...
!  where for j = k, coefficients >= a_k+2 are zero
!  
  freqsq = (j+m)*(j+m + 2./gamm1) - m**2

  akprev = 1.0  ! this is a_0 which is the amplitude
  etar = akprev
  !!print*,'mode = ',j,m,' nu^2 = ',freqsq,' a_0 = ',akprev
!
!--the co-efficients for the terms above a_0 are calculated using
!  the recurrence relation between the a_k's
!
  do k = 2,j,2
     kprev = k-2
     denom = real((kprev + 2 + m)**2 - m**2)
     ak = akprev*(kprev**2 + 2.*kprev*m + 2.*(kprev+m)/gamm1 - freqsq)/denom
     !!print*,'coeff ',k,' = ',ak,k**2,2.*k/gamm1
     etar = etar + ak*rad**k
     akprev = ak
  enddo
  
  etar = etar * rad**m

end function etar

!
!--function that evaluates the polynomial for v(r/re) for a given radial mode
!  (from the power series solution to the 2nd order ODE)
!
real function detadr(j,m,rad,gamma)
  implicit none
  integer :: j,m,k,kprev   ! j is the radial mode, m is the theta mode
  real :: rad,gamma,denom,term1,term2
  real :: ak,akprev,gamm1,freqsq
!
!--this solution is for arbitrary gamma
!
  gamm1 = gamma - 1.
  if (gamm1.lt.1.e-3) then
     print*,'error gamma -1 <= 0'
     detadr = 0.
     return
  endif
!
!--the solution is of the form
!  drhor = a_0 + a_2 (r/re)**2 + a_4 (r/re)**4 + ...
!  where for j = k, coefficients >= a_k+2 are zero
!  
  freqsq = (j+m)*(j+m + 2./gamm1) - m**2

  detadr = 0.
  akprev = 1.0  ! this is a_0 which is the amplitude
  term1 = akprev
  term2 = 0.
!  print*,'mode = ',j,m,' nu^2 = ',freqsq,' a_0 = ',akprev
!
!--the co-efficients for the terms above a_0 are calculated using
!  the recurrence relation between the a_k's
!
  do k = 2,j,2
     kprev = k-2
     denom = real((kprev + 2 + m)**2 - m**2)
     ak = akprev*(kprev**2 + 2.*kprev*m + 2.*(kprev+m)/gamm1 - freqsq)/denom
     !!print*,'coeff ',k,' = ',ak,k*ak,rad,(k-1)
     term1 = term1 + ak*rad**k
     term2 = term2 + k*ak*rad**(k-1)
     akprev = ak
  enddo
  
  if (m.eq.0) then
     detadr = term2
  else
     detadr = m*rad**(m-1)*term1 + rad**m*term2
  endif
  
end function detadr

end subroutine modify_dump
