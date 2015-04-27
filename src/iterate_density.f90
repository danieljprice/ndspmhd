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

!!------------------------------------------------------------------------
!! Iterates the density calculation so that rho and h are
!! calculated self-consistently. 
!!
!! rho_a is a function of h_a and vice versa
!! via the summation rho_a = sum_b m_b W_ab(h_a)
!! we assume h \propto 1/rho^(1/ndim)
!!
!! we determine the change in h and recalculate the density only 
!! for those particles where h changes significantly after the density
!! is calculated.
!!
!! For details see Price and Monaghan 2004, MNRAS 348, 139
!!
!! Number density option has been added Aug 2005.
!! IMPORTANT: if number density formalism is used, must change the
!!            relevant lines in rates also.
!!
!!------------------------------------------------------------------------

subroutine iterate_density
  use dimen_mhd
  use debug
  use loguns
  
  use bound
  use hterms
  use linklist, only:numneigh
  use options, only:ikernav,ihvar,ibound,maxdensits,tolh,usenumdens
  use part,    only:npart,ntotal,itype,x,pmass,hh,vel,rho,rhoalt,itypebnd
  use setup_params
  use density_summations
  use rates, only:dhdt,drhodt
  !!use linklist, only:numneigh
!
!--define local variables
!
  implicit none
  integer :: i,j,itsdensitymax
  integer :: ncalc,ncalcprev,ncalctotal,nrhosmall,nwarn
  integer, dimension(npart) :: redolist, redolistprev
  real :: hnew,func,dfdh
  real :: rhoi,dhdrhoi,omegai,densnumi,dhdni,dwdhsumi,d2hdrho2i
  real, dimension(2*size(rho)) :: hhin,dndt ! allow space in case of ghost reallocation
  logical :: converged,redolink
  
!!  integer :: itest
!
!--allow for tracing flow
!      
  if (trace) write(iprint,*) ' Entering subroutine iterate_density' 
!
!--set maximum number of iterations to perform
! 
  if ((ikernav.eq.3).and.(ihvar.ne.0)) then
     itsdensitymax = maxdensits  ! perform 1 fixed point iteration
  else
     itsdensitymax = 0   ! no iterations
  endif
!
!--Loop to find rho and h self-consistently (if using Springel/Hernquist)
!
  itsdensity = 0
  ncalctotal = 0
  ncalc = npart   ! number of particles to calculate density on
  redolink = .false.
  ncalcprev = 0
  where (itype /= itypebnd)
     gradh = 0.
     gradhn = 0.
     gradsoft = 0.
     gradgradh = 0.
     drhodt = 0.
     dhdt = 0.
  end where
  hhin(1:npart) = hh(1:npart)
  if (any(hh(1:npart).le.tiny(hh))) then
     write(iprint,*) 'error: h <= 0 in density call'
     call quit
  endif
  dhdni = 0.
  dhdrhoi = 0.

!!  itest = 416
  
  if (ncalc.eq.npart) then
     do j=1,npart     
        redolist(j) = j
     enddo
  endif
  if (any(pmass(1:npart).ne.pmass(1))) then
     rhomin = 0. !!minval(rho(1:npart))
  else
     rhomin = 0.
  endif
  
  iterate: do while ((ncalc.gt.0).and.(itsdensity.le.itsdensitymax))     
     itsdensity = itsdensity + 1

     if (redolink) then
        if (any(ibound.gt.1)) call set_ghost_particles
        if (idebug(1:3).eq.'den') write(iprint,*) 'relinking...'
        call set_linklist
     endif
!
!--calculate the density (using h), for either all the particles or
!  only on a partial list
!     
     if (ncalc.eq.npart) then
        call density(x,pmass,hh,vel,rho,drhodt,rhoalt,dndt,gradh,gradhn,gradsoft,gradgradh,npart,ntotal) ! symmetric for particle pairs
     else
        call density_partial(x,pmass,hh,vel,rho,drhodt,rhoalt,dndt,gradh,gradhn,gradsoft,gradgradh,ntotal,ncalc,redolist)
     endif
!
!--rhoalt is the number density times the particle mass
!  (also computed with ikernelalt instead of ikernel)
!
     ncalctotal = ncalctotal + ncalc
     ncalcprev = ncalc
     redolistprev(1:ncalcprev) = redolist(1:ncalcprev)
     ncalc = 0
     redolink = .false.
     nrhosmall = 0
     nwarn = 0
     
     if (ikernav.eq.3) then
        do j=1,ncalcprev
           i = redolistprev(j)

           if (itype(i).ne.itypebnd) then
              if (rho(i).le.1.e-6) then
                 if (rho(i).le.0.) then
                    write(iprint,*) 'error: rho(',i,') = ',rho(i),hh(i),pmass(i)
                    call quit
                 else
                    nrhosmall = nrhosmall + 1
                    !!write(iprint,*) 'Warning : rho < 1e-6 '
                 endif
              endif
              
              if (usenumdens) then
                 densnumi = (hfact/hh(i))**ndim                ! this is the number density compatible with the old h
!                 dhdni = -hh(i)/(ndim*(densnumi + rhomin))          ! deriv of this
                 dhdni = -hh(i)/(ndim*(rhoalt(i) + rhomin))          ! deriv of this
                 omegai =  1. - dhdni*gradhn(i)           
                 if (omegai.lt.1.e-5) then
                    !print*,'warning: omega < 1.e-5 ',i,omegai
                    if (abs(omegai).eq.0.) omegai = 1. ! call quit
                 endif
                 gradhn(i) = 1./omegai
                 gradh(i) = gradh(i)*dhdni
                 func = densnumi - rhoalt(i)
                 dfdh = omegai/dhdni
                 gradsoft(i) = gradsoft(i)*dhdni
              else
                 rhoi = pmass(i)/((hh(i) - h_min)/hfact)**ndim - rhomin ! this is the rho compatible with the old h
                 dhdrhoi = -(hh(i) - h_min)/(ndim*(rho(i) + rhomin))          ! deriv of this
!                 dhdrhoi = -hh(i)/(ndim*(rhoi + rhomin))          ! deriv of this
                 dwdhsumi = gradh(i)
                 omegai =  1. - dhdrhoi*gradh(i)
                 if (omegai.lt.1.e-5) then
                    !print*,'warning: omega < 1.e-5 ',i,omegai
                    if (abs(omegai).eq.0.) omegai = 1. !call quit
                 endif
                 gradh(i) = 1./omegai   ! this is what *multiplies* the kernel gradient in rates etc
                 func = rhoi - rho(i)
                 dfdh = omegai/dhdrhoi
                 gradsoft(i) = gradsoft(i)*dhdrhoi
                 !--gradgradhi is the "zeta" term in Price (2010)
                 d2hdrho2i = hh(i)*(ndim+1)/(rho(i)*ndim)**2
                 gradgradh(i) = rho(i)*(d2hdrho2i*dwdhsumi + dhdrhoi**2*gradgradh(i))
              endif
!
!--perform Newton-Raphson iteration to get new h
!                    
              hnew = hh(i) - func/dfdh
              if (hnew.gt.1.2*hh(i)) then
                 hnew = 1.2*hh(i)
              elseif (hnew.lt.0.8*hh(i)) then
                 hnew = 0.8*hh(i)
              endif
!
!--overwrite if iterations are going wrong
!
              if (usenumdens .and. (hnew.le.0 .or. gradhn(i).le.0)) then
                 nwarn = nwarn + 1
!                 !print*,' warning: h or omega < 0 in iterations',i,hnew,gradhn(i)
                 hnew = hfact*(1./rhoalt(i))**dndim   ! ie h proportional to 1/n^dimen
              elseif (.not. usenumdens .and. (hnew.le.0 .or. gradh(i).le.tiny(gradh))) then
!                 nwarn = nwarn + 1
                 hnew = hfact*(pmass(i)/(rho(i)+rhomin))**dndim   ! ie h proportional to 1/rho^dimen
              elseif (itsdensity.gt.100) then
                 hnew = hfact*(pmass(i)/(rho(i)+rhomin))**dndim   ! ie h proportional to 1/rho^dimen              
              endif
              if (numneigh(i).le.1) then
                 nwarn = nwarn + 1
                 !write(iprint,*) ' WARNING: particle ',i,' has no neighbours, increasing h'
                 hnew = hh(i) + psep
                 !write(iprint,*) 'NO NEIGHBOURS : rho,h = ',rho(i),hh(i),'setting h = ',hnew
                 redolink = .true.
              endif
!
!--if this particle is not converged, add to list of particles to recalculate
!              
              !!converged = abs(func)/rhoalt(i) < tolh .and. omegai > 0.
              converged = (abs((hnew-hh(i))/hhin(i)) < tolh .and. omegai > 0.) &
                          .or. itsdensitymax.eq.0
              
              testconvergence: if (.not.converged) then
                 ncalc = ncalc + 1
                 redolist(ncalc) = i
!            PRINT*,'not converged',i,abs(hnew-hh(i))/hh(i),rho(i),   &
!                 ncalc,redolist(ncalc)
!
!--update smoothing length only if taking another iteration
!
                 if (itsdensity.le.itsdensitymax .and. itype(i).ne.itypebnd) then
                    !print*,'hh new, old ',i,' = ',hnew,hh(i),abs((hnew-hh(i))/hh(i))
                    hh(i) = hnew
                 elseif (itsdensity.eq.itsdensitymax .and. .not.converged) then
                    write(iprint,*) 'ERROR: density not converged'
                    write(iprint,*) 'particle ',i,' h = ',hh(i),' hnew = ',hnew,'(rho -rhoi)/rho =',abs(func)/hhin(i)
                 endif
                 if (hnew.gt.hhmax) then
                    redolink = .true.
                 endif
              else
!
!--normalise arrays
!                 
                 if (usenumdens) then
                    dndt(i) = dndt(i)*gradhn(i)
                    dhdt(i) = dhdni*dndt(i)
                    drhodt(i) = drhodt(i) + gradh(i)*dndt(i)
                 else
                    if (ikernav.eq.3) drhodt(i) = drhodt(i)*gradh(i)
                    dhdt(i) = dhdrhoi*drhodt(i)
                 endif
              endif testconvergence
              
           endif   ! itype .NE. 1
        enddo
        if (nwarn.gt.0) then
           !write(iprint,*) ' WARNING: h or omega < 0 or no neighbours in iterations ',nwarn,' times'
        endif
        if (nrhosmall.gt.0) then
           write(iprint,"(a,i3,a,i8,a)") ' WARNING: iteration ',itsdensity,': rho < 1.e-6 on ',nrhosmall,' particles'
        endif

        if ((idebug(1:3).eq.'den').and.(ncalc.gt.0)) then
           write(iprint,*) ' density, iteration ',itsdensity,' ncalc = ',ncalc,':',redolist(1:ncalc)
        endif
     elseif (ihvar.gt.0) then
        do i=1,npart
           dhdrhoi = -hh(i)/(ndim*(rho(i) + rhomin))          ! deriv of this
           dhdt(i) = dhdrhoi*drhodt(i)
        enddo  
        if (ndim.gt.1) write(iprint,*) ' min. nneigh = ',minval(numneigh(1:npart)), &
                        ' max nneigh = ',maxval(numneigh(1:npart))
     else ! if ihvar = 0
!
!--overwrite gradh if not using variable h
! 
        do i=1,npart
           dhdt(i) = 0.
           gradh(i) = 1.
           gradhn(i) = 0.
           gradgradh(i) = 0.
        enddo   
     endif
!
!--write over boundary particles
!     
     if (any(ibound.eq.1)) then
        do i=1,npart      ! update fixed parts and ghosts
           if (itype(i).eq.itypebnd) then
              j = ireal(i)
              if (j > 0) then
                 rho(i) = rho(j)
                 rhoalt(i) = rhoalt(j)
                 drhodt(i) = drhodt(j)
                 dhdt(i) = dhdt(j)
                 hh(i) = hh(j)
                 gradh(i) = gradh(j)
                 gradhn(i) = gradhn(j)
                 gradsoft(i) = gradsoft(j)
              else
                 !!rho(i) = rhoin(i)
                 !write(iprint,*) 'Warning: ireal not set for fixed parts'
              endif
           endif
        enddo
     endif
     if (any(ibound > 1)) then   ! update ghosts
        do i=npart+1,ntotal
           j = ireal(i)
           if (j > 0) then
              rho(i) = rho(j)
              rhoalt(i) = rhoalt(j)
              drhodt(i) = drhodt(j)
              dhdt(i) = dhdt(j)
              hh(i) = hh(j)
              gradh(i) = gradh(j)
              gradhn(i) = gradhn(j)
              gradsoft(i) = gradsoft(j)
           endif
        enddo
     endif

  enddo iterate
  

!--NB: itsdensity is also used in step  
  if (itsdensity.gt.itsdensitymax .and. itsdensitymax.gt.0) then
     write(iprint,*) ' ERROR: DENSITY NOT CONVERGED ON ',ncalc,' PARTICLES'
     call quit
  elseif (itsdensity > 5 .or. usenumdens) then
     write(iprint,"(a,i2,a,f6.3,a,i3,a,i5)") &
      ' Density, its = ',itsdensity,' mean: ',ncalctotal/real(npart),&
      ' neigh min: ',minval(numneigh(1:npart)), &
      ' max: ',maxval(numneigh(1:npart))
  endif

  return
end subroutine iterate_density
