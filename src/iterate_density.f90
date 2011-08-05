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
!!------------------------------------------------------------------------

subroutine iterate_density
  use dimen_mhd
  use debug
  use loguns
  
  use bound
  use hterms
  use options
  use part
  use setup_params
  use density_summations
  use linklist, only:numneigh
!
!--define local variables
!
  implicit none
  integer :: i,j,itsdensitymax
  integer :: ncalc,ncalcprev,ncalctotal,isize
  integer, dimension(npart) :: redolist, redolistprev
  real :: tol,hnew,htemp,func,dfdh,deltarho
  real :: rhoi,dhdrhoi,omegai
  real, dimension(npart) :: rho_old
  logical :: converged,redolink
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
  tol = 1.e-4
  ncalctotal = 0
  ncalc = npart   ! number of particles to calculate density on
  redolink = .false.
  ncalcprev = 0
  gradh = 0.
  if (ncalc.eq.npart) then
     do j=1,npart     
        redolist(j) = j
     enddo
  endif
  
  iterate: do while ((ncalc.gt.0).and.(itsdensity.le.itsdensitymax)) !!    &
    !!!!   .and. ncalc.ne.ncalcprev)
     
     itsdensity = itsdensity + 1
     rho_old(1:npart) = rho(1:npart)

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
        call density(x,pmass,hh,rho,gradh,npart) ! symmetric for particle pairs
        !call output(0.0,1)
     else
        call density_partial(x,pmass,hh,rho,gradh,npart,ncalc,redolist)
     endif
     
     ncalctotal = ncalctotal + ncalc
     ncalcprev = ncalc
     redolistprev(1:ncalcprev) = redolist(1:ncalcprev)
     ncalc = 0
     redolink = .false.
     
     if (ihvar.ne.0) then
        do j=1,ncalcprev
           i = redolistprev(j)

           if (itype(i).ne.1) then
              if (rho(i).le.1.e-3) then
                 if (rho(i).le.0.) then
                    write(iprint,*) 'rho: rho -ve, dying rho(',i,') = ',rho(i),hh(i),pmass(i)
                    call quit
                 else
                    write(iprint,*) 'Warning : rho < 1e-3 '
                 endif
              endif
              
              rhoi = pmass(i)/(hh(i)/hfact)**ndim  ! this is the rho compatible with the old h
              dhdrhoi = -hh(i)/(ndim*rhoi)          ! deriv of this
              omegai =  1. - dhdrhoi*gradh(i)
              if (omegai.lt.1.e-5) then
                 print*,'warning: omega < 1.e-5 ',i,omegai
                 if (abs(omegai).eq.0.) call quit
              endif
              gradh(i) = 1./omegai   ! this is what *multiplies* the kernel gradient in rates etc
!
!--perform Newton-Raphson iteration to get new h
!      
              func = rhoi - rho(i)
              dfdh = omegai/dhdrhoi
              
              hnew = hh(i) - func/dfdh
!
!--overwrite if iterations are going wrong
!
              if (hnew.le.0. .or. gradh(i).le.0.) then
                 print*,' warning: h or omega < 0 in iterations ',i,hnew,gradh(i)
                 hnew = hfact*(pmass(i)/rho(i))**dndim   ! ie h proportional to 1/rho^dimen
              endif
              !if (numneigh(i).le.1) then
              !   print*,'NO NEIGHBOURS : rho = ',rho(i),' h = ',hnew,hh(i)
              !   print*,' ERROR: particle has no neighbours, increasing h'
              !   hnew = max(hh(i),hnew) + psep
              !endif
!
!--if this particle is not converged, add to list of particles to recalculate
!
!             PRINT*,'hnew - hh(i) = ',abs(hnew-hh(i))/hh(i)
              
              !!converged = abs((rho(i)-rhoi)/rho(i)) < tol .and. omegai > 0.
              converged = abs((hnew-hh(i))/hh(i)) < tol .and. omegai > 0.
              if (.not.converged) then
                 ncalc = ncalc + 1
                 redolist(ncalc) = i
!            PRINT*,'not converged',i,abs(hnew-hh(i))/hh(i),rho(i),   &
!                 ncalc,redolist(ncalc)
!
!--update smoothing length only if taking another iteration
!
                 if (itsdensity.lt.itsdensitymax .and. itype(i).ne.1) then
                    !print*,'hh new, old ',i,' = ',hnew,hh(i),abs((hnew-hh(i))/hh(i))
                    hh(i) = hnew
                 elseif (itsdensity.eq.itsdensitymax) then
                    write(iprint,*) 'ERROR: density not converged'
                    write(iprint,*) 'particle ',i,' h = ',hh(i),' hnew = ',hnew,'(h - hnew) /h = ',abs(hnew-hh(i))/hh(i)
                 endif
                 if (hnew.gt.hhmax) then
                    redolink = .true.
                 endif
              endif
           endif   ! itype .NE. 1
        enddo
        
        if ((idebug(1:3).eq.'den').and.(ncalc.gt.0)) then
           write(iprint,*) ' density, iteration ',itsdensity,' ncalc = ',ncalc,':',redolist(1:ncalc)
        endif
        
     endif
!
!--write over boundary particles
!     
     if (any(ibound.eq.1)) then
        do i=1,npart      ! update fixed parts and ghosts
           if (itype(i).eq.1) then
              j = ireal(i)
              if (j.ne.0) then
                 rho(i) = rho(j)
                 hh(i) = hh(j)
                 gradh(i) = gradh(j)
              else
                 rho(i) = rho_old(i)
                 !!write(iprint,*) 'Warning: ireal not set for fixed parts'
              endif
           endif
        enddo
     endif
     if (any(ibound.gt.1)) then   ! update ghosts
        do i=npart+1,ntotal
           j = ireal(i)
           rho(i) = rho(j)
           hh(i) = hh(j)
           gradh(i) = gradh(j)
        enddo
     endif
     
     !!call output(0.0,1)
     !!read*
  enddo iterate

!--NB: itsdensity is also used in step  
  if (itsdensity.gt.2 .or. (itsdensity.gt.1 .and. ndim.ge.2)) then
     write(iprint,*) ' Finished density, iterations = ',itsdensity, ncalctotal
  endif
  
  return
  
 contains
 
!------------------------------------------------------------
! these functions define the relationship between h and rho
!------------------------------------------------------------
  real function hrho(hfacti,pmassi,rhoi)
   use dimen_mhd
   use loguns
   implicit none
   real :: pmassi,rhoi,hfacti
   
   if (rhoi.le.0.) write(iprint,*) '*** ERROR: rho < 0 in hrho'
   hrho = hfacti*(pmassi/rhoi)**dndim
    
  end function hrho
  
  real function rhoh(hfacti,pmassi,hi)
   use dimen_mhd
   implicit none
   real :: hfacti, pmassi, hi
   
   rhoh = pmassi*(hfacti/hi)**dndim
  
  end function rhoh
  
  real function dhdrho(hi,rhoi)
   use dimen_mhd
   implicit none
   real :: hi, rhoi
  
   dhdrho = -hi/(ndim*rhoi)
  
  end function dhdrho
  
  real function drhodh(rhoi,hi)
   use dimen_mhd
   implicit none
   real :: rhoi, hi
   
   drhodh = -ndim*rhoi/hi
   
  end function drhodh
  
end subroutine iterate_density
