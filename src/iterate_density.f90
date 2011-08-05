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
  use part_in	! for rhoin on fixed particles
  use setup_params
!
!--define local variables
!
  implicit none
  integer :: i,j,itsdensity,itsdensitymax
  integer :: ncalc,ncalcprev,ncalctotal,isize
  integer, dimension(:), allocatable :: redolist, redolistprev
  real :: tol,hnew
  real, dimension(:), allocatable :: rho_old,gradh_old
  logical :: converged,redolink
!
!--allow for tracing flow
!      
  if (trace) write(iprint,*) ' Entering subroutine iterate_density' 
!
!--allocate memory for local array copies
!
  isize = size(rho)
  allocate( rho_old(isize), gradh_old(isize) )
  allocate( redolist(ntotal), redolistprev(ntotal) )
!
!--set maximum number of iterations to perform
! 
  if ((ikernav.eq.3).and.(ihvar.ne.0)) then
     itsdensitymax = 20	! perform 1 fixed point iteration
  else
     itsdensitymax = 0	! no iterations
  endif
!
!--Loop to find rho and h self-consistently (if using Springel/Hernquist)
!
  itsdensity = 0
  tol = 1.e-2
  ncalctotal = 0
  ncalc = npart	! number of particles to calculate density on
  redolink = .false.
  ncalcprev = 0
  if (ncalc.eq.npart) then
     do j=1,npart     
        redolist(j) = j
     enddo
  endif
  
  iterate: do while ((ncalc.gt.0).and.(itsdensity.le.itsdensitymax)) !!    &
    !!!!   .and. ncalc.ne.ncalcprev)
     
     itsdensity = itsdensity + 1
     if (isize.ne.size(rho)) then
        deallocate(rho_old,gradh_old)
        isize = size(rho)
        allocate(rho_old(isize),gradh_old(isize))
     endif
     rho_old = rho
     gradh_old = gradh
     if (redolink) then
        if (any(ibound.gt.1)) call set_ghost_particles
        if (idebug(1:3).eq.'den') write(iprint,*) 'relinking...'
        call link
     endif
     
     if (ncalc.eq.npart) then	! calculate density on all particles
        call density		! do this symmetrically
     else	                ! calculate density on partial list of particles              
        call density_partial(ncalc,redolist)
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
                    write(iprint,*) 'rho: rho -ve, dying rho(',i,') = ',rho(i)
                    call quit
                 else
                    write(iprint,*) 'Warning : rho < 1e-3 '
                 endif
              endif
              
              gradh(i) = gradh(i)/rho(i)    ! now that rho is known
              if (abs(1.-gradh(i)).lt.1.e-5) print*,'eek'
!
!--perform Newton-Raphson iteration on rho
!      
!             rho(i) = rho_old(i) - (rho_old(i) - rho(i))/(1.-gradh(i))
!
!--work out new smoothing length h
!
              hnew = hfact*(pmass(i)/rho(i))**hpower	! ie h proportional to 1/rho^dimen
!
!--if this particle is not converged, add to list of particles to recalculate
!
!             PRINT*,'hnew - hh(i) = ',abs(hnew-hh(i))/hh(i)
              
              converged = abs((hnew-hh(i))/hh(i)) < tol	  
              if (.not.converged) then
                 ncalc = ncalc + 1
                 redolist(ncalc) = i
!	         PRINT*,'not converged',i,abs(hnew-hh(i))/hh(i),rho(i),	&
!     	         ncalc,redolist(ncalc)
!
!--update smoothing length only if taking another iteration
!
                 if (itsdensity.lt.itsdensitymax) hh(i) = hnew
                 if (hnew.gt.hhmax) then
                    redolink = .true.
                 endif
              endif
           endif	! itype .NE. 1
        enddo
        
        if ((idebug(1:3).eq.'den').and.(ncalc.gt.0)) then
           write(iprint,*) ' density, iteration ',itsdensity,' ncalc = ',ncalc,':',redolist(1:ncalc)
        endif
        
     endif
     
     if (any(ibound.gt.1)) then
        do i=npart+1,ntotal		! update ghosts
           j = ireal(i)
           rho(i) = rho(j)
           hh(i) = hh(j)
           gradh(i) = gradh(j)       
        enddo
     elseif (any(ibound.eq.1)) then		! update fixed particles
        where (itype(:).eq.1)
           rho(:) = rhoin(:)
	   hh(:) = hhin(:)
        end where
     endif
     
  enddo iterate
  
  if (itsdensity.gt.1) write(iprint,*) ' Finished density, iterations = ', &
                                       itsdensity, ncalctotal
  
  return
end subroutine iterate_density
