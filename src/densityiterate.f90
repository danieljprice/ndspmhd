!!------------------------------------------------------------------------
!! Iterates the density calculation so that rho and h are
!! calculated self-consistently.
!!
!! In this version the iterations are performed in an inner loop
!! for each particle.
!!
!! Uses my complicated dependence of h on rho to better maintain the
!! exact neighbour number
!!
!!------------------------------------------------------------------------

subroutine densityiterate
 use dimen_mhd, only:ndim
 use debug, only:trace
 use loguns, only:iprint
 
 use bound, only:hhmax
 use hterms, only:gradh,gradsoft
 use linklist, only:iamincell
 use options, only:maxdensits,tolh,ibound
 use part, only:x,pmass,hh,rho,npart,ntotal
 use rates, only:dhdt,drhodt
 use setup_params, only:hfact,psep
 use get_neighbour_lists, only:get_neighbour_list_partial
 use kernels, only:interpolate_kernel,radkern2
 implicit none
 real, parameter :: pi = 3.1415926536
 integer :: i,j,n,icell,icellprev,nneigh,nneighi,itsdensity
 integer :: minneigh,maxneigh,ncalc,nitsmax
 integer, dimension(ntotal) :: listneigh
 real :: const,sum1,sum2,dwdhterm1,hi,hi1,hi21,hnew,dwdhi
 real :: func0,func,dfunc,dfdh,dhdrhoi,omegai,rhomin,rhoi
 real, dimension(ndim) :: dx
 real :: rij2,q2i,pmassj,wabi,grkerni !,grgrkerni
 logical :: converged
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' Entering subroutine densityiterate' 
!
!--set "constants"
!
 select case(ndim)
 case(1)
    const = radkern2/6.
 case(2)
    const = 0.25*radkern2
 case(3)
    const = 0.3*radkern2
 end select
 const = 0.
 rhomin = 0.
 icellprev = 0
 minneigh = huge(minneigh)
 maxneigh = 0
!
!--loop over all particles
!
 ncalc = 0
 nitsmax = 0
 
 do i=1,npart
!
!--set sums to zero, also set initial values of h, rho etc
!
    hi = hh(i)
    func0 = pmass(i)*hfact**ndim
    itsdensity = 0
!
!--begin iteration loop
!
    converged = .false.
    iterate: do while (.not.converged)    
       itsdensity = itsdensity + 1
       ncalc = ncalc + 1

!
!--get neighbour list for current particle
! 
       icell = iamincell(i)
!
!--if different to previous cell used, get the list of neighbours for this cell
!  (common to all particles in the cell)
!
       if (hh(i).gt.hhmax) then
          if (any(ibound.gt.1)) call set_ghost_particles
          write(iprint,*) 'relinking...hi = ',hi,hhmax
          call set_linklist
          icell = iamincell(i)
          call get_neighbour_list_partial(icell,listneigh,nneigh)
          write(iprint,*) 'new nneigh = ',nneigh
       elseif (icell.NE.icellprev) then
          call get_neighbour_list_partial(icell,listneigh,nneigh)
       endif
       icellprev = icell
       sum1 = 0.
       sum2 = 0.
       dwdhterm1 = 0.
       hi = hh(i)
       hi1 = 1./hi
       hi21 = hi1**2
       nneighi = 0
!
!--loop over current particle's neighbours
!
       loop_over_neighbours: do n = 1,nneigh
          j = listneigh(n)
          dx(:) = x(:,i) - x(:,j)
          rij2 = dot_product(dx,dx)
          q2i = rij2*hi21
!      
!--do interaction if r/h < compact support size
!
          if (q2i.LT.radkern2) then
             nneighi = nneighi + 1
             pmassj = pmass(j)
!      
!--interpolate from kernel table (using hi)
!
             call interpolate_kernel(q2i,wabi,grkerni)             
             sum1 = sum1 + pmassj*wabi
             dwdhterm1 = dwdhterm1 - pmassj*grkerni*sqrt(q2i)
!             sum2 = sum2 + pmassj*grgrkerni
!             if (i.ne.j) sum2 = sum2 + pmassj*(ndim-1)*grkerni/sqrt(q2i)
!                print*,'q2i,w = ',q2i,wabi,grgrkerni
          endif
       enddo loop_over_neighbours
!
!--work out sums, see if h is converged
!
       func = sum1 !!+ const*sum2
       dfdh = dwdhterm1*hi1
       print*,i,'iteration',itsdensity,' rho = ',sum1*hi1**ndim,'h = ',hi
       print*,'  nneigh = ',nneighi, 'dfdh = ',dfdh
       dfunc = func - func0
       
       converged = (abs(dfunc/func) < tolh) .or. itsdensity.ge.maxdensits

          rho(i) = sum1*hi1**ndim
          rhoi = pmass(i)/(hh(i)/hfact)**ndim - rhomin ! this is the rho compatible with the old h
          dwdhi = dwdhterm1*hi1**(ndim+1) - ndim*hi1*rhoi
          dhdrhoi = -hi/(ndim*(rhoi + rhomin))
          omegai = 1. - dhdrhoi*dwdhi
          gradh(i) = 1./omegai   ! this is what *multiplies* the kernel gradient in rates etc
          func = rhoi - rho(i)
          dfdh = omegai/dhdrhoi
          hnew = hh(i) - func/dfdh
          converged = (abs((hnew-hh(i))/hi) < tolh .and. omegai > 0.)

       if (nneighi.le.1) then
          !print*,'NO NEIGHBOURS : rho = ',rho(i),' h = ',hnew,hh(i)
          print*,' WARNING: particle ',i,' has no neighbours, increasing h'
          hnew = hh(i) + psep
          converged = .false.
       endif
!
!--set h for next iterations if not converged
!
       if (.not.converged) then
!
!--perform Newton-Raphson iteration to get new h
!

!          hnew = hi - func/dfdh
          !--fixed point iteration
          print*,'h new = ',hnew,' func/dfdh = ',func/dfdh
!          if (hnew.lt.0.) then
!             write(iprint,*) 'taking fixed point, h = ',hfact*(pmass(i)/sum1)**dndim*hi,hi - func/dfdh
!             hnew = hfact*(pmass(i)/sum1)**dndim*hi
!          endif
          hh(i) = hnew
!          read*
       endif


    enddo iterate

    nitsmax = max(nitsmax,itsdensity)
    if (itsdensity.eq.maxdensits) write(iprint,*) 'ERROR: not converged ',i
!
!--having finished iterations, now calculate extra crap
!
    rhoi = sum1*hi1**ndim
    rho(i) = sum1*hi1**ndim
    hh(i) = hi
    dwdhi = dwdhterm1*hi1**(ndim+1) - ndim*hi1*rhoi
    dhdrhoi = -hi/(ndim*(rho(i) + rhomin))
    omegai = 1. - dhdrhoi*dwdhi
    gradh(i) = 1./omegai
    gradsoft(i) = 0.    
    drhodt(i) = 0.
    dhdt(i) = 0.
    minneigh = min(minneigh,nneighi)
    maxneigh = max(maxneigh,nneighi)
!    print*,i,'its = ',itsdensity,'neigh = ',nneighi,' rho = ',rho(i)
!    read*

 enddo
 write(iprint,*) ' min. nneigh = ',minneigh,' max nneigh = ',maxneigh
 write(iprint,*) ' iterations = ',nitsmax,ncalc

end subroutine densityiterate
