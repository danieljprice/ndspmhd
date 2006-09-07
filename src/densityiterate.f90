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
 use dimen_mhd, only:ndim,dndim
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
 use kernels, only:interpolate_kernel_dens,radkern,radkern2
 implicit none
 real, parameter :: pi = 3.1415926536
 integer :: i,j,n,icell,icellprev,nneigh,nneighi,itsdensity
 integer :: minneigh,maxneigh,ncalc,nitsmax
 integer, dimension(ntotal) :: listneigh
 real :: const,sum1,sum2,hi,hi21,hnew
 real :: func0,func,dfunc
 real, dimension(ndim) :: dx
 real :: rij2,q2i,pmassj,wabi,grkerni,grgrkerni
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
 print*,'const = ',const,'hhmax=',hhmax
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
       if (hi.gt.hhmax) then
          if (any(ibound.gt.1)) call set_ghost_particles
          write(iprint,*) 'relinking...'
          call set_linklist
          call get_neighbour_list_partial(icell,listneigh,nneigh)
       elseif (icell.NE.icellprev) then
          call get_neighbour_list_partial(icell,listneigh,nneigh)
       endif
       icellprev = icell
       sum1 = 0.
       sum2 = 0.
       hi21 = 1./hi**2
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
             call interpolate_kernel_dens(q2i,wabi,grkerni,grgrkerni)             
             sum1 = sum1 + pmassj*wabi
             sum2 = sum2 + pmassj*grgrkerni
             if (i.ne.j) sum2 = sum2 + pmassj*(ndim-1)*grkerni/sqrt(q2i)
!                print*,'q2i,w = ',q2i,wabi,grgrkerni
          endif
       enddo loop_over_neighbours
!
!--work out sums, see if h is converged
!
       func = sum1 + const*sum2
!       print*,'sum1, sum2 = ',sum1,const*sum2
       dfunc = func - func0
       
       converged = (abs(dfunc/func) < tolh) .or. itsdensity.ge.maxdensits
!
!--set h for next iterations if not converged
!
       if (.not.converged) then
          !--fixed point iteration
!          print*,i,'hold = ',hi,' neigh = ',nneighi,'rho = ',sum1/(hi**ndim)
!          print*,'del2rho =',sum2/(hi**(ndim+2))
!          print*,'fact=',hfact*(pmass(i)/(sum1 + const*sum2))**dndim
          hi = hfact*(pmass(i)/(sum1 + const*sum2))**dndim*hi
!          select case(ndim)
!          case(1)
!          print*,'Neigh 1 =',2.*radkern*sum1/pmass(i)
!          print*,'Neigh 2 =',(2.*radkern*sum1 + radkern**3*2./3.*sum2)/pmass(i)
!          case(3)
!          print*,'Neigh 1 =',4./3.*pi*radkern**3*sum1/pmass(i)
!          print*,'Neigh 2 =',(4./3.*pi*radkern**3*sum1 + 2.*pi/5.*radkern**5*sum2)/pmass(i)
!          end select
!          read*
       endif
       if (nneighi.le.1) then
          !print*,'NO NEIGHBOURS : rho = ',rho(i),' h = ',hnew,hh(i)
          print*,' WARNING: particle ',i,' has no neighbours, increasing h'
          hi = hi + psep
          converged = .false.
       endif

    enddo iterate

    nitsmax = max(nitsmax,itsdensity)
    if (itsdensity.eq.maxdensits) write(iprint,*) 'ERROR: not converged ',i
!
!--having finished iterations, now calculate extra crap
!
    rho(i) = sum1/(hi**ndim)
    hh(i) = hi
    drhodt(i) = 0.
    dhdt(i) = 0.
    gradh(i) = 1.
    gradsoft(i) = 0.
    minneigh = min(minneigh,nneighi)
    maxneigh = max(maxneigh,nneighi)
!    print*,i,'its = ',itsdensity,'neigh = ',nneighi,' rho = ',rho(i)
!    read*

 enddo
 write(iprint,*) ' min. nneigh = ',minneigh,' max nneigh = ',maxneigh
 write(iprint,*) ' iterations = ',nitsmax,ncalc

end subroutine densityiterate
