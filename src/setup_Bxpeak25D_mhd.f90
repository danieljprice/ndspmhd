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

!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Setup for the Bx peak advection problem in Dedner et al JCP 175, 645  !!
!!  Bx = r(x^2 + y^2)/sqrt(4pi) (ie div B .ne. 0)                         !!
!!  Basically to see how an initially non-zero div B propagates           !!
!!                                                                        !!
!!------------------------------------------------------------------------!!

subroutine setup
!
!--include relevant global variables
!
 use dimen_mhd
 use debug
 use loguns
 use bound
 use eos
 use options
 use part
 use setup_params
 
 use uniform_distributions
!
!--define local variables
!      
 implicit none
 integer :: ipart,ibump,ierr
 real :: denszero,przero
 real :: pri,rbump,rr,rbump2
 real :: totmass,gam1,massp,const
 real, dimension(ndim) :: xorigin, dx
 real, dimension(ndimv) :: Bzero
!
!--check number of dimensions is right
!
 if (ndimv.ne.3) stop ' need ndimv=3 for this problem'
!
!--set boundaries
!                        
 ibound = 3     ! reflective ghosts (boundaries not important in this problem)
 nbpts = 0      ! no fixed particles
 xmin(:) = -0.5 ! unit square
 xmax(:) = 1.5
 const = sqrt(4.*pi) 
!
!--setup parameters for the problem
! 
 xorigin(:) = 0.0 ! co-ordinates of the centre of the initial blast
 rbump = 1./sqrt(8.)        ! radius of the initial bump
! read(rootname(6:7),*,iostat=ierr) ibump
! if (ierr.eq.0) then
!    print*,'ibump = ',ibump
! else
!    read(rootname(6:6),*,iostat=ierr) ibump
!    if (ierr.eq.0) then
!       print*,'ibump = ',ibump
!    else
       ibump = 1
!    endif
! endif
 rbump = rbump/ibump  !!/16.
 rbump2 = rbump*rbump
 Bzero(:) = 0.
 if (imhd.ne.0) Bzero(3) = 1.0/const        ! uniform field in bz direction
 przero = 6.0                ! initial pressure
 denszero = 1.0                ! ambient density
 
 gam1 = gamma - 1.

 write(iprint,*) 'Two dimensional div B advection problem '
 write(iprint,10) denszero,rbump,Bzero(3),przero
10 format(/,' density  = ',f10.3,', size of bump = ',f6.3,/, &
            ' initial Bz   = ',f6.3,', pressure = ',f6.3,/)
!
!--setup uniform density grid of particles (2D) 
!  (determines particle number and allocates memory)
!
 call set_uniform_cartesian(2,psep,xmin,xmax,fill=.true.)        ! 2 = close packed arrangement

 ntotal = npart
!
!--determine particle mass in ambient medium
!
 totmass = denszero*product(xmax(:)-xmin(:))
 massp = totmass/float(ntotal) ! average particle mass
!
!--now assign particle properties
! 
 do ipart=1,ntotal
    dx(:) = x(:,ipart)-xorigin(:) 
    rr = dot_product(dx,dx)
    Bfield(:,ipart) = bzero(:)
    if (rr.le.rbump2) then
       Bfield(1,ipart) = ((rr/rbump2)**4 - 2.*(rr/rbump2)**2 + 1.)/const
    endif  
    pmass(ipart) = massp
    dens(ipart) = denszero
    vel(:,ipart) = 0.
    vel(1:ndim,ipart) = 1.0
    pri = przero 
    uu(ipart) = pri/(gam1*denszero)
 enddo
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
            
 return
end

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
