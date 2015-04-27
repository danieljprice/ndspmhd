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

subroutine set_fixedbound
  use dimen_mhd
  use debug
  use loguns
  use bound
  use options
  use hterms
  use part
  implicit none
  integer :: i,j,nstart,npart1,iboundold(ndim)
  real :: rr, rmin
  real, dimension(ndim) :: dx
  
  if (trace) write(iprint,*) 'Entering subroutine set_fixedbound'
  
  if (all(itype.eq.0)) then
!
!--in 1D setup can just specify the number of particles to hold fixed
!  this is for backwards compatibility of the code
!
     if (ndim.EQ.1 .AND. nbpts.GT.0) then
        write(iprint,*) 'fixing first and last ',nbpts,' particles'
        itype(1:nbpts) = itypebnd
        nstart = npart - nbpts + 1
        itype(nstart:ntotal) = itypebnd
        !!--fix the values of rho, h equal to those just outside the fixed zone
        ireal(1:nbpts) = nbpts+1
        rho(1:nbpts) = rho(nbpts+1)
        hh(1:nbpts) = hh(nbpts+1)
        gradh(1:nbpts) = gradh(nbpts+1)
        ireal(nstart:npart) = npart-nbpts
        rho(nstart:npart) = rho(npart-nbpts)
        hh(nstart:npart) = hh(npart-nbpts)
        gradh(nstart:npart) = gradh(npart-nbpts)
     elseif (nbpts.eq.0) then
!
!--in >1D we setup ghost particles initially but then fix them
! 
        write(iprint,10) ' Setting up fixed particles from ghosts'
        iboundold = ibound
        ibound = 0
        where (iboundold.eq.1) ibound = 2
        call set_ghost_particles  ! setup ghost particles
        !--now fix the ghost particles that have been set
        npart1 = npart + 1
        itype(npart1:ntotal) = itypebnd  ! set all these particles to be fixed
        npart = ntotal              ! no ghosts
        do i=1,npart
           if (itype(i).eq.itypebnd .or. itype(i).eq.itypebnd2) nbpts = nbpts + 1
        enddo
        write(iprint,*) nbpts,' fixed particles set'
        ibound = iboundold
     endif
  else
!
!--if fixed particles have been set in the setup routine, assign the ireal
!  array so that the density can be copied from particles just outside the 
!  fixed zone
!    
     nbpts = 0
     do i=1,npart
        if (itype(i).eq.itypebnd .or. itype(i).eq.itypebnd2) nbpts = nbpts + 1
     enddo
     write(iprint,*) nbpts,' fixed particles set: finding nearest real parts' 
     do i=1,npart
        if (itype(i).eq.itypebnd) then
           rmin = huge(rmin)
           do j=1,npart
              if (j.ne.i) then   ! find closest real particle
                 select case(itype(j))
                 case(itypegas,itypegas1,itypegas2)
                    dx = x(:,i) - x(:,j)
                    rr = DOT_PRODUCT(dx,dx)
                    if (rr.lt.rmin) then
                       rmin = rr
                       ireal(i) = j
                    endif
                 end select
              endif
           enddo
           !ireal(i) = 0
           !idebug = 'fixed'
           !if (ireal(i).eq.0) stop 'error finding nearest particle to fixed part'
           if (idebug(1:5).eq.'fixed') write(iprint,*) ' particle ',i,' copied from ',ireal(i)
        endif
     enddo
  endif

10 format(/,a,/)  
  
end subroutine set_fixedbound
