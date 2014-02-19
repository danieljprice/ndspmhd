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

!
!--checks (user defined) setup of primitive quantities for errors
!
subroutine check_setup
  use dimen_mhd, only:ndim
  use debug
  use loguns
  use part
  use options, only:igravity,ibound
  use bound,   only:xmin,xmax
  implicit none
  integer :: i,j
  real, parameter :: vbig = 1.e12
  real, dimension(ndim) :: xcentre
  real :: dxmean,dxij
  
  write(iprint,5) ' Checking setup... '
5 format(/,a)

  if (npart.le.0) then
     write(iprint,10) 'no particles (npart = 0)'
     stop
  endif
  if (.not.allocated(rho)) then
     write(iprint,10) 'memory not allocated'
     stop
  endif

  xcentre(:) = 0.
  do i=1,npart
!
!--check for negative densities/thermal energies
!
     if (dens(i).lt.tiny(dens)) then
        write(iprint,20) 'density <= 0 ',dens(i),i
        stop
!     elseif (isnan(dens(i))) then
!        write(iprint,20) 'NaNs in density',dens(i),i
!        stop
     endif
     if (pmass(i).lt.tiny(pmass)) then
        write(iprint,20) 'pmass <= 0 ',pmass(i),i
        stop
     elseif (pmass(i).gt.vbig) then
        write(iprint,20) 'huge pmass values ',pmass(i),i
        stop
     endif     
     if (uu(i).lt.0.) then
        write(iprint,20) 'uu < 0 ',uu(i),i
        stop
     endif
     if (hh(i).lt.0.) then
        write(iprint,20) 'h <= 0',hh(i),i
     endif
!
!--check for huge velocities
!
     if (any(abs(vel(1:ndimV,i)).gt.vbig)) then
        write(iprint,10) ' contains huge velocities!! '
        stop
!     else
!        do j=1,ndimV
!           if (isnan(vel(j,i))) then
!              write(iprint,20) 'NaNs in velocities ',vel(j,i),i
!              stop
!           endif
!        enddo
     endif
     
     xcentre(:) = xcentre(:) + pmass(i)*x(:,i)
!
!--check for particles outside the boundary
!
     do j=1,ndim
        if (ibound(j).gt.0) then
           if (x(j,i).lt.xmin(j) .or. x(j,i).gt.xmax(j)) then
              write(iprint,20) 'particle outside boundary',x(1,i),i
              stop
           endif
        endif
     enddo
     
     if (itype(i).lt.0 .or. itype(i).gt.1000) then
        write(iprint,30) 'particle type out of range ',itype(i),i
        stop
     endif
  enddo
  
  if (all((abs(pmass(1:ntotal)-pmass(1)).lt.epsilon(0.)))) then
     write(iprint,"(a)") ' -> equal mass particles'
  else
     write(iprint,"(a)") ' WARNING: Unequal mass particles used' 
  endif
!
!--check that no particles are on top of each other
!
!  do i=1,npart
!     do j=i+1,npart
!        xsep = dot_product(x(:,i)-x(:,j),x(:,i)-x(:,j))
!        if (xsep.lt.tiny(xsep)) then
!           write(iprint,20) 'non-unique particle position, x= ',x(1,i),i
!           write(iprint,20) 'non-unique particle position, x= ',x(1,j),j
!           stop
!        endif
!     enddo
!  enddo
!
!--find mean particle separation
! 
 if (igravity.eq.1 .or. igravity.eq.2) then
    dxmean = 0.
    do i=1,ntotal
       do j=i+1,ntotal
          dxij = sqrt(dot_product(x(:,i)-x(:,j),x(:,i)-x(:,j)))
          dxmean = dxmean + dxij
       enddo
    enddo
    dxmean = dxmean/real((ntotal**2 - ntotal)/2)
    write(iprint,*) 'mean particle spacing = ',dxmean
    write(iprint,*) 'suggested softening h = ',dxmean/40.,' to ',dxmean/35.
 endif
!
!--warnings only
! 
  if (any(uu(1:npart).eq.0.)) then
     write(iprint,*) ' WARNING: uu = 0 on some particles'
  endif
  
  xcentre = xcentre / SUM(pmass(1:npart))

10 format(/,' ERROR IN PARTICLE SETUP: ',a,/)
20 format(/,' ERROR IN PARTICLE SETUP: ',a,1pe10.3,' particle ',i5,/)
30 format(/,' ERROR IN PARTICLE SETUP: ',a,i10,' particle ',i5,/)

   write(iprint,*) '-> centre of mass is at : ',xcentre(:)

   write(iprint,*) '-> OK'

end subroutine check_setup
