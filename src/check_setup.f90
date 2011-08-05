!
!--checks (user defined) setup of primitive quantities for errors
!
subroutine check_setup
  use dimen_mhd, only:ndim
  use debug
  use loguns
  use part
  implicit none
  integer :: i
  real, parameter :: vbig = 1.e12
  real, dimension(ndim) :: xcentre
  
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
!
!--check for huge velocities
!
     if (any(abs(vel(1:ndimV,i)).gt.vbig)) then
        write(iprint,10) ' contains huge velocities!! '
        stop
     endif
     
     xcentre(:) = xcentre(:) + pmass(i)*x(:,i)
     
  enddo

!
!--warnings only
! 
  if (any(uu(1:npart).eq.0.)) then
     write(iprint,*) ' WARNING: uu = 0 on some particles'
  endif
  
  xcentre = xcentre / SUM(pmass(1:npart))

10 format(/,' ERROR IN PARTICLE SETUP: ',a,/)
20 format(/,' ERROR IN PARTICLE SETUP: ',a,1pe10.3,' particle ',i5,/)

   write(iprint,*) '-> centre of mass is at : ',xcentre(:)

   write(iprint,*) '-> OK'

end subroutine check_setup
