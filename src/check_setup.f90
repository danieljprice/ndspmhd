!
!--checks (user defined) setup of primitive quantities for errors
!
subroutine check_setup
  use debug
  use loguns
  use part
  implicit none
  integer :: i
  
  write(iprint,5,ADVANCE='NO') ' Checking setup... '
5 format(/,a)

  if (npart.le.0) then
     write(iprint,10) 'no particles (npart = 0)'
     stop
  endif
  if (.not.allocated(rho)) then
     write(iprint,10) 'memory not allocated'
     stop
  endif

  do i=1,npart
!
!--check for negative densities/thermal energies
!
     if (dens(i).le.0.) then
        write(iprint,20) 'density <= 0 ', dens(i)
        stop
     endif
     if (uu(i).le.0.) then
        write(iprint,20) 'uu <= 0', uu(i)
        stop
     endif
  enddo

10 format(/,'ERROR IN PARTICLE SETUP:',a,/)
20 format(/,'ERROR IN PARTICLE SETUP:',a,1pe10.3,/)

   write(iprint,*) 'OK'

end subroutine check_setup
