subroutine set_fixedbound
  use dimen_mhd
  use debug
  use loguns
  use options
  use hterms
  use part
  implicit none
  integer :: i,nstart,npart1
  
  if (all(itype.eq.0)) then
!
!--in 1D setup can just specify the number of particles to hold fixed
!  this is for backwards compatibility of the code
!
     if (ndim.EQ.1 .AND. nbpts.GT.0) then
        write(iprint,*) 'fixing first and last ',nbpts,' particles'
        itype(1:nbpts) = 1
        nstart = npart - nbpts + 1
        itype(nstart:ntotal) = 1
        !!--fix the values of rho, h equal to those just outside the fixed zone
        rho(1:nbpts) = rho(nbpts+1)
        hh(1:nbpts) = hh(nbpts+1)
        gradh(1:nbpts) = gradh(nbpts+1)
        rho(nstart:npart) = rho(npart-nbpts)
        hh(nstart:npart) = hh(npart-nbpts)
        gradh(nstart:npart) = gradh(npart-nbpts)
     elseif (nbpts.eq.0) then
!
!--in >1D we setup ghost particles initially but then fix them
! 
        write(iprint,10) ' Setting up fixed particles'
        call set_ghost_particles  ! setup ghost particles  
        !--now fix the ghost particles that have been set
        npart1 = npart + 1
        itype(npart1:ntotal) = 1  ! set all these particles to be fixed
        npart = ntotal		  ! no ghosts
     endif
  endif
10 format(/,a,/)  
  
end subroutine set_fixedbound
