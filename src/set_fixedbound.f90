subroutine set_fixedbound
  use dimen_mhd
  use debug
  use loguns
  use bound
  use options
  use hterms
  use part
  implicit none
  integer :: i,j,nstart,npart1
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
        itype(1:nbpts) = 1
        nstart = npart - nbpts + 1
        itype(nstart:ntotal) = 1
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
        call set_ghost_particles  ! setup ghost particles  
        !--now fix the ghost particles that have been set
        npart1 = npart + 1
        itype(npart1:ntotal) = 1  ! set all these particles to be fixed
        npart = ntotal              ! no ghosts
     endif
  else
!
!--if fixed particles have been set in the setup routine, assign the ireal
!  array so that the density can be copied from particles just outside the 
!  fixed zone
!    
     nbpts = 0
     do i=1,npart
        if (itype(i).eq.1 .or. itype(i).eq.2) nbpts = nbpts + 1
     enddo
     write(iprint,*) nbpts,' fixed particles set: finding nearest real parts' 
     do i=1,npart
        if (itype(i).eq.1) then
           rmin = 1.e10
           do j=1,npart
              if (j.ne.i .and. itype(j).ne.1) then   ! find closest real particle
                 dx = x(:,i) - x(:,j)
                 rr = DOT_PRODUCT(dx,dx)
                 if (rr.lt.rmin) then
                    rmin = rr
                    ireal(i) = j
                 endif
              endif
           enddo
           !ireal(i) = 0
           if (ireal(i).eq.0) stop 'error finding nearest particle to fixed part'
           if (idebug(1:5).eq.'fixed') write(iprint,*) ' particle ',i,' copied from ',ireal(i)
        endif
     enddo
  endif

10 format(/,a,/)  
  
end subroutine set_fixedbound
