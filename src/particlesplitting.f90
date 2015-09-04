!!------------------------------------------------------------------------!!
!! Particle splitting module (VERY EXPERIMENTAL)
!!------------------------------------------------------------------------!!
module particlesplit
 implicit none
  
contains

subroutine particle_splitting(nsplit)
 use dimen_mhd
 use debug
 use loguns
 use options, only:rhocrit
 use part, only:npart,ntotal,rho,pmass
 implicit none
 integer, intent(out) :: nsplit
 integer :: npartnew,i
 real, parameter :: pmassinit = 1./5000.
!
!--allow for tracing flow
!
!
!--initialise quantities
!
 npartnew = npart
 do i=1,npart
    if (rho(i).gt.rhocrit .and. pmass(i).ge.pmassinit) then
       print*,'splitting ',i,'rho, rhocrit = ',rho(i),rhocrit
       call splitpart(i,npartnew)
    endif
 enddo

 npart = npartnew
 ntotal = npart

 nsplit = npartnew - npart

 return
end subroutine particle_splitting

subroutine splitpart(i,npartnew)
 use part
 use mem_allocation, only:alloc
 use random,         only:ran1
 implicit none
 integer, intent(in) :: i
 integer, intent(inout) :: npartnew
 integer :: j,npartold
 integer, parameter :: nmake = 5
 real :: pmassnew
 integer :: iseed = -123
 save iseed
 

 npartold = npartnew
 npartnew = npartold + nmake
 !print*,'iseed = ',iseed
 
 pmassnew = pmass(i)/real(nmake)
 !--create new particle(s)
 if (npartnew.gt.size(x(1,:))) call alloc(npart+100*nmake)
 do j=npartold+1,npartnew
    call copy_particle(j,i)
    x(:,j) = x(:,i)*(1. + 0.001*hh(i)*(ran1(iseed)-0.5))
    !print*,j,' xnew = ',x(:,j),' xold = ',x(:,j)
    vel(:,j) = vel(:,i)
    pmass(j) = pmassnew
 enddo
 pmass(i) = pmassnew
 
end subroutine splitpart

end module particlesplit
