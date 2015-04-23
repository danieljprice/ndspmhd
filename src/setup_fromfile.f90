!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Setup by reading a dump from a file                                   !!
!!  Should be compatible with the output format given in output           !!
!!                                                                        !!
!!------------------------------------------------------------------------!!

SUBROUTINE setup
!
!--include relevant global variables
!
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE bound
 USE derivB
 USE eos
 USE options
 USE part
 USE setup_params
 USE rates, only:force
 use cons2prim, only:primitive2conservative
!
!--define local variables
!      
 IMPLICIT NONE
 INTEGER :: i,ierr
 CHARACTER(LEN=30) :: setupfile
 LOGICAL :: iexist
 REAL :: dummy,frad,rhat(3),omegai,rpart,v2onr

!
!--set name of setup file
!
 setupfile = 'start.out'

 WRITE(iprint,*) 'Reading setup from file: ',TRIM(setupfile)
 INQUIRE (file=setupfile,exist=iexist)
 IF (.not.iexist) THEN
    WRITE(iprint,*) 'Can''t find setup file: ',TRIM(setupfile)    
    STOP
 ENDIF
!
!--open setup file
!
 OPEN(UNIT=ireadf,FILE=setupfile,ERR=666,STATUS='old',FORM='formatted')
!
!--read header line
!
 npart = 15000
 ntotal = npart
 
 CALL alloc(ntotal)
!
!--read one step (only) from file
!
 i = 0
 ierr = 0
 do while (i < ntotal .and. ierr .eq. 0)
    i = i + 1
    if (ndim.eq.3) then
       read(ireadf,*,iostat=ierr) x(:,i),dummy,vel(:,i),pmass(i),dens(i),spsound(i)
    else
       read(ireadf,*,iostat=ierr) x(1:2,i),dummy,dummy,vel(:,i),pmass(i),dens(i),spsound(i)
    endif
!    print*,i,x(:,i),dummy,vel(:,i),pmass(i)
    if (ierr .eq. 0) then
       uu(i) = spsound(i)**2/(gamma*(gamma-1.))
       pr(i) = (gamma-1.)*uu(i)*dens(i)
       polyk = pr(i)/dens(i)**gamma
       !!print*,'polyk = ',polyk
       if (imhd.ne.0) then
          bfield(:,i) = 0.
          divb(i) = 0.
          curlb(:,i) = 0.
       endif
    endif
 enddo
 npart = i-1
 ntotal = npart
 print*,' npart = ',npart
 if (ierr > 0) then
    print*,'ERROR reading '//trim(setupfile)
    stop
 endif
!
!--close the file
!
 close(unit=ireadf)
 write(iprint,*) 'Finished reading setup file: everything is AOK'
!
!--set bounds of setup
!
 xmax = 0.   ! irrelevant
 xmin = 0.
 ibound = 0	! no boundaries 
 iexternal_force = 2	! use central potential
!
!--now change things according to the specific setup required
!
 vel = 0.
 where (rho > 0.)
    hh = hfact*(pmass(i)/rho(:))**dndim
 end where
 
 call primitive2conservative
 !
 !--balance pressure forces with centripedal acceleration
 !
 
 do i=1,npart
    rpart = sqrt(dot_product(x(1:2,i),x(1:2,i)))
    rhat(1:2) = x(1:2,i)/rpart
    rhat(3) = 0.
    frad = abs(dot_product(force(:,i),rhat(:)))
    omegai = sqrt(frad/rpart)
    vel(1,i) = -omegai*x(2,i)/sqrt(rpart)
    vel(2,i) = omegai*x(1,i)/sqrt(rpart)
    if (ndimV.ge.3) vel(3,i) = 0.
    !!print*,'forcez = ',force(3,i),force(2,i)
    !!v2onr = dot_product(vel(1:2,i),vel(1:2,i))/rpart
    !!print*,'v2onr = ',v2onr*rhat(2),force(2,i)
 enddo
 
! zz(:) = 0.
! zz(3) = 1.
! do i=1,npart
!    rpart = sqrt(dot_product(x(:,i),x(:,i)))
!    rhat(:) = x(:,i)/rpart
!    frad = abs(dot_product(force(:,i),rhat(:)))
!    omegai = sqrt(frad/rpart)
!    rxy = sqrt(dot_product(x(1:2,i),x(1:2,i)))
!    theta = ACOS(rxy/rpart)
!    !--compute zhat (unit vector perpendicular to spherical r)
!    zhat(:) = zz(:) - SIN(theta)*rhat(:)
!    !--vel = omega*zhat cross r
!    vel(1,i) = omegai*(zhat(2)*x(3,i) - zhat(3)*x(2,i))
!    vel(2,i) = omegai*(zhat(3)*x(1,i) - zhat(1)*x(3,i))
!    vel(3,i) = omegai*(zhat(1)*x(2,i) - zhat(2)*x(1,i))
! enddo

 return

666 stop 'error opening setup file'

END SUBROUTINE setup

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
