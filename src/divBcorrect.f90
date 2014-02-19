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

!!--------------------------------------------------------------------------
!! Interface subroutine for the div B correction by the projection method
!!
!!--------------------------------------------------------------------------

SUBROUTINE divBcorrect(npts,ntot)
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE bound
 USE derivB
 USE options
 USE part
 USE timestep

 use getdivB
 use getcurl
 
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: npts, ntot
 REAL, PARAMETER :: pi = 3.14159265358979
 INTEGER :: i,icall
 REAL :: phi,dpi,ecrap,momcrap
 REAL, DIMENSION(npts) :: source
 REAL, DIMENSION(size(rho)) :: divBonrho
 REAL, DIMENSION(ndimV,npts) :: sourcevec,curlA
 REAL, DIMENSION(ndim,npts) :: gradphi
 LOGICAL :: debugging
 SAVE icall
!
!--allow for tracing flow
!      
 debugging = .false.
 IF (trace) WRITE(iprint,*) ' Entering subroutine divBcorrect'
 IF (idebug(1:4).EQ.'divB') debugging = .true.
 
 IF (npts.ne.npart) STOP 'npts should be = npart for memory allocation'
 SELECT CASE(ndim)
  CASE(2)
     dpi = 1./(2.*pi)
  CASE(3)
     dpi = 1./(4.*pi)
  CASE default
     dpi = 1.
     !WRITE(iprint,*) 'divergence correction can''t be done in 1D'
     !RETURN
 END SELECT

 icall = icall + 1

 SELECT CASE(idivBzero)
!----------------------------------------------------
! Standard projection method (scalar potential)
!----------------------------------------------------
    CASE(10)
!
!--get neighbours and calculate density if required
!
       if (debugging) write(iprint,*) ' linking ...'
       if (any(ibound.gt.1)) call set_ghost_particles
       call set_linklist
       if (debugging) write(iprint,*) ' calculating density...'
       if (icty.le.0) call iterate_density
!
!--calculate div B source term for poisson equation
!       
       if (debugging) write(iprint,*) ' calculating div B before correction'
       call get_divB(divBonrho(1:ntot),ntot)
       if (ntot.gt.npart) divBonrho(npart+1:ntot) = 0.
       divB(1:ntot) = rho(1:ntot)*divBonrho(1:ntot)

       if (debugging) then
          call output(time,nsteps)   ! output div B and Bfield before correction
          call evwrite(real(icall),ecrap,momcrap)
       endif
!
!--specify the source term for the Poisson equation
!    
       source(1:npart) = pmass(1:npart)*divBonrho(1:npart)*dpi
!
!--calculate the correction to the magnetic field
! 
       write(iprint,"(a)",ADVANCE='NO') ' div B correction by scalar projection step...'
       CALL direct_sum_poisson(x(:,1:npart),source,phi,gradphi,npart)
!
!--correct the magnetic field
!              
       Bfield(1:ndim,1:npart) = Bfield(1:ndim,1:npart) - gradphi(1:ndim,1:npart)
       
       write(iprint,*) 'done'
       
       IF (imhd.GE.11) THEN
          Bevol(1:ndim,1:npart) = Bfield(1:ndim,1:npart)
       ELSE
          DO i=1,npart
             Bevol(1:ndim,i) = Bfield(1:ndim,i)/rho(i)
          ENDDO
       ENDIF
       if (ntot.gt.npart) then
          do i=npart+1,ntot
             call copy_particle(i,ireal(i))
          enddo
       endif
       !!CALL primitive2conservative ! so Bfield -> Bevol
       if (debugging) then
          write(iprint,*) ' calculating div B after correction'
          CALL get_divB(divBonrho,ntot)
          !
          !--set divB to zero on ghosts and on real counterparts
          !  this is to avoid problems at the boundary
          !
          do i=npart+1,ntot
             divBonrho(i) = 0.
             divBonrho(ireal(i)) = 0.
          enddo
          divB(1:ntot) = rho(1:ntot)*divBonrho(1:ntot)
          call output(time,nsteps)   ! output div B and Bfield after correction
          call evwrite(real(icall),ecrap,momcrap)
       endif

!----------------------------------------------------
! Current projection method (vector potential)
!----------------------------------------------------       
    CASE(11)
!
!--get neighbours and calculate density if required
!
       if (debugging) write(iprint,*) ' linking ...'
       if (any(ibound.gt.1)) call set_ghost_particles
       call set_linklist
       if (debugging) write(iprint,*) ' calculating density...'
       if (icty.le.0) call iterate_density
!
!--calculate div B source term for poisson equation
!       
       if (debugging) write(iprint,*) ' calculating current before correction'
       call get_curl(1,ntot,x,pmass,rho,hh,Bfield,curlB)
       if (ntot.gt.npart) curlB(:,npart+1:ntot) = 0.
       
       if (debugging) then
          write(iprint,*) ' calculating div B before correction'
          call get_divB(divBonrho,ntot)
          if (ntot.gt.npart) divBonrho(npart+1:ntot) = 0.
          divB(1:ntot) = rho(1:ntot)*divBonrho(1:ntot)

          call output(time,nsteps)   ! output curl B and Bfield before correction
          call evwrite(real(icall),ecrap,momcrap)
       endif
!
!--specify the source term for the Poisson equation
!    
       sourcevec = 0.
       do i=1,npart
          sourcevec(:,i) = pmass(i)*curlB(:,i)/rho(i)*dpi
       enddo
!
!--calculate the correction to the magnetic field
! 
       write(iprint,"(a)",ADVANCE='NO') ' div B correction by vector projection step...'
       CALL direct_sum_poisson_vec(x(:,1:npart),sourcevec,curlA,npts)
!
!--correct the magnetic field
!              
       !!print*,'Bfield = ',Bfield(:,3436),curlA(:,3436)
       Bfield(1:ndim,1:npart) = curlA(1:ndim,1:npart)
       
       write(iprint,*) 'done'

       IF (imhd.GE.11) THEN
          Bevol(1:ndim,1:npart) = Bfield(1:ndim,1:npart)
       ELSE
          DO i=1,npart
             Bevol(1:ndim,i) = Bfield(1:ndim,i)/rho(i)
          ENDDO
       ENDIF
       if (ntot.gt.npart) then
          !
          !--overwrite B near boundary
          !
          do i=npart+1,ntot
             call copy_particle(ireal(i),i)
          enddo
       endif
       
       if (debugging) then
          write(iprint,*) ' calculating div/curl B after correction'
          call get_divB(divBonrho,ntot)
          !
          !--set divB to zero on ghosts and on real counterparts
          !  this is to avoid problems at the boundary
          !
          do i=npart+1,ntot
             divBonrho(i) = 0.
             divBonrho(ireal(i)) = 0.
          enddo
          divB(1:ntot) = rho(1:ntot)*divBonrho(1:ntot)
          
          call get_curl(1,ntot,x,pmass,rho,hh,Bfield,curlB)
          if (ntot.gt.npart) curlB(:,npart+1:ntot) = 0.
          
          call output(time,nsteps)   ! output div B and Bfield after correction
          call evwrite(real(icall),ecrap,momcrap)
       endif      
       
 END SELECT
 
 RETURN
END SUBROUTINE divBcorrect
