!!--------------------------------------------------------------------------
!! Interface subroutine for the div B correction by the projection method
!!
!!--------------------------------------------------------------------------

SUBROUTINE divBcorrect(npts,ntot)
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE bound
 USE part
 USE derivB
 USE options
 USE timestep
 IMPLICIT NONE
 REAL, PARAMETER :: pi = 3.14159265358979
 INTEGER :: i,ntot,npts,icall
 REAL :: phi,dpi,ecrap,momcrap
 REAL, DIMENSION(npts) :: source
 REAL, DIMENSION(ntot) :: divBonrho
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
     WRITE(iprint,*) 'divergence correction can''t be done in 1D'
     RETURN
 END SELECT

 icall = icall + 1

 SELECT CASE(idivBzero)
 
    CASE(10)
!
!--get neighbours and calculate density if required
!
       if (debugging) write(iprint,*) ' linking ...'
       call set_linklist
       if (debugging) write(iprint,*) ' calculating density...'
       if (icty.le.0) call iterate_density
!
!--calculate div B source term for poisson equation
!       
       if (debugging) write(iprint,*) ' calculating div B before correction'
       call get_divB(divBonrho,ntot)
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
       write(iprint,"(a)",ADVANCE='NO') ' div B correction by projection step...'
       CALL direct_sum_poisson(x(:,1:npart),source,phi,gradphi,npart)
!
!--correct the magnetic field
!              
       Bfield(1:ndim,1:npart) = Bfield(1:ndim,1:npart) - gradphi(1:ndim,1:npart)
       
       write(iprint,*) 'done'
       
       if (debugging) then
          write(iprint,*) ' calculating div B after correction'
          CALL get_divB(divBonrho,ntot)
	  if (ntot.gt.npart) divBonrho(npart+1:ntot) = 0.
          divB(1:ntot) = rho(1:ntot)*divBonrho(1:ntot)
       endif
       
       IF (imhd.GE.11) THEN
          Bcons(1:ndim,1:npart) = Bfield(1:ndim,1:npart)
       ELSE
          DO i=1,npart
             Bcons(1:ndim,i) = Bfield(1:ndim,i)/rho(i)
	  ENDDO
       ENDIF
       if (ntot.gt.npart) then
	  do i=npart+1,ntot
	     call copy_particle(i,ireal(i))
          enddo
       endif
       !!CALL primitive2conservative ! so Bfield -> Bcons
       
       if (debugging) then
          call output(time,nsteps)   ! output div B and Bfield after correction
          call evwrite(real(icall),ecrap,momcrap)
       endif
       
    CASE(11)
!       CALL direct_sum_current
       
 END SELECT
 
 RETURN
END SUBROUTINE divBcorrect
