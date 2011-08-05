!!--------------------------------------------------------------------------
!! Interface subroutine for the div B correction by the projection method
!!
!!--------------------------------------------------------------------------

SUBROUTINE divBcorrect(npts,ntot)
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE part
 USE derivB
 USE options
 USE timestep
 IMPLICIT NONE
 REAL, PARAMETER :: pi = 3.14159265358979
 INTEGER :: ntot,npts
 REAL :: phi,dpi
 REAL, DIMENSION(npts) :: source
 REAL, DIMENSION(ntot) :: divBonrho
 REAL, DIMENSION(ndim,npts) :: gradphi
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine divBcorrect'
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

 SELECT CASE(idivBzero)
 
    CASE(10)
       print*,' linking ...'
       CALL link
       print*,'calculating density...'
       CALL density
       print*,' calculating div B before correction'
       CALL get_divB(divBonrho,ntot)
       if (ntot.gt.npart) divBonrho(npart+1:ntot) = 0.
       divB(1:ntot) = rho(1:ntot)*divBonrho(1:ntot)

       CALL output(time,nsteps)   ! output div B and Bfield before correction
!
!--specify the source term for the Poisson equation
!    
       source(1:npart) = pmass(1:npart)*divBonrho(1:npart)*dpi
!
!--calculate the correction to the magnetic field
! 
       print*,' calculating correction to magnetic field...'
       CALL direct_sum_poisson(x(:,1:npart),source,phi,gradphi,npart)
!
!--correct the magnetic field
!              
       Bfield(1:ndim,1:npart) = Bfield(1:ndim,1:npart) - gradphi(1:ndim,1:npart)
       print*,' calculating div B after correction'
       CALL get_divB(divBonrho,ntot)
       if (ntot.gt.npart) divBonrho(npart+1:ntot) = 0.
       divB(1:ntot) = rho(1:ntot)*divBonrho(1:ntot)
       CALL primitive2conservative ! so Bfield -> Bcons
       CALL output(time,nsteps)   ! output div B and Bfield after correction
       
    CASE(2)
!       CALL direct_sum_current
       
 END SELECT
 
 RETURN
END SUBROUTINE divBcorrect
