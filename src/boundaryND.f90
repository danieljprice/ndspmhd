!!-------------------------------------------------------------------------
!! update of particles at end of step for periodic (particles crossing
!! domain) and inflow/outflow (adds / subtracts particles from the domain)
!!
!! for ND case only periodic boundaries implemented
!!-------------------------------------------------------------------------
	 
SUBROUTINE boundary
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE bound
 USE options
 USE part
 USE part_in
! USE setup_params
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i,j,nnew,nsub,jdim
 REAL :: psepleft,psepright
 LOGICAL :: debugging
!
!--allow for tracing flow
!      
 IF (trace) THEN
    WRITE(iprint,*) ' Entering subroutine boundary'
    debugging = .true.
 ENDIF   
 
!---------------------------------------------------------------------------
!  inflow/outflow when fixed particles are used.
!---------------------------------------------------------------------------     
 IF (ANY(ibound.EQ.1)) THEN
    IF (ndim.GT.1) WRITE(iprint,*) ' inflow/outflow not implemented in ND'
 ENDIF
!----------------------------------------------------------------------
! periodic boundary conditions - allow particles to cross the domain
!----------------------------------------------------------------------

 IF (ANY(ibound.EQ.3)) THEN
    
    DO i=1,npart

       DO jdim=1,ndim
          IF (ibound(jdim).EQ.3) THEN
	     IF (x(jdim,i).GT.xmax(jdim)) THEN
!	        print*,'ss xold,xmax,xnew = ',jdim,x(jdim,i),xmax(jdim),xmin(jdim) + x(jdim,i) - xmax(jdim)
	        x(jdim,i) = xmin(jdim) + x(jdim,i) - xmax(jdim)
	        xin(jdim,i) = x(jdim,i)
	     ELSEIF(x(jdim,i).LT.xmin(jdim)) THEN
!	        print*,'ss xold,xmin,xnew = ',jdim,x(jdim,i),xmin(jdim),xmax(jdim) + x(jdim,i) - xmin(jdim)	     
!	        read*
                x(jdim,i) = xmax(jdim) - (xmin(jdim) - x(jdim,i))
	        xin(jdim,i) = x(jdim,i)
             ENDIF	  
	  ENDIF
       ENDDO
    
    ENDDO   

 ENDIF

 RETURN 
END SUBROUTINE boundary
