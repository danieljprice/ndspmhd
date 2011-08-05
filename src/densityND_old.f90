!!------------------------------------------------------------------------
!! Computes the density by direct summation over the particles neighbours
!! ie. rho_a = sum_b m_b W_ab
!! *** grad h terms not right in ND
!!------------------------------------------------------------------------

SUBROUTINE density
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE bound
 USE hterms
 USE kernel
 USE linklist
 USE options
 USE part
 USE setup_params
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i,j,n
 INTEGER :: icell,icellloop,ipart,iprev,ncell,nneigh,index,index1
 INTEGER :: listneigh(nlistdim) ! up to 10% of particles in each cell
 INTEGER :: idone
 INTEGER, DIMENSION(3*2**(ndim-1) - 1) :: neighcell
!
!  (particle properties - local copies)
!      
 REAL :: rij,rij2,rho1i
 REAL :: hi,hi1,h2,h3,hav,hj,hj1,hi2,hj2,hi3,hj3
 REAL :: hfacwab,hfacwabi,hfacwabj,hfacgrkerni,hfacgrkernj
 REAL, DIMENSION(:), ALLOCATABLE :: rho_old,hh_old
 REAL, DIMENSION(ndim) :: dx
!
!  (kernel quantities)
!
 REAL :: q2,q2i,q2j      
 REAL :: wab,wabi,wabj,dwdx,dxx,weight      
 REAL :: grkerni,grkernj,dgrwdx
 REAL :: ddq2table
!
!  (grad h terms)
!
 INTEGER :: itsdensity,itsdensitymax
 REAL :: dwdhi,dwdhj,rhoerr,vol,tol	!,dhdrhoi,dhdrhoj

!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine density'  
!
!--allocate memory for local array copies
!
 ALLOCATE( rho_old(ntotal), hh_old(ntotal) )
 
 IF ((ikernav.EQ.3).AND.(ihvar.NE.0)) THEN
    itsdensitymax = 5	! perform 1 fixed point iteration
 ELSE
    itsdensitymax = 0	! no iterations
 ENDIF
!
!--Loop to find rho and h self-consistently (if using Springel/Hernquist)
!
 itsdensity = 0
 rhoerr = 1.e6
 tol = 1.e-3
 iterate_density: DO WHILE ((rhoerr.GT.tol).AND.(itsdensity.LE.itsdensitymax))

 itsdensity = itsdensity + 1
!
!--initialise quantities
!
 DO i=1,npart
    rho_old(i) = rho(i)	! don't need these
    hh_old(i) = hh(i)	! " " just to show convergence
    rho(i) = 0.
    gradh(i) = 0.
 ENDDO
!
!--Loop over all the link-list cells
!
 loop_over_cells: DO icell=1,ncellsloop		! step through all cells
!
!--get the list of neighbours for this cell 
!  (common to all particles in the cell)
!
    CALL get_neighbour_list(icell,neighcell,listneigh,nneigh)
!
!--now loop over all particles in the current cell
!
    i = ifirstincell(icell)
    idone = -1	! note density summation includes current particle
    IF (i.NE.-1) iprev = i

    loop_over_cell_particles: DO WHILE (i.NE.-1)		! loop over home cell particles

!       PRINT*,'Doing particle ',i,nneigh,' neighbours',pmass(i)
       idone = idone + 1
       hi = hh(i)
       hi2 = hi*hi
       hi3 = hi*hi2
       hi1 = 1./hi
       
       hfacwabi = hi1**ndim
       hfacgrkerni = hi1**(ndim+1)    
!
!--for each particle in the current cell, loop over its neighbours
!
       loop_over_neighbours: DO n = idone+1,nneigh
	  j = listneigh(n)
	  dx(:) = x(:,i) - x(:,j)
	  hj = hh(j)
	  hj2 = hj*hj
	  hj3 = hj*hj2
	  hj1 = 1./hj
!
!--calculate averages of smoothing length if using this averaging
!			 
	  hav = 0.5*(hi + hj)
	  h2 = hav*hav
	  hfacwab = 1./hav**ndim
	  hfacwabj = hj1**ndim
	  hfacgrkernj = hj1**(ndim+1)
	  
	  rij2 = DOT_PRODUCT(dx,dx)
	  rij = SQRT(rij2)
	  q2 = rij2/h2
	  q2i = rij2/hi2
	  q2j = rij2/hj2	
!          PRINT*,' neighbour,r/h,dx,hi,hj ',j,SQRT(q2),dx,hi,hj
!	
!--do interaction if r/h < compact support size
!  don't calculate interactions between ghost particles
!
	 IF (((q2.LT.radkern2).OR.(q2i.LT.radkern2).OR.(q2j.LT.radkern2))  &
	     .AND. .NOT.(i.GT.npart.AND.j.GT.npart)) THEN
!	
!--interpolate from kernel table		
!  (use either average h or average kernel gradient)
!
!          IF (i.EQ.100.OR.j.EQ.100) THEN
!	  PRINT*,' neighbour,r/h,dx,hi,hj ',i,j,SQRT(q2),dx,hi,hj
!	  ENDIF
            ddq2table = 1./dq2table
	    IF (ikernav.EQ.1) THEN		
	       index = INT(q2*ddq2table)	! nearest index in table
	       index1 = index + 1
	       IF (index.GT.ikern) index = ikern
	       IF (index1.GT.ikern) index1 = ikern
	       dxx = q2 - index*dq2table 	! increment along from index pt		
	       dwdx =  (wij(index1)-wij(index))*ddq2table ! slope
	       wab = (wij(index)+ dwdx*dxx)*hfacwab	! divide by h in 1D		
            ELSE
!  (using hi)
	       index = INT(q2i*ddq2table)	! nearest index in table
	       index1 = index + 1
	       IF (index.GT.ikern) index = ikern
	       IF (index1.GT.ikern) index1 = ikern
	       dxx = q2i - index*dq2table 	! increment along from index pt
	       dwdx =  (wij(index1)-wij(index))*ddq2table ! slope
	       wabi = (wij(index)+ dwdx*dxx)*hfacwabi	! divide by h in 1D		
	       dgrwdx =  (grwij(index1)-grwij(index))*ddq2table ! slope
	       grkerni = (grwij(index)+ dgrwdx*dxx)*hfacgrkerni	! divide by h**3. in 1D		
!  (using hj)
	       index = INT(q2j*ddq2table)	! nearest index in table
	       index1 = index + 1
	       IF (index.GT.ikern) index = ikern
	       IF (index1.GT.ikern) index1 = ikern
	       dxx = q2j - index*dq2table 	! increment along from index pt
	       dwdx =  (wij(index1)-wij(index))*ddq2table ! slope
	       wabj = (wij(index)+ dwdx*dxx)*hfacwabj	! divide by h in 1D		
	       dgrwdx =  (grwij(index1)-grwij(index))*ddq2table ! slope
	       grkernj = (grwij(index)+ dgrwdx*dxx)*hfacgrkernj	! divide by h**3. in 1D		
!  (calculate average)  		
	       wab = 0.5*(wabi + wabj)
!
!--derivative w.r.t. h for grad h correction terms (and dhdrho)
!	       
	       dwdhi = -rij*grkerni*hi1 - wabi*hi1
!               dhdrhoi = -hi	! divide by  should use rho(i)
	       dwdhj = -rij*grkernj*hj1 - wabj*hj1
!	        dhdrhoj = -hj/rho_old(j)	! see above
		   				  
	    ENDIF
!
!--calculate density
!
	    weight = 1.0
	    IF (j.EQ.i) weight = 0.5
	    IF (ikernav.EQ.3) THEN
	       rho(i) = rho(i) + pmass(j)*wabi*weight
	       rho(j) = rho(j) + pmass(i)*wabj*weight
!
!--correction term for variable smoothing lengths
!  this is the small bit that should be 1-gradh
!  need to divide by rho once rho is known

	       gradh(i) = gradh(i) - ndim*hi*pmass(j)*weight*dwdhi
	       gradh(j) = gradh(j) - ndim*hj*pmass(i)*weight*dwdhj
		  
	    ELSE
	       rho(i) = rho(i) + pmass(j)*wab*weight		
	       rho(j) = rho(j) + pmass(i)*wab*weight	
	    ENDIF
!	 ELSE
!	    PRINT*,' r/h > 2 '      
	 
	 ENDIF
	    
       ENDDO loop_over_neighbours
	    
       iprev = i
       IF (iprev.NE.-1) i = ll(i)		! possibly should be only IF (iprev.NE.-1)
    ENDDO loop_over_cell_particles
            
 ENDDO loop_over_cells

! PRINT*,'finished density, iteration ',itsdensity
! PRINT*,' rho(100) = ',rho(100)
! READ* 

 rhoerr = 0.
 vol = 0.
 IF (ihvar.NE.0) THEN
    DO i=nbpts+1,npart-nbpts		! reset h to new value of density
       IF (rho(i).LE.0.) THEN
	  WRITE(iprint,*) 'rho: rho -ve, dying rho(',i,') = ',rho(i)
	  CALL quit
       ENDIF

       gradh(i) = gradh(i)/rho_old(i)    ! now that rho is known
       IF (abs(1.-gradh(i)).LT.1.e-5) PRINT*,'eek'
!
!--perform Newton-Raphson iteration
!      
       rho(i) = rho_old(i) - (rho_old(i) - rho(i))/(1.-gradh(i))
!
!       IF (ikernav.EQ.3) THEN
       rho1i = 1./rho(i)
       rhoerr = rhoerr + pmass(i)*abs(rho(i) - rho_old(i))*rho1i
       vol = vol + pmass(i)
!       gradh(i) = gradh(i)*rho1i    ! now that rho is known
       hh(i) = hfact*(pmass(i)*rho1i)**hpower	! ie h proportional to 1/rho^dimen
!       ENDIF
    ENDDO

    rhoerr = rhoerr/vol
    IF ((rhoerr.LT.tol).AND.(itsdensity.GT.1))    &
       WRITE(iprint,*) 'finished density, iteration ',itsdensity,' error = ',rhoerr

 ENDIF

!
!--update ghost particles (copy values from real particles)
!
 DO i=npart+1,ntotal
    j = ireal(i)
    rho(i) = rho(j)
    hh(i) = hh(j)
    gradh(i) = gradh(j)	! same
 ENDDO 
      
!      IF (itsdensity.GT.1) THEN
!         PRINT*,'rho_old(200) ,rho(200),hh_old(200),hh(200)=',
!     &          rho_old(100),rho(100),hh_old(100),hh(100)
!      ENDIF
      
 ENDDO iterate_density

 RETURN
END SUBROUTINE density
      
