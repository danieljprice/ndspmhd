!!------------------------------------------------------------------------
!! Computes a smoothed estimate of the magnetic field
!! by summation over the particles neighbours
!! ie. B_a = sum_b m_b B_b/rho_b W_ab
!! also computes the unity function
!! B_a = sum_b m_b /rho_b Wab
!! requires rho to be previously known
!!------------------------------------------------------------------------

SUBROUTINE get_smoothB
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE bound
 USE kernel
 USE linklist
 USE options
 USE part
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i,j,n
 INTEGER :: icell,ipart,iprev,ncell,nneigh,index,index1
 INTEGER :: listneigh(nlistdim) ! up to 10% of particles in each cell
 INTEGER :: neighcell(3),idone
!
!  (particle properties - local copies)
!      
 REAL :: dx,rij,rij2,dvx
 REAL :: hi,h2,h3,hav,hj,hi2,hj2,hi3,hj3
 REAL, ALLOCATABLE, DIMENSION(:) :: unity
!
!  (kernel quantities)
!
 REAL :: q2,q2i,q2j      
 REAL :: wab,wabi,wabj,dwdx,dxx,weight      
 REAL :: grkerni,grkernj,dgrwdx
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine get_smoothB'  
!
!--initialise quantities
!
 ALLOCATE(unity(ntotal))

 DO i=1,ntotal
    Bsmooth(:,i) = 0.
    unity(i) = 0.
 ENDDO
!
!--Loop over all the link-list cells
!
 DO icell=1,ncells		! step through all cells
    i = ifirstincell(icell)
!
!--get the list of neighbours for this cell 
!  (common to all particles in the cell)
!
    CALL get_neighbour_list(icell,neighcell,listneigh,nneigh)
!
!--now loop over all particles in the current cell
!
    idone = -1	! note density summation includes current particle
    iprev = i
    DO WHILE (ll(iprev).NE.-1)		! loop over home cell particles
!      PRINT*,'Doing particle ',i
       idone = idone + 1
       hi = hh(i)
       hi2 = hi*hi
       hi3 = hi*hi2		    
!
!--for each particle in the current cell, loop over its neighbours
!
       DO n = idone+1,nneigh
	  j = listneigh(n)
!         print*,' ... neighbour, h=',j,0.5*(hi + hh(j))
	  dx = x(i) - x(j)		
	  hj = hh(j)
	  hj2 = hj*hj
	  hj3 = hj*hj2
!
!--calculate averages of smoothing length if using this averaging
!			 
	  hav = 0.5*(hi + hj)
	  h2 = hav*hav
		
	  rij2 = dx*dx
	  rij = SQRT(rij2)
	  q2 = rij2/h2
	  q2i = rij2/hi2
	  q2j = rij2/hj2	
!	
!--interpolate from kernel table		
!
	 IF ((q2.LT.radkern2).OR.(q2i.LT.radkern2).OR.(q2j.LT.radkern2)) THEN
!		
!--use either average h or average kernel gradient
!
	    IF (ikernav.EQ.1) THEN		
	       index = INT(q2/dq2table)	! nearest index in table
	       index1 = index + 1
	       IF (index.GT.ikern) index = ikern
	       IF (index1.GT.ikern) index1 = ikern
	       dxx = q2 - index*dq2table 	! increment along from index pt		
	       dwdx =  (wij(index1)-wij(index))/dq2table ! slope
	       wab = (wij(index)+ dwdx*dxx)/hav	! divide by h in 1D		
            ELSE
!  (using hi)
	       index = INT(q2i/dq2table)	! nearest index in table
	       index1 = index + 1
	       IF (index.GT.ikern) index = ikern
	       IF (index1.GT.ikern) index1 = ikern
	       dxx = q2i - index*dq2table 	! increment along from index pt
	       dwdx =  (wij(index1)-wij(index))/dq2table ! slope
	       wabi = (wij(index)+ dwdx*dxx)/hi	! divide by h in 1D		
	       dgrwdx =  (grwij(index1)-grwij(index))/dq2table ! slope
	       grkerni = (grwij(index)+ dgrwdx*dxx)/hi2	! divide by h**3. in 1D		
!  (using hj)
	       index = INT(q2j/dq2table)	! nearest index in table
	       index1 = index + 1
	       IF (index.GT.ikern) index = ikern
	       IF (index1.GT.ikern) index1 = ikern
	       dxx = q2j - index*dq2table 	! increment along from index pt
	       dwdx =  (wij(index1)-wij(index))/dq2table ! slope
	       wabj = (wij(index)+ dwdx*dxx)/hj	! divide by h in 1D		
	       dgrwdx =  (grwij(index1)-grwij(index))/dq2table ! slope
	       grkernj = (grwij(index)+ dgrwdx*dxx)/hj2	! divide by h**3. in 1D		
!  (calculate average)  		
	       wab = 0.5*(wabi + wabj)
	    ENDIF
!
!--calculate quantities
!
	    weight = 1.0
	    IF (j.EQ.i) weight = 0.5
	    IF (ikernav.EQ.3) THEN
	       Bsmooth(:,i) = Bsmooth(:,i) 	&
	                    + pmass(j)*Bfield(:,j)/rho(j)*wabi*weight
	       Bsmooth(:,j) = Bsmooth(:,j) 	&
	       	            + pmass(i)*Bfield(:,i)/rho(i)*wabj*weight
	       unity(i) = unity(i) + pmass(j)/rho(j)*wabi*weight
	       unity(j) = unity(j) + pmass(i)/rho(i)*wabj*weight	       		 
	    ELSE
	       Bsmooth(:,i) = Bsmooth(:,i) 	&
	                    + pmass(j)*Bfield(:,j)/rho(j)*wab*weight
	       Bsmooth(:,j) = Bsmooth(:,j) 	&
	       	            + pmass(i)*Bfield(:,i)/rho(i)*wab*weight
	       unity(i) = unity(i) + pmass(j)/rho(j)*wab*weight
	       unity(j) = unity(j) + pmass(i)/rho(i)*wab*weight	       		 
	    ENDIF
	       
	 ENDIF
	    
       ENDDO	! over neighbours
	    
       iprev = i
       IF (iprev.NE.-1) i = ll(i)
    ENDDO	! over home cell particles
            
 ENDDO	! over cells

 DO i=1,npart
    Bsmooth(:,i) = Bsmooth(:,i)/unity(i)
 ENDDO 

 RETURN
END SUBROUTINE get_smoothB
      
