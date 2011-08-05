!!------------------------------------------------------------------------
!! Computes an SPH estimate of the curl of a vector quantity
!! This version computes the curl for all particles
!! and therefore only does each pairwise interaction once
!!------------------------------------------------------------------------

SUBROUTINE get_curl(curlBonrho,ntot)
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE bound
 USE kernel
 USE linklist
 USE options
 USE part
 USE setup_params
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: ntot
 INTEGER :: i,j,n
 INTEGER :: icell,iprev,nneigh
 INTEGER, DIMENSION(ntot) :: listneigh ! neighbour list
 INTEGER :: idone
 INTEGER, DIMENSION(3**ndim) :: neighcell
!
!  (particle properties - local copies)
!      
 REAL :: rij,rij2
 REAL :: hi,hi1,hav,hav1,hj,hj1,h2,hi2,hj2
 REAL :: hfacwab,hfacwabi,hfacwabj
 REAL :: rho21i, rho21j, pmassi
 REAL, DIMENSION(ndim) :: dx
 REAL, DIMENSION(ndimV) :: dr, dB, Bi, Bj, curlBi
 REAL, DIMENSION(ndimV,ntot), INTENT(OUT) :: curlBonrho
!
!  (kernel quantities)
!
 REAL :: q2,q2i,q2j      
 REAL :: wab,wabi,wabj
 REAL :: grkern,grkerni,grkernj
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine get_curl, ntot=',ntot  
!
!--initialise quantities
!
 listneigh = 0
 curlBonrho = 0.
!
!--Loop over all the link-list cells
!
 loop_over_cells: DO icell=1,ncellsloop ! step through all cells
!
!--get the list of neighbours for this cell 
!  (common to all particles in the cell)
!
    CALL get_neighbour_list(icell,neighcell,listneigh,nneigh)
!
!--now loop over all particles in the current cell
!
    i = ifirstincell(icell)
    idone = -1     ! note density summation includes current particle
    IF (i.NE.-1) iprev = i

    loop_over_cell_particles: DO WHILE (i.NE.-1)          ! loop over home cell particles

       !       PRINT*,'Doing particle ',i,nneigh,' neighbours',hh(i)
       idone = idone + 1
       hi = hh(i)
       hi1 = 1./hi
       hi2 = hi*hi
       hfacwabi = hi1**ndim
       rho21i = 1./rho(i)**2
       pmassi = pmass(i)
       if (imhd.ge.11) then      ! if mag field variable is B
          Bi(:) = Bevol(:,i)
       elseif (imhd.gt.0) then      ! if mag field variable is B/rho
          Bi(:) = Bfield(:,i)
       endif
!
!--for each particle in the current cell, loop over its neighbours
!
       loop_over_neighbours: DO n = idone+1,nneigh
          j = listneigh(n)
          IF ((j.NE.i).AND..NOT.(j.GT.npart .AND. i.GT.npart)) THEN
             ! don't count particle with itself       
             dx(:) = x(:,i) - x(:,j)
             hj = hh(j)
             hj1 = 1./hj
             hj2 = hj*hj
!
!--calculate averages of smoothing length if using this averaging
!                
             hav = 0.5*(hi + hj)
             hav1 = 1./hav
             h2 = hav*hav
             hfacwab = hav1**ndim
             hfacwabj = hj1**ndim
             rho21j = 1./rho(j)**2
             
             rij2 = DOT_PRODUCT(dx,dx)
             rij = SQRT(rij2)
             q2 = rij2/h2
             q2i = rij2/hi2
             q2j = rij2/hj2     
             dr(1:ndim) = dx(1:ndim)/rij  ! unit vector
             if (ndimV.gt.ndim) dr(ndim+1:ndimV) = 0. 

             !          PRINT*,' neighbour,r/h,dx,hi,hj ',j,SQRT(q2),dx,hi,hj
!     
!--do interaction if r/h < compact support size
!  don't calculate interactions between ghost particles
!
             IF (((q2i.LT.radkern2).OR.(q2j.LT.radkern2))  &
                  .AND. .NOT.(i.GT.npart.AND.j.GT.npart)) THEN
!     
!--interpolate from kernel table          
!  (use either average h or average kernel gradient)
!
                !       PRINT*,' neighbour,r/h,dx,hi,hj ',i,j,SQRT(q2),dx,hi,hj
                IF (ikernav.EQ.1) THEN          
                   CALL interpolate_kernel(q2,wab,grkern)
                   wab = wab*hfacwab
                   grkern = grkern*hfacwab*hj1
                ELSE
                   !  (using hi)
                   CALL interpolate_kernel(q2i,wabi,grkerni)
                   wabi = wabi*hfacwabi
                   grkerni = grkerni*hfacwabi*hi1
                   !  (using hj)
                   CALL interpolate_kernel(q2j,wabj,grkernj)
                   wabj = wabj*hfacwabj
                   grkernj = grkernj*hfacwabj*hj1
                   !  (calculate average)            
                   wab = 0.5*(wabi + wabj)                  
                   grkern = 0.5*(grkerni + grkernj)            
                ENDIF

                if (imhd.ge.11) then      ! if B is mag field variable
                   Bj(:) = Bevol(:,j)
                elseif (imhd.ne.0) then      ! if B/rho is mag field variable
                   Bj(:) = Bfield(:,j)                                
                endif
!
!--calculate div B
!
                dB = Bi(:) - Bj(:)
                if (ndimV.eq.3) then
                   curlBi(1) = dB(2)*dr(3) - dB(3)*dr(2)
                   curlBi(2) = dB(3)*dr(1) - dB(1)*dr(3)
                   curlBi(3) = dB(1)*dr(2) - dB(2)*dr(1)
                elseif (ndimV.eq.2) then  ! just Jz in 2D
                   curlBi(1) = dB(1)*dr(2) - dB(2)*dr(1)
                   curlBi(2) = 0.
                endif
                !
                !--compute rho * current density J
                !
                curlBonrho(:,i) = curlBonrho(:,i) - pmass(j)*curlBi(:)*grkern
                curlBonrho(:,j) = curlBonrho(:,j) - pmassi*curlBi(:)*grkern
                !      ELSE
                !         PRINT*,' r/h > 2 '      
                
             ENDIF
          ENDIF! j .ne. i   
       ENDDO loop_over_neighbours
       
       iprev = i
       IF (iprev.NE.-1) i = ll(i)          ! possibly should be only IF (iprev.NE.-1)
    ENDDO loop_over_cell_particles
    
 ENDDO loop_over_cells

 do i=1,ntot
    curlBonrho(:,i) = curlBonrho(:,i)/rho(i)**2
 enddo

 RETURN
END SUBROUTINE get_curl
      
