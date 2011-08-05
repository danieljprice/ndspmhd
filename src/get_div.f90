!!------------------------------------------------------------------------
!! Computes an SPH estimate of div B given a vector B
!!
!! This version computes div B on all particles
!! and therefore only does each pairwise interaction once
!!
!! particle quantities are explicitly passed rather than through the 
!! global modules
!!
!! assumes there has been a previous call to density
!!------------------------------------------------------------------------

SUBROUTINE get_divBgradpsi(divB,Bin,x,hh,pmass,rho,npart,ntot)
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE bound
 USE kernel
 USE linklist
 USE options
 USE get_neighbour_lists
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: npart,ntot
 REAL, DIMENSION(ndim,ntot), INTENT(IN) :: x
 REAL, DIMENSION(ntot), INTENT(IN) :: hh,pmass,rho
 REAL, DIMENSION(ndimV,ntot), INTENT(IN) :: Bin
 REAL, DIMENSION(ntot), INTENT(OUT) :: divB
!
!--define local variables
!
 INTEGER :: i,j,n
 INTEGER :: icell,iprev,nneigh
 INTEGER, DIMENSION(ntot) :: listneigh ! neighbour list
 INTEGER :: idone
!
!  (particle properties - local copies)
!      
 REAL :: rij,rij2
 REAL :: hi,hi1,hav,hav1,hj,hj1,h2,hi2,hj2
 REAL :: hfacwab,hfacwabi,hfacwabj
 REAL :: pmassi,pmassj,projdB
 REAL, DIMENSION(ndim) :: dx
 REAL, DIMENSION(ndimV) :: dr
!
!  (kernel quantities)
!
 REAL :: q2,q2i,q2j      
 REAL :: wab,wabi,wabj
 REAL :: grkern,grkerni,grkernj,weight
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine get_divBgradpsi, npart,ntot=',npart,ntot 
!
!--initialise quantities
!
 listneigh = 0
 divB = 0.
 rho = 0.
!
!--Loop over all the link-list cells
!
 loop_over_cells: DO icell=1,ncellsloop          ! step through all cells
!
!--get the list of neighbours for this cell 
!  (common to all particles in the cell)
!
    CALL get_neighbour_list(icell,listneigh,nneigh)
!
!--now loop over all particles in the current cell
!
    i = ifirstincell(icell)
    idone = -1     ! note density summation includes current particle
    IF (i.NE.-1) iprev = i

    loop_over_cell_particles: DO WHILE (i.NE.-1) ! loop over home cell particles

       !PRINT*,'Doing particle ',i,nneigh,' neighbours',hh(i),Bin(:,i)
       idone = idone + 1
       hi = hh(i)
       hi1 = 1./hi
       hi2 = hi*hi
       hfacwabi = hi1**ndim
       pmassi = pmass(i)
!
!--for each particle in the current cell, loop over its neighbours
!
       loop_over_neighbours: DO n = idone+1,nneigh
          j = listneigh(n)
          IF (.NOT.(j.GT.npart .AND. i.GT.npart)) THEN
             ! do count particle with itself       
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
             pmassj = pmass(j)
             
             rij2 = DOT_PRODUCT(dx,dx)
             rij = SQRT(rij2)
             q2 = rij2/h2
             q2i = rij2/hi2
             q2j = rij2/hj2
             dr = 0.
             IF (j.NE.i) THEN
                dr(1:ndim) = dx(1:ndim)/rij  ! unit vector
             ENDIF
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
                !PRINT*,' neighbour,r/h,dx,hi,hj ',i,j,SQRT(q2),dx,hi,hj
                IF (ikernav.EQ.1) THEN          
                   CALL interpolate_kernel(q2,wab,grkern)
                   wab = wab*hfacwab
                   grkern = grkern*hfacwab*hj1
                   grkerni = grkern
                   grkernj = grkern
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
                   if (ikernav.eq.2) then
                      wabi = wab
                      wabj = wab
                      grkerni = grkern
                      grkernj = grkern
                   endif           
                ENDIF
!
!--calculate div B and grad psi
!
                projdB = DOT_PRODUCT(Bin(:,i)-Bin(:,j),dr)

                divB(i) = divB(i) - pmassj*projdB*grkerni
                divB(j) = divB(j) - pmassi*projdB*grkernj            
                                
             ENDIF

         ENDIF! .not. j>npart .and. i>npart   
       ENDDO loop_over_neighbours
       
       iprev = i
       IF (iprev.NE.-1) i = ll(i)          ! possibly should be only IF (iprev.NE.-1)
    ENDDO loop_over_cell_particles
    
 ENDDO loop_over_cells

!
!--do divisions by rho
!  note that grad psi returns the correct term appropriate to evolving either
!  B or B/rho (ie. grad psi or grad psi / rho respectively). This should be
!  the same in rates.
!
 DO i=1,npart
    if (rho(i).ge.0.) then
       divB(i) = divB(i)/(gradh(i)*rho(i))
    else
       write(*,*) 'ERROR: get_divBgradpsi: rho < 0.'
    endif
    !print*,'divB, rho = ',divB(i),rho(i),Bin(:,i)
    !if (mod(i,20).eq.0) read*
 ENDDO
 
 DO i=npart+1,ntot
    j = ireal(i)
    divB(i) = divB(j)
 ENDDO
! read*

 RETURN
END SUBROUTINE get_divBgradpsi
      
