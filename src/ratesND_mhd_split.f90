!!--------------------------------------------------------------------
!! Computes the rates of change of everything (forces, energy eqn etc)
!! This is the core of the SPH algorithm
!!--------------------------------------------------------------------

SUBROUTINE get_rates
! USE dimen_mhd
 USE debug
 USE loguns
 
 USE artvi
 USE eos
 USE hterms
 USE kernel
 USE linklist
 USE options
 USE part
 USE rates
 USE timestep
 USE xsph
 USE anticlumping
 USE setup_params	! for hfact

 USE derivB
 USE fmagarray
!
!--variables shared only between rates and subroutines
!
 USE rates_local
 USE rates_mhdlocal
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i,j,n
 INTEGER :: icell,ipart,iprev,ncell,nneigh,index,index1
 INTEGER :: listneigh(nlistdim) 
 INTEGER :: idone
 INTEGER, DIMENSION(3*2**(ndim-1) -1) :: neighcell
!
!  (kernel-related quantities - wab,hfacwab,grkern's shared)
!
 REAL :: hi,hj,hav,hi2,hj2,h2
 REAL :: hfacwabi,hfacwabj,hfacgrkern,hfacgrkerni,hfacgrkernj
 REAL, DIMENSION(ndim) :: dx
 REAL :: q2,q2i,q2j
 REAL :: dgrwdx,dxx
 REAL :: wabi,wabj,dwdx
!
!  (joe's mhd fix - Rjoe,wabjoei shared)
!
 INTEGER :: indexjoe,indexjoe1
 REAL :: wabjoe,dwdxjoe,dxxjoe,wabjoej,deltap2,q2joe
!
! (other particle related)
!
 REAL :: rhoi5,rhoj5,rho1i,rho1j
 REAL :: valfveni,valfvenj,valfven2i,valfven2j
!
! (my av switch)
!
 REAL :: qab,vdotri,vdotrj
 REAL, ALLOCATABLE, DIMENSION(:) :: avsource
!
!  (artificial viscosity quantities - vsig,viss shared)
!      
 REAL :: vsigi,vsigj,vsigav,vsig2i,vsig2j,vsigproji,vsigprojj
 REAL :: spsoundi,spsoundj
 REAL :: source,tdecay    
!
!  (time step criteria)
!      
 REAL :: vsigdtc,vmag
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine get_rates'
!
!--initialise quantities
!      
 dtcourant = 1.e6   
 ALLOCATE (avsource(ntotal))	! should deallocate at end of subroutine

 DO i=1,ntotal	! using ntotal just makes sure they are zero for ghosts
  force(:,i) = 0.0
  drhodt(i) = 0.0
  dudt(i) = 0.0
  dendt(i) = 0.0
  dBfielddt(:,i) = 0.0
  fmag(:,i) = 0.0
  divB(i) = 0.0
  curlB(:,i) = 0.0
  xsphterm(:,i) = 0.0
  avsource(i) = 0.0
 ENDDO
!
!--calculate kernel for Joe's correction term
!
 IF (ianticlump.EQ.1) THEN
  q2joe = (1./hfact)**2	! 1/hfact is initial particle spacing in units of h
!  print*,'q2joe = ',q2joe,q2joe/dq2table
  indexjoe = q2joe/dq2table
  dxxjoe = q2joe - indexjoe*dq2table 
  indexjoe1 = indexjoe + 1
  IF (indexjoe.GT.ikern) indexjoe = ikern
  IF (indexjoe1.GT.ikern) indexjoe1 = ikern
  dwdxjoe = (wij(indexjoe1) - wij(indexjoe))/dq2table
  wabjoei = (wij(indexjoe) + dwdxjoe*dxxjoe) ! need to divide by hav later       
 ENDIF
! print*,'wabjoe = ',wabjoei
 
!
!--Loop over all the link-list cells
!
 loop_over_cells: DO icell=1,ncellsloop		! step through all cells

!    PRINT*,'> doing cell ',icell
!
!--get the list of neighbours for this cell 
!  (common to all particles in the cell)
!
    CALL get_neighbour_list(icell,neighcell,listneigh,nneigh)
	 
    i = ifirstincell(icell)	! start with first particle in cell
    idone = -1
    IF (i.NE.-1) iprev = i

    loop_over_cell_particles: DO WHILE (i.NE.-1)	! loop over home cell particles

!    PRINT*,'Doing particle ',i,x(:,i),valfven2i,rhoi,hh(i)
       idone = idone + 1
       rhoi = rho(i)
       rho1i = 1./rhoi
       rho2i = rhoi*rhoi
       rhoi5 = SQRT(rhoi)
       pri = pr(i)
       Prho2i = pr(i)/rho2i
       spsoundi = spsound(i)
       veli(:) = vel(:,i)
       pmassi = pmass(i)
       alphai = alpha(i)
       IF (imhd.GE.11) THEN	! if mag field variable is B
          Bi(:) = Bfield(:,i)
          Brhoi(:) = Bi(:)/rhoi
       ELSEIF (imhd.GT.0) THEN	! if mag field variable is B/rho
          Brhoi(:) = Bfield(:,i)
          Bi(:) = Brhoi(:)*rhoi
       ENDIF
! mhd definitions
       Brho2i = DOT_PRODUCT(Brhoi,Brhoi)
       valfven2i = Brho2i*rhoi
       valfveni = SQRT(valfven2i)
              
       gradhi = 1. - gradh(i)
       IF (gradhi.LE.0. .OR. gradhi.GT.2) THEN
          WRITE(iprint,*) 'Error in grad h terms, part ',i,gradhi
       ENDIF   
       hi = hh(i)
       hi2 = hi*hi		    
       IF (hi.LE.0.) THEN
          WRITE(iprint,*) ' rates: h <= 0 particle',i,hi
	  CALL quit
       ENDIF
       hfacwabi = 1./hi**ndim
       hfacgrkerni = 1./hi**(ndim+1)
!
!--for each particle in the current cell, loop over its neighbours
!
       loop_over_neighbours: DO n = idone+1,nneigh
       
          j = listneigh(n)
	  IF ((j.NE.i).AND..NOT.(j.GT.npart .AND. i.GT.npart)) THEN		! don't count particle with itself
!	     print*,' ... neighbour, h=',j,hh(j),rho(j),x(:,j)
	     dx(:) = x(:,i) - x(:,j)
	     hj = hh(j)
	     hj2 = hj*hj
!
!--calculate averages of smoothing length if using this averaging
!			 
	     hav = 0.5*(hi + hj)
	     h2 = hav*hav
	     hfacwab = 1./hav**ndim
	     hfacwabj = 1./hj**ndim
	     hfacgrkern = 1./hav**(ndim+1)
	     hfacgrkernj = 1./hj**(ndim+1)
	     
	     rij2 = DOT_PRODUCT(dx,dx)
	     rij = SQRT(rij2)
	     IF (rij.EQ.0.) THEN
	        WRITE(iprint,*) 'rates: dx = 0 i,j,dx,hi,hj=',i,j,dx,hi,hj
                CALL quit
	     ENDIF	
	     q2 = rij2/h2
	     q2i = rij2/hi2
	     q2j = rij2/hj2
	     dr(1:ndim) = dx(1:ndim)/rij
	     IF (ndimV.GT.ndim) dr(ndim+1:ndimV) = 0.
!
!--linear interpolation from kernel table to get kernel and its derivative
!
	     IF ((q2.LT.radkern2).OR.(q2i.LT.radkern2)	& ! if within 2h
                                 .OR.(q2j.LT.radkern2)) THEN
!
!--use either average h, average kernel gradient or Springel/Hernquist type
!
		IF (ikernav.EQ.1) THEN		
		   index = INT(q2/dq2table)	! nearest index in table
		   index1 = index + 1
		   IF (index.GT.ikern) index = ikern
		   IF (index1.GT.ikern) index1 = ikern
		   dxx = q2 - index*dq2table 	! increment along from index pt
		   dgrwdx =  (grwij(index1)-grwij(index))/dq2table ! slope
		   grkern = (grwij(index)+ dgrwdx*dxx)*hfacgrkern	! divide by h**3. in 1D		
		   grkerni = grkern	! ie no variable h terms
		   grkernj = grkern	!  "  "    "      "   "
		   dwdx =  (wij(index1)-wij(index))/dq2table ! slope
		   wab = (wij(index)+ dwdx*dxx)*hfacwab	! divide by h in 1D	
                ELSE
!  (using hi)
		   index = INT(q2i/dq2table)	! nearest index in table
		   index1 = index + 1
		   IF (index.GT.ikern) index = ikern
		   IF (index1.GT.ikern) index1 = ikern
		   dxx = q2i - index*dq2table 	! increment along from index pt
		   dgrwdx =  (grwij(index1)-grwij(index))/dq2table ! slope
		   grkerni = (grwij(index)+ dgrwdx*dxx)*hfacgrkerni	! divide by h**3. in 1D		
		   dwdx =  (wij(index1)-wij(index))/dq2table ! slope
		   wabi = (wij(index)+ dwdx*dxx)*hfacwabi	! divide by h in 1D		
!  (using hj)
		   index = INT(q2j/dq2table)	! nearest index in table
		   index1 = index + 1
		   IF (index.GT.ikern) index = ikern
		   IF (index1.GT.ikern) index1 = ikern
		   dxx = q2j - index*dq2table 	! increment along from index pt
		   dgrwdx =  (grwij(index1)-grwij(index))/dq2table ! slope
		   grkernj = (grwij(index)+ dgrwdx*dxx)*hfacgrkernj	! divide by h**3. in 1D		
		   dwdx =  (wij(index1)-wij(index))/dq2table ! slope
		   wabj = (wij(index)+ dwdx*dxx)*hfacwabj	! divide by h in 1D		
!  (calculate average)
		   grkern = 0.5*(grkerni + grkernj)
		   wab = 0.5*(wabi + wabj)
		   
		   IF (ikernav.NE.3) THEN  ! if not using grad h correction
		      grkerni = grkern
		      grkernj = grkern
		   ENDIF
		ENDIF
!
!--define local copies of quantities (these are shared in module rates_local)
!
	        velj(:) = vel(:,j)		
	        dvel(:) = veli(:) - velj(:)
		dvdotr = DOT_PRODUCT(dvel,dr)
	        rhoj = rho(j)
		rho1j = 1./rhoj
		rho2j = rhoj*rhoj
	        rhoj5 = SQRT(rhoj)
		rhoij = rhoi*rhoj
		rhoav = 0.5*(rhoi+rhoj)
		prj = pr(j)
	        Prho2j = pr(j)/rho2j
		spsoundj = spsound(j)
	        pmassj = pmass(j)
		gradhj = 1. - gradh(j)
!
!--mhd definitions (these are shared in rates_mhdlocal)
!
  	        IF (imhd.GE.11) THEN	! if B is mag field variable
		 Bj(:) = Bfield(:,j)
		 Brhoj(:) = Bj(:)/rhoj
		ELSEIF (imhd.NE.0) THEN	! if B/rho is mag field variable
		 Brhoj(:) = Bfield(:,j)
		 Bj(:) = Brhoj(:)*rhoj
		ELSE			! no MHD
		 Brhoj(:) = 0.
		 Bj(:) = 0.
		ENDIF
		IF (imhd.NE.0) THEN
		 dB(:) = Bi(:) - Bj(:)
                 projBi = DOT_PRODUCT(Bi,dr)
		 projBj = DOT_PRODUCT(Bj,dr)
		 projdB = DOT_PRODUCT(dB,dr)
		 projBrhoi = DOT_PRODUCT(Brhoi,dr)
		 projBrhoj = DOT_PRODUCT(Brhoj,dr)
		 Brho2j = DOT_PRODUCT(Brhoj,Brhoj)
		 valfven2j = Brho2j*rhoj
		 valfvenj = SQRT(valfven2j)
                ENDIF
!
!--calculate the maximum signal velocity (used in AV and timestep control)
!
!--maximum velocity for timestep control
		vmag = SQRT(DOT_PRODUCT(dvel,dvel))
!		vsigdtc = vmag + spsoundi + spsoundj  	&
!			       + valfveni + valfvenj
		vsig = 0.
		viss = 0.
		IF (vmag.EQ.0.) vmag = 1.e-6
		vunit(:) = abs(dvel(:))/vmag
	        IF ((iavlim.EQ.1).AND.(ndimV.GT.ndim)) THEN
		   avsource(i) = avsource(i) + pmassj*(DOT_PRODUCT(dvel,vunit))*grkern
		   avsource(j) = avsource(j) + pmassi*(DOT_PRODUCT(dvel,vunit))*grkern
		ENDIF
!
!--calculate maximum signal velocity (this is used for the timestep control
!  and also in the artificial viscosity)
!
		IF (dvdotr.LT.0.) THEN
		   vunit(:) = -dvel(:)/vmag	! unit vector in direction of velocity
		   IF (ndimV.GT.ndim) THEN	! if using 1.5D or 2.5D
		      viss = abs(DOT_PRODUCT(dvel,vunit))
		   ELSE				! normal
		      viss = abs(dvdotr)		    
		   ENDIF
		ENDIF   
!
!--max signal velocity (Joe's)
!
!		vsigi = 0.5*(SQRT(spsoundi**2 + valfven2i	 	&
!      		    - 2.*spsoundi*projBi/rhoi5)	&
!                             +SQRT(spsoundi**2 + valfven2i 		&
!      		    + 2.*spsoundi*projBi/rhoi5))
!		vsigj = 0.5*(SQRT(spsoundj**2 + valfven2j		&
!             	    - 2.*spsoundj*projBj/rhoj5)	&
!                            +SQRT(spsoundj**2 + valfven2j		&
!     		    + 2.*spsoundj*projBj/rhoj5))
!
!--max signal velocity
!
                vsig2i = spsoundi**2 + valfven2i
		vsig2j = spsoundj**2 + valfven2j				
		
		vsigproji = vsig2i**2 - 4*(spsoundi*projBi)**2*rho1i
		vsigprojj = vsig2j**2 - 4*(spsoundj*projBj)**2*rho1j
		
		vsigi = SQRT(0.5*(vsig2i + SQRT(vsigproji)))
		vsigj = SQRT(0.5*(vsig2j + SQRT(vsigprojj)))
!
!--magnetosonic speed
!		vsigi = SQRT(spsoundi**2+valfven2i+beta*viss*viss)
!		vsigj = SQRT(spsoundj**2+valfven2j+beta*viss*viss)

		vsig = vsigi + vsigj + beta*viss
	
		vsigdtc = vsig		! this is Joe's sigma
!
!--compute artificial dissipation terms if particles are approaching
!		
                IF ((iav.NE.0).AND.(dvdotr.LT.0)) CALL artvis_terms(i,j)
!
!--time step control (courant and viscous)
!
		vsigdtc = vsigi + vsigj + beta*abs(dvdotr)
	        dtcourant = min(dtcourant,0.8*hav/vsigdtc)
!
!--time derivative of density (continuity equation)
!	        
		IF (icty.EQ.2) THEN	! good for rho discontinuous
		   drhodti = rhoi*pmassj*dvdotr*grkern/rhoj
		   drhodtj = rhoj*pmassi*dvdotr*grkern/rhoi				   
		ELSE ! compute this even if direct sum - gives divv for art vis.
		   drhodti = pmassj*dvdotr*grkerni/gradhi	! gradhi = 1 if not set
		   drhodtj = pmassi*dvdotr*grkernj/gradhj
		ENDIF
				
		drhodt(i) = drhodt(i) + drhodti
		drhodt(j) = drhodt(j) + drhodtj
!
!--pressure term (choice of several)
!
	        IF (iprterm.EQ.2) THEN	! standard alternative form
		   prterm = (pri*grkerni + prj*grkernj)/rhoij
		ELSEIF (iprterm.EQ.3) THEN	! Hernquist/Katz
		   prterm = (2.*SQRT(pri*prj)/rhoij)*grkern		   
		ELSEIF (iprterm.NE.0) THEN	! default
		   prterm = Prho2i*grkerni/gradhi 		&
     		          + Prho2j*grkernj/gradhj		   
		ENDIF

		force(:,i) = force(:,i) - pmassj*prterm*dr(:)
		force(:,j) = force(:,j) + pmassi*prterm*dr(:)		
!
!--subtract external forces
!	
		IF (iexternal_force.EQ.1) THEN	! linear force
		   force(1:ndim,i) = force(1:ndim,i) - x(1:ndim,i)
		   force(1:ndim,j) = force(1:ndim,j) - x(1:ndim,j)
		ENDIF
!		
!--Lorentz force and time derivative of B terms
!
                IF (imhd.NE.0) CALL mhd_terms(i,j)
!
!--energy equation
!
                IF (iener.NE.0) CALL energy_terms(i,j)
!
!--calculate XSPH term for moving the particles
!
	        IF (ixsph.EQ.1) THEN
		   xsphterm(:,i) = xsphterm(:,i) - pmassj*dvel(:)/rhoav*wab
		   xsphterm(:,j) = xsphterm(:,j) + pmassi*dvel(:)/rhoav*wab
		ENDIF

	     ELSE	! if outside 2h
!		  PRINT*,'outside 2h, not calculated, r/h=',sqrt(q2)
	     ENDIF	! r < 2h or not

	   ENDIF		! j.ne.i
        
	  ENDDO loop_over_neighbours

       iprev = i
       IF (iprev.NE.-1) i = ll(i)		! possibly should be only IF (iprev.NE.-1)
    
    ENDDO loop_over_cell_particles
            
 ENDDO loop_over_cells
!
!--do the divisions by rho etc (this is for speed - so calculations are not
!  done multiple times within the loop)
!
 DO i=1,npart
    rhoi = rho(i)
    IF (imhd.NE.0) THEN
       divB(i) = divB(i)/rhoi		!*rhoi
       IF ((imhd.EQ.1).OR.(imhd.EQ.3).OR.(imhd.EQ.11)) THEN
          dBfielddt(:,i) = dBfielddt(:,i)/rhoi
       ENDIF
    ENDIF
!
!--calculate time derivative of alpha (artificial dissipation coefficient)
!  see Morris and Monaghan (1997)
!     
    IF (iavlim.EQ.1) THEN
!       IF (ndimV.GT.ndim) THEN
!          print*,' source ',i,' = ',avsource(i)/rhoi,drhodt(i)/rhoi
!          source = max(avsource(i)/rhoi,0.0)
!       ELSE    
          source = MAX(drhodt(i)/rhoi,0.0)    ! source term is div v
!       ENDIF
       IF (imhd.GE.11) THEN
          valfven2i = DOT_PRODUCT(Bfield(:,i),Bfield(:,i))/rhoi
       ELSE
          valfven2i = DOT_PRODUCT(Bfield(:,i),Bfield(:,i))*rhoi       
       ENDIF
       vsig = SQRT(spsound(i)**2. + valfven2i) 	! approximate vsig only
       tdecay = hh(i)/(avdecayconst*vsig)	! decay time (use vsig)
       daldt(i) = (alphamin - alpha(i))/tdecay + avfact*source
    ENDIF
      
 ENDDO
     
 IF (trace) WRITE(iprint,*) ' Exiting subroutine get_rates'
      
 RETURN
 END SUBROUTINE get_rates
