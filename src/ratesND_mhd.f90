!!--------------------------------------------------------------------
!! Computes the rates of change of everything (forces, energy eqn etc)
!! This is the core of the SPH algorithm
!!
!! Changes log:
!! v3.6:
!! 15/10/03 bug fix in artificial stress
!!--------------------------------------------------------------------

SUBROUTINE get_rates
! USE dimen_mhd
 USE debug
 USE loguns
 
 USE artvi
 USE eos
 USE gravity
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
  
 USE fmagarray
 USE derivB
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i,j,n
 INTEGER :: icell,ipart,iprev,ncell,nneigh,index,index1
 INTEGER, ALLOCATABLE, DIMENSION(:) :: listneigh
 INTEGER :: idone
 INTEGER, DIMENSION(3*2**(ndim-1) -1) :: neighcell
!
!  (particle properties - local copies and composites)
!
 REAL :: rij,rij2
 REAL :: rhoi,rho1i,rho2i,rho21i,rhoj,rho1j,rho2j,rho21j,rhoav,rhoav1,rhoij
 REAL :: pmassi,pmassj
 REAL :: pri,prj,pr2i,pr2j,Prho2i,Prho2j,prterm
 REAL :: hi,hi1,hj,hj1,hi2,hi21,hj2,hj21,hi3,hj3,hav,h2,h3,hav1,h21
 REAL :: drhodti,drhodtj
 REAL :: hfacwab,hfacwabi,hfacwabj,hfacgrkern,hfacgrkerni,hfacgrkernj
 REAL, DIMENSION(ndim) :: dx
!
!  (velocity)
!      
 REAL, DIMENSION(ndimV) :: veli,velj,dvel,vrho2i,vrho2j
 REAL, DIMENSION(ndimV) :: vterm,prvterm
 REAL, DIMENSION(ndimV) :: dr
 REAL :: dvdotr,projprv,projvterm,v2i,v2j
!
!  (mhd)
!     
 REAL, DIMENSION(ndimB) :: Brhoi,Brhoj,Bi,Bj,dB
 REAL, DIMENSION(ndimB) :: faniso,fmagi,fmagj,fcorr
 REAL, DIMENSION(ndimB) :: curlBi,divBvec
 REAL :: fiso,divBonrho,divBvec
 REAL :: valfveni,valfvenj,valfven2i,valfven2j
 REAL :: BidotdB,BjdotdB,Brho2i,Brho2j
 REAL :: projBrhoi,projBrhoj,projBi,projBj,projdB
 REAL :: prvaniso,prvaniso2  
 REAL :: Bidotvj,Bjdotvi
!
!  (mhd art. vis)
!
 REAL, DIMENSION(ndimB) :: Bav,vunit,Bvisc,dBdtvisc
 REAL :: Bab,vsigmag,rhoi5,rhoj5,B2i,B2j,ediffB
 REAL :: vsig2i,vsig2j,vsigproji,vsigprojj
 REAL :: vissv,vissB,vissu
!
! (my av switch)
!
 REAL :: qab,vdotri,vdotrj
 REAL, ALLOCATABLE, DIMENSION(:) :: avsource
!
!  (kernel related quantities)
!
 REAL :: q2,q2i,q2j
 REAL :: grkern,grkerni,grkernj,dgrwdx,dxx
 REAL :: wab,wabi,wabj,dwdx
!
!  (joe's mhd fix)
!
 INTEGER :: indexjoe,indexjoe1,ndim1
 REAL :: wabjoe,dwdxjoe,dxxjoe,wabjoei,wabjoej,deltap2
 REAL :: Rjoe,q2joe,prvanisoi,prvanisoj,Bidotvi,Bjdotvj
!
!  (artificial viscosity quantities)
!      
 REAL :: vsig,vsigi,vsigj,vsigav,viss,eni,enj,ediff,qdiff
 REAL :: spsoundi,spsoundj,visc,envisc,uvisc
 REAL :: alphai,alphaav,source,tdecay    
!
!  (time step criteria)
!      
 REAL :: vsigdtc,vmag
!
!  (variable smoothing length terms)
!
 REAL :: gradhi,gradhj
 INTEGER :: ierr
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine get_rates'
!
!--initialise quantities
!      
 dtcourant = 1.e6  

 nlistdim = ntotal
 ALLOCATE ( listneigh(nlistdim),STAT=ierr )
 ALLOCATE (avsource(ntotal), STAT=ierr )	! should deallocate at end of subroutine
 IF (ierr.NE.0) WRITE(iprint,*) ' Error allocating avsource, ierr = ',ierr

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
 
!  CALL interpolate_kernel(q2joe,wabjoei,grkerni)
 
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
       rho2i = rhoi*rhoi
       rhoi5 = SQRT(rhoi)
       rho1i = 1./rhoi
       rho21i = rho1i*rho1i       
       pri = pr(i)
       pr2i = pr(i)*pr(i)
       Prho2i = pr(i)*rho21i
       spsoundi = spsound(i)
       veli(:) = vel(:,i)
       vrho2i(:) = veli(:)*rho21i
       pmassi = pmass(i)
       alphai = alpha(i)
       IF (imhd.GE.11) THEN	! if mag field variable is B
          Bi(:) = Bfield(:,i)
          Brhoi(:) = Bi(:)*rho1i
       ELSEIF (imhd.GT.0) THEN	! if mag field variable is B/rho
          Brhoi(:) = Bfield(:,i)
          Bi(:) = Brhoi(:)*rhoi
       ENDIF
! mhd definitions
       Brho2i = DOT_PRODUCT(Brhoi,Brhoi)
       valfven2i = Brho2i*rhoi
       valfveni = SQRT(valfven2i)
              
       gradhi = 1./(1. - gradh(i))
       IF (gradhi.LE.0.5) THEN
          WRITE(iprint,*) 'Error in grad h terms, part ',i,gradhi
       ENDIF   
       hi = hh(i)
       hi1 = 1./hi
       hi21 = hi1*hi1
       hi2 = hi*hi		    
       hi3 = hi*hi2
       IF (hi.LE.0.) THEN
          WRITE(iprint,*) ' rates: h <= 0 particle',i,hi
	  CALL quit
       ENDIF
       hfacwabi = hi1**ndim
       hfacgrkerni = hfacwabi*hi1
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
	     hj3 = hj*hj2
	     hj1 = 1./hj
	     hj21 = hj1*hj1
!
!--calculate averages of smoothing length if using this averaging
!			 
	     hav = 0.5*(hi + hj)
	     h2 = hav*hav
	     h3 = hav*h2
	     hav1 = 1./hav
	     h21 = hav1*hav1
	     hfacwab = hav1**ndim
	     hfacwabj = hj1**ndim
	     hfacgrkern = hfacwab*hav1
	     hfacgrkernj = hfacwabj*hj1
	     
	     rij2 = DOT_PRODUCT(dx,dx)
	     rij = SQRT(rij2)
	     IF (rij.EQ.0.) THEN
	        WRITE(iprint,*) 'rates: dx = 0 i,j,dx,hi,hj=',i,j,dx,hi,hj
                CALL quit
	     ENDIF	
	     q2 = rij2*h21
	     q2i = rij2*hi21
	     q2j = rij2*hj21
	     dr(1:ndim) = dx(1:ndim)/rij	! unit vector
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
	           CALL interpolate_kernel(q2,wab,grkern)
	           wab = wab*hfacwab
	           grkern = grkern*hfacgrkern
		   grkerni = grkern
		   grkernj = grkern
                ELSE
!  (using hi)
                   CALL interpolate_kernel(q2i,wabi,grkerni)
	           wabi = wabi*hfacwabi
	           grkerni = grkerni*hfacgrkerni
!  (using hj)
	           CALL interpolate_kernel(q2j,wabj,grkernj)
                   wabj = wabj*hfacwabj
	           grkernj = grkernj*hfacgrkernj
!  (calculate average)  		
                   grkern = 0.5*(grkerni + grkernj)
	           wab = 0.5*(wabi + wabj)
!  (grad h terms)  
		   IF (ikernav.EQ.3) THEN  ! if using grad h correction
 		      gradhj = 1./(1. - gradh(j))
		      grkerni = grkerni*gradhi
		      grkernj = grkernj*gradhj
		   ELSE  ! if not using grad h correction		   
		      grkerni = grkern
		      grkernj = grkern
		   ENDIF
		ENDIF
!
!--define local copies of quantities
!
	        velj(:) = vel(:,j)		
	        dvel(:) = veli(:) - velj(:)
		dvdotr = DOT_PRODUCT(dvel,dr)
	        rhoj = rho(j)
		rho1j = 1./rhoj
		rho2j = rhoj*rhoj
		rho21j = rho1j*rho1j
	        rhoj5 = SQRT(rhoj)
		rhoij = rhoi*rhoj
		prj = pr(j)
		pr2j = pr(j)*pr(j)		
	        Prho2j = pr(j)*rho21j
		spsoundj = spsound(j)
		vrho2j(:) = velj(:)*rho21j
	        pmassj = pmass(j)
		alphaav = 0.5*(alphai + alpha(j))		
!  (mhd definitions)
  	        IF (imhd.GE.11) THEN	! if B is mag field variable
		 Bj(:) = Bfield(:,j)
		 Brhoj(:) = Bj(:)*rho1j
		ELSEIF (imhd.NE.0) THEN	! if B/rho is mag field variable
		 Brhoj(:) = Bfield(:,j)
		 Bj(:) = Brhoj(:)*rhoj		    	          
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
!--artificial viscosity terms
!
	        rhoav = 0.5*(rhoi + rhoj)
		rhoav1 = 1./rhoav
!--maximum velocity for timestep control
		vmag = SQRT(DOT_PRODUCT(dvel,dvel))
!		vsigdtc = vmag + spsoundi + spsoundj  	&
!			       + valfveni + valfvenj
		vsig = 0.
		viss = 0.
		visc = 0.
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
		vsigi = 0.5*(SQRT(spsoundi**2 + valfven2i	 	&
      		    - 2.*spsoundi*projBi/rhoi5)	&
                             +SQRT(spsoundi**2 + valfven2i 		&
      		    + 2.*spsoundi*projBi/rhoi5))
		vsigj = 0.5*(SQRT(spsoundj**2 + valfven2j		&
             	    - 2.*spsoundj*projBj/rhoj5)	&
                            +SQRT(spsoundj**2 + valfven2j		&
     		    + 2.*spsoundj*projBj/rhoj5))
!
!--max signal velocity (my version)
!
!                vsig2i = spsoundi**2 + valfven2i
!		vsig2j = spsoundj**2 + valfven2j				
!		
!		vsigproji = vsig2i**2 - 4.*(spsoundi*projBi)**2*rho1i
!		vsigprojj = vsig2j**2 - 4.*(spsoundj*projBj)**2*rho1j
		
!		IF (vsigproji.LT.0. .OR. vsigprojj.LT.0.) THEN
!		   WRITE(iprint,*) ' rates: vsig det < 0 ',   &
!		     'i: ',vsigproji,vsig2i**2,  	&
!		     4*(spsoundi*projBi)**2*rho1i,	&
!		     spsoundi,projBi,rho1i,		&
!		     'j: ',vsigprojj,vsig2j**2,   	&
!		     4*(spsoundj*projBj)**2*rho1j,	&
!		     spsoundj,projBj,rho1j
!		   CALL quit  
!		ENDIF
		
!		vsigi = SQRT(0.5*(vsig2i + SQRT(vsigproji)))
!		vsigj = SQRT(0.5*(vsig2j + SQRT(vsigprojj)))
!
!--magnetosonic speed
!		vsigi = SQRT(spsoundi**2+valfven2i+beta*viss*viss)
!		vsigj = SQRT(spsoundj**2+valfven2j+beta*viss*viss)
!--crap version		    
!		vsigi = MAX(spsoundi,SQRT(valfven2i))
!		vsigj = MAX(spsoundj,SQRT(valfven2j))

		vsig = vsigi + vsigj + beta*viss
	
!		vsigdtc = vsig		! this is Joe's sigma
                vsigdtc = vsigi + vsigj + beta*abs(dvdotr)
		
		IF (dvdotr.LT.0. .AND. iav.NE.0) THEN	! only for approaching particles
!--viscosity (kinetic energy term)
                   visc = 0.5*alphaav*vsig*viss*rhoav1
		   IF (ndimV.GT.ndim) THEN
		      vissv = -0.5*(DOT_PRODUCT(veli,vunit)	&
		                  - DOT_PRODUCT(velj,vunit))**2
!		      vissv = -0.5*DOT_PRODUCT(dvel,dvel)
     		   ELSE
		      vissv = -0.5*(DOT_PRODUCT(veli,dr) 	&
		                  - DOT_PRODUCT(velj,dr))**2		    
		   ENDIF
!--ohmic dissipation (magnetic energy term)
		   Bvisc(:) = dB(:) - dr(:)*projdB
		   dBdtvisc(:) = 0.5*alphaav*vsig*Bvisc(:)*rhoav1**2
		   vissB = -0.5*(DOT_PRODUCT(dB,dB)-projdB**2)*rhoav1
!--vissu is the dissipation energy from thermal conductivity
		   vissu = gconst*(uu(i) - uu(j))
!--envisc is the total contribution to the thermal energy equation
		   envisc = 0.5*alphaav*vsig*(vissv+vissB)*rhoav1*grkern!*abs(rx)
	           uvisc = 0.5*alphaav*vsig*vissu*rhoav1*grkern		   
		ELSE
		   visc = 0.0
		   dBdtvisc(:) = 0.
		   envisc = 0.
		   uvisc = 0.
!     		   vissB = -0.5*(DOT_PRODUCT(dB,dB)-projdB**2)
!		   dBdtvisc(:) = 0.5*alphaav*vsig*Bvisc(:)*rhoav1**2
!		   envisc = 0.5*alphaav*vsig*(vissB)*rhoav1*grkern!*abs(rx)
!		   print*,'Bvisc, envisc = ',dBdtvisc,envisc
	        ENDIF		
!
!--time step control (courant and viscous)
!
!		IF (0.8*hav/vsigdtc .LT. dtcourant) THEN
!		   PRINT*,' new dtcourant = ',0.8*hav/vsigdtc,hav,vsigdtc,i,j
!		ENDIF
	        dtcourant = min(dtcourant,0.8*hav/vsigdtc)
!
!--time derivative of density (continuity equation)
!	        
		IF (icty.EQ.2) THEN	! good for rho discontinuous
		   drhodti = rhoi*pmassj*dvdotr*grkern*rho1j
		   drhodtj = rhoj*pmassi*dvdotr*grkern*rho1i				   
		ELSE ! compute this even if direct sum - gives divv for art vis.
		   drhodti = pmassj*dvdotr*grkerni	! gradhi = 1 if not set
		   drhodtj = pmassi*dvdotr*grkernj
		ENDIF
				
		drhodt(i) = drhodt(i) + drhodti
		drhodt(j) = drhodt(j) + drhodtj
!
!--pressure term (choice of several)
!
	        IF (iprterm.EQ.1) THEN	! standard alternative form
		   prterm = (pri*grkerni + prj*grkernj)/rhoij + visc*grkern       		   	      
		ELSEIF (iprterm.EQ.2) THEN	! Hernquist/Katz
		   prterm = (2.*SQRT(pri*prj)/rhoij + visc)*grkern		   
		ELSEIF (iprterm.EQ.4) THEN	! use in 1.5D
		   prterm = Prho2i*grkerni 		&
     		          + Prho2j*grkernj		   
		ELSE		! default
		   prterm = Prho2i*grkerni 		&
     		          + Prho2j*grkernj + visc*grkern       		   
		ENDIF
!
!--add pressure and viscosity terms to force (equation of motion)
!
		IF (iprterm.EQ.4) THEN
		  force(:,i) = force(:,i)	&
		             - pmassj*(prterm*dr(:)+visc*vunit(:)*grkern)
		  force(:,j) = force(:,j)	&
		             + pmassi*(prterm*dr(:)+visc*vunit(:)*grkern)  
!		  force(1,i) = force(1,i)	&
!		             - pmassj*(prterm*dr(1)+visc*dr(1)*grkern)
!		  force(1,j) = force(1,j)	&
!		             + pmassi*(prterm*dr(1)+visc*dr(1)*grkern)  
!		  force(2,i) = force(2,i)	&
!		             - pmassj*(viscy*dr(1)*grkern)
!		  force(2,j) = force(2,j)	&
!		             + pmassi*(viscy*dr(1)*grkern)  

		ELSE
		  force(:,i) = force(:,i) - pmassj*prterm*dr(:)
		  force(:,j) = force(:,j) + pmassi*prterm*dr(:)		
                ENDIF	
!		
!--Lorentz force and time derivative of B terms
!
                IF (imhd.NE.0) THEN
!
!--calculate curl B for current (only if field is 3D)
!
		  IF (ndimB.EQ.3) THEN
		     curlBi(1) = dB(2)*dr(3) - dB(3)*dr(2)
		     curlBi(2) = dB(3)*dr(1) - dB(1)*dr(3)
		     curlBi(3) = dB(1)*dr(2) - dB(2)*dr(1)
		  ENDIF

		  IF (imagforce.EQ.1) THEN	! vector form (dot products)
		  		  
		  BidotdB = DOT_PRODUCT(Brhoi,dB)
		  BjdotdB = DOT_PRODUCT(Brhoj,dB)
		  fmagi(:) = grkern*(dr(:)*BidotdB - projBrhoi*dB(:))		  
		  fmagj(:) = grkern*(dr(:)*BjdotdB - projBrhoj*dB(:))
		  		  			
		  ELSEIF (imagforce.EQ.2) THEN	! tensor formalism
!
!--isotropic mag force (pressure gradient)
!		  
	          fiso = 0.5*(Brho2i*grkerni + Brho2j*grkernj)
!
!--anisotropic force (including variable smoothing length terms)
!
		  faniso(:) = Brhoi(:)*projBrhoi*grkerni  &
		            + Brhoj(:)*projBrhoj*grkernj		   
!
!--Joe's correction term (preserves momentum conservation)
!  *** in 1D only do this in the x-direction?
		  IF (ianticlump.EQ.1) THEN
		     wabjoe = wabjoei*hfacwab
		     Rjoe = 0.5*eps*(wab/wabjoe)**neps
!		     IF (Rjoe.GT.0.1) PRINT*,'Rjoe = ',Rjoe,i,j,wab,wabjoe
!		     faniso(:) = faniso(:) - Rjoe*faniso(:)  
		     faniso(1:ndim) = faniso(1:ndim) - Rjoe*faniso(1:ndim)  
		  ENDIF
!
!--add contributions to magnetic force
!	  
	          fmagi(:) = faniso(:) - fiso*dr(:)
		  
		  ELSEIF (imagforce.EQ.3) THEN	! D.Price version

		  fiso = grkern*0.5*(Brho2i + Brho2j)
		  faniso(:) = grkern*(Brhoi(:)*projBrhoj + Brhoj(:)*projBrhoi)     
	          fmagi(:) = faniso(:) - fiso*dr(:)
		  		  
		  ELSEIF (imagforce.EQ.4) THEN ! form that goes with alt. time deriv

		  rhoij = rhoi*rhoj
		  fiso = grkern*0.5*(Brho2i + Brho2j)
		  faniso(:) = grkern*(Bi(:)*projBi + Bj(:)*projBj)/rhoij
	          fmagi(:) = faniso(:) - fiso*dr(:)
		  	  
		  ELSEIF (imagforce.EQ.5) THEN	! Morris' Hybrid form

		  rhoij = rhoi*rhoj
		  fiso = grkern*0.5*(Brho2i + Brho2j)
		  faniso(:) = grkern*(Bj(:)*projBj - Bi(:)*projBi)/rhoij
	          fmagi(:) = faniso(:) - fiso*dr(:)		  
		  
		  ENDIF
!
!--compute divergence of B
!
		  divB(i) = divB(i) - pmassj*DOT_PRODUCT(dB,dr)*grkern
		  divB(j) = divB(j) - pmassi*DOT_PRODUCT(dB,dr)*grkern
!
!--compute current density J
!
		  curlB(:,i) = curlB(:,i) - pmassj*curlBi(:)*grkern
		  curlB(:,j) = curlB(:,j) - pmassi*curlBi(:)*grkern
!
!--stabilising correction term from Borve et al. (2001)
!
                  IF (idivBzero.EQ.3) THEN	          
		     divBvec(:) = Brhoi(:)*rho1i + Brhoj(:)*rho1j
		     divBonrho = grkern*(DOT_PRODUCT(divBvec,dr))
		     fcorr(:) = Bi(:)*divBonrho
		     fmagi(:) = fmagi(:) - fcorr(:)
		  ENDIF		   
!
!--add Lorentz force to total force
!        	
		  IF (imagforce.EQ.1) THEN
		     fmag(:,i) = fmag(:,i) + pmassj*fmagi(:)/rho2i
		     fmag(:,j) = fmag(:,j) - pmassi*fmagj(:)/rho2j
		     force(:,i) = force(:,i) + pmassj*fmagi(:)/rho2i
		     force(:,j) = force(:,j) - pmassi*fmagj(:)/rho2j
		  ELSEIF (imagforce.EQ.5) THEN	! Morris' Hybrid force
		     fmag(:,i) = fmag(:,i) + pmassj*(faniso(:)-fiso*dr(:))
		     fmag(:,j) = fmag(:,j) + pmassi*(faniso(:)+fiso*dr(:))
		     force(:,i) = force(:,i) + pmassj*(faniso(:)-fiso*dr(:))
		     force(:,j) = force(:,j) + pmassi*(faniso(:)+fiso*dr(:)) 			        
		  ELSE 	! symmetric forces fmagxi = -fmagxj
		     fmag(:,i) = fmag(:,i) + pmassj*fmagi(:)
		     fmag(:,j) = fmag(:,j) - pmassi*fmagi(:)
		     force(:,i) = force(:,i) + pmassj*fmagi(:)
		     force(:,j) = force(:,j) - pmassi*fmagi(:) 
		  ENDIF
!
!--time derivative of magnetic field (divide by rho later)
!
!  (evolving B/rho)
	          IF (imhd.EQ.1) THEN
!		     IF (ianticlump.EQ.1 .AND. imagforce.EQ.2) THEN
!		     dBfielddt(:,i) = dBfielddt(:,i)		&
!                   - pmassj*(dvel(:)*projBrhoi)*grkerni*(1.-Rjoe) 
!		     dBfielddt(:,j) = dBfielddt(:,j) 		&
!     		   - pmassi*(dvel(:)*projBrhoj)*grkernj*(1.-Rjoe)	     
!		     ELSE
		     dBfielddt(:,i) = dBfielddt(:,i)		&
                   - pmassj*(dvel(:)*projBrhoi)*grkerni 
		     dBfielddt(:,j) = dBfielddt(:,j) 		&
     		   - pmassi*(dvel(:)*projBrhoj)*grkernj
!		     ENDIF

		     IF (iav.EQ.2) THEN		! add dissipative term ***check
	                dBfielddt(:,i) = dBfielddt(:,i)		&
                         + rhoi*pmassj*dBdtvisc(:)*grkern		   
		        dBfielddt(:,j) = dBfielddt(:,j)		&
                         - rhoj*pmassi*dBdtvisc(:)*grkern
	             ENDIF

                  ELSEIF (imhd.EQ.2) THEN	! good for rho discts
		     dBfielddt(:,i) = dBfielddt(:,i)		&
                      - pmassj*(dvel(:)*projBrhoi*rho1j)*grkern		   
		     dBfielddt(:,j) = dBfielddt(:,j) 		&
     		      - pmassi*(dvel(:)*projBrhoj*rho1i)*grkern
		  ELSEIF (imhd.EQ.3) THEN	! add v*div B term
		     dBfielddt(:,i) = dBfielddt(:,i)		&
                     - pmassj*(dvel(:)*projBrhoi		&
                              +veli(:)*DOT_PRODUCT(dB,dr))*grkern		   
		     dBfielddt(:,j) = dBfielddt(:,j) 		&
      		      - pmassi*(dvel(:)*projBrhoj		&
                                +velj(:)*DOT_PRODUCT(dB,dr))*grkern
!   (evolving B)
		  ELSEIF (imhd.EQ.11) THEN	! note divided by rho later		  
		     dBfielddt(:,i) = dBfielddt(:,i)		&
            + pmassj*(Bi(:)*dvdotr - dvel(:)*projBi)*grkerni   
		     dBfielddt(:,j) = dBfielddt(:,j) 		&
            + pmassi*(Bj(:)*dvdotr - dvel(:)*projBj)*grkernj
		  
		     IF (iav.EQ.2) THEN		! add dissipative term
	                dBfielddt(:,i) = dBfielddt(:,i)		&
                         + rho2i*pmassj*dBdtvisc(:)*grkern		   
		        dBfielddt(:,j) = dBfielddt(:,j)		&
                         - rho2j*pmassi*dBdtvisc(:)*grkern
	             ENDIF
		  
		  ENDIF
		
		ELSE	! if no mhd
		   Brho2j = 0.
		   Brho2i = 0.  
	           Brhoi(:) = 0.
		   Brhoj(:) = 0.
		   projBrhoi = 0.
		   projBrhoj = 0.
		ENDIF

!
!--energy equation (choice of several)
!               
	        IF (iener.EQ.2) THEN	! 1st law of thermodynamics
                    dudt(i) = dudt(i) + Prho2i*drhodti + pmassj*(envisc + uvisc)
     		    dudt(j) = dudt(j) + Prho2j*drhodtj + pmassi*(envisc - uvisc)
		ELSEIF (iener.GE.3) THEN	! total energy/particle
		   Bidotvj = DOT_PRODUCT(Brhoi,velj)
		   Bjdotvi = DOT_PRODUCT(Brhoj,veli)
! (isotropic stress)
		   prvterm(:) = (Prho2i+0.5*Brho2i)*velj(:)*grkerni &
                              + (Prho2j+0.5*Brho2j)*veli(:)*grkernj
! (anisotropic stress)
		   prvaniso =  - Bidotvj*projBrhoi*grkerni	& 
		               - Bjdotvi*projBrhoj*grkernj
		   projprv = DOT_PRODUCT(prvterm,dr)		   
! (Joe's anticlumping term)		   
		   IF (ianticlump.EQ.1 .AND. imagforce.EQ.2) THEN
!		      prvaniso = prvaniso - Rjoe*prvaniso
! (if applied in x-direction only)
                      Bidotvi = DOT_PRODUCT(Brhoi(1:ndim),veli(1:ndim))
		      Bjdotvi = DOT_PRODUCT(Brhoj(1:ndim),veli(1:ndim))
		      
                      Bidotvj = DOT_PRODUCT(Brhoi(1:ndim),velj(1:ndim))
		      Bjdotvj = DOT_PRODUCT(Brhoj(1:ndim),velj(1:ndim))
		      
! KNOWN BUG HERE (not quite consistent derivation)		      
!	    	      prvaniso = prvaniso 		&
!		      - Rjoe*(- Bidotvj*projBrhoi*grkerni	&
!		              - Bjdotvi*projBrhoj*grkernj)

 		     prvanisoi = -Rjoe*(-Bidotvi*projBrhoi*grkerni &
		                        -Bjdotvi*projBrhoj*grkernj)
		     prvanisoj = -Rjoe*(-Bidotvj*projBrhoi*grkerni &
		     			-Bjdotvj*projBrhoj*grkernj)
		   ENDIF	   
! (dissipation term)
		   IF (ndimV.GT.ndim) THEN
!		      v2i = DOT_PRODUCT(veli,veli)	! total kinetic energy
!		      v2j = DOT_PRODUCT(velj,velj)
		      v2i = DOT_PRODUCT(veli,vunit)**2.	! along velocity line
		      v2j = DOT_PRODUCT(velj,vunit)**2.	! (use in 1.5D)
		   ELSE
		      v2i = DOT_PRODUCT(veli,dr)**2	! energy along line
		      v2j = DOT_PRODUCT(velj,dr)**2	! of sight		   
		   ENDIF
		   B2i = (DOT_PRODUCT(Bi,Bi) - DOT_PRODUCT(Bi,dr)**2)	!  "   " 
		   B2j = (DOT_PRODUCT(Bj,Bj) - DOT_PRODUCT(Bj,dr)**2)
   		   enj = gconst*uu(j) + 0.5*v2j + 0.5*B2j*rhoav1
		   eni = gconst*uu(i) + 0.5*v2i + 0.5*B2i*rhoav1
		   ediff = eni - enj 
		   IF (dvdotr.LT.0) THEN ! bug was in line below
		   qdiff = -0.5*alphaav*vsig*ediff*grkern*rhoav1
		   ELSE	! this is just the MHD bits if applying Bvisc everywhere
		   qdiff = 0.
		   ENDIF

		   dendt(i) = dendt(i) - pmassj*(projprv+prvaniso+prvanisoi+qdiff)
		   dendt(j) = dendt(j) + pmassi*(projprv+prvaniso+prvanisoj+qdiff)     		
		ELSE	! default is 'symmetric' form in Monaghan (1992)
	           dudt(i) = dudt(i) + 0.5*pmassj*prterm*dvdotr
	           dudt(j) = dudt(j) + 0.5*pmassi*prterm*dvdotr 
		ENDIF
!
!--calculate XSPH term for moving the particles
!
	        IF (ixsph.EQ.1) THEN
		   xsphterm(:,i) = xsphterm(:,i) - pmassj*dvel(:)*rhoav1*wab
		   xsphterm(:,j) = xsphterm(:,j) + pmassi*dvel(:)*rhoav1*wab
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
!--calculate gravitational force on all the particles
!
! IF (igravity.NE.0) CALL direct_sum_poisson(x,pmass,poten,fgrav,ntotal)

 DO i=1,npart
!
!--subtract external forces
!
    IF (itoystar.EQ.1) force(1:ndim,i) = force(1:ndim,i) - x(1:ndim,i)
!
!--add gravitational force
!
    IF (igravity.NE.0) force(1:ndim,i) = force(1:ndim,i) + fgrav(:,i)
!
!--damp force if appropriate
!
    IF (damp.NE.0.) force(:,i) = force(:,i) - damp*vel(:,i)
!
!--do the divisions by rho etc (this is for speed - so calculations are not
!  done multiple times within the loop)
!
    rho1i = 1./rho(i)
    IF (imhd.NE.0) THEN
       divB(i) = divB(i)*rho1i		!*rhoi
       IF ((imhd.EQ.1).OR.(imhd.EQ.3).OR.(imhd.EQ.11)) THEN
          dBfielddt(:,i) = dBfielddt(:,i)*rho1i
       ENDIF
    ENDIF    
!
!--calculate time derivative of the smoothing length
!
    IF (ihvar.EQ.2) THEN
       dhdt(i) = -hh(i)/(ndim*rho(i))*drhodt(i)
    ELSE
       dhdt(i) = 0.    
    ENDIF
!
!--calculate time derivative of alpha (artificial dissipation coefficient)
!  see Morris and Monaghan (1997)
!     
    IF (iavlim.EQ.1) THEN
!       IF (ndimV.GT.ndim) THEN
!          print*,' source ',i,' = ',avsource(i)*rho1i,drhodt(i)*rho1i
!          source = max(avsource(i)*rho1i,0.0)
!       ELSE    
          source = MAX(drhodt(i)*rho1i,0.0)    ! source term is div v
!         source = MAX(drhodt(i)*rho1i*(1.5-alpha(i)),0.0)    ! source term is div v
!       ENDIF
       valfven2i = 0.
       IF (imhd.GE.11) THEN
          valfven2i = DOT_PRODUCT(Bfield(:,i),Bfield(:,i))*rho1i
       ELSEIF (imhd.NE.0) THEN
          valfven2i = DOT_PRODUCT(Bfield(:,i),Bfield(:,i))*rho(i)       
       ENDIF
       vsig = SQRT(spsound(i)**2. + valfven2i) 	! approximate vsig only
       tdecay = hh(i)/(avconst*vsig)	! decay time (use vsig)
       daldt(i) = (alphamin - alpha(i))/tdecay + avfact*source
    ENDIF
      
 ENDDO
     
 DEALLOCATE ( avsource )
 IF (trace) WRITE(iprint,*) ' Exiting subroutine get_rates'
      
 RETURN
 END SUBROUTINE get_rates
