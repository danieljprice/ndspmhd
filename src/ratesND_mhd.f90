!!--------------------------------------------------------------------
!! Computes the rates of change of the conserved variables
!! (forces, energy etc)
!! This is the core of the SPH algorithm
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
 INTEGER, DIMENSION(3**ndim) :: neighcell
!
!  (particle properties - local copies and composites)
!
 REAL :: rij,rij2
 REAL :: rhoi,rho1i,rho2i,rho21i,rhoj,rho1j,rho2j,rho21j,rhoav,rhoav1,rhoij
 REAL :: pmassi,pmassj
 REAL :: pri,prj,Prho2i,Prho2j,prterm
 REAL :: hi,hi1,hj,hj1,hi2,hi21,hj2,hj21,hi3,hj3,hav,h2,h3,hav1,h21
 REAL :: drhodti,drhodtj
 REAL :: hfacwab,hfacwabi,hfacwabj,hfacgrkern,hfacgrkerni,hfacgrkernj
 REAL, DIMENSION(ndim) :: dx
!
!  (velocity)
!      
 REAL, DIMENSION(ndimV) :: veli,velj,dvel
 REAL, DIMENSION(ndimV) :: vterm,prvterm
 REAL, DIMENSION(ndimV) :: dr
 REAL :: dvdotr,projprv,projvterm,v2i,v2j
!
!  (mhd)
!     
 REAL, DIMENSION(ndimB) :: Brhoi,Brhoj,Bi,Bj,dB
 REAL, DIMENSION(ndimB) :: faniso,fmagi,fmagj,fcorr
 REAL, DIMENSION(ndimB) :: curlBi
 REAL :: fiso,divBonrho
 REAL :: valfven2i,valfven2j
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
 REAL :: vissv,vissB,vissu,vsigii,vsigjj
 REAL :: avterm
!
!  (alternative forms)
!
 REAL, DIMENSION(:), ALLOCATABLE :: phi
 REAL :: phii,phii1,phii_on_phij,phij_on_phii
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
 REAL :: alphai,alphaav,source,tdecay1    
!
!  (time step criteria)
!      
 REAL :: vsigdtc,vmag,zero,fhmax, fonh, forcemag
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
!--allocate memory for local arrays
!
 nlistdim = ntotal
 ALLOCATE ( listneigh(nlistdim),STAT=ierr )
 IF (ierr.NE.0) WRITE(iprint,*) ' Error allocating neighbour list, ierr = ',ierr
 ALLOCATE ( phi(ntotal), STAT=ierr )
 IF (ierr.NE.0) WRITE(iprint,*) ' Error allocating phi, ierr = ',ierr 
!
!--initialise quantities
!      
 dtcourant = 1.e6  
 zero = 1.e-10

 DO i=1,ntotal	! using ntotal just makes sure they are zero for ghosts
  force(:,i) = 0.0
  drhodt(i) = 0.0
  dudt(i) = 0.0
  dendt(i) = 0.0
  dBconsdt(:,i) = 0.0
  fmag(:,i) = 0.0
  divB(i) = 0.0
  curlB(:,i) = 0.0
  xsphterm(:,i) = 0.0
 ENDDO
!
!  set alternative forms for the SPH equations here
!  phi can be any scalar variable  
!
 SELECT CASE(ialtform)
    CASE(1)		! this gives the (P_a + P_b) / (rho_a rho_b) form
       phi(1:ntotal) = rho(1:ntotal)
    CASE(2)		! this gives the HK89 form SQRT(Pa Pb)/rhoa rhob
       phi(1:ntotal) = SQRT(pr(1:ntotal))/rho(1:ntotal)    
    CASE DEFAULT	! this gives the usual continuity, momentum and induction eqns
       phi(1:ntotal) = 1.0 
 END SELECT
!
!--calculate kernel for the MHD anticlumping term
!
 IF (ianticlump.EQ.1) THEN
  q2joe = (1./hfact)**2	! 1/hfact is initial particle spacing in units of h 
  CALL interpolate_kernel(q2joe,wabjoei,grkerni)
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

!    PRINT*,'Doing particle ',i,x(:,i),valfven2i,rho(i),hh(i)
       idone = idone + 1
       rhoi = rho(i)
       rho2i = rhoi*rhoi
       rhoi5 = SQRT(rhoi)
       rho1i = 1./rhoi
       rho21i = rho1i*rho1i       
       pri = pr(i)
       Prho2i = pr(i)*rho21i
       spsoundi = spsound(i)
       veli(:) = vel(:,i)
       pmassi = pmass(i)
       alphai = alpha(i)
       phii = phi(i)
       phii1 = 1./phii
       IF (imhd.GE.11) THEN	! if mag field variable is B
          Bi(:) = Bcons(:,i)
          Brhoi(:) = Bi(:)*rho1i
       ELSEIF (imhd.GT.0) THEN	! if mag field variable is B/rho
          Brhoi(:) = Bcons(:,i)
          Bi(:) = Bfield(:,i)
       ENDIF
! mhd definitions
       Brho2i = DOT_PRODUCT(Brhoi,Brhoi)
       valfven2i = Brho2i*rhoi
              
       gradhi = 1./(1. - gradh(i))
!       IF (gradhi.LE.0.5) THEN
!          WRITE(iprint,*) 'Error in grad h terms, part ',i,gradhi
!       ENDIF   
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
!	     print*,' ... neighbour, h=',j,hh(j),rho(j),x(:,j)
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
	        Prho2j = pr(j)*rho21j
		spsoundj = spsound(j)
	        pmassj = pmass(j)
		alphaav = 0.5*(alphai + alpha(j))
		phii_on_phij = phii/phi(j)
		phij_on_phii = phi(j)*phii1 		
!  (mhd definitions)
  	        IF (imhd.GE.11) THEN	! if B is mag field variable
		 Bj(:) = Bcons(:,j)
		 Brhoj(:) = Bj(:)*rho1j
		ELSEIF (imhd.NE.0) THEN	! if B/rho is mag field variable
		 Brhoj(:) = Bcons(:,j)
		 Bj(:) = Bfield(:,j)		    	          
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
                ENDIF
!
!--artificial viscosity terms
!
	        rhoav = 0.5*(rhoi + rhoj)
		rhoav1 = 1./rhoav
!--maximum velocity for timestep control
		vmag = SQRT(DOT_PRODUCT(dvel,dvel))
!		vsigdtc = vmag + spsoundi + spsoundj  	&
!			       + sqrt(valfven2i) + sqrt(valfven2j)
		vsig = 0.
		viss = 0.
		visc = 0.
		IF (vmag.EQ.0.) vmag = 1.e-6
		vunit(:) = abs(dvel(:))/vmag
!
!--calculate maximum signal velocity (this is used for the timestep control
!  and also in the artificial viscosity)
!
		IF (dvdotr.LT.0.) THEN
		   vunit(:) = -dvel(:)/vmag	! unit vector in direction of velocity
		   !vunit = dr
		   viss = abs(dvdotr)		    
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
		vsig = vsigi + vsigj + beta*viss
	
!		vsigdtc is the signal velocity used in the timestep control
                vsigdtc = vsigi + vsigj + beta*abs(dvdotr)

		IF (dvdotr.LT.0. .AND. iav.NE.0) THEN	! only for approaching particles
                   IF (iav.EQ.1) THEN
		   avterm = 0.5*alphaav*vsig*rhoav1
!--viscosity (kinetic energy term)
                   visc = avterm*viss*grkern
	           vissv = -0.5*(DOT_PRODUCT(veli,dr) 	&
		               - DOT_PRODUCT(velj,dr))**2		    
!--ohmic dissipation (magnetic energy term)
		   Bvisc(:) = (dB(:) - dr(:)*projdB)*rhoav1
		   dBdtvisc(:) = avterm*Bdiss_frac*Bvisc(:)
		   vissB = -0.5*Bdiss_frac*(DOT_PRODUCT(dB,dB)-projdB**2)*rhoav1
!--vissu is the dissipation energy from thermal conductivity
		   vissu = udiss_frac*(uu(i) - uu(j))
!--envisc is the total contribution to the thermal energy equation
		   envisc = avterm*(vissv+vissB)*grkern!*abs(rx)
	           uvisc = avterm*vissu*grkern		   
		   ELSEIF (iav.EQ.2) THEN
		   STOP 'crap av choice'
!--viscosity (kinetic energy term)
		   vsigii = vsigi + 0.5*beta*viss
		   vsigjj = vsigj + 0.5*beta*viss   
		   vsig = vsigii/rhoi*grkerni + vsigjj/rhoj*grkernj                
		   visc = 0.5*alphaav*vsig*viss
		   vissv = -0.5*(DOT_PRODUCT(veli,dr) 	&
		              - DOT_PRODUCT(velj,dr))**2		    
!--ohmic dissipation (magnetic energy term)
		   Bvisc(:) = dB(:) - dr(:)*projdB
		   dBdtvisc(:) = 0.5*alphaav*vsig*Bvisc(:)*rhoav1
		   vissB = -0.5*(DOT_PRODUCT(dB,dB)-projdB**2)*rhoav1
!--vissu is the dissipation energy from thermal conductivity
		   vissu = udiss_frac*(uu(i) - uu(j))
!--envisc is the total contribution to the thermal energy equation
		   envisc = 0.5*alphaav*vsig*(vissv+vissB) !*abs(rx)
	           uvisc = 0.5*alphaav*vsig*vissu		   
		   
		   ENDIF
		ELSE
		   visc = 0.0
		   dBdtvisc(:) = 0.
		   envisc = 0.
		   uvisc = 0.
	        ENDIF	
!
!--time step control (courant and viscous)
!
	        IF (vsigdtc.GT.zero) dtcourant = min(dtcourant,0.8*hav/vsigdtc)
!
!--pressure term (generalised form)
!
       	        IF (ialtform.GE.0) THEN
                prterm = phii_on_phij*Prho2i*grkerni 		&
     		       + phij_on_phii*Prho2j*grkernj	   	   
                ENDIF
!
!--add pressure and viscosity terms to force (equation of motion)
!
		IF (ndimV.GT.ndim) THEN 
		  force(:,i) = force(:,i)	&
		             - pmassj*(prterm*dr(:)+visc*vunit(:))
		  force(:,j) = force(:,j)	&
		             + pmassi*(prterm*dr(:)+visc*vunit(:))  

		ELSE
		  force(:,i) = force(:,i) - pmassj*(prterm+visc)*dr(:)
		  force(:,j) = force(:,j) + pmassi*(prterm+visc)*dr(:)		
                ENDIF	
!
!--time derivative of density (continuity equation) in generalised form
!  compute this even if direct sum - gives divv for art vis.
!
!	        drhodti = phii_on_phij*pmassj*dvdotr*grkerni
!		drhodtj = phij_on_phii*pmassi*dvdotr*grkernj
	        drhodti = pmassj*dvdotr*grkerni
		drhodtj = pmassi*dvdotr*grkernj

				
		drhodt(i) = drhodt(i) + drhodti
		drhodt(j) = drhodt(j) + drhodtj
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
		  		  			
		  ELSEIF (imagforce.EQ.2) THEN	! tensor formalism in generalised form
!
!--isotropic mag force (pressure gradient)
!		  
	          fiso = 0.5*(Brho2i*phij_on_phii*grkerni &
		            + Brho2j*phii_on_phij*grkernj)
!
!--anisotropic force (including variable smoothing length terms)
!
		  faniso(:) = Brhoi(:)*projBrhoi*phij_on_phii*grkerni  &
		            + Brhoj(:)*projBrhoj*phii_on_phij*grkernj		   
!
!--Joe's correction term (preserves momentum conservation)
!  *** in 1D only do this in the x-direction?
		  IF (ianticlump.EQ.1) THEN
		     wabjoe = wabjoei*hfacwab
		     Rjoe = 0.5*eps*(wab/wabjoe)**neps
!		     IF (Rjoe.GT.0.1) PRINT*,'Rjoe = ',Rjoe,i,j,wab,wabjoe
		     faniso(:) = faniso(:) - Rjoe*faniso(:)  
!		     faniso(1:ndim) = faniso(1:ndim) - Rjoe*faniso(1:ndim)  
		  ENDIF
!
!--add contributions to magnetic force
!	  
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
!--time derivative of magnetic field (divide by rho later) - in generalised form
!
!  (evolving B/rho)
		  IF (imhd.EQ.1) THEN
		     dBconsdt(:,i) = dBconsdt(:,i)		&
                   - phii_on_phij*pmassj*(dvel(:)*projBrhoi)*grkerni 
		     dBconsdt(:,j) = dBconsdt(:,j) 		&
     		   - phij_on_phii*pmassi*(dvel(:)*projBrhoj)*grkernj

		     IF (iav.NE.0) THEN		! add dissipative term ***check
	                dBconsdt(:,i) = dBconsdt(:,i)		&
                         + rhoi*pmassj*dBdtvisc(:)*grkern		   
		        dBconsdt(:,j) = dBconsdt(:,j)		&
                         - rhoj*pmassi*dBdtvisc(:)*grkern
	             ENDIF
!   (evolving B)
		  ELSEIF (imhd.EQ.11) THEN	! note divided by rho later		  
		     dBconsdt(:,i) = dBconsdt(:,i)		&
            + pmassj*(Bi(:)*dvdotr - dvel(:)*projBi)*grkerni   
		     dBconsdt(:,j) = dBconsdt(:,j) 		&
            + pmassi*(Bj(:)*dvdotr - dvel(:)*projBj)*grkernj
		  
		     IF (iav.NE.0) THEN		! add dissipative term
	                dBconsdt(:,i) = dBconsdt(:,i)		&
                         + rho2i*pmassj*dBdtvisc(:)*grkern		   
		        dBconsdt(:,j) = dBconsdt(:,j)		&
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
!--energy equation (total or thermal)
!               
		IF (iener.GE.3) THEN	! total energy/particle (generalised form)
		   Bidotvj = DOT_PRODUCT(Brhoi,velj)
		   Bjdotvi = DOT_PRODUCT(Brhoj,veli)
! (isotropic stress)
		   prvterm(:) = (Prho2i+0.5*Brho2i)*phii_on_phij*velj(:)*grkerni &
                              + (Prho2j+0.5*Brho2j)*phij_on_phii*veli(:)*grkernj
! (anisotropic stress)
		   prvaniso =  - Bidotvj*projBrhoi*phii_on_phij*grkerni	& 
		               - Bjdotvi*projBrhoj*phij_on_phii*grkernj
		   projprv = DOT_PRODUCT(prvterm,dr)		   
! (add source term for anticlumping term)		   
		   IF (ianticlump.EQ.1 .AND. imagforce.EQ.2) THEN
!		      prvaniso = prvaniso - Rjoe*prvaniso
! (if applied in x-direction only)
                      Bidotvi = DOT_PRODUCT(Brhoi(1:ndim),veli(1:ndim))
		      Bjdotvi = DOT_PRODUCT(Brhoj(1:ndim),veli(1:ndim))
		      
                      Bidotvj = DOT_PRODUCT(Brhoi(1:ndim),velj(1:ndim))
		      Bjdotvj = DOT_PRODUCT(Brhoj(1:ndim),velj(1:ndim))
		      
 		      prvanisoi = -Rjoe*(-Bidotvi*projBrhoi*phii_on_phij*grkerni &
		                         -Bjdotvi*projBrhoj*phij_on_phii*grkernj)
		      prvanisoj = -Rjoe*(-Bidotvj*projBrhoi*phii_on_phij*grkerni &
		     			 -Bjdotvj*projBrhoj*phij_on_phii*grkernj)
		   ENDIF	   
! (dissipation term)
		   v2i = DOT_PRODUCT(veli,dr)**2	! energy along line
		   v2j = DOT_PRODUCT(velj,dr)**2	! of sight		   
		   B2i = (DOT_PRODUCT(Bi,Bi) - DOT_PRODUCT(Bi,dr)**2)	!  "   " 
		   B2j = (DOT_PRODUCT(Bj,Bj) - DOT_PRODUCT(Bj,dr)**2)
   		   enj = udiss_frac*uu(j) + 0.5*v2j + Bdiss_frac*0.5*B2j*rhoav1
		   eni = udiss_frac*uu(i) + 0.5*v2i + Bdiss_frac*0.5*B2i*rhoav1
		   ediff = eni - enj 
		   IF (dvdotr.LT.0) THEN ! bug was in line below
		    IF (iav.EQ.1) THEN
		       qdiff = -avterm*ediff*grkern
		    ELSEIF (iav.EQ.2) THEN
		       qdiff = -0.5*alphaav*vsig*ediff
		    ENDIF
		   ELSE	! this is just the MHD bits if applying Bvisc everywhere
		   qdiff = 0.
		   ENDIF

		   dendt(i) = dendt(i) - pmassj*(projprv+prvaniso+prvanisoi+qdiff)
		   dendt(j) = dendt(j) + pmassi*(projprv+prvaniso+prvanisoj+qdiff)     		
		ENDIF
!
!--compute thermal energy derivative even if using total energy equation
!		
		IF (iener.NE.0) THEN	!  1st law of thermodynamics
                   dudt(i) = dudt(i) + Prho2i*drhodti + pmassj*(envisc + uvisc)
     		   dudt(j) = dudt(j) + Prho2j*drhodtj + pmassi*(envisc - uvisc)
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

 fhmax = 0.0
 dtforce = 1.e6
 
 DO i=1,npart
!
!--subtract external forces
!
    IF (iexternal_force.NE.0) CALL external_forces(i,iexternal_force)
!
!--add self-gravity force
!
    IF (igravity.NE.0) force(1:ndim,i) = force(1:ndim,i) + fgrav(1:ndim,i)
!
!--damp force if appropriate
!
    IF (damp.GT.0.) force(:,i) = force(:,i) - damp*vel(:,i)
!
!--do the divisions by rho etc (this is for speed - so calculations are not
!  done multiple times within the loop)
!
    rho1i = 1./rho(i)
    IF (imhd.NE.0) THEN
       divB(i) = divB(i)*rho1i		!*rhoi
       dBconsdt(:,i) = dBconsdt(:,i)*rho1i
    ENDIF    
!
!--calculate time derivative of the smoothing length
!
    IF (ihvar.EQ.2 .OR. ihvar.EQ.3) THEN
       dhdt(i) = -hh(i)/(ndim*rho(i))*drhodt(i)
    ELSE
       dhdt(i) = 0.    
    ENDIF
!
!--if using the thermal energy equation, set the energy derivative
!
    IF (iener.NE.3 .AND. iener.NE.0) THEN
       dendt(i) = dudt(i)
    ENDIF
!
!--calculate maximum force/h for the timestep condition
!  also check for errors in the force
!
    IF ( ANY(force(:,i).GT.1.e8)) THEN
       WRITE(iprint,*) 'rates: force ridiculous ',force(:,i)
       CALL quit
    ENDIF
    forcemag = SQRT(DOT_PRODUCT(force(:,i),force(:,i)))   
    fonh = forcemag/hh(i)
    IF (fonh.GT.fhmax) fhmax = fonh
!
!--calculate time derivative of alpha (artificial dissipation coefficient)
!  see Morris and Monaghan (1997)
!     
    IF (iavlim.NE.0) THEN
       IF (iavlim.EQ.2) THEN
          IF (idivBzero.EQ.1) THEN
             source = MAX(drhodt(i)*rho1i*(2.0-alpha(i)),abs(divB(i))*SQRT(rho1i),0.0)    ! source term is div v	 
	  ELSE	  
	     source = MAX(drhodt(i)*rho1i*(2.0-alpha(i)),0.0)    ! source term is div v
          ENDIF
       ELSE
          IF (idivBzero.EQ.1) THEN
             source = MAX(drhodt(i)*rho1i,abs(divB(i))*SQRT(rho1i),0.0)    ! source term is div v	 
	  ELSE          
	     source = MAX(drhodt(i)*rho1i,0.0)    ! source term is div v
          ENDIF
       ENDIF
       valfven2i = 0.
       IF (imhd.NE.0) valfven2i = DOT_PRODUCT(Bfield(:,i),Bfield(:,i))*rho1i
       vsig = SQRT(spsound(i)**2. + valfven2i) 	! approximate vsig only
       tdecay1 = (avconst*vsig)/hh(i)	! 1/decay time (use vsig)
       daldt(i) = (alphamin - alpha(i))*tdecay1 + avfact*source
    ELSE
       daldt(i) = 0.
    ENDIF
      
 ENDDO
!
!--calculate timestep constraint from the forces
!  dtforce is returned together with dtcourant to the main timestepping loop
!
 IF (fhmax.LE.0.) THEN
    WRITE(iprint,*) 'rates: fhmax <=0 :',fhmax
    CALL quit
 ENDIF    
 IF (dtforce.GT.0.0) dtforce = SQRT(1./fhmax)
     
 IF (ALLOCATED(listneigh)) DEALLOCATE(listneigh)
 IF (ALLOCATED(phi)) DEALLOCATE(phi)
 IF (trace) WRITE(iprint,*) ' Exiting subroutine get_rates'
      
 RETURN
 END SUBROUTINE get_rates
