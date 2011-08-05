cc----------------------------------------------------
cc   A Blast wave is set up as in Colella and Woodward
cc   and calculated by Marti et al phys.rev.D 1991
cc   Density is constant so particles are equispaced
cc   Initidist velocity is zero. Pressure jumps.
cc---------------------------------------------------

      SUBROUTINE setup
c
c--include relevant global variables
c
      INCLUDE 'COMMONS/dimen_mhd'
      INCLUDE 'COMMONS/bound'
      INCLUDE 'COMMONS/eos'      
      INCLUDE 'COMMONS/options'
      INCLUDE 'COMMONS/part'      
      INCLUDE 'COMMONS/part_in'      
      INCLUDE 'COMMONS/setup_params'      
c
c--define local variables
c      
      REAL dist,xcentre,delta,vxi,term,term2,rhozero,massp,gam1
      REAL exx,uuleft,uuright
c      REAL xmin,xmax	! local since they are repositioned in initialise

      ibound = 2
c      nbpts = 10
      xmin = -0.5
      xmax = 0.5
      xcentre  = 0.0
      gam1 = gamma - 1.
      dist  = 0.5*psep
      rhozero =  1.0   
      uuright = 0.01/(rhozero*gam1)        ! thermal energy on the right
      uuleft = 1000./(rhozero*gam1)       ! thermal energy on the left
      massp = rhozero*psep                 ! the mass per sph particle
      hfact = hfact*massp	! this is so h = hfact/rho (used in sum_density)
      npart = 0
c
c--mhd
c
c      Bin(:,:) = 0. 
      i = 1
      xin(1) = xmin + 0.5*psep       

      DO WHILE (xin(i).lt.xmax)       
         npart = npart + 1
         xin(i) = (i-1.)*psep - 0.5
         delta = (xin(i) - xcentre)/dist
         rhoin(i) = rhozero
         
	 IF (delta.gt.20.) THEN
           enin(i) = uuright
         ELSEIF (delta.lt.-20.) THEN
           enin(i) = uuleft
         ELSE
           exx = exp(delta)
           enin(i) = (uuleft + uuright*exx)/(1 + exx)
         ENDIF
	 
         hhin(i) = hfact/rhozero
         velin(:,i) = 0.
         uuin(i) = enin(i) - 0.5*DOT_PRODUCT(velin(:,i),velin(:,i))
         pmass(i) = massp
         i = i + 1
         xin(i) = (i-1)*psep - 0.5
      ENDDO
      
      ntotal = npart
              
      RETURN
      END
