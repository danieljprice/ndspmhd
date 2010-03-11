      SUBROUTINE cp_distribute(rmin, rmax, d0, npart, x, y, z, idim)
c************************************************************
c                                                           *
c  This subroutine positions particles on a close packed    *
c  lattice in three dimensions                              *
c  The lattice is centred on the origin and trimmed to      *
c  give a spherical distribution rmin < r < rmax            *
c                                                           *
c  D. Price University of Exeter 2005                       *
c  hacked from an original routine by Matthew Bate          *
c                                                           *
c************************************************************

      IMPLICIT NONE ! those sweet, sweet words
      INTEGER np,npart,idim
      REAL rmin,rmax,d0
      REAL x(idim),y(idim),z(idim)

      INTEGER i,k,l,m,nx,ny,nz
      INTEGER jy,jz
      REAL xmin,xmax,ymin,ymax,zmin,zmax
      REAL facx,facy,facz,delx,dely,stepx,stepy,stepz
      REAL deltax,deltay,deltaz
      REAL xstart,ystart,zstart,xi,yi,zi,rr2
      REAL xcentre,ycentre,zcentre

c
c--Set uniform particle distribution, centred at the origin
c
      xcentre = 0.
      ycentre = 0.
      zcentre = 0.
      xmin = -rmax
      ymin = -rmax
      zmin = -rmax
      xmax = rmax
      ymax = rmax
      zmax = rmax

      facx = 1.
      facy = SQRT(3./4.)
      facz = SQRT(6.)/3.
      delx = 0.5*d0
      dely = d0*SQRT(3.)/6.

      deltax = xmax - xmin
      deltay = ymax - ymin
      deltaz = zmax - zmin
      nx = INT(deltax/(facx*d0)) + 1
      ny = INT(deltay/(facy*d0)) + 2
      nz = INT(deltaz/(facz*d0)) + 1
      np = nx*ny*nz
      WRITE(*,*) 'nx,ny,nz = ',nx,ny,nz
      
      stepx = d0*facx
      stepy = d0*facy
      stepz = d0*facz
c
c--set the limits so that the particles are 
c  exactly centred on the origin
c
      xmin = xcentre - (nx-1)/2*stepx
      xmax = xcentre + (nx-1)/2*stepx
      ymin = ycentre - (ny-1)/2*stepy
      ymax = ycentre + (ny-1)/2*stepy
      zmin = zcentre - (nz-1)/2*stepz
      zmax = zcentre + (nz-1)/2*stepz

      k = 0
      l = 1
      m = 1
      npart = 0
      DO i = 1, np
         k = k + 1
         IF (k.GT.nx) THEN
            k = 1
            l = l + 1
            IF (l.GT.ny) THEN
               l = 1
               m = m + 1
               IF (m.GT.nz) THEN
                  k = 1
                  l = 1
                  nz = nz + 1
               ENDIF
            ENDIF
         ENDIF
            
         xstart = xmin
	 ystart = ymin
	 zstart = zmin
         jy = MOD(l, 2)
         jz = MOD(m, 3)
	 IF (jz.EQ.0) THEN	! 3rd layer
	    ystart = ystart + 2.*dely
	    IF (jy.EQ.0) xstart = xstart + delx		   
	 ELSEIF (jz.EQ.2) THEN	! 2nd layer	  
	    ystart = ystart + dely
	    IF (jy.EQ.1) xstart = xstart + delx
	 ELSEIF (jy.EQ.0) THEN    ! first layer	  
	    xstart = xstart + delx
	 ENDIF

         xi = xstart + FLOAT(k - 1)*stepx
         yi = ystart + FLOAT(l - 1)*stepy
         zi = zstart + FLOAT(m - 1)*stepz 
c
c--trim to fit radius. Do not allow particles to have *exactly* rmax
c  (this stops round-off error from giving non-zero centre of mass)
c
         rr2 = xi*xi + yi*yi + zi*zi
         IF (rr2.le.(rmax**2- 0.01*d0**2) .AND.
     &       rr2.ge.rmin**2) THEN
            npart = npart + 1
            IF (npart.GT.idim) THEN
               WRITE(*,*) 'ERROR! arrays too small in cp_distribute!'
               STOP
            ENDIF
            x(npart) = xi
            y(npart) = yi
            z(npart) = zi
         ENDIF

      ENDDO

      RETURN
      END
