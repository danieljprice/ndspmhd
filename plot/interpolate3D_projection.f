!--------------------------------------------------------------------------
!     program to interpolate from particle data to even grid of pixels
!
!     The data is smoothed using the SPH summation interpolant,
!     that is, we compute the smoothed array according to
!
!     datsmooth(pixel) = sum_b m_b dat_b/rho_b W(r-r_b, h_b)
! 
!     where _b is the quantity at the neighbouring particle b and
!     W is the smoothing kernel, for which we use the usual cubic spline
!
!     ** In this version 3D data is interpolated to a 2D grid by use of an
!     ** integrated form of the kernel (that is W_ab in this case is
!     ** the integral through the 3D kernel to give a 2D kernel)
!     ** This results in a column density map of the interpolated quantity
!     ** From a similar routine by Matthew Bate.
!
!     Input: particle coordinates  : x,y   (npart)
!            particle masses       : pmass (npart)
!            density on particles  : rho   (npart) - must be computed separately
!            smoothing lengths     : hh    (npart) - could be computed from density
!            scalar data to smooth : dat   (npart)
!
!     Output: smoothed data 	   : datsmooth (npixx,npixy)
!
!     Written by Daniel Price September 2003
!--------------------------------------------------------------------------

      SUBROUTINE interpolate3D_projection(
     &    x,y,pmass,rho,hh,dat,npart,
     &    xmin,ymin,datsmooth,npixx,npixy,pixwidth)
      
      USE column
      IMPLICIT NONE
      REAL, PARAMETER :: pi = 3.1415926536      
      INTEGER, INTENT(IN) :: npart,npixx,npixy
      REAL, INTENT(IN), DIMENSION(npart) :: x,y,pmass,rho,hh,dat
      REAL, INTENT(IN) :: xmin,ymin,pixwidth
      REAL, INTENT(OUT), DIMENSION(npixx,npixy) :: datsmooth

      INTEGER :: i,j,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
      INTEGER :: index, index1
      REAL :: hi,hi1,h2,radkern,qq,wab,rab,const
      REAL :: term,dx,dy,xpix,ypix
      REAL :: dxx,dwdx
      REAL :: dmaxcoltable

      datsmooth = 0.
      term = 0.
      dmaxcoltable = 1./REAL(maxcoltable)
      PRINT*,'projecting from particles to pixels...'
!
!--loop over particles
!      
      DO i=1,npart
!
!--set kernel related quantities
!
         hi = hh(i)
	 hi1 = 1./hi
	 h2 = hi*hi
	 radkern = 2.*hi	! radius of the smoothing kernel
!         const = 10./(7.*pi*h2)	! normalisation constant
         const = hi1*hi1
	 IF (rho(i).NE.0.) term = const*pmass(i)*dat(i)/rho(i) 
!
!--for each particle work out which pixels it contributes to
!               
	 ipixmin = INT((x(i) - radkern - xmin)/pixwidth)
	 jpixmin = INT((y(i) - radkern - ymin)/pixwidth)
	 ipixmax = INT((x(i) + radkern - xmin)/pixwidth)
	 jpixmax = INT((y(i) + radkern - ymin)/pixwidth)
	 
         IF (ipixmin.LT.1) ipixmin = 1	! make sure they only contribute
	 IF (jpixmin.LT.1) jpixmin = 1  ! to pixels in the image
	 IF (ipixmax.GT.npixx) ipixmax = npixx
	 IF (jpixmax.GT.npixy) jpixmax = npixy
!
!--loop over pixels, adding the contribution from this particle
!
	    DO jpix = jpixmin,jpixmax
	       DO ipix = ipixmin,ipixmax

		  ypix = ymin + (jpix)*pixwidth - 0.5*pixwidth
		  xpix = xmin + (ipix)*pixwidth - 0.5*pixwidth
		  dy = ypix - y(i)
		  dx = xpix - x(i)
		  rab = SQRT(dx**2 + dy**2)
		  qq = rab*hi1
!
!--SPH kernel - integral through cubic spline
!  interpolate from a pre-calculated table
!		     
!!--find nearest index in table
                  index = INT(0.5*qq*maxcoltable) + 1
		  IF (index.LT.maxcoltable) THEN
		   index1 = index + 1
!--find increment along from this index
  		   dxx = qq - index*maxcoltable
!--find gradient		  
		   dwdx = (coltable(index1) - coltable(index))*dmaxcoltable
!--compute value of integrated kernel
		   wab = coltable(index) + dwdx*dxx
!
!--calculate data value at this pixel using the summation interpolant
!		  
		   datsmooth(ipix,jpix) = datsmooth(ipix,jpix) + term*wab		          
                  ENDIF
      
               ENDDO
	    ENDDO

      ENDDO
      
      RETURN
      
      END SUBROUTINE interpolate3D_projection
