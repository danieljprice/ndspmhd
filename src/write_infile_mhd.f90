!!-----------------------------------------------------------------
!! Writes an input file
!!-----------------------------------------------------------------

SUBROUTINE write_infile(infile)
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE artvi
 USE eos
 USE options
 USE polyconst
 USE setup_params
 USE timestep
 USE xsph
 USE anticlumping
!
!--define local variables
!      
 IMPLICIT NONE
 CHARACTER(LEN=*), INTENT(IN) :: infile   
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine write_infile'
              
 OPEN(UNIT=iread,ERR=999,FILE=infile,STATUS='replace',FORM='formatted')
  WRITE(iread,10) psep
  WRITE(iread,20) tmax,tout,nmax,nout
  WRITE(iread,30) gamma
  WRITE(iread,40) iener,udiss_frac,alphaBmin,polyk
  WRITE(iread,50) icty,ndirect
  WRITE(iread,60) iprterm
  WRITE(iread,70) iav,alphamin,beta
  WRITE(iread,80) iavlim,avdecayconst
  WRITE(iread,90) ikernav
  WRITE(iread,100) ihvar,hfact
  WRITE(iread,110) idumpghost
  WRITE(iread,120) imhd,imagforce
  WRITE(iread,130) idivBzero,psidecayfact
  WRITE(iread,140) ianticlump,eps,neps
  WRITE(iread,150) ixsph,xsphfac
  WRITE(iread,160) igravity
  WRITE(iread,170) damp
  WRITE(iread,180) ikernel
 CLOSE(UNIT=iread)

10 FORMAT(f14.10,22x,'! particle separation')
20 FORMAT(f7.3,2x,f7.3,2x,i9,1x,i5,3x,'! tmax, tout, nmax, nout')
30 FORMAT(f14.12,22x,'! gamma ')
40 FORMAT(i1,2x,f5.3,2x,f5.3,2x,f5.3,21x,'! type of energy equation, udiss_frac, alphaBmin(for iner=3), polyk(for iener=0)')
50 FORMAT(i1,2x,i9,24x,'! type of cty equation (0:direct sum 1:time deriv 2:alt time deriv)')
60 FORMAT(i1,35x,'! type of pressure term (0:normal 1:Pa+Pb/rhoa*rhob 2:Hernquist/Katz )')
70 FORMAT(i1,2x,f5.3,2x,f5.3,21x,'! artificial viscosity type, alpha(min), beta(0:off 1:M92 2:M97)')
80 FORMAT(i1,2x,f5.3,28x,'! use Morris & Monaghan av limiter, constant for this(0.1-0.2)')
90 FORMAT(i1,35x,'! type of kernel averaging (1:average h, 2:average grad Wab 3:Springel/Hernquist)')
100 FORMAT(i1,2x,f5.3,28x,'! hvariable (0:no 1:yes) initial h factor')
110 FORMAT(i1,35x,'! dump ghost particles? (0: no 1: yes)')
120 FORMAT(i2,4x,i1,29x,'! Magnetic field (0:off 1:on) and force algorithm(1:vector 2:tensor)')
130 FORMAT(i2,2x,f5.3,28x,'! divergence correction method (0:none 1:projection 2: kernel 3: Borve 4: Joe 5: dans kernel)')
140 FORMAT(i1,2x,f5.3,2x,i2,24x,'! anticlumping term (0:off 1:on), eps, power')
150 FORMAT(i1,2x,f5.3,28x,'! use XSPH, parameter')
160 FORMAT(i1,35x,'! self-gravity')
170 FORMAT(f7.4,28x,'! artificial damping (0.0 or few percent)')
180 FORMAT(i1,35x,'! kernel type (0: cubic spline, 3:quintic)')

 WRITE(iprint,200) infile
200 FORMAT (' Input file ',a20,' created successfully')

 RETURN
      
999   STOP 'Error creating input file, exiting...'
      
END SUBROUTINE write_infile
