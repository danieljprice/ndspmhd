!!-----------------------------------------------------------------
!! Reads parameters for the run from the input file
!!-----------------------------------------------------------------

SUBROUTINE read_infile(infile)
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE artvi
 USE eos
 USE options
 USE setup_params
 USE timestep
 USE xsph
 USE anticlumping
!
!--define local variables
!      
 IMPLICIT NONE
 CHARACTER(LEN=*), INTENT(IN) :: infile 
 CHARACTER(LEN=1) :: ians  
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine read_infile'
              
 OPEN(UNIT=iread,ERR=999,FILE=infile,STATUS='old',FORM='formatted')
  READ(iread,*) psep
  READ(iread,*) tmax,tout,nmax,nout
  READ(iread,*) gamma
  READ(iread,*) iener,gconst,polyk
  READ(iread,*) icty,ndirect
  READ(iread,*) iprterm
  READ(iread,*) iav,alphamin,beta
  READ(iread,*) iavlim,avdecayconst
  READ(iread,*) ikernav
  READ(iread,*) ihvar,hfact
  READ(iread,*) idumpghost
  READ(iread,*) imhd,imagforce
  READ(iread,*) idivBzero
  READ(iread,*) ianticlump,eps,neps
  READ(iread,*) ixsph,xsphfac
 CLOSE(UNIT=iread)
!
!--check options for possible errors
!      
 IF (psep.LT.1.e-5) WRITE(iprint,100) 'psep < 1.e-5'
 IF (tout.GT.tmax) WRITE(iprint,100) 'no output tout > tmax'
 IF (nout.GT.nmax) WRITE(iprint,100) 'no output nout > nmax'
 IF (gamma.LT.1.) WRITE(iprint,100) 'gamma < 1.0 '
 IF ((iener.EQ.3).AND.(gconst.LT.0.)) THEN
    WRITE(iprint,100) 'gconst < 0.'
 ELSEIF ((iener.EQ.0).AND.(polyk.LT.0.)) THEN
    WRITE(iprint,100) 'polyk < 0.'      
 ENDIF
 IF ((iav.NE.0).AND.(alphamin.LT.0. .OR. beta.LT.0.) ) THEN
    WRITE(iprint,100) 'AV alpha or beta < 0.'      
 ENDIF
 IF ((iavlim.GT.0).AND.(alphamin.GE.1.)) THEN
    WRITE(iprint,100) 'using AV limiter, but alphamin set > 1.0'
 ENDIF
 IF ((iavlim.GT.0).AND.((avdecayconst.LE.0.01).OR.(avdecayconst.GT.0.5))) THEN
    WRITE(iprint,100) 'AV decay constant not in range 0.01-0.5'
 ENDIF     
 IF ((ikernav.LE.0).OR.(ikernav.GT.3)) THEN
    WRITE(iprint,100) 'kernel averaging not set (ikernav)'
 ENDIF
 IF ((hfact.LE.1.0).OR.(hfact.GT.2.0)) THEN
    WRITE(iprint,100) 'hfact too low/high (1.0 < hfact < 2.0)'
 ENDIF

100   FORMAT(/' read_infile: warning: ',a)
 RETURN
            
999   WRITE(iprint,1000) infile      
1000  FORMAT (' Input file ',a20,' not found')
      WRITE(iprint,*) ' Would you like to create one with default options?'
      READ*,ians
      IF (ians.EQ.'y'.OR.ians.EQ.'Y') THEN
!
!--set default options
!         
	 psep = 0.01
	 tmax = 10.0
	 tout = 1.0
	 nmax = 1000000
	 nout = -1
	 gamma = 5./3.
	 iener = 2
	 gconst = 0.0
	 polyk = 0.6
	 icty = 0
	 ndirect = nmax
	 iprterm = 0
	 iav = 2
	 alphamin = 0.1
	 beta = 2.0
	 iavlim = 1
	 avdecayconst = 0.1
	 ikernav = 2
	 ihvar = 2
	 hfact = 1.5
	 idumpghost = 0
	 imhd = 11
	 imagforce = 2
	 idivBzero = 0
	 ianticlump = 1
	 eps = 0.8
	 neps = 6
	 ixsph = 0
	 xsphfac = 0.5
	 
	 CALL write_infile(infile)
      ENDIF

      STOP 'exiting...'
      
END SUBROUTINE read_infile
