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
 CHARACTER(LEN=1) :: ians  
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine read_infile'
              
 OPEN(UNIT=iread,ERR=999,FILE=infile,STATUS='old',FORM='formatted')
  READ(iread,*,ERR=50) psep
  READ(iread,*,ERR=50) tmax,tout,nmax,nout
  READ(iread,*,ERR=50) gamma
  READ(iread,*,ERR=50) iener,polyk
  READ(iread,*,ERR=50) icty,ndirect,maxdensits
  READ(iread,*,ERR=50) iprterm
  READ(iread,*,ERR=50) iav,alphamin,alphaumin,alphaBmin,beta
  READ(iread,*,ERR=50) iavlim(:),avdecayconst
  READ(iread,*,ERR=50) ikernav
  READ(iread,*,ERR=50) ihvar,hfact
  READ(iread,*,ERR=50) idumpghost
  READ(iread,*,ERR=50) imhd,imagforce
  READ(iread,*,ERR=50) idivBzero,psidecayfact
  READ(iread,*,ERR=50) ianticlump,eps,neps
  READ(iread,*,ERR=50) ixsph,xsphfac
  READ(iread,*,ERR=50) igravity
  READ(iread,*,ERR=50) damp
  READ(iread,*,ERR=50) ikernel
 CLOSE(UNIT=iread)

 GOTO 55
50 WRITE(iprint,*) 'Error reading infile: re-writing with current options'
   ians = 'y'
   GOTO 1001
55 CONTINUE
!
!--check options for possible errors
!      
 IF (psep.LT.1.e-5) WRITE(iprint,100) 'psep < 1.e-5'
 IF (tout.GT.tmax) WRITE(iprint,100) 'no output tout > tmax'
 IF (nout.GT.nmax) WRITE(iprint,100) 'no output nout > nmax'
 IF (nout.EQ.0) STOP 'error in input: nout = 0'
 IF (gamma.LT.1.) WRITE(iprint,100) 'gamma < 1.0 '
 IF (abs(gamma-1.).lt.1.e-3 .AND. iener.NE.0) STOP 'must use iener = 0 for isothermal eos'
 IF ((iener.GT.0).AND.(alphaumin.LT.0.).OR.(alphaBmin.LT.0.)) THEN
    WRITE(iprint,100) 'alphaumin or alphaBmin < 0.'
 ELSEIF ((iener.EQ.0).AND.(polyk.LT.0.)) THEN
    WRITE(iprint,100) 'polyk < 0.'      
 ENDIF
 IF ((iav.NE.0).AND.(alphamin.LT.0. .OR. beta.LT.0.) ) THEN
    WRITE(iprint,100) 'AV alpha or beta < 0.'      
 ENDIF
 IF ((iavlim(1).GT.0).AND.(alphamin.GE.1.)) THEN
    WRITE(iprint,100) 'using AV limiter, but alphamin set > 1.0'
 ENDIF
 IF (ANY(iavlim.GT.0).AND.((avdecayconst.LE.0.01).OR.(avdecayconst.GT.0.5))) THEN
    WRITE(iprint,100) 'AV decay constant not in range 0.01-0.5'
 ENDIF     
 IF ((ikernav.LE.0).OR.(ikernav.GT.3)) THEN
    WRITE(iprint,100) 'kernel averaging not set (ikernav)'
 ENDIF
 IF ((hfact.LE.1.0).OR.(hfact.GT.2.0)) THEN
    WRITE(iprint,100) 'hfact too low/high (1.0 < hfact < 2.0)'
 ENDIF
 IF (psidecayfact.LT.0.0) THEN
    WRITE(iprint,100) 'psidecayfact < 0.0'
 ENDIF
 
100   FORMAT(/' read_infile: warning: ',a)
 RETURN
            
999   WRITE(iprint,1000) infile      
1000  FORMAT (' Input file ',a20,' not found')
      WRITE(*,*) ' Would you like to create one with default options?'
      READ*,ians

1001  CONTINUE
      IF (ians.EQ.'y'.OR.ians.EQ.'Y') CALL write_infile(infile)

      STOP 'exiting...'
      
END SUBROUTINE read_infile
