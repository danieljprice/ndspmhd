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
  READ(iread,*) psep
  READ(iread,*) tmax,tout,nmax,nout
  READ(iread,*) gamma
  READ(iread,*) iener,gconst,polyk
  READ(iread,*) icty,ndirect
  READ(iread,*) ialtform
  READ(iread,*) iav,alphamin,beta
  READ(iread,*) iavlim,avconst
  READ(iread,*) ikernav
  READ(iread,*) ihvar,hfact
  READ(iread,*) idumpghost
  READ(iread,*) imhd,imagforce
  READ(iread,*) idivBzero
  READ(iread,*) ianticlump,eps,neps
  READ(iread,*) ixsph,xsphfac
  READ(iread,*) igravity
  READ(iread,*) damp
  READ(iread,*) ikernel
 CLOSE(UNIT=iread)
!
!--check options for possible errors
!      
 IF (psep.LT.1.e-5) WRITE(iprint,100) 'psep < 1.e-5'
 IF (tout.GT.tmax) WRITE(iprint,100) 'no output tout > tmax'
 IF (nout.GT.nmax) WRITE(iprint,100) 'no output nout > nmax'
 IF (nout.EQ.0) STOP 'error in input: nout = 0'
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
 IF ((iavlim.GT.0).AND.((avconst.LE.0.01).OR.(avconst.GT.0.5))) THEN
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
      IF (ians.EQ.'y'.OR.ians.EQ.'Y') CALL write_infile(infile)

      STOP 'exiting...'
      
END SUBROUTINE read_infile
