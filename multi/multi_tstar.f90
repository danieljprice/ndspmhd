!!-----------------------------------------------------------------------
!! computes multiple infiles where one or several parameters are varied
!!
!! this one does toy star oscillations
!!-----------------------------------------------------------------------
PROGRAM multirun
 USE dimen_mhd
 USE loguns
 USE artvi
 USE eos
 USE options
 USE polyconst
 USE setup_params
 USE timestep
 USE xsph
 USE anticlumping
 IMPLICIT NONE
 INTEGER :: i,j,nruns,iener_default,norder
 INTEGER :: int_from_string
 CHARACTER :: filename*15,infile*20,tstarfile*20,filenum*2,charnruns*3
 REAL :: H,C,A,sigma,omega,period
 
 nruns = 0
 iread = 11
 CALL getarg(1,filename)
 CALL getarg(2,charnruns)
 
 nruns = int_from_string(charnruns)
 
 IF (filename.EQ.'') THEN
  PRINT*,' Enter filename for runs'
  READ*,filename
 ENDIF
 
 IF (nruns.LE.0) THEN
  PRINT*,' Enter number of runs:'
  READ*,nruns
 ENDIF
 
 PRINT*,' filename = ',filename,' nruns = ',nruns 
  
! PRINT*,' Enter particle separation run 1 '
! READ*,psep
!
!--set all options (either set them all here or read from an infile)
! 
 tmax = 10.0
 tout = 1.0
 nmax = 1000000
 nout = -1
 gamma = 1.6666666666666666666666667
 iener = 2
 iener_default = 2
 gconst = 0.0
 polyk = 0.6
 icty = 0
 ndirect = nmax
 iprterm = 4
 iav = 2
 alphamin = 0.1
 beta = 2.0
 iavlim = 1
 avconst = 0.1
 ikernav = 3
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
!
!-- or read the generic input file
! 
 PRINT*,' reading multirun.in... '
 CALL read_infile('multirun.in')
 PRINT*,' initial psep = ',psep
 iener_default = iener
  
 DO i=0,nruns
    IF (i.GE.10) THEN
       filenum = ACHAR(48+i/10)//ACHAR(48+mod(i,10))
       infile = filename(1:LEN_TRIM(filename))//filenum//'.in'
       tstarfile = filename(1:LEN_TRIM(filename))//filenum//'.tstar'
    ELSE
       filenum(1:1) = ACHAR(48+mod(i,10))
       infile = filename(1:LEN_TRIM(filename))//filenum(1:1)//'.in'
       tstarfile = filename(1:LEN_TRIM(filename))//filenum(1:1)//'.tstar'
    ENDIF
    PRINT*,' writing input file ',infile, ' psep = ',psep
    
    H = 1.0
    C = 1.0
    A = 0.05/SQRT(2.)
    sigma = 1./SQRT(2.)
    norder = i
!
!--run for 10 periods with output every period
!    
    omega = SQRT(0.5*(norder+1)*(norder+2))
    period = 2.*3.1415926536/omega
    tout = 0.25*period
    tmax = 10.*period
    
    PRINT*,'Writing ',tstarfile,' with initial left/right states'
    OPEN(UNIT=11,FILE=tstarfile,STATUS='replace',FORM='formatted')
       WRITE(11,*) H,C,A
       WRITE(11,*) sigma
       WRITE(11,*) norder
    CLOSE(UNIT=11)
    
    CALL write_infile(infile)

 ENDDO
 
END PROGRAM multirun
