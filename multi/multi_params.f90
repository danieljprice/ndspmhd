!!-----------------------------------------------------------------------
!! computes multiple infiles where one or several parameters are varied
!!
!! this one does Brio/Wu shock test, varying the stress parameters
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
 INTEGER :: i,j,nruns,iener_default
 INTEGER :: int_from_string
 CHARACTER :: filename*15,infile*20,shkfile*20,filenum*2,charnruns*3
 REAL :: rholeft,rhoright,prleft,prright
 REAL :: vxleft,vxright,vyleft,vyright,vzleft,vzright
 REAL :: Bxinit,Byleft,Byright,Bzleft,Bzright,const
 
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
 ialtform = 4
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
 
 const = 1./SQRT(4.*pi)
 eps = 0.2	! start very low
 neps = 3
 
 DO i=1,nruns
    IF (i.GE.10) THEN
       filenum = ACHAR(48+i/10)//ACHAR(48+mod(i,10))
       infile = TRIM(filename)//filenum//'.in'
       shkfile = TRIM(filename)//filenum//'.shk'
    ELSE
       filenum(1:1) = ACHAR(48+mod(i,10))
       infile = filename(1:LEN_TRIM(filename))//filenum(1:1)//'.in'
       shkfile = filename(1:LEN_TRIM(filename))//filenum(1:1)//'.shk'
    ENDIF
    PRINT*,' writing input file ',infile, ' psep = ',psep
    
!
!--brio/wu problem
!
       psep = 0.0056
       gamma = 2.0
       tmax = 0.1
       tout = 0.05
       
       eps = eps + 0.1
       IF (MOD(i,20).EQ.0) THEN
          neps = neps + 1
	  eps = 0.2
       ENDIF

       PRINT*,'neps = ',neps
       PRINT*,'eps = ',eps

       rholeft = 1.0
       rhoright = 0.125
       prleft = 1.0
       prright = 0.1
       vxleft = 0.
       vxright = 0.
       vyleft = 0.
       vyright = 0.
       vzleft = 0.
       vzright = 0.
       Bxinit = 0.75
       Byleft = 1.0
       Byright = -1.0
       Bzleft = 0.
       Bzright = 0.     
    
    PRINT*,'Writing ',shkfile,' with initial left/right states'
    OPEN(UNIT=11,FILE=shkfile,STATUS='replace',FORM='formatted')
       WRITE(11,*) rholeft,rhoright
       WRITE(11,*) prleft,prright
       WRITE(11,*) vxleft,vxright
       WRITE(11,*) vyleft,vyright              
       WRITE(11,*) vzleft,vzright
       WRITE(11,*) Bxinit
       WRITE(11,*) Byleft,Byright
       WRITE(11,*) Bzleft,Bzright       
    CLOSE(UNIT=11)
    
    CALL write_infile(infile)

 ENDDO
 
END PROGRAM multirun
