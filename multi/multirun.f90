!!-----------------------------------------------------------------------
!! computes multiple infiles where one or several parameters are varied
!!
!! this varies h and the div B cleaning parameters
!!-----------------------------------------------------------------------
PROGRAM multirun
 USE dimen_mhd
 USE loguns
 USE artvi
 USE eos
 USE options
 USE setup_params
 USE timestep
 USE xsph
 USE anticlumping
 IMPLICIT NONE
 INTEGER :: i,j,nruns,iener_default
 INTEGER :: ierr
 CHARACTER :: filename*15,infile*20,filenum*2,charnruns*3,shkfile*20
 REAL :: rholeft,rhoright,prleft,prright
 REAL :: vxleft,vxright,vyleft,vyright,vzleft,vzright
 REAL :: Bxinit,Byleft,Byright,Bzleft,Bzright,const
 
 nruns = 0
 iread = 11
 CALL getarg(1,filename)
 CALL getarg(2,charnruns)
 
 read(charnruns,*,iostat=ierr) nruns
 
 IF (filename.EQ.'' .OR. ierr.NE.0) THEN
  PRINT*,'Usage: multirun filename nruns'
  STOP
 ENDIF
 
 PRINT*,' filename = ',filename,' nruns = ',nruns 
  
!
!--set default options
!
 CALL set_default_options
!
!-- or read the generic input file
! 
 PRINT*,' reading multirun.in... '
 CALL read_infile('multirun.in')
 PRINT*,' initial psep = ',psep
 iener_default = iener
 
 psidecayfact = 1.0
 !hfact = 1.0
 j = 0

 const = 1./SQRT(4.*pi)

       psep = 0.0025
       gamma = 1.66666666666667
       tmax = 0.08
       tout = 0.04
       
       rholeft = 1.0
       rhoright = 1.0
       prleft = 20.0
       prright = 1.0
       vxleft = 10.0
       vxright = -10.0
       vyleft = 0.
       vyright = 0.
       vzleft = 0.
       vzright = 0.
       Bxinit = 5.*const
       Byleft = 5.*const
       Byright = 5.*const
       Bzleft = 0.
       Bzright = 0.
        
 DO i=1,nruns
    IF (i.GE.10) THEN
       filenum = ACHAR(48+i/10)//ACHAR(48+mod(i,10))
       infile = TRIM(filename)//filenum//'.in'
       shkfile = TRIM(filename)//filenum//'.shk'
    ELSE
       filenum(1:1) = ACHAR(48+mod(i,10))
       infile = TRIM(filename)//filenum(1:1)//'.in'
       shkfile = TRIM(filename)//filenum(1:1)//'.shk'
    ENDIF
    
    j = j + 1
    !IF (j.GT.9) THEN
    !   j = 0
    !   psep = 0.5*psep
    !   psidecayfact = -0.1
    !ENDIF
    !psidecayfact = psidecayfact + 0.1
    !hfact = hfact + 0.2

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
        
    PRINT*,'Writing input file ',infile, ' psep = ',psep,' sigma = ',psidecayfact,' hfact = ',hfact
    
    CALL write_infile(infile)

 ENDDO
 
END PROGRAM multirun
