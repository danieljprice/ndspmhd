!!-----------------------------------------------------------------------
!! writes multiple infiles where one or several parameters are varied
!! this version changes equation types etc
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
 INTEGER :: i,j,nruns
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
 
 IF (filename.EQ.'' .OR. nruns.LE.0) THEN
  PRINT*,'Usage: multirun filename nruns'
  STOP
 ENDIF
 
 PRINT*,' filename = ',filename,' nruns = ',nruns 
!
!--set default options
!
 CALL set_default_options
!
!--setup for shock tube
!
 !psep = 0.001485
 const = 1./SQRT(4.*pi)
    
       
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
!
!--read the generic input file
! 
 PRINT*,' reading multirun.in... '
 CALL read_infile('multirun.in')
 PRINT*,' initial psep = ',psep

 gamma = 1.66666666666667
 tmax = 0.5
 tout = 0.05
 
 DO i=1,nruns
    SELECT CASE(i)
       CASE(1)
        imhd = 11
        hfact = 1.5
        iav = 1
       CASE(2)
        iav = 2
       CASE(3)
        iav = 3
       CASE(4)
        hfact = 1.5
        imhd = 1
       CASE(5)
        iav = 2
       CASE(6)
        iav = 3
       CASE(7)
        imhd = 11
        hfact = 1.2	
	iav = 1
       CASE(8)
	iav = 2
       CASE(9)
	iav = 3        	
       CASE(10)
        hfact = 1.5
        imhd = 11
	iav = 2
        idivBzero = 2
	psidecayfact = 0.0
       CASE(11)
        idivBzero = 2
	psidecayfact = 0.1
       CASE(12)
        idivBzero = 2
	psidecayfact = 0.2
       CASE(13)
        idivBzero = 2
	psidecayfact = 0.4
    END SELECT
    IF (i.GE.10) THEN
       filenum = ACHAR(48+i/10)//ACHAR(48+mod(i,10))
       infile = TRIM(filename)//filenum//'.in'
       shkfile = TRIM(filename)//filenum//'.shk'
    ELSE
       filenum(1:1) = ACHAR(48+mod(i,10))
       infile = TRIM(filename)//filenum(1:1)//'.in'
       shkfile = TRIM(filename)//filenum(1:1)//'.shk'
    ENDIF
    PRINT*,' writing input file ',infile
    CALL write_infile(infile)
    
!    PRINT*,'Writing ',shkfile,' with initial left/right states'
!    OPEN(UNIT=11,FILE=shkfile,STATUS='replace',FORM='formatted')
!       WRITE(11,*) rholeft,rhoright
!       WRITE(11,*) prleft,prright
!       WRITE(11,*) vxleft,vxright
!       WRITE(11,*) vyleft,vyright              
!       WRITE(11,*) vzleft,vzright
!       WRITE(11,*) Bxinit
!       WRITE(11,*) Byleft,Byright
!       WRITE(11,*) Bzleft,Bzright       
!    CLOSE(UNIT=11)
    
    ENDDO
 
END PROGRAM multirun
