!!-----------------------------------------------------------------------
!! computes multiple infiles where one or several parameters are varied
!!
!! this one does all 7 MHD shock tube tests
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
 
 USE infiles
 IMPLICIT NONE
 INTEGER :: i,j,nruns,iener_default
 INTEGER :: ierr
 CHARACTER :: filename*15,infile*20,shkfile*20,filenum*2,charnruns*3
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
!--read the generic input file
! 
 PRINT*,' reading multirun.in... '
 CALL read_infile('multirun.in')
 PRINT*,' initial psep = ',psep
 iener_default = iener
 
 const = 1./SQRT(4.*pi)
 
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
    PRINT*,' writing input file ',infile, ' psep = ',psep
    
    IF (i.EQ.1) THEN	! Brio/Wu
       psep = 0.0056
       gamma = 2.0
       tmax = 0.1
       tout = 0.05

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
    ELSEIF (i.EQ.2) THEN	! Balsara 1
       psep = 0.00375
       gamma = 1.6666666666667
       tmax = 0.15
       tout = 0.075

       rholeft = 1.0
       rhoright = 0.2
       prleft = 1.0
       prright = 0.1
       vxleft = 0.
       vxright = 0.
       vyleft = 0.
       vyright = 0.
       vzleft = 0.
       vzright = 0.
       Bxinit = 1.
       Byleft = 1.0
       Byright = 0.0
       Bzleft = 0.
       Bzright = 0.     

    ELSEIF (i.EQ.3) THEN	! 7 discont.
       psep = 0.001485
       gamma = 1.6666666666667
       tmax = 0.2
       tout = 0.1

       rholeft = 1.08
       rhoright = 1.0
       prleft = 0.95
       prright = 1.0
       vxleft = 1.2
       vxright = 0.
       vyleft = 0.01
       vyright = 0.
       vzleft = 0.5
       vzright = 0.
       Bxinit = 2.*const
       Byleft = 3.6*const
       Byright = 4.0*const
       Bzleft = 2.*const
       Bzright = 2.*const         
    ELSEIF (i.EQ.4) THEN		! isothermal 6 discont.
       psep = 0.001485
       gamma = 1.0
       tmax = 0.2
       tout = 0.1
       iener = 0	! isothermal
       polyk = 1.0	! isothermal sound speed

       rholeft = 1.08
       rhoright = 1.0
       prleft = 0.95
       prright = 1.0
       vxleft = 1.2
       vxright = 0.
       vyleft = 0.01
       vyright = 0.
       vzleft = 0.5
       vzright = 0.
       Bxinit = 2.*const
       Byleft = 3.6*const
       Byright = 4.0*const
       Bzleft = 2.*const
       Bzright = 2.*const         

    ELSEIF (i.EQ.5) THEN	! rarefaction
       psep = 0.002
       gamma = 1.66666666666666667
       tmax = 0.1
       tout = 0.05
       iener = iener_default	! back from isothermal

       rholeft = 1.0
       rhoright = 1.0
       prleft = 1.0
       prright = 1.0
       vxleft = -1.0
       vxright = 1.
       vyleft = 0.
       vyright = 0.
       vzleft = 0.
       vzright = 0.
       Bxinit = 0.
       Byleft = 1.
       Byright = 1.
       Bzleft = 0.
       Bzright = 0.
        
    ELSEIF (i.EQ.6) THEN        ! two fast shocks
       psep = 0.0025
       gamma = 1.66666666666667
       tmax = 0.03
       tout = 0.015
       
       rholeft = 1.0
       rhoright = 1.0
       prleft = 1.0
       prright = 1.0
       vxleft = 36.87
       vxright = -36.87
       vyleft = -0.1546
       vyright = 0.
       vzleft = -0.03864
       vzright = 0.
       Bxinit = 4.*const
       Byleft = 4.*const
       Byright = 4.*const
       Bzleft = 1.*const
       Bzright = 1.*const   
    
    ELSEIF (i.EQ.7) THEN   ! 2D shock problem
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
    ENDIF
    
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
