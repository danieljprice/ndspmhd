!!-----------------------------------------------------------------------
!! computes multiple infiles where one or several parameters are varied
!!
!! this one varies hfact from 1.0 to 2.0
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
 CHARACTER :: filename*15,infile*20,rfile*20,filenum*2,charnruns*3
 REAL :: rstress
 
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
!-- or read the generic input file
! 
 PRINT*,' reading multirun.in... '
 CALL read_infile('multirun.in')
  
 rstress = 1.0
 
 DO i=1,nruns
    IF (i.GE.10) THEN
       filenum = ACHAR(48+i/10)//ACHAR(48+mod(i,10))
       infile = TRIM(filename)//filenum//'.in'
       rfile = TRIM(filename)//filenum//'.rstress'
    ELSE
       filenum(1:1) = ACHAR(48+mod(i,10))
       infile = TRIM(filename)//filenum(1:1)//'.in'
       rfile = TRIM(filename)//filenum(1:1)//'.rstress'
    ENDIF
    PRINT*,' writing input file ',infile
    PRINT*,' rstress = ',rstress
    
    OPEN(unit=20,file=rfile,status='replace',form='formatted')
     WRITE(20,*) rstress
    CLOSE(unit=20)
    
    CALL write_infile(infile)
    
    rstress = rstress - 1.0

 ENDDO
 
END PROGRAM multirun
