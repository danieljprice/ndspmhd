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
 USE setup_params
 USE timestep
 USE xsph
 USE anticlumping
 
 USE infiles
 IMPLICIT NONE
 INTEGER :: i,j,nruns
 INTEGER :: ierr
 CHARACTER :: filename*15,infile*20,filenum*2,charnruns*3
 
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
 
 hfact = 1.0
 PRINT*,' initial hfact = ',hfact
 
 DO i=1,nruns
    IF (i.GE.10) THEN
       filenum = ACHAR(48+i/10)//ACHAR(48+mod(i,10))
       infile = TRIM(filename)//filenum//'.in'
    ELSE
       filenum(1:1) = ACHAR(48+mod(i,10))
       infile = TRIM(filename)//filenum(1:1)//'.in'
    ENDIF
    PRINT*,' writing input file ',infile, ' hfact = ',hfact
    
    CALL write_infile(infile)
    
    hfact = hfact + 0.05

 ENDDO
 
END PROGRAM multirun
