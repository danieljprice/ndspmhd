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
 
 USE infiles
 IMPLICIT NONE
 INTEGER :: i,j,nruns,iener_default
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
 PRINT*,' initial psep = ',psep
 iener_default = iener
 
 psidecayfact = -0.1
 hfact = 1.0
 j = 0
 
 DO i=1,nruns
    IF (i.GE.10) THEN
       filenum = ACHAR(48+i/10)//ACHAR(48+mod(i,10))
       infile = TRIM(filename)//filenum//'.in'
    ELSE
       filenum(1:1) = ACHAR(48+mod(i,10))
       infile = TRIM(filename)//filenum(1:1)//'.in'
    ENDIF
    
    j = j + 1
    IF (j.GT.9) THEN
       j = 0
       psep = 0.5*psep
       psidecayfact = -0.1
    ENDIF
    psidecayfact = psidecayfact + 0.1
    
    PRINT*,'Writing input file ',infile, ' psep = ',psep,' sigma = ',psidecayfact
    
    CALL write_infile(infile)

 ENDDO
 
END PROGRAM multirun
