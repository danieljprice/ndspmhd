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
 USE setup_params
 USE timestep
 USE xsph
 USE anticlumping
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
!--read the generic input file
! 
 PRINT*,' reading multirun.in... '
 CALL read_infile('multirun.in')
 PRINT*,' initial psep = ',psep
 
 DO i=1,nruns
    SELECT CASE(i)
       CASE(1)
        iener = 2
       CASE(2)
        imhd = 1
       CASE(3)
        iener = 3
        imhd = 11
       CASE(4)
        iener = 3
	imhd = 1
       CASE(5)
        ikernav = 1
       CASE(6)
        ikernav = 2
       CASE(7)
        ikernav = 3
       CASE(8)
        iav = 1
       CASE(9)
        iav = 3
       CASE(10)
        iav = 2
       CASE(11)
        idivBzero = 2
	psidecayfact = 0.1
       CASE(12)
        idivBzero = 2
	psidecayfact = 0.2
       CASE(13)
        idivBzero = 2
	psidecayfact = 0.4
       CASE(14)
        idivBzero = 2
	psidecayfact = 0.8
       CASE(15)
        idivBzero = 2
	psidecayfact = 1.6
    END SELECT
    IF (i.GE.10) THEN
       filenum = ACHAR(48+i/10)//ACHAR(48+mod(i,10))
       infile = TRIM(filename)//filenum//'.in'
    ELSE
       filenum(1:1) = ACHAR(48+mod(i,10))
       infile = TRIM(filename)//filenum(1:1)//'.in'
    ENDIF
    PRINT*,' writing input file ',infile
    CALL write_infile(infile)
 ENDDO
 
END PROGRAM multirun
