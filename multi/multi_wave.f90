!!-----------------------------------------------------------------------
!! computes multiple infiles where one or several parameters are varied
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
 CHARACTER :: filename*15,infile*20,filenum*2,charnruns*3
 
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
 
 DO i=1,nruns
    IF (i.GE.10) THEN
       filenum = ACHAR(48+i/10)//ACHAR(48+mod(i,10))
       infile = filename(1:LEN_TRIM(filename))//filenum//'.in'
    ELSE
       filenum(1:1) = ACHAR(48+mod(i,10))
       infile = filename(1:LEN_TRIM(filename))//filenum(1:1)//'.in'
    ENDIF
    PRINT*,' writing input file ',infile, ' psep = ',psep
    CALL write_infile(infile)
    psep = psep/2.		! ie. decrease resolution
 ENDDO
 
END PROGRAM multirun
