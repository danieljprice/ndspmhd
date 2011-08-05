!!-----------------------------------------------------------------------
!! writes multiple infiles where one or several parameters are varied
!! this version changes equation types etc
!!-----------------------------------------------------------------------
PROGRAM multirun
 USE dimen_mhd
 USE loguns
 USE artvi
 USE options
 USE setup_params
 USE timestep
 USE xsph
 
 USE infiles
 IMPLICIT NONE
 INTEGER :: i,j,nruns,nres
 INTEGER :: ierr
 CHARACTER :: filename*15,infile*20,filenum*2,charnruns*3,res*4
 
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
 
 iavlim(1) = 2
 alphamin = 0.
 alphaumin = 0.
 psep = 0.04
 nres = 3
 DO i=1,nres
    psep = psep*0.5
    select case(i)
    case(1)
      res = 'lres'
    case(2)
      res = 'sres'
    case(3)
      res = 'hres'    
    end select

    DO j = 1,4
       select case(j)
       case(1)
         iener = 2
         iprterm = 0
         iavlim(2) = 0
         infile = 'he'//res//'.in' 
       case(2)
         iener = 2
         iprterm = 0
         iavlim(2) = 1
         infile = 'hecond'//res//'.in' 
       case(3)
         iener = 10
         iprterm = 10
         iavlim(2) = 0
         infile = 'hert'//res//'.in' 
       case(4)
         iener = 10
         iprterm = 10
         iavlim(2) = 1
         infile = 'hertcond'//res//'.in' 
       end select    
       PRINT*,' writing input file ',infile
       CALL write_infile(infile)
    ENDDO
 ENDDO
 
END PROGRAM multirun
