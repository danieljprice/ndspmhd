!---------------------------------------------------------
! Program takes n data files with concurrent timestep data
! and makes new data files with data at simultaneous times
!
! modified 10/10/03 so it doesn't hog so much memory
! also re-ordered dat array for optimisation
! now does infinite number of steps (so long as the character
! string is long enough!)
!---------------------------------------------------------
PROGRAM combinedat
 IMPLICIT NONE
 INTEGER :: i,j,k
 INTEGER, PARAMETER :: ncolmax=25
 INTEGER, PARAMETER :: nfilesmax=200
 INTEGER, PARAMETER :: npartmax = 2000
 INTEGER :: ifile,nfiles,infile,ioutfile
 INTEGER :: nsteps,ncol
 INTEGER :: int_from_string
 CHARACTER(len=24), DIMENSION(nfilesmax) :: filename
 CHARACTER(len=20) :: rootname
 CHARACTER(len=2) :: fnum
 CHARACTER(len=6) :: charnsteps,charnfiles
 REAL, DIMENSION(nfilesmax) :: time,gamma,hfact
 INTEGER, DIMENSION(nfilesmax) :: npart,nprint,ndim,ndimV,ndata
 REAL, DIMENSION(ncolmax,npartmax,nfilesmax) :: dat
 
 nfiles = 0
 CALL getarg(1,rootname)	! get rootname off command line
 CALL getarg(2,charnfiles)	! get nfiles off command line
 CALL getarg(3,charnsteps)	! get nsteps off command line
 
 nfiles = int_from_string(charnfiles)
 nsteps = int_from_string(charnsteps)
 
 PRINT*,' runname = ',TRIM(rootname)
 PRINT*,' nfiles = ',nfiles
 PRINT*,' nsteps = ',nsteps
  
 IF ((nfiles.LE.0).OR.(nfiles.GT.nfilesmax)) THEN
    PRINT*,'Enter number of files to read:'
    READ*,nfiles
 ENDIF
 
 IF (nfiles.GT.nfilesmax) THEN
    nfiles = nfilesmax
    PRINT*,'nfiles > nfilesmax, nfiles = ',nfiles
 ENDIF
 
 IF ((nsteps.LE.0).OR.(nsteps.GT.1000)) THEN
    PRINT*,'Enter number of steps to read:'
    READ*,nsteps
 ENDIF   
!
!--open all data files
! 
 DO i=1,nfiles
    IF (rootname(1:1).EQ.' ') THEN
       PRINT*,' Enter filename number ',i
       READ*,filename(i)
    ELSE
       IF (i.LT.10) THEN
       filename(i) = TRIM(rootname)//ACHAR(48+MOD(i,10))//'.dat'
       ELSE
       filename(i) = TRIM(rootname)//ACHAR(48+i/10)//ACHAR(48+MOD(i,10))//'.dat'
       ENDIF
    ENDIF
    PRINT*,' opening ',filename(i)
    infile = 10+i  ! logical unit number for input files (*must* be same as below)
    OPEN(UNIT=infile,ERR=101,FILE=filename(i),STATUS='old',FORM='formatted')
    GOTO 102
101 PRINT*,'Error opening ',filename(i)
102 CONTINUE    
 ENDDO
!
!--read each time step
! 
 DO j=1,nsteps
!
!--open new data file for this time step
!
    ioutfile = 1	! logical unit number for output file (re-use each time)
    fnum = ACHAR(48+j/10)//ACHAR(48+mod(j,10))
    PRINT*,'---------------------------------------------'
    PRINT*,' Creating data file ',TRIM(filename(1))//fnum
!    READ*
    OPEN(UNIT=ioutfile,FILE=TRIM(filename(1))//fnum,	&
         STATUS='replace',FORM='formatted')
    
    DO ifile=1,nfiles
       infile = 10+ifile ! logical unit number for input files
       PRINT*,'reading ',filename(ifile)
       READ(infile,*,END=221) time(ifile),npart(ifile),nprint(ifile), &
                              gamma(ifile),hfact(ifile),ndim(ifile),  &
                              ndimV(ifile),ndata(ifile)
       ncol = ndata(ifile)			      			      
       READ(infile,*,END=222) (dat(1:ncol,k,ifile), k=1,nprint(ifile))
       PRINT*,' t = ',time(ifile)
!       PRINT*,' first line = ',dat(ifile,1:ncol,1)

      GOTO 223
221   CONTINUE	! timestep not there at all
         PRINT*,' *** timestep missing... skipping '
	 GOTO 223
222   CONTINUE	! catch the error if data incomplete
         PRINT*,' *** timestep incomplete, want ',nprint(ifile),' found ',k-1
         nprint(ifile) = k-1
223   CONTINUE     
    ENDDO
    
    PRINT*,' Writing to data file... '
    DO ifile=1,nfiles
       PRINT*,' t = ',time(ifile),npart(ifile),nprint(ifile)
       WRITE(ioutfile,*) time(ifile),npart(ifile),nprint(ifile), &
                         gamma(ifile),hfact(ifile),ndim(ifile),  &
			 ndimV(ifile),ndata(ifile)
       DO k=1,nprint(ifile)
          WRITE(ioutfile,10) dat(1:ncol,k,ifile)
       ENDDO  
10     FORMAT (16(1pe14.6,1x))	! make sure the format statement has >/=	
					! max number of columns in the write statement
    ENDDO
    CLOSE(UNIT=ioutfile)
    
 ENDDO	! over timesteps
 
 DO i=1,nfiles
    ifile = 10+i
    CLOSE(UNIT=ifile)
 ENDDO
      
END PROGRAM combinedat

!---------------------------------------------------------------------
!  function to convert a character string representing a number
!  to an integer value
!
!  Daniel Price, Institute of Astronomy, Cambridge, Oct 2003
!---------------------------------------------------------------------

FUNCTION int_from_string(string)
 IMPLICIT NONE
 INTEGER :: int_from_string,idigit,i,ipower,maxdigits
 CHARACTER(LEN=*) :: string
 
 ipower = -1
 maxdigits = LEN(string)
 int_from_string = 0

 DO i=maxdigits,1,-1	! down through the characters
    IF (string(i:i).NE.' ') THEN
       ipower = ipower + 1
       idigit = IACHAR(string(i:i)) - 48
       IF ((idigit.LT.0).OR.(idigit.GT.10)) THEN
          PRINT*,'Error: non-numeric character in input string'
	  int_from_string = 0
	  RETURN
       ELSE
          int_from_string = int_from_string + idigit*10**ipower
       ENDIF
    ENDIF
 ENDDO   
 
END FUNCTION int_from_string
