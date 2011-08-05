!---------------------------------------------------------------------------
! module containing various utility subroutines 
! related to reading from ascii files
!
! written by Daniel Price, University of Exeter, dprice@astro.ex.ac.uk
! 24th April '07
!
! this is a standalone module with no dependencies
!---------------------------------------------------------------------------
module asciiutils
 implicit none
 public :: read_asciifile,get_ncolumns,ncolumnsline,safename
 
 private

contains

!---------------------------------------------------------------------------
! Generic subroutine to read all lines of an ascii file
! returns array of character strings (one per line)
! up to a maximum corresponding to the size of the array
!---------------------------------------------------------------------------
subroutine read_asciifile(filename,nlinesread,charline)
 implicit none
 character(len=*), intent(in) :: filename
 integer, intent(out) :: nlinesread
 character(len=*), dimension(:), intent(out) :: charline
 integer, parameter :: iunit = 66 ! logical unit number for read operation
 integer :: ierr,i,maxlines
 logical :: iexist
 
 nlinesread = 0
 !--if file does not exist, do nothing and return
 inquire(file=filename,exist=iexist)
 if (.not.iexist) return
 
 open(unit=iunit,file=filename,status='old',form='formatted',iostat=ierr)
 !--error opening file (but file does exist)
 if (ierr /= 0) then
    print "(a)",' ERROR opening '//trim(filename)
    return
 endif
 
 maxlines = size(charline)
 do i=1,maxlines
    read(iunit,"(a)",err=66,end=99) charline(i)
 enddo
 !--end of array limits
 print "(a)",' WARNING: array limits reached reading '//trim(filename)//' max = ',maxlines
 nlinesread = maxlines
 close(unit=iunit)
 return

 !--error encountered
66 continue
  print "(a,i6)",' ERROR reading '//trim(filename)//' at line ',i-1
  nlinesread = i-1
  close(unit=iunit)
  return

 !--reached end of file (the expected behaviour)
99 continue
  nlinesread = i-1
  close(unit=iunit)
  return

end subroutine read_asciifile

!---------------------------------------------------------------------------
! utility to work out number of columns of real numbers
! in an ascii file
!
! file must already be open and at the start
! slightly ad-hoc but its the best way I could think of!
!---------------------------------------------------------------------------
subroutine get_ncolumns(lunit,ncolumns,nheaderlines)
 implicit none
 integer, intent(in) :: lunit
 integer, intent(out) :: ncolumns,nheaderlines
 integer :: ierr,ncolprev,ncolsthisline
 character(len=2000) :: line
 logical :: nansinfile,infsinfile

 nheaderlines = 0
 line = ' '
 ierr = 0
 ncolumns = 0
 ncolprev = 666
 ncolsthisline = 0
 nansinfile = .false.
 infsinfile = .false.
!
!--loop until we find two consecutive lines with the same number of columns (but non zero)
!
 do while ((len_trim(line).eq.0 .or. ncolsthisline.ne.ncolprev .or. ncolumns.le.0) .and. ierr.eq.0)
    ncolprev = ncolumns
    read(lunit,"(a)",iostat=ierr) line
    if (index(line,'NaN').gt.0) nansinfile = .true.
    if (index(line,'Inf').gt.0) infsinfile = .true.
    if (ierr.eq.0) ncolsthisline = ncolumnsline(line)
    if (ncolsthisline.ne.0) nheaderlines = nheaderlines + 1
    ncolumns = ncolsthisline
 enddo
 !--subtract 2 from the header line count (the last two lines which were the same)
 nheaderlines = max(nheaderlines - 2,0)
 if (ierr .gt.0 .or. ncolumns.le.0) then
    ncolumns = 0
 elseif (ierr .lt. 0) then
    print*,ncolumns,ncolprev
 endif
 if (nansinfile) print "(a)",' INDIAN BREAD WARNING!! NaNs in file!!'
 if (infsinfile) print "(a)",' WARNING!! Infs in file!!'
 rewind(lunit)

 if (ncolumns.eq.0) then
    print "(a)",' ERROR: no columns of real numbers found'
 else
    print "(a,i3)",' number of data columns = ',ncolumns
 endif
 
end subroutine get_ncolumns

!---------------------------------------------------------------------------
!
! function returning the number of columns of real numbers from a given line
!
!---------------------------------------------------------------------------
integer function ncolumnsline(line)
 implicit none
 character(len=*), intent(in) :: line
 real :: dummyreal(100)
 integer :: ierr,i

 dummyreal = -666.0
 
 ierr = 0
 read(line,*,iostat=ierr) (dummyreal(i),i=1,size(dummyreal))
 if (ierr .gt. 0) then
    ncolumnsline = -1
    return
 endif

 i = 1
 ncolumnsline = 0
 do while(abs(dummyreal(i)+666.).gt.1.e-10)
    ncolumnsline = ncolumnsline + 1
    i = i + 1
    if (i.gt.size(dummyreal)) then
       print "(a)",'*** ERROR: too many columns in file'
       ncolumnsline = size(dummyreal)
       return
    endif
 enddo

end function ncolumnsline

!---------------------------------------------------------------------------
!
! function stripping '/', '\' and spaces out of filenames
!
!---------------------------------------------------------------------------
function safename(string)
 implicit none
 character(len=*), intent(in) :: string
 character(len=len(string)) :: safename
 integer :: ipos
 
 safename = string

 !--remove forward slashes which can be mistaken for directories: replace with '_'
 ipos = index(safename,'/')
 do while (ipos.ne.0)
    safename = safename(1:ipos-1)//'_'//safename(ipos+1:len_trim(safename))
    ipos = index(safename,'/')
 enddo
 
 !--remove spaces: replace with '_'
 ipos = index(trim(safename),' ')
 do while (ipos.ne.0)
    safename = safename(1:ipos-1)//'_'//safename(ipos+1:len_trim(safename))
    ipos = index(trim(safename),' ')
 enddo
 
 !--remove escape sequences: remove '\' and position following
 ipos = index(trim(safename),'\')
 do while (ipos.ne.0)
    safename = safename(1:ipos-1)//safename(ipos+2:len_trim(safename))
    ipos = index(trim(safename),'\')
 enddo

end function safename

end module asciiutils
