!------------------------------------------------------------------------
!
! module containing subroutines for applying transformations to data
! prior to plotting.
!
! the transformations are:
!
!  1) log10(x)
!  2) |x|
!  3) 1/x
!  4) sqrt(x)     
!  5) x^2
!
! * combinations of transformations are done when the 
!   input number is > 10 (e.g. 321 means 1/x, then abs, then log10)
!
! subroutines contained within this module are the following:
!
!   transform        : applies transformation to a one dimensional array
!   transform2       : applies the transformation to a two dimensional array
!   transform_limits : transforms the plot limits appropriately
!   transform_label  : changes the plot label appropriately
!
! Written for use in supersphplot by
! Daniel Price, Institute of Astronomy, Cambridge, UK 2002-2004
!               University of Exeter, UK              2004-
! dprice@astro.ex.ac.uk
!
!------------------------------------------------------------------------
module transforms
 integer, parameter, public :: ntrans = 5  ! this is the number of different transformations

contains
!------------------------------------------------------------------------
!
!  subroutine returns log, 1/x of a given array
!
!  * can specify up to 9 individual operations to perform
!  * combinations of transformations are done when the 
!    input number is > 10 (e.g. 321 means 1/x, then abs, then log10)
!
!------------------------------------------------------------------------
subroutine transform(array,itrans)
  implicit none
  integer, intent(in) :: itrans
  real, dimension(:), intent(inout) :: array
  real, dimension(size(array)) :: arraytemp
  character(len=20) :: string
  integer :: i
  !
  !--extract the digits from the input number 
  !
  if (itrans.gt.0) then      

     write(string,*) itrans
     !
     !--do a transformation for each digit
     !
     arraytemp = array

     do i=1,len_trim(string)
        !
        !--perform transformation appropriate to this digit     
        !
        select case(string(i:i))
        case('1')
           where (arraytemp > 0)
              arraytemp = log10(arraytemp)
           elsewhere
              arraytemp = 0.
           end where
        case('2')
           arraytemp = abs(arraytemp)    
        case('3')
           where (arraytemp .ne. 0)
              arraytemp = 1./arraytemp
           elsewhere
              arraytemp = 0.
           end where
        case('4') 
           where (arraytemp .gt. 0)
              arraytemp = sqrt(arraytemp)
           elsewhere
              arraytemp = 0.
           end where
        case('5')
           arraytemp = arraytemp**2             
        end select
     enddo

     array = arraytemp
  endif

end subroutine transform

!------------------------------------------------------------------------
!
!  same as transform but for a two dimensional array
!  applies the transformation to the same array as was input
!
!------------------------------------------------------------------------
subroutine transform2(array,itrans,isizex,isizey)
  implicit none
  integer, intent(in) :: itrans,isizex,isizey
  real, dimension(isizex,isizey), intent(inout) :: array
  real, dimension(isizex,isizey) :: arraytemp
  character(len=20) :: string
  integer :: i
  !
  !--extract the digits from the input number 
  !
  if (itrans.gt.0) then      

     write(string,*) itrans
     !
     !--do a transformation for each digit     
     !
     arraytemp = array

     do i=1,len_trim(string)
        !
        !--perform transformation appropriate to this digit     
        !
        select case(string(i:i))
        case('1')
           where (arraytemp > 0)
              arraytemp = log10(arraytemp)
           elsewhere
              arraytemp = 0.
           end where
        case('2')
           arraytemp = abs(arraytemp)    
        case('3')
           where (arraytemp .ne. 0)
              arraytemp = 1./arraytemp
           elsewhere
              arraytemp = 0.
           end where
        case('4') 
           where (arraytemp .gt. 0)
              arraytemp = sqrt(arraytemp)
           elsewhere
              arraytemp = 0.
           end where
        case('5')
           arraytemp = arraytemp**2         
        end select
     enddo

     array = arraytemp
 
  endif

end subroutine transform2

!------------------------------------------------------------------------
!
!  same as transform but for the plot limits
!  (min can become max and vice versa)
!
!------------------------------------------------------------------------
subroutine transform_limits(xmin,xmax,itrans)
  implicit none
  integer, intent(in) :: itrans
  real, intent(inout) :: xmin,xmax
  real :: xmintemp,xmaxtemp
  character(len=20) :: string
  integer :: i
  !
  !--extract the digits from the input number 
  !
  if (itrans.gt.0) then      

     write(string,*) itrans
     !
     !--do a transformation for each digit     
     !
     xmintemp = xmin
     xmaxtemp = xmax

     do i=1,len_trim(string)
        !
        !--perform transformation appropriate to this digit     
        !
        select case(string(i:i))
        case('1')
           if (xmintemp > 0) then
              xmintemp = log10(xmintemp)
           elseif (xmintemp.eq.0) then
              print*,' log10(xmin = 0): min set to 10-12'
              xmintemp = -12.
           endif
           if (xmaxtemp > 0) then
              xmaxtemp = log10(xmaxtemp)
           elseif (xmaxtemp.eq.0) then
              print*,' log10(xmax = 0): max set to 10-12'
              xmaxtemp = -12.
           endif
        case('2')
           if ((xmintemp.lt.0. .and. xmaxtemp.gt.0.) &
           .or.(xmaxtemp.lt.0. .and. xmintemp.gt.0.)) then
           !
           !--minimum is zero if limits have opposite signs
           !
              xmaxtemp = max(abs(xmintemp),abs(xmaxtemp))
              xmintemp = 0.
           else
           !
           !--or just take magnitude
           !
              xmintemp = abs(xmintemp)
              xmaxtemp = abs(xmaxtemp)
           endif
        case('3')
           if (xmintemp .ne. 0) then
              xmintemp = 1./xmintemp
           else
              xmintemp = 0.
           endif
           if (xmaxtemp .ne. 0) then
              xmaxtemp = 1./xmaxtemp
           else
              xmaxtemp = 0.             
           endif
        case('4') 
           if (xmintemp .ge. 0) then
              xmintemp = sqrt(xmintemp)
           else
              xmintemp = 0.
           endif
           if (xmaxtemp .ge. 0) then
              xmaxtemp = sqrt(xmaxtemp)
           else
              xmaxtemp = 0.             
           endif
        case('5')
           xmintemp = xmintemp**2
           xmaxtemp = xmaxtemp**2
        end select
     enddo

     xmin = min(xmintemp,xmaxtemp)
     xmax = max(xmintemp,xmaxtemp)

  endif

end subroutine transform_limits

!------------------------------------------------------------------------
!
!  function to adjust the label of a plot if log, 1/x etc
!
!  Note: *cannot* put print or write statements into this function
!        as it is used in the middle of write or print statements
!        this means that finding the digits is a bit trickier
!      
!------------------------------------------------------------------------
function transform_label(label,itrans)
  implicit none
  integer :: itrans,itransmulti,i,ndigits
  integer, dimension(20) :: digit
  character(len=*) :: label
  character(len=len(label)+20) :: transform_label
  character(len=len(label)+20) :: temp_label      
  !
  !--extract the digits from the input number 
  !            
  if (itrans.gt.0) then      
     call get_digits(itrans,digit,ndigits)         
     temp_label = label      
     !
     !--do a transformation for each digit     
     !
     do i=1,ndigits
        itransmulti = digit(i)
        !
        !--perform transformation appropriate to this digit     
        !
        select case(itransmulti)
        case(1)
           temp_label = 'log\d10\u'//trim(temp_label)
        case(2)
           temp_label = '|'//trim(temp_label)//'|'
        case(3)
           temp_label = '1/'//trim(temp_label)
        case(4)
           temp_label = 'SQRT('//trim(temp_label)//')'
        case(5)
           temp_label = trim(temp_label)//'\u2\d'
        case default
           temp_label = trim(temp_label)
        end select
     enddo

     transform_label = temp_label
  else
     transform_label = label     
  endif

end function transform_label

!------------------------------------------------------------------------
!     get_digits: for an integer i returns number of digits it contains
!     and a list of these *without* using write statements
!
!     i            : integer to split into digits
!     nmax           : dimensions of digits array
!     digits(nmax) : array of digits
!     ndigits      : number of digits in i
!------------------------------------------------------------------------

subroutine get_digits(i,digits,ndigits)
  implicit none
  integer, intent(in) :: i
  integer, intent(out) :: ndigits
  integer, intent(out), dimension(*) :: digits
  integer :: j,isubtract,idigit

  ndigits = 0
  
  isubtract = 0      

  do j=size(digits),0,-1
     if (i.ge.10**j) then 
        ndigits = ndigits + 1
        idigit = (i - isubtract)/10**j
        digits(ndigits) = idigit
        isubtract = isubtract + digits(ndigits)*10**j
     endif
  enddo

end subroutine get_digits

end module transforms
