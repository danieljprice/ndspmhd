subroutine sobolsequence(ndim,x)
  implicit none
  integer, intent(in) :: ndim
  real, intent(out), dimension(abs(ndim)) :: x
  integer, parameter :: maxbit = 30, maxdim = 6
  !
  ! subroutine to generate a quasi-random sequence of points in n
  ! dimensions using Sobol sequences. Algorithm is from numerical 
  ! recipes but I have tried to make sense of the program and get
  ! rid of horrible things like EQUIVALENCE statements and GOTO's.
  !
  ! initial call should be with initialise = .true.
  ! this initialises a set of maxbit direction numbers for 
  ! up to maxdim different Sobol sequences.
  !
  ! with initialise = .false. the vector x(1..ndim) is returned containing
  ! the next values from ndim of these sequences.
  !
  ! Daniel Price, Institute of Astronomy, Cambridge UK, 27/7/04
  ! dprice@ast.cam.ac.uk
  !
  ! How the algorithm works is described in (e.g.) Numerical Recipes and
  ! in Antonov & Saleev (1979), USSR Computational Mathematics and 
  ! Computational Physics, vol 19, no.1, pp. 252-256
  ! 
  ! Each sequence (ie. we use one sequence per dimension) is based on 
  ! the next in the sequence of primitive polynomials modulo 2.
  ! This program works for up to 6 dimensions, for which the 
  ! first 6 primitive polynomials modulo 2 are
  ! x + 1             (ipolydeg = 1)
  ! x^2 + x + 1       (ipolydeg = 2) 
  ! x^3 + x + 1       (ipolydeg = 3)
  ! x^3 + x^2 + 1     (ipolydeg = 3)
  ! x^4 + x + 1       (ipolydeg = 4)
  ! x^4 + x^3 + 1     (ipolydeg = 4)
  !
  ! The coefficients of these polynomials are stored in a slightly obscure
  ! manner. Since they only be zero or one and since the coefficient of
  ! one and also the co-efficient of the highest order are always one, only
  ! the middle coefficients are stored. These middle coefficients are written
  ! as the bits of a binary number which are stored as a decimal.
  ! For example, for the last polynomial the coefficients are:
  !
  ! x^4 + 1*x^3 + 0*x^2 + 0*x + 1
  !
  ! so taking only the middle coefficients, the binary number is 100
  ! which in decimal is 4. This decimal is stored for each polynomial
  ! in the array ipolycoeff.
  !
  integer, dimension(maxdim), parameter :: ipolydeg = (/1,2,3,3,4,4/)
  integer, dimension(maxdim), parameter :: ipolycoeff = (/0,1,1,2,1,4/)
  integer, dimension(maxdim) :: ix
  integer, dimension(maxdim,maxbit) :: istart, ibitsequence
  integer :: icoeff,inuminseq
  integer :: i,j,k,iseq,ibit,index
  real :: factor
  save factor,inuminseq,ix,ibitsequence

  !
  !--istart contains the starting values for the recurrence (ie. the first
  !  four bits) 
  !  these can be arbitrary odd integers less than 2**deg where deg is the 
  !  degree of the starting polynomial for each sequence
  !  However, in the initialisations below, not all are freely specifiable
  !  as some are fixed by the recurrence (e.g. for degree 1, 1st bit is 1
  !  so 2nd, 3rd, 4th go 3,5,15 as set by the recurrence)
  !
  ibitsequence(:,1) = 1                  ! 1st bit for degree 1, so < 2
  ibitsequence(:,2) = (/3,1,3,3,1,1/)    ! 2nd bit for degree 2, so < 4
  ibitsequence(:,3) = (/5,7,7,3,3,5/)    ! 3rd bit for degree 3, so < 8
  ibitsequence(:,4) = (/15,1,5,15,13,9/) ! 4th bit for degree 4, so < 16 
  ibitsequence(:,5:maxbit) = 0

  !print*,'factor = ',factor,' number in sequence =',inuminseq
  !print*,'ibitsequence = ',ibitsequence

  if (ndim.le.0) then
     ix(:) = 0
     inuminseq = 0
     if (ibitsequence(1,1).ne.1) return
     factor = 1./2.**maxbit
     !print*,'fac = ',factor,maxbit,2**maxbit
     !
     !--initialise the maxdim Sobol sequences (nothing returned in x)
     !  this means calculating all of the direction numbers
     !
     do iseq = 1,maxdim
        !
        !--normalise the starting values of each sequence
        !  (these are those initialised above) to give the direction numbers
        !
        do j=1,ipolydeg(iseq)
           ibitsequence(iseq,j) = ibitsequence(iseq,j)*2**(maxbit-j)
        enddo
!        print*,'iseq = ',iseq,' istart = ',ibitsequence(iseq,:)
        !
        !--then calculate the other starting values (ie from degree->maxbit)
        !  using the recurrence relation
        !
        do j=ipolydeg(iseq)+1,maxbit
           icoeff = ipolycoeff(iseq)
           i = ibitsequence(iseq,j-ipolydeg(iseq))
           i = IEOR(i,i/2**ipolydeg(iseq))
           do k = ipolydeg(iseq)-1,1,-1
              if (IAND(icoeff,1).ne.0) i = IEOR(i,ibitsequence(iseq,j-1))
              icoeff = icoeff/2
           enddo
           ibitsequence(iseq,j) = i
        enddo
!        print*,'iseq = ',iseq,' istart = ',ibitsequence(iseq,:)
        
     enddo
  else
     !
     !--calculate the next vector in the sequence
     !
     ibit = inuminseq
     print*,'number in sequence = ',inuminseq
     !--find rightmost zero bit in inuminseq
     j = 1
     do while(IAND(ibit,1).NE.0 .and. j.lt.maxbit)
        ibit = ibit/2
        j = j + 1
     enddo
     print*,'rightmost zero bit = ',j-1,(j-1)*maxdim
     if (j.eq.maxbit) print*,'ERROR: maxbit too small in sobol sequence'
     ibit = (j-1)*maxdim
     !
     !--the inuminseq+1 th number in the sequence is obtained by an XOR of
     !  the previous number in the sequence with the direction number
     !  V_i with the position of the rightmost zero bit in inuminseq.
     !  (this is given by ibit)
     !
     !  the direction numbers are in the ibitsequence array. The sequence goes
     !  from 1->maxdim for every bit, so the 'i'th direction number is 
     !  given by adding (iseq-1)*maxdim + ibit
     !
     do iseq=1,min(ndim,maxdim)
        index = ibit + iseq
        print*,'index = ',index,' = ',index/maxdim,index-iseq*index/maxdim
        ix(iseq) = IEOR(ix(iseq),ibitsequence(index/maxdim,ibit))
        print*,'int x = ',index,ix(iseq)
!
!--convert the integer number to a floating point number for output
!
        x(iseq) = ix(iseq)*factor
     enddo
     inuminseq = inuminseq + 1

  endif
  return
end subroutine sobolsequence
