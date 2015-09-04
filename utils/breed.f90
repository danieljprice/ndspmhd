!------------------------------------------
! Dan's super kernel-breeding program
! Written September 2015
! daniel.price@monash.edu
!------------------------------------------
! Implements a genetic algorithm to
! "breed" smoothing kernels towards
! some desired goal
!
! Code assumes an external script that
! produces the 20 tabulated "child"
! kernels in tables in the specified
! directory:
!  01.dat
!  02.dat
!  03.dat
!  ...
!
! Also in this directory should be a
! file called "ranked.list" consisting of
! a list of kernel numbers ranked by score
!
! By convention the kernel tables are 
! equally spaced in q2 for use in ndspmhd
!------------------------------------------
module genalg
 implicit none
 integer, parameter :: dp = 8
 real(dp), parameter :: mut_prob = 0.4
 real(dp), parameter :: mut_ampl = 1.0
 
 
 real(dp), parameter :: radkern = 2.
 real(dp), parameter :: radkern2 = radkern*radkern
 real(dp), parameter :: pi = 3.1415926536
 
contains

 !------------------------------------------------------
 ! join two kernels to make a new "child"
 ! we take a linear combination with random coefficient
 !------------------------------------------------------
 subroutine mate(ikern,wkern1,wkern2,wkern,iseed)
  use random, only:ran2
  integer, intent(in)  :: ikern,iseed
  real(dp),    intent(in)  :: wkern1(0:ikern),wkern2(0:ikern)
  real(dp),    intent(out) :: wkern(0:ikern)
  integer :: i
  real(dp) :: f
  
  f = ran2(iseed)
  !print*,'f = ',f
  do i=0,ikern
     wkern(i) = f*wkern1(i) + (1.-f)*wkern2(i)
  enddo

 end subroutine mate

 !------------------------------------------------------
 ! Mutation method 1: perturb W with Gaussian of random
 !  amplitude and wavelength
 !------------------------------------------------------
 subroutine mutate(ikern,wkern,iseed,iseed2)
  use random, only:ran2,rayleigh_deviate
  integer, intent(in)    :: ikern
  real(dp),intent(inout) :: wkern(0:ikern)
  integer, intent(inout) :: iseed,iseed2
  integer :: i
  real(dp) :: q0,h,ampl,q2,q,dq2table,f
  
  f = ran2(iseed)
  !print*,'f = ',f,iseed,iseed2
  if (f < mut_prob) then
     q0 = ran2(iseed)*radkern
     h = ran2(iseed)*radkern
     ampl = mut_ampl*rayleigh_deviate(iseed2)
     dq2table = radkern2/real(ikern,kind=dp)
     !print*,'MUTATING = ',q0,' h = ',h,' ampl = ',ampl

     do i=0,ikern
        q2 = i*dq2table
        q  = sqrt(q2)
        wkern(i) = wkern(i)*(1. + ampl*exp(-(q - q0)**2/(2.*h**2)))
     enddo
  endif

 end subroutine mutate

 !------------------------------------------------------
 ! Return Gaussian random deviate
 ! https://en.wikipedia.org/wiki/Box-Muller_transform
 !------------------------------------------------------
 real(dp) function gauss_dev(iseed)
  use random, only:ran2
  integer, intent(inout) :: iseed
  real(dp), parameter :: pi = 3.1415936536
  real(dp) :: u1,u2
  
  u1 = ran2(iseed)
  u2 = ran2(iseed)
  gauss_dev = sqrt(-2.*log(u1))*cos(2.*pi*u2)
  
 end function gauss_dev

 !------------------------------------------------------
 ! Mutation method 2: Mutate on 2nd derivative
 !  We then integrate to get effect on kernel
 !------------------------------------------------------
 subroutine mutate2(ikern,wkern,grkern,grgrwkern,iseed,iseed2)
  use random, only:ran2,rayleigh_deviate
  integer, intent(in)  :: ikern
  real(dp),intent(inout) :: wkern(0:ikern),grkern(0:ikern),grgrwkern(0:ikern)
  integer, intent(inout) :: iseed,iseed2
  integer :: i,mode
  real(dp), parameter :: pindex = 1.5
  real(dp) :: h,ampl,q2,q,dq2table,f,phi,kx,pow,k0,q0
  
  f = ran2(iseed)
  !print*,'f = ',f,iseed,iseed2
  if (f < mut_prob) then
     !q0 = ran2(iseed)*radkern
     phi = -pi + 2.*pi*ran2(iseed)
     h = ran2(iseed)*radkern
     kx = 2.*pi/h
     k0 = 2.*pi/radkern
     !mode = int(radkern/h) + 1
     !kx = mode*k0
     ! power law decay in amplitudes
     pow = 1./(kx/k0)**pindex
     ampl = mut_ampl*rayleigh_deviate(iseed2)*pow
     !! use gaussian deviate as mutation can be either + or -
     !ampl = mut_ampl*gauss_dev(iseed2)*pow
     dq2table = radkern2/real(ikern,kind=dp)
     !print*,'MUTATE = ',phi,'mode=',mode,' lambda = ',h,2.*pi/kx,' ampl = ',ampl
     do i=0,ikern
        q2 = i*dq2table
        q  = sqrt(q2)
        !grgrwkern(i) = grgrwkern(i)*(1. + ampl*exp(-(q - q0)**2/(2.*h**2)))
        grgrwkern(i) = grgrwkern(i)*(1. + ampl*sin(kx*q + phi))
     enddo
     !print*,' INTEGRATING...'
     ! integrate to get gradient
     call integrate(ikern,grgrwkern,grkern)
     ! integrate to get kernel
     call integrate(ikern,grkern,wkern)
  endif

 end subroutine mutate2

 !-----------------------------------------
 ! routine to normalise the kernel table
 ! appropriately in 1, 2 and 3D
 !-----------------------------------------
 subroutine normalise(ikern,wkern,c)
  integer, intent(in)    :: ikern
  real(dp),    intent(inout) :: wkern(0:ikern)
  real(dp),    intent(out)   :: c(3)
  real(dp) :: dq2table,q2,q,f(3),qprev,dq
  integer :: i
 
  dq2table = radkern2/real(ikern,kind=dp)

  ! trapezoidal rule to get integral under kernel
  f = 0.
  qprev = 0.
  do i=1,ikern
     q2 = i*dq2table
     q  = sqrt(q2)
     dq = q - qprev
     f(1) = f(1) + 0.5*(wkern(i) + wkern(i-1))*dq
     f(2) = f(2) + 0.5*(q*wkern(i) + qprev*wkern(i-1))*dq
     f(3) = f(3) + 0.5*(q2*wkern(i) + qprev**2*wkern(i-1))*dq
     qprev = q
  enddo
  f(1) = 2.*f(1)
  f(2) = 2.*pi*f(2)
  f(3) = 4.*pi*f(3)
  !print*,' integral is ',f(:)
  c(:) = 1./f(:)
  !print*,' normalisation is ',c(:)
 
 end subroutine normalise

 !-----------------------------------------------
 ! differentiate the kernel table to get
 ! tables for the derivative and 2nd deriv
 ! We use 2nd-order accurate finite differencing
 ! on non-uniform meshes
 !-----------------------------------------------
 subroutine diff(ikern,wkern,grkern,grgrkern)
  integer, intent(in)    :: ikern
  real(dp),    intent(in)    :: wkern(0:ikern)
  real(dp),    intent(out)   :: grkern(0:ikern),grgrkern(0:ikern)
  real(dp) :: dq2table,q2,q,qm1,qm2,qm3,qp1,qp2,qp3
  real(dp) :: dq,dqm1,dqm2,dqm3,dqp1,dqp2,a,b,c,d
  integer :: i
 
  dq2table = radkern2/real(ikern,kind=dp)

  ! finite differencing to get kernel derivatives
  do i=0,ikern
     q2 = i*dq2table
     q  = sqrt(q2)
     if (i==0 .or. i==1) then
     ! forward diff
        qp1 = sqrt((i+1)*dq2table)
        qp2 = sqrt((i+2)*dq2table)
        qp3 = sqrt((i+3)*dq2table)
        dq = qp1 - q
        dqp1 = qp2 - qp1
        dqp2 = qp3 - qp2
        a = -(2.*dq + dqp1)/(dq*(dq + dqp1))
        b = (dq + dqp1)/(dq*dqp1)
        c = -dq/(dqp1*(dq + dqp1))      
        grkern(i) = a*wkern(i) + b*wkern(i+1) + c*wkern(i+2)
     ! second deriv, forward diff
        a = (6.*dq + 4.*dqp1 + 2.*dqp2)/(dq*(dq + dqp1)*(dq + dqp1 + dqp2))
        b = -(4.*(dq + dqp1) + 2.*dqp2)/(dq*dqp1*(dqp1 + dqp2))
        c = (4.*dq + 2.*(dqp1 + dqp2))/((dqp1 + dq)*dqp1*dqp2)
        d = -(4.*dq + 2.*dqp1)/((dq + dqp1 + dqp2)*(dqp2 + dqp1)*dqp2)
        grgrkern(i) = a*wkern(i) + b*wkern(i+1) + c*wkern(i+2) + d*wkern(i+3)     
     elseif (i==ikern .or. i==ikern-1) then
     ! backward diff
        qm1 = sqrt((i-1)*dq2table)
        qm2 = sqrt((i-2)*dq2table)
        qm3 = sqrt((i-3)*dq2table)
        dqm2 = qm1 - qm2
        dqm1 = q - qm1
        dqm3 = qm2 - qm3
        a = dqm1/(dqm2*(dqm1+dqm2))
        b = -(dqm1 + dqm2)/(dqm1*dqm2)
        c = (2.*dqm1 + dqm2)/(dqm1*(dqm1 + dqm2))
        grkern(i) = a*wkern(i-2) + b*wkern(i-1) + c*wkern(i)
     ! second deriv, backwards
        a = -(4.*dqm1 + 2.*dqm2)/(dqm3*(dqm3 + dqm2)*(dqm3 + dqm2 + dqm1))
        b = (4.*dqm1 + 2.*(dqm2 + dqm3))/(dqm3*dqm2*(dqm2 + dqm1))
        c = -(4.*(dqm1 + dqm2) + 2.*dqm2)/((dqm2 + dqm3)*dqm2*dqm1)
        d = (6.*dqm1 + 4.*dqm2 + 2.*dqm3)/((dqm1 + dqm2 + dqm3)*(dqm1 + dqm2)*dqm1)
        grgrkern(i) = a*wkern(i-3) + b*wkern(i-2) + c*wkern(i-1) + d*wkern(i)
     else
        qp2 = sqrt((i+2)*dq2table)
        qp1 = sqrt((i+1)*dq2table)
        qm1 = sqrt((i-1)*dq2table)
        dqm1 = q - qm1
        dq   = qp1 - q
        dqp1 = qp2 - qp1
     ! 2nd order unequal grid finite diff
        grkern(i) = -dq/(dqm1*(dq+dqm1))*wkern(i-1) &
                   + (dq-dqm1)/(dq*dqm1)*wkern(i) &
                   + dqm1/(dq*(dq + dqm1))*wkern(i+1)
     ! 2nd deriv - BUT THIS IS ONLY FIRST ORDER ACCURATE
        !a = 2./(dqm1*(dq + dqm1))
        !b = -2./(dqm1*dq)
        !c = 2./(dq*(dq + dqm1))
        !grgrkern(i) = a*wkern(i-1) + b*wkern(i) + c*wkern(i+1)
     ! two nodes forward, one back 2nd derivative estimate
        a = 2.*(2.*dq + dqp1)/(dqm1*(dqm1 + dq)*(dqm1 + dq + dqp1))
        b = -2.*(2.*dq + dqp1 - dqm1)/(dqm1*dq*(dq + dqp1))
        c =  2.*(dq + dqp1 - dqm1)/((dqm1 + dq)*dq*dqp1)
        d = -2.*(dq - dqm1)/((dqm1 + dq + dqp1)*(dq + dqp1)*dqp1)
        grgrkern(i) = a*wkern(i-1) + b*wkern(i) + c*wkern(i+1) + d*wkern(i+2)
     endif
  enddo
 
 end subroutine diff

 !-----------------------------------------------
 ! Integrate the kernel table to get
 ! the kernel from the kernel derivative
 ! or kernel derivative from 2nd derivative.
 ! Called twice this does opposite of diff
 !-----------------------------------------------
 subroutine integrate(ikern,df,f)
  integer,     intent(in)    :: ikern
  real(dp),    intent(in)    :: df(0:ikern)
  real(dp),    intent(out)   :: f(0:ikern)
  real(dp) :: dq2table,q,qprev,dq
  integer :: i
  
  dq2table = radkern2/real(ikern,kind=dp)

  ! integrate to get kernel gradient
  f(ikern) = 0.
  qprev = sqrt(ikern*dq2table)
  do i=ikern-1,0,-1
     q  = sqrt(i*dq2table)
     dq = q - qprev
     f(i) = f(i+1) + 0.5*(df(i) + df(i+1))*dq
     qprev = q
  enddo

 end subroutine integrate
 
 !-----------------------------
 ! read kernel table from file
 !-----------------------------
 subroutine read_kernel(filename,ikern,wkern,ierr)
  character(len=*), intent(in) :: filename
  integer, intent(in)  :: ikern
  real(dp),    intent(out) :: wkern(0:ikern)
  integer, intent(out) :: ierr
  integer :: lu,i
  real(dp) :: c(3),dum
  
  open(newunit=lu,file=filename,status='old',iostat=ierr)
  if (ierr /= 0) return  
  read(lu,*) c(1:3)
  do i=0,ikern
     read(lu,*) dum,wkern(i)
  enddo
  close(lu)

 end subroutine read_kernel
 
 !----------------------------
 ! write kernel table to file
 !----------------------------
 subroutine write_kernel(filename,ikern,c,wkern,grkern,grgrkern,ierr)
  character(len=*), intent(in) :: filename
  integer, intent(in)  :: ikern
  real(dp),    intent(in)  :: wkern(0:ikern),grkern(0:ikern),grgrkern(0:ikern),c(3)
  integer, intent(out) :: ierr
  integer :: lu,i
  real(dp) :: q,dq2table

  dq2table = radkern2/real(ikern,kind=dp)  
  open(newunit=lu,file=filename,status='replace',iostat=ierr)
  if (ierr /= 0) return  
  write(lu,*) c(1:3)
  do i=0,ikern
     q = sqrt(i*dq2table)
     write(lu,*) q,wkern(i),grkern(i),grgrkern(i)
  enddo
  close(lu)

 end subroutine write_kernel
 
 !-------------------------------------------------------
 ! read two kernels and breed them to create a "child"
 ! call the mutation routine to give a possible mutation
 !-------------------------------------------------------
 subroutine breed_pair(file1,file2,fileout,iseed,iseed2)
  integer, parameter :: ikern = 4000
  character(len=*), intent(in)    :: file1,file2,fileout
  integer,          intent(inout) :: iseed,iseed2
  real(dp), dimension(0:ikern)    :: wkern1,wkern2,wkern,grkern,grgrkern
  real(dp) :: c(3)
  integer :: ierr1,ierr2,ierrw

  call read_kernel(file1,ikern,wkern1,ierr1)
  call read_kernel(file2,ikern,wkern2,ierr2)
  if (ierr1 /= 0 .or. ierr2 /= 0) stop 'error reading kernel files'

  call mate(ikern,wkern1,wkern2,wkern,iseed)
  !call mutate(ikern,wkern,iseed,iseed2)
  call diff(ikern,wkern,grkern,grgrkern)
  call mutate2(ikern,wkern,grkern,grgrkern,iseed,iseed2)

  call normalise(ikern,wkern,c)
  wkern = wkern*c(2)
  c(:) = c(:)/c(2)
  call diff(ikern,wkern,grkern,grgrkern)

  call write_kernel(fileout,ikern,c,wkern,grkern,grgrkern,ierrw)
  if (ierrw /= 0) print*,' ERROR during write'

 end subroutine breed_pair

 !--------------------------------------------
 ! procedure just to copy a kernel unchanged
 ! into the next generation
 !--------------------------------------------
 subroutine save_child(filein,fileout,ierr1)
  integer, parameter :: ikern = 4000
  character(len=*), intent(in)    :: filein,fileout
  real(dp), dimension(0:ikern)    :: wkern,grkern,grgrkern
  real(dp) :: c(3)
  integer :: ierr1,ierrw

  call read_kernel(filein,ikern,wkern,ierr1)
  if (ierr1 /= 0) stop 'error reading kernel files'

  call normalise(ikern,wkern,c)
  wkern = wkern*c(2)
  c(:) = c(:)/c(2)
  call diff(ikern,wkern,grkern,grgrkern)
  print "(3a)",trim(filein),' -> ',trim(fileout)

  call write_kernel(fileout,ikern,c,wkern,grkern,grgrkern,ierrw)
  if (ierrw /= 0) print*,' ERROR during write'

 end subroutine save_child

 !----------------------------------------------------------
 ! select pairs from the ranked list to use for breeding
 ! we do this probabilistically, giving higher weight
 ! to higher ranked kernels (currently prob \propto 1/rank)
 !----------------------------------------------------------
 subroutine select_pairs(dir,iseed,nkids,pairs,ierr,best)
  use random, only:ran2
  character(len=*), intent(in) :: dir
  integer, intent(inout) :: iseed
  integer, intent(in)    :: nkids
  integer, intent(out)   :: pairs(2,nkids),best(2)
  integer, intent(out)   :: ierr
  integer :: lu,i,i1,i2,npairs
  integer  :: label(nkids)
  real(dp) :: err(nkids),score(nkids),tot,myran,prob
  
  open(newunit=lu,file=trim(dir)//'/ranked.list',iostat=ierr)
  do i=1,nkids
     read(lu,*,iostat=ierr) label(i),err(i)
  enddo
  best(1) = label(1)
  best(2) = label(2)
  close(lu)
  if (ierr /= 0) return
  
  ! score is 1/err, normalise to get
  ! probability of breeding
  tot = 0.
  do i=1,nkids
     score(i) = 1./real(i) !1./exp(log10(err(i))) !1./real(i)**2 !1./err(i)**4
     tot      = tot + score(i)
  enddo
  score(:) = score(:)/tot
  
  !do i=1,nkids
  !   print "(i2,a,i2,a,f10.4)",i,' kernel = ',label(i),' prob = ',score(i)
  !enddo
  !print*,' total prob = ',sum(score)
  
  ! now roll the dice to see who breeds with who
  npairs = 0
  do while(npairs < nkids)
     ! pick two random partners
     i1 = int(nkids*ran2(iseed))+1
     i2 = int(nkids*ran2(iseed))+1
     myran = ran2(iseed)
     prob  = sqrt(score(i1)*score(i2))
     !print*,'trying ',i1,i2,myran,prob
     if (myran < prob .and. i1 /= i2) then
        npairs = npairs + 1
        pairs(1,npairs) = label(i1)
        pairs(2,npairs) = label(i2)
        !print*,' pair ',npairs,' breeding ',i1,i2
     endif
  enddo
  
 end subroutine select_pairs

end module genalg

!----------------------------------------------
! driver program for breeding procedure
! operates in 3 modes, depending on the number
! of arguments supplied on the command line:
!
!  1) given a single directory as argument
!     breeds all kernels in dir to give
!     the next generation
!  2) given a single file "kernel.dat"
!     just copies this to a new table
!     with appropriate normalisation applied
!  3) given two files e.g. k1.dat k2.dat
!     breeds a single pair to "kernel-new.dat"
!----------------------------------------------
program breed
 use genalg
 implicit none
 integer :: nargs,iseed,iseed2,lu,i,ierr
 character(len=120) :: file1,file2,fileout,dir
 logical  :: iexist
 integer, parameter :: nkids = 20
 integer :: pairs(2,nkids),best(2)
 
 iseed = -2345
 iseed2 = -125
 inquire(file='seed',exist=iexist)
 if (iexist) then
    open(newunit=lu,file='seed',status='old')
    read(lu,*) iseed,iseed2
    close(lu)
    !print*,' GOT ISEED = ',iseed,iseed2
 endif

 nargs = command_argument_count()
 if (nargs == 2) then
    call get_command_argument(1,file1)
    call get_command_argument(2,file2)
    fileout = 'kernel-new.dat'
    call breed_pair(file1,file2,fileout,iseed,iseed2)
    print "(5a)",trim(file1),' + ',trim(file2),' -> ',trim(fileout)
 elseif (nargs /= 1) then
    stop 'usage: breed dir'
 else
    call get_command_argument(1,dir)

    if (index(dir,'kernel')/=0) then
       ! just normalise an existing kernel
       call save_child(dir,'kernel-new.dat',ierr)
    else
       call select_pairs(dir,iseed,nkids,pairs,ierr,best)
       if (ierr /= 0) then
          print*,'ERROR: could not select pairs in '//trim(dir)
          stop
       endif

       do i=1,nkids-2
          write(file1,"(a,'/',i2.2,'.dat')") trim(dir),pairs(1,i)
          write(file2,"(a,'/',i2.2,'.dat')") trim(dir),pairs(2,i)
          write(fileout,"(i2.2,'.dat')") i
          !print "(5a)",trim(file1),' + ',trim(file2),' -> ',trim(fileout)
          call breed_pair(file1,file2,fileout,iseed,iseed2)
       enddo

       write(file1,"(a,'/',i2.2,'.dat')") trim(dir),best(1)
       write(fileout,"(i2.2,'.dat')") nkids-1
       call save_child(file1,fileout,ierr)

       write(file2,"(a,'/',i2.2,'.dat')") trim(dir),best(2)
       write(fileout,"(i2.2,'.dat')") nkids
       call save_child(file2,fileout,ierr)

       ! print*,' WRITING ',-abs(iseed),-abs(iseed2)

    endif
 endif

 ! update random seed
 open(newunit=lu,file='seed',status='replace')
 write(lu,*) -abs(iseed),-abs(iseed2)
 close(lu)

end program breed
