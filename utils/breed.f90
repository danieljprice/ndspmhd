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
 use kernel_utils, only:normalise,differentiate,diff,integrate
 implicit none
 integer, parameter :: dp = 8
 real(dp), parameter :: mut_prob = 0.4
 real(dp), parameter :: mut_ampl = 0.5
 
 
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
 subroutine mutate2(ikern,wkern,grkern,grgrwkern,iseed,iseed2,imethod)
  use random, only:ran2,rayleigh_deviate
  integer, intent(in)  :: ikern
  real(dp),intent(inout) :: wkern(0:ikern),grkern(0:ikern),grgrwkern(0:ikern)
  integer, intent(inout) :: iseed,iseed2
  integer :: i,imethod
  real(dp), parameter :: pindex = 1.5
  real(dp) :: h,ampl,q2,q,dq2table,f,phi,kx,pow,k0
  
  imethod = 0
  f = ran2(iseed)
  !print*,'f = ',f,iseed,iseed2
  if (f < mut_prob) then
     phi = -pi + 2.*pi*ran2(iseed)
     h = ran2(iseed)*radkern
     kx = 2.*pi/h
     k0 = 2.*pi/radkern
     !mode = int(radkern/h) + 1
     !kx = mode*k0
     ! power law decay in amplitudes
     pow = 1./(kx/k0)**pindex
     imethod = int(ran2(iseed)*2.) + 1
     select case(imethod)
     case(2)
        call smooth(ikern,grgrwkern)
     case default
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
     end select
     ! integrate to get gradient
     call integrate(ikern,grgrwkern,grkern,radkern2)
     ! integrate to get kernel
     call integrate(ikern,grkern,wkern,radkern2)
  endif

 end subroutine mutate2

 !------------------------------------------------------
 ! Mutation method 3: Mutate on 3rd derivative
 !  We then integrate to get effect on kernel
 !------------------------------------------------------
 subroutine mutate3(ikern,wkern,grkern,grgrwkern,iseed,iseed2,imethod,ampl)
  use random, only:ran2,rayleigh_deviate
  integer, intent(in)  :: ikern
  real(dp),intent(inout) :: wkern(0:ikern),grkern(0:ikern),grgrwkern(0:ikern)
  integer, intent(inout) :: iseed,iseed2
  integer, intent(out)   :: imethod
  real(dp), intent(out)  :: ampl
  integer :: i
  real(dp), parameter :: pindex = 1.0
  real(dp) :: h,q2,q,dq2table,f,phi,kx,pow,k0
  real(dp) :: gr3w(0:ikern),gr4w(0:ikern)
  
  f = ran2(iseed)
  imethod = 0
  ampl = 0.
  !print*,'f = ',f,iseed,iseed2
  if (f < mut_prob) then
     phi = -pi + 2.*pi*ran2(iseed)
     h = ran2(iseed)*radkern
     kx = 2.*pi/h
     k0 = 2.*pi/radkern
     !mode = int(radkern/h) + 1
     !kx = mode*k0
     ! power law decay in amplitudes
     pow = 1./(kx/k0)**pindex
     ampl = mut_ampl*rayleigh_deviate(iseed2)*pow

     dq2table = radkern2/real(ikern,kind=dp)
     ! get 3rd deriv
     call differentiate(ikern,grgrwkern,gr3w,radkern2)
     !call differentiate(ikern,gr3w,gr4w,radkern2)

!     call diff(ikern,grgrwkern,gr3w,gr4w,radkern2)    
     !print*,'MUTATE = ',phi,'mode=',mode,' lambda = ',h,2.*pi/kx,' ampl = ',ampl
     imethod = int(ran2(iseed)*4.) + 1
     
     if (imethod==4) then
        call smooth(ikern,grgrwkern)
     else
        do i=0,ikern
           q2 = i*dq2table
           q  = sqrt(q2)
           select case(imethod)
           case(3)
              grgrwkern(i) = grgrwkern(i)*(1 + 0.1*ampl*sin(kx*q + phi))
   !        case(2)
   !           gr4w(i) = gr4w(i)*(1 + 0.01*ampl*sin(kx*q + phi))
           case(2)
              gr3w(i) = gr3w(i) + ampl*sin(kx*q + phi)*(1. - q2/radkern2)**3
           case default
              gr3w(i) = gr3w(i)*(1 + 0.1*ampl*sin(kx*q + phi))
           end select
        enddo
        ! integrate to get 3rd deriv
        !if (imethod == 2) call integrate(ikern,gr4w,gr3w,radkern2)
        ! integrate to get 2nd deriv
        if (imethod < 3) call integrate(ikern,gr3w,grgrwkern,radkern2)
     endif
        ! integrate to get gradient
        call integrate(ikern,grgrwkern,grkern,radkern2)
        ! integrate to get kernel
        call integrate(ikern,grkern,wkern,radkern2)
  endif

 end subroutine mutate3
 
 !-----------------------------------------
 ! Apply a smoothing operator to a kernel
 ! We do this by solving a diffusion
 ! equation, evolving for several steps
 ! to eliminate high frequency noise
 !-----------------------------------------
 subroutine smooth(ikern,wkern)
  integer,  intent(in)    :: ikern
  real(dp), intent(inout) :: wkern(0:ikern)
  real(dp) :: dq2table,wnew(0:ikern) !,ddq22
  integer :: i,j

  dq2table = radkern2/real(ikern,kind=dp)
  !ddq22 = dq2table
  !dt = deltax^2/kappa, but coeff is 0.5*dt*kappa/deltax^2
  
  do j=1,10
     do i=0,ikern
        if (i==0) then
           wnew(i) = wkern(i+1) !+ 0.1*(2.*wkern(i) - 5.*wkern(i+1) + 4.*wkern(i+2) - wkern(i+3))
        elseif (i==ikern) then
           wnew(i) = wkern(i-1) ! + 0.1*(wkern(i-3) + 4.*wkern(i-2) - 5.*wkern(i-1) + 2.*wkern(i))     
        else
           wnew(i) = wkern(i) + 0.25*(wkern(i-1) - 2.*wkern(i) + wkern(i+1))
        endif
     enddo
     wkern = wnew
  enddo
  
 end subroutine smooth
 
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
  real(dp) :: q,dq2table,d2W

  dq2table = radkern2/real(ikern,kind=dp)  
  open(newunit=lu,file=filename,status='replace',iostat=ierr)
  if (ierr /= 0) return
  write(lu,*) c(1:3)
  do i=0,ikern
     q = sqrt(i*dq2table)
     if (q > 0.) then
        d2W = grgrkern(i) + 2.*grgrkern(i)/q
     else
        d2W = 0.
     endif
     write(lu,*) q,wkern(i),grkern(i),grgrkern(i),d2W
     !write(lu,*) q,wkern(i),grkern(i),grgrkern(i),d2W
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
  real(dp) :: c(3),ampl
  integer :: ierr1,ierr2,ierrw,imethod

  call read_kernel(file1,ikern,wkern1,ierr1)
  call read_kernel(file2,ikern,wkern2,ierr2)
  if (ierr1 /= 0 .or. ierr2 /= 0) stop 'error reading kernel files'

  call mate(ikern,wkern1,wkern2,wkern,iseed)
  !call mutate(ikern,wkern,iseed,iseed2)
  call diff(ikern,wkern,grkern,grgrkern,radkern2)
  ampl = 0.
  call mutate2(ikern,wkern,grkern,grgrkern,iseed,iseed2,imethod)
  !call mutate3(ikern,wkern,grkern,grgrkern,iseed,iseed2,imethod,ampl)

  call normalise(ikern,wkern,c,radkern2)
  wkern = wkern*c(3)
  c(:) = c(:)/c(3)
  call diff(ikern,wkern,grkern,grgrkern,radkern2)

  call write_kernel(fileout,ikern,c,wkern,grkern,grgrkern,ierrw)
  open(unit=100,file=trim(fileout)//'.mut',status='replace',form='formatted')
  write(100,*) imethod,ampl
  close(100)
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

  call normalise(ikern,wkern,c,radkern2)
  wkern = wkern*c(3)
  c(:) = c(:)/c(3)
!  call differentiate(ikern,wkern,grkern,radkern2)
!  call differentiate(ikern,grkern,grgrkern,radkern2)
  call diff(ikern,wkern,grkern,grgrkern,radkern2)
  !call diffu(ikern,wkern,grkern,grgrkern,radkern2)
  print "(3a)",trim(filein),' save-> ',trim(fileout)

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
    print*, 'usage: breed dir           (to breed all)'
    print*, '   or: breed kernel.dat    (to normalise a kernel table)'
    stop
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
