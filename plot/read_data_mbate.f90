!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR READING UNFORMATTED OUTPUT FROM MATTHEW BATE'S CODE
! (ie. STRAIGHT FROM THE DATA DUMP)
!
! *** CONVERTS TO SINGLE PRECISION ***
!
! the data is stored in the global array dat
!
! >> this subroutine must return values for the following: <<
!
! ncolumns    : number of data columns
! ndim, ndimV : number of spatial, velocity dimensions
! nstepsread  : number of steps read from this file
!
! maxplot,maxpart,maxstep      : dimensions of main data array
! dat(maxplot,maxpart,maxstep) : main data array
!
! npartoftype(1:6,maxstep) : number of particles of each type in each timestep
!
! time(maxstep)       : time at each step
! gamma(maxstep)      : gamma at each step 
!
! most of these values are stored in global arrays 
! in the module 'particle_data'
!-------------------------------------------------------------------------

subroutine read_data(rootname,indexstart,nstepsread)
  use particle_data
  use params
  use settings_data, only:ndim,ndimV,ncolumns,ncalc
  use mem_allocation
  implicit none
  integer, intent(in) :: indexstart
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname
  integer, parameter :: maxptmass = 1000
  real, parameter :: pi=3.141592653589
  integer :: i,j,ifile,ierr
  integer :: npart_max,nstep_max,ncolstep
  logical :: iexist
    
  character(len=3) :: fileno
  character(len=len(rootname)+10) :: dumpfile
  integer :: nprint, nghosti, n1, n2, nptmass
  integer, dimension(:), allocatable :: isteps, iphase
  integer, dimension(maxptmass) :: listpm
  
  !--use these lines if dump is double precision
  real(doub_prec), dimension(:,:), allocatable :: dattemp
  real(doub_prec), dimension(:), allocatable :: dummy
  real(doub_prec) :: udisti,umassi,utimei
  real(doub_prec) :: timei, gammai
  real(doub_prec) :: rhozero, RK2
  real(doub_prec) :: escap,tkin,tgrav,tterm
  real(doub_prec) :: dtmax, tcomp

  !--use these lines for single precision
  !real, dimension(:,:), allocatable :: dattemp
  !real, dimension(:), allocatable :: dummy
  !real(doub_prec) :: udisti,umassi,utimei
  !real :: timei, gammai
  !real :: rhozero, RK2
  !real :: escap,tkin,tgrav,tterm
  !real :: dtmax,tcomp

  nstepsread = 0
  nstep_max = 0
  npart_max = maxpart
  ifile = 1
  !
  !--for rootnames without the '00', read all files starting at #1
  !
  if (len_trim(rootname).lt.7) then
     ifile = 1
     if (len_trim(rootname).eq.4) then
        write(fileno,"(i1,i1,i1)") ifile/100,mod(ifile,100)/10,mod(ifile,10)
        dumpfile = rootname(1:4)//fileno 
     elseif (len_trim(rootname).eq.5) then
        write(fileno,"(i1,i1)") ifile/10,mod(ifile,10)
        dumpfile = rootname(1:5)//trim(fileno)     
     endif
  else
     dumpfile = trim(rootname)   
  endif
  !
  !--check if first data file exists
  !
  inquire(file=dumpfile,exist=iexist)
  if (.not.iexist) then
     print "(a)",' *** error: ',trim(dumpfile),' file not found ***'    
     return
  endif
  !
  !--fix number of spatial dimensions
  !
  ndim = 3
  ndimV = 3
  ncolstep = 10  ! number of columns in file
  ncolumns = ncolstep
  !
  !--allocate memory initially
  !
  nstep_max = max(nstep_max,indexstart,4)

  j = indexstart
  nstepsread = 0
  
  do while (iexist)
     write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)
     !
     !--open the (unformatted) binary file and read the number of particles
     !
     open(unit=15,iostat=ierr,file=dumpfile,status='old',form='unformatted')
     if (ierr /= 0) then
        print*,'*** ERROR OPENING ',trim(dumpfile),' ***'
     else
        !
        !--read the number of particles in the first step,
        !  allocate memory and rewind
        !
        read(15,end=55,iostat=ierr) udisti,umassi,utimei,nprint 
        if (.not.allocated(dat) .or. nprint.gt.npart_max) then
           npart_max = max(npart_max,INT(1.1*nprint))
           call alloc(npart_max,nstep_max,ncolstep+ncalc)
        endif
        rewind(15)
     endif
     if (ierr /= 0) then
        print*,'*** ERROR READING TIMESTEP HEADER ***'
     else
!
!--loop over the timesteps in this file
!     
     over_steps_in_file: do     
        npart_max = max(npart_max,nprint)
!
!--allocate/reallocate memory if j > maxstep
!
        if (j.gt.maxstep) then
           call alloc(maxpart,j+2*nstepsread,maxcol)
        endif
!
!--allocate a temporary array for double precision variables
!
        if (allocated(dattemp)) deallocate(dattemp)
        allocate(dattemp(npart_max,ncolstep),stat=ierr)
        if (ierr /= 0) print*,'not enough memory in read_data'
!
!--allocate a dummy arrays for data I want to throw away
!
        if (allocated(dummy)) deallocate(dummy)
        allocate(dummy(npart_max),stat=ierr)
        if (ierr /= 0) print*,'not enough memory in read_data'

        if (allocated(isteps)) deallocate(isteps)
        allocate(isteps(npart_max),stat=ierr)
        if (ierr /= 0) print*,'not enough memory in read_data'

        if (allocated(iphase)) deallocate(iphase)
        allocate(iphase(npart_max),stat=ierr)
        if (ierr /= 0) print*,'not enough memory in read_data'

!
!--now read the timestep data in the dumpfile
!
        write(*,"(a,i5,a)",advance="no") '| step ',j,': '

        read(15,end=55,iostat=ierr) udisti, umassi, utimei, &
             nprint, n1, n2, timei, gammai, rhozero, RK2, &
             (dattemp(i,7), i=1, nprint),escap, tkin, tgrav, tterm, &
             (dattemp(i,1), i=1, nprint), (dattemp(i,2), i=1, nprint), &
             (dattemp(i,3), i=1, nprint), (dattemp(i,4), i=1, nprint), &
             (dattemp(i,5), i=1, nprint), (dattemp(i,6), i=1, nprint), &
             (dattemp(i,8), i=1, nprint), (dattemp(i,9), i=1, nprint), &
             (dattemp(i,10), i=1, nprint), (dummy(i),i=1,nprint), &
             dtmax, (isteps(i), i=1,nprint), (iphase(i),i=1,nprint), &
             nptmass, (listpm(i), i=1,nptmass)
        
        if (ierr /= 0) then
           print "(a)",'*** INCOMPLETE DATA ON LAST TIMESTEP ***'
           nstepsread = nstepsread + 1
           exit over_steps_in_file
        else
           nstepsread = nstepsread + 1
        endif
!
!--convert to single precision
!
        print *,'t = ',timei,' ntotal = ',nprint
        if (nptmass.gt.0) then
           print *,'WARNING: nptmasses = ',nptmass,' but nothing done in read_data'
        endif
        print "(a)",'| converting to single precision... '
        dat(1:nprint,1:ncolstep,j) = real(dattemp(1:nprint,1:ncolstep))

!
!--convert to physical units
!
        dat(1:nprint,10,j) = dat(1:nprint,10,j)*real(umassi/udisti**3)
        if (allocated(dattemp)) deallocate(dattemp)
        if (allocated(dummy)) deallocate(dummy)
        if (allocated(isteps)) deallocate(isteps)
        if (allocated(iphase)) deallocate(iphase)

        npartoftype(1,j) = nprint-nghosti
        npartoftype(2,j) = nghosti

        gamma(j) = real(gammai)
        tcomp = sqrt((3.*pi)/(32*rhozero))
        time(j) = real(timei)/tcomp
        j = j + 1

     enddo over_steps_in_file
     
     endif

55 continue
  !
  !--reached end of file
  !
  close(15)
  if (j-1 .gt. 0) then
     print*,'>> end of dump file: nsteps =',j-1,'ntot = ', &
           sum(npartoftype(:,j-1)),'nghost=',npartoftype(2,j-1)
  endif
     !
     !--if just the rootname has been input, 
     !  set next filename and see if it exists
     !
  ifile = ifile + 1
  if (len_trim(rootname).eq.4) then
     write(fileno,"(i1,i1,i1)") ifile/100,mod(ifile,100)/10,mod(ifile,10)
     dumpfile = rootname(1:4)//fileno 
     inquire(file=dumpfile,exist=iexist)
     elseif (len_trim(rootname).eq.5) then
     write(fileno,"(i1,i1)") ifile/10,mod(ifile,10)
     dumpfile = rootname(1:5)//trim(fileno)     
     inquire(file=dumpfile,exist=iexist)
  else
     iexist = .false. ! exit loop
  endif
enddo
   
return
                    
end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
  use labels
  use params
  use settings_data
  use geometry, only:labelcoord
  implicit none
  integer :: i
  
  if (ndim.le.0 .or. ndim.gt.3) then
     print*,'*** ERROR: ndim = ',ndim,' in set_labels ***'
     return
  endif
  if (ndimV.le.0 .or. ndimV.gt.3) then
     print*,'*** ERROR: ndimV = ',ndimV,' in set_labels ***'
     return
  endif
    
  do i=1,ndim
     ix(i) = i
  enddo
  ivx = 4
  ih = 7        !  smoothing length
  iutherm = 8  !  thermal energy
  ipmass = 9   !  particle mass      
  irho = 10     ! location of rho in data array
  
  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  do i=1,ndimV
     label(ivx+i-1) = 'v\d'//labelcoord(i,1)
  enddo
  label(irho) = 'density (g/cm\u3\d)'
  label(iutherm) = 'u'
  label(ih) = 'h       '
  label(ipmass) = 'particle mass'     

    !
  !--set labels for vector quantities
  !
  iamvec(ivx:ivx+ndimV-1) = ivx
  labelvec(ivx:ivx+ndimV-1) = 'v'
  do i=1,ndimV
     label(ivx+i-1) = trim(labelvec(ivx))//'\d'//labelcoord(i,1)
  enddo
  !
  !--set labels for each particle type
  !
  ntypes = 2  !!maxparttypes
  labeltype(1) = 'gas'
  !!labeltype(2) = 'ghost'
  labeltype(2) = 'sink'
 
!-----------------------------------------------------------

  return 
end subroutine set_labels
