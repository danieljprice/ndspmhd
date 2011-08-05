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
! ivegotdata  : flag which indicates successful data read
!
! maxplot,maxpart,maxstep      : dimensions of main data array
! dat(maxplot,maxpart,maxstep) : main data array
!
! npartoftype(1:6,maxstep) : number of particles of each type in each timestep
! ntot(maxstep)       : total number of particles in each timestep
! iam(maxpart,maxstep): integer identification of particle type
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
  use settings_data
  implicit none
  integer, intent(IN) :: indexstart
  integer, intent(OUT) :: nstepsread
  character(LEN=*), intent(IN) :: rootname
  integer, parameter :: maxptmass = 1000
  integer :: i,j,ifile,ierr
  integer :: npart_max,nstep_max
  logical :: iexist
    
  character(LEN=3) :: fileno
  character(LEN=LEN(rootname)+10) :: dumpfile
  integer :: nprint, nghosti, n1, n2, nptmass
  integer, dimension(:), allocatable :: isteps, iphase
  integer, dimension(maxptmass) :: listpm
  
  real(doub_prec), dimension(:,:), allocatable :: dattemp
  real(doub_prec), dimension(:), allocatable :: dummy
  real(doub_prec) :: udisti,umassi,utimei, umagfdi, timei, gammai
  real(doub_prec) :: rhozero, RK2
  real(doub_prec) :: escap,tkin,tgrav,tterm,tmag
  real(doub_prec) :: dtmax

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
  ivegotdata = .false.
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
  ncolumns = 19  ! number of columns in file  
  !
  !--allocate memory initially
  !
  nstep_max = max(nstep_max,indexstart,2)

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
        read(15,end=55,iostat=ierr) udisti,umassi,utimei,umagfdi,nprint 
        if (.not.allocated(dat) .or. nprint.gt.npart_max) then
           npart_max = max(npart_max,INT(1.1*nprint))
           call alloc(npart_max,nstep_max,ncolumns)
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
           call alloc(maxpart,j+1,maxcol)
        endif
!
!--allocate a temporary array for double precision variables
!
        if (allocated(dattemp)) deallocate(dattemp)
        allocate(dattemp(npart_max,ncolumns))
!
!--allocate a dummy arrays for data I want to throw away
!
        if (allocated(dummy)) deallocate(dummy)
        allocate(dummy(npart_max))
        if (allocated(isteps)) deallocate(isteps)
        allocate(isteps(npart_max))
        if (allocated(iphase)) deallocate(iphase)
        allocate(iphase(npart_max))
!
!--now read the timestep data in the dumpfile
!
        read(15,end=55,iostat=ierr) udisti, umassi, utimei, umagfdi,  &
             nprint, nghosti, n1, n2, timei, gammai, rhozero, RK2, &
             (dattemp(i,7), i=1, nprint), (dattemp(i,8), i=1,nprint), &
             escap, tkin, tgrav, tterm, tmag, &
             (dattemp(i,1), i=1, nprint), (dattemp(i,2), i=1, nprint), &
             (dattemp(i,3), i=1, nprint), (dattemp(i,4), i=1, nprint), &
             (dattemp(i,5), i=1, nprint), (dattemp(i,6), i=1, nprint), &
             (dattemp(i,9), i=1, nprint), (dattemp(i,10), i=1, nprint), &
             (dattemp(i,11), i=1, nprint), (dattemp(i,12), i=1, nprint), &  
             (dattemp(i,13), i=1, nprint), (dattemp(i,14), i=1, nprint), &
             (dattemp(i,15), i=1, nprint), (dattemp(i,16), i=1, nprint), &
             (dattemp(i,17), i=1, nprint), (dattemp(i,18), i=1, nprint), &
             (dattemp(i,19), i=1, nprint), (dummy(i),i=1,nprint), &
             dtmax, (isteps(i), i=1,nprint), (iphase(i),i=1,nprint), &
             nptmass, (listpm(i), i=1,nptmass)
        
        if (ierr /= 0) then
           print "(a)",'|*** ERROR READING TIMESTEP ***'
           cycle over_steps_in_file
        else
           nstepsread = nstepsread + 1
        endif
!
!--convert to single precision
!
        print "(a,i5,a,i8)",'| step ',j,': ntotal = ',nprint
        print "(a)",'| converting to single precision... '
        dat(1:nprint,1:ncolumns,j) = real(dattemp(1:nprint,1:ncolumns))

        iam(1:nprint,j) = iphase(1:nprint)
        deallocate(dattemp)
        deallocate(dummy)
        deallocate(isteps)
        deallocate(iphase)

        ntot(j) = nprint
        npartoftype(1,j) = nprint-nghosti
        npartoftype(2,j) = nghosti

        gamma(j) = real(gammai)
        time(j) = real(timei)
        if (npart_max.eq.130000) npart_max = nprint
        j = j + 1

     enddo over_steps_in_file
     
     endif

55 continue
  !
  !--reached end of file
  !
  close(15)

  print*,'>> end of dump file: nsteps =',j-1,'ntot = ',ntot(j-1),'nghost=',npartoftype(2,j-1)

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
   
ivegotdata = .true.
return
                    
end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
  use labels
  use params
  use settings_data
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
  irho = 18     ! location of rho in data array
  iutherm = 16  !  thermal energy
  ih = 7        !  smoothing length
  ipmass = 17   !  particle mass      
  
  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  do i=1,ndimV
     label(ivx+i-1) = 'v\d'//labelcoord(i,1)
  enddo
  label(irho) = '\gr'      
  label(iutherm) = 'u'
  label(ih) = 'h       '
  label(ipmass) = 'particle mass'     
  label(8) = 'alpha'
  label(19) = 'psi' 
  
  label(ndim + ndimV+5) = '\ga'
  if (ncolumns.gt.11) then
     iBfirst = 9    ! location of Bx
     do i=1,ndimV
        label(iBfirst + i-1) = 'B\d'//labelcoord(i,1) !' (x10\u-3\d)'	!//'/rho'
     enddo
     idivB = 12
     label(idivB) = 'div B'
     do i=1,ndimV
        label(13 + i-1) = 'J'//labelcoord(i,1)
     enddo
  else
     iBfirst = 0
  endif
  !
  !--set labels for vector quantities
  !
  iamvec(ivx:ivx+ndimV-1) = ivx
  labelvec(ivx:ivx+ndimV-1) = 'v'
  do i=1,ndimV
     label(ivx+i-1) = trim(labelvec(ivx))//'\d'//labelcoord(i,1)
  enddo
  !--mag field
  if (iBfirst.gt.0) then
     iamvec(iBfirst:iBfirst+ndimV-1) = iBfirst
     labelvec(iBfirst:iBfirst+ndimV-1) = 'B'
     do i=1,ndimV
        label(iBfirst+i-1) = trim(labelvec(iBfirst))//'\d'//labelcoord(i,1)
     enddo
  endif
  !--current density
  iamvec(13:13+ndimV-1) = 13
  labelvec(13:13+ndimV-1) = 'J'
  do i=1,ndimV
     label(13+i-1) = trim(labelvec(13))//'\d'//labelcoord(i,1)
  enddo
  
  !
  !--set labels for each particle type
  !
  ntypes = 3  !!maxparttypes
  labeltype(1) = 'gas'
  labeltype(2) = 'ghost'
  labeltype(3) = 'sink'
 
!-----------------------------------------------------------

  return 
end subroutine set_labels
