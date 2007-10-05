!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR GENERAL ASCII DATA FORMATS
!
! SOME CHOICES FOR THIS FORMAT CAN BE SET USING THE FOLLOWING
!  ENVIRONMENT VARIABLES:
!
! ASPLASH_NCOLUMNS can be used to override the automatic ncolumns choice
!
! e.g. setenv ASPLASH_NCOLUMNS=10
!
! the data is stored in the global array dat
!
! >> this subroutine must return values for the following: <<
!
! ncolumns    : number of data columns
! ndim, ndimV : number of spatial, velocity dimensions
! nstepsread  : number of steps read from this file
!
! dat(maxplot,maxpart,maxstep) : main data array
!
! npartoftype(1:6,maxstep) : number of particles of each type in each timestep
! ntot(maxstep)       : total number of particles in each timestep
!
! time(maxstep)       : time at each step
! gamma(maxstep)      : gamma at each step 
!
! most of these values are stored in global arrays 
! in the module 'particle_data'
!-------------------------------------------------------------------------

subroutine read_data(rootname,indexstart,nstepsread)
  use particle_data, only:dat,npartoftype,time,gamma,maxpart,maxcol,maxstep
  use params
  use settings_data, only:ndim,ndimV,ncolumns,ncalc
  use mem_allocation, only:alloc
  use asciiutils, only:get_ncolumns
  use system_utils, only:ienvironment
  implicit none
  integer, intent(in) :: indexstart
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname
  integer :: i,j,ierr,iunit,ncolstep,ncolenv
  integer :: nprint,npart_max,nstep_max,icol,nheaderlines
  logical :: iexist,timeset,gammaset
  real :: dummyreal
  character(len=len(rootname)+4) :: dumpfile

  nstepsread = 0
  nstep_max = 0
  npart_max = maxpart
  iunit = 15  ! logical unit number for input

  dumpfile = trim(rootname)
  !
  !--check if first data file exists
  !
  inquire(file=dumpfile,exist=iexist)
  if (.not.iexist) then
     print "(a)",' *** error: '//trim(dumpfile)//': file not found ***'    
     return
  endif
  !
  !--fix number of spatial dimensions (0 means no particle coords)
  !
  ndim = 0
  ndimV = 0

  j = indexstart
  nstepsread = 0
  
  write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)
  !
  !--open the file and read the number of particles
  !
  open(unit=iunit,iostat=ierr,file=dumpfile,status='old',form='formatted')
  if (ierr /= 0) then
     print "(a)",'*** ERROR OPENING '//trim(dumpfile)//' ***'
  else
     call get_ncolumns(iunit,ncolstep,nheaderlines)
     !--override columns setting with environment variable
     ncolenv = ienvironment('ASPLASH_NCOLUMNS',-1)
     if (ncolenv.gt.0) then
        print "(a,i3,a)",' setting ncolumns = ',ncolenv,' from ASPLASH_NCOLUMNS environment variable'
        ncolstep = ncolenv
     endif
     if (ncolstep.le.0) then
        print "(a)",'*** ERROR: zero/undetermined number of columns in file ***'
        return
     endif
     !
     !--allocate memory initially
     !
     nprint = 101
     nstep_max = max(nstep_max,indexstart,1)
     if (.not.allocated(dat) .or. (nprint.gt.npart_max) .or. (ncolstep+ncalc).gt.maxcol) then
        npart_max = max(npart_max,INT(1.1*(nprint)))
        call alloc(npart_max,nstep_max,ncolstep+ncalc)
     endif
  endif

  npart_max = max(npart_max,nprint)
  ncolumns = ncolstep
!
!--allocate/reallocate memory if j > maxstep
!
  if (j.gt.maxstep) then
     call alloc(maxpart,j+1,maxcol)
  endif
!
!--read header lines, try to use it to set time
!
  timeset = .false.
  gammaset = .false.
  do i=1,nheaderlines
     read(iunit,*,iostat=ierr) dummyreal
     if (timeset .and. .not.gammaset .and. ierr.eq.0 &
        .and. dummyreal.gt.0.999999 .and. dummyreal.lt.2.000001) then
        print*,'setting gamma = ',dummyreal,' from header line ',i
        gamma(j) = dummyreal
        gammaset = .true.
     endif
     if (ierr.eq.0 .and. .not. timeset) then
        time(j) = dummyreal
        timeset = .true.
        print*,'setting time = ',dummyreal,' from header line ',i
     endif
  enddo
!
!--now read the timestep data in the dumpfile
!
  i = 0
  ierr = 0
  overparts: do while (ierr == 0)
     i = i + 1
     if (i.gt.npart_max) then ! reallocate memory if necessary
        npart_max = 10*npart_max
        call alloc(npart_max,nstep_max,ncolstep+ncalc)
     endif
     read(iunit,*,iostat=ierr) (dat(i,icol,j),icol = 1,ncolstep)
  enddo overparts

  nprint = i - 1
  nstepsread = nstepsread + 1
  if (ierr < 0) then
     print*,'end of file: npts = ',nprint
  elseif (ierr > 0) then
     print*,' *** error reading file, npts = ',nprint,' ***'
  endif


  npartoftype(:,j) = 0
  npartoftype(1,j) = nprint

  close(iunit)
     
return
end subroutine read_data

!!-------------------------------------------------------------------
!! set labels for each column of data
!!
!! read these from a file called 'columns' in the current directory
!! then take sensible guesses as to which quantities are which
!! from the column labels
!!
!!-------------------------------------------------------------------

subroutine set_labels
  use labels, only:label,labeltype,ix,irho,ipmass,ih,iutherm,ipr,ivx,iBfirst,iamvec,labelvec
  use params
  use settings_data, only:ncolumns,ntypes,ndim,ndimV,UseTypeInRenderings
  use geometry, only:labelcoord
  implicit none
  integer :: i,ierr,ndimVtemp
!
!--read column labels from the columns file if it exists
!  
  open(unit=51,file='columns',status='old',iostat=ierr)
  if (ierr /=0) then
     print "(3(/,a))",' WARNING: columns file not found: using default labels',&
                    ' To change the labels, create a file called ''columns'' ',&
                    '  in the current directory with one label per line'
  else
     overcols: do i=1,ncolumns
        read(51,"(a)",iostat=ierr) label(i)
!
!--guess positions of various quantities from the column labels
!
        if (label(i)(1:3).eq.'den' .or. label(i)(1:3).eq.'rho') then
           irho = i
        elseif (label(i)(1:5).eq.'pmass' .or. label(i)(1:13).eq.'particle mass' &
                .or. index(label(i),'mass').ne.0) then
           ipmass = i
        !--use first column labelled h as smoothing length
        elseif (ih.eq.0 .and. (label(i)(1:1).eq.'h' &
                .or. label(i)(1:6).eq.'smooth')) then
           ih = i
        elseif (trim(label(i)).eq.'u'.or.label(i)(1:6).eq.'utherm' &
            .or.trim(label(i)).eq.'internal energy') then
           iutherm = i
        elseif (label(i)(1:2).eq.'pr') then
           ipr = i
        elseif (ivx.eq.0 .and. label(i)(1:1).eq.'v') then
           ivx = i
           ndimV = 1
        endif
        !--set ndimV as number of columns with v as label
        if (ivx.gt.0 .and. i.gt.ivx .and. i.le.ivx+2) then
           if (label(i)(1:1).eq.'v') ndimV = i - ivx + 1
        endif
        if (iBfirst.eq.0 .and. (label(i)(1:2).eq.'bx' .or. label(i)(1:2).eq.'Bx')) then
           iBfirst = i
        endif
        !--set ndimV as number of columns with v as label
        if (iBfirst.gt.0 .and. i.gt.iBfirst .and. i.le.iBfirst+2) then
           if (label(i)(1:1).eq.'b') then
              ndimVtemp = i - iBfirst + 1
              if (ndimV.gt.0 .and. ndimVtemp.gt.ndimV) then
                 print "(a)",' WARNING: possible confusion with vector dimensions'
                 ndimV = ndimVtemp
              endif
           endif
        endif        
        if (ierr < 0) then
           print "(a,i3)",' ERROR: end of file in columns file: read to column ',i-1
           exit overcols
        elseif (ierr > 0) then
           print "(a)",' *** error reading from columns file ***'
           exit overcols
        endif
     enddo overcols
     close(unit=51)
  endif
  
  if (label(1)(1:1).eq.'x') then
     ndim = 1
     ix(1) = 1
     if (label(2)(1:1).eq.'y') then
        ndim = 2
        ix(2) = 2
        if (label(3)(1:1).eq.'z') then
           ndim = 3
           ix(3) = 3
        endif
     endif
  endif
  if (ndim.gt.0) print "(a,i1)",' Assuming number of dimensions = ',ndim
  if (ndimV.gt.0) print "(a,i1)",' Assuming vectors have dimension = ',ndimV
  if (irho.gt.0) print "(a,i2)",' Assuming density in column ',irho
  if (ipmass.gt.0) print "(a,i2)",' Assuming particle mass in column ',ipmass
  if (ih.gt.0) print "(a,i2)",' Assuming smoothing length in column ',ih
  if (iutherm.gt.0) print "(a,i2)",' Assuming thermal energy in column ',iutherm
  if (ipr.gt.0) print "(a,i2)",' Assuming pressure in column ',ipr
  if (ivx.gt.0) then
     if (ndimV.gt.1) then
        print "(a,i2,a,i2)",' Assuming velocity in columns ',ivx,' to ',ivx+ndimV-1     
     else
        print "(a,i2)",' Assuming velocity in column ',ivx
     endif
  endif
  if (ndim.eq.0 .or. irho.eq.0 .or. ipmass.eq.0 .or. ih.eq.0) then
     print "(4(/,a))",' NOTE: Rendering capabilities cannot be enabled', &
                 '  until positions of density, smoothing length and particle', &
                 '  mass are known (for the ascii read the simplest way is to ', &
                 '  label the relevant columns appropriately in the columns file)'
  endif
  
  if (ivx.gt.0) then
     iamvec(ivx:ivx+ndimV-1) = ivx
     labelvec(ivx:ivx+ndimV-1) = 'v'
     do i=1,ndimV
       label(ivx+i-1) = 'v\d'//labelcoord(i,1)
     enddo
  endif
  if (iBfirst.gt.0) then
     iamvec(iBfirst:iBfirst+ndimV-1) = ivx
     labelvec(iBfirst:iBfirst+ndimV-1) = 'B'
     do i=1,ndimV
       label(iBfirst+i-1) = 'B\d'//labelcoord(i,1)
     enddo
  endif
  !
  !--set labels for each particle type
  !
  ntypes = 1 !!maxparttypes
  labeltype(1) = 'gas'
  UseTypeInRenderings(1) = .true.
  
 
!-----------------------------------------------------------

  return 
end subroutine set_labels
