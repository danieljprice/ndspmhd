module dumpfiles
 implicit none
 
contains
!!--------------------------------------------------------------------
!!  subroutine to write output to data file
!!  (change this to change format of output)
!!-------------------------------------------------------------------

subroutine write_dump(t,dumpfile)
 use dimen_mhd
 use loguns
 
 use eos
 use options
 use part
 use setup_params
 use derivb
 use rates
!
!--define local variables
!
 implicit none
 real, intent(in) :: t
 character(len=*), intent(in) :: dumpfile
 integer :: nprint,ndata
 integer :: i,ierr,iformat

 if (idumpghost.eq.1) then
    nprint = ntotal
 else
    nprint = npart
 endif
!
!--open dumpfile
!
 open(unit=idatfile,file=dumpfile,status='replace',form='unformatted',iostat=ierr)
 if (ierr /= 0) then
    write(iprint,*) 'error: can''t create new dumpfile ',trim(dumpfile)
    stop
 endif
!
!--write timestep header to data file
!
 if (imhd.ne.0) then
    ndata = ndim + 11 + 3*ndimV ! number of columns apart from co-ords
 else
    ndata = ndim + 8 + ndimV
    if (igeom.ne.0) ndata = ndata + 2 + ndimV
 endif
 write(idatfile,iostat=ierr) t,npart,nprint,gamma,hfact,ndim,ndimV,ndata,igeom,iformat
 if (ierr /= 0) then
    write(iprint,*) '*** error writing timestep header to dumpfile ',trim(dumpfile)
 endif
!
!--write the data (primitive variables) to the data file
!  data is written in two blocks :
! 
! 1) essential variables needed to restart the code
! 2) variables that are output for information only
!
! MHD variables are written after the hydro ones in 
! each case
!
  if (imhd.ne.0) then
     !--essential variables
     write(idatfile) x(:,1:nprint)
     write(idatfile) vel(:,1:nprint)
     write(idatfile) hh(1:nprint)
     write(idatfile) dens(1:nprint)
     write(idatfile) uu(1:nprint)
     write(idatfile) pmass(1:nprint)
     write(idatfile) alpha(1:3,1:nprint)
     write(idatfile) Bfield(:,1:nprint)
     write(idatfile) psi(1:nprint)
     !--info only
     write(idatfile) pr(1:nprint)
     write(idatfile) -drhodt(1:nprint)/rho(1:nprint)
     write(idatfile) divB(1:nprint)
     write(idatfile) curlB(:,1:nprint)
  else
     !--essential variables
     write(idatfile) x(:,1:nprint)
     write(idatfile) vel(:,1:nprint)
     write(idatfile) hh(1:nprint)
     write(idatfile) dens(1:nprint)
     write(idatfile) uu(1:nprint)
     write(idatfile) pmass(1:nprint)
     write(idatfile) alpha(1:2,1:nprint)
     !--info only
     write(idatfile) pr(1:nprint)
     write(idatfile) -drhodt(1:nprint)/rho(1:nprint)
     if (igeom.ne.0) then
        write(idatfile) rho(1:nprint)
        write(idatfile) sqrtg(1:nprint)
        write(idatfile) pmom(:,1:nprint)
     endif
  endif

 close(unit=idatfile)

 return
end subroutine write_dump

!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  setup by reading a dump from a file                                   !!
!!  should be compatible with the output format given in write_dump       !!
!!                                                                        !!
!!------------------------------------------------------------------------!!

subroutine read_dump(dumpfile,tfile)
!
!--include relevant global variables
!
 use dimen_mhd
 use debug
 use loguns

 use bound
 use derivb
 use eos
 use options
 use part
 use setup_params
 use geometry
!
!--define local variables
!      
 implicit none
 character(len=*), intent(in) :: dumpfile
 real, intent(out) :: tfile
 integer :: i,j,ndimfile,ndimvfile,npartfile,nprintfile,ncolumns
 integer :: ierr, iformat,igeomfile
 real :: gammafile,hfactfile
 logical :: iexist, mhdfile
 real, dimension(ndim) :: xnew

 igeomfile = 0
!
!--set name of setup file
!
 write(iprint,*) 'reading setup from file: ',trim(dumpfile)
 inquire (file=dumpfile,exist=iexist)
 if (.not.iexist) then
    write(iprint,*) 'can''t find dump file: ',trim(dumpfile)    
    stop
 endif
!
!--open setup file
!
 open(unit=ireadf,file=dumpfile,status='old',form='unformatted',iostat=ierr)
 if (ierr /= 0) then
    write(iprint,*) 'error opening dump file: ',trim(dumpfile)
    stop
 endif
!
!--read header line
!
! read(ireadf,iostat=ierr) tfile,npartfile,nprintfile,gammafile,hfactfile, &
!      ndimfile,ndimvfile,ncolumns,igeomfile,iformat
 read(ireadf,iostat=ierr) tfile,npartfile,nprintfile,gammafile,hfactfile, &
      ndimfile,ndimvfile,ncolumns
 if (ierr /= 0) then
    write(iprint,*) 'error reading header from dump file: ',trim(dumpfile)
    stop 
 endif
 write(iprint,*) 'time = ',tfile,' in dump file '
!
!--check for compatibility with current settings
!
 if (ndimfile.ne.ndim) then
    write(iprint,*) '***ERROR: x dimensions not equal between dump file and code',ndim,ndimfile
    if (ndimfile.le.0 .or. ndimfile.ge.3) then
       write(iprint,*) 'This could be because the code was compiled in double precision'
       write(iprint,*) ' whilst moddump has been compiled in single (or vice versa)' 
    endif
    stop
 endif
 if (ndimVfile.ne.ndimV) stop 'v dimensions not equal between dump file and code'
 if (abs(gammafile-gamma).gt.1.e-3) write(iprint,10) 'gamma',gammafile,gamma
 if (abs(hfactfile-hfact).gt.1.e-3) write(iprint,10) 'hfact',hfactfile,hfact
10 format('warning: ',a,' changed from original setup: old = ',f9.6,' new = ',f9.6)

 if ((nprintfile.ne.npartfile).and.(all(ibound.le.1))) then
    write(iprint,*) 'warning: setup file contains ghosts, but none set'
 endif

 if (ncolumns.ge.(ndim+6+ndimV)) then
    if (ncolumns.lt.ndim+2*ndimV+8 .or. igeomfile.ne.0) then
       mhdfile = .false.
       write(iprint,*) '(non-mhd input file)'
       if (imhd.gt.0) write(iprint,*) 'warning: reading non-mhd file, but mhd is on'
    else
       mhdfile = .true.
       write(iprint,*) '(mhd input file)'
       if (imhd.le.0) write(iprint,*) 'warning: mhd input file, but mhd is off'       
    endif
 else
    mhdfile = .false.
    write(iprint,*) 'WARNING: header suggests not enough data columns in input file'
 endif
!
!--allocate memory for the number of particles
!
 npart = npartfile
 ntotal = npartfile

 call alloc(nprintfile)
!
!--read data from file (only bits needed to restart the run - do not read ghosts)
!
 read(ireadf,iostat=ierr) (x(1:ndim,i),i=1,npart)
 if (ierr /= 0) then
    write(iprint,*) '*** error reading data from dumpfile ***'
    stop
 endif
 read(ireadf,iostat=ierr) vel(1:ndimV,1:npart)
 if (ierr /= 0) stop 'error reading dumpfile'
 read(ireadf,iostat=ierr) hh(1:npart)
 if (ierr /= 0) stop 'error reading dumpfile'
 read(ireadf,iostat=ierr) dens(1:npart)
 if (ierr /= 0) stop 'error reading dumpfile'
 read(ireadf,iostat=ierr) uu(1:npart)
 if (ierr /= 0) stop 'error reading dumpfile'
 read(ireadf,iostat=ierr) pmass(1:npart)
 if (ierr /= 0) stop 'error reading dumpfile'
 if (mhdfile .and. imhd.ne.0) then
    read(ireadf,iostat=ierr) alpha(1:3,1:npart)
    if (ierr /= 0) then
       write(iprint,*) '*** error in dumpfile : non-MHD file ***'
    endif
    read(ireadf,iostat=ierr) Bfield(1:ndimV,1:npart)
    read(ireadf,iostat=ierr) psi(1:npart)
 elseif (mhdfile) then
    !--read alpha in MHD format
    read(ireadf,iostat=ierr) alpha(1:3,1:npart)    
 else
    read(ireadf,iostat=ierr) alpha(1:2,1:npart)
 endif
!
!--close the file
!
 close(unit=ireadf)
 
!
!--convert to appropriate coordinate system
!
 if (igeomfile.ne.igeom) then
    write(iprint,*) 'CONVERTING file from coord system ',igeomfile,' to ',igeom
    do i=1,npart
       call coord_transform(x(:,i),ndim,igeomfile,xnew(:),ndim,igeom)
       x(:,i) = xnew(:)
    enddo
    write(iprint,*) 'WARNING: vector transform not yet implemented'
 endif
 
 
 write(iprint,*) 'finished reading setup file: everything is aok'

!
!--set bounds of setup
!                   
 xmax = 0.   ! irrelevant
 xmin = 0.
 ibound = 0	! no boundaries 
 iexternal_force = 1	! use toy star force

 if (igeom.eq.2) then
    ibound(1) = 2 ! reflective in r
    xmin(1) = 0.0
    if (ndim.ge.2) then 
       ibound(2) = 3 ! periodic in phi
       xmin(2) = -pi ! phi min
       xmax(2) = pi ! phi max
    endif
 endif
!
!--now change things according to the specific setup required
!
! vel(:,:) = 0.

 return

end subroutine read_dump

end module dumpfiles

