!------------------------------------------------------------------------------!
! NDSPMHD: A Smoothed Particle (Magneto)Hydrodynamics code for (astrophysical) !
! fluid dynamics simulations in 1, 2 and 3 spatial dimensions.                 !
!                                                                              !
! (c) 2002-2014 Daniel Price                                                   !
!                                                                              !
! http://users.monash.edu.au/~dprice/ndspmhd                                   !
! daniel.price@monash.edu -or- dprice@cantab.net (forwards to current address) !
!                                                                              !
!  NDSPMHD comes with ABSOLUTELY NO WARRANTY.                                  !
!  This is free software; and you are welcome to redistribute                  !
!  it under the terms of the GNU General Public License                        !
!  (see LICENSE file for details) and the provision that                       !
!  this notice remains intact. If you modify this file, please                 !
!  note section 2a) of the GPLv2 states that:                                  !
!                                                                              !
!  a) You must cause the modified files to carry prominent notices             !
!     stating that you changed the files and the date of any change.           !
!                                                                              !
!  ChangeLog:                                                                  !
!------------------------------------------------------------------------------!

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
 
 use bound
 use eos
 use options
 use part
 use setup_params
 use derivb
 use rates
 use hterms, only:gradh
!
!--define local variables
!
 implicit none
 real, intent(in) :: t
 character(len=*), intent(in) :: dumpfile
 integer :: i,nprint,ncolumns
 integer :: ierr,iformat

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
 ncolumns = ndim + 2*ndimV + 4
 if (imhd.ne.0) then
    ncolumns = ncolumns + 8 + 2*ndimV ! number of columns
    iformat = 2
    if (imhd.lt.0) ncolumns = ncolumns + ndimV
 else
    ncolumns = ncolumns + 5
    iformat = 1
    if (igravity.ne.0) ncolumns = ncolumns + 1
    if (allocated(del2v)) ncolumns = ncolumns + 1
 endif
 if (geom(1:4).ne.'cart') then
    ncolumns = ncolumns + 2 + ndimV
    iformat = iformat + 2
 endif
 if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
    iformat = 5
    ncolumns = ncolumns + ndimV + 1
 endif
 
 write(idatfile,iostat=ierr) t,npart,nprint,gamma,hfact,ndim,ndimV, &
      ncolumns,iformat,ibound,xmin(1:ndim),xmax(1:ndim),len(geom),geom
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
  !--essential variables
  do i=1,ndim
     write(idatfile) x(i,1:nprint)
  enddo
  do i=1,ndimV
     write(idatfile) vel(i,1:nprint)
  enddo
  write(idatfile) hh(1:nprint)
  if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
     write(idatfile) rho(1:nprint)  
  else
     write(idatfile) dens(1:nprint)
  endif
  write(idatfile) uu(1:nprint)
  write(idatfile) pmass(1:nprint)
  
  if (imhd.ne.0) then
     do i=1,3
        write(idatfile) alpha(i,1:nprint)
     enddo
     do i=1,ndimV
        write(idatfile) Bfield(i,1:nprint)
     enddo
     write(idatfile) psi(1:nprint)
     !--info only
     write(idatfile) pr(1:nprint)
     write(idatfile) -drhodt(1:nprint)/rho(1:nprint)
     write(idatfile) divB(1:nprint)
     do i=1,ndimV
        write(idatfile) curlB(i,1:nprint)
     enddo
     write(idatfile) gradh(1:nprint)
     do i=1,ndimV
        write(idatfile) force(i,1:nprint)
     enddo
     if (imhd.lt.0) then
        do i=1,ndimV
           write(idatfile) Bevol(i,1:nprint)
        enddo
        write(idatfile) Bconst(:)
     endif
  else
     do i=1,2
        write(idatfile) alpha(i,1:nprint)
     enddo
     !--info only
     write(idatfile) pr(1:nprint)
     write(idatfile) -drhodt(1:nprint)/rho(1:nprint)
     write(idatfile) gradh(1:nprint)
     do i=1,ndimV
        write(idatfile) force(i,1:nprint)
     enddo
     if (igravity.ne.0) write(idatfile) poten(1:nprint)
     if (allocated(del2v)) write(idatfile) del2v(1:nprint)
  endif
  if (geom(1:4).ne.'cart') then
     write(idatfile) rho(1:nprint)
     write(idatfile) sqrtg(1:nprint)
     do i=1,ndimV
        write(idatfile) pmom(i,1:nprint)
     enddo
  endif
  if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
     write(idatfile) dustfrac(1:nprint)
     do i=1,ndimV
        write(idatfile) deltav(i,1:nprint)
     enddo
  endif
  write(idatfile) itype(1:nprint)

 close(unit=idatfile)

 return
end subroutine write_dump

!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  setup by reading a dump from a file                                   !!
!!  should be compatible with the output format given in write_dump       !!
!!                                                                        !!
!!------------------------------------------------------------------------!!

subroutine read_dump(dumpfile,tfile,copysetup)
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
 use mem_allocation, only:alloc
!
!--define local variables
!      
 implicit none
 character(len=*), intent(in) :: dumpfile
 real, intent(out) :: tfile
 logical, optional, intent(in) :: copysetup
 integer :: i,ndimfile,ndimvfile,npartfile,nprintfile,ncolumns
 integer :: ierr, iformat,lengeom,nread
 character(len=len(geom)) :: geomfile
 real :: gammafile,hfactfile
 logical :: iexist

 geomfile = 'cartes'
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
                          ndimfile,ndimvfile,ncolumns,iformat, &
                          ibound(1:ndimfile),xmin(1:ndimfile),xmax(1:ndimfile),lengeom,geomfile(1:lengeom)

 if (ierr /= 0) then
    write(iprint,*) 'error reading header from dump file: ',trim(dumpfile)
    stop 
 endif
 write(iprint,*) 'time = ',tfile,' in dump file '
 write(iprint,*) 'xmin = ',xmin(1:ndimfile),' xmax = ',xmax(1:ndimfile)
 write(iprint,*) 'ibound = ',ibound(1:ndimfile),' npart in file = ',npartfile
!
!--copy setup options from header if desired
!
 if (present(copysetup)) then
    if (copysetup) then
       gamma = gammafile
       hfact = hfactfile
    endif
 endif
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

 if (imhd.eq.0 .and. (iformat.eq.2 .or. iformat.eq.4)) then
    write(iprint,*) 'warning: mhd input file, but MHD is off'
 elseif (imhd.ne.0 .and. (iformat.ne.2 .and. iformat.ne.4)) then
    write(iprint,*) 'WARNING: non-mhd infile but MHD is on (Bfield set to 0)'
 elseif (imhd.lt.0 .and. (ncolumns.lt.(ndim + 5*ndimV + 12) .or. iformat.ne.2)) then
    write(iprint,*) 'ERROR: cannot re-start with vector potential from this file'
    stop
 elseif ((idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) .and. (iformat.ne.5)) then
    write(iprint,*) 'ERROR: idust=1 but dump file does not contain dustfrac or deltav arrays'
    !stop
 endif
!
!--switch current geometry to that of the file if not convertible
!  (need to know whether we have a non-cartesian geometry before
!   memory allocation)
!
 geomsetup = geomfile
 if (geomsetup(1:6).ne.geom(1:6)) then
    select case(geomsetup(1:6))
       case('sphrpt','cylrpz','cartes')
          write(iprint,*) 'file geometry = ',geomfile
       case('cylrzp')
          if (ndim.eq.2) then
             write(iprint,*) '=> Using original co-ordinate system'
             geom = geomsetup
          else
             write(iprint,*) 'file geometry = ',geomfile          
          endif
       case default
         write(iprint,*) '=> Using original co-ordinate system'
         geom = geomsetup
    end select
 endif
!
!--allocate memory for the number of particles
!
 npart = npartfile
 ntotal = npartfile

 call alloc(nprintfile)
!
!--zero quantities which may not be explicitly read
!
 Bfield = 0.
 alpha = 0.
!
!--read data from file (only bits needed to restart the run - do not read ghosts)
!
 do i=1,ndim
    read(ireadf,iostat=ierr) x(i,1:npart)
 enddo
 if (ierr /= 0) then
    write(iprint,*) '*** error reading data from dumpfile ***'
    stop
 endif
 do i=1,ndimV
    read(ireadf,iostat=ierr) vel(i,1:npart)
 enddo
 if (ierr /= 0) stop 'error reading dumpfile'
 read(ireadf,iostat=ierr) hh(1:npart)
 if (ierr /= 0) stop 'error reading dumpfile'
 read(ireadf,iostat=ierr) dens(1:npart)
 if (ierr /= 0) stop 'error reading dumpfile'
 read(ireadf,iostat=ierr) uu(1:npart)
 if (ierr /= 0) stop 'error reading dumpfile'
 read(ireadf,iostat=ierr) pmass(1:npart)
 if (ierr /= 0) stop 'error reading dumpfile'
 nread = ndim + ndimV + 4
 
 if (iformat.eq.2 .or. iformat.eq.4 .and. imhd.ne.0) then
    do i=1,3
       read(ireadf,iostat=ierr) alpha(i,1:npart)
    enddo
    if (ierr /= 0) then
       write(iprint,*) '*** error in dumpfile : non-MHD file ***'
    endif
    do i=1,ndimV
       read(ireadf,iostat=ierr) Bfield(i,1:npart)
    enddo
    read(ireadf,iostat=ierr) psi(1:npart)
    nread = nread + 4 + ndimV
    
    !--read vector/euler potentials if required
    if (imhd.lt.0) then
       !--skip other quantities
       do i=1,4+2*ndimV
          read(ireadf,iostat=ierr)    
          if (ierr /= 0) stop 'readdump: error skipping columns for vector potential'
       enddo
       do i=1,ndimV
          read(ireadf,iostat=ierr) Bevol(i,1:npart)
          if (ierr /= 0) stop 'readdump: error reading vector potential'
       enddo
       read(ireadf,iostat=ierr) Bconst(1:ndimV)
       if (ierr /= 0) then
          write(iprint,*) '*** error reading Bconst from file ***'
          write(iprint,*) 'please enter values for Bconst:'
          read*,Bconst(1:ndimV)
       endif
       nread = nread + 5 + 3*ndimV
    endif
 elseif (iformat.eq.2) then
    !--read alpha in MHD format
    do i=1,3
       read(ireadf,iostat=ierr) alpha(i,1:npart)    
    enddo
    nread = nread + 3
 else
    do i=1,2
       read(ireadf,iostat=ierr) alpha(i,1:npart)
    enddo
    nread = nread + 2
 endif
 if ((idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) .and. iformat.eq.5) then
    read(ireadf,iostat=ierr) dustfrac(1:npart)
    do i=1,ndimV
       read(ireadf,iostat=ierr) deltav(i,1:npart)
    enddo
 endif
 !--skip extra quantities that come after alpha but before itype
 do i=1,ncolumns-nread
    read(ireadf,iostat=ierr)
 enddo

 read(ireadf,iostat=ierr) itype(1:npart)
 if (ierr /= 0) then
    write(iprint,*) 'WARNING: itype array not present in file'
 endif
!
!--close the file
!
 close(unit=ireadf)
 
 write(iprint,*) 'finished reading setup file: everything is aok'

 return

end subroutine read_dump

end module dumpfiles

