!-----------------------------------------------------------------
!     module containing output routines for writing SPH data
!     as read by SPLASH to output file in various formats
!
!     (c) D. Price 22/01/08
!-----------------------------------------------------------------
module write_sphdata
 public :: issphformat,write_sphdump
 private

contains

!-----------------------------------------------------------------
! utility to check if a format selection is valid
!-----------------------------------------------------------------
logical function issphformat(string)
 implicit none
 character(len=*), intent(in) :: string

 issphformat = .false.
 select case(trim(string))
 case('ascii','ASCII')
     issphformat = .true.
 case('binary','BINARY')
     issphformat = .true.
 case('rsph','RSPH')
     issphformat = .true.
 case('phantom','PHANTOM')
     issphformat = .true.
 end select
 
 if (.not.issphformat) then
    print "(a)",' convert mode ("splash to X dumpfiles"): '
    print "(a,/)",' splash to ascii   : convert SPH data to ascii file dumpfile.ascii'
    print "(a)",  '        to binary  : convert SPH data to simple unformatted binary dumpfile.binary '
    print "(a)",  '                      write(1) time,npart,ncolumns'
    print "(a)",  '                      do i=1,npart'
    print "(a)",  '                         write(1) dat(1:ncolumns),itype'
    print "(a)",  '                      enddo'
    print "(a)",  '        to phantom : convert SPH data to binary dump file for PHANTOM'
 endif
 
 return
end function issphformat

subroutine write_sphdump(time,gamma,dat,npart,ntypes,npartoftype,masstype,itype,ncolumns,filename,outformat)
 use labels,         only:labeltype,label,irho,ipmass,ix
 use settings_units, only:unitslabel,units
 use settings_data,  only:ndim
 use params,         only:int1
 use write_data_phantom, only:write_sphdata_phantom
 implicit none
 integer, intent(in)                          :: npart,ntypes,ncolumns
 integer, intent(in), dimension(:)            :: npartoftype
 integer(kind=int1), intent(in), dimension(:) :: itype
 real, intent(in)                             :: time,gamma
 real, intent(in), dimension(npart,ncolumns)  :: dat
 real, intent(in), dimension(:)               :: masstype
 character(len=*), intent(in)                 :: filename,outformat
 integer, parameter :: iunit = 83
 integer, parameter :: maxline = 1000
 integer            :: ierr,i,idim,i1,i2
 character(len=40)  :: fmtstring,fmtstring2,fmtstringlab,outfile

 select case(trim(outformat))
 case ('ascii','ASCII')
    print "(/,5('-'),'>',a,i2,a,1x,'<',5('-'),/)",' WRITING TO FILE '//trim(filename)//'.ascii WITH ',ncolumns,' COLUMNS'

    !--format the header lines to go in the ascii file 
    write(fmtstring,"(i10,a)") ncolumns,'(1pe15.7,1x)'
    fmtstring2 = '('//trim(adjustl(fmtstring))//',i1)'
    fmtstring = '('//trim(adjustl(fmtstring))//')'

    write(fmtstringlab,"(i10,a)") ncolumns,'(a15,1x),a'
    fmtstringlab = '(''#'',1x,'//trim(adjustl(fmtstringlab))//')'

    open(unit=iunit,file=trim(filename)//'.ascii',status='replace',form='formatted',iostat=ierr)
       if (ierr /= 0) then
          print "(a)",' ERROR OPENING FILE FOR WRITING'
          return
       endif
       write(iunit,"(a)",iostat=ierr) '# '//trim(filename)// &
                          '.ascii created by SPLASH, a visualisation tool for SPH data (c) 2008 Daniel Price'
       write(iunit,"('#')")
       write(iunit,"('#',1x,'time:',13x,'time unit (',a,')')",iostat=ierr) trim(unitslabel(0))
       write(iunit,"('#',2(1x,1pe15.7))",iostat=ierr) time,units(0)
       write(iunit,"('#')")
       write(iunit,"('#',1x,'npart:',6(1x,a12))",iostat=ierr) (trim(labeltype(i)),i=1,ntypes)
       write(iunit,"('#',7x,6(1x,i12))",iostat=ierr) npartoftype(1:ntypes)
       write(iunit,"('# units:')")
       write(iunit,"('#'"//fmtstring(2:),iostat=ierr) units(1:ncolumns)
       write(iunit,fmtstringlab,iostat=ierr) unitslabel(1:ncolumns)
       write(iunit,"('#')")
       !
       !--write body
       !
       if (size(itype).gt.1) then
          write(iunit,fmtstringlab,iostat=ierr) label(1:ncolumns),'itype'
          do i=1,npart
             write(iunit,fmtstring2,err=100) dat(i,1:ncolumns),itype(i)
          enddo
       else
          write(iunit,fmtstringlab,iostat=ierr) label(1:ncolumns)
          do i=1,npart
             write(iunit,fmtstring,err=100) dat(i,1:ncolumns)
          enddo       
       endif
    close(iunit)

    return
100 continue
    close(iunit)
    print*,'ERROR WRITING ASCII FILE'
    return

 case ('binary','BINARY')
!
!--This is the most basic binary (ie. unformatted) file format I could think of,
!  as an alternative to ascii for large files.
!
    print "(/,5('-'),'>',a,i2,a,1x,'<',5('-'),/)",' WRITING TO FILE '//trim(filename)//'.binary WITH ',ncolumns,' COLUMNS'
    open(unit=iunit,file=trim(filename)//'.binary',status='replace',form='unformatted',iostat=ierr)
       if (ierr /= 0) then
          print "(a)",' ERROR OPENING FILE FOR WRITING'
          return
       endif
       write(iunit,iostat=ierr) time,npart,ncolumns
       if (ierr /= 0) then
          print "(a)",' ERROR WRITING HEADER LINE TO BINARY FILE '
       endif
       !
       !--write body
       !
       if (size(itype).gt.1) then
          do i=1,npart
             write(iunit,err=200) dat(i,1:ncolumns),int(itype(i))
          enddo
       else
          do i=1,npart
             write(iunit,err=200) dat(i,1:ncolumns)
          enddo       
       endif
    close(iunit)

    return
200 continue
    close(iunit)
    print*,'ERROR WRITING BINARY FILE'
    return
 case ('rsph','RSPH')
!
!--Files for Steinar Borve's RSPH format
!
    if (ndim.lt.2) then
       print "(a)",' ERROR: cannot write RSPH format for < 2D'
       return
    endif
    outfile = 'rsph2D_pos.dat'
    print "(a)",' writing to '//trim(outfile)
    open(unit=iunit,file=outfile,status='replace',form='formatted',iostat=ierr)
    if (ierr /= 0) then
       print "(a)",' ERROR OPENING '//trim(outfile)//' FOR WRITING'
       return
    endif
    write(iunit,"(i1)") ndim
    write(iunit,"(a)") 'position'
    write(iunit,"(i4)") maxline
    write(iunit,*) (minval(dat(1:npart,ix(i))),i=1,ndim)
    write(iunit,*) (maxval(dat(1:npart,ix(i))),i=1,ndim)
    write(iunit,*) time
    write(iunit,*) npart
    do idim=1,2
       i1 = 1
       i2 = 0
       do while (i2 < npart)
          i2 = min(i2 + maxline,npart)
          write(iunit,*) dat(i1:i2,ix(idim))
          i1 = i2 + 1
       enddo
    enddo
    close(unit=iunit)

    outfile = 'rsph2D_rho.dat'
    print "(a)",' writing to '//trim(outfile)
    open(unit=iunit,file=outfile,status='replace',form='formatted',iostat=ierr)
    if (ierr /= 0) then
       print "(a)",' ERROR OPENING FILE FOR WRITING'
       return
    endif
    write(iunit,"(i1)") ndim
    write(iunit,"(a)") 'density'
    write(iunit,"(i4)") maxline
    write(iunit,*) (minval(dat(1:npart,ix(i))),i=1,ndim)
    write(iunit,*) (maxval(dat(1:npart,ix(i))),i=1,ndim)
    i1 = 1
    i2 = 0
    do while (i2 < npart)
       i2 =  min(i2 + maxline,npart)
       write(iunit,*) dat(i1:i2,irho)
       i1 = i2 + 1
    enddo
    close(unit=iunit)

    outfile = 'rsph2D_siz.dat'
    print "(a)",' writing to '//trim(outfile)
    open(unit=iunit,file=outfile,status='replace',form='formatted',iostat=ierr)
    if (ierr /= 0) then
       print "(a)",' ERROR OPENING FILE FOR WRITING'
       return
    endif
    write(iunit,"(i1)") ndim
    write(iunit,"(a)") 'size'
    write(iunit,"(i4)") maxline
    write(iunit,*) (minval(dat(1:npart,ix(i))),i=1,ndim)
    write(iunit,*) (maxval(dat(1:npart,ix(i))),i=1,ndim)
    i1 = 1
    i2 = 0
    do while (i2 < npart)
       i2 =  min(i2 + maxline,npart)
       write(iunit,*) ((dat(i,ipmass)/dat(i,irho))**(1./ndim),i=i1,i2)
       i1 = i2 + 1
    enddo
    close(unit=iunit)

 case('phantom','PHANTOM')  
     if (size(itype).gt.1 .and. sum(npartoftype(2:)).gt.0) then
        print "(a)",' ERROR: writing of PHANTOM dumps not implemented with mixed types'
     else
        call write_sphdata_phantom(time,gamma,dat,npart,ntypes,npartoftype,&
                                   masstype,ncolumns,filename)
     endif
 case default
    print "(a)",' ERROR: unknown output format '''//trim(outformat)//''' in write_sphdump'
    return
 end select
 
end subroutine write_sphdump

end module write_sphdata
