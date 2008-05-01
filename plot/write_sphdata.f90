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
 case('ascii')
     issphformat = .true.
 end select
 
 if (.not.issphformat) then
    print "(a)",' possible formats for convert mode ("splash to X"): '
    print "(a)",' splash to ascii : convert SPH data to ascii file'
 endif
 
 return
end function issphformat

subroutine write_sphdump(time,dat,npart,ntypes,npartoftype,itype,ncolumns,filename,outformat)
 use labels, only:labeltype,label
 use settings_units, only:unitslabel,units
 use params, only:int1
 implicit none
 integer, intent(in) :: npart,ntypes,ncolumns
 integer, intent(in), dimension(:) :: npartoftype
 integer(kind=int1), intent(in), dimension(:) :: itype
 real, intent(in) :: time
 real, intent(in), dimension(npart,ncolumns) :: dat
 character(len=*), intent(in) :: filename,outformat
 integer, parameter :: iunit = 83
 integer :: ierr,i
 character(len=40) :: fmtstring,fmtstring2,fmtstringlab
 
 write(fmtstring,"(i10,a)") ncolumns,'(1pe15.7,1x)'
 fmtstring2 = '('//trim(adjustl(fmtstring))//',i1)'
 fmtstring = '('//trim(adjustl(fmtstring))//')'

 write(fmtstringlab,"(i10,a)") ncolumns,'(a15,1x),a'
 fmtstringlab = '(''#'',1x,'//trim(adjustl(fmtstringlab))//')'

 select case(trim(outformat))
 case ('ascii')
    print "(/,5('-'),'>',a,i2,a,1x,'<',5('-'),/)",' WRITING TO FILE '//trim(filename)//'.ascii WITH ',ncolumns,' COLUMNS'
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
 case default
    print "(a)",' ERROR: unknown output format '''//trim(outformat)//''' in write_sphdump'
    return
 end select
 
end subroutine write_sphdump

end module write_sphdata
