!-----------------------------------------------
! reads default options from file
! uses namelist input to group the options
! these are specified in the modules
!-----------------------------------------------
subroutine defaults_read
 use filenames
 use multiplot
 use settings_data
 use settings_part
 use settings_page
 use settings_render
 use settings_vecplot
 use settings_xsecrot
 use settings_powerspec
 use exact
 implicit none
 logical :: iexist
 integer :: ierr,i
 
 inquire (exist=iexist, file='defaults')
 if (iexist) then
    open(unit=1,file='defaults',status='old',form='formatted')
    
    ierr = 0
    read(1,NML=dataopts,end=77,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading data options from defaults'    
    
    ierr = 0
    read(1,NML=plotopts,end=77,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading plot options from defaults'

    ierr = 0
    read(1,NML=pageopts,end=77,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading page options from defaults'

    ierr = 0
    read(1,NML=renderopts,end=77,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading render options from defaults'

    ierr = 0
    read(1,NML=vectoropts,end=77,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading vector plot options from defaults'

    ierr = 0
    read(1,NML=xsecrotopts,end=77,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading xsec/rotation options from defaults'

    ierr = 0
    read(1,NML=powerspecopts,end=77,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading power spectrum options from defaults'

    ierr = 0
    read(1,NML=exactparams,end=77,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading exact solution parameters from defaults'    
  
    ierr = 0
    read(1,NML=multi,end=77,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading multiplot options from defaults'

    do i=1,maxfile
       read(1,*,end=66,iostat=ierr) rootname(i)
    enddo
66  continue

    close(unit=1)
    print*,'read default options from file '
    return
 else
    print*,'defaults file not found: using program settings'
    return
 endif
 
77 continue
 print*,'**** warning: end of file in defaults ****'
 close(unit=1)

 return
end subroutine defaults_read
