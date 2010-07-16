!
!--Wrapper for output of dumps : sets file name and writes to log
!
subroutine output(t,nstep)
 use debug
 use loguns
 use part, only:npart,ntotal
 use options, only:imhd,ibound
 
 use dumpfiles
 use cons2prim, only:conservative2primitive
 implicit none
 integer, intent(in) :: nstep
 real, intent(in) :: t
 character(len=len(rootname)+10) :: dumpfile
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' entering subroutine output' 
!
!--calculate primitive variables from conservatives for output
!  nstep <0 means do not do this as we are on a quit dump
!
 if (nstep.ge.0) then
    if (imhd.lt.0) then
       if (any(ibound.ge.2)) call set_ghost_particles
       call set_linklist
    endif
    call conservative2primitive  ! also calls equation of state
 endif
!
!--create new dumpfile if rootname contains numbers
!
 ifile = ifile + 1
 write(dumpfile,"(a,'_',i5.5,'.dat')") trim(rootname),ifile

!
!--write timestep info to log file
!      
!! write(iprint,"('| ',a73,'|')") trim(dumpfile)
 write (iprint,10) t,abs(nstep),npart,ntotal-npart,trim(dumpfile)
10 format('| t = ',f7.3,' | step = ',i8,' | npart = ',i7,        &
             ' | nghost = ',i5,' |',a)
!! write(iprint,"(76('-'))")

 call write_dump(t,dumpfile)

end subroutine output
