!----------------------------------------------------------------
!
! module for probability distribution function calculation
! and/or plotting on the particles
!
!----------------------------------------------------------------
module pdfs
 implicit none
 public :: pdfcalc,pdfplot,write_pdf

contains

!-----------------------------------------------------------------
!
! subroutine bins particles into x, works out number in each bin,
! calculates normalisation for PDF.
!
!-----------------------------------------------------------------
subroutine pdfcalc(npart,xpart,xminplot,xmaxplot,nbins,xbin,pdf,pdfmin,pdfmax,itransx,itransy,labelx,icolours,rhopart)
 use transforms, only:transform,transform_inverse,transform_limits
 implicit none
 integer, intent(in) :: npart,nbins
 real, dimension(:), intent(in) :: xpart
 real, intent(in) :: xminplot,xmaxplot
 real, intent(out), dimension(nbins) :: xbin,pdf
 real, intent(out) :: pdfmin,pdfmax
 integer, intent(in) :: itransx,itransy
 character(len=*), intent(in) :: labelx
 integer, intent(in), dimension(:) :: icolours
 real, intent(in), dimension(:), optional :: rhopart
 integer :: ibin,i
 real :: dx,totprob,fi,fprev,xbini,xbinprev,dxprev
 real :: xmin,xmax
 logical :: volweighted
 
 if (present(rhopart)) then
    print "(a,i3,a)",' calculating (volume weighted) PDF using ',nbins,' bins'
    volweighted = .true.
 else
    print "(a,i3,a)",' calculating (density weighted) PDF using ',nbins,' bins'
 endif
 !
 !--set bins in PDF: must always use all the particles
 !  note that xpart will already have been transformed
 !  so these are min and max in transformed space
 !
 xmin = minval(xpart(1:npart))
 xmax = maxval(xpart(1:npart))
 
 dx = (xmax - xmin)/real(nbins - 1)
 do ibin=1,nbins
    xbin(ibin) = xmin + (ibin-0.5)*dx
 enddo
!
!--now calculate probability of finding a particle at each x
!
 pdf(:) = 0.
 do i=1,npart
    !--do not use hidden particles
    if (icolours(i).ge.0) then
       ibin = int((xpart(i) - xmin)/dx) + 1
       if (ibin.lt.1) ibin = 1
       if (ibin.gt.nbins) ibin = nbins

       if (volweighted) then
          pdf(ibin) = pdf(ibin) + 1./rhopart(i)
       else
          pdf(ibin) = pdf(ibin) + 1.
       endif
    endif
 enddo
!
!--get total area under pdf by trapezoidal rule
!
 totprob = 0.
 fprev = 0.
 xbinprev = xmin
 if (itransx.gt.0) call transform_inverse(xbinprev,itransx)
 dxprev = 0.
 
 do ibin=1,nbins
    fi = pdf(ibin)
    if (itransx.gt.0) then
       xbini = xmin + ibin*dx
       call transform_inverse(xbini,itransx)
    ! == not used == \int pdf dx = ln(10)*\int x*pdf d(log_10 x)
    ! instead just use \int pdf dx with dx varying
       totprob = totprob + 0.5*((xbini - xbinprev)*fi + dxprev*fprev)
    else
       totprob = totprob + 0.5*dx*(fi + fprev)
    endif
    dxprev = xbini - xbinprev
    fprev = fi
 enddo
!
!--normalise pdf so total area is unity
!
 print*,'normalisation factor = ',totprob ! =npart*dx for equispaced
 pdf(1:nbins) = pdf(1:nbins)/totprob
 
 call write_pdf(nbins,xbin,pdf,labelx,itransx,volweighted)
!
!--return min and max for adaptive plot limit setting
!  (exclude zero as min)
! 
 pdfmin = minval(pdf(1:nbins),mask=(pdf(1:nbins).gt.0.))
 pdfmax = maxval(pdf(1:nbins))
!
!--apply transformations to y data
!
 if (itransy.gt.0) then
    call transform(pdf,itransy)
    call transform_limits(pdfmin,pdfmax,itransy)
 endif
 
end subroutine pdfcalc

!-----------------------------------------------------------------
! interface which controls plotting of PDF
! (so can easily change properties of PDF plotting,
!  e.g. histogram vs. line)
!-----------------------------------------------------------------
subroutine pdfplot(nbins,xbin,pb)
 use plotutils, only:plotline,plotbins
 implicit none
 integer, intent(in) :: nbins
 real, dimension(:), intent(in) :: xbin,pb

!
!--plot as line segment, with blanking at zero
! 
 call plotline(nbins,xbin,pb,blank=0.)
!
!--plot as histogram, with blanking of zero
!
! call plotbins(nbins,xbin,pb,blank=0.)

end subroutine pdfplot

!-----------------------------------------------------------------
! routine to write pdf to file
!-----------------------------------------------------------------
subroutine write_pdf(nbins,xbin,pb,labelx,itransx,volweighted)
 use filenames, only:rootname,ifileopen
 use transforms, only:transform_label,transform_inverse
 implicit none
 character(len=*), intent(in) :: labelx
 integer, intent(in) :: nbins,itransx
 real, intent(in), dimension(nbins) :: xbin,pb
 logical, intent(in) :: volweighted
 real, dimension(nbins) :: xbintemp
 integer :: i,ierr
 integer, parameter :: iunit = 86
 logical :: warned
 
 print "(a)",' writing to '//trim(rootname(ifileopen))//'_pdf_'//trim(labelx)//'.dat'
 open(unit=iunit,file=trim(rootname(ifileopen))//'_pdf_'//trim(labelx)//'.dat', &
      form='formatted',status='replace',iostat=ierr)
 if (ierr /= 0) then
    print "(a)",'ERROR: could not open file: no output'
 endif
 
 if (volweighted) then
    write(iunit,"(a)",iostat=ierr) '# volume weighted PDF: calculated using SPLASH (c)2008 Daniel Price '
 else
    write(iunit,"(a)",iostat=ierr) '# density weighted PDF: calculated using SPLASH (c)2008 Daniel Price '
 endif
 if (ierr /= 0) print "(a)",' ERROR writing header line'
 write(iunit,"(a,i5,a)",iostat=ierr) '# ',nbins,' bins evenly spaced in '//trim(transform_label(labelx,itransx))

 warned = .false.
 !--output x bins in un-transformed space
 xbintemp = xbin
 if (itransx.gt.0) call transform_inverse(xbintemp,itransx)
 !--dump bins to file
 do i=1,nbins
    write(iunit,*,iostat=ierr) xbintemp(i),pb(i)
    if (ierr /= 0 .and. .not.warned) then
       print "(a)",' ERRORS during write'
       warned = .true.
    endif
 enddo
 close(iunit)
 
 return
end subroutine write_pdf

end module pdfs
