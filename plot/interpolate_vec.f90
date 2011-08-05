!--------------------------------------------------------------------------
!
!  Module with some auxiliary routines related to vector interpolation
!
!--------------------------------------------------------------------------
module interpolate_vec
 implicit none
 public :: interpolate_vec_average
 public :: mask_vectors
 
 private
 
contains

!--------------------------------------------------------------------------
!
!  Hides vector arrows (sets them to zero) where there are no particles
!  contained in the pixel (as opposed to merely contributing to the pixel)
!
!  means you can avoid funny looking plots with arrows in otherwise
!  empty regions
!
!  Daniel Price 26/3/07
!
!--------------------------------------------------------------------------
subroutine mask_vectors(xplot,yplot,itype,npart,xmin,xmax,ymin,ymax, &
                        vecpixx,vecpixy,npixvecx,npixvecy,minincell)
 implicit none
 integer, intent(in) :: npart,npixvecx,npixvecy,minincell
 integer, dimension(npart), intent(in) :: itype
 real, dimension(npart), intent(in) :: xplot,yplot
 real, intent(in) :: xmin,xmax,ymin,ymax
 real, dimension(npixvecx,npixvecy), intent(inout) :: vecpixx,vecpixy
 integer, dimension(npixvecx,npixvecy) :: nincell
 integer :: icellx,icelly,j
 real :: dxcell1,dycell1
 character(len=16) :: chmin
 
 !--write nice, neat information line
 write(chmin,"(g10.0)") minincell
 print "(2x,a)",'(hiding arrows where there are < '//trim(adjustl(chmin))//' particles in pixel cell)'

 dxcell1 = (npixvecx - 1)/(xmax-xmin + tiny(xmin))
 dycell1 = (npixvecy - 1)/(ymax-ymin + tiny(ymin))
!
!--count particles which fall into each pixel ("cell")
!
 nincell(:,:) = 0
 do j=1,npart
    if (itype(j).ge.0) then ! exclude not-plotted particles
       icellx = int((xplot(j) - xmin)*dxcell1) + 1
       icelly = int((yplot(j) - ymin)*dycell1) + 1
       !--count number of particles in each cell
       if (icellx.gt.0 .and. icellx.le.npixvecx &
          .and. icelly.gt.0 .and. icelly.le.npixvecy) then
          nincell(icellx,icelly) = nincell(icellx,icelly) + 1
       endif
    endif
 enddo 
!
!--set vector arrow lengths to zero in cells where there are no particles
!
 where (nincell.lt.minincell)
    vecpixx = 0.
    vecpixy = 0. 
 end where
 
 return
end subroutine mask_vectors

!--------------------------------------------------------------------------
!    Interpolates vector quantity from particles to even grid of pixels
!
!    This version just does a simple averaging by binning particles
!    and taking the average of vx,vy in the cell to give a vector for
!    that cell. This is because the interpolation of a vector quantity is
!    usually to a *coarser* grid than the particles.
!
!    Input: particle coordinates  : x,y   (npart)
!           vector data to smooth : vecx  (npart)
!                                   vecy  (npart)
!           grid setup : xmin, ymin, dx
!
!     Output: smoothed vector field   : vecpixx (npixx,npixy)
!                                     : vecpixy (npixx,npixy)
!
!     Daniel Price, Institute of Astronomy, Cambridge, 20/8/04
!--------------------------------------------------------------------------

subroutine interpolate_vec_average(x,y,vecx,vecy,itype, &
     xmin,ymin,dx,vecpixx,vecpixy,npart,npixx,npixy)

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,vecx,vecy
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,dx
  real, intent(out), dimension(npixx,npixy) :: vecpixx, vecpixy
  integer :: i,j,k,ix,iy
  integer, dimension(npixx,npixy) :: ihoc,numcell
  integer, dimension(npart) :: ll

  print*,'averaging vector field onto pixels...'
  if (dx.le.0.) then
     print*,'interpolate_vec: error: pixel width <= 0'
     return
  endif
  !
  !--interpolation is to a coarser grid, so just average
  !  bin particles into cells using a link list 
  !     
  ihoc(:,:) = -1   ! head of chain
  numcell(:,:) = 0
  do i=1,npart
     if (itype(i).ge.0) then
        ix = int((x(i)-xmin)/dx)+1
        iy = int((y(i)-ymin)/dx)+1
        if ((ix.ge.1).and.(ix.le.npixx).and.(iy.ge.1).and.(iy.le.npixy)) then
           ll(i)=ihoc(ix,iy)   ! set link list of this particle to old head of list
           ihoc(ix,iy) = i            ! set head of chain to this particle
        endif
     endif
  enddo
  !
  !--add up total vx,vy in each cell
  !
  vecpixx(:,:) = 0.
  vecpixy(:,:) = 0.
  do j=1,npixy
     do i=1,npixx
        k = ihoc(i,j)
        do while (k.ne.-1)
           vecpixx(i,j) = vecpixx(i,j) + vecx(k)
           vecpixy(i,j) = vecpixy(i,j) + vecy(k)
           numcell(i,j) = numcell(i,j) + 1
           k = ll(k)
        enddo
     enddo
  enddo
  !
  !--divide by number of particles in that cell to get average vx,vy
  !
  do j=1,npixy
     do i=1,npixx
        if (numcell(i,j).ne.0) then
           vecpixx(i,j) = vecpixx(i,j)/float(numcell(i,j))
           vecpixy(i,j) = vecpixy(i,j)/float(numcell(i,j))         
        endif
     enddo
  enddo
  return
  
end subroutine interpolate_vec_average

end module interpolate_vec
