module ppmutils
 implicit none

contains

subroutine get_colours_from_ppm(rgbtable,ierr)
 implicit none
 integer, intent(out) :: ierr
 real, dimension(:,:), intent(out) :: rgbtable
 character(len=100) :: filename
 character(len=2) :: ptype
 integer :: npixx,npixy,ncolours,i
 
 write(*,"(1x,a)",ADVANCE='NO') 'enter ppm file to read colour table from: '
 read*,filename
 open(unit=1,file=filename,status='old',form='formatted',iostat=ierr)
 if (ierr /= 0) then
    print*,'cannot open ',trim(filename)
    return
 endif
 read(1,*,iostat=ierr) ptype
 if (ierr /= 0) then
    print*,'error: cannot read format type in ',trim(filename)
    return
 elseif (ptype(1:1).ne.'P') then
    print*,trim(filename),' is not pnm format -> convert first using pbm tools'
    ierr = 1
    return
 elseif (ptype.ne.'P3') then
    print*,trim(filename),' is not plain pnm -> convert first using pnmtoplainpnm'
    ierr = 2
    return
 endif
 print*,'-> reading ',trim(filename)
 read(1,*) npixx,npixy,ncolours
 print*,'ppm image size is ',npixx,'x',npixy
 ncolours = ncolours + 1
 print*,'ppm image contains ',ncolours,' colours'
 if (ncolours.gt.size(rgbtable(1,:))) then
    print*,'error: too many colours in ppm -> use pnmquant to get 256 colours only'
    ierr = 3
    return
 endif

 read(1,*,iostat=ierr) (rgbtable(:,i),i=1,ncolours+1)
 if (ierr /= 0) then
    print*,'error encountered whilst reading colour table'
    return
 else
    print*,'colour table successfully read '
    do i=1,ncolours
       rgbtable(:,i) = rgbtable(:,i)/real(ncolours-1)
       print*,'r,g,b = ',rgbtable(:,i)
    enddo
 endif

 close(1)
 return
end subroutine get_colours_from_ppm

end module ppmutils
