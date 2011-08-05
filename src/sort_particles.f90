subroutine sort_particles
 use loguns, only:iprint
 use part, only:rho,npart
 use mem_allocation, only:alloc
 use linklist, only:ncellsloop,ifirstincell,ll
 use get_neighbour_lists, only:get_neighbour_list
 implicit none
 integer, dimension(size(rho)) :: ilist,listneigh
 logical, dimension(size(rho)) :: sorted
 integer :: i,icell,ipart,nneigh,n,j,idone
 real :: t1,t2
 
 write(iprint,*) 'sorting particles...'
 call cpu_time(t1)
 
 ipart = 0
 sorted = .false.
 do icell=1,ncellsloop
    call get_neighbour_list(icell,listneigh,nneigh)
 
    i = ifirstincell(icell)
    idone = -1
    do while (i.ne.-1 .and. ipart.lt.npart)
       idone = idone + 1
       if (i.le.npart .and. .not.sorted(i)) then
          ipart = ipart + 1
          ilist(ipart) = i
          sorted(i) = .true.
          if (ipart.lt.npart) then
             do n=idone+1,nneigh
                j = listneigh(n)
                if (j.le.npart .and. .not.sorted(j)) then
                   ipart = ipart + 1
                   ilist(ipart) = j
                   sorted(j) = .true.
                endif
             enddo
          endif
       endif
       i = ll(i)
    enddo
 enddo
 print*,'ncellsloop = ',ncellsloop,' ipart = ',ipart, 'npart= ',npart
 !
 !--do not re-order ghosts (but ireal will be wrong)
 !
! do i=1,npart
!    ilist(i) = i
! enddo
 do i=npart+1,size(ilist)
    ilist(i) = i
 enddo
 
 call alloc(size(rho),sortlist=ilist)
 call cpu_time(t2)
 write(iprint,*) 'completed in ',t2-t1,'s'

end subroutine sort_particles

subroutine indexx(n, arr, indx)
!************************************************************
!                                                           *
!  This is INDEXX using the quicksort algorithm.            *
!                                                           *
!************************************************************
 implicit none
 integer, parameter :: m=7, nstack=500
 integer, intent(in) :: n
 real, dimension(n), intent(in) :: arr
 integer, dimension(n), intent(out) :: indx

 integer :: i,j,k,l,ir,jstack,indxt,itemp
 integer, dimension(nstack) :: istack
 real :: a

 do j = 1, n
    indx(j) = j
 enddo
 jstack = 0
 l = 1
 ir = n

1 if (ir - l.lt.m) then
   do j = l + 1, ir
      indxt = indx(j)
      a = arr(indxt)
      do i = j - 1, 1, -1
         if (arr(indx(i)).le.a) goto 2
         indx(i + 1) = indx(i)
      end do
      i = 0
2     indx(i + 1) = indxt
   end do
   if (jstack.eq.0) return
   ir = istack(jstack)
   l = istack(jstack - 1)
   jstack = jstack - 2
  else
   k = (l + ir)/2
   itemp = indx(k)
   indx(k) = indx(l + 1)
   indx(l + 1) = itemp
   if (arr(indx(l + 1)).gt.arr(indx(ir))) then
      itemp = indx(l + 1)
      indx(l + 1) = indx(ir)
      indx(ir) = itemp
   endif
   if (arr(indx(l)).gt.arr(indx(ir))) then
      itemp = indx(l)
      indx(l) = indx(ir)
      indx(ir) = itemp
   endif
   if (arr(indx(l + 1)).gt.arr(indx(l))) then
      itemp = indx(l + 1)
      indx(l + 1) = indx(l)
      indx(l) = itemp
   endif
   i = l + 1
   j = ir
   indxt = indx(l)
   a = arr(indxt)

3  continue
   i = i + 1
   if (arr(indx(i)).lt.a) goto 3
4  continue
   j = j - 1
   if (arr(indx(j)).gt.a) goto 4
   if (j.lt.i) goto 5
   itemp = indx(i)
   indx(i) = indx(j)
   indx(j) = itemp
   goto 3

5  indx(l) = indx(j)
   indx(j) = indxt
   jstack = jstack + 2
   if (jstack.gt.nstack) then
      print*,'fatal error!!! stacksize exceeded in sort'
      print*,'need to set parameter nstack higher in subroutine indexx '
      stop
   endif
   if (ir - i + 1.ge.j - l) then
      istack(jstack) = ir
      istack(jstack - 1) = i
      ir = j - 1
   else
      istack(jstack) = j - 1
      istack(jstack - 1) = l
      l = i
   endif
 endif

goto 1
end subroutine indexx
