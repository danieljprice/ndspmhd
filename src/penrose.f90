module penrosetile

contains

subroutine penrose(np,x)
 implicit none
 integer, intent(out) :: np
 real, dimension(:,:), intent(out) :: x
 
!******************************************************************
!  This program calculate the vertex positions of 2D quasiperiodic tilings
!  using the generalized dual method. It is just a simple program
!  from were you can start other applications like the 3D version (it is just
!  very easy to do,from the present code).
!  First, you give the star vectors e_{i},
!  i.e. a set of unitary vectors that have the required symmetry. For example,
!  in the Penrose it is required 5 fold symmetry, thus, the star is a set
!  of 5 vectors that point to the vertex of a pentagon. Then, all intersections
!  of a family of lines perpendicular to each vector is calculated. Following
!  the formulas of G.G. Naumis, J.L. Aragón,"Analytic expressions for the 
!  vertex coordinates of quasiperiodic lattices",Z.  Kristallogr. 218, 397 (2003).
!  This program has been made by Gerardo G. Naumis and José Luis Aragon,
!  in 2002, at the Instituto de Fisica, National University of Mexico.
!  To draw the tiling, just make a program that joins all points separated
!  by a distance "1". See our mathematica version of the program in the
!  same webpage:  www.fisica.unam.mx/naumis
!  If you have comments, contact: G. Naumis naumis@fisica.unam.mx
!*********************************************************************
!  Variable declaration
 real ex(5),ey(5),ux(5),uy(5),gamma(5),rx,ry,arx0,ary0,aux,al,pi
 real aij,arx,ary
! common h,ntp,u,v,pi,n0,ntg,n,nfr
 real    :: arxp,aryp,co1,co2,ri
 integer :: n,n1
 integer :: i,j,l,m,ipart
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!    Parameters
 pi=4.0*atan(1.0)
 al=2.0*pi/5.0
!*******************************************************************
!  Following the notation of Steinhardt and Levine, on needs to define
!  phases that shifts the family of lines on each vector. 
!  Different isomorphisms (different types of tilings) are obtained 
!  if these phases are changed. If the sum of these phases is 0, 
!  one obtain the Penrose isomorphism class, which is the most famous.
!  see: J Socolar, Steinhardt, Phys. Rev. B 34, 617 (1986).
!  However, many people forget to tell you that you must be careful 
!  with singular tilings,i.e., if these shifts are so that
!  one has points where three of more lines cross, the tiles are no longer
!  rhombuses, one can have hexagons (duals of three ´lines that meet in a point)
!  octagons (four lines), decagon (5 lines). This last one ocuurs at the origin
!  if all phases are set to zero.
!  A FAQ is "why I don´t get the right Penrose?. First, you will always
!  get different portions of the Penrose tiling unless you know which are the
!  exact parameters of the tile that you want to compare. Usually,instead one
!  needs to check that both tilings are in the same local isomorphism class,
!  that means that each patch of a given radius, is observed in both tilings.
!  To check the Penrose, as first check, verify the proportions of each 
!  vertex configuration and second neighbour. A table was published in:
!   G.G. Naumis, "Density of states and first spectral moments of a 
!  quasiperiodic lattice", J. Phys: Condens. Matter 11 (1999) 7 143-7 154.
!  I an exotic symmetry is required, change the number of star vectors
!************************************************************************  
 gamma(1)=0.0
 gamma(2)=-0.50
 gamma(3)=0.25
 gamma(4)=0.25
 gamma(5)=0.44
!************************************************************************
!      Definition of star (basis vectors), and generation of auxiliary
!      vectors that allow to get the intersection of each family of lines
!************************************************************************
 do i=1,5
    ri=real(i)-1.0
    ex(i)=cos(al*ri)
    ey(i)=sin(al*ri)
    ux(i)=-sin(al*ri)
    uy(i)=cos(al*ri)
    write(*,*) ex(i),ey(i)
 enddo
 
 ipart = 0
!**************************************************************************
! A loop is carried over I and J which are couples of vectors that generate
! all the intersections. For each copule I and J, we have a family of lines
! indexed by m and n, which are the "numeration" of each line over the directions
! I and J.
!**************************************************************************
 do i=1,5
    do j=i+1,5
       if (i.ne.j) then 
          aij=((ex(i)*ey(j))-(ey(i)*ex(j)))

          do n=-10,10
             do m=-10,10
                rx=(real(n)*ex(i))+(real(m)*ex(j))
                ry=(real(n)*ey(i))+(real(m)*ey(j))
                arx0=(-(real(m)+gamma(j))*ey(i)+(real(n)+gamma(i))*ey(j))/aij
                ary0=(((real(m)+gamma(j))*ex(i)-(real(n)+gamma(i))*ex(j)))/aij
                co1=(arx0*ex(i))+(ary0*ey(i))-real(n)-gamma(i)
                co2=(arx0*ex(j))+(ary0*ey(j))-real(m)-gamma(j)
!***************************************************************************
! Write the family and interseection coordinates of two lines,
! which in fact correspond to the center of a rhombus.
!***************************************************************************
!         write(*,*) 'N',N,' M',M,' INT=',arx0,ary0
                if ((ABS(co1).GT.0.0001).OR.(ABS(co2).GT.0.0001)) then
                   write(*,*) 'error',I,J,n,m,co1,co2
                endif
!         write(*,*) 'I',I,' J',J,' N',N,' M',M,co1
!****************************************************************************
! Dual transformation 
! Calculate ordinal coordinates with respect to each family of lines
! of the generated intersection. Of course, coordinates I and J are known
! already, since they were used to generate the point. Thus, such step
! is not performed for lines of family I and J.
!****************************************************************************
               arxp=rx
               aryp=ry
               overl: do l=1,5
                  if ((l.ne.i).and.(l.ne.j)) then
                     aux=(arx0*ex(l))+(ary0*ey(l))-gamma(l)
                     if (aux.ge.0.0) then
                        n1=aint(aux)
                     else
                        n1=aint(aux)-1
                     endif
!******************************************************************
! Perform in a recursive way the dual transformation by making a linear
! combination with the right coefficients of the star vectors.
!******************************************************************
                     arx=rx+(real(n1)*ex(l))
                     ary=ry+(real(n1)*ey(l))
                     arxp=arxp+((aux-0.5)*ex(l))
                     aryp=aryp+((aux-0.5)*ey(l))
                     rx=arx
                     ry=ary
                  endif
               enddo overl
!  Once a point is found in the tiling, the corresponding other 3 points of
!  a rhombus are generated by adding the proper vectors
               if ((abs(rx).lt.9.0).and.(abs(ry).lt.9.0)) then
!                 write(6,*) arx0,ary0
                  ipart = ipart + 1
                  x(1,ipart) = rx
                  x(2,ipart) = ry
                  ipart = ipart + 1
                  x(1,ipart) = rx-ex(i)
                  x(2,ipart) = ry-ey(i)
                  ipart = ipart + 1
                  x(1,ipart) = rx-ex(j)
                  x(2,ipart) = ry-ey(j)
                  ipart = ipart + 1
                  x(1,ipart) = rx-ex(i)-ex(j)
                  x(2,ipart) = ry-ey(i)-ey(j)
                  
                 ! write(4,*) rx,ry,1.0
                 ! write(4,*) rx-ex(i),ry-ey(i),1.0
                 ! write(4,*) rx-ex(j),ry-ey(j),1.0
                 ! write(4,*) rx-ex(i)-ex(j),ry-ey(i)-ey(j),1.0
!                 write(2,*) arxp,aryp,1.0 
!                 write(2,*) arxp-EX(i),aryp-EY(i),1.0  
!                 write(2,*) arxp-EX(j),aryp-EY(j),1.0  
!                 write(2,*) arxp-EX(i)-EX(j),aryp-EY(i)-EY(j),1.0     
               endif
            enddo
         enddo
      endif
   enddo
 enddo
 
 np = ipart
 m = 0
 do ipart=1,np
    rx = x(1,ipart)
    ry = x(2,ipart)
    do j=1,ipart-1
       if (abs(rx - x(1,j)).lt.tiny(0.) .and. &
           abs(ry - x(2,j)).lt.tiny(0.)) then
          print*,ipart,' is same position as ',j,rx,ry
          m = m + 1
       endif
    enddo
 enddo
 print*,' same position found ',m,' times'

end subroutine penrose

end module penrosetile
