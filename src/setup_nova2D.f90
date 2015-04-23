!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Setup for the blob evaporation problem                                !!
!!                                                                        !!
!!  dense disc of fluid at origin                                         !!
!!                                                                        !!
!!------------------------------------------------------------------------!!

subroutine setup
!
!--include relevant global variables
!
 use dimen_mhd
 use debug
 use loguns
 use bound
 use eos
 use options
 use part
 use setup_params
 
 use uniform_distributions
 use cons2prim, only:primitive2conservative
!
!--define local variables
!      
 implicit none
 integer :: ipart
 real :: rs,rhoc,rhomid,rhozero,vxc,vxmid,vxzero,prc,prmid,przero,gamma1,gamma2
 real :: xi,yi,totmass,totvol,psepc,psepmid,massp,pri,gam1,Rfunc,pmassi
 real, dimension(ndim) :: xtemp
 real, dimension(ndimV) :: Bzero
 logical, parameter :: equalmass = .true.

 Rfunc(yi) = -50. + 10.*cos(2.*pi*yi/100.)

!
!--set boundaries
!                        
 ibound(1) = 1  ! fixed
 ibound(2) = 3  ! periodic
 nbpts = 0        ! no fixed particles initially
 xmin(1) = -300.0
 xmax(1) = 100.0
 xmin(2) = -50.0
 xmax(2) = 50.0
 rs= 0.0
 rhoc = 0.12
 rhomid = 1.7
 rhozero = 5.88713
 vxc = 59.44
 vxmid = 59.44
 vxzero = 24.03
 prc = 10.0
 prmid = 10.0
 przero = 3007
 Bzero(:) = 0.
 gamma1 = 1.45
 gamma2 = 1.8
 iprterm = 11
!
!--setup parameters for the problem
! 
 psepc = psep*(rhozero/rhoc)**(1./ndim)
 psepmid = psep*(rhozero/rhomid)**(1./ndim)

 totvol = product(xmax - xmin)
 totmass = rhozero*totvol
 gam1 = gamma - 1.
 !pext = przero

 write(iprint,*) 'NOVA laser test '
 write(iprint,10) rs,gamma1,gamma2
10 format(/,' rs  = ',f10.3,', gamma = ',f6.3,f6.3/)
!
!--get mass by fake uniform setup
!
 call set_uniform_cartesian(1,psep,xmin,xmax)
 massp = totmass/real(npart)

 if (equalmass) then
    npart = 0
!
!--setup 2D particle distribution using masks
!  (determines particle number and allocates memory)
!
    xtemp(:) = xmax(:)
    xtemp(1) = rs
    call set_uniform_cartesian(1,psepc,xmin,xtemp,mask=-4,fill=.true.)
    call set_uniform_cartesian(1,psepmid,xmin,xtemp,mask=4,fill=.true.)

    xtemp(:) = xmin(:)
    xtemp(1) = rs
    call set_uniform_cartesian(1,psep,xtemp,xmax,fill=.true.)
 endif
 ntotal = npart
!
!--now assign particle properties
!
 do ipart=1,ntotal
    xi = x(1,ipart)
    yi = x(2,ipart)
    vel(:,ipart) = 0.
    if (xi.lt.Rfunc(yi) .and. xi.lt.rs) then
       vel(1,ipart) = vxc
       dens(ipart) = rhoc
       pri = prc
       psi(ipart) = gamma1
       pmassi = massp*rhoc/rhozero
    elseif (xi.lt.rs) then
       vel(1,ipart) = vxmid
       dens(ipart) = rhomid
       pri = prmid
       psi(ipart) = gamma2
       pmassi = massp*rhomid/rhozero
    else
       vel(1,ipart) = vxzero
       dens(ipart) = rhozero
       pri = przero
       psi(ipart) = gamma2
       pmassi = massp
    endif
    if (xi.lt.(xmin(1)+3.*hfact*psep) .or. xi.gt.xmax(1)-3.*hfact*psep) then
       itype(ipart) = 1
    endif
    hh(ipart) = 1.2*psep
    if (equalmass) then
       pmass(ipart) = massp
    else
       pmass(ipart) = pmassi
    endif
    if (imhd.ne.0) then
       Bfield(:,ipart) = Bzero(:)
    endif
    uu(ipart) = pri/((psi(ipart)-1.)*dens(ipart))
 enddo

 Bconst(:) = Bzero(:)
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
            
 return
end subroutine setup

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
