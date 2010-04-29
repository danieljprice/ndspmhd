!----------------------------------------------------------------
!     Set up a uniform density cartesian grid of particles in ND
!----------------------------------------------------------------

subroutine setup
!
!--include relevant global variables
!
 use dimen_mhd
 use debug
 use loguns
 use bound
 use options
 use part
 use setup_params
 use eos, only:gamma
 use artvi, only:alphamin
 use uniform_distributions
!
!--define local variables
!            
 implicit none
 integer :: i
 real :: massp,volume,totmass,spsoundi,uuzero,hzero,kfac,dxbnd
 real :: denszero,przero,betazero,bzero,dy,valfveni,viscnu
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup (2D Tearing Mode)'
!
!--set boundaries
!
 ibound(:) = 3  ! periodic in y,z
 ibound(1) = 1  ! fixed in x
 nbpts = 0      ! use ghosts not fixed
 xmin(1) = -40. ! set position of boundaries
 xmax(1) = 40.
 xmin(2) = 0.5*xmin(1) ! give aspect ratio
 xmax(2) = 0.5*xmax(1)

 dxbnd = 5.*hfact*psep
 xmin(1) = xmin(1) - dxbnd
 xmax(1) = xmax(1) + dxbnd
!
!--set up the uniform density grid
!
 call set_uniform_cartesian(1,psep,xmin,xmax,fill=.true.)
 npart = ntotal
 print*,'npart =',npart
!
!--determine particle mass
!
 denszero = 1.0
 volume = product(xmax(:)-xmin(:))
 totmass = denszero*volume
 massp = totmass/float(ntotal) ! average particle mass
 
 !--Set the sound speed, and determine P and u accordingly
 spsoundi = 20.0
 przero = (denszero*spsoundi**2)/gamma
 uuzero = przero/((gamma-1.)*denszero)
 
 !--Quentin's parameters:
 !uuzero = denszero**(2./3.)
 !przero = (gamma-1.)*denszero*uuzero
 !spsoundi = sqrt(gamma*przero/denszero)
 
 hzero  = hfact*psep
 
 !--use these lines to set beta and determine b
 !betazero = 10.0
 !bzero = sqrt(2.*przero/betazero) 

 !--use these lines to set b and determine beta
 bzero = 1.
 betazero = przero/(0.5*bzero**2)

 valfveni = bzero/sqrt(denszero)
 if (iresist.eq.0) print*,'WARNING: iresist should be non-zero for this test'
 write(iprint,*) ' initial sound speed = ',spsoundi,' przero = ',przero
 write(iprint,*) ' beta = ',betazero,' Bzero = ',bzero
 if (ndim.eq.2) then
   viscnu = 1./8.*alphamin*spsoundi*hzero
 elseif (ndim.eq.3) then
   viscnu = 1./10.*alphamin*spsoundi*hzero
 else
   viscnu = tiny(viscnu)
 endif
 write(iprint,*) ' viscous nu = ',viscnu,' resistive eta  = ',etamhd
 dy = xmax(2) - xmin(2)
 write(iprint,*) ' t (sound crossing)     = ',dy/spsoundi
 write(iprint,*) ' t(alfven crossing)     = ',dy/valfveni
 write(iprint,*) ' t(viscous diffusion)   = ',dy**2/viscnu
 write(iprint,*) ' t(resistive diffusion) = ',dy**2/etamhd
!
!--now assign particle properties
!
 kfac = 20.*(xmax(1)-xmin(1))
 do i=1,ntotal
    if (abs(x(1,i)).gt.(xmax(1)-dxbnd)) itype(i) = 1
    dens(i) = denszero
    pmass(i) = massp
    vel(:,i) = 0.
    uu(i) = uuzero
    Bfield(:,i) = 0.
    Bfield(2,i) = bzero*tanh(x(1,i)*kfac)
 enddo
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
  
 return
end
