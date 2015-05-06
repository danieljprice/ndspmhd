!------------------------------------------------------------------------------!
! NDSPMHD: A Smoothed Particle (Magneto)Hydrodynamics code for (astrophysical) !
! fluid dynamics simulations in 1, 2 and 3 spatial dimensions.                 !
!                                                                              !
! (c) 2002-2015 Daniel Price                                                   !
!                                                                              !
! http://users.monash.edu.au/~dprice/ndspmhd                                   !
! daniel.price@monash.edu -or- dprice@cantab.net (forwards to current address) !
!                                                                              !
!  NDSPMHD comes with ABSOLUTELY NO WARRANTY.                                  !
!  This is free software; and you are welcome to redistribute                  !
!  it under the terms of the GNU General Public License                        !
!  (see LICENSE file for details) and the provision that                       !
!  this notice remains intact. If you modify this file, please                 !
!  note section 2a) of the GPLv2 states that:                                  !
!                                                                              !
!  a) You must cause the modified files to carry prominent notices             !
!     stating that you changed the files and the date of any change.           !
!                                                                              !
!  ChangeLog:                                                                  !
!------------------------------------------------------------------------------!

!!-------------------------------------------------------------------------
!! update of particles at end of step for periodic (particles crossing
!! domain) and inflow/outflow (adds / subtracts particles from the domain)
!!
!! for ND case only periodic boundaries implemented
!! NOTE: AT THE MOMENT THIS DOES *NOT* DO MEMORY REALLOCATION IF npart > ARRAY SIZE
!!-------------------------------------------------------------------------
       
SUBROUTINE boundary
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE bound
 USE options
 USE part
 USE part_in

! USE setup_params
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i,nnew,nsub,nsubtemp,jdim,npartin
 REAL :: psepleft,psepright
 LOGICAL :: debugging
!
!--allow for tracing flow
!      
 debugging = .true.
 IF (trace) THEN
    WRITE(iprint,*) ' Entering subroutine boundary'
    debugging = .true.
 ENDIF   
 
!---------------------------------------------------------------------------
!  inflow/outflow when fixed particles are used.
!---------------------------------------------------------------------------     
 IF (ANY(ibound.EQ.1)) THEN

    IF (ndim.GT.1) THEN    
!       WRITE(iprint,*) 'warning: inflow/outflow not implemented in ND'
    ELSE      
    
    npartin = npart
!
!--this is just pulled straight from the 1D code - needs to be
!  modified if applied in 2 or 3D since it assumes that particles
!  are ordered from left to right in the 1D domain

!
!--determine particle separation at left and right boundaries
!
    psepleft = x(1,2) - x(1,1)
    psepright = x(1,npart) - x(1,npart-1)
    
    IF (psepright.LT.0.) WRITE(iprint,*) 'Boundary: warning psepright < 0', &
                               psepright,x(1,npart-1),x(1,npart),npart
    IF (psepleft.LT.0.) WRITE(iprint,*) 'Boundary: warning psepleft < 0',psepleft,x(1,1),x(1,2)
!
!--inflow to left boundary
!
    IF (x(1,1).GT.(xmin(1)+psepleft)) THEN
       IF (debugging) WRITE(iprint,*) 'inflow to left boundary npart = ',npart+1
       IF ((x(1,1)-xmin(1)).GT.2*psepleft) THEN
          WRITE(iprint,*)' Need more than one particle',psepleft,(x(1,1)-xmin(1))/psepleft
          nnew = INT((x(1,1)-xmin(1))/psepleft)
          WRITE(iprint,*)' Number of particles needed = ',nnew
       ENDIF        
       DO i=npart,1,-1            ! relabel particles
!          print*,' particle ',i+1,' = ',i
!--copy both primitive and conservative variables
          x(1,i+1) = x(1,i)            ! (must copy all particle properties)
          itype(i+1) = itype(i)        
          call copy_particle(i+1,i)
       ENDDO       
!
!--now make new particle number 1
!  (particle properties are automatically set in the next iteration of step)
!
       npart = npart + 1      ! add new particle
       x(1,1) = x(1,1) - psepleft
       call copy_particle(1,2)  ! copy quantities from particle 2
       itype(1:nbpts) = 1
       itype(nbpts+1) = 0
!
!       PRINT*,' new particle x(',1,') = ',x(1)
!
!--outflow from left boundary
!       
    ELSEIF (x(1,1).LT.xmin(1)) THEN
       nsub = INT((xmin(1)-x(1,1))/psepleft) !+ 1
       IF (debugging .and. nsub.gt.0)       &
          WRITE(iprint,*) 'outflow from left boundary npart = ',npart-nsub,nsub
       IF (nsub.GT.0) THEN
          npart = npart - nsub      ! subtract particle(s)
          DO i=1,npart            ! relabel particles         
             x(1,i) = x(1,i+nsub)             ! (must copy all particle properties)
             call copy_particle(i,i+nsub)
          ENDDO
       ENDIF
              
    ENDIF
!
!--right boundary
!    
    IF (x(1,npart).GT.xmax(1)) THEN            ! outflow from right boundary
       nsub = INT((x(1,npart)-xmax(1))/psepright) !+ 1
       IF (debugging) WRITE(iprint,*) 'outflow from right boundary npart = ',npart-nsub,nsub       
!--if more than one particle, check other particles close to boundary
       IF (nsub.GE.2) THEN
          WRITE(iprint,*) 'need more than one particle ',nsub
          nsubtemp = 0
          DO i=0,nsub
             IF (x(1,npart-i).GT.xmax(1)) THEN
                nsubtemp = nsubtemp + 1
             ENDIF
          ENDDO
          nsub = nsubtemp
          WRITE(iprint,*) 'really need more than one particle ',nsub
       ENDIF        
       npart = npart - nsub      ! subtract particle(s)
    ELSEIF (x(1,npart).LT.(xmax(1)-psepright)) THEN      ! inflow from right boundary   
       IF (debugging) WRITE(iprint,*) 'inflow to right boundary npart = ',npart+1
       IF ((x(1,npart)-xmax(1)).GT.2*psepright) THEN
          WRITE(iprint,*)' Need more than one particle',psepright,(x(1,npart)-xmax(1))/psepright
          nnew = INT((x(1,npart)-xmax(1))/psepright)
          WRITE(iprint,*)' Number of particles needed = ',nnew
       ENDIF        

       npart = npart + 1      ! add only one new particle
!
!--create new particle number npart
!  (force etc calculated next step - doesn't matter so long as its a fixed
!   particle to start with)
!
       x(1,npart) = x(1,npart-1) + psepright
       call copy_particle(npart,npart-1)
       itype(npart) = 1
       itype(npart-nbpts) = 0
       
    ENDIF
!
!  adjust ireal and itype if particles have been added/subtracted
!    
    IF (npartin.ne.npart) THEN
       ireal = 0
       itype = 0
       itype(1:nbpts) = 1
       ireal(1:nbpts) = nbpts+1
       itype(npart-nbpts+1:npart) = 1
       ireal(npart-nbpts+1:npart) = npart-nbpts
!       WRITE(iprint,*) ' New number of particles = ',npart
    ENDIF

    ntotal = npart
!
!--make sure all the initial quantities are the same
!
    DO i=1,npart
       xin(1,i) = x(1,i)
       velin(:,i) = vel(:,i)
       rhoin(i) = rho(i)
       hhin(i) = hh(i)
       enin(i) = en(i)
       Bevolin(:,i) = Bevol(:,i)
       alphain(:,i) = alpha(:,i)
!
!--check to see if particles are ordered left to right
!       
!       IF (i.LT.npart.AND.x(1,i).GT.x(1,i+1)) THEN
!          WRITE(iprint,*) 'Warning: particles crossed ',i,x(1,i),i+1,x(1,i+1)       
!       ENDIF
    ENDDO
        
    ENDIF

 ENDIF
!----------------------------------------------------------------------
! periodic boundary conditions - allow particles to cross the domain
!----------------------------------------------------------------------

 IF (ANY(ibound.EQ.3)) THEN
    
    DO i=1,npart
! -- this is very f77 but I can't get it right using where statements        
       DO jdim=1,ndim
          IF (ibound(jdim).EQ.3) THEN
             IF (x(jdim,i).GT.xmax(jdim)) THEN
!               print*,'ss xold,xmax,xnew = ',jdim,x(jdim,i),xmax(jdim),xmin(jdim) + x(jdim,i) - xmax(jdim)
                x(jdim,i) = xmin(jdim) + x(jdim,i) - xmax(jdim)
             ELSEIF(x(jdim,i).LT.xmin(jdim)) THEN
!              print*,'ss xold,xmin,xnew = ',jdim,x(jdim,i),xmin(jdim),xmax(jdim) + x(jdim,i) - xmin(jdim)           
!                read*
                x(jdim,i) = xmax(jdim) - (xmin(jdim) - x(jdim,i))
             ENDIF        
        ENDIF
       ENDDO
    
    ENDDO   

 ENDIF

 RETURN 
END SUBROUTINE boundary
