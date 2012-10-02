module getvbary
 use get_neighbour_lists
 implicit none
 
contains

!!------------------------------------------------------------------------
!! Computes the barycentric velocity in two steps
!! Indicate which one is computed
!! This routine is used for calculating high drag regimes
!!
!!istep = 1: first part of the Kdrag barycentric
!!istep = 2: second part of the Kdrag barycentric
!!istep = 3: Monaghan-like type of barycentric
!!------------------------------------------------------------------------

subroutine get_vbary(istep,vbary_in,itype,x,hh,pmass,rho,vbary_out,dubarydt_out)
 use dimen_mhd, only:ndim,ndimV
 use debug,     only:trace
 use loguns,    only:iprint

 use kernels,   only:interpolate_kerneldrag
 use options,   only:idrag_nature,idrag_structure,Kdrag
 use part,      only:npart,ntotal,spsound,itypegas
 use eos,       only:gamma
 use linklist

 implicit none
 integer, intent(in)                        :: istep
 real, dimension(ndimV,ntotal), intent(in)  :: vbary_in
 integer, dimension(ntotal), intent(in)     :: itype
 real, dimension(ndim,ntotal), intent(in)   :: x
 real, dimension(ntotal), intent(in)        :: hh,pmass,rho
 real, dimension(ndimV,ntotal), intent(out) :: vbary_out
 real, dimension(ntotal), intent(out)       :: dubarydt_out
!
!--define local variables
!
 integer :: i,j,n
 integer :: icell,iprev,nneigh
 integer :: idone
 integer, dimension(ntotal) :: listneigh ! neighbour list
!
!  (drag quantities)
!
 integer :: itypei,itypej
 logical :: iskip_drag
 real    :: coeff_gei_1,coeff_dq_1,coeff_dq_4
 real    :: pmassi,rhoi,dubarydt_outi
 real    :: pmassj,pmassij,rhoj,rhoij
 real    :: dv2,vij,V0,f,dragcoeff,dragterm
 real    :: s2_over_m,spsoundgas
 real    :: hi,hi1,hi21,hfacwabi
 real    :: hj,hj1,hj21,hfacwabj
 real    :: q2i,q2j,rij,rij1,rij2
 real    :: wabi,wabj,wab
 real, dimension(ndim)  :: xi,xj,dx
 real, dimension(ndimV) :: vbary_ini,vbary_inj,dvel,dr,drdrag,vbary_outi
 real, parameter :: pow_drag_exp = 0.4
 real, parameter :: a2_set       = 0.5
 real, parameter :: a3_set       = 0.5
 real, parameter :: pi  = 3.141592653589

!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' entering subroutine get_vbary,npart,ntotal=' &
                            ,npart,ntotal
!
!--sanity check
!
 if ((istep.ne.1).and.(istep.ne.2).and.((istep.ne.3))) then
    print*,'in get_vbary, istep ne to 1 or 2 or 3',istep
    stop
 endif
!
!--initialise quantities
!
 listneigh     = 0
 iskip_drag    = .false.
 coeff_gei_1   = 4./3.*sqrt(8.*pi/gamma)
 coeff_dq_1    = sqrt(0.5*gamma)
 coeff_dq_4    = 9.*pi*gamma/128.    
 s2_over_m     = 1./coeff_gei_1
 f             = 0.
 V0            = 0.
 dragcoeff     = 0.
 dr(:)         = 0.
 do i=1,ntotal !--clean the globally cumulative quantities
    vbary_out(:,i)  = 0.
    dubarydt_out(i) = 0.
 enddo
!
!--loop over all the link-list cells
!
 loop_over_cells: do icell=1,ncellsloop          ! step through all cells
!
!--get the list of neighbours for this cell 
!  (common to all particles in the cell)
!
    call get_neighbour_list(icell,listneigh,nneigh)
!
!--now loop over all particles in the current cell
!
    i = ifirstincell(icell)
    idone = -1     ! note density summation includes current particle
    if (i.ne.-1) iprev = i

    loop_over_cell_particles: do while (i.ne.-1) ! loop over home cell particles

       idone = idone + 1
       xi(:)        = x(:,i)
       vbary_ini(:) = vbary_in(:,i)
       vbary_outi(:) = 0.
       dubarydt_outi  = 0. 
       itypei       = itype(i)
       pmassi       = pmass(i)
       rhoi         = rho(i)
       hi           = hh(i)
       hi1          = 1./hi
       hi21         = hi1*hi1
       hfacwabi     = hi1**ndim
!
!--for each particle in the current cell, loop over its neighbours
!
       loop_over_neighbours: do n = idone+1,nneigh
          j      = listneigh(n)
          itypej = itype(j)
!          if ((j.ne.i) .and. .not.(j.gt.npart .and. i.gt.npart) .and. (itypei.ne.itypej)) then
          if ((j.ne.i).and. (itypei.ne.itypej)) then
             ! do count particle with itself
             xj(:)    = x(:,j)
             dx(:)    = xi(:) - xj(:)
             rij2     = dot_product(dx,dx)
             rij      = sqrt(rij2)
             rij1     = 1./rij !eventually infinity
             dr(1:ndim) = rij1*dx(1:ndim)
             hj       = hh(j)
             hj1      = 1./hj !!1./hj
             hj21     = hj1*hj1
             hfacwabj = hj1**ndim
             q2i = rij2*hi21
             q2j = rij2*hj21
            
!
!--get informations on the differential velocities
!
              vbary_inj(:) = vbary_in(:,j)
              if ((istep.eq.1).or.(istep.eq.3)) then
                 dvel(:)   = vbary_ini(:) - vbary_inj(:)
              else
                 dvel(:)   = 0.5*(vbary_ini(:) - vbary_inj(:))
              endif
              dv2          = dot_product(dvel,dvel) 

              if (rij.le.epsilon(rij)) then !two particles at the same position
                 if (dv2.le.epsilon(dv2)) then !no differential velocity => no drag
                     drdrag(:)      = 0.
                     iskip_drag     = .true.
                 else ! Change dr so that the drag is colinear to the differential velocity
                     vij            = sqrt(dv2)
                     drdrag(:)      = dvel(:)/vij
                 endif
              else ! dr = drdrag
                 drdrag(:) = dr(:)
              endif

!---start the drag calculation
   !if (.not.iskip_drag) then
!
!--get the j particle extra properties
!
                pmassj    = pmass(j)
                pmassij   = pmassi + pmassj
                rhoj      = rho(j)
                rhoij     = rhoi+rhoj !Careful rhoij .ne. rhoi*rhoj here
                if (itypei.eq.itypegas) then
                   spsoundgas = spsound(i)
                else
                   spsoundgas = spsound(j)
                endif
!
!--calculate the kernel(s)
! 
                call interpolate_kerneldrag(q2i,wabi)
                wabi     = wabi*hfacwabi
                call interpolate_kerneldrag(q2j,wabj)
                wabj     = wabj*hfacwabj
                if (itypei.eq.itypegas) then
                   wab = wabi
                else
                   wab = wabj
                endif
!
!--calculate the quantities needed for the drag force
!
                V0 = dot_product(dvel,drdrag)
                select case(idrag_nature)
                    case(1) !--constant drag
                       dragcoeff = Kdrag/(rhoi*rhoj)
                    case(2) !--Epstein regime
                       select case(idrag_structure)
                          case(1,4,5) !--linear, third order, PM expression
                             dragcoeff = coeff_gei_1*s2_over_m*spsoundgas
                          case(3)
                             dragcoeff = pi*s2_over_m
                          case default
                             print*,'ERROR drag calculation: wrong idrag_structure'
                       end select
                end select

                select case(idrag_structure)
                   case(1) !--linear regime
                      f = 1.
                   case(2) !--power law
                      f = dv2**(0.5*pow_drag_exp)
                   case(3) !--quadratic
                      f = sqrt(dv2)
                   case(4) !--cubic expansion
                      select case(idrag_nature)
                         case(1) !--forced drag
                            f = 1. + a3_set*dv2
                         case(2) !--Epstein
                            f = 1. + 0.2*coeff_dq_1*coeff_dq_1*dv2/(spsoundgas*spsoundgas)
                         case default
                            print*,'ERROR drag calculation cubic: wrong drag nature' 
                      end select
                   case(5) !--PM expression
                      select case(idrag_nature)
                         case(1) !--forced drag
                            f = sqrt(1. + a2_set*dv2)
                         case(2) !--Epstein
                            f = sqrt(1. + coeff_dq_4*dv2/(spsoundgas*spsoundgas))
                         case default
                            print*,'ERROR drag calculation PM: wrong drag nature'
                      end select
                   case default
                      print*,'this value for idrag_structure does not exist'
                end select
                

!
!--update the force and the energy
!   
                if (istep.eq.1) then
                   dragterm = ndim*wab*f*V0/dragcoeff
                elseif(istep.eq.2) then             
                   dragterm = - ndim*wab*dragcoeff*f*V0
                elseif (istep.eq.3) then
                   dragterm = ndim*wab*V0/rhoij
                else
                  print*,'no other type of istep in vbary'
                  stop
                endif
                
                if (istep.eq.1) then
                    vbary_outi(:)  = vbary_outi(:) - & 
                                    pmassj/(rhoj*rhoij)*dragterm*drdrag(:)
                    vbary_out(:,j) = vbary_out(:,j) + &
                                    pmassi/(rhoi*rhoij)*dragterm*drdrag(:)
                elseif(istep.eq.2) then
                    vbary_outi(:)  = vbary_outi(:) - pmassj*dragterm*drdrag(:)
                    vbary_out(:,j) = vbary_out(:,j) + pmassi*dragterm*drdrag(:)
                elseif (istep.eq.3) then
                    vbary_outi(:)  = vbary_outi(:) - & 
                                     pmassj*dragterm*drdrag(:)
                    vbary_out(:,j) = vbary_out(:,j) + &
                                     pmassi*dragterm*drdrag(:)
                else
                   print*,'no other type of istep in vbary'
                   stop
                endif
!----------------------------------------------------------------------
                if (itypei.eq.itypegas) then
                   dubarydt_outi   = dubarydt_outi + pmassj/(rhoj*rhoij)*dragterm
                endif
                if (itypej.eq.itypegas) then
                   dubarydt_out(j) = dubarydt_out(j) + pmassi/(rhoi*rhoij)*dragterm
                endif               

         endif! .not. j>npart .and. i>npart   
       enddo loop_over_neighbours

       vbary_out(:,i)  = vbary_out(:,i)  + vbary_outi(:)
       dubarydt_out(i) = dubarydt_out(i) + dubarydt_outi
       
       iprev = i
       if (iprev.ne.-1) i = ll(i)          ! possibly should be only if (iprev.ne.-1)
    enddo loop_over_cell_particles
 enddo loop_over_cells
 return
end subroutine get_vbary
end module getvbary
