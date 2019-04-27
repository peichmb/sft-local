! ******************************************************************************
! ****                          MODULE UPDATE                               ****
! ******************************************************************************
! 
! Contains the kernel of the code.
!
! PROCEDURES
! -- update_positions: updates the position of all the features in the
!      linked list.
! -- update umaps: updates the additional velocity maps (ux and uy) with
!      the shear velocity and the aggregate inflow velocity
! -- update_rwu: updates the random walk velocities of the features
! -- cancel: marks the features that are to be cancelled
! -- delete: deletes the marked features

module update
    
    use params
    use ustencil
    use vars
    use grid
    use lsflows
    use random
    implicit none

    INTERFACE 
       FUNCTION ZBQLU01(DUMMY)
         DOUBLE PRECISION :: ZBQLU01
         DOUBLE PRECISION, INTENT(IN) :: DUMMY
       END FUNCTION ZBQLU01
    END INTERFACE        

contains    

    subroutine update_positions(dt,boundary)
        
        logical :: boundary
        integer :: i        
        double precision :: dt, ssdx, ssdy ! dt, small scale dx and dy
        type(feature), pointer :: curr
        
        curr => feature_head
        
        do while (associated(curr%next))
            curr => curr%next
            
            ssdx = 0.d0
            ssdy = 0.d0
            
            if (inc_inflows) then
                ssdx = ssdx + ux(curr%xi,curr%yi)*dt
                ssdy = ssdy + uy(curr%xi,curr%yi)*dt
            end if
            if (inc_rw) then
                ssdx = ssdx + curr%urwx*dt
                ssdy = ssdy + curr%urwy*dt
            end if
            
            curr%x = curr%x + ssdx
            curr%y = curr%y + ssdy
            
            curr%ssx = curr%ssx + ssdx
            curr%ssy = curr%ssy + ssdy
            
            if (inc_mflow) then
                curr%y = curr%y+uy_mf(curr%y,lat0)*dt
            end if
            if (inc_shear) then
                curr%x = curr%x+ux_sh(curr%y,lat0)*dt
            end if
            
            curr%xi = int((curr%x-xmin+dx/2.)/dx)
            curr%yi = int((curr%y-ymin+dy/2.)/dy)
            
            if (curr%xi < 1 .or. curr%xi > nx .or. &
              & curr%yi < 1 .or. curr%yi > ny) then
                boundary = .true.
                print *, curr%x, curr%y, curr%xi, curr%yi
            end if             
        end do
        
    end subroutine update_positions
    
    subroutine update_umaps()

        integer :: i
        type(feature), pointer :: curr
               
!        if (inc_shear) then
!            ux = ushear
!        else
!            ux = 0.d0
!        end if
!        if (inc_mflow) then
!            uy = umflow
!        else
!            uy = 0.d0
!        end if
        
!        if (inc_inflows) then
            ! This aggregates the stencil of inflow velocities centered on each
            !   of the features.
            
            ux = 0.d0
            uy = 0.d0
            curr => feature_head
            do while (associated(curr%next))
                curr => curr%next                
                ux(:,:) = ux(:,:) &
                  & + uxst(nx-curr%xi+1:2*nx-curr%xi,ny-curr%yi+1:2*ny-curr%yi)
                uy(:,:) = uy(:,:) &
                  & + uyst(nx-curr%xi+1:2*nx-curr%xi,ny-curr%yi+1:2*ny-curr%yi)                   
            end do
!        end if
        
!        print *,'-----',sqrt(maxval(ux**2+uy**2))
!        open(32,file='lll')
!        do i=1,nx
!            write(32,*) xx(i),ux(i,ny/2)
!        end do
!        stop
        
    end subroutine update_umaps
   
    subroutine update_rwu(nsteps, nrw)
    
        integer :: i, nsteps, nrw
        double precision :: r1,t
        type(feature), pointer :: curr
        curr => feature_head    
        if (sim_rwc) then
        
            if( mod(nsteps,nrw) == 0) then        
                do while (associated(curr%next))
                    curr => curr%next
                    r1 = ZBQLU01(DUMMY)
                    curr%urwx = dx_rw*cos(2.d0*pi*r1)/dt_rw
                    curr%urwy = dx_rw*sin(2.d0*pi*r1)/dt_rw
                end do
            else
                return
            end if
        
        else                            
            ! The tokens are delivered in the initconds module. These are to
            !   avoid the change of all the random velocities happening at the
            !   same time. Note that this is useless if we have set sim_rw to
            !   true in the parameters file.
            do while (associated(curr%next))
                curr => curr%next
                if ( mod(nsteps-curr%rwtoken,nrw) == 0 ) then
                    r1 = ZBQLU01(DUMMY)
                    curr%urwx = dx_rw*cos(2.d0*pi*r1)/dt_rw
                    curr%urwy = dx_rw*sin(2.d0*pi*r1)/dt_rw
                end if
            end do
        
        end if
            
    end subroutine update_rwu
    
    subroutine cancel()
    
        integer :: i
        double precision :: x1, y1, x2, y2, d2
        logical :: call_delete
        type(feature), pointer :: curr1, curr2
        
        if (.not. associated(feature_head%next)) return
        
        call_delete = .false.
        curr1 => feature_head%next
        
        do while (associated(curr1%next))
            if (curr1%cancel) then
                curr1 => curr1%next
            else
                curr2 => curr1%next
                do
                    if (curr2%cancel) then
                        if (associated(curr2%next)) then
                            curr2 => curr2%next                        
                        else
                            curr1 => curr1%next
                            exit
                        end if
                    else
                        x1 = curr1%x
                        y1 = curr1%y
                        x2 = curr2%x
                        y2 = curr2%y
                        d2 = (x2-x1)**2 + (y2-y1)**2 
                        if (d2 < ctres2 .and. curr1%pol /= curr2%pol) then
                            call_delete = .true.
                            curr1%cancel = .true.
                            curr2%cancel = .true.
                            curr1 => curr1%next
                            exit
                        else
                            if (associated(curr2%next)) then
                                curr2 => curr2%next
                            else
                                curr1 => curr1%next
                                exit
                            end if
                        end if
                    end if
                end do
            end if
        end do
                        
        if (call_delete) call delete()            
            
    end subroutine cancel
        
    subroutine delete()
            
        type(feature), pointer :: curr1, curr2
        
        curr1 => feature_head
        
        if (.not. associated(curr1%next)) then
            print *, 'Empty list'
            return
        else   
            do while (associated(curr1%next))
                curr2 => curr1%next
                if (curr2%cancel) then
                    if (.not. associated(curr2%next)) then
                        ! End of list
                        nullify(curr1%next)
                        deallocate(curr2)                      
                    else
                        curr1%next => curr2%next                       
                        deallocate(curr2)
                    end if                    
                else
                    if (.not. associated(curr2%next)) then
                        ! End of list
                        exit                    
                    else
                        curr1 => curr1%next
                    end if
                end if
            end do
        end if
    
    end subroutine delete

end module update
