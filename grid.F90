! ******************************************************************************
! ****                               MODULE GRID                            ****
! ******************************************************************************

! This module contains the grid-related variables.
! The grid is used to add up all the inflow and shear velocities; with this
!   approach we save computations.
!
! VARIABLES:
! -- xx, yy: Coordinates
! -- ux, uy: x- and y- components of the velocity
! -- dx, dy: x- and y- size of cells
!
! PROCEDURES:
! -- subroutine init_grid: Initializes the grid variables

module grid
    
    use params
    implicit none
    double precision :: dx, dy
    double precision, allocatable, dimension(:) :: xx, yy
    double precision, allocatable, dimension(:,:) :: ux, uy
    
contains
    
    subroutine init_grid()
    
        integer :: i

        dx = (xmax-xmin)/(nx-1)
        dy = (ymax-ymin)/(ny-1)
        
        allocate(xx(nx))
        allocate(yy(ny))
        allocate(ux(nx,ny))
        allocate(uy(nx,ny))  
        
        ux = 0.d0
        uy = 0.d0      
        
        do i=1,nx
            xx(i)=xmin + dx*(i-1)
        end do
        do i = 1, ny
            yy(i)=ymin + dy*(i-1)
        end do
        
        print *
        print *, '** init_grid **'
        print *, ' - dx = ', dx
        print *, ' - dy = ', dy
        print *, ' - |dx-dy|/dx (%) = ', abs(dx-dy)/dx*100.d0
        
    end subroutine init_grid
    
end module grid
