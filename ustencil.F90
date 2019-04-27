! ******************************************************************************
! ****                          MODULE USTENCIL                             ****
! ******************************************************************************
!
! Contains the stencils of inflow velocities, ie,  maps of velocities that are
!   calculated once at the beginning of the program and then placed on the 
!   velocity maps (ux, uy). The central cell of the stencil placed right on top
!   of the map cell corresponding to the position of the feature (xi, yi).
!

module ustencil
    
    use params
    use grid
    implicit none
    
    integer :: nxst, nyst
    double precision :: lst, xmin_st, xmax_st, ymin_st, ymax_st
    double precision, allocatable, dimension(:) :: xxst, yyst
    double precision, allocatable, dimension(:,:) :: uxst
    double precision, allocatable, dimension(:,:) :: uyst
    
contains

    subroutine init_stencil()

        integer :: i, j
        double precision, external :: uinfx, uinfy
                
        nxst = 2*nx - 1
        nyst = 2*ny - 1
        
        xmin_st = -(nxst/2)*dx
        xmax_st = (nxst/2)*dx
        ymin_st = -(nyst/2)*dy
        ymax_st = (nyst/2)*dy
        
        if (.not. inc_inflows) return

        allocate(xxst(nxst))
        allocate(yyst(nyst))
        allocate(uxst(nxst,nyst))
        allocate(uyst(nxst,nyst))
        
        do i=1,nxst
            xxst(i)=xmin_st + dx*(i-1)
        end do
        do j=1,nyst
            yyst(j)=ymin_st + dy*(j-1)
        end do
        
        do j=1,nyst
        do i=1,nxst
            call uinf(xxst(i), yyst(j), uxst(i,j), uyst(i,j))
        end do
        end do
        
        print *
        print *, '** init_stencil **'
        print *, ' - Number of points in stencil: ', nxst, ' x ', nyst
        print *, ' - Number of points within 2sigma: ', 2.d0*sigma/dx
        print *, ' - xmin_st = ', xmin_st
        print *, ' - xmax_st = ', xmax_st
        print *, ' - ymin_st = ', ymin_st
        print *, ' - ymax_st = ', ymax_st
        print *, ' - xmid = ', xxst(nxst/2+1)
        print *, ' - ymid = ', yyst(nyst/2+1)
          
    end subroutine init_stencil
    
    function bfield(r)

        double precision :: r, bfield
        bfield = exp(-r**2/2.d0/sigma**2)
                
        return
        
    end function bfield
    
    subroutine uinf(x,y,uix,uiy)
    
        double precision :: x, y, uix, uiy
        double precision :: r, phi, bz
        
        r = sqrt(x**2+y**2)
        phi = atan2(y,x)
        bz = bfield(r)
        
        uix = - ifac*alpha*beta*r/sigma**2*bz**beta*cos(phi)
        uiy = - ifac*alpha*beta*r/sigma**2*bz**beta*sin(phi)
        
    end subroutine uinf
    
end module ustencil
    



