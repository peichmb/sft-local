! ******************************************************************************
! ****                          MODULE VARS                                 ****
! ******************************************************************************
!
! Contains the definition of the basic linked list cell
!   and the head of the list.

module vars

    use params
    use grid
    use lsflows
    implicit none
    
    type feature
        double precision :: x,y ! position (x and y coordinates)
        double precision :: urwx, urwy ! random walk velocity (x and y comps.)
        double precision :: ssx, ssy ! Small scale displacement (i.e. without large scale flows)
        integer :: xi, yi ! Map coordinates (i.e. grid cell numbers)
        integer :: pol, rwtoken, id ! Polarity, token for rw vel. change and id.
        logical :: cancel ! Marker to know whether to delete the feature
        type(feature), pointer :: next ! Pointer to the next cell of the list.
    end type feature
    
    type(feature), pointer :: feature_head
    
contains

    subroutine init_vars()
        
        allocate(feature_head)
        nullify(feature_head%next)
        feature_head%x = 0.d0
        feature_head%y = 0.d0
        feature_head%id = 0
        feature_head%xi = 0
        feature_head%yi = 0
        feature_head%urwx = 0.d0
        feature_head%urwy = 0.d0
        feature_head%ssx = 0.d0
        feature_head%ssy = 0.d0
        feature_head%pol = 0
        feature_head%rwtoken = 0
        feature_head%cancel = .false.
        
    end subroutine init_vars
    
end module vars
