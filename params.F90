! ******************************************************************************
! ****                          MODULE PARAMS                               ****
! ******************************************************************************

! Contains the problem parameters read from the input file
!   and some derived quantities

module params
    
    implicit none
    
    logical :: inc_shear, inc_rw, inc_inflows, verbose, dtopt, sim_rwc, rwcc, &
             & inc_mflow
    integer :: nft, nx, ny, rinit, rseed, idum
    double precision :: xmin, xmax, ymin, ymax
    double precision :: lat0, radius, sigma, alpha, beta, ifac
    double precision :: tmax, eta, dx_rw, dx_ref, dt_rw, twrite, ctres, ctres2
    double precision :: ic1, ic2, ic3 ! Additional parameters for initial conds.
    double precision :: pi = 3.1415926536d0
    double precision :: DUMMY = 1.d0 ! Parameters for Marsaglia's RNG
    character(100) :: path, scenario
    
contains

    subroutine init_params()
     
        open(10, file='parameters')
        read(10,*) nx
        read(10,*) xmax
        read(10,*) ymax
        read(10,*)
        read(10,*) scenario
        read(10,*) ic1
        read(10,*) ic2
        read(10,*) ic3
        read(10,*) nft
        read(10,*) lat0
        read(10,*) radius
        read(10,*) ctres
        read(10,*) sigma
        read(10,*) alpha
        read(10,*) beta
        read(10,*) ifac
        read(10,*)
        read(10,*) tmax
        read(10,*) twrite
        read(10,*)
        read(10,*) eta
        read(10,*) dx_rw
        read(10,*) dtopt
        read(10,*) dx_ref
        read(10,*) sim_rwc
        read(10,*) rwcc ! Random Walk Constrained Cancellation
        read(10,*)
        read(10,*) inc_shear
        read(10,*) inc_mflow
        read(10,*) inc_inflows
        read(10,*) inc_rw
        read(10,*)
        read(10,*) path
        read(10,*) verbose
        read(10,*)
        read(10,*) rseed
        close(10)
        
        ! Initialization of some other parameters
        dt_rw = dx_rw**2/eta/1.d6/4.d0
        if (twrite == 0.e0) twrite = dt_rw
        xmin = -xmax
        ymin = -ymax
        ny = nint(ymax/xmax*nx)
        ctres2 = ctres**2
        idum = -rseed ! Initialization for NR ran2 subroutine must be negative
        call ZBQLINI(rseed)
        
        if (verbose) then       
            print *
            print *, 'STENCIL CODE'
            print *, '------------'
            print *
            print *, 'x-Size of grid (points): ', nx
            print *, 'y-Size of grid (points): ', ny
            print *, 'xmin = ', xmin
            print *, 'xmax = ', xmax
            print *, 'ymin = ', ymin
            print *, 'ymax = ', ymax
            print *, 'Scenario: '//trim(scenario)
            print *, 'I. conds. parameter #1: ', ic1
            print *, 'I. conds. parameter #2: ', ic2
            print *, 'I. conds. parameter #3: ', ic3
            print *, 'Number of features: ', nft
            print *, 'AR latitude (deg): ', lat0
            print *, 'AR radius (m): ', radius
            print *, 'Cancel threshold: ', ctres
            print *, 'sigma (m) = ', sigma
            print *, 'alpha (m/s) = ', alpha
            print *, 'beta = ', beta
            print *, 'Inflow mult. factor = ',ifac
            print *, 'Max. time: ', tmax
            print *, 'Writing interval: ', twrite
            print *, 'eta (km2/s) = ', eta
            print *, 'RW step size (m) = ', dx_rw
            print *, 'Time step optimization?', dtopt
            print *, 'Reference stepsize: ', dx_ref
            print *, 'Simultaneous RW change? ', sim_rwc
            print *, 'Random walk constrained cancellation?', rwcc
            print *, 'Include shear? ', inc_shear
            print *, 'Include meridional flow? ', inc_mflow
            print *, 'Include inflows? ', inc_inflows
            print *, 'Include random walk? ', inc_rw
            print *, 'Random seed: ', rseed
            print *, 'Output directory: '//trim(path)
            print *, 'Verbose output?', verbose        
        end if
        
        ! Convert latitude to radians
        lat0 = lat0*pi/180.d0
        
    end subroutine init_params   
    
end module params
