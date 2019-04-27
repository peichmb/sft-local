program stencil

    use params
    use grid
    use ustencil
    use lsflows
    use vars
    use iobelda, only: int2str
    use update
    use initconds
    use timestep
    use random
    implicit none
    
    logical :: end_loop=.false., wrfile=.false.
    logical :: ss, sm
    integer :: i, nsteps=0, nwrite=1, fnamelen=50, fnumlen=4, nrw
    double precision :: t=0.d0, dt0, dt, dt20, tac=0.d0, miv
    character(4) :: fnum
    character(100) :: fname
        
    ! Initialize
    call init()

    ! Maximum inflow velocity value
    if (inc_inflows) then
        call update_umaps()
        miv = maxval(sqrt(ux**2+uy**2)) 
    end if
    print *
    print *,'Maximum inflow velocity: ',miv
    
    ! Calculate stepsize
    call calc_dt(dt0, nrw)
    dt = dt0
    
    ! Deliver tokens for non simultaneous change of RW velocities
    if (.not. sim_rwc) call tokens(nrw)
    
    ! Write initial snapshot
    fname=trim(path)//'/shot0000'
    print *
    call wrout(trim(fname), len(trim(fname)))
    
    ! ---------
    ! Main loop
    ! ---------
    
    do while (t<tmax)
        
        ! Check if it is writing time
        if (tac+dt0 >= twrite) then
            dt = twrite-tac
            wrfile = .true.
        end if
        
        ! Cancel out nearby opposite polarity elements
        if (rwcc) then
            if ( mod(nsteps, nrw) == 0 ) then
!                print *,'CANCELLATION NOW ', t
                call cancel()
            end if
        else
            call cancel()
        end if
        
        ! Update positions and pos. indices
        call update_positions(dt, end_loop)
        if (end_loop) exit
        
        ! Update velocity maps
        if (inc_inflows) call update_umaps()
        
        ! Update random walk velocities
        if (inc_rw) call update_rwu(nsteps, nrw)
        
        ! Write if it is time
        if (wrfile) then
            call int2str(nwrite, fnumlen, fnum)
            fname=trim(path)//'/shot'//fnum
            call wrout(trim(fname), len(trim(fname)))
            dt = dt0
            tac = 0.d0
            wrfile = .false.
            nwrite = nwrite+1
        end if
        
        ! Update time and step number
        t = t+dt
        tac = tac+dt
        nsteps = nsteps + 1
        
    end do
        
    ! ---------------
    ! Main loop - end
    ! ---------------    

    print *
    print *, 't_final = ',t
    print *, 'number of steps = ',nsteps
    if (end_loop) then
        print *, 'At least one particle escaped.'
    else
        print *, 'No elements escaped.'
    end if
    
contains

    subroutine init()
        
        implicit none
        type(feature), pointer :: curr
        double precision :: r
        
        call init_params()
        call init_grid()
        call init_stencil()        
        call init_vars()
        call init_initconds()       
                
    end subroutine init
    
    subroutine wrout(fname, fnamelen)
    
        implicit none
        
        integer :: fnamelen, i, j
        double precision :: ssr2
        character(fnamelen) :: fname
        type(feature), pointer :: curr
        
        curr => feature_head
        
        if (verbose) print *,'Writing output at time t = ',t
        open(10,file=fname)
        write(10,*) t
        do while (associated(curr%next))
            curr => curr%next
            ssr2 = ( curr%ssx**2+curr%ssy**2 )/1.d6 ! en km2
            write(10,'(1i6,2f20.6,1i4,3f20.6)') curr%id, curr%x, curr%y, curr%pol, curr%ssx, curr%ssy, ssr2/4.d0/t
        end do
        close(10)
                
    end subroutine wrout
    
end program stencil



