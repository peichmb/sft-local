! ******************************************************************************
! ****                          MODULE TIMESTEP                             ****
! ******************************************************************************
!
! Contains the subroutine calc_dt
!
! PROCEDURES:
! -- calc_dt: calculates the stepsize. There are two possibilities:
! ---- Compute it from a reference random walk step. I.e., the timestep of
!    a random walk with the reference DX and a diffusion coefficient eta.
! ---- Compute it as one hundredth of the DT corresponding to DX = 20Mm and 
!    a diffusion coefficient eta.
module timestep

use params

contains

    subroutine calc_dt(dt0, nrw)
    
        double precision :: dt20, dt0
        integer nrw
        
        if (dtopt) then
            dt20 = 1.d6/eta ! Reference stepsize.
            nrw = nint(dt_rw/dt20)
            if (nrw == 0) nrw = 1
            dt0 = dt_rw/float(nrw)
        else
            dt0 = dx_ref**2/4.d0/1.d6/eta
            nrw = nint(dt_rw/dt0)
        end if

        print *
        print *, '** calc_dt **'
        print *, ' - Timestep sizes:'
        print *, ' - dt(rw) = ', dt_rw
        print *, ' - dt = ', dt0
        print *, ' - dt(rw)/dt = ', nrw
        
    end subroutine calc_dt
    
end module timestep
