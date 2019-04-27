! ******************************************************************************
! ****                          MODULE INITCONDS                            ****
! ******************************************************************************
!
! This module sets up the linked list containing the initial conditions of
!   each magnetic feature.

module initconds
        
    use params
    use grid    
    use vars
    use random
    implicit none
    
    INTERFACE 
      FUNCTION ZBQLU01(DUMMY)
        DOUBLE PRECISION :: ZBQLU01
        DOUBLE PRECISION, INTENT(IN) :: DUMMY
      END FUNCTION ZBQLU01
   END INTERFACE
   
   INTERFACE 
      FUNCTION ZBQLNOR(MU,SIGMA)
        DOUBLE PRECISION :: ZBQLNOR
        DOUBLE PRECISION, INTENT(IN) :: MU, SIGMA
      END FUNCTION ZBQLNOR
    END INTERFACE 
   
contains

    subroutine init_initconds()
        
        integer :: i
        double precision :: r1, r2, alpha, beta, w
        type(feature), pointer :: temp, curr
    
        if (trim(scenario) == 'monopolar') then    
            ! Build linked list according to initial scenario
            curr => feature_head
            do i=1,nft
                r1 = ZBQLU01(DUMMY)
                r2 = ZBQLU01(DUMMY)
                allocate(temp)
                temp%x = radius*sqrt(r2)*cos(2.d0*pi*r1)           
                temp%y = radius*sqrt(r2)*sin(2.d0*pi*r1)
                temp%xi = int((temp%x-xmin+dx/2.)/dx)
                temp%yi = int((temp%y-ymin+dy/2.)/dy)
                temp%pol = 1
                temp%ssx = 0.d0
                temp%ssy = 0.d0
                ! Initial random walk velocities
                if (inc_rw) then
                    r1 = ZBQLU01(DUMMY)
                    temp%urwx = dx_rw*cos(2.d0*pi*r1)/dt_rw
                    temp%urwy = dx_rw*sin(2.d0*pi*r1)/dt_rw
                else
                    temp%urwx = 0.d0
                    temp%urwy = 0.d0
                end if                
                nullify(temp%next)
                curr%next => temp
                curr => temp
            end do
        
        else if (trim(scenario) == 'monogauss') then
            
            ! Build linked list according to initial scenario
            curr => feature_head
            do i=1,nft
                allocate(temp)
                ! Box-Muller algorithm for Gaussian RN dist.
                !do
                !    r1 = 2.d0*ran2(idum)-1.d0
                !    r2 = 2.d0*ran2(idum)-1.d0
                !    w = r1*r1 + r2*r2
                !    if (w < 1.d0) exit
                !end do
                !w = sqrt((-2.d0*log(w))/w)
                !temp%x = r1*w*radius
                !temp%y = r2*w*radius                
                temp%x = ZBQLNOR(0.d0,radius)
                temp%y = ZBQLNOR(0.d0,radius)
                temp%xi = int((temp%x-xmin+dx/2.)/dx)
                temp%yi = int((temp%y-ymin+dy/2.)/dy)
                temp%pol = 1
                ! Initial random walk velocities
                if (inc_rw) then
                    r1 = ran2(idum)
                    temp%urwx = dx_rw*cos(2.d0*pi*r1)/dt_rw
                    temp%urwy = dx_rw*sin(2.d0*pi*r1)/dt_rw
                else
                    temp%urwx = 0.d0
                    temp%urwy = 0.d0
                end if                
                nullify(temp%next)
                curr%next => temp
                curr => temp
            end do        
        
        
        
                    
        else if (trim(scenario) == 'bipolar') then    
            ! Build linked list according to initial scenario
            curr => feature_head
            beta = pi/2.d0 - ic1/180.d0*pi ! tilt angle (0 - 180)
            do i=1,nft
                r1 = ran2(idum)
                r2 = ran2(idum)
                alpha = 2.d0*pi*r1
                allocate(temp)
                temp%x = radius*sqrt(r2)*cos(alpha)           
                temp%y = radius*sqrt(r2)*sin(alpha)
                temp%xi = int((temp%x-xmin+dx/2.)/dx)
                temp%yi = int((temp%y-ymin+dy/2.)/dy)
                if (alpha >= beta .and. alpha < beta+pi) then
                    temp%pol = 1
                else
                    temp%pol = -1
                end if
                ! Initial random walk velocities
                if (inc_rw) then
                    r1 = ran2(idum)
                    temp%urwx = dx_rw*cos(2.d0*pi*r1)/dt_rw
                    temp%urwy = dx_rw*sin(2.d0*pi*r1)/dt_rw
                else
                    temp%urwx = 0.d0
                    temp%urwy = 0.d0
                end if                
                nullify(temp%next)
                curr%next => temp
                curr => temp
            end do
        else if (trim(scenario) == 'bipolar2') then    
            ! Build linked list according to initial scenario
            curr => feature_head
            beta = ic1/180.d0*pi ! tilt angle (0 - 180)
            do i=1,nft/2 ! Left spot
                r1 = ZBQLU01(DUMMY)
                r2 = ZBQLU01(DUMMY)
                alpha = 2.d0*pi*r1
                allocate(temp)
                temp%x = 0.5d0*radius*(sqrt(r2)*cos(alpha)-cos(beta))
                temp%y = 0.5d0*radius*(sqrt(r2)*sin(alpha)+sin(beta))
                temp%xi = int((temp%x-xmin+dx/2.)/dx)
                temp%yi = int((temp%y-ymin+dy/2.)/dy)
                temp%pol = 1
                temp%id = i
                temp%ssx = 0.d0
                temp%ssy = 0.d0
                ! Initial random walk velocities
!                if (inc_rw) then
!                    r1 = ran2(idum)
!                    temp%urwx = dx_rw*cos(2.d0*pi*r1)/dt_rw
!                    temp%urwy = dx_rw*sin(2.d0*pi*r1)/dt_rw
!                else
                    temp%urwx = 0.d0
                    temp%urwy = 0.d0
!                end if                
                nullify(temp%next)
                curr%next => temp
                curr => temp
            end do
            do i=1,nft/2 ! Right spot
                r1 = ZBQLU01(DUMMY)
                r2 = ZBQLU01(DUMMY)
                alpha = 2.d0*pi*r1
                allocate(temp)
                temp%x = 0.5d0*radius*(sqrt(r2)*cos(alpha)+cos(beta))
                temp%y = 0.5d0*radius*(sqrt(r2)*sin(alpha)-sin(beta))
                temp%xi = int((temp%x-xmin+dx/2.)/dx)
                temp%yi = int((temp%y-ymin+dy/2.)/dy)
                temp%pol = -1
                temp%id = i+nft/2
                ! Initial random walk velocities
!                if (inc_rw) then
!                    r1 = ran2(idum)
!                    temp%urwx = dx_rw*cos(2.d0*pi*r1)/dt_rw
!                    temp%urwy = dx_rw*sin(2.d0*pi*r1)/dt_rw
!                else
                    temp%urwx = 0.d0
                    temp%urwy = 0.d0
!                end if                
                nullify(temp%next)
                curr%next => temp
                curr => temp
            end do
        else if (trim(scenario) == 'ctest') then    
            ! Cancel test
            curr => feature_head
            nft = 4
            do i=1,nft
                allocate(temp)
                if (i == 1) then
                    temp%x = -50.d6
                    temp%y = -0.45d6
                    temp%xi = int((temp%x-xmin+dx/2.)/dx)
                    temp%yi = int((temp%y-ymin+dy/2.)/dy)
                else if (i == 2) then
                    temp%x = 50.1d6
                    temp%y = 0.45d6
                    temp%xi = int((temp%x-xmin+dx/2.)/dx)
                    temp%yi = int((temp%y-ymin+dy/2.)/dy)
                else if (i == 3) then
                    temp%x = 50.d6
                    temp%y = 0.45d6
                    temp%xi = int((temp%x-xmin+dx/2.)/dx)
                    temp%yi = int((temp%y-ymin+dy/2.)/dy)
                else
                    temp%x = 0.d6
                    temp%y = 100.4d6
                    temp%xi = int((temp%x-xmin+dx/2.)/dx)
                    temp%yi = int((temp%y-ymin+dy/2.)/dy)
                end if
                
                if (temp%x >= 0.) then
                    temp%pol = 1
                else
                    temp%pol = -1
                end if
                ! Initial random walk velocities
                if (inc_rw) then
                    r1 = ran2(idum)
                    temp%urwx = dx_rw*cos(2.d0*pi*r1)/dt_rw
                    temp%urwy = dx_rw*sin(2.d0*pi*r1)/dt_rw
                else
                    temp%urwx = 0.d0
                    temp%urwy = 0.d0
                end if                
                nullify(temp%next)
                curr%next => temp
                curr => temp
                print *, i
            end do
        else
            print *, 'Initial scenario not implemented.'
            stop
        end if
                
    end subroutine init_initconds
    
    subroutine tokens(nrw)
    
        integer :: itk, nrw
        double precision :: tk
        type(feature), pointer :: curr
    
        curr => feature_head
        do while (associated(curr%next))
            curr => curr%next
            tk = ran2(idum)
            itk = int(tk*float(nrw))
            curr%rwtoken = itk+1
        end do
    
    end subroutine tokens
    
end module initconds
