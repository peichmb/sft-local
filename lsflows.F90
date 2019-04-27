module lsflows

    use params
    use grid
    implicit none
    
contains
       
    function ux_sh(y,lat0)
        
        double precision :: y, lat0, vxsh0, ux_sh
        double precision :: lat, rsun = 695500.d3, omega
        double precision :: ddtrs, omega0 !deg/day to rad/s
        
        lat = y/rsun + lat0 ! Latitudes are in radians

        ddtrs = pi/180.d0/3600.d0/24.d0

        omega0 = ddtrs*(13.38d0-2.30d0*sin(lat0)**2-1.62d0*sin(lat0)**4)
        omega =  ddtrs*(13.38d0-2.30d0*sin(lat)**2-1.62d0*sin(lat)**4)
        
        ux_sh = rsun*(omega-omega0)
        
    end function ux_sh
    
    function uy_mf(y,lat0)
    
        double precision :: y, lat0, lat, rsun = 695500.d3, uy_mf
        lat = y/rsun + lat0 ! Latitudes are in radians
        
        uy_mf = 11.d0*sin(2.4d0*lat)
    
    end function uy_mf
    
end module lsflows
