! ******************************************************************************
! ****                          MODULE IOBELDA                            ****
! ******************************************************************************
!
! Several I/O functions and subroutines
 
module iobelda

contains

! ---- SUBROUTINE int2str ----!
! Returns a string of an integer number right-justified with 0's
! INPUT
! -- n -> integer number
! -- l -> string length
! OUTPUT
! -- filenum -> string
!
! Example: n=32, l=5 => filenum will contain the characters '00032'

    subroutine int2str(n,l,filenum)
    
        implicit none
        integer :: n, nc, l, i, q
        character(l) :: filenum
        
        nc = n ! I don't want to modify n
                
        do i = l-1,0,-1
            q = nc/10**i
            filenum(l-i:l-i) = char(48+q)
            nc = nc - q*10**i
        end do
        
    end subroutine int2str
    
end module iobelda
