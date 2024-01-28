program test

use constants, only : dp, pi, pio2, twpi, p_zp

use mt_serial,  only :  sgrnd,  &
                        grnd

use coordinate, only :  coo_rot_v,      &
                        coo_rot_cv

implicit none

integer :: ii, ij, k, n
integer, dimension(3) :: iarr1, iarr2
integer, dimension(100) :: ilist
real(kind=dp) :: ri, rj, rk, rl, rm, rn
real(kind=dp), dimension(:,:), allocatable :: coo
real(kind=dp), dimension(100) :: rarr
character(LEN=15),dimension(3) :: carr

allocate (coo(3,10))

coo(:,1) = [1d-16, 1d0, -1d0]
coo(:,2) = [-1d-15, 0d0, 1d0]

write (*,*) coo(:,1), coo(:,2)

call coo_rot_v(2, p_zp, coo(:,1:2), pio2)

write (*,*) coo(:,1), coo(:,2)




end program
