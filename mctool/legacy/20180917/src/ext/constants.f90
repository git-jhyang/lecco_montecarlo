module constants
implicit none
save

integer,                     parameter :: dp      = kind(1d0)
integer,                     parameter :: VARLEN  = 64
!integer,                     parameter :: VARLEN  = 128
integer,                     parameter :: PATHLEN = 256
integer,                     parameter :: LINELEN = 256
integer,                     parameter :: CHKOK   = 0
integer,                     parameter :: CHKERR  = 1
integer,                     parameter :: CHKEOF  = -1
real(kind=dp), dimension(3), parameter :: p_o     = [ 0.0_dp,  0.0_dp,  0.0_dp]
real(kind=dp), dimension(3), parameter :: p_xp    = [ 1.0_dp,  0.0_dp,  0.0_dp]
real(kind=dp), dimension(3), parameter :: p_xn    = [-1.0_dp,  0.0_dp,  0.0_dp]
real(kind=dp), dimension(3), parameter :: p_yp    = [ 0.0_dp,  1.0_dp,  0.0_dp]
real(kind=dp), dimension(3), parameter :: p_yn    = [ 0.0_dp, -1.0_dp,  0.0_dp]
real(kind=dp), dimension(3), parameter :: p_zp    = [ 0.0_dp,  0.0_dp,  1.0_dp]
real(kind=dp), dimension(3), parameter :: p_zn    = [ 0.0_dp,  0.0_dp, -1.0_dp]
real(kind=dp),               parameter :: sqrt2   = 1.41421356237309504880_dp
real(kind=dp),               parameter :: euler   = 0.57721566490153286061_dp
real(kind=dp),               parameter :: pi      = 3.14159265358979323846_dp
real(kind=dp),               parameter :: pio2    = 1.57079632679489661923_dp
real(kind=dp),               parameter :: twpi    = 6.28318530717958647693_dp
real(kind=dp),               parameter :: kb_ev   = 8.6173303d-5
real(kind=dp),               parameter :: kb_j    = 1.38064852d-23
real(kind=dp),               parameter :: n_avo   = 6.022140857d+23
real(kind=dp),               parameter :: h_ev    = 4.135667662d-15
real(kind=dp),               parameter :: hb_ev   = 6.582119514d-16
real(kind=dp),               parameter :: q_at    = 1.6021766208d-19

end module constants
