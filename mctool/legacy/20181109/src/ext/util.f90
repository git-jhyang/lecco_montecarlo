module util

use constants,  only :  dp, pi, twpi, pio2, &
                        p_xp, p_zp

use iotool,     only :  io_printerr

implicit none

interface assert_data
   module procedure     assert_data_i0, assert_data_i1, &
                        assert_data_i2, assert_data_i3, &
                        assert_data_r0, assert_data_r1, &
                        assert_data_r2, assert_data_r3, &
                        assert_data_d0, assert_data_d1, &
                        assert_data_d2, assert_data_d3, &
                        assert_data_l0, assert_data_l1, &
                        assert_data_l2, assert_data_l3, &
                        assert_data_c0, assert_data_c1
end interface
 
contains
   
 !-----------------------------------------------------------------------------!
   
 function angle(vec1, vec2) result (ang)
   
   ! calculate angle between two vectors, vec1 and vec2, using intrinsic acos
   ! function
   
   implicit none
   
   real(kind=dp), dimension(3), intent(in) :: vec1, vec2  ! input two vectors
   real(kind=dp)                           :: ang
   
   real(kind=dp) :: br, len1, len2

   real(kind=dp), dimension(3,2) :: vecs

   integer :: i

   len1 = sqrt(dot_product(vec1, vec1))       ! length of vector
   len2 = sqrt(dot_product(vec2, vec2))

   if (len1 .ne. 0d0) then           ! change to unit vector
      vecs(:,1) = vec1(:)/len1
   else
      ang = 0d0
      return
   endif

   if (len2 .ne. 0d0) then           ! change to unit vector
      vecs(:,2) = vec2(:)/len2
   else
      ang = 0d0
      return
   endif

   br = dot_product(vecs(:,1), vecs(:,2))  ! dot product of vectors

   !write (*,*) br
   ang = acos(br)
!   if (abs(br) .lt.  1.0d-16) angle = pio2       
                ! if dot_product is too small, angle is 90
!   if (    br  .gt.  1.0d0 - 1.0d-16) angle = 0d0          
                ! if dot_product is near 1, angle is 0
!   if (    br  .lt. -1.0d0 + 1.0d-16) angle = pi

 end function angle

!-----------------------------------------------------------------------------!

 function gaussian(sigma, x, u) result(val)

   ! calculate gaussian value at x
   ! sigma : width of function
   ! u : center of function
   
   implicit none

   real(kind=dp), intent(in) :: sigma, x, u
   real(kind=dp) :: val
   
   val = 1.0d0/(sigma*sqrt(2*pi))*exp(-5d-1*((x-u)/sigma)**2)
   
 end function gaussian

!-----------------------------------------------------------------------------!
   
 function cross_product(vec1, vec2) result(vecout)
   
   implicit none
   
   real(kind=dp), dimension(3), intent(in) :: vec1, vec2
   real(kind=dp), dimension(3)             :: vecout
   
   vecout(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
   vecout(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
   vecout(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
   
 end function cross_product
 
!-----------------------------------------------------------------------------!

 subroutine rot_x(c1, c2, theta) 
   
   ! rotate c1 vector with theta around x axis, and gives c2 vector
   ! c1 : input vector
   ! c2 : output vector
   ! theta : rotaion angle
   
   implicit none
   
   real(kind=dp), dimension(3), intent(in)  :: c1
   real(kind=dp), dimension(3), intent(out) :: c2
   real(kind=dp),               intent(in)  :: theta
   
   real(kind=dp) :: rmat(3,3)
   integer :: i, j
   
   rmat(:,:) = 0d0      ! rotational matrix
   c2(:) = 0d0
   
   rmat(2,2) = cos(theta)
   rmat(2,3) = -sin(theta)
   rmat(3,2) = sin(theta)
   rmat(3,3) = cos(theta)
   
   rmat(1,1) = 1.0d0
   
   do i = 1, 3
   do j = 1, 3
      if (abs(rmat(i,j)) .lt. 1.0d-14) rmat(i,j) = 0d0
      c2(i) = c2(i) + rmat(i,j)*c1(j)
   enddo
   enddo
   
 end subroutine rot_x
   
!-----------------------------------------------------------------------------!

 subroutine rot_y(c1, c2, theta)
   
   ! same as rot_x
   implicit none
   
   real(kind=dp), dimension(3), intent(in)  :: c1
   real(kind=dp), dimension(3), intent(out) :: c2
   real(kind=dp),               intent(in)  :: theta
   
   real(kind=dp) :: rmat(3,3)
   integer :: i, j
   
   rmat(:,:) = 0d0
   c2(:) = 0d0
   
   rmat(1,1) = cos(theta)
   rmat(1,3) = sin(theta)
   rmat(3,1) = -sin(theta)
   rmat(3,3) = cos(theta)
   
   rmat(2,2) = 1.0d0
   
   do i = 1, 3
   do j = 1, 3
      if (abs(rmat(i,j)) .lt. 1.0d-14) rmat(i,j) = 0d0
      c2(i) = c2(i) + rmat(i,j)*c1(j)
   enddo
   enddo
   
 end subroutine rot_y
   
!-----------------------------------------------------------------------------!

 subroutine rot_z(c1, c2, theta)
   
   ! same as rot_x
   implicit none
   
   real(kind=dp), dimension(3), intent(in)  :: c1
   real(kind=dp), dimension(3), intent(out) :: c2
   real(kind=dp),               intent(in)  :: theta
   
   real(kind=dp) :: rmat(3,3)
   integer :: i, j
   
   rmat(:,:) = 0d0
   c2(:) = 0d0
   
   rmat(1,1) = cos(theta)
   rmat(1,2) = -sin(theta)
   rmat(2,1) = sin(theta)
   rmat(2,2) = cos(theta)
   rmat(3,3) = 1.0d0
   
   do i = 1, 3
   do j = 1, 3
      if (abs(rmat(i,j)) .lt. 1.0d-14) rmat(i,j) = 0d0
      c2(i) = c2(i) + rmat(i,j)*c1(j)
   enddo
   enddo
   
 end subroutine rot_z

!-----------------------------------------------------------------------------!
   
 subroutine rot_vec(pnt1, pnt2, vec, th)
   
   ! rotate 'vec' by 'th' around the vector goes through from
   ! 'pnt1' to 'pnt2'.
   
   implicit none
   
   real(kind=dp), dimension(3), intent(in)    :: pnt1, pnt2
   real(kind=dp), dimension(3), intent(inout) :: vec
   real(kind=dp),               intent(in)    :: th

   real(kind=dp) :: a, b, c, r, u, v, w
   real(kind=dp), dimension(3)   :: vecout
   real(kind=dp), dimension(3,4) :: rmat
   integer :: ii, ij
   
   a = pnt1(1)             ! origin point
   b = pnt1(2)
   c = pnt1(3)
   
   u = pnt2(1) - pnt1(1)   ! direction of vector
   v = pnt2(2) - pnt1(2)
   w = pnt2(3) - pnt1(3)
   
   r = sqrt(u*u + v*v + w*w)       ! vector length
   
   u = u/r                 ! normalization
   v = v/r
   w = w/r
   
   ! Source : http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
   
   rmat(1,1) = u*u + (v*v + w*w)*cos(th)
   rmat(1,2) = u*v*(1.0d0-cos(th)) - w*sin(th)
   rmat(1,3) = u*w*(1.0d0-cos(th)) + v*sin(th)
   rmat(1,4) = (a*(v*v + w*w) - u*(b*v + c*w))*(1.0d0-cos(th)) + (b*w - c*v)*sin(th)
   
   rmat(2,1) = u*v*(1.0d0-cos(th)) + w*sin(th)
   rmat(2,2) = v*v + (u*u + w*w)*cos(th)
   rmat(2,3) = v*w*(1.0d0-cos(th)) - u*sin(th)
   rmat(2,4) = (b*(u*u + w*w) - v*(a*u + c*w))*(1.0d0-cos(th)) + (c*u - a*w)*sin(th)
   
   rmat(3,1) = u*w*(1.0d0-cos(th)) - v*sin(th)
   rmat(3,2) = v*w*(1.0d0-cos(th)) + u*sin(th)
   rmat(3,3) = w*w + (u*u + v*v)*cos(th)
   rmat(3,4) = (c*(u*u + v*v) - w*(a*u + b*v))*(1.0d0-cos(th)) + (a*v - b*u)*sin(th)
   
   vecout(:) = 0d0
   do ii = 1, 3
      do ij = 1, 3
         vecout(ii) = vecout(ii) + rmat(ii,ij)*vec(ij)
      enddo
      vecout(ii) = vecout(ii) + rmat(ii,4)
   enddo
   
   vec = vecout
   
 end subroutine rot_vec
   
!-----------------------------------------------------------------------------!
   
 function gen_rseed() result(seed)
   
   ! generate 1 random integer based on cluster time
   
   implicit none
   
   integer :: seed, ii
   integer, dimension(8) :: date
   
   call date_and_time(values = date)
   
   ! write (*,*) date
   seed = 1
   do ii = 1, 6
      seed = seed + date(ii)
   enddo
   
   do ii = 7, 8
      if (mod(date(ii), 2) .eq. 0) then
         seed = seed + date(ii) 
      else 
         seed = seed*date(ii)
      endif
   enddo
   
 end function  
   
!-----------------------------------------------------------------------------!

 subroutine dir_to_cart(LattVec, coo)     ! fractional to cartesian
   implicit none
   
   real(kind=dp), dimension(3,3), intent(in)    :: LattVec ! lattice vector
   real(kind=dp), dimension(3),   intent(inout) :: coo     ! xyz inform.
   
   real(kind=dp) :: vol
   real(kind=dp), dimension(3)   :: LenVec, ang, coo_out
   real(kind=dp), dimension(3,3) :: conv_mat    
   integer :: i
   
   do i = 1, 3
      LenVec(i) = sqrt(dot_product(LattVec(:,i),LattVec(:,i)))
   end do
   
   ! calculate lattice angle
   ang(1) = angle(LattVec(:,2), LattVec(:,3))    ! alpha
   ang(2) = angle(LattVec(:,3), LattVec(:,1))    ! beta
   ang(3) = angle(LattVec(:,1), LattVec(:,2))    ! gamma
   
   ! conversion matrix comes from following address
   ! http://www.ruppweb.org/Xray/tutorial/Coordinate%20system%20transformation.htm
   
   vol = lenvec(1)*lenvec(2)*lenvec(3)*sqrt(1.0d0 - cos(ang(1))**2 - &
   &     cos(ang(2))**2 - cos(ang(3))**2 + 2*cos(ang(1))*cos(ang(2))*cos(ang(3)))
   
   ! conv_mat is inverse of 'M'
   
   conv_mat(1,1) = lenvec(1)
   conv_mat(1,2) = lenvec(2)*cos(ang(3))
   conv_mat(1,3) = lenvec(3)*cos(ang(2))
   
   conv_mat(2,1) = 0d0
   conv_mat(2,2) = lenvec(2)*sin(ang(3))
   conv_mat(2,3) = lenvec(3)*(cos(ang(1)) - cos(ang(2))*cos(ang(3)))/sin(ang(3))
   
   conv_mat(3,1) = 0d0
   conv_mat(3,2) = 0d0
   conv_mat(3,3) = vol/(lenvec(1)*lenvec(2)*sin(ang(3)))
   
   do i = 1, 3
      do while (coo(i) .lt. 0d0)         ! if fractional is negative, 
         coo(i) = coo(i) + 1.0d0
      enddo
      do while (coo(i) .ge. 1.0d0)    ! if fractional is larger than 1, 
         coo(i) = coo(i) - 1.0d0
      enddo
   enddo
   
   coo_out(:) = 0d0
   do i = 1, 3     ! make cartesian
      coo_out(:) = coo_out(:) + conv_mat(:,i)*coo(i)
   enddo
   
   coo(:) = coo_out(:)
   
 end subroutine dir_to_cart
   
!-----------------------------------------------------------------------------!

 subroutine cart_to_dir(LattVec, coo)     ! cartesian to fractional
   implicit none
   
   real(kind=dp), dimension(3,3), intent(in)    :: LattVec
   real(kind=dp), dimension(3),   intent(inout) :: coo
   
   real(kind=dp) :: vol
   real(kind=dp), dimension(3)   :: lenvec, ang, coo_out
   real(kind=dp), dimension(3,3) :: conv_mat

   integer :: i
   
   do i = 1, 3
      lenvec(i) = sqrt(dot_product(LattVec(:,i), LattVec(:,i)))
   end do
   
   ! write (*,*) lenvec(:)
   
   ! calculate angles
   ang(1) = angle(LattVec(:,2), LattVec(:,3))    ! alpha
   ang(2) = angle(LattVec(:,3), LattVec(:,1))    ! beta
   ang(3) = angle(LattVec(:,1), LattVec(:,2))    ! gamma
   
   ! write (*,*) ang(:)
   
   ! matrix is 'M' from same reference shown above
   
   vol = lenvec(1)*lenvec(2)*lenvec(3)*sqrt(1.0d0 - cos(ang(1))**2 - &
         cos(ang(2))**2 - cos(ang(3))**2 + 2*cos(ang(1))*cos(ang(2))*cos(ang(3)))
   
   conv_mat(1,1) = 1.0d0/lenvec(1)
   conv_mat(1,2) = -cos(ang(3))/(lenvec(1)*sin(ang(3)))
   conv_mat(1,3) = lenvec(2)*lenvec(3)*(cos(ang(3))*cos(ang(1)) - & 
                   cos(ang(2)))/(sin(ang(3))*vol)
   
   conv_mat(2,1) = 0d0
   conv_mat(2,2) = 1.0d0/(lenvec(2)*sin(ang(3)))
   conv_mat(2,3) = lenvec(1)*lenvec(3)*(cos(ang(2))*cos(ang(3)) - &
                   cos(ang(1)))/(sin(ang(3))*vol)
   
   conv_mat(3,1) = 0d0
   conv_mat(3,2) = 0d0
   conv_mat(3,3) = lenvec(1)*lenvec(2)*sin(ang(3))/vol
   
   coo_out(:) = 0d0
   
   do i = 1, 3
      coo_out(:) = coo_out(:) + conv_mat(:,i)*coo(i)
   enddo
   
   coo = coo_out

   do i = 1, 3
      do while (coo(i) .lt. 0d0)
         coo(i) = coo(i) + 1.0d0
      enddo
      do while (coo(i) .ge. 1.0d0 )
         coo(i) = coo(i) - 1.0d0
      enddo
   enddo
   
 end subroutine cart_to_dir
   
!-----------------------------------------------------------------------------!

 subroutine cart_to_sph(co_c, co_s)
   
   ! convert cartesian coordinate to spherical coordinate
   implicit none
   
   real(kind=dp), dimension(3), intent(in)  :: co_c
   real(kind=dp), dimension(3), intent(out) :: co_s
   
   co_s(1) = sqrt(dot_product(co_c, co_c))
   co_s(2) = angle(co_c, p_zp)
   co_s(3) = angle((/co_c(1),co_c(2),0d0/), p_xp) 
   if (co_c(2) .lt. 0) co_s(3) = twpi - co_s(3)
   
 end subroutine cart_to_sph

!-----------------------------------------------------------------------------!
   
 subroutine sph_to_cart(co_s, co_c)
   ! convert spherical coordinate to cartesian coordinate
   implicit none
   
   real(kind=dp), dimension(3), intent(in)  :: co_s
   real(kind=dp), dimension(3), intent(out) :: co_c
   
   co_c(1) = co_s(1)*sin(co_s(2))*cos(co_s(3))
   co_c(2) = co_s(1)*sin(co_s(2))*sin(co_s(3))
   co_c(3) = co_s(1)*cos(co_s(2))
   
 end subroutine sph_to_cart
   
!-----------------------------------------------------------------------------!

 subroutine latt_refin(LattVec)
   implicit none
   
   real(kind=dp), dimension(3,3), intent(inout) :: LattVec
   
   integer :: ii, ij

   real(kind=dp), dimension(3) :: lenvec, ang
   
   lenvec(:) = 0d0

   do ii = 1, 3
      do ij = 1, 3
         lenvec(ii) = lenvec(ii) + LattVec(ij,ii)**2
      enddo
      lenvec(ii) = sqrt(lenvec(ii))
   enddo
   
   ang(1) = angle(LattVec(:,2), LattVec(:,3))    ! alpha
   ang(2) = angle(LattVec(:,3), LattVec(:,1))    ! beta
   ang(3) = angle(LattVec(:,1), LattVec(:,2))    ! gamma
   
   LattVec(1,1) = lenvec(1)
   LattVec(2,1) = 0d0
   LattVec(3,1) = 0d0   
   
   LattVec(1,2) = lenvec(2)*cos(ang(3))
   LattVec(2,2) = sqrt(lenvec(2)**2 - LattVec(1,2)**2)
   LattVec(3,2) = 0d0   
   
   LattVec(1,3) = lenvec(3)*cos(ang(2))
   LattVec(2,3) = (cos(ang(1))*lenvec(2)*lenvec(3) - LattVec(1,2)*LattVec(1,3))/LattVec(2,2)
   LattVec(3,3) = sqrt(lenvec(3)**2 - LattVec(1,3)**2 - LattVec(2,3)**2)
   
   do ii = 1, 3
      do ij = 1, 3
         if (abs(LattVec(ij,ii)) .lt. 1.0d-14) LattVec(ij,ii) = 0d0
      enddo
   enddo
   
 end subroutine latt_refin

!-----------------------------------------------------------------------------!

 subroutine time_formatting(time, tForm)
   implicit none
   real(kind=dp),    intent(in)  :: time
   character(LEN=*), intent(out) :: tForm
   
   integer :: t_h, t_m, t_s, t_ms

   if (len(tForm) .lt. 14) then
      call io_printerr("internal error - time_formatting","input string is short")
   endif

   t_h  = floor(time/3600)
   t_m  = floor(time/60 - t_h*60)
   t_s  = floor(time - t_m*60 - t_h*3600)
   t_ms = floor((time - floor(time))*100)
   if (t_h .ne. 0) then
      write (tForm,'(I4,"h",I2.2,"m",I2.2,".",I2.2,"s")') t_h, t_m, t_s, t_ms
   elseif (t_m .ne. 0) then
      write (tForm,'(5x,I2,"m",I2.2,".",I2.2,"s")') t_m, t_s, t_ms
   else
      write (tForm,'(8x,I2,".",I2.2,"s")') t_s, t_ms
   endif

 end subroutine time_formatting

!-----------------------------------------------------------------------------!

 function assert_data_i0(dat1, dat2) result (isErr)
   implicit none
   integer, intent(in) :: dat1, dat2
   logical             :: isErr

   isErr = .false.
   if (dat1 .ne. dat2) isErr = .True.

 end function assert_data_i0

 function assert_data_i1(dat1, dat2, a) result (isErr)
   implicit none
   integer,               intent(in) :: a
   integer, dimension(a), intent(in) :: dat1, dat2
   logical                           :: isErr

   integer :: ia

   isErr = .false.
   do ia = 1, a
      if (dat1(a) .ne. dat2(a)) then
         isErr = .True.
         return
      endif
   enddo

 end function assert_data_i1

 function assert_data_i2(dat1, dat2, a, b) result (isErr)
   implicit none
   integer,                 intent(in) :: a, b
   integer, dimension(a,b), intent(in) :: dat1, dat2
   logical                             :: isErr

   integer :: ia, ib

   isErr = .false.
   do ib = 1, b
   do ia = 1, a
      if (dat1(ia,ib) .ne. dat2(ia,ib)) then
         isErr = .True.
         return
      endif
   enddo
   enddo

 end function assert_data_i2

 function assert_data_i3(dat1, dat2, a, b, c) result (isErr)
   implicit none
   integer,                   intent(in) :: a, b, c
   integer, dimension(a,b,c), intent(in) :: dat1, dat2
   logical                               :: isErr

   integer :: ia, ib, ic

   isErr = .false.
   do ic = 1, c
   do ib = 1, b
   do ia = 1, a
      if (dat1(ia,ib,ic) .ne. dat2(ia,ib,ic)) then
         isErr = .True.
         return
      endif
   enddo
   enddo
   enddo

 end function assert_data_i3

!-----------------------------------------------------------------------------!

 function assert_data_r0(dat1, dat2) result (isErr)
   implicit none
   real,    intent(in) :: dat1, dat2
   logical             :: isErr

   isErr = .false.
   if (dat1 .ne. dat2) isErr = .True.

 end function assert_data_r0

 function assert_data_r1(dat1, dat2, a) result (isErr)
   implicit none
   integer,               intent(in) :: a
   real,    dimension(a), intent(in) :: dat1, dat2
   logical                           :: isErr

   integer :: ia

   isErr = .false.
   do ia = 1, a
      if (dat1(a) .ne. dat2(a)) then
         isErr = .True.
         return
      endif
   enddo

 end function assert_data_r1

 function assert_data_r2(dat1, dat2, a, b) result (isErr)
   implicit none
   integer,                 intent(in) :: a, b
   real,    dimension(a,b), intent(in) :: dat1, dat2
   logical                             :: isErr

   integer :: ia, ib

   isErr = .false.
   do ib = 1, b
   do ia = 1, a
      if (dat1(ia,ib) .ne. dat2(ia,ib)) then
         isErr = .True.
         return
      endif
   enddo
   enddo

 end function assert_data_r2

 function assert_data_r3(dat1, dat2, a, b, c) result (isErr)
   implicit none
   integer,                   intent(in) :: a, b, c
   real,    dimension(a,b,c), intent(in) :: dat1, dat2
   logical                               :: isErr

   integer :: ia, ib, ic

   isErr = .false.
   do ic = 1, c
   do ib = 1, b
   do ia = 1, a
      if (dat1(ia,ib,ic) .ne. dat2(ia,ib,ic)) then
         isErr = .True.
         return
      endif
   enddo
   enddo
   enddo

 end function assert_data_r3

!-----------------------------------------------------------------------------!

 function assert_data_d0(dat1, dat2) result (isErr)
   implicit none
   real(kind=dp), intent(in) :: dat1, dat2
   logical                   :: isErr

   isErr = .false.
   if (dat1 .ne. dat2) isErr = .True.

 end function assert_data_d0

 function assert_data_d1(dat1, dat2, a) result (isErr)
   implicit none
   integer,                     intent(in) :: a
   real(kind=dp), dimension(a), intent(in) :: dat1, dat2
   logical                                 :: isErr

   integer :: ia

   isErr = .false.
   do ia = 1, a
      if (dat1(a) .ne. dat2(a)) isErr = .True.
   enddo

 end function assert_data_d1

 function assert_data_d2(dat1, dat2, a, b) result (isErr)
   implicit none
   integer,                       intent(in) :: a, b
   real(kind=dp), dimension(a,b), intent(in) :: dat1, dat2
   logical                                   :: isErr

   integer :: ia, ib

   isErr = .false.
   do ib = 1, b
   do ia = 1, a
      if (dat1(ia,ib) .ne. dat2(ia,ib)) then
         isErr = .True.
         return
      endif
   enddo
   enddo

 end function assert_data_d2

 function assert_data_d3(dat1, dat2, a, b, c) result (isErr)
   implicit none
   integer,                         intent(in) :: a, b, c
   real(kind=dp), dimension(a,b,c), intent(in) :: dat1, dat2
   logical                                     :: isErr

   integer :: ia, ib, ic

   isErr = .false.
   do ic = 1, c
   do ib = 1, b
   do ia = 1, a
      if (dat1(ia,ib,ic) .ne. dat2(ia,ib,ic)) then
         isErr = .True.
         return
      endif
   enddo
   enddo
   enddo

 end function assert_data_d3

!-----------------------------------------------------------------------------!

 function assert_data_l0(dat1, dat2) result (isErr)
   implicit none
   logical, intent(in) :: dat1, dat2
   logical             :: isErr

   isErr = .false.
   if (dat1 .ne. dat2) isErr = .True.

 end function assert_data_l0

 function assert_data_l1(dat1, dat2, a) result (isErr)
   implicit none
   integer,               intent(in) :: a
   logical, dimension(a), intent(in) :: dat1, dat2
   logical                           :: isErr

   integer :: ia

   isErr = .false.
   do ia = 1, a
      if (dat1(a) .ne. dat2(a)) then
         isErr = .True.
         return
      endif
   enddo

 end function assert_data_l1

 function assert_data_l2(dat1, dat2, a, b) result (isErr)
   implicit none
   integer,                 intent(in) :: a, b
   logical, dimension(a,b), intent(in) :: dat1, dat2
   logical                             :: isErr

   integer :: ia, ib

   isErr = .false.
   do ib = 1, b
   do ia = 1, a
      if (dat1(ia,ib) .ne. dat2(ia,ib)) then
         isErr = .True.
         return
      endif
   enddo
   enddo

 end function assert_data_l2

 function assert_data_l3(dat1, dat2, a, b, c) result (isErr)
   implicit none
   integer,                   intent(in) :: a, b, c
   logical, dimension(a,b,c), intent(in) :: dat1, dat2
   logical                               :: isErr

   integer :: ia, ib, ic

   isErr = .false.
   do ic = 1, c
   do ib = 1, b
   do ia = 1, a
      if (dat1(ia,ib,ic) .ne. dat2(ia,ib,ic)) then
         isErr = .True.
         return
      endif
   enddo
   enddo
   enddo

 end function assert_data_l3

!-----------------------------------------------------------------------------!

 function assert_data_c0(dat1, dat2) result (isErr)
   implicit none
   character(LEN=*), intent(in) :: dat1, dat2
   logical                      :: isErr

   isErr = .false.
   if (dat1 .ne. dat2) isErr = .True.

 end function assert_data_c0

 function assert_data_c1(dat1, dat2, a) result (isErr)
   implicit none
   integer,                        intent(in) :: a
   character(LEN=*), dimension(a), intent(in) :: dat1, dat2
   logical                                    :: isErr

   integer :: ia

   isErr = .false.
   do ia = 1, a
      if (dat1(a) .ne. dat2(a)) then
         isErr = .True.
         return
      endif
   enddo

 end function assert_data_c1

 function assert_data_c2(dat1, dat2, a, b) result (isErr)
   implicit none
   integer,                          intent(in) :: a, b
   character(LEN=*), dimension(a,b), intent(in) :: dat1, dat2
   logical                                      :: isErr

   integer :: ia, ib

   isErr = .false.
   do ib = 1, b
   do ia = 1, a
      if (dat1(ia,ib) .ne. dat2(ia,ib)) then
         isErr = .True.
         return
      endif
   enddo
   enddo

 end function assert_data_c2

 function assert_data_c3(dat1, dat2, a, b, c) result (isErr)
   implicit none
   integer,                            intent(in) :: a, b, c
   character(LEN=*), dimension(a,b,c), intent(in) :: dat1, dat2
   logical                                        :: isErr

   integer :: ia, ib, ic

   isErr = .false.
   do ic = 1, c
   do ib = 1, b
   do ia = 1, a
      if (dat1(ia,ib,ic) .ne. dat2(ia,ib,ic)) then
         isErr = .True.
         return
      endif
   enddo
   enddo
   enddo

 end function assert_data_c3

end module
