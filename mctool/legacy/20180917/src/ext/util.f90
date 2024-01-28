module util

use constants,  only :  dp, pi, twpi, pio2, &
                        p_xp, p_zp

implicit none
 
contains
   
 !-----------------------------------------------------------------------------!
   
 function angle(vec1, vec2)
   
   ! calculate angle between two vectors, vec1 and vec2, using intrinsic acos
   ! function
   
   implicit none
   
   real(kind=dp) :: vec1(3), vec2(3)        ! input two vectors
   real(kind=dp) :: angle
   
   real(kind=dp) :: br, len1, len2, vecs(3,2)

   integer :: i

   len1 = sqrt(dot_product(vec1, vec1))       ! length of vector
   len2 = sqrt(dot_product(vec2, vec2))

   if (len1 .ne. 0d0) then           ! change to unit vector
      vecs(:,1) = vec1(:)/len1
   else
      angle = 0d0
      return
   endif

   if (len2 .ne. 0d0) then           ! change to unit vector
      vecs(:,2) = vec2(:)/len2
   else
      angle = 0d0
      return
   endif

   br = dot_product(vecs(:,1), vecs(:,2))  ! dot product of vectors

   !write (*,*) br
   angle = acos(br)
!   if (abs(br) .lt.  1.0d-16) angle = pio2       
                ! if dot_product is too small, angle is 90
!   if (    br  .gt.  1d0 - 1.0d-16) angle = 0d0          
                ! if dot_product is near 1, angle is 0
!   if (    br  .lt. -1d0 + 1.0d-16) angle = pi

 end function angle

!-----------------------------------------------------------------------------!

 function gaussian(sigma, x, u)

   ! calculate gaussian value at x
   ! sigma : width of function
   ! u : center of function
   
   implicit none
   real(kind=dp), intent(in) :: sigma, x, u
   real(kind=dp) :: gaussian
   
   gaussian = 1d0/(sigma*sqrt(2*pi))*exp(-5d-1*((x-u)/sigma)**2)
   
 end function gaussian

!-----------------------------------------------------------------------------!

 function get_memuse() result(memuse)
   implicit none
   integer :: memuse
   
 end function
   
!-----------------------------------------------------------------------------!
 ! subroutines
   
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
   
   real(kind=dp), intent(in) :: c1(3), theta
   real(kind=dp), intent(out) :: c2(3)
   
   real(kind=dp) :: rmat(3,3)
   integer :: i, j
   
   rmat(:,:) = 0d0      ! rotational matrix
   c2(:) = 0d0
   
   rmat(2,2) = cos(theta)
   rmat(2,3) = -sin(theta)
   rmat(3,2) = sin(theta)
   rmat(3,3) = cos(theta)
   
   rmat(1,1) = 1d0
   
   do i = 1, 3
   do j = 1, 3
      if (abs(rmat(i,j)) .lt. 1d-14) rmat(i,j) = 0d0
      c2(i) = c2(i) + rmat(i,j)*c1(j)
   enddo
   enddo
   
 end subroutine rot_x
   
!-----------------------------------------------------------------------------!

 subroutine rot_y(c1, c2, theta)
   
   ! same as rot_x
   implicit none
   
   real(kind=dp), intent(in) :: c1(3), theta
   real(kind=dp), intent(out) :: c2(3)
   
   real(kind=dp) :: rmat(3,3)
   integer :: i, j
   
   rmat(:,:) = 0d0
   c2(:) = 0d0
   
   rmat(1,1) = cos(theta)
   rmat(1,3) = sin(theta)
   rmat(3,1) = -sin(theta)
   rmat(3,3) = cos(theta)
   
   rmat(2,2) = 1d0
   
   do i = 1, 3
   do j = 1, 3
      if (abs(rmat(i,j)) .lt. 1d-14) rmat(i,j) = 0d0
      c2(i) = c2(i) + rmat(i,j)*c1(j)
   enddo
   enddo
   
 end subroutine rot_y
   
!-----------------------------------------------------------------------------!

 subroutine rot_z(c1, c2, theta)
   
   ! same as rot_x
   implicit none
   
   real(kind=dp), intent(in) :: c1(3), theta
   real(kind=dp), intent(out) :: c2(3)
   
   real(kind=dp) :: rmat(3,3)
   integer :: i, j
   
   rmat(:,:) = 0d0
   c2(:) = 0d0
   
   rmat(1,1) = cos(theta)
   rmat(1,2) = -sin(theta)
   rmat(2,1) = sin(theta)
   rmat(2,2) = cos(theta)
   rmat(3,3) = 1d0
   
   do i = 1, 3
   do j = 1, 3
      if (abs(rmat(i,j)) .lt. 1d-14) rmat(i,j) = 0d0
      c2(i) = c2(i) + rmat(i,j)*c1(j)
   enddo
   enddo
   
 end subroutine rot_z

!-----------------------------------------------------------------------------!
   
 subroutine rot_vec(pnt1, pnt2, vec, th)
   
   ! rotate 'vec' by 'th' around the vector goes through from 'pnt1' to 'pnt2'.
   
   implicit none
   
   real(kind=dp), intent(in) :: pnt1(3), pnt2(3), th
   real(kind=dp), intent(inout) :: vec(3)
   
   real(kind=dp) :: a, b, c, r, u, v, w, vecout(3), rmat(3,4) 
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
   rmat(1,2) = u*v*(1d0-cos(th)) - w*sin(th)
   rmat(1,3) = u*w*(1d0-cos(th)) + v*sin(th)
   rmat(1,4) = (a*(v*v + w*w) - u*(b*v + c*w))*(1d0-cos(th)) + (b*w - c*v)*sin(th)
   
   rmat(2,1) = u*v*(1d0-cos(th)) + w*sin(th)
   rmat(2,2) = v*v + (u*u + w*w)*cos(th)
   rmat(2,3) = v*w*(1d0-cos(th)) - u*sin(th)
   rmat(2,4) = (b*(u*u + w*w) - v*(a*u + c*w))*(1d0-cos(th)) + (c*u - a*w)*sin(th)
   
   rmat(3,1) = u*w*(1d0-cos(th)) - v*sin(th)
   rmat(3,2) = v*w*(1d0-cos(th)) + u*sin(th)
   rmat(3,3) = w*w + (u*u + v*v)*cos(th)
   rmat(3,4) = (c*(u*u + v*v) - w*(a*u + b*v))*(1d0-cos(th)) + (a*v - b*u)*sin(th)
   
   vecout(:) = 0d0
   do ii = 1, 3
      do ij = 1, 3
         vecout(ii) = vecout(ii) + rmat(ii,ij)*vec(ij)
      enddo
      vecout(ii) = vecout(ii) + rmat(ii,4)
      if (abs(vecout(ii)) .lt. 1.0d-14) vecout(ii) = 0.0d0
   enddo
   
   vec = vecout
   
 end subroutine rot_vec
   
!-----------------------------------------------------------------------------!
   
 function gen_rseed() result(seed)
   
   ! generate 1 random integer based on cluster time
   
   implicit none
   
   integer :: seed
   integer :: date(8), ii
   
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

 subroutine dir_to_cart(lvec, coord)     ! fractional to cartesian
   implicit none
   
   real(kind=dp), intent(in) :: lvec(3,3)       ! lattice vector
   real(kind=dp), intent(inout) :: coord(3)     ! xyz inform.
   
   real(kind=dp) :: conv_mat(3,3), lenvec(3), ang(3), vol, coord_new(3)
   integer :: i
   
   do i = 1, 3
      lenvec(i) = sqrt(dot_product(lvec(:,i),lvec(:,i)))
   end do
   
   ! calculate lattice angle
   ang(1) = angle(lvec(:,2), lvec(:,3))    ! alpha
   ang(2) = angle(lvec(:,3), lvec(:,1))    ! beta
   ang(3) = angle(lvec(:,1), lvec(:,2))    ! gamma
   
   ! conversion matrix comes from following address
   ! http://www.ruppweb.org/Xray/tutorial/Coordinate%20system%20transformation.htm
   
   vol = lenvec(1)*lenvec(2)*lenvec(3)*sqrt(1d0 - cos(ang(1))**2 - &
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
   
   coord_new(:) = 0d0
   
   do i = 1, 3
      do while ( coord(i) .lt. 0 )         ! if fractional is negative, 
         coord(i) = coord(i) + 1d0
      enddo
      do while ( coord(i) .ge. 1d0 )    ! if fractional is larger than 1, 
         coord(i) = coord(i) - 1d0
      enddo
   enddo
   
   do i = 1, 3     ! make cartesian
      coord_new(1) = coord_new(1) + conv_mat(1,i)*coord(i)
      coord_new(2) = coord_new(2) + conv_mat(2,i)*coord(i)
      coord_new(3) = coord_new(3) + conv_mat(3,i)*coord(i)
   enddo
   
   coord(:) = coord_new(:)
   
 end subroutine dir_to_cart
   
!-----------------------------------------------------------------------------!

 subroutine cart_to_dir(lvec, coord)     ! cartesian to fractional
   implicit none
   
   real(kind=dp), intent(in) :: lvec(3,3)
   real(kind=dp), intent(inout) :: coord(3)
   
   real(kind=dp) :: conv_mat(3,3), lenvec(3), ang(3), vol, coord_new(3)
   integer :: i
   
   do i = 1, 3
      lenvec(i) = sqrt(dot_product(lvec(:,i), lvec(:,i)))
   end do
   
   ! write (*,*) lenvec(:)
   
   ! calculate angles
   ang(1) = angle(lvec(:,2), lvec(:,3))    ! alpha
   ang(2) = angle(lvec(:,3), lvec(:,1))    ! beta
   ang(3) = angle(lvec(:,1), lvec(:,2))    ! gamma
   
   ! write (*,*) ang(:)
   
   ! matrix is 'M' from same reference shown above
   
   vol = lenvec(1)*lenvec(2)*lenvec(3)*sqrt(1d0 - cos(ang(1))**2 - &
         cos(ang(2))**2 - cos(ang(3))**2 + 2*cos(ang(1))*cos(ang(2))*cos(ang(3)))
   
   conv_mat(1,1) = 1d0/lenvec(1)
   conv_mat(1,2) = -cos(ang(3))/(lenvec(1)*sin(ang(3)))
   conv_mat(1,3) = lenvec(2)*lenvec(3)*(cos(ang(3))*cos(ang(1)) - & 
                   cos(ang(2)))/(sin(ang(3))*vol)
   
   conv_mat(2,1) = 0d0
   conv_mat(2,2) = 1d0/(lenvec(2)*sin(ang(3)))
   conv_mat(2,3) = lenvec(1)*lenvec(3)*(cos(ang(2))*cos(ang(3)) - &
                   cos(ang(1)))/(sin(ang(3))*vol)
   
   conv_mat(3,1) = 0d0
   conv_mat(3,2) = 0d0
   conv_mat(3,3) = lenvec(1)*lenvec(2)*sin(ang(3))/vol
   
   coord_new(:) = 0d0
   
   do i = 1, 3
      coord_new(1) = coord_new(1) + conv_mat(1,i)*coord(i)
      coord_new(2) = coord_new(2) + conv_mat(2,i)*coord(i)
      coord_new(3) = coord_new(3) + conv_mat(3,i)*coord(i)
   enddo
   
   do i = 1, 3
      do while ( coord_new(i) .lt. 0 )
         coord_new(i) = coord_new(i) + 1d0
      enddo
      do while ( coord_new(i) .ge. 1d0 )
         coord_new(i) = coord_new(i) - 1d0
      enddo
   enddo
   
   coord(:) = coord_new(:)
   
 end subroutine cart_to_dir
   
!-----------------------------------------------------------------------------!

 subroutine cart_to_sph(co_c, co_s)
   
   ! convert cartesian coordinate to spherical coordinate
   implicit none
   
   real(kind=dp), intent(in) :: co_c(3)
   real(kind=dp), intent(out) :: co_s(3)
   
   co_s(1) = sqrt(dot_product(co_c, co_c))
   co_s(2) = angle(co_c, p_zp)
   co_s(3) = angle((/co_c(1),co_c(2),0d0/), p_xp) 
   if (co_c(2) .lt. 0) co_s(3) = twpi - co_s(3)
   
 end subroutine cart_to_sph

!-----------------------------------------------------------------------------!
   
 subroutine sph_to_cart(co_s, co_c)
   ! convert spherical coordinate to cartesian coordinate
   implicit none
   
   real(kind=dp), intent(in) :: co_s(3)
   real(kind=dp), intent(out) :: co_c(3)
   
   co_c(1) = co_s(1)*sin(co_s(2))*cos(co_s(3))
   co_c(2) = co_s(1)*sin(co_s(2))*sin(co_s(3))
   co_c(3) = co_s(1)*cos(co_s(2))
   
 end subroutine sph_to_cart
   
!-----------------------------------------------------------------------------!

 subroutine latt_refin(lvec)
   implicit none
   
   real(kind=dp), intent(inout) :: lvec(3,3)
   
   integer :: ii, ij
   real(kind=dp) :: lvec_tmp(3,3), lenvec(3), ang(3)
   
   lenvec(:) = 0d0
   do ii = 1, 3
      do ij = 1, 3
         lenvec(ii) = lenvec(ii) + lvec(ij,ii)**2
      enddo
      lenvec(ii) = sqrt(lenvec(ii))
   enddo
   
   ang(1) = angle(lvec(:,2), lvec(:,3))    ! alpha
   ang(2) = angle(lvec(:,3), lvec(:,1))    ! beta
   ang(3) = angle(lvec(:,1), lvec(:,2))    ! gamma
   
   lvec(1,1) = lenvec(1)
   lvec(2,1) = 0d0
   lvec(3,1) = 0d0   
   
   lvec(1,2) = lenvec(2)*cos(ang(3))
   lvec(2,2) = sqrt(lenvec(2)**2 - lvec(1,2)**2)
   lvec(3,2) = 0d0   
   
   lvec(1,3) = lenvec(3)*cos(ang(2))
   lvec(2,3) = (cos(ang(1))*lenvec(2)*lenvec(3) - lvec(1,2)*lvec(1,3))/lvec(2,2)
   lvec(3,3) = sqrt(lenvec(3)**2 - lvec(1,3)**2 - lvec(2,3)**2)
   
   do ii = 1, 3
   do ij = 1, 3
      if (abs(lvec(ij,ii)) .lt. 1.0d-14) lvec(ij,ii) = 0d0
   enddo
   enddo
   
 end subroutine latt_refin
   
end module
