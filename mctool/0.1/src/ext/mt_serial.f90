! A C-program for MT19937: Real number version
!   genrand() generates one pseudorandom real number (double)
! which is uniformly distributed on [0,1]-interval, for each
! call. sgenrand(seed) set initial values to the working area
! of 624 words. Before genrand(), sgenrand(seed) must be
! called once. (seed is any 32-bit integer except for 0).
! Integer generator is obtained by modifying two lines.
!   Coded by Takuji Nishimura, considering the suggestions by
! Topher Cooper and Marc Rieffel in July-Aug. 1997.
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Library General Public
! License as published by the Free Software Foundation; either
! version 2 of the License, or (at your option) any later
! version.
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Library General Public License for more details.
! You should have received a copy of the GNU Library General
! Public License along with this library; if not, write to the
! Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
! 02111-1307  USA
!
! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
! When you use this, send an email to: matumoto@math.keio.ac.jp
! with an appropriate reference to your work.
!
!***********************************************************************
! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
!
!   genrand()      -> double precision function grnd()
!   sgenrand(seed) -> subroutine sgrnd(seed)
!                     integer seed
!
! This program uses the following non-standard intrinsics.
!   ishft(i,n): If n>0, shifts bits in i by n positions to left.
!               If n<0, shifts bits in i by n positions to right.
!   iand (i,j): Performs logical AND on corresponding bits of i and j.
!   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
!   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
!
!***********************************************************************
! this main() outputs first 1000 generated numbers
module mt_serial

! example
! this main() outputs first 1000 generated numbers
!      program main
!
!      implicit integer(i-n)
!      implicit double precision(a-h,o-z)
!
!      parameter(no=1000)
!      dimension r(0:7)
!
!      call sgrnd(4357)
!                         any nonzero integer can be used as a seed
!      do j=0,no-1
!         r(mod(j,8))=grnd()
!         if(mod(j,8).eq.7) then
!            write(*,'(8(f8.6,'' ''))') (r(k),k=0,7)
!         else if(j.eq.no-1) then
!            write(*,'(8(f8.6,'' ''))') (r(k),k=0,mod(no-1,8))
!         endif
!      end do
!      end program

contains 

!------------------------------------------------------------------!
! grnd() : double precision with range of [0,1)
! grni() : integer          with range of [0,2^32-1]
!------------------------------------------------------------------!

 subroutine sgrnd(seed)

   implicit integer(a-z)

! Period parameters
   integer, parameter :: N = 624
   integer :: mt(0:N-1)    ! the array for the state vector

   common /block/mti,mt
   save /block/
 
!      setting initial seeds to mt[N] using
!      the generator Line 25 of Table 1 in
!      [KNUTH 1981, The Art of Computer Programming
!         Vol. 2 (2nd Ed.), pp102]

   mt(0) = iand(seed,-1)
   do mti = 1, N-1
      mt(mti) = iand(69069 * mt(mti-1),-1)
   enddo
      
 end subroutine

!------------------------------------------------------------------!

 double precision function grnd()

   implicit none

   integer :: y

   y = grni()

   if (y .lt. 0) then
      grnd = (dble(y)+2.0d0**32)/(2.0d0**32)
   else
      grnd = dble(y)/(2.0d0**32)
   endif

 end function grnd

!------------------------------------------------------------------!

 integer function grni()

   implicit integer(a-z)

! Period parameters
   integer, parameter :: N     =  624
   integer, parameter :: N1    =  N+1
   integer, parameter :: M     =  397
   integer, parameter :: MATA  = -1727483681
!                                    constant vector a
   integer, parameter :: UMASK = -2147483648
!                                    most significant w-r bits
   integer, parameter :: LMASK =  2147483647
!                                    least significant r bits

! Tempering parameters
   integer, parameter :: TMASKB= -1658038656
   integer, parameter :: TMASKC= -272236544

   integer :: mt(0:N-1)
   !                     the array for the state vector
   common /block/mti,mt
   save   /block/
   data   mti/N1/
   !                     mti==N+1 means mt[N] is not initialized

   dimension mag01(0:1)
   data mag01/0, MATA/
   save mag01
   !                        mag01(x) = x * MATA for x=0,1

   TSHFTU(y)=ishft(y,-11)
   TSHFTS(y)=ishft(y,7)
   TSHFTT(y)=ishft(y,15)
   TSHFTL(y)=ishft(y,-18)
   
   if (mti .ge. N) then            ! generate N words at one time
      if (mti .eq. N+1) then       ! if sgrnd() has not been called,
         call sgrnd(4357)          ! a default initial seed is used
      endif
      do kk = 0, N-M-1
         y      = ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
         mt(kk) = ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
      enddo
   
      do kk = N-M, N-2
         y      = ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
         mt(kk) = ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
      enddo

      y       = ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
      mt(N-1) = ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
      mti     = 0
   endif
   
   y   = mt(mti)
   mti = mti+1
   y   = ieor(y,TSHFTU(y))
   y   = ieor(y,iand(TSHFTS(y),TMASKB))
   y   = ieor(y,iand(TSHFTT(y),TMASKC))
   y   = ieor(y,TSHFTL(y))
   
   grni = y
   
 end function grni

end module mt_serial
