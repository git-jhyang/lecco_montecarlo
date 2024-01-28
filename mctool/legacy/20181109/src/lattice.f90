module lattice

use constants,  only :  CHKOK, CHKERR,  &
                        VARLEN

use mt_serial,  only :  grnd

use iotool,     only :  io_printhead, &
                        io_printerr,  &
                        io_isodd

use parallel,   only :  pMast,          &
                        nrnk,           &
                        pp_bcast

use input,      only :  Var_Input

implicit none
save

!-------------------------------------------------------------------!
!       Var_Latt        : variables of lattice
!-------------------------------------------------------------------!
!       occ      : number of occupied site. can be differed from
!                  'nSite' of 'vinp' if there is defect
!       conn     : if 'lconn(0)' of 'vinp' is .true., conn = 2, otherwise 1
!       ldim     : dimension of lattice, 'nx', 'ny', and 'nz' of 'vinp'
!       num2lidx : convert n-th occputation number to lattice index
!-------------------------------------------------------------------!

type, public :: Var_Latt
   integer                              :: nSite
   integer, dimension(:,:), allocatable :: num2lidx
   logical                              :: iserr = .false.
   logical                              :: pp    = .false.
end type Var_Latt

!-------------------------------------------------------------------!
!       LattStructure : lattice structure, saving occpuation number
!       typ     : type lattice structure
!-------------------------------------------------------------------!

integer,                                private :: nNear
integer,                                private :: connectivity
integer, dimension(3),                  private :: LattDim
integer, dimension(:,:,:), allocatable, private :: LattStructure
logical,                                private :: lat_isinit = .false.
logical,                                private :: lat_iserr  = .false.
integer,                                private :: LattType

!-------------------------------------------------------------------!
!       functions and subroutines
!-------------------------------------------------------------------!
!       = functions = 
!       lat_load(vinp) : occupying lattice and return 'vlat'
!       lat_growth(vlat) : growing lattice based on diffusion-limited
!                       aggregation (DLA) model. if there is a defect,
!                       return .true.
!       lat_findnear(vlat,lidx,lidx_nears) : find nesrest occupied sites
!                       of 'lidx' and save in 'lidx_nears'. if there is 
!                       more than one occupied site, resturn .true.
!       lidx2num(lidx) : get index 'lidx' and return occupation number 'i'
!      
!       = subroutines=
!       pp_bcast_vlat(vlat) : broadcast 'vlat' to all processors
!       lat_origin(vlat,lidx) : get origin 'lidx' of lattice 
!       lat_pntgen(vlat,lidx) : get random initial position of particle 
!                       at outside of clustor
!       lat_pntmove(vlat,lidx) : move point particle to random direction
!       lat_print(vlat) : print lattice
!-------------------------------------------------------------------!

public ::       lat_init,       lat_fin

private ::      lat_bcast,      lat_growth,     &
                lat_gen_ori,    lat_gen_pnt,    &
                lat_move_pnt,   lat_find_near,  &
                lidx2num,       lat_writelatt,  &
                lat_shift,      lat_print

contains

!-------------------------------------------------------------------!
!         load lattice
!-------------------------------------------------------------------!

 function lat_init(vinp) result(vlat)
   implicit none
   type(Var_Input), intent(in) :: vinp
   type(Var_Latt)              :: vlat

   logical        :: isdefect

   if (lat_isinit) then
      if (pMast) call io_printerr('lattice initialization','lattice is already initialized')
      vlat%iserr = .true.
      return
   endif

   if (pMast) then
      vlat%nSite   = vinp%nSite
      connectivity = 1
      if (vinp%isConn(0)) connectivity = 2
      LattType     = vinp%LattType
      LattDim      = [vinp%nx, vinp%ny, vinp%nz]

      select case (vinp%mode)
      case ('OPT')
         select case (LattType)
         case (0)
            nNear = 6
         case (1)
            nNear = 8
         case (2,3)
            nNear = 12
         end select
         isdefect = lat_growth(vlat)
         if (isdefect) then
            write (*,200) ' Vacnacy detected !!!-----------------------------------------------------'
            write (*,201) '   Total number of lattice sites : ', vinp%nSiteTotal
            write (*,201) '   Number of occupied sites      : ', vlat%nSite
            write (*,202) '   Lattice occupation ratio      : ', float(vlat%nSite)/float(vinp%nSiteTotal)
            write (*,*)
         endif
         write (*,200) ' DLA growth finished -----------------------------------------------------'
         write (*,*)
      case ('GEN')
         allocate (vlat%num2lidx(3,vlat%nSite))
         vlat%num2lidx(:,1) = [2, 2, 2]
         if (vinp%isConn(0)) vlat%num2lidx(:,2) = [2,3,2]
      end select
   endif

   call lat_bcast(vlat)

200 format(1x,A)
201 format(1x,A,I8)
202 format(1x,A,F13.4)

   lat_isinit = .true.

 end function lat_init

!-------------------------------------------------------------------!

 subroutine lat_fin(vlat)
   implicit none
   type(Var_Latt), intent(inout) :: vlat

   if (.not. lat_isinit) then
      if (pMast) call io_printerr('lattice finalize','lattice is not initialized')
      vlat%iserr = .true. 
      return
   endif

   if (allocated (vlat%num2lidx)) deallocate (vlat%num2lidx)
   if (allocated (LattStructure)) deallocate (LattStructure)

   vlat%pp    = .false.
   lat_isinit = .false.

 end subroutine lat_fin

!-------------------------------------------------------------------!
!          broadcast lattice variables
!-------------------------------------------------------------------!

 subroutine lat_bcast(vlat)
   implicit none
   type(Var_Latt), intent(inout) :: vlat

   call pp_bcast(vlat%iserr)
   if (vlat%iserr) return
   if (vlat%pp) then
      if (pMast) call io_printerr('lattice bcast','lattice is already bcasted')
      vlat%iserr = .true.
      return
   endif

   call pp_bcast(vlat%nSite)
   if (.not. allocated (vlat%num2lidx)) allocate (vlat%num2lidx(3, vlat%nSite))
   call pp_bcast(vlat%num2lidx, 3, vlat%nSite)

   vlat%pp = .true.

 end subroutine lat_bcast

!-------------------------------------------------------------------!
!         growing lattice with DLA model
!-------------------------------------------------------------------!

 function lat_growth(vlat) result(isdefect)
   implicit none
   type(Var_Latt), intent(inout) :: vlat
   logical                       :: isdefect

   logical :: iscontact

   integer :: iConn, iSite, ixyz, iCnt

   integer, dimension(3,2)              :: lidx_near
   integer, dimension(:,:), allocatable :: num2lidx_tmp

   select case (LattType)
   case (0) ! S - simple
      allocate (LattStructure(LattDim(1), LattDim(2), LattDim(3)))
   case (1) ! BOC - body cendtered
   case (2) ! FC - face centered
      allocate (LattStructure(LattDim(1)*2, LattDim(2)*2, LattDim(3)*2))
   case (3) ! HCP - hexagonal closed packing
      allocate (LattStructure(LattDim(1), LattDim(2), LattDim(3)*2))
   end select

   allocate (num2lidx_tmp(3,vlat%nSite))

   LattStructure(:,:,:) = 0
   num2lidx_tmp(:,:)    = 0
   isdefect             = .false.

   lidx_near(:,1) = lat_gen_ori()
   if (connectivity .eq. 2) lidx_near(:,2) = lat_move_pnt(lidx_near(:,1))
   do iConn = 1, connectivity
      call lat_writelatt(lidx_near(:,iConn), iConn)
      num2lidx_tmp(:,iConn) = lidx_near(:,iConn)
  enddo

   iSite = connectivity + 1
   iCnt  = 0

   L_GROWTH: do while (iSite .le. vlat%nSite)
      lidx_near(:,1) = lat_gen_pnt()
      iscontact      = lat_find_near(lidx_near(:,1))

      L_DIFFUS: do while (.not. iscontact)
         lidx_near(:,1) = lat_move_pnt(lidx_near(:,1))
         iscontact    = lat_find_near(lidx_near(:,1))
      enddo L_DIFFUS
      call lat_writelatt(lidx_near(:,1), iSite)
      num2lidx_tmp(:,iSite) = lidx_near(:,1)

      if (connectivity .eq. 2) then
         lidx_near(:,2) = lat_move_pnt(lidx_near(:,1))
         isdefect = .true.
         do ixyz = 1, 3
            if (lidx_near(ixyz,1) .ne. lidx_near(ixyz,2)) isdefect = .false.
         enddo
      endif

      if (.not. isdefect) then
         iSite = iSite + 1
         if (connectivity .eq. 2) then
            call lat_writelatt(lidx_near(:,2), iSite)
            num2lidx_tmp(:,iSite) = lidx_near(:,2)
            iSite = iSite + 1
         endif
      else
         iCnt = iCnt + 1
         call lat_writelatt(lidx_near(:,1), 0)
      endif

      if (iCnt .gt. 50*vlat%nSite) then
         vlat%nSite = iSite - 1
         exit L_GROWTH
      endif
   enddo L_GROWTH

   allocate (vlat%num2lidx(3, vlat%nSite))
   vlat%num2lidx(:,:) = 0
   vlat%num2lidx(:,1:vlat%nSite) = num2lidx_tmp(:,1:vlat%nSite)
   deallocate (num2lidx_tmp, LattStructure)

 end function lat_growth

!-------------------------------------------------------------------!

 function lat_gen_ori() result(lidx)
   implicit none
   integer, dimension(3) :: lidx

   logical :: ix_odd, iy_odd, iz_odd

   select case(LattType)
   case (0)
      lidx(1) = floor(LattDim(1)*0.5d0) + 1
      lidx(2) = floor(LattDim(2)*0.5d0) + 1
      lidx(3) = floor(LattDim(3)*0.5d0) + 1
   case (2)
      lidx(1) = LattDim(1) + 1
      lidx(2) = LattDim(2) + 1
      lidx(3) = LattDim(3) + 1
      ix_odd = io_isodd(lidx(1))
      iy_odd = io_isodd(lidx(2))
      iz_odd = io_isodd(lidx(3))
      if     ((ix_odd .and. iy_odd) .or. ((.not. ix_odd) .and. (.not. iy_odd))) then 
         if (.not. iz_odd) lidx(3) = lidx(3) - 1
      elseif (((.not. ix_odd) .and. iy_odd) .or. (ix_odd .and. (.not. iy_odd))) then
         if (      iz_odd) lidx(3) = lidx(3) - 1
      endif 
   case (3)
      lidx(1) = floor(LattDim(1)*0.5d0) + 1
      lidx(2) = floor(LattDim(2)*0.5d0) + 1
      lidx(3) = LattDim(3) + 1
   end select

 end function lat_gen_ori

!-------------------------------------------------------------------!

 function lat_gen_pnt() result(lidx)
   implicit none
   integer, dimension(3) :: lidx
   
   logical :: ix_odd, iy_odd, iz_odd

   find_empty : do
      select case(LattType)
      case(0)
         lidx(1) = floor(grnd()*LattDim(1)) + 1
         lidx(2) = floor(grnd()*LattDim(2)) + 1
         lidx(3) = floor(grnd()*LattDim(3)) + 1
      case(2)
         lidx(1) = floor(grnd()*LattDim(1)*2) + 1
         lidx(2) = floor(grnd()*LattDim(2)*2) + 1
         lidx(3) = floor(grnd()*LattDim(3)*2) + 1
         ix_odd = io_isodd(lidx(1))
         iy_odd = io_isodd(lidx(2))
         iz_odd = io_isodd(lidx(3))
         if     ((ix_odd .and. iy_odd) .or. ((.not. ix_odd) .and. (.not. iy_odd))) then
            if (.not. iz_odd) lidx(3) = lidx(3) - 1
         elseif (((.not. ix_odd) .and. iy_odd) .or. (ix_odd .and. (.not. iy_odd))) then
            if (      iz_odd) lidx(3) = lidx(3) + 1
         endif
      case (3)
         lidx(1) = floor(grnd()*LattDim(1)) + 1
         lidx(2) = floor(grnd()*LattDim(2)) + 1
         lidx(3) = floor(grnd()*LattDim(3)*2) + 1
      end select
      if (lidx2num(lidx) .eq. 0) exit find_empty
   enddo find_empty

 end function lat_gen_pnt

!-------------------------------------------------------------------!

 function lat_move_pnt(lidx_in) result(lidx_out)
   implicit none
   integer, dimension(3), intent(in) :: lidx_in
   integer, dimension(3)             :: lidx_out

   integer :: ixyz, iNear, cnt

   logical, dimension(0:nNear-1) :: isTried

   isTried(:) = .false.

   cnt = 0
   LFE : do ! loop for find empty
      do
         iNear = floor(grnd()*nNear)
         if (.not. isTried(iNear)) exit
      enddo
      lidx_out = lidx_in + lat_shift(iNear, io_isodd(lidx_in(3)))
      do ixyz = 1, 3
         if (lidx_out(ixyz) .gt. size(LattStructure,ixyz)) lidx_out(ixyz) = 1
         if (lidx_out(ixyz) .eq. 0) lidx_out(ixyz) = size(LattStructure,ixyz)
      enddo

      if (lidx2num(lidx_out) .eq. 0) return ! find empty site

      isTried(iNear) = .true.
      do iNear = 0, nNear-1
         if (.not. isTried(iNear)) cycle LFE
      enddo
      lidx_out = lidx_in
      return ! fail to find empty site

   enddo LFE

 end function lat_move_pnt 

!-------------------------------------------------------------------!

 function lat_find_near(lidx_in, lidx_nears) result(isContact)
   implicit none
   integer, dimension(3),                intent(in)  :: lidx_in
   integer, dimension(3,0:12), optional, intent(out) :: lidx_nears
   logical                                           :: isContact

   integer :: ixyz, iNear, iCnt

   integer, dimension(3) :: lidx
   
   if (present(lidx_nears)) lidx_nears(:,:) = 0

   iscontact = .false.   
   iCnt      = 0

   L1NEAR : do iNear = 0, nNear-1
      lidx = lidx_in + lat_shift(iNear, io_isodd(lidx_in(3)))
      do ixyz = 1, 3
         if (lidx(ixyz) .gt. size(LattStructure,ixyz)) lidx(ixyz) = 1
         if (lidx(ixyz) .eq. 0) lidx(ixyz) = size(LattStructure,ixyz)
      enddo
      if (lidx2num(lidx) .ne. 0) then
         iCnt = iCnt + 1
         iscontact = .true.
         if (present(lidx_nears)) lidx_nears(:,iCnt) = lidx
         if (.not. present(lidx_nears)) return
      endif
   enddo L1NEAR

   if (present(lidx_nears)) lidx_nears(1,0) = iCnt

 end function lat_find_near

!-------------------------------------------------------------------!

 function lidx2num(lidx) result (iSite)
   implicit none
   integer, dimension(3), intent(in) :: lidx
   integer                           :: iSite

   iSite = LattStructure(lidx(1),lidx(2),lidx(3))

 end function lidx2num

!-------------------------------------------------------------------!

 subroutine lat_writelatt(lidx, iSite)
   implicit none
   integer, dimension(3), intent(in) :: lidx
   integer,               intent(in) :: iSite
 
   LattStructure(lidx(1),lidx(2),lidx(3)) = iSite

 end subroutine lat_writelatt

!-------------------------------------------------------------------!

 function lat_shift(iNear, isoddZ) result (shift)
   implicit none
   integer, intent(in)   :: iNear
   logical, intent(in)   :: isoddZ
   integer, dimension(3) :: shift

   select case (LattType)
   case (0)
      if (iNear .eq. 0) shift = [ 1,  0,  0]
      if (iNear .eq. 1) shift = [-1,  0,  0]
      if (iNear .eq. 2) shift = [ 0,  1,  0]
      if (iNear .eq. 3) shift = [ 0, -1,  0]
      if (iNear .eq. 4) shift = [ 0,  0,  1]
      if (iNear .eq. 5) shift = [ 0,  0, -1]
   case (1)
      if (iNear .eq. 0) shift = [ 1,  1,  1]
      if (iNear .eq. 1) shift = [ 1,  1, -1]
      if (iNear .eq. 2) shift = [ 1, -1,  1]
      if (iNear .eq. 3) shift = [ 1, -1, -1]
      if (iNear .eq. 4) shift = [-1,  1,  1]
      if (iNear .eq. 5) shift = [-1,  1, -1]
      if (iNear .eq. 6) shift = [-1, -1,  1]
      if (iNear .eq. 7) shift = [-1, -1, -1]
   case (2)
      if (iNear .eq.  0) shift = [ 1,  0, -1]
      if (iNear .eq.  1) shift = [-1,  0, -1]
      if (iNear .eq.  2) shift = [ 0,  1, -1]
      if (iNear .eq.  3) shift = [ 0, -1, -1]
      if (iNear .eq.  4) shift = [ 1,  1,  0]
      if (iNear .eq.  5) shift = [ 1, -1,  0]
      if (iNear .eq.  6) shift = [-1,  1,  0]
      if (iNear .eq.  7) shift = [-1, -1,  0]
      if (iNear .eq.  8) shift = [ 1,  0,  1]
      if (iNear .eq.  9) shift = [-1,  0,  1]
      if (iNear .eq. 10) shift = [ 0,  1,  1]
      if (iNear .eq. 11) shift = [ 0, -1,  1]
   case (3)
      if (iNear .eq.  0) shift = [ 0,  1,  0]
      if (iNear .eq.  1) shift = [ 0, -1,  0]
      if (iNear .eq.  2) shift = [ 1,  0,  0]
      if (iNear .eq.  3) shift = [ 1, -1,  0]
      if (iNear .eq.  4) shift = [-1,  0,  0]
      if (iNear .eq.  5) shift = [-1,  1,  0]
      if (iNear .eq.  6) shift = [ 0,  0,  1]
      if (iNear .eq.  7) shift = [ 0,  0, -1]
     if (isoddZ) then
      if (iNear .eq.  8) shift = [ 0, -1,  1]
      if (iNear .eq.  9) shift = [-1,  0,  1]
      if (iNear .eq. 10) shift = [ 0, -1, -1]
      if (iNear .eq. 11) shift = [-1,  0, -1]
     else
      if (iNear .eq.  8) shift = [ 0,  1,  1]
      if (iNear .eq.  9) shift = [ 1,  0,  1]
      if (iNear .eq. 10) shift = [ 0,  1, -1]
      if (iNear .eq. 11) shift = [ 1,  0, -1]
     endif
   end select
 end function lat_shift 

!-------------------------------------------------------------------!

 subroutine lat_print()
   implicit none

   integer :: ix, iy, iz
   logical :: ix_odd, iy_odd, iz_odd, ix_eve, iy_eve, iz_eve

   select case (LattType)
   case (2)
      do iz = 1, LattDim(3)*2
         iz_odd = io_isodd(iz)
         iz_eve = .not. iz_odd
         do ix = 1, LattDim(1)*2
            ix_odd = io_isodd(ix)
            ix_eve = .not. ix_odd
            do iy = 1, LattDim(2)*2
               iy_odd = io_isodd(iy)
               iy_eve = .not. iy_odd
               if (iz_odd) then
                  if ((ix_odd .and. iy_odd) .or. (ix_eve .and. iy_eve)) then
                     if (LattStructure(ix,iy,iz) .eq. 0) write (*,100,advance='no') 'O'
                     if (LattStructure(ix,iy,iz) .ne. 0) write (*,100,advance='no') 'X'
                  else
                     write (*,100,advance='no') ' '
                  endif
               else
                  if ((ix_odd .and. iy_odd) .or. (ix_eve .and. iy_eve)) then
                     write (*,100,advance='no') ' '
                  else
                     if (LattStructure(ix,iy,iz) .eq. 0) write (*,100,advance='no') 'O'
                     if (LattStructure(ix,iy,iz) .ne. 0) write (*,100,advance='no') 'X'
                  endif
               endif
            enddo
            write (*,*)
         enddo
         write (*,*)
         write (*,*) '-------------------------------------'
         write (*,*)
      enddo  
   end select

100 format(A2)

 end subroutine lat_print

!-------------------------------------------------------------------!

end module lattice
