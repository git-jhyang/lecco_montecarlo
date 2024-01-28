module lattice

use constants,  only :  CHKOK, CHKERR,  &
                        VARLEN

use mt_serial,  only :  grnd

use mcio,       only :  mcio_printhead, &
                        mcio_printerr,  &
                        mcio_isodd

use parallel,   only :  pMast,          &
                        nrnk,           &
                        pp_bcast

use input,      only :  Var_Input

implicit none
save

!-------------------------------------------------------------------!
!       Var_Latt        : variables of lattice
!-------------------------------------------------------------------!
!       occ     : number of occupied site. can be differed from
!                 'nSite' of 'vinp' if there is defect
!       conn    : if 'lconn(0)' of 'vinp' is .true., conn = 2, otherwise 1
!       ldim    : dimension of lattice, 'nx', 'ny', and 'nz' of 'vinp'
!       num2id  : convert n-th occputation number to lattice index
!-------------------------------------------------------------------!

type, public :: Var_Latt
   integer                              :: nSite
   integer, dimension(:,:), allocatable :: num2id
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
character(LEN=VARLEN),                  private :: LattType

!-------------------------------------------------------------------!
!       functions and subroutines
!-------------------------------------------------------------------!
!       = functions = 
!       lat_load(vinp) : occupying lattice and return 'vlat'
!       lat_growth(vlat) : growing lattice based on diffusion-limited
!                       aggregation (DLA) model. if there is a defect,
!                       return .true.
!       lat_findnear(vlat,id,id_nears) : find nesrest occupied sites
!                       of 'id' and save in 'id_nears'. if there is 
!                       more than one occupied site, resturn .true.
!       id2num(id) : get index 'id' and return occupation number 'i'
!      
!       = subroutines=
!       pp_bcast_vlat(vlat) : broadcast 'vlat' to all processors
!       lat_origin(vlat,id) : get origin 'id' of lattice 
!       lat_pntgen(vlat,id) : get random initial position of particle 
!                       at outside of clustor
!       lat_pntmove(vlat,id) : move point particle to random direction
!       lat_print(vlat) : print lattice
!-------------------------------------------------------------------!

public ::       lat_init,         lat_fin

private ::      lat_bcast,      lat_growth,     &
                lat_gen_ori,    lat_gen_pnt,    &
                lat_move_pnt,   lat_find_near,  &
                id2num,         lat_writelatt,  &
                lat_shift,      lat_print

contains

!-------------------------------------------------------------------!
!         load lattice
!-------------------------------------------------------------------!

 function lat_init(vinp) result(vlat)
   implicit none
   type(Var_Input), intent(in)    :: vinp

   type(Var_Latt) :: vlat
   logical        :: isdefect

   if (lat_isinit) then
      if (pMast) call mcio_printerr('lattice initialization','lattice is already initialized')
      vlat%iserr = .true.
      return
   endif

   if (pMast) then
      vlat%nSite   = vinp%nSite
      connectivity = 1
      if (vinp%isConnect(0)) connectivity = 2
      LattType     = vinp%LattType
      LattDim      = [vinp%nx, vinp%ny, vinp%nz]

      select case (vinp%mode)
      case ('OPT')
         select case (LattType)
         case ('S')
            nNear = 6
         case ('FC','HCP')
            nNear = 12
         case ('BC')
            nNear = 8
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
         allocate (vlat%num2id(3,vlat%nSite))
         vlat%num2id(:,1) = [2, 2, 2]
         if (vinp%isConnect(0)) vlat%num2id(:,2) = [2,3,2]
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
      if (pMast) call mcio_printerr('lattice finalize','lattice is not initialized')
      vlat%iserr = .true. 
      return
   endif

   if (allocated (vlat%num2id)) deallocate (vlat%num2id)
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
      if (pMast) call mcio_printerr('lattice bcast','lattice is already bcasted')
      vlat%iserr = .true.
      return
   endif

   call pp_bcast(vlat%nSite)
   if (.not. allocated (vlat%num2id)) allocate (vlat%num2id(3, vlat%nSite))
   call pp_bcast(vlat%num2id, 3, vlat%nSite)

   vlat%pp = .true.

 end subroutine lat_bcast

!-------------------------------------------------------------------!
!         growing lattice with DLA model
!-------------------------------------------------------------------!

 function lat_growth(vlat) result(isdefect)
   implicit none
   type(Var_Latt), intent(inout)        :: vlat

   logical                              :: iscontact, isdefect
   integer                              :: iConn, iSite, ixyz, iCnt
   integer, dimension(3,2)              :: id_near
   integer, dimension(:,:), allocatable :: num2id_tmp

   select case (LattType)
   case ('S')
      allocate (LattStructure(LattDim(1), LattDim(2), LattDim(3)))
   case ('FC')
      allocate (LattStructure(LattDim(1)*2, LattDim(2)*2, LattDim(3)*2))
   end select

   allocate (num2id_tmp(3,vlat%nSite))

   LattStructure(:,:,:) = 0
   num2id_tmp(:,:)      = 0
   isdefect             = .false.

   id_near(:,1) = lat_gen_ori()
   if (connectivity .eq. 2) id_near(:,2) = lat_move_pnt(id_near(:,1))
   do iConn = 1, connectivity
      call lat_writelatt(id_near(:,iConn), iConn)
      num2id_tmp(:,iConn) = id_near(:,iConn)
  enddo

   iSite = connectivity + 1
   iCnt  = 0

   loop_growth: do while (iSite .le. vlat%nSite)
      id_near(:,1) = lat_gen_pnt()
      iscontact    = lat_find_near(id_near(:,1))

      loop_diffuse: do while (.not. iscontact)
         id_near(:,1) = lat_move_pnt(id_near(:,1))
         iscontact    = lat_find_near(id_near(:,1))
      enddo loop_diffuse
      call lat_writelatt(id_near(:,1), iSite)
      num2id_tmp(:,iSite) = id_near(:,1)

      if (connectivity .eq. 2) then
         id_near(:,2) = lat_move_pnt(id_near(:,1))
         isdefect = .true.
         do ixyz = 1, 3
            if (id_near(ixyz,1) .ne. id_near(ixyz,2)) isdefect = .false.
         enddo
      endif

      if (.not. isdefect) then
         iSite = iSite + 1
         if (connectivity .eq. 2) then
            call lat_writelatt(id_near(:,2), iSite)
            num2id_tmp(:,iSite) = id_near(:,2)
            iSite = iSite + 1
         endif
      else
         iCnt = iCnt + 1
         call lat_writelatt(id_near(:,1), 0)
      endif

      if (iCnt .gt. 50*vlat%nSite) then
         vlat%nSite = iSite - 1
         exit loop_growth
      endif
   enddo loop_growth

   allocate (vlat%num2id(3, vlat%nSite))
   vlat%num2id(:,:) = 0
   vlat%num2id(:,1:vlat%nSite) = num2id_tmp(:,1:vlat%nSite)
   deallocate (num2id_tmp, LattStructure)

 end function lat_growth

!-------------------------------------------------------------------!

 function lat_gen_ori() result(id)
   implicit none
   integer, dimension(3) :: id

   logical :: ix_odd, iy_odd, iz_odd

   select case(LattType)
   case('S')
      id(1) = LattDim(1)/2 + 1
      id(2) = LattDim(2)/2 + 1
      id(3) = LattDim(3)/2 + 1
   case('FC')
      id(1) = LattDim(1) + 1
      id(2) = LattDim(2) + 1
      id(3) = LattDim(3) + 1
      ix_odd = mcio_isodd(id(1))
      iy_odd = mcio_isodd(id(2))
      iz_odd = mcio_isodd(id(3))
      if     ((ix_odd .and. iy_odd) .or. ((.not. ix_odd) .and. (.not. iy_odd))) then 
         if (.not. iz_odd) id(3) = id(3) - 1
      elseif (((.not. ix_odd) .and. iy_odd) .or. (ix_odd .and. (.not. iy_odd))) then
         if (      iz_odd) id(3) = id(3) - 1
      endif 
   end select

 end function lat_gen_ori

!-------------------------------------------------------------------!

 function lat_gen_pnt() result(id)
   implicit none
   integer, dimension(3) :: id
   
   logical :: ix_odd, iy_odd, iz_odd

   find_empty : do
      select case(LattType)
      case('S')
         id(1) = floor(grnd()*LattDim(1)) + 1
         id(2) = floor(grnd()*LattDim(2)) + 1
         id(3) = floor(grnd()*LattDim(3)) + 1
      case('FC')
         id(1) = floor(grnd()*LattDim(1)*2) + 1
         id(2) = floor(grnd()*LattDim(2)*2) + 1
         id(3) = floor(grnd()*LattDim(3)*2) + 1
         ix_odd = mcio_isodd(id(1))
         iy_odd = mcio_isodd(id(2))
         iz_odd = mcio_isodd(id(3))
         if     ((ix_odd .and. iy_odd) .or. ((.not. ix_odd) .and. (.not. iy_odd))) then
            if (.not. iz_odd) id(3) = id(3) - 1
         elseif (((.not. ix_odd) .and. iy_odd) .or. (ix_odd .and. (.not. iy_odd))) then
            if (      iz_odd) id(3) = id(3) + 1
         endif
      end select
      if (id2num(id) .eq. 0) exit find_empty
   enddo find_empty

 end function lat_gen_pnt

!-------------------------------------------------------------------!

 function lat_move_pnt(id) result(id_new)
   implicit none
   integer, dimension(3), intent(in) :: id

   integer :: ixyz, iNear, cnt
   integer, dimension(3)     :: id_new
   logical, dimension(0:nNear-1) :: isEmpty

   isEmpty(:) = .true.

   cnt = 0
   LFE : do ! loop for find empty
      do
         iNear = floor(grnd()*nNear)
         if (isEmpty(iNear)) exit
      enddo
      id_new = id + lat_shift(iNear)
      do ixyz = 1, 3
         if (id_new(ixyz) .gt. size(LattStructure,ixyz)) id_new(ixyz) = 1
         if (id_new(ixyz) .eq. 0) id_new(ixyz) = size(LattStructure,ixyz)
      enddo

      if (id2num(id_new) .eq. 0) return ! find empty site

      isEmpty(iNear) = .false.
      do iNear = 0, nNear-1
         if (isEmpty(iNear)) cycle LFE
      enddo
      id_new = id
      return ! fail to find empty site

   enddo LFE

 end function lat_move_pnt 

!-------------------------------------------------------------------!

 function lat_find_near(id, id_nears) result(isContact)
   implicit none
   integer, dimension(3),                intent(in)  :: id
   integer, dimension(3,0:12), optional, intent(out) :: id_nears

   logical               :: isContact
   integer               :: ixyz, iNear, iCnt
   integer, dimension(3) :: id_new
   
   if (present(id_nears)) id_nears(:,:) = 0
   iscontact                            = .false.   

   iCnt = 0

   loop_1st_near : do iNear = 0, nNear-1
      id_new = id + lat_shift(iNear)
      do ixyz = 1, 3
         if (id_new(ixyz) .gt. size(LattStructure,ixyz)) id_new(ixyz) = 1
         if (id_new(ixyz) .eq. 0) id_new(ixyz) = size(LattStructure,ixyz)
      enddo
      if (id2num(id_new) .ne. 0) then
         iCnt = iCnt + 1
         if (present(id_nears)) id_nears(:,iCnt) = id_new
         if (.not. present(id_nears)) exit loop_1st_near
      endif
   enddo loop_1st_near

   if (iCnt .ne. 0) then
      iscontact = .true.
      if (present(id_nears)) id_nears(1,0) = iCnt
   endif

 end function lat_find_near

!-------------------------------------------------------------------!

 function id2num(id) result(num)
   implicit none
   integer, dimension(3), intent(in) :: id

   integer :: num

   num = LattStructure(id(1),id(2),id(3))

 end function id2num

!-------------------------------------------------------------------!

 subroutine lat_writelatt(id, num)
   implicit none
   integer, dimension(3), intent(in) :: id
   integer,               intent(in) :: num
 
   LattStructure(id(1), id(2), id(3)) = num

 end subroutine lat_writelatt

!-------------------------------------------------------------------!

 function lat_shift(iNear) result(shift)
   implicit none
   integer, intent(in) :: iNear

   integer, dimension(3) :: shift

   select case (LattType)
   case ('S')
      if (iNear .eq. 0) shift = [ 1,  0,  0]
      if (iNear .eq. 1) shift = [-1,  0,  0]
      if (iNear .eq. 2) shift = [ 0,  1,  0]
      if (iNear .eq. 3) shift = [ 0, -1,  0]
      if (iNear .eq. 4) shift = [ 0,  0,  1]
      if (iNear .eq. 5) shift = [ 0,  0, -1]
   case ('FC')
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
   case ('BC')
      if (iNear .eq. 0) shift = [ 1,  1,  1]
      if (iNear .eq. 1) shift = [ 1,  1, -1]
      if (iNear .eq. 2) shift = [ 1, -1,  1]
      if (iNear .eq. 3) shift = [ 1, -1, -1]
      if (iNear .eq. 4) shift = [-1,  1,  1]
      if (iNear .eq. 5) shift = [-1,  1, -1]
      if (iNear .eq. 6) shift = [-1, -1,  1]
      if (iNear .eq. 7) shift = [-1, -1, -1]
   end select
 end function lat_shift 

!-------------------------------------------------------------------!

 subroutine lat_print()
   implicit none

   integer :: ix, iy, iz
   logical :: ix_odd, iy_odd, iz_odd, ix_eve, iy_eve, iz_eve

   select case (LattType)
   case ('FC')
      do iz = 1, LattDim(3)*2
         iz_odd = mcio_isodd(iz)
         iz_eve = .not. iz_odd
         do ix = 1, LattDim(1)*2
            ix_odd = mcio_isodd(ix)
            ix_eve = .not. ix_odd
            do iy = 1, LattDim(2)*2
               iy_odd = mcio_isodd(iy)
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
