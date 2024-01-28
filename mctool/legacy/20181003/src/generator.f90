module generator

use constants,  only :  VARLEN, dp,     &
                        p_xp,   p_yp,   &
                        p_zp,   pi

use mcio,       only :  mcio_isodd,     &
                        mcio_printerr

use input,      only :  Var_Input

use coordinate, only :  Var_Coord,      &
                        coo_trans,      &
                        coo_rot_v,      &
                        coo_rot_cv,     &
                        coo_align_v

 implicit none

 type, public :: Var_Gen
   integer                                    :: nAtom
   integer,       dimension(:),   allocatable :: nAtomSpec
   integer                                    :: iGrid
   integer                                    :: Cnt_all
   integer                                    :: Cnt_fail
   logical                                    :: iserr = .false.
   logical                                    :: isfin = .false.
   real(kind=dp), dimension(3,3)              :: LattVec
   real(kind=dp), dimension(:,:), allocatable :: coo
 end type

 integer,                                      private :: nSpec
 integer,                                      private :: nUnit
 integer,                                      private :: nAtoms_max
 integer,                                      private :: iTransY
 integer,                                      private :: iTransZ
 integer,                                      private :: iRotL
 integer,                                      private :: iRotG
 integer,                                      private :: nGridRotG
 integer,                                      private :: nGridRotL
 integer,       dimension(3),                  private :: nGridRotL_all
 integer,       dimension(2),                  private :: nGridTrans
 integer,       dimension(:),     allocatable, private :: nAtoms
 integer,       dimension(:),     allocatable, private :: iTypeB
 integer,       dimension(:,:),   allocatable, private :: nAtomSpec_Unit
 integer,       dimension(:,:),   allocatable, private :: OrderSpec_Unit

 logical,                                      private :: gen_isinit = .false.
 logical,                                      private :: gen_iserr  = .false.
 logical,       dimension(:,:),   allocatable, private :: isMove

 real(kind=dp),                                private :: r_cut
 real(kind=dp),                                private :: AngRotG_max
 real(kind=dp), dimension(3),                  private :: AngRotL_max
 real(kind=dp), dimension(3),                  private :: TransLatt
 real(kind=dp), dimension(2,2),                private :: DisTrans_range
 real(kind=dp), dimension(:,:,:), allocatable, private :: coo_orig
 real(kind=dp), dimension(:,:,:), allocatable, private :: coo_B
 real(kind=dp), dimension(:,:,:), allocatable, private :: coo_shift

 contains

!-------------------------------------------------------------------!

 function gen_init(vinp, vcoo) result(vgen)
   implicit none
   type(Var_Input), intent(in) :: vinp
   type(Var_Coord), intent(in) :: vcoo
   type(Var_Gen)               :: vgen

   integer :: iUnit, ixyz, iSpec, iSpec_Unit, id01, id02

   if (gen_isinit) then
      call mcio_printerr('generator initialize','generator is already initialized')
      vgen%iserr = .true.
      return
   endif

   nSpec      = vcoo%nSpec
   nUnit      = vinp%nUnit
   allocate (nAtoms(nUnit))
   nAtoms     = vinp%nAtoms
   nAtoms_max = maxval(nAtoms)
   vgen%nAtom = sum(nAtoms(:))

   allocate (nAtomSpec_Unit(size(vcoo%nAtomSpec_Unit,1),nUnit), &
             OrderSpec_Unit(0:size(vcoo%nAtomSpec_Unit,1),nUnit))

   nAtomSpec_Unit = vcoo%nAtomSpec_Unit
   OrderSpec_Unit = vcoo%nOrderSPec_Unit

   allocate (vgen%nAtomSpec(nSpec)) 
   vgen%nAtomSpec(:) = 0
   do iSpec = 1, nSpec
      do iUnit = 1, nUnit
      do iSpec_Unit = 1, OrderSpec_Unit(0,iUnit)
         if (OrderSpec_Unit(iSpec_Unit,iUnit) .ne. iSpec) cycle
         vgen%nAtomSpec(iSpec) = vgen%nAtomSpec(iSpec) + &
                                 nAtomSpec_Unit(iSpec_Unit,iUnit)
      enddo
      enddo
   enddo
  
   allocate (iTypeB(nAtoms(2)))
   do iSpec = 1, OrderSpec_Unit(0,2)
      id02 = sum(nAtomSpec_Unit(1:iSpec,2))
      id01 = id02 - nAtomSpec_Unit(iSpec,2) + 1
      iTypeB(id01:id02) = iSpec
   enddo

   allocate (coo_orig(3,nAtoms_max,nUnit), &
             coo_shift(3,nAtoms_max,nUnit))

   coo_orig       = vcoo%coo_unit

   nGridRotG      = vinp%nGridRotG
   nGridRotL_all  = vinp%nGridRotL
   nGridTrans     = vinp%nGridTrans
   AngRotG_max    = vinp%AngRotG
   AngRotL_max    = vinp%AngRotL
   DisTrans_range = vinp%DisTrans

   r_cut          = vinp%r_cut

   call gen_coo_rotB(coo_orig(:,1:nAtoms(2),2))
   call gen_calc_LattVec(vgen%LattVec)

   vgen%iGrid    = 0
   vgen%Cnt_all  = 0
   vgen%Cnt_fail = 0

   iTransY    = 1
   iTransZ    = 1
   iRotG      = 1
   iRotL      = 1

   allocate (vgen%coo(3,vgen%nAtom))
   gen_isinit = .true.

 end function gen_init

!-------------------------------------------------------------------!

 subroutine gen_load_coo(vgen)
   implicit none
   type(Var_Gen), intent(inout) :: vgen

   real(kind=dp) :: AngRotG, DisTransY, DisTransZ, r_min
   real(kind=dp), dimension(3) :: TransVec

   if (.not. gen_isinit) then
      call mcio_printerr('generator coord. loader','generator is not initialized')
      vgen%iserr = .true.
      return
   endif

   vgen%iGrid  = vgen%iGrid + 1

   loop_find : do
      if (vgen%isfin) exit loop_find

      vgen%Cnt_all = vgen%Cnt_all + 1

      AngRotG     = dble(iRotG-1)*AngRotG_max/dble(nGridRotG-1)
      TransVec(1) = 0d0
      TransVec(2) = DisTrans_range(1,1) + dble(iTransY-1)* &
                   (DisTrans_range(2,1) - DisTrans_range(1,1))/dble(nGridTrans(1)-1)
      TransVec(3) = DisTrans_range(1,2) + dble(iTransZ-1)* &
                   (DisTrans_range(2,2) - DisTrans_range(1,2))/dble(nGridTrans(2)-1)

      coo_shift(:,1:nAtoms(1),1) = coo_orig(:,1:nAtoms(1),1)
      coo_shift(:,1:nAtoms(2),2) = coo_B(:,:,iRotL)

      call coo_rot_cv(nAtoms(1), p_zp, coo_shift(:,1:nAtoms(1),1), AngRotG)
      call coo_trans(nAtoms(2), TransVec,  coo_shift(:,1:nAtoms(2),2))

      call coo_trans(nAtoms(1), TransLatt, coo_shift(:,1:nAtoms(1),1))
      call coo_trans(nAtoms(2), TransLatt, coo_shift(:,1:nAtoms(2),2))

      iRotL = iRotL + 1
      if (iRotL .gt. nGridRotL) then
         iTransY = iTransY + 1
         iRotL   = 1
      endif
      if (iTransY .gt. nGridTrans(1)) then
         iTransZ = iTransZ + 1
         iTransY = 1
      endif
      if (iTransZ .gt. nGridTrans(2)) then
         iRotG   = iRotG + 1
         iTransZ = 1
      endif
      if (iRotG .gt. nGridRotG) vgen%isfin = .true.

      r_min = gen_calc_min_distance()
      if (r_min .le. r_cut) then
         vgen%Cnt_fail = vgen%Cnt_fail + 1
         cycle loop_find
      endif
      call gen_coo_gather(vgen)
      exit loop_find
   enddo loop_find

 end subroutine gen_load_coo

!-------------------------------------------------------------------!

 subroutine gen_coo_rotB(coo_in)
   implicit none
   real(kind=dp), dimension(3,nAtoms(2)), intent(in) :: coo_in

   integer       :: iRotX, iRotY, iRotZ, iAtom1, iAtom2
   integer       :: iGridRot1, iGridRot2, nGridRotL_re
   logical       :: isDuplicate
   real(kind=dp) :: AngRotX, AngRotY, AngRotZ

   real(kind=dp),     dimension(3,nAtoms(2),3)      :: coo_tmp
   real(kind=dp),     dimension(:,:,:), allocatable :: coo_B_tmp
   character(LEN=50), dimension(2)                  :: info

   nGridRotL = nGridRotL_all(1)*nGridRotL_all(2)*nGridRotL_all(3)

   if (nGridRotL .le. 1) then
      if (nGridRotL .eq. 0) nGridRotL = 1
      allocate (coo_B(3,nAtoms(2),1))
      coo_B(1:3,1:nAtoms(2),1) = coo_in
      return
   endif

   allocate (coo_B(3,nAtoms(2),nGridRotL), &
             coo_B_tmp(3,nAtoms(2),nGridRotL))

   iGridRot1 = 0

   do iRotZ = 1, nGridRotL_all(3)
      AngRotZ = dble(iRotZ-1)*AngRotL_max(3)/dble(nGridRotL_all(3)-1)
      coo_tmp(:,:,3) = coo_in
      call coo_rot_cv(nAtoms(2), p_zp, coo_tmp(:,:,3), AngRotZ)
      do iRotY = 1, nGridRotL_all(2)
         AngRotY = dble(iRotY-1)*AngRotL_max(2)/dble(nGridRotL_all(2)-1)
         coo_tmp(:,:,2) = coo_tmp(:,:,3)
         call coo_rot_cv(nAtoms(2), p_yp, coo_tmp(:,:,2), AngRotY)
         do iRotX = 1, nGridRotL_all(1)
            AngRotX = dble(iRotX-1)*AngRotL_max(1)/dble(nGridRotL_all(1)-1)
            iGridRot1 = iGridRot1 + 1
            coo_tmp(:,:,1) = coo_tmp(:,:,2)
            call coo_rot_cv(nAtoms(2), p_xp, coo_tmp(:,:,1), AngRotX)
            coo_B(:,:,iGridRot1) = coo_tmp(:,:,1)
         enddo
      enddo
   enddo

   coo_B_tmp(:,:,1) = coo_B(:,:,1)
   nGridRotL_re      = 1

   loop_grid_cnt : do iGridRot1 = 2, nGridRotL
      loop_grid_ref : do iGridRot2 = 1, nGridRotL_re
         loop_atom_cnt : do iAtom1 = 1, nAtoms(2)
            write (info(1),'(I2,3F15.8)') iTypeB(iAtom1), coo_B(:,iAtom1,iGridRot1)
            loop_atom_ref : do iAtom2 = 1, nAtoms(2)
               write (info(2),'(I2,3F15.8)') iTypeB(iAtom2), coo_B_tmp(:,iAtom2,iGridRot2)
               if (info(1) .eq. info(2)) cycle loop_atom_cnt
            enddo loop_atom_ref
            cycle loop_grid_ref
         enddo loop_atom_cnt
         cycle loop_grid_cnt
      enddo loop_grid_ref
      nGridRotL_re = nGridRotL_re + 1
      coo_B_tmp(:,:,nGridRotL_re) = coo_B(:,:,iGridRot1)
   enddo loop_grid_cnt

   nGridRotL = nGridRotL_re
   deallocate (coo_B)
   allocate (coo_B(3,nAtoms(2),nGridRotL))
   coo_B(1:3,1:nAtoms(2),1:nGridRotL) = coo_B_tmp(1:3,1:nAtoms(2),1:nGridRotL)
   deallocate (coo_B_tmp)

 end subroutine gen_coo_rotB

!-------------------------------------------------------------------!

 subroutine gen_calc_LattVec(LattVec)
   implicit none
   real(kind=dp), dimension(3,3), intent(out) :: LattVec
 
   integer :: iUnit, iAtom, ixyz

   real(kind=dp) :: r

   real(kind=dp), dimension(3)       :: vec
   real(kind=dp), dimension(3,nUnit) :: centroid, r_max


   LattVec(:,:)  = 0d0
   centroid(:,:) = 0d0
   r_max(:,:)    = 0d0

   do iUnit = 1, nUnit
      do iAtom = 1, nAtoms(iUnit)
         centroid(:,iUnit) = centroid(:,iUnit) + coo_orig(:,iAtom,iUnit)
      enddo
      centroid(:,iUnit) = centroid(:,iUnit)/dble(nAtoms(iUnit))
      do iAtom = 1, nAtoms(iUnit)
         vec = coo_orig(:,iAtom,iUnit) - centroid(:,iUnit)
         do ixyz = 1, 3
            r = abs(vec(ixyz))
            if (r .gt. r_max(ixyz,iUnit)) r_max(ixyz,iUnit) = r
         enddo
      enddo
   enddo

   vec = centroid(:,2) - centroid(:,1)
   r_max(:,2) = maxval(r_max(:,2))

   do ixyz = 1, 3
      r = abs(vec(ixyz)) + minval(r_max(ixyz,:))
      if (r .lt. maxval(r_max(ixyz,:))) r = maxval(r_max(ixyz,:))
      LattVec(ixyz,ixyz) = r + maxval(r_max(ixyz,:)) + 1.8d1
      TransLatt(ixyz)    = LattVec(ixyz,ixyz)*0.5d0
   enddo

 end subroutine gen_calc_LattVec

!-------------------------------------------------------------------!

 function gen_calc_min_distance() result(r_min)
   implicit none
   real(kind=dp) :: r_min

   integer       :: iAtom1, iAtom2
   real(kind=dp) :: r
   real(kind=dp), dimension(3) :: vec

   r_min = 1d4

   do iAtom1 = 1, nAtoms(1)
   do iAtom2 = 1, nAtoms(2)
      vec = coo_shift(:,iAtom1,1) - coo_shift(:,iAtom2,2)
      r   = dot_product(vec, vec)
      if (r .lt. r_min) r_min = r
   enddo
   enddo

   r_min = sqrt(r_min)

 end function gen_calc_min_distance

!-------------------------------------------------------------------!

 subroutine gen_coo_gather(vgen)
   implicit none
   type(Var_Gen), intent(inout) :: vgen

   integer :: iSpec, iUnit, iSpecUnit
   integer :: id01, id02, id03, id04, nAtomSpecStack

   loop_spec : do iSpec = 1, nSpec
      nAtomSpecStack = 0
      loop_unit : do iUnit = 1, nUnit
      loop_order : do iSpecUnit = 1, OrderSpec_Unit(0,iUnit)
         if (OrderSpec_Unit(iSpecUnit,iUnit) .ne. iSpec) cycle loop_order
         id01 = sum(vgen%nAtomSpec(1:iSpec)) - vgen%nAtomSpec(iSpec) + &
                nAtomSpecStack + 1
         nAtomSpecStack = nAtomSpecStack + nAtomSpec_Unit(iSpecUnit,iUnit)
         id02 = id01 + nAtomSpec_Unit(iSpecUnit,iUnit) - 1
         id04 = sum(nAtomSpec_Unit(1:iSpecUnit,iUnit))
         id03 = id04 - nAtomSpec_Unit(iSpecUnit,iUnit) + 1
         vgen%coo(1:3,id01:id02) = coo_shift(1:3,id03:id04,iUnit)
      enddo loop_order
      enddo loop_unit
   enddo loop_spec

 end subroutine gen_coo_gather

!-------------------------------------------------------------------!

end module generator
