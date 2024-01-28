module optimizer

use constants,  only :  VARLEN, dp, p_zp, twpi

use mt_serial,  only :  grnd

use mcio,       only :  mcio_isodd,     &
                        mcio_printerr

use parallel,   only :  pMast,  &
                        pp_bcast

use input,      only :  Var_Input

use coordinate, only :  Var_Coord,      &
                        coo_trans,      &
                        coo_rot_v,      &
                        coo_rot_cv,     &
                        coo_align_v

use lattice,    only :  Var_Latt

 implicit none

 type, public :: Var_Opt
   integer                                    :: nAtom
   integer,       dimension(:),   allocatable :: nAtomSpec
   integer,       dimension(:),   allocatable :: TypeAtoms
   logical                                    :: pp    = .false.
   logical                                    :: iserr = .false.
   real(kind=dp), dimension(3,3)              :: LattVec
   real(kind=dp), dimension(:,:), allocatable :: coo
 end type

 integer,                                        private :: nSpec
 integer,                                        private :: nUnit
 integer,                                        private :: nSite
 integer,                                        private :: nAtoms_max
 integer,       dimension(:),       allocatable, private :: nAtoms
 integer,       dimension(:,:),     allocatable, private :: num2id

 logical,                                        private :: opt_isinit = .false.
 logical,                                        private :: opt_iserr  = .false.
 logical,       dimension(:),       allocatable, private :: isConnect
 logical,       dimension(:,:),     allocatable, private :: isMove
 
 character(LEN=VARLEN),                          private :: LattType

 real(kind=dp), dimension(3,3),                  private :: LattVec
 real(kind=dp), dimension(:,:,:,:), allocatable, private :: coo_orig
 real(kind=dp), dimension(:,:,:,:), allocatable, private :: coo_ran
 real(kind=dp), dimension(:,:,:,:), allocatable, private :: coo_shift

contains

!-------------------------------------------------------------------!
 
 function opt_init(vinp, vcoo, vlat) result(vopt)
   implicit none
   type(Var_Input), intent(in) :: vinp
   type(Var_Coord), intent(in) :: vcoo
   type(Var_Latt),  intent(in) :: vlat

   type(Var_Opt) :: vopt
   integer       :: iSite, iUnit, ixyz, iSpec, iSpec_Unit, ii, id01, id02

   if (opt_isinit) then
      if (pMast) call mcio_printerr('optimizer initialize','optimizer is already initialized')
      vopt%iserr = .true.
      return
   endif

   nSpec      = vcoo%nSpec
   nUnit      = vinp%nUnit
   nSite      = vlat%nSite
   allocate (nAtoms(nUnit))
   nAtoms     = vinp%nAtoms
   allocate (num2id(3, nSite))
   num2id     = vlat%num2id
   nAtoms_max = maxval(nAtoms)
   allocate (isConnect(0:nUnit))
   isConnect  = vinp%isConnect
   allocate (isMove(nUnit, 0:1))
   isMove     = vcoo%isMove

   LattType   = vinp%LattType
   LattVec    = vcoo%LattVec

   if (pMast) then
      vopt%LattVec(:,1) = LattVec(:,1)*vinp%nx
      vopt%LattVec(:,2) = LattVec(:,2)*vinp%ny
      vopt%LattVec(:,3) = LattVec(:,3)*vinp%nz

      vopt%nAtom  = sum(nAtoms)*nSite
      if (isConnect(0)) vopt%nAtom = vopt%nAtom/2

      allocate (vopt%nAtomSpec(nSpec),      &
                vopt%TypeAtoms(vopt%nAtom), &
                vopt%coo(3,vopt%nAtom))

      vopt%nAtomSpec(:) = 0
      do iSpec = 1, nSpec
         do iUnit = 1, nUnit
         do iSpec_Unit = 1, vcoo%nOrderSpec_unit(0,iUnit)
            if (vcoo%nOrderSpec_unit(iSpec_Unit,iUnit) .ne. iSpec) cycle
            do ii = 1, 0, -1
               if (.not. isMove(iUnit, ii)) cycle
               vopt%nAtomSpec(iSpec) = vopt%nAtomSpec(iSpec) +                 &
                                       vcoo%nAtomSpec_unit(iSpec_Unit,iUnit) * &
                                       floor(nSite*0.5d0)
            enddo
            if (mcio_isodd(nSite)) vopt%nAtomSpec(iSpec) = vopt%nAtomSpec(iSpec) + 1
         enddo
         enddo
         id01 = sum(vopt%nAtomSpec(1:iSpec)) - vopt%nAtomSpec(iSpec) + 1
         id02 = sum(vopt%nAtomSpec(1:iSpec))
         vopt%TypeAtoms(id01:id02) = iSpec
      enddo
      call opt_coo_init(vcoo%coo_unit)
      do iSite = 1, nSite
         call opt_coo_gather(vcoo, vopt, iSite)
      enddo
      vopt%iserr = opt_iserr
   endif

   call pp_bcast(vopt%iserr)
   call opt_bcast(vopt)
 
  if (vopt%iserr) return

   opt_isinit = .true.

 end function opt_init
 
!-------------------------------------------------------------------!

 subroutine opt_fin(vopt)
   implicit none
   type(Var_Opt), intent(inout) :: vopt

   if (.not. opt_isinit) then
      if (pMast) call mcio_printerr('optimizer finalize','optimizer is not initialized')
      vopt%iserr = .true.
      return
   endif

   if (allocated (vopt%nAtomSpec)) deallocate (vopt%nAtomSpec)
   if (allocated (vopt%TypeAtoms)) deallocate (vopt%TypeAtoms)
   if (allocated (vopt%coo)) deallocate (vopt%coo)
   if (allocated (nAtoms)) deallocate (nAtoms)
   if (allocated (num2id)) deallocate (num2id)
   if (allocated (isConnect)) deallocate (isConnect)
   if (allocated (isMove)) deallocate (isMove)
   if (allocated (coo_orig)) deallocate (coo_orig)
   if (allocated (coo_ran)) deallocate (coo_ran)
   if (allocated (coo_shift)) deallocate (coo_shift)

   vopt%pp    = .false.
   opt_isinit = .false.

 end subroutine opt_fin

!-------------------------------------------------------------------!

 subroutine opt_ranwalk(rate_dis, rate_ang, iSite_in, iUnit_in)
   implicit none
   real(kind=dp),           intent(in) :: rate_dis
   real(kind=dp),           intent(in) :: rate_ang
   integer,                 intent(in) :: iSite_in
   integer,       optional, intent(in) :: iUnit_in

   integer :: iUnit

   call opt_coo_ranwalk(iSite_in, rate_dis, rate_ang, &
                        coo_orig(:,:,:,iSite_in), coo_ran(:,:,:,iSite_in))
   if (present(iUnit_in)) then
      do iUnit = 1, nUnit
         if (iUnit .eq. iUnit_in) cycle
         coo_ran(:,:,iUnit,iSite_in) = coo_orig(:,:,iUnit,iSite_in)
      enddo
   endif
   call opt_coo_shift(iSite_in, coo_ran(:,:,:,iSite_in), coo_shift(:,:,:,iSite_in))

 end subroutine opt_ranwalk

!-------------------------------------------------------------------!

 subroutine opt_coo_update(iSite_in)
   implicit none
   integer, intent(in) :: iSite_in

   coo_orig(:,:,:,iSite_in) = coo_ran(:,:,:,iSite_in)

 end subroutine

!-------------------------------------------------------------------!

 subroutine opt_coo_rollback(iSite_in)
   implicit none
   integer, intent(in) :: iSite_in

   call opt_coo_shift(iSite_in, coo_orig(:,:,:,iSite_in), coo_shift(:,:,:,iSite_in))

 end subroutine

!-------------------------------------------------------------------!

 subroutine opt_bcast(vopt)
   type(Var_Opt), intent(inout) :: vopt

   call pp_bcast(vopt%iserr)
   if (vopt%iserr) return
   if (vopt%pp) then
      if (pMast) call mcio_printerr('optimizer bcast','optimizer is already bcasted')
      vopt%iserr = .true.
      return
   endif

   call pp_bcast(vopt%nAtom)
   call pp_bcast(vopt%LattVec, 3, 3)
   if (.not. allocated (vopt%nAtomSpec)) allocate (vopt%nAtomSpec(nSpec))
   call pp_bcast(vopt%nAtomSpec, nSpec)
   if (.not. allocated (vopt%TypeAtoms)) allocate (vopt%TypeAtoms(vopt%nAtom))
   call pp_bcast(vopt%TypeAtoms, vopt%nAtom)
   if (.not. allocated (vopt%coo)) allocate (vopt%coo(3, vopt%nAtom))
   call pp_bcast(vopt%coo, 3, vopt%nAtom)

   vopt%pp = .true.

 end subroutine opt_bcast

!-------------------------------------------------------------------!
 
 subroutine opt_coo_init(coo_in)
   implicit none
   real(kind=dp), dimension(3,nAtoms_max,nUnit), intent(in) :: coo_in

   integer :: iSite1, iSite2, iUnit1, iUnit2, iAtom, nAtom
   integer :: cnt_stop

   real(kind=dp)                                :: rad, dis, val
   real(kind=dp)                                :: rate_dis, rate_ang, rate_tol
   real(kind=dp), dimension(3)                  :: vec
   real(kind=dp), dimension(nUnit)              :: radius
   real(kind=dp), dimension(3,nUnit,nSite)      :: centroid
   real(kind=dp), dimension(3,nAtoms_max,nUnit) :: coo_tmp

   allocate (coo_orig (3, nAtoms_max, nUnit, nSite), &
             coo_ran  (3, nAtoms_max, nUnit, nSite), &
             coo_shift(3, nAtoms_max, nUnit, nSite))

   do iUnit1 = 1, nUnit
      radius(:) = 0d0
      centroid(:,iUnit1,1) = 0d0
      nAtom = nAtoms(iUnit1)
      do iAtom = 1, nAtom
         centroid(:,iUnit1,1) = centroid(:,iUnit1,1) + &
                                coo_in(:,iAtom,iUnit1)/dble(nAtom)
      enddo
      do iAtom = 1, nAtom
         vec = coo_in(:,iAtom,iUnit1) - centroid(:,iUnit1,1)
         rad = sqrt(dot_product(vec, vec))
         if (rad .gt. radius(iUnit1)) radius(iUnit1) = rad
      enddo
      radius(iUnit1) = radius(iUnit1) + 1.0d0
      if (isConnect(iUnit1)) radius(iUnit1) = radius(iUnit1)*0.9d0
   enddo

   do iSite1 = 1, nSite
      coo_orig(:,:,:,iSite1) = coo_in
      rate_dis = 1d0
      rate_ang = 1d0
      rate_tol = 1d0
      cnt_stop = 0
      loop_find : do
         cnt_stop = cnt_stop + 1
         if (mod(cnt_stop,1000) .eq. 0) rate_tol = rate_tol*0.99d0
         if (rate_tol .lt. 0.85d0) then
            write (*,*) ' !!! WARNING : fail to locate UNIT in free volume : ',iSite1
            write (*,*) '               recommend to increase your cell volume'
            exit loop_find
         endif
         coo_tmp = coo_orig(:,:,:,iSite1)
         do iUnit1 = 1, nUnit
            if (.not. isMove(iUnit1,mod(iSite1,2))) cycle
            nAtom = nAtoms(iUnit1)
            val   = rate_ang*twpi*grnd()
            call coo_rot_v(nAtom, p_zp, coo_tmp(:,1:nAtom,iUnit1), val)
         enddo
         call opt_coo_ranwalk(iSite1, rate_dis, rate_ang, coo_tmp, coo_ran(:,:,:,iSite1))
         call opt_coo_shift(iSite1, coo_ran(:,:,:,iSite1), coo_shift(:,:,:,iSite1))
         do iUnit1 = 1, nUnit
            if (.not. isMove(iUnit1,mod(iSite1,2))) cycle
            nAtom = nAtoms(iUnit1)
            centroid(:,iUnit1,iSite1) = 0d0
            do iAtom = 1, nAtom
               centroid(:,iUnit1,iSite1) = centroid(:,iUnit1,iSite1) + &
                                           coo_shift(:,iAtom,iUnit1,iSite1)/dble(nAtom)
            enddo
            do iUnit2 = 1, iUnit1
               if (iUnit2 .eq. iUnit1) exit
               if (.not. isMove(iUnit2,mod(iSite1,2))) cycle
               vec = centroid(:,iUnit1,iSite1) - centroid(:,iUnit2,iSite1)
               dis = sqrt(dot_product(vec, vec))
               rad = (radius(iUnit1) + radius(iUnit2))*rate_tol
               if (dis .lt. rad) cycle loop_find
            enddo
         enddo
         if (iSite1 .eq. 1) then
            call opt_coo_update(iSite1)
            exit loop_find
         endif
         loop_site : do iSite2 = 1, iSite1-1
            do iUnit2 = 1, nUnit
            if (.not. isMove(iUnit2,mod(iSite2,2))) cycle
            do iUnit1 = 1, nUnit
            if (.not. isMove(iUnit1,mod(iSite1,2))) cycle
               vec = centroid(:,iUnit1,iSite1) - centroid(:,iUnit2,iSite2)
               dis = sqrt(dot_product(vec, vec))
               rad = (radius(iUnit1) + radius(iUnit2))*rate_tol
               val = (rad - dis)/rad
               if (val .gt. 0d0) then
                  rate_dis = 1d0
                  rate_ang = 1d0
                  if (val .gt. 0.1d0) cycle loop_find
                  rate_dis = rad - dis
                  rate_ang = val
                  coo_orig(:,:,iUnit1,iSite1) = coo_ran(:,:,iUnit1,iSite1)
                  cycle loop_find
               endif
            enddo
            enddo
         enddo loop_site
         call opt_coo_update(iSite1)
         exit loop_find
      enddo loop_find
   enddo

 end subroutine opt_coo_init

!-------------------------------------------------------------------!

 subroutine opt_coo_ranwalk(iSite_in, rate_dis, rate_ang, coo_in, coo_out)
   implicit none
   integer,                                      intent(in)  :: iSite_in
   real(kind=dp),                                intent(in)  :: rate_dis
   real(kind=dp),                                intent(in)  :: rate_ang
   real(kind=dp), dimension(3,nAtoms_max,nUnit), intent(in)  :: coo_in
   real(kind=dp), dimension(3,nAtoms_max,nUnit), intent(out) :: coo_out

   integer :: iUnit, ixyz, nAtom

   real(kind=dp)               :: len_vec, rot_ang
   real(kind=dp), dimension(3) :: dis_vec, rot_vec

   do iUnit = 1, nUnit
      if (.not. isMove(iUnit,mod(iSite_in,2))) cycle
      nAtom   = nAtoms(iUnit)
      rot_ang = rate_ang*twpi*grnd()
      coo_out(1:3,1:nAtom,iUnit) = coo_in(1:3,1:nAtom,iUnit)
      if (isConnect(iUnit)) then
         call coo_rot_v(nAtom, p_zp, coo_out(:,1:nAtom,iUnit), rot_ang)
      else
         do ixyz = 1, 3
            rot_vec(ixyz) = grnd()-0.5d0
            dis_vec(ixyz) = grnd()-0.5d0
         enddo
         call coo_rot_cv(nAtom, rot_vec, coo_out(:,1:nAtom,iUnit), rot_ang)
         len_vec = sqrt(dot_product(dis_vec, dis_vec))
         dis_vec(:) = rate_dis*dis_vec(:)/len_vec
         call coo_trans(nAtom, dis_vec, coo_out(:,1:nAtom,iUnit))
      endif
   enddo

 end subroutine opt_coo_ranwalk

!-------------------------------------------------------------------!

 subroutine opt_coo_shift(iSite_in, coo_in, coo_out)
   implicit none
   integer,                                      intent(in)  :: iSite_in
   real(kind=dp), dimension(3,nAtoms_max,nUnit), intent(in)  :: coo_in
   real(kind=dp), dimension(3,nAtoms_max,nUnit), intent(out) :: coo_out

   integer                   :: iUnit, ixyz, iAtom, icnt, nAtom
   integer, dimension(3,2)   :: lidx
   integer, dimension(3,0:1) :: lidx_pbc

   real(kind=dp), dimension(3)                  :: dir_conn
   real(kind=dp), dimension(3,nUnit)            :: LattLoc
   real(kind=dp), dimension(3,2)                :: M

   lidx = opt_lidx_get(iSite_in)

   if (isConnect(0)) then
      lidx_pbc(:,1) = lidx(:,1)
      lidx_pbc(:,0) = lidx(:,2)
      do ixyz = 1, 3
         if ((lidx(ixyz,1) - lidx(ixyz,2)) .gt.  1) lidx_pbc(ixyz,0) = lidx(ixyz,1) + 1
         if ((lidx(ixyz,1) - lidx(ixyz,2)) .lt. -1) lidx_pbc(ixyz,1) = lidx(ixyz,2) + 1
      enddo
   endif

   icnt = 0
   LattLoc(:,:) = 0d0

   select case(LattType)
   case('S')
      M(:,:) = 1d0
      if (isConnect(0)) M(:,1) = 0.5d0
   case('FC')
      M(:,:) = 0.5d0
      if (isConnect(0)) M(:,1) = 0.25d0
   end select

   if (isConnect(0)) then
      do iUnit = 1, nUnit
         if (.not. isConnect(iUnit)) icnt = icnt + 1
         if (.not. isMove(iUnit,mod(iSite_in,2))) cycle
         do ixyz = 1, 3
            if (isConnect(iUnit)) then
               LattLoc(:,iUnit) = LattLoc(:,iUnit) + LattVec(:,ixyz) * M(ixyz,1) * &
                                  dble(lidx_pbc(ixyz,1) + lidx_pbc(ixyz,0) - 2)
            else
               LattLoc(:,iUnit) = LattLoc(:,iUnit) + LattVec(:,iXYZ) * M(ixyz,2) * &
                                  dble(lidx_pbc(ixyz,(mod(icnt,2))) - 1)
            endif
         enddo
      enddo
      dir_conn(:) = lidx_pbc(:,0) - lidx_pbc(:,1)
   elseif (.not. isConnect(0)) then
      do ixyz = 1, 3
         LattLoc(:,1) = LattLoc(:,1) + LattVec(:,ixyz) * M(ixyz,1) * dble(lidx(ixyz,1)-1)
      enddo
   endif

   do iUnit = 1, nUnit
      if (.not. isMove(iUnit,mod(iSite_in,2))) cycle
      nAtom = nAtoms(iUnit)
      coo_out(:,1:nAtom,iUnit) = coo_in(:,1:nAtom,iUnit)
      if (      isConnect(0)) call coo_align_v(nAtom, p_zp, dir_conn, coo_out(:,1:nAtom,iUnit))
      if (.not. isConnect(0)) LattLoc(:,iUnit) = LattLoc(:,1)
      call coo_trans(nAtom, LattLoc(:,iUnit), coo_out(:,1:nAtom,iUnit))
   enddo

 end subroutine opt_coo_shift

!-------------------------------------------------------------------!

 subroutine opt_coo_gather(vcoo, vopt, iSite_in)
   implicit none
   type(Var_Coord), intent(in)    :: vcoo
   type(Var_Opt),   intent(inout) :: vopt
   integer,         intent(in)    :: iSite_in

   integer :: iSpec, iUnit, iSpecUnit, ii, stack, id01, id02, id03, id04
   integer, dimension(2) :: nSpAtom

   stack = iSite_in - 1
   if (isConnect(0)) stack = floor(stack*0.5d0)*2
   loop_spec : do iSpec = 1, nSpec
      nSpAtom(1) = vopt%nAtomSpec(iSpec)*stack/nSite
      nSpAtom(2) = 0
      loop_unit : do iUnit = 1, nUnit
         loop_order : do iSpecUnit = 1, vcoo%nOrderSpec_unit(0,iUnit)
            if (vcoo%nOrderSpec_unit(iSpecUnit,iUnit) .ne. iSpec) cycle loop_order
            ii = 2
            do while (.true.)
               ii = ii - 1
               if (ii .lt. 0) exit
               if (.not. isMove(iUnit,ii)) cycle
               if (.not. isConnect(0)) ii = ii - 2
               id01 = sum(vopt%nAtomSpec(1:iSpec)) - vopt%nAtomSpec(iSpec) + sum(nSpAtom) + 1
               nSpAtom(2) = nSpAtom(2) + vcoo%nAtomSpec_unit(iSpecUnit,iUnit)
               if (isConnect(0) .and. (mod(iSite_in,2) .ne. ii)) cycle
               id02 = id01 + vcoo%nAtomSpec_unit(iSpecUnit,iUnit) - 1
               id04 = sum(vcoo%nAtomSpec_unit(1:iSpecUnit,iUnit))
               id03 = id04 - vcoo%nAtomSpec_unit(iSpecUnit,iUnit) + 1
               vopt%coo(1:3,id01:id02) = coo_shift(1:3,id03:id04,iUnit,iSite_in)
! write (*,'(A,2I6,A,4I4)') vcoo%NameAtom(iSpec), id01, id02, ' // ', id03, id04, iUnit, iSite_in
            enddo
         enddo loop_order
      enddo loop_unit
   enddo loop_spec

 end subroutine opt_coo_gather

!-------------------------------------------------------------------!

 function opt_lidx_get(iSite_in) result (lidx_out)
   implicit none
   integer, intent(in) :: iSite_in

   integer, dimension(3,2) :: lidx_out

   if (.not. isConnect(0)) then
      lidx_out(:,1) = num2id(:,iSite_in)
      lidx_out(:,2) = lidx_out(:,1)
   elseif (mcio_isodd(iSite_in)) then
      lidx_out(:,1) = num2id(:,iSite_in)
      lidx_out(:,2) = num2id(:,iSite_in + 1)
   else
      lidx_out(:,1) = num2id(:,iSite_in - 1)
      lidx_out(:,2) = num2id(:,iSite_in)
   endif

 end function opt_lidx_get
   
!-------------------------------------------------------------------!
 
end module optimizer
