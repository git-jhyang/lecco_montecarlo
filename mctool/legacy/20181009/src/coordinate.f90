module coordinate

use constants,  only :  dp, LINELEN,    &
                        p_o, p_zp,      &
                        CHKOK,          &
                        CHKERR,         &
                        CHKEOF

use util,       only :  dir_to_cart,    &
                        cart_to_dir,    &
                        cart_to_sph,    &
                        sph_to_cart,    &
                        cross_product,  &
                        assert_data,    &
                        angle, rot_vec

use iotool,     only :  io_getunit,     &
                        io_isodd,       &
                        io_printhead,   &
                        io_printerr

use vasp,       only :  poscar,         &
                        poscar_read,    &
                        poscar_write

use parallel,   only :  pMast,          &
                        pp_bcast

use input,      only :  Var_Input

implicit none
save

!-------------------------------------------------------------------!
!       Var_Coord       : variables of coordination
!-------------------------------------------------------------------!
!       nSpec          : number of total atom species w/o duplication
!       nAtomSpec_unit      : number of atoms of each species of each unit
!       nOrderSpec_unit     : sorting order of each atom block
!       
!       NameAtom         : name of each atom species
!
!       iserr           : .true. if there is an error
!       ismove(n,0:1)   : .false. is dummy coordinate. not counted
!       
!       lvec            : lattice vector
!       coo_unit       : cartesian coordinate of each unit
!       chg_units       : point charge of each unit (for amber)
!-------------------------------------------------------------------!

type, public :: Var_Coord
   integer                                          :: nSpec
   integer,          dimension(:,:),    allocatable :: nAtomSpec_unit
   integer,          dimension(:,:),    allocatable :: nOrderSpec_unit

   character(LEN=3), dimension(:),      allocatable :: NameAtom

   logical                                          :: iserr = .false.
   logical                                          :: pp    = .false.
   logical,          dimension(:,:),    allocatable :: isMove

   real(kind=dp),    dimension(3,3)                 :: LattVec
   real(kind=dp),    dimension(:,:,:),  allocatable :: coo_unit
   real(kind=dp),    dimension(:,:,:),  allocatable :: force_unit
end type Var_Coord 

integer,                            private :: nUnit
integer,                            private :: nSpecs_max
integer,                            private :: nAtoms_max
integer, dimension(:), allocatable, private :: nSpecs
integer, dimension(:), allocatable, private :: nAtoms

logical, private :: coo_isinit = .false.
logical, private :: coo_iserr  = .false.

!-------------------------------------------------------------------!
!       functions and subroutines
!-------------------------------------------------------------------!
!       = functions = 
!       coo_read(vinp)  : read coodination file using 'vinp' and return
!                         'vcoo'
!       coo_chkdim(line,dim)  : when read line, check number of variable
!                       error returns 'true'
!       coo_chkcart(u,line,iserr)  : check the poscar is whether 'dirct'
!                       or 'cartesian'. if cart, return 'true'.
!
!       = subroutines=
!       pp_bcast_vcoo(vinp,vcoo)  : broadcast 'vcoo' to all processor
!       coo_alloc(vcoo,n,nSpec,nAtoms)  : allocate 'vcoo'
!                       variables using nUnit 'n', 'nSpec', 
!                       and 'nAtoms'
!       coo_getposcar(f, n, m, lvec, names, nums, coo, chg)  : 
!                       read vasp pocscar with filename 'f' and return
!                       nSpec 'n', tnum_at 'm', 'lvec', 'names', 'nums'
!                       'coo', and 'chg'
!       coo_getorder(vcoo, n, m, name_units)  : get 'nOrderSpec_unit'
!                       using nSpecs_max 'n', nUnit 'm', and 'name_unit'
!       coo_printinfo(vinp, vcoo)  : print information about 'vcoo'
!       coo_trans(n,vec,coo)  : translate 'n' atoms coordinate 'coo'
!                       with vector 'vec'
!       coo_rot_v(n,vec,coo,ang)  : rotate 'n' atoms coordinate 'coo'
!                       about vector 'vec' from origin with angle 'ang'
!       coo_rot_cv(n,vec,coo,ang)  : rotate 'n' atoms coordinate 'coo'
!                       about vector 'vec' from its centroid with angle 'ang'
!       coo_align_v(n,vec_in,vec_out,coo)  : align 'n' atoms coordinate 
!                       'coo' from 'vec_in' to 'vec_out'
!-------------------------------------------------------------------!

public ::       coo_init,       coo_fin,        &
                coo_getposcar,  coo_writeposcar,&
                coo_trans,      coo_rot_v,      &
                coo_rot_cv,     coo_align_v

private ::      coo_bcast,      coo_defaults,   &
                coo_print_info, coo_getorder

contains 

!-------------------------------------------------------------------!

 function coo_init(vinp) result(vcoo)
   implicit none
   type(Var_Input), intent(in) :: vinp
   type(Var_Coord)             :: vcoo

   integer :: iUnit, iCnt

   real(kind=dp),    dimension(:,:,:), allocatable :: LattVec_unit
   character(LEN=2), dimension(:,:),   allocatable :: NameAtom_unit

   if (coo_isinit) then
      if (pMast) call io_printerr('coord initalize','coord is already initialized')
      vcoo%iserr = .true.
      return
   endif

   nUnit      = vinp%nUnit
   allocate (nSpecs(nUnit), nAtoms(nUnit))
   nSpecs     = vinp%nSpecs
   nSpecs_max = maxval(nSpecs)
   nAtoms     = vinp%nAtoms
   nAtoms_max = maxval(nAtoms)

   allocate (vcoo%nAtomSpec_unit(nSpecs_max, nUnit),    &
             vcoo%nOrderSpec_unit(0:nSpecs_max, nUnit), &
             vcoo%isMove(nUnit, 0:1))

   if (pMast) then
      allocate (vcoo%coo_unit(3, nAtoms_max, nUnit), &
                vcoo%force_unit(3, nAtoms_max, nUnit))
      call coo_defaults(vcoo)
      allocate(LattVec_unit(3, 3, nUnit),       &
               NameAtom_unit(nSpecs_max, nUnit))
      NameAtom_unit(:,:) = 'xx'
      do iUnit = 1, nUnit
         call coo_getposcar(vinp%fname_coo_in(iUnit), nSpecs(iUnit),    &
                            nAtoms(iUnit), LattVec_unit(:,:,iUnit),     &
                            NameAtom_unit(1:nSpecs(iUnit),iUnit),       &
                            vcoo%nAtomSpec_Unit(1:nSpecs(iUnit),iUnit), &
                            vcoo%coo_unit(:,1:nAtoms(iUnit),iUnit),     &
                            vcoo%force_unit(:,1:nAtoms(iUnit),iUnit))
         call coo_trans(nAtoms(iUnit), -vinp%unit_cen(:,iUnit), &
                        vcoo%coo_unit(:,1:nAtoms(iUnit),iUnit))
         call coo_align_v(nAtoms(iUnit), vinp%unit_dir(:,iUnit), p_zp, &
                          vcoo%coo_unit(:,1:nAtoms(iUnit),iUnit))
      enddo
      if (vinp%isConnect(0)) then
         iCnt = 0
         do iUnit = 1, nUnit
            if (.not. vinp%isConnect(iUnit)) then
               iCnt = iCnt + 1
               if (.not. io_isodd(iCnt)) vcoo%isMove(iUnit,1) = .false.
            endif
            vcoo%isMove(iUnit,0) = .not. vcoo%isMove(iUnit,1)
         enddo
      endif
      vcoo%nOrderSpec_unit(0,1:nUnit) = nSpecs(1:nUnit)
      call coo_getorder(vcoo, NameAtom_unit)
      vcoo%iserr = coo_iserr
      vcoo%LattVec = LattVec_unit(:,:,1)
      deallocate (LattVec_unit, NameAtom_unit)
   endif

   call coo_bcast(vcoo)
   if (pMast) call coo_print_info(vinp, vcoo)
   coo_isinit = .true.

 end function coo_init

!-------------------------------------------------------------------!

 subroutine coo_fin(vcoo)
   implicit none
   type(Var_Coord), intent(inout) :: vcoo

   if (.not. coo_isinit) then
      if (pMast) call io_printerr('coord finalize','coord is not initalized')
      vcoo%iserr = .true.
      return
   endif

   if (allocated (vcoo%nAtomSpec_unit)) deallocate (vcoo%nAtomSpec_unit)
   if (allocated (vcoo%nOrderSpec_unit)) deallocate (vcoo%nOrderSpec_unit)
   if (allocated (vcoo%NameAtom)) deallocate (vcoo%NameAtom)
   if (allocated (vcoo%isMove)) deallocate (vcoo%isMove)
   if (allocated (vcoo%coo_unit)) deallocate (vcoo%coo_unit)
   if (allocated (vcoo%force_unit)) deallocate (vcoo%force_unit)
   if (allocated (nSpecs)) deallocate (nSpecs)
   if (allocated (nAtoms)) deallocate (nAtoms)

   coo_isinit = .false.
   vcoo%pp    = .false.

 end subroutine coo_fin

!-------------------------------------------------------------------!
!          broadcast Coordinate variables
!-------------------------------------------------------------------!

 subroutine coo_bcast(vcoo)
   implicit none
   type(Var_Coord), intent(inout) :: vcoo

   integer :: ii

   call pp_bcast(vcoo%iserr) 
   if (vcoo%iserr) return
   if (vcoo%pp) then
      if (pMast) call io_printerr('coord bcast','coord is already bcasted')
      vcoo%iserr = .true.
      return
   endif

   call pp_bcast(vcoo%nSpec)
   if (.not. allocated (vcoo%nAtomSpec_unit)) allocate (vcoo%nAtomSpec_unit(nSpecs_max, nUnit))
   call pp_bcast(vcoo%nAtomSpec_unit, nSpecs_max, nUnit)
   if (.not. allocated (vcoo%nOrderSpec_unit)) allocate (vcoo%nOrderSpec_unit(0:nSpecs_max, nUnit))
   call pp_bcast(vcoo%nOrderSpec_unit(0:nSpecs_max,1:nUnit), nSpecs_max+1, nUnit)
   if (.not. allocated (vcoo%NameAtom)) allocate (vcoo%NameAtom(vcoo%nSpec))
   call pp_bcast(vcoo%NameAtom, vcoo%nSpec)
   if (.not. allocated (vcoo%isMove)) allocate (vcoo%isMove(nUnit, 0:1))
   do ii = 0, 1
      call pp_bcast(vcoo%isMove(:,ii), nUnit)
   enddo
   call pp_bcast(vcoo%LattVec, 3, 3)
!   if (.not. allocated(vcoo%coo_unit)) &
!        allocate(vcoo%coo_unit(3, maxval(vinp%nAtoms), vinp%nUnit))
!   call pp_bcast(vcoo%coo_unit, 3, maxval(vinp%nAtoms), vinp%nUnit)
!   if (.not. allocated(vcoo%chg_units)) &
!        allocate(vcoo%chg_units(maxval(vinp%nAtoms), vinp%nUnit))
!   call pp_bcast(vcoo%chg_units, maxval(vinp%nAtoms), vinp%nUnit)

   vcoo%pp = .true.

 end subroutine coo_bcast

!-------------------------------------------------------------------!

 subroutine coo_defaults(vcoo)
   implicit none
   
   type(Var_Coord), intent(inout) :: vcoo

   vcoo%nAtomSpec_unit  = 0
   vcoo%nOrderSpec_unit = 0
   vcoo%isMove          = .true.
   vcoo%coo_unit        = 0d0
!   vcoo%force_unit        = 0d0
   
 end subroutine coo_defaults

!-------------------------------------------------------------------!

 subroutine coo_print_info(vinp, vcoo)
   implicit none
   type(Var_Input), intent(in) :: vinp
   type(Var_Coord), intent(in) :: vcoo

   integer :: iUnit, iSpec, iAtom, iCnt1, iCnt2

   if (coo_iserr) return

   iCnt1 = 0
   iCnt2 = 0
   call io_printhead('Coordinate information')
   do iUnit = 1, nUnit
      write (*,210) ' Coordination No. ',iUnit,' -----------------------------------------------------'
      write (*,205) '   Coordination filename         : ', trim(adjustl(vinp%fname_coo_in(iUnit)))
      write (*,230) '   Connects two latt.points      : ', vinp%isConnect(iUnit)
      write (*,225) '   Number of atom species        : ', nSpecs(iUnit)
      write (*,200,advance = 'no') '   Atom types                    : '
      do iSpec = 1, nSpecs(iUnit)
         write (*,202,advance = 'no') adjustr(vcoo%NameAtom(vcoo%nOrderSpec_unit(iSpec,iUnit)))
      enddo
      write (*,*)
      write (*,200,advance = 'no') '   Number of atoms               : '
      do iSpec = 1, nSpecs(iUnit)
         write (*,220,advance = 'no') vcoo%nAtomSpec_unit(iSpec,iUnit)
      enddo
      write (*,*)
      write (*,200) '   Cartesian coordinates         : '
      do iAtom = 1, nAtoms(iUnit)
         write (*,240) vcoo%coo_unit(:,iAtom,iUnit)!,' ( ',vcoo%force_unit(ij,ii),' ) '
      enddo
      write (*,*)
   enddo

200 format(1x,A)
202 format(A4)
205 format(1x,A,A)
210 format(1x,A,I2,A)
220 format(I4)
225 format(1x,A,I4)
230 format(1x,A,L4)
240 format(4x,3F15.7)!,A,F6.3,A)

 end subroutine coo_print_info

!-------------------------------------------------------------------!

 subroutine coo_getorder(vcoo, NameAtom_unit)
   implicit none
   type(Var_Coord),                               intent(inout) :: vcoo
   character(LEN=2), dimension(nSpecs_max,nUnit), intent(in)    :: NameAtom_unit

   logical :: isovlp
   integer :: iUnit, iSpec, iOrderSpec, iCnt
   character(LEN=2), dimension(:), allocatable :: NameAtom_tmp

   if (coo_iserr) return

   vcoo%nSpec = vcoo%nOrderSpec_unit(0,1)
   allocate (vcoo%NameAtom(vcoo%nSpec))
   
   vcoo%NameAtom = NameAtom_unit(1:vcoo%nSpec,1)
   vcoo%nOrderSpec_unit(1:vcoo%nSpec,1) = [1:vcoo%nSpec]

   if (nUnit .eq. 1) return

   do iUnit = 2, nUnit
      do iOrderSpec = 1, vcoo%nOrderSpec_unit(0,iUnit)
         isovlp = .false.
         do iSpec = 1, vcoo%nSpec
            if (vcoo%NameAtom(iSpec) .eq. NameAtom_unit(iOrderSpec,iUnit)) then
               isovlp = .true.
               exit
            endif
         enddo
         if (isovlp) then
            vcoo%nOrderSpec_unit(iOrderSpec,iUnit) = iSpec
         else
            vcoo%nOrderSpec_unit(iOrderSpec,iUnit) = vcoo%nSpec + 1
            allocate (NameAtom_tmp(vcoo%nSpec))
            NameAtom_tmp = vcoo%NameAtom
            deallocate (vcoo%NameAtom)
            allocate (vcoo%NameAtom(vcoo%nSpec+1))
            vcoo%NameAtom(1:vcoo%nSpec) = NameAtom_tmp
            vcoo%NameAtom(vcoo%nSpec+1) = NameAtom_unit(iOrderSpec,iUnit)
            deallocate (NameAtom_tmp)
            vcoo%nSpec = vcoo%nSpec + 1
         endif
      enddo
   enddo

 end subroutine coo_getorder

!-------------------------------------------------------------------!

 subroutine coo_getposcar(f_name, nSpec, nAtom, LattVec, NameAtom, &
                          nAtomSpec, coo_out, force_out)
   implicit none
   character(LEN=*),                               intent(in)  :: f_name
   integer,                                        intent(in)  :: nSpec
   integer,                                        intent(in)  :: nAtom
   real(kind=dp),    dimension(3,3),               intent(out) :: LattVec
   character(LEN=2), dimension(nSpec),             intent(out) :: NameAtom
   integer,          dimension(nSpec),             intent(out) :: nAtomSpec
   real(kind=dp),    dimension(3,nAtom),           intent(out) :: coo_out
   real(kind=dp),    dimension(3,nAtom), optional, intent(out) :: force_out

   type(poscar) :: DatPos
   integer      :: u, iAtom

   u = io_getunit()

   open(unit = u, file=f_name, status='old', action='read')

   call poscar_read(u, DatPos)

   if (.not. DatPos%isCart) then
      do iAtom = 1, nAtom
         call dir_to_cart(DatPos%LattVec, DatPos%coo(:,iAtom))
      enddo
   endif

   isErr = assert_data(nSpec, DatPos%nSpec)
   if (isErr) call io_printerr('coordinate - input error','nSpec')
   isErr = assert_data(nAtom, DatPos%nAtom)
   if (isErr) call io_printerr('coordinate - input error','nAtom')

   LattVec   = DatPos%LattVec
   NameAtom  = DatPos%NameAtoms
   nAtomSpec = DatPos%nAtomSpec
   coo_out   = DatPos%coo
   force_out = DatPos%force

   deallocate (DatPos%NameAtom, DatPos%nAtomSpec, &
               DatPos%coo,      DatPos%force)
   close(u)

 end subroutine coo_getposcar

!-------------------------------------------------------------------!

 subroutine coo_writeposcar(f_name, nSpec, nAtom, f_comm, LattVec, &
                            NameAtom, nAtomSpec, isCart, coo_in,   &
                            force_in, islong_in)
   implicit none
   character(LEN=*),                               intent(in) :: f_name
   integer,                                        intent(in) :: nSpec
   integer,                                        intent(in) :: nAtom
   character(LEN=*),                               intent(in) :: f_comm
   real(kind=dp),    dimension(3,3),               intent(in) :: LattVec
   character(LEN=*), dimension(nSpec),             intent(in) :: NameAtom
   integer,          dimension(nSpec),             intent(in) :: nAtomSpec
   logical,                              optional, intent(in) :: isCart
   real(kind=dp),    dimension(3,nAtom),           intent(in) :: coo_in
   real(kind=dp),    dimension(3,nAtom), optional, intent(in) :: force_in
   logical,                              optional, intent(in) :: islong_in

   type(poscar) :: DatPos
   integer      :: u

   u = io_getunit()
   open (unit=u, file=f_name, action='write')

   allocate (DatPos%NameAtoms(nSpec), DatPos%nAtomSpec(nSpec), &
             DatPos%coo(3,nAtom),     DatPos%force(3,nAtom))

   DatPos%Des        = f_comm
   DatPos%LattVec    = LattVec
   DatPos%nSpec      = nSpec
   DatPos%nAtom      = nAtom
   DatPos%NameAtoms  = NameAtom
   DatPos%nAtomSpec  = nAtomSpec
   DatPos%isCart     = .true.
   if (present(isCart)) DatPos%isCart = isCart
   DatPos%coo        = coo_in
   DatPos%force(:,:) = 0d0
   if (present(force_in)) DatPos%force  = force_in 

   call poscar_write(u, DatPos, islong_in=islong_in)

   deallocate (DatPos%NameAtom, DatPos%nAtomSpec, &
               DatPos%coo,      DatPos%force)
   close (u)

 end subroutine coo_writeposcar

!-------------------------------------------------------------------!

 subroutine coo_trans(nAtom, vec, coo)
   implicit none
   integer,                           intent(in)    :: nAtom
   real(kind=dp), dimension(3),       intent(in)    :: vec
   real(kind=dp), dimension(3,nAtom), intent(inout) :: coo

   integer :: iAtom

   do iAtom = 1, nAtom
      coo(:,iAtom) = coo(:,iAtom) + vec
   enddo

 end subroutine coo_trans

!-------------------------------------------------------------------!

 subroutine coo_rot_v(nAtom, vec, coo, ang)
   implicit none
   integer,                           intent(in)    :: nAtom
   real(kind=dp),                     intent(in)    :: ang
   real(kind=dp), dimension(3),       intent(in)    :: vec
   real(kind=dp), dimension(3,nAtom), intent(inout) :: coo

   integer :: iAtom, ixyz

   do iAtom = 1, nAtom
      call rot_vec(p_o, vec, coo(:,iAtom), ang)
      do ixyz = 1, 3
         if (abs(coo(ixyz,iAtom)) .lt. 1.0d-14) coo(ixyz,iAtom) = 0d0
      enddo
   enddo

 end subroutine coo_rot_v

!-------------------------------------------------------------------!

 subroutine coo_rot_cv(nAtom, vec, coo, ang)
   implicit none
   integer,                           intent(in)    :: nAtom
   real(kind=dp),                     intent(in)    :: ang
   real(kind=dp), dimension(3),       intent(in)    :: vec
   real(kind=dp), dimension(3,nAtom), intent(inout) :: coo

   integer                     :: iAtom, ixyz
   real(kind=dp), dimension(3) :: centroid, coo_tmp
   real(kind=dp), dimension(3) :: coo_sph

   centroid(:) = 0d0

   do iAtom = 1, nAtom
      centroid(:) = centroid(:) + coo(:,iAtom)
   enddo
   centroid(:) = centroid(:)/dble(nAtom)

   do iAtom = 1, nAtom
      coo_tmp(:) = coo(:,iAtom) - centroid(:)
      call rot_vec(p_o, vec, coo_tmp, ang)
      coo(:,iAtom) = coo_tmp(:) + centroid(:)
      do ixyz = 1, 3
         if (abs(coo(ixyz,iAtom)) .lt. 1.0d-14) coo(ixyz,iAtom) = 0d0
      enddo
   enddo

 end subroutine coo_rot_cv

!-------------------------------------------------------------------!

 subroutine coo_align_v(n, vec_in, vec_out, coo)
   implicit none
   integer,                       intent(in)    :: n ! tnum_at
   real(kind=dp), dimension(3),   intent(in)    :: vec_in
   real(kind=dp), dimension(3),   intent(in)    :: vec_out
   real(kind=dp), dimension(3,n), intent(inout) :: coo

   integer                     :: ii
   real(kind=dp)               :: ang, len_vec
   real(kind=dp), dimension(3) :: vec_ver,   vec_tmp

   vec_ver = cross_product(vec_in, vec_out)
   ang     = angle(vec_in, vec_out)
   len_vec = sqrt(dot_product(vec_ver, vec_ver))

   if (len_vec .ne. 0d0) vec_ver(:) = vec_ver(:)/len_vec
   if (len_vec .eq. 0d0) then
      if (ang .eq. 0d0) then
         return
      else
         coo(:,:) = -coo(:,:)
         return
      endif
   endif

   do ii = 1, n
      vec_tmp = coo(:,ii)
      call rot_vec(p_o, vec_ver, vec_tmp, ang)
      coo(:,ii) = vec_tmp
   enddo

 end subroutine coo_align_v

end module coordinate
