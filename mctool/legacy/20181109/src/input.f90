module input

use constants,  only :  dp, pi,         &
                        VARLEN,         &
                        PATHLEN,        &
                        LINELEN,        &
                        CHKOK,          &
                        CHKERR,         &
                        CHKEOF  

use iotool,     only :  io_getunit,     &
                        io_getcapital,  &
                        io_assert_file, &
                        io_findline,    &
                        io_findstr,     &
                        io_printhead,   &
                        io_printerr,    &
                        io_checkdim,    &
                        io_linesplit

use parallel,   only :  pMast,   &
                        pp_bcast

implicit none
save

!-------------------------------------------------------------------!
!       Var_Input       : variables of input
!-------------------------------------------------------------------!
!       = integer =
!       num_units       : 'NUM.UNITS'   : number of units
!       nOptSite        : 'OPT.SITE'    : maximum number of step 
!                                         for site selection
!       nOptMig         : 'OPT.MIG'     : maximum number of step 
!                                         for migration
!       nOptTry         : 'OPT.TRY'     : maximum number of step 
!                                         for trial coordinate
!       nx, ny, and nz  : 'LATT.DIM'    : dimension of super cell
!       nSiteTotal      : not given     : number of total possible 
!                                         occupation sites
!       nSite           : not given     : number of total occupied 
!                                         sites
!       seed            : 'SEED.MASTER' : input seed number for mt_stream
!       nSpecs(n)       : 'N.SPECIES'   : number of atom species for each unit
!       nAtoms(n)       : 'N.ATOM'      : number of total atoms for each unit
!
!       = logical = 
!       iserr           : not given     : 'true' if error occured
!       lparam          : 'L.PARAM'     : if 'true', rewrite parameter file
!       isTrace         : 'L.COORD'     : if 'true', write trajectory
!       lseed           : not given     : if 'SEED.MASTER' given, become 'true'
!       isConn(n)       : 'L.CONNECT'   : if 'true', connect unit 'n' 
!                                         connects two lattice points
!       isConn(0)       : not given     : 'true' if one of isConn(n) is 'true' 
!
!       = character =
!       method          : 'ENE.METHOD'  : energy calculation method
!                                         'AENET' for AENET
!       LattType        : 'LATT.TYPE'   : lattice type
!                                         'S' for simple lattice, default
!                                         'BC', 'BCL', 'BCC' for base centered
!                                         'FC', 'FCL', 'FCC' for face centered
!                                         'HCP', 'CP' for hexagonal closed packing
!       fname_input     : given at run  : input filename, 'opt.in' is default
!       fname_coo_out   : 'FILE.COO.OUT': coordinate trajectory and output filename 
!       fname_pot       : 'FILE.POT'    : potential filename
!       fname_coo_in(n) : 'FILE.COO.IN' : coordinate input filename
!
!       = real =
!       temperature     : 'TEMPERATURE' : temperature for possibility calc.
!       rate_occ        : 'RATE.OCC'    : lattice occupation ratio
!       unit_cen(3,n)   : 'UNIT.CENTER' : center of eachunit in cartesian
!       unit_dir(3,n)   : 'UNIT.DIRECTION' : main axis direction of 
!                                            each unit 
!-------------------------------------------------------------------!

type :: Var_Input
! common
   integer                                             :: nUnit
   integer                                             :: seed = -999
   integer,                dimension(:),   allocatable :: nSpecs
   integer,                dimension(:),   allocatable :: nAtoms
   logical                                             :: iserr = .false.
   logical                                             :: lseed
   logical,                dimension(:),   allocatable :: isConn
   character(LEN=VARLEN)                               :: mode
   character(LEN=PATHLEN)                              :: fname_in
   character(LEN=PATHLEN)                              :: fname_coo_out
   character(LEN=PATHLEN), dimension(:),   allocatable :: fname_coo_in
   real(kind=dp),          dimension(3,3)              :: LattVec
   real(kind=dp),          dimension(:,:), allocatable :: unit_cen
   real(kind=dp),          dimension(:,:), allocatable :: unit_dir

! optimizer
   integer                :: nOptSite
   integer                :: nOptMig 
   integer                :: nOptTry
   integer                :: nx, ny, nz
   integer                :: nSiteTotal
   integer                :: nSite
   integer                :: LattType
   logical                :: pp    = .false.
   logical                :: isTrace
   logical                :: selective_site
   logical                :: selective_unit
   logical                :: ispbc
   logical                :: lforce
   character(LEN=VARLEN)  :: method
   character(LEN=PATHLEN) :: fname_amber
   real(kind=dp)          :: temperature
   real(kind=dp)          :: rate_occ
   real(kind=dp)          :: rate_dis
   real(kind=dp)          :: rate_ang

! generator
   integer                       :: nGridRotG
   integer,       dimension(3)   :: nGridRotL
   integer,       dimension(3)   :: nGridTrans
   real(kind=dp)                 :: AngRotG
   real(kind=dp), dimension(3)   :: AngRotL
   real(kind=dp), dimension(2,3) :: DisTrans
 end type Var_Input

logical, private :: inp_isinit = .false.
logical, private :: inp_iserr  = .false.

!-------------------------------------------------------------------!
!       functions and subroutines
!-------------------------------------------------------------------!
!       = functions = 
!       inp_read(f)     : read input file 'f' and return 'vinp'
!      
!       = subroutines=
!       inp_readnetwork(n,names,networks,stat) : 
!                       read AENET network files and match to 'name_at'
!                       after coordination loading. if error occured, 
!                       'stat' is none-zero
!       pp_bcast_vinp(vinp)     : broadcast 'vinp' to all processor
!       inp_defaults(vinp)      : load default values of 'vinp'
!       inp_printinfo(vinp)     : print information about 'vinp'
!       inp_extend(vinp)        : when connected and non-connected 
!                       molecules exist at the same time, case of 
!                       isConn(0) and .not. isConn(n), duplicate 
!                       non-connected molecule with their variables
!       inp_read_fnames(u, chr, f, n, isfile)  : 
!                       read 'n' filenames 'f' of input coordinates.
!                       if one file does not exists, isfile become false
!       inp_read_vecs(u, chr, vec, n)       : read 'n' vectors 'vec'
!-------------------------------------------------------------------!

public ::       inp_init,       inp_fin,  &
                inp_read_network

private ::      inp_defaults,      &
                inp_read_opt,      &
                inp_read_gen,      &
                inp_bcast,         &
                inp_printinfo_opt, &
                inp_printinfo_gen, &
                inp_extend,        &
                inp_read_fnames,   &
                inp_read_vecs

!-------------------------------------------------------------------!
!       inp_read_val    : read input value
!-------------------------------------------------------------------!
!       inp_read_i1(u,k,v)      : read 1 integer 'v' with keyword 'k'
!                                 in file indexed by 'u'
!       inp_read_in(u,k,v,n)    : read 'n'-dim integer array of 'v'
!                                 with keyword 'k' in file indexed by 'u'
!       inp_read_r1(u,k,v)      : read 1 real 'v' with keyword 'k'
!                                 in file indexed by 'u'
!       inp_read_rn(u,k,v,n)    : read 'n'-dim real array of 'v'
!                                 with keyword 'k' in file indexed by 'u'
!       inp_read_c1(u,k,v)      : read 1 character 'v' with keyword 'k'
!                                 in file indexed by 'u'
!       inp_read_cn(u,k,v,n)    : read 'n'-dim character array of 'v'
!                                 with keyword 'k' in file indexed by 'u'
!       inp_read_l1(u,k,v)      : read 1 logical 'v' with keyword 'k'
!                                 in file indexed by 'u'
!       inp_read_ln(u,k,v,n)    : read 'n'-dim logical array of 'v'
!                                 with keyword 'k' in file indexed by 'u'
!-------------------------------------------------------------------!

interface inp_read_val
   module procedure     inp_read_val_i1, inp_read_val_in, &
                        inp_read_val_r1, inp_read_val_rn, &
                        inp_read_val_c1, inp_read_val_cn, &
                        inp_read_val_l1, inp_read_val_ln
end interface

private ::      inp_read_val_i1, inp_read_val_in, &
                inp_read_val_r1, inp_read_val_rn, &
                inp_read_val_c1, inp_read_val_cn, &
                inp_read_val_l1, inp_read_val_ln

contains

!-------------------------------------------------------------------!
!       read input function
!-------------------------------------------------------------------!

 function inp_init(mode) result(vinp)
   implicit none
   character(LEN=*), intent(in) :: mode
   type(Var_Input)              :: vinp

   logical         :: isfile

   if (inp_isinit) then
      if (pMast) call io_printerr('input initialize','input is already initialized')
      vinp%iserr = .true.
      return
   endif

   if (pMast) then
      select case(trim(adjustl(mode)))
      case ('OPT')
         isfile = io_assert_file('opt.in')
         inp_iserr = .not. isfile
         if (isfile) vinp = inp_read_opt('opt.in')
      case ('GEN')
         isfile = io_assert_file('gen.in')
         inp_iserr = .not. isfile
         if (isfile) vinp = inp_read_gen('gen.in')
      end select
      vinp%mode = trim(adjustl(mode))
      vinp%iserr = inp_iserr
   endif

   call inp_bcast(vinp)
   if (vinp%iserr) return

   if (pMast) then
      select case (trim(adjustl(mode)))
      case ('OPT')
         call inp_printinfo_opt(vinp)
      case ('GEN')
         call inp_printinfo_gen(vinp)
      end select
   endif

   inp_isinit = .true.

 end function inp_init

!-------------------------------------------------------------------!

 subroutine inp_fin(vinp)
   implicit none
   type(Var_Input), intent(inout) :: vinp

   if (.not. inp_isinit) then
      if (pMast) call io_printerr('input finalize','input is not initialized')
      vinp%iserr = .true.
      return
   endif

   if (allocated (vinp%nSpecs)) deallocate (vinp%nSpecs)
   if (allocated (vinp%nAtoms)) deallocate (vinp%nAtoms)
   if (allocated (vinp%isConn)) deallocate (vinp%isConn)
   if (allocated (vinp%fname_coo_in)) deallocate (vinp%fname_coo_in)
   if (allocated (vinp%unit_cen)) deallocate (vinp%unit_cen)
   if (allocated (vinp%unit_dir)) deallocate (vinp%unit_dir)

   inp_isinit = .false.
   vinp%pp    = .false.

 end subroutine inp_fin

!-------------------------------------------------------------------!

 function inp_read_opt(fname) result(vinp)
   implicit none
   character(LEN=*), intent(in) :: fname
   type(Var_Input)              :: vinp

   integer               :: u, nUnit, iUnit
   logical               :: isfolder
   character(LEN=VARLEN) :: LattTypeChr
   integer, dimension(3) :: ldim

   u = io_getunit()

   vinp%fname_in = fname
   vinp%mode     = 'OPT'
   call inp_defaults(vinp)

   open(u, file=fname, status = 'old', action = 'read')
   call inp_read_val(u, 'NUM.UNITS', nUnit)
   vinp%nUnit = nUnit

   allocate (vinp%nSpecs(nUnit))
   allocate (vinp%nAtoms(nUnit))
   allocate (vinp%isConn(0:nUnit))
   allocate (vinp%fname_coo_in(nUnit))
   allocate (vinp%unit_cen(3,nUnit))
   allocate (vinp%unit_dir(3,nUnit))

   call inp_read_val(u, 'OPT.SITE',    vinp%nOptSite, .true.)
   call inp_read_val(u, 'OPT.MIG',     vinp%nOptMig,  .true.)
   call inp_read_val(u, 'OPT.TRY',     vinp%nOptTry,  .true.)
   call inp_read_val(u, 'SEED.MASTER', vinp%seed,     .true.)

   if (vinp%seed .ne. -999) vinp%lseed = .true.

   call inp_read_val(u, 'LATT.DIM',    ldim, 3)

   vinp%nx = ldim(1)
   vinp%ny = ldim(2)
   vinp%nz = ldim(3)

   call inp_read_val(u, 'N.SPECIES',    vinp%nSpecs, nUnit)
   call inp_read_val(u, 'N.ATOM',       vinp%nAtoms, nUnit)
   call inp_read_val(u, 'L.TRAJECT',    vinp%isTrace, .true.)
   call inp_read_val(u, 'L.CONNECT',    vinp%isConn(1:nUnit), nUnit)
   call inp_read_val(u, 'L.SELEC.SITE', vinp%selective_site, .true.)
   call inp_read_val(u, 'L.SELEC.UNIT', vinp%selective_unit, .true.)

   if (vinp%selective_unit) vinp%selective_site = .true.
   
   vinp%isConn(0) = .false.
   do iUnit = 1, nUnit
      if (vinp%isConn(iUnit)) vinp%isConn(0) = .true.
   enddo

   call inp_read_val(u, 'ENE.METHOD', vinp%method)
   call inp_read_val(u, 'LATT.TYPE',  LattTypeChr, .true.)

   call inp_read_val(u, 'FILE.COO.OUT', vinp%fname_coo_out)
   call inp_read_val(u, 'FILE.AMBER',   vinp%fname_amber)

   call inp_read_val(u, 'RATE.OCC',    vinp%rate_occ,    .true.)
   call inp_read_val(u, 'TEMPERATURE', vinp%temperature, .true.)
   if (vinp%temperature .lt. 1d-15) vinp%temperature = 1d-15

   call inp_read_fnames(u, 'FILE.COO.IN',    vinp%fname_coo_in, nUnit)
   call inp_read_vecs  (u, 'UNIT.CENTER',    vinp%unit_cen,     nUnit)
   call inp_read_vecs  (u, 'UNIT.DIRECTION', vinp%unit_dir,     nUnit)

   close(u)

   if (vinp%isTrace) then
      inquire(directory = 'trjs', exist = isfolder)
      if (.not. isfolder) call system('mkdir trjs')
   endif

   LattTypeChr = io_getcapital(trim(adjustl(LattTypeChr)))

   select case (trim(adjustl(LattTypeChr)))
   case default 
      vinp%LattType = 0
      vinp%nSiteTotal = vinp%nx*vinp%ny*vinp%nz
   case ('S','SIMPLE')
      vinp%LattType = 0
      vinp%nSiteTotal = vinp%nx*vinp%ny*vinp%nz
   case ('BOC','BOCL')
      vinp%LattType = 1
      vinp%nSiteTotal = vinp%nx*vinp%ny*vinp%nz*2
   case ('FC', 'FCL')
      vinp%LattType = 2
      vinp%nSiteTotal = vinp%nx*vinp%ny*vinp%nz*4
   case ('HCP', 'CP')
      vinp%LattType = 3
      vinp%nSiteTotal = vinp%nx*vinp%ny*vinp%nz*4
   case ('BAC','BACL')
      vinp%LattType = 4
      vinp%nSiteTotal = vinp%nx*vinp%ny*vinp%nz*4
   end select

   vinp%nSite = floor(vinp%nSiteTotal*vinp%rate_occ)
   if (vinp%nSite .eq. 0) vinp%nSite = 1
   if (vinp%isConn(0) .and. (mod(vinp%nSite,2) .eq. 1)) vinp%nSite = vinp%nSite + 1
   if (vinp%isConn(0)) call inp_extend(vinp)

 endfunction inp_read_opt

!-------------------------------------------------------------------!

 function inp_read_gen(fname) result(vinp)
   implicit none
   character(LEN=*), intent(in) :: fname
   type(Var_Input)              :: vinp

   integer :: u, nUnit, iUnit
   logical :: isfolder
   
   u = io_getunit()

   vinp%fname_in = fname
   vinp%mode     = 'GEN'
   call inp_defaults(vinp)

   open(u, file=fname, status='old', action='read')
   call inp_read_val(u, 'NUM.UNITS', nUnit)
   if (nUnit .ne. 2) then
      call io_printerr('input parameter - NUM.UNITS','3 or more units are not supported yet')
      return
   endif
   vinp%nUnit = nUnit

   allocate (vinp%nSpecs(nUnit))
   allocate (vinp%nAtoms(nUnit))
   allocate (vinp%isConn(0:nUnit))
   allocate (vinp%fname_coo_in(nUnit))
   allocate (vinp%unit_cen(3,nUnit))
   allocate (vinp%unit_dir(3,nUnit))

   call inp_read_val(u, 'N.SPECIES',   vinp%nSpecs, nUnit)
   call inp_read_val(u, 'N.ATOM',      vinp%nAtoms, nUnit)

   call inp_read_val(u, 'FILE.COO.OUT', vinp%fname_coo_out)

   call inp_read_fnames(u, 'FILE.COO.IN',    vinp%fname_coo_in, nUnit)
   call inp_read_vecs  (u, 'UNIT.CENTER',    vinp%unit_cen,     nUnit)
   call inp_read_vecs  (u, 'UNIT.DIRECTION', vinp%unit_dir,     nUnit)

   call inp_read_val(u, 'ROT.GLOBAL.ANG',   vinp%AngRotG)
   call inp_read_val(u, 'ROT.GLOBAL.GRID',  vinp%nGridRotG)
   call inp_read_val(u, 'ROT.LOCAL.ANG.X',  vinp%AngRotL(1))
   call inp_read_val(u, 'ROT.LOCAL.ANG.Y',  vinp%AngRotL(2))
   call inp_read_val(u, 'ROT.LOCAL.ANG.Z',  vinp%AngRotL(3))
   call inp_read_val(u, 'ROT.LOCAL.GRID.X', vinp%nGridRotL(1))
   call inp_read_val(u, 'ROT.LOCAL.GRID.Y', vinp%nGridRotL(2))
   call inp_read_val(u, 'ROT.LOCAL.GRID.Z', vinp%nGridRotL(3))
   call inp_read_val(u, 'TRANS.DIS.X',  vinp%DisTrans(:,1), 2)
   call inp_read_val(u, 'TRANS.DIS.Y',  vinp%DisTrans(:,2), 2)
   call inp_read_val(u, 'TRANS.DIS.Z',  vinp%DisTrans(:,3), 2)
   call inp_read_val(u, 'TRANS.GRID.X', vinp%nGridTrans(1))
   call inp_read_val(u, 'TRANS.GRID.Y', vinp%nGridTrans(2))
   call inp_read_val(u, 'TRANS.GRID.Z', vinp%nGridTrans(3))

   vinp%AngRotG = vinp%AngRotG*pi/1.8d2
   vinp%AngRotL(:) = vinp%AngRotL(:)*pi/1.8d2

   close(u)

   inquire (directory='outputs', exist=isfolder)
   if (.not. isfolder) call system('mkdir outputs')

 end function inp_read_gen

!-------------------------------------------------------------------!

 subroutine inp_read_network(fname_in, nSpec, NameAtom, networks, iserr)
   implicit none
   character(LEN=*),                   intent(in)  :: fname_in
   integer,                            intent(in)  :: nSpec
   character(LEN=*), dimension(nSpec), intent(in)  :: NameAtom
   character(LEN=*), dimension(nSpec), intent(out) :: networks
   logical,                            intent(out) :: iserr

   integer                :: iSpec1, iSpec2, u, CHK, ii, idx
   character(LEN=LINELEN) :: line

   character(LEN=3),       dimension(nSpec) :: NameAtom_net
   character(LEN=PATHLEN), dimension(nSpec) :: f_names

   iserr       = .false.
   networks(:) = 'NONE'
   u = io_getunit()
   open(u, file=fname_in, status = 'old', action = 'read')

   call io_findline(u, 'AENET.NETWORKS', line, CHK)
   do iSpec1 = 1, nSpec
      read (u, '(A)') line
      write (line, '(A)') trim(adjustl(line))
      read (line, *) NameAtom_net(iSpec1), f_names(iSpec1)
   enddo
   close(u)

   loop_spec_1 : do iSpec1 = 1, nSpec
      loop_spec_2 : do iSpec2 = 1, nSpec
         if (trim(NameAtom_net(iSpec2)) .eq. trim(NameAtom(iSpec1))) then
            networks(iSpec1) = trim(adjustl(f_names(iSpec2)))
            exit loop_spec_2
         endif
      enddo loop_spec_2
      if (networks(iSpec1) .eq. 'NONE') then
         call io_printerr('aenet network not found with "AENET.NETWORKS"', NameAtom(iSpec1))
         iserr = .true.
      endif
   enddo loop_spec_1

end subroutine inp_read_network

!-------------------------------------------------------------------!
!           broadcast input variables
!-------------------------------------------------------------------!

 subroutine inp_bcast(vinp)
   implicit none
   type(Var_Input), intent(inout) :: vinp

   if (vinp%mode .eq. 'GEN') return
   call pp_bcast(vinp%iserr) 
   if (vinp%iserr) return
   if (vinp%pp) then
      if (pMast) call io_printerr('input bcast','input is already bcasted')
      vinp%iserr = .true.
      return
   endif

   call pp_bcast(vinp%mode)
   call pp_bcast(vinp%nUnit)
   call pp_bcast(vinp%seed)
   call pp_bcast(vinp%lseed)
   if (.not. allocated (vinp%nSpecs)) allocate (vinp%nSpecs(vinp%nUnit))
   call pp_bcast(vinp%nSpecs, vinp%nUnit)
   if (.not. allocated (vinp%nAtoms)) allocate (vinp%nAtoms(vinp%nUnit))
   call pp_bcast(vinp%nAtoms, vinp%nUnit)

   call pp_bcast(vinp%nOptSite)
   call pp_bcast(vinp%nOptMig)
   call pp_bcast(vinp%nOptTry)
   call pp_bcast(vinp%nx)
   call pp_bcast(vinp%ny)
   call pp_bcast(vinp%nz)
   call pp_bcast(vinp%nSiteTotal)
   call pp_bcast(vinp%nSite)
   call pp_bcast(vinp%ispbc)
   if (.not. allocated (vinp%isConn)) allocate (vinp%isConn(0:vinp%nUnit))
   call pp_bcast(vinp%isConn, vinp%nUnit+1)
   call pp_bcast(vinp%method)
   call pp_bcast(vinp%LattType)

   vinp%pp = .true.

 end subroutine inp_bcast

!-------------------------------------------------------------------!
!           default input variables
!-------------------------------------------------------------------!

 subroutine inp_defaults(vinp)
   implicit none
   type(Var_Input), intent(inout) :: vinp

   vinp%nUnit           = 1

   vinp%nSpecs(:)       = 0
   vinp%nAtoms(:)       = 0

   vinp%isConn(:)       = .false.

   vinp%fname_coo_out   = 'CONTCAR'
   vinp%fname_coo_in(:) = 'NONE'

   vinp%unit_cen(:,:)   = 0d0
   vinp%unit_dir(:,:)   = 0d0
   vinp%unit_dir(3,:)   = 1d0

   select case (vinp%mode)
   case ('OPT')
      vinp%nOptSite       = 1500
      vinp%nOptMig        = 2000
      vinp%nOptTry        = 200
      vinp%nx             = 1
      vinp%ny             = 1
      vinp%nz             = 1
      vinp%seed           = -999

      vinp%isTrace        = .true.
      vinp%selective_site = .false.
      vinp%selective_unit = .false.
      vinp%lseed          = .false.
      vinp%ispbc          = .true.
      vinp%lforce         = .false.

      vinp%method         = 'NONE'
      vinp%LattType       = 0

      vinp%fname_amber    = 'NONE'

      vinp%temperature    = 3d2
      vinp%rate_occ       = 7d-1
      vinp%rate_dis       = 7d-1
      vinp%rate_ang       = 3d-3
   case ('GEN')
      vinp%AngRotG        = 1.8d2
      vinp%nGridRotG      = 5
      vinp%AngRotL(:)     = 1.8d2
      vinp%nGridRotL(:)   = 5
      vinp%DisTrans(1,:)  = -3d0
      vinp%DisTrans(2,:)  = 3d0
      vinp%nGridTrans(:)  = 5
   end select

 end subroutine inp_defaults

!-------------------------------------------------------------------!
!           print input variables
!-------------------------------------------------------------------!

 subroutine inp_printinfo_opt(vinp)
   implicit none
   type (Var_Input), intent(in) :: vinp

   integer               :: iUnit
   character(LEN=VARLEN) :: LattType_long

   select case(vinp%LattType)
   case (0)
      LattType_long = 'Simple Lattice'
   case(1)
      LattType_long = 'Body Centered Lattice'
   case(2)
      LattType_long = 'Face Centered Lattice'
   case(3)
      LattType_long = 'Hexagonal Closed Packing'
   end select

   call io_printhead('Input parameters')
   write (*,201) ' Input file                      : ', 'opt.in'
   write (*,*) 
   write (*,200) ' Parameters for geometry optimization ------------------------------------'
   write (*,203) '   Max.Num. of site selection    : ', vinp%nOptSite
   write (*,203) '   Max.Num. of migration step    : ', vinp%nOptMig 
   write (*,203) '   Number of trial structure     : ', vinp%nOptTry
   write (*,205) '   Temperature                   : ', vinp%temperature
   write (*,*)
   write (*,200) ' Lattice information -----------------------------------------------------'
   write (*,201) '   Lattice type                  : ', trim(LattType_long)
   write (*,204) '   Lattice dimension (x,y,z)     : ', vinp%nx, vinp%ny, vinp%nz
   write (*,203) '   Total number of lattice sites : ', vinp%nSiteTotal
   write (*,203) '   Number of occupied sites      : ', vinp%nSite 
   write (*,205) '   Lattice occupation ratio      : ', float(vinp%nSite)/float(vinp%nSiteTotal)
   write (*,202) '   Connection between latt.point : ', vinp%isConn(0)
   write (*,*)     
   write (*,200) ' File management --------------------------------------------------------'
   write (*,203) '   Number of input coordination  : ', vinp%nUnit
  do iUnit = 1, vinp%nUnit
   write (*,201) '   Input coordination filename   : ', trim(adjustl(vinp%fname_coo_in(iUnit)))
   write (*,206) '     Center of coordinate        : ', vinp%unit_cen(:,iUnit)
   write (*,206) '     Direction of coordinate     : ', vinp%unit_dir(:,iUnit)
  enddo
   write (*,201) '   Output coordination filename  : ', trim(adjustl(vinp%fname_coo_out))
   write (*,*)

200 format(1x,A)
201 format(1x,A,A)
202 format(1x,A,L8)
203 format(1x,A,I8)
204 format(1x,A,3I4)
205 format(1x,A,F13.4)
206 format(1x,A,3F10.4)

 end subroutine inp_printinfo_opt

!-------------------------------------------------------------------!

 subroutine inp_printinfo_gen(vinp)
   implicit none
   type (Var_Input), intent(in) :: vinp

   integer               :: iUnit

   call io_printhead('Input parameters')
   write (*,201) ' Input file                      : ', 'gen.in'
   write (*,*) 
   write (*,200) ' Parameters for geometry generation --------------------------------------'
   write (*,202) '   Global rotation angle         : ', vinp%AngRotG
   write (*,203) '   Global rotation grid          : ', vinp%nGridRotG
   write (*,204) '   Local rotation angles (x,y,z) : ', vinp%AngRotL(1:3)
   write (*,205) '   Local rotation grid   (x,y,z) : ', vinp%nGridRotL(1:3)
   write (*,206) '   Translation min/max (y)       : ', vinp%DisTrans(:,1)
   write (*,207) '   Translation grid    (y)       : ', vinp%nGridTrans(1)
   write (*,206) '   Translation min/max (z)       : ', vinp%DisTrans(:,2)
   write (*,207) '   Translation grid    (z)       : ', vinp%nGridTrans(2)
   write (*,*)
   write (*,200) ' File management --------------------------------------------------------'
   write (*,203) '   Number of input coordination  : ', vinp%nUnit
  do iUnit = 1, vinp%nUnit
   write (*,201) '   Input coordination filename   : ', trim(adjustl(vinp%fname_coo_in(iUnit)))
   write (*,204) '     Center of coordinate        : ', vinp%unit_cen(:,iUnit)
   write (*,204) '     Direction of coordinate     : ', vinp%unit_dir(:,iUnit)
  enddo
   write (*,201) '   Output coordination filename  : ', trim(adjustl(vinp%fname_coo_out))
   write (*,*)

200 format(1x,A)
201 format(1x,A,A)
202 format(1x,A,F12.3)
203 format(1x,A,I8)
204 format(1x,A,3F10.3)
205 format(1x,A,I6,4x,I6,4x,I6)
206 format(1x,A,2F10.3)
207 format(1x,A,I6,4x,I6)

 end subroutine inp_printinfo_gen

!-------------------------------------------------------------------!
!           duplicate non-connected molecules in connected lattice
!-------------------------------------------------------------------!

 subroutine inp_extend(vinp)
   implicit none

   type(Var_Input), intent(inout) :: vinp

   integer,                dimension(:),   allocatable :: nSpecs_ext(:)
   integer,                dimension(:),   allocatable :: nAtoms_ext(:)
   logical,                dimension(:),   allocatable :: isConn_ext(:)
   character(LEN=PATHLEN), dimension(:),   allocatable :: fname_in_ext(:)
   real(kind=dp),          dimension(:,:), allocatable :: unit_cen_ext(:,:)
   real(kind=dp),          dimension(:,:), allocatable :: unit_dir_ext(:,:)
   
   integer :: iUnit, iCnt, nUnit

   nUnit = vinp%nUnit
   do iUnit = 1, vinp%nUnit
      if (.not. vinp%isConn(iUnit)) nUnit = nUnit + 1
   enddo

   allocate (nSpecs_ext(nUnit))
   allocate (nAtoms_ext(nUnit))
   allocate (isConn_ext(0:nUnit))
   allocate (fname_in_ext(nUnit))
   allocate (unit_cen_ext(3,nUnit))
   allocate (unit_dir_ext(3,nUnit))

   isConn_ext(0) = vinp%isConn(0)
   iCnt = 0
   do iUnit = 1, vinp%nUnit
      iCnt = iCnt + 1
      nSpecs_ext(iCnt)     = vinp%nSpecs(iUnit)
      nAtoms_ext(iCnt)     = vinp%nAtoms(iUnit)
      isConn_ext(iCnt)     = vinp%isConn(iUnit)
      fname_in_ext(iCnt)   = vinp%fname_coo_in(iUnit)
      unit_cen_ext(:,iCnt) = vinp%unit_cen(:,iUnit)
      unit_dir_ext(:,iCnt) = vinp%unit_dir(:,iUnit)
      if (.not. vinp%isConn(iUnit)) then
         iCnt = iCnt + 1
         nSpecs_ext(iCnt)     = vinp%nSpecs(iUnit)
         nAtoms_ext(iCnt)     = vinp%nAtoms(iUnit)
         isConn_ext(iCnt)     = vinp%isConn(iUnit)
         fname_in_ext(iCnt)   = vinp%fname_coo_in(iUnit)
         unit_cen_ext(:,iCnt) = vinp%unit_cen(:,iUnit)
         unit_dir_ext(:,iCnt) = -vinp%unit_dir(:,iUnit)
      endif
   enddo

   vinp%nUnit = nUnit

   deallocate (vinp%nSpecs)
   deallocate (vinp%nAtoms)
   deallocate (vinp%isConn)
   deallocate (vinp%fname_coo_in)
   deallocate (vinp%unit_cen)
   deallocate (vinp%unit_dir)

   allocate (vinp%nSpecs(nUnit))
   allocate (vinp%nAtoms(nUnit))
   allocate (vinp%isConn(0:nUnit))
   allocate (vinp%fname_coo_in(nUnit))
   allocate (vinp%unit_cen(3,nUnit))
   allocate (vinp%unit_dir(3,nUnit))

   vinp%nSpecs       = nSpecs_ext 
   vinp%nAtoms       = nAtoms_ext
   vinp%isConn       = isConn_ext
   vinp%fname_coo_in = fname_in_ext
   vinp%unit_cen     = unit_cen_ext
   vinp%unit_dir     = unit_dir_ext
   
   deallocate (nSpecs_ext)
   deallocate (nAtoms_ext)
   deallocate (isConn_ext)
   deallocate (fname_in_ext)
   deallocate (unit_cen_ext)
   deallocate (unit_dir_ext)
   
 end subroutine inp_extend

!===================================================================!
!           read input variables
!===================================================================!
!-------------------------------------------------------------------!
!           read filenames
!-------------------------------------------------------------------!

 subroutine inp_read_fnames(u, chr_in, fnames, nUnit)
   implicit none
   integer,                            intent(in)  :: u
   integer,                            intent(in)  :: nUnit
   character(LEN=*),                   intent(in)  :: chr_in
   character(LEN=*), dimension(nUnit), intent(out) :: fnames

   character(LEN=LINELEN) :: line
   integer                :: iUnit, IOCHK, CHK
   logical                :: isFile

   isfile = .true.

   call io_findline(u, chr_in, line, CHK)
   if (CHK .ne. CHKOK) then
      call io_printerr('input error - '//chr_in,'KEYWORD is missing - '//chr_in)
      goto 100
   endif

   do iUnit = 1, nUnit
        ! loop for get filenames
      read(u,'(A)', iostat=IOCHK, end=100) line
      call io_checkdim(line, 1, CHK)
      if (CHK .ne. CHKOK) then
         write (*,*) CHK
         call io_printerr('input error - '//chr_in,'line is empty',line)
         goto 100
      endif
      call io_linesplit(line, fnames(iUnit), CHK)
      if (CHK .ne. CHKOK) goto 100
      if (trim(fnames(iUnit)) .eq. 'END') then
         call io_printerr('input error - '//chr_in,'needs more file')
         goto 100
      endif
      isFile = isFile .and. io_assert_file(fnames(iUnit))
   enddo
   if (.not. isfile) goto 100
   return

100 continue
   if ((CHK .le. CHKEOF) .or. (IOCHK .le. CHKEOF)) &
        call io_printerr('input error - '//chr_in,'touched EoF')

   inp_iserr = .true.

 end subroutine inp_read_fnames

!-------------------------------------------------------------------!
!           read vectors
!-------------------------------------------------------------------!

 subroutine inp_read_vecs(u, chr_in, vecs, nUnit)
   implicit none
   integer,                           intent(in)  :: u
   integer,                           intent(in)  :: nUnit
   character(LEN=*),                  intent(in)  :: chr_in
   real(kind=dp), dimension(3,nUnit), intent(out) :: vecs

   character(LEN=LINELEN) :: line
   integer                :: iUnit, CHK, IOCHK 

   call io_findline(u, chr_in, line, CHK)
   if (CHK .ne. CHKOK) then
      call io_printerr('input error - '//chr_in,'KEYWORD is missing - '//chr_in)
      goto 100
   endif

   do iUnit = 1, nUnit
      read(u, '(A)', iostat=IOCHK, end=100) line
      call io_checkdim(line, 3, CHK)
      if (CHK .ne. CHKOK) then
         call io_printerr('input error - '//chr_in,'need 3 values')
         goto 100
      endif
      call io_linesplit(line, 3, vecs(:,iUnit), CHK)
      if (CHK .ne. CHKOK) goto 100
   enddo
   return

100 continue
   if ((CHK .le. CHKEOF) .or. (IOCHK .le. CHKEOF)) write (*,*) 'touched end of file'
   inp_iserr = .true.

 end subroutine inp_read_vecs

!-------------------------------------------------------------------!
!           integer
!-------------------------------------------------------------------!

 subroutine inp_read_val_i1(u, chr_in, val_out, skip_err_in) 
   implicit none
   integer,           intent(in)    :: u
   character(LEN=*),  intent(in)    :: chr_in
   integer,           intent(inout) :: val_out
   logical, optional, intent(in)    :: skip_err_in

   character(LEN=LINELEN) :: line
   integer :: CHK
   logical :: skip_err = .false.

   if (present(skip_err_in)) skip_err = skip_err_in

   call io_findline(u, chr_in, line, CHK)
   if (CHK .eq. CHKOK) then
      call io_checkdim(line, 1, CHK)
      if (CHK .ne. CHKOK) goto 100
      call io_linesplit(line, val_out, CHK)
      if (CHK .ne. CHKOK) goto 100
   elseif (CHK .eq. CHKERR) then
      goto 100  
   endif
   return       

100 continue
   if (skip_err) return
        ! 1. error occured due to 'incorrectly assigned varible'
        ! 2. error during 'io_findline'
   inp_iserr = .true.
   call io_printerr('input reading - '//chr_in,"KEYWORD not found")

 end subroutine inp_read_val_i1

!-------------------------------------------------------------------!

 subroutine inp_read_val_in(u, chr_in, val_out, n, skip_err_in)
   implicit none
   integer,                        intent(in)    :: u
   integer,                        intent(in)    :: n
   character(LEN=*),               intent(in)    :: chr_in
   integer,          dimension(n), intent(inout) :: val_out
   logical,          optional,     intent(in)    :: skip_err_in

   character(LEN=LINELEN) :: line
   integer :: CHK
   logical :: skip_err = .false.

   if (present(skip_err_in)) skip_err = skip_err_in

   call io_findline(u, chr_in, line, CHK)
   if (CHK .eq. CHKOK) then
      call io_checkdim(line, n, CHK)
      if (CHK .ne. CHKOK) goto 100
      call io_linesplit(line, n, val_out, CHK)
      if (CHK .ne. CHKOK) goto 100
   elseif (CHK .eq. CHKERR) then
      goto 100
   endif      
   return

100 continue
   if (skip_err) return
        ! 1. error occured due to 'not properly assigned varible'
        ! 2. error during 'io_findline'
   inp_iserr = .true.
   call io_printerr('input reading - '//chr_in,"KEYWORD not found")

 end subroutine inp_read_val_in

!-------------------------------------------------------------------!
!           real
!-------------------------------------------------------------------!

 subroutine inp_read_val_r1(u, chr_in, val_out, skip_err_in)
   implicit none
   integer,           intent(in)    :: u
   character(LEN=*),  intent(in)    :: chr_in
   real(kind=dp),     intent(inout) :: val_out
   logical, optional, intent(in)    :: skip_err_in

   character(LEN=LINELEN) :: line
   integer                :: CHK
   logical :: skip_err = .false.

   if (present(skip_err_in)) skip_err = skip_err_in


   call io_findline(u, chr_in, line, CHK)
   if (CHK .eq. CHKOK) then
      call io_checkdim(line, 1, CHK)
      if (CHK .ne. CHKOK) goto 100
      call io_linesplit(line, val_out, CHK)
      if (CHK .ne. CHKOK) goto 100
   elseif (CHK .eq. CHKERR) then
      goto 100
   endif      
   return

100 continue
   if (skip_err) return
        ! 1. error occured due to 'not properly assigned varible'
   inp_iserr = .true.
   call io_printerr('input reading - '//chr_in,"KEYWORD not found")

 end subroutine inp_read_val_r1

!-------------------------------------------------------------------!

 subroutine inp_read_val_rn(u, chr_in, val_out, n, skip_err_in)
   implicit none
   integer,                     intent(in)    :: u
   integer,                     intent(in)    :: n
   character(LEN=*),            intent(in)    :: chr_in
   real(kind=dp), dimension(n), intent(inout) :: val_out
   logical,       optional,     intent(in)    :: skip_err_in

   character(LEN=LINELEN) :: line
   integer                :: CHK
   logical :: skip_err = .false.

   if (present(skip_err_in)) skip_err = skip_err_in

   call io_findline(u, chr_in, line, CHK)
   if (CHK .eq. CHKOK) then
      call io_checkdim(line, n, CHK)
      if (CHK .ne. CHKOK) goto 100
      call io_linesplit(line, n, val_out, CHK)
      if (CHK .ne. CHKOK) goto 100
   elseif (CHK .eq. CHKERR) then
      goto 100
   endif      
   return

100 continue
   if (skip_err) return
        ! 1. error occured due to 'not properly assigned varible'
        ! 2. error during 'io_findline'
   inp_iserr = .true.
   call io_printerr('input reading - '//chr_in,"KEYWORD not found")

 end subroutine inp_read_val_rn

!-------------------------------------------------------------------!
!           character
!-------------------------------------------------------------------!

 subroutine inp_read_val_c1(u, chr_in, val_out, skip_err_in)
   implicit none
   integer,           intent(in)    :: u
   character(LEN=*),  intent(in)    :: chr_in
   character(LEN=*),  intent(inout) :: val_out
   logical, optional, intent(in)    :: skip_err_in

   character(LEN=LINELEN) :: line
   integer                :: CHK
   logical :: skip_err = .false.

   if (present(skip_err_in)) skip_err = skip_err_in

   call io_findline(u, chr_in, line, CHK)
   if (CHK .eq. CHKOK) then
      call io_checkdim(line, 1, CHK)
      if (CHK .ne. CHKOK) goto 100
      call io_linesplit(line, val_out, CHK)
      if (CHK .ne. CHKOK) goto 100
   elseif (CHK .eq. CHKERR) then
      goto 100
   endif      
   return

100 continue
   if (skip_err) return
        ! 1. error occured due to 'not properly assigned varible'
        ! 2. error during 'io_findline'
   inp_iserr = .true.
   call io_printerr('input reading - '//chr_in,"KEYWORD not found")

 end subroutine inp_read_val_c1

!-------------------------------------------------------------------!

 subroutine inp_read_val_cn(u, chr_in, val_out, n, skip_err_in)
   implicit none
   integer,                        intent(in)    :: u
   integer,                        intent(in)    :: n
   character(LEN=*),               intent(in)    :: chr_in
   character(LEN=*), dimension(n), intent(inout) :: val_out
   logical,          optional,     intent(in)    :: skip_err_in

   character(LEN=LINELEN) :: line
   integer                :: CHK
   logical :: skip_err = .false.

   if (present(skip_err_in)) skip_err = skip_err_in

   call io_findline(u, chr_in, line, CHK)
   if (CHK .eq. CHKOK) then
      call io_checkdim(line, n, CHK)
      if (CHK .ne. CHKOK) goto 100
      call io_linesplit(line, n, val_out, CHK)
      if (CHK .ne. CHKOK) goto 100
   elseif (CHK .eq. CHKERR) then
      goto 100
   endif      
   return

100 continue
   if (skip_err) return
        ! 1. error occured due to 'not properly assigned varible'
        ! 2. error during 'io_findline'
   inp_iserr = .true.
   call io_printerr('input reading - '//chr_in,"KEYWORD not found")

 end subroutine inp_read_val_cn

!-------------------------------------------------------------------!
!           logical
!-------------------------------------------------------------------!

 subroutine inp_read_val_l1(u, chr_in, val_out, skip_err_in)
   implicit none
   integer,           intent(in)    :: u
   character(LEN=*),  intent(in)    :: chr_in
   logical,           intent(inout) :: val_out
   logical, optional, intent(in)    :: skip_err_in

   character(LEN=LINELEN) :: line
   integer                :: CHK
   logical :: skip_err = .false.

   if (present(skip_err_in)) skip_err = skip_err_in

   call io_findline(u, chr_in, line, CHK)
   if (CHK .eq. CHKOK) then
      call io_checkdim(line, 1, CHK)
      if (CHK .ne. CHKOK) goto 100
      call io_linesplit(line, val_out, CHK)
      if (CHK .ne. CHKOK) goto 100
   elseif (CHK .eq. CHKERR) then
      goto 100
   endif      
   return       ! if failed to find 'keyword (chr_in)' in file, return defaults

100 continue
   if (skip_err) return
        ! 1. error occured due to 'not properly assigned varible'
        ! 2. error during 'io_findline'
   inp_iserr = .true.
   call io_printerr('input reading - '//chr_in,"KEYWORD not found")

 end subroutine inp_read_val_l1

!-------------------------------------------------------------------!

 subroutine inp_read_val_ln(u, chr_in, val_out, n, skip_err_in)
   implicit none
   integer,                        intent(in)    :: u
   integer,                        intent(in)    :: n
   character(LEN=*),               intent(in)    :: chr_in
   logical,          dimension(n), intent(inout) :: val_out
   logical,          optional,     intent(in)    :: skip_err_in

   character(LEN=LINELEN) :: line
   integer                :: CHK
   logical :: skip_err = .false.

   if (present(skip_err_in)) skip_err = skip_err_in

   call io_findline(u, chr_in, line, CHK)
   if (CHK .eq. CHKOK) then
      call io_checkdim(line, n, CHK)
      if (CHK .ne. CHKOK) goto 100
      call io_linesplit(line, n, val_out, CHK)
      if (CHK .ne. CHKOK) goto 100
   elseif (CHK .eq. CHKERR) then
      goto 100
   endif      
   return       ! if failed to find 'keyword (chr_in)' in file, return defaults

100 continue
   if (skip_err) return
        ! 1. error occured due to 'not properly assigned varible'
        ! 2. error during 'io_findline'
   inp_iserr = .true.
   call io_printerr('input reading - '//chr_in,"KEYWORD not found")

 end subroutine inp_read_val_ln

!-------------------------------------------------------------------!

end module input
