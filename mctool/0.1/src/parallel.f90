module parallel

  use constants,        only :  dp
!  use mt_stream,        only :  mt_state
  use util,             only :  gen_rseed
  use iotool,         only :  io_printhead

implicit none
save

include 'mpif.h'

!-------------------------------------------------------------------!
!       mpi variables
!-------------------------------------------------------------------!
!       mts     : mt_stream input
!       pMast   : if nrnk is root, return .true.
!       nrnk    : rank of processor
!       nprc    : number of processor
!       status  : mpi status
!       
!       pp_vinp : if Var_Input is broadcasted, return .true.
!       pp_vcoo : if Var_Coord is broadcasted, return .true.
!       pp_vlat : if Var_Latt is broadcasted, return .true.
!
!       ppInit  : if mpi is initialized, return .true.
!       IERR    : MPI error code
!       root    : rank of root processor
!-------------------------------------------------------------------!

!type(mt_state),                      public :: mts
integer,                             public :: root   = 0
integer,                             public :: nrnk
integer,                             public :: nprc
integer, dimension(MPI_STATUS_SIZE), public :: status

logical,                             public :: pMast

integer,                             private :: IERR
logical,                             private :: ppInit = .false.

!-------------------------------------------------------------------!
!       subroutines for parallel run
!-------------------------------------------------------------------!
!       pp_init()       : initialize parallel run
!       pp_fin()        : finilize parallel run
!       pp_info()       : print parallel run information
!       pp_mtserial(n)  : initialize mt_19937 random number generator
!                         for root processor with seed number 'n'
!       pp_mtstream(n)  : initialize mt_stream parallel random number
!                         generator with seed number 'n'
!       pp_barrier()    : mpi_barrier for MPI_COMM_WORLD
!-------------------------------------------------------------------!

private ::      pp_info
public ::       pp_init,        &
                pp_fin,         &
                pp_mtserial,    &
!                pp_mtstream,    &
                pp_barrier

!-------------------------------------------------------------------!
!       pp_bcast     : broadcast variable with CALL MPI_BCAST 
!-------------------------------------------------------------------!
!       pp_bcast_i   : one integer
!       pp_bcast_i1d : integer array 
!       pp_bcast_i2d : 2-dimensional integer matrix
!       pp_bcast_i3d : 3-dimensional integer matrix
!       pp_bcast_l   : on logical variable
!       pp_bcast_l1d : logical array
!       pp_bcast_d   : one real variable
!       pp_bcast_d1d : double precision array
!       pp_bcast_d2d : 2-dim double precision matrix
!       pp_bcast_d3d : 3-dim double precision matrix
!       pp_bcast_d4d : 3-dim double precision matrix
!       pp_bcast_c   : one string
!       pp_bcast_c1d : string array
!-------------------------------------------------------------------!

interface pp_bcast
   module procedure     pp_bcast_i,     pp_bcast_i1d,   &
                        pp_bcast_i2d,   pp_bcast_i3d,   &
                        pp_bcast_l,     pp_bcast_l1d,   &
                        pp_bcast_d,     pp_bcast_d1d,   &
                        pp_bcast_d2d,   pp_bcast_d3d,   &
                        pp_bcast_d4d,                   &
                        pp_bcast_c,     pp_bcast_c1d
end interface pp_bcast

private ::      pp_bcast_i,     pp_bcast_i1d,   &
                pp_bcast_i2d,   pp_bcast_i3d,   &
                pp_bcast_l,     pp_bcast_l1d,   &
                pp_bcast_d,     pp_bcast_d1d,   &
                pp_bcast_d2d,   pp_bcast_d3d,   &
                pp_bcast_d4d,                   &
                pp_bcast_c,     pp_bcast_c1d

public ::       pp_bcast

!-------------------------------------------------------------------!
!       pp_send     : send variable to 'rank' via CALL MPI_SEND 
!-------------------------------------------------------------------!
!       pp_send_d(val, rank) : one real variable 'val'
!       pp_send_d1d(val, rank, a) : 1-dim real matrix 'val'
!       pp_send_d4d(val, rank, a, b, c, d) : 4-dim real matrix 'val' 
!-------------------------------------------------------------------!

interface pp_send
   module procedure     pp_send_d,      pp_send_d1d,   &
                        pp_send_d4d
end interface pp_send

private ::      pp_send_d,      pp_send_d1d,    &
                pp_send_d4d

public ::       pp_send

!-------------------------------------------------------------------!
!       pp_recv     : receive variable from 'rank' via CALL MPI_RECV 
!-------------------------------------------------------------------!
!       pp_recv_d(val, rank) : one real variable
!       pp_recv_d1d(val, rank, a) : 1-dim real matrix 'val'
!       pp_recv_d4d(val, rank, a, b, c, d) : 4-dim real matrix 'val'
!-------------------------------------------------------------------!

interface pp_recv
   module procedure     pp_recv_d,      pp_recv_d1d,    &
                        pp_recv_d4d
end interface pp_recv

private ::      pp_recv_d,      pp_recv_d1d,    &
                pp_recv_d4d

public ::       pp_recv

contains

!===================================================================!
!       start parallel run
!===================================================================!

 subroutine pp_init()

   implicit none

   if (ppInit) return

   call mpi_init(IERR)

   call mpi_comm_size(MPI_COMM_WORLD, nprc, IERR)
   call mpi_comm_rank(MPI_COMM_WORLD, nrnk, IERR)

   if (nrnk .eq. root) then
      pMast = .true.
   else
      pMast = .false.
   endif

   ppInit = .true.

   call pp_info()

 end subroutine pp_init

!===================================================================!
!       finish parallel run
!===================================================================!

 subroutine pp_fin()
   implicit none

   if (.not. ppInit) return
   call mpi_finalize(IERR)

   ppInit = .false.
 
 end subroutine pp_fin

!===================================================================!
!       print parallel machine information
!===================================================================!

 subroutine pp_info()
   implicit none

   if (.not. ppInit) return
   if (.not. pMast) return

   if (nprc .eq. 1) call io_printhead('Serial run')
   if (nprc .ne. 1) call io_printhead('Parallel run')

   write (*,'(1X,A,I5)') ' Number of processors            : ',nprc 
   write (*,*)
 end subroutine pp_info

!===================================================================!
!       initialize Mersenne Twister
!===================================================================!

 subroutine pp_mtserial(seed_in)
   use mt_serial,       only :  sgrnd
   implicit none
   integer, optional, intent(in) :: seed_in

   integer                       :: seed, ii

   if (.not. ppInit) return
   if (.not. pMast) return

   if (.not. present(seed_in)) seed = gen_rseed()
   if (      present(seed_in)) seed = seed_in

   if (.not. present(seed_in)) then
      write (*,200) ' Mersenne Twister seed number --------------------------------------------'
      write (*,201) '   Master seed number            : ', seed
      write (*,200) '        note - same number of core and master seed number with input   '
      write (*,200) '               file were required for reproducing previous calculation.'
      write (*,*)
   else
      write (*,200) ' Input SEED NUMBER detected ----------------------------------------------'
      write (*,201) '   Master seed number from input : ', seed
      write (*,200) '        note - reproducing previous calculatioin. please check your inputs'
      write (*,*)
   endif
   call sgrnd(seed)

200 format(1x,A)
201 format(1x,A,I15)
202 format(1x,A,I3,A,I15)
 end subroutine pp_mtserial

!-------------------------------------------------------------------!

!  subroutine pp_mtstream(seed_in)
!    use mt_stream,       only :  set_mt19937,    &
!                                 new,            &
!                                 init,           &
!                                 create_stream
!    implicit none
!    integer, optional, intent(in) :: seed_in
! 
!    type(mt_state)                :: mts_seed
!    integer                       :: seed, ii
!    integer, dimension(nprc)      :: seeds
! 
!    if (.not. ppInit) return
! 
!    if (pMast) then
!       if (.not. present(seed_in)) seed = gen_rseed()
!       if (      present(seed_in)) seed = seed_in
! 
!       if (.not. present(seed_in)) then
!          write (*,200) ' Mersenne Twister seed number --------------------------------------------'
!          write (*,201) '   Master seed number            : ', seed
!          write (*,200) '        note - same number of core and master seed number with input   '
!          write (*,200) '               file were required for reproducing previous calculation.'
!          write (*,*)
!       else
!          write (*,200) ' Input SEED NUMBER detected ----------------------------------------------'
!          write (*,201) '   Master seed number from input : ', seed
!          write (*,200) '        note - reproducing previous calculatioin. please check your inputs'
!          write (*,*)
!       endif
!    endif
!    call pp_bcast(seed) 
! 
! ! set mt_stream
! 
!    call set_mt19937
!    call new(mts_seed)
!    call init(mts_seed, seed)
!    call create_stream(mts_seed, mts, nrnk)
! 
! 200 format(1x,A)
! 201 format(1x,A,I15)
! 202 format(1x,A,I3,A,I15)
!  end subroutine pp_mtstream

!===================================================================!
!       mpi_bcast of basic variables
!===================================================================!
!-------------------------------------------------------------------!
!       integer
!-------------------------------------------------------------------!

 subroutine pp_bcast_i(val, rank)
   implicit none
   integer,           intent(inout) :: val
   integer, optional, intent(in)    :: rank

   if (present(rank)) then
      call mpi_bcast(val, 1, MPI_INTEGER, rank, MPI_COMM_WORLD, IERR)
   else
      call mpi_bcast(val, 1, MPI_INTEGER, root, MPI_COMM_WORLD, IERR)
   endif

 end subroutine pp_bcast_i

!-------------------------------------------------------------------!

 subroutine pp_bcast_i1d(val, n, rank)
   implicit none
   integer,               intent(in)    :: n
   integer, dimension(n), intent(inout) :: val
   integer, optional,     intent(in)    :: rank

   if (present(rank)) then
      call mpi_bcast(val, n, MPI_INTEGER, rank, MPI_COMM_WORLD, IERR)
   else
      call mpi_bcast(val, n, MPI_INTEGER, root, MPI_COMM_WORLD, IERR)
   endif

 end subroutine pp_bcast_i1d

!-------------------------------------------------------------------!

 subroutine pp_bcast_i2d(val, n, m, rank)
   implicit none
   integer,                 intent(in)    :: n, m
   integer, dimension(n,m), intent(inout) :: val
   integer, optional,       intent(in)    :: rank

   if (present(rank)) then
      call mpi_bcast(val, n*m, MPI_INTEGER, rank, MPI_COMM_WORLD, IERR)
   else
      call mpi_bcast(val, n*m, MPI_INTEGER, root, MPI_COMM_WORLD, IERR)
   endif

 end subroutine pp_bcast_i2d

!-------------------------------------------------------------------!

 subroutine pp_bcast_i3d(val, n, m, l, rank)
   implicit none
   integer,                   intent(in)    :: n, m, l
   integer, dimension(n,m,l), intent(inout) :: val
   integer, optional,         intent(in)    :: rank

   if (present(rank)) then
      call mpi_bcast(val, n*m*l, MPI_INTEGER, rank, MPI_COMM_WORLD, IERR)
   else
      call mpi_bcast(val, n*m*l, MPI_INTEGER, root, MPI_COMM_WORLD, IERR)
   endif

 end subroutine pp_bcast_i3d

!-------------------------------------------------------------------!
!       logical
!-------------------------------------------------------------------!

 subroutine pp_bcast_l(val, rank)
   implicit none
   logical,           intent(inout) :: val
   integer, optional, intent(in)    :: rank

   if (present(rank)) then
      call mpi_bcast(val, 1, MPI_LOGICAL, rank, MPI_COMM_WORLD, IERR)
   else
      call mpi_bcast(val, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, IERR)
   endif

 end subroutine pp_bcast_l

!-------------------------------------------------------------------!

 subroutine pp_bcast_l1d(val, n, rank)
   implicit none
   integer,               intent(in)    :: n
   logical, dimension(n), intent(inout) :: val
   integer, optional,     intent(in)    :: rank

   if (present(rank)) then
      call mpi_bcast(val, n, MPI_LOGICAL, rank, MPI_COMM_WORLD, IERR)
   else
      call mpi_bcast(val, n, MPI_LOGICAL, root, MPI_COMM_WORLD, IERR)
   endif

 end subroutine pp_bcast_l1d

!-------------------------------------------------------------------!
!       real(kind=dp)
!-------------------------------------------------------------------!

 subroutine pp_bcast_d(val, rank)
   implicit none
   real(kind=dp),           intent(inout) :: val
   integer,       optional, intent(in)    :: rank

   if (present(rank)) then
      call mpi_bcast(val, 1, MPI_DOUBLE_PRECISION, rank, MPI_COMM_WORLD, IERR)
   else
      call mpi_bcast(val, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, IERR)
   endif

 end subroutine pp_bcast_d

!-------------------------------------------------------------------!

 subroutine pp_bcast_d1d(val, n, rank)
   implicit none
   integer,                     intent(in)    :: n
   real(kind=dp), dimension(n), intent(inout) :: val
   integer,       optional,     intent(in)    :: rank

   if (present(rank)) then
      call mpi_bcast(val, n, MPI_DOUBLE_PRECISION, rank, MPI_COMM_WORLD, IERR)
   else
      call mpi_bcast(val, n, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, IERR)
   endif

 end subroutine pp_bcast_d1d

!-------------------------------------------------------------------!

 subroutine pp_bcast_d2d(val, a, b, rank)
   implicit none
   integer,                       intent(in)    :: a, b
   real(kind=dp), dimension(a,b), intent(inout) :: val
   integer,       optional,       intent(in)    :: rank

   if (present(rank)) then
      call mpi_bcast(val, a*b, MPI_DOUBLE_PRECISION, rank, MPI_COMM_WORLD, IERR)
   else
      call mpi_bcast(val, a*b, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, IERR)
   endif

 end subroutine pp_bcast_d2d

!-------------------------------------------------------------------!

 subroutine pp_bcast_d3d(val, a, b, c, rank)
   implicit none
   integer,                         intent(in)    :: a, b, c
   real(kind=dp), dimension(a,b,c), intent(inout) :: val
   integer,       optional,         intent(in)    :: rank

   integer :: ic

   do ic = 1, c
   if (present(rank)) then
      call mpi_bcast(val(:,:,ic), a*b, MPI_DOUBLE_PRECISION, rank, MPI_COMM_WORLD, IERR)
   else
      call mpi_bcast(val(:,:,ic), a*b, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, IERR)
   endif
   enddo

 end subroutine pp_bcast_d3d

!-------------------------------------------------------------------!

 subroutine pp_bcast_d4d(val, a, b, c, d, rank)
   implicit none
   integer,                           intent(in)    :: a, b, c, d
   real(kind=dp), dimension(a,b,c,d), intent(inout) :: val
   integer,       optional,           intent(in)    :: rank

   integer :: ic, id

   do id = 1, d
   do ic = 1, c
   if (present(rank)) then
      call mpi_bcast(val(:,:,ic,id), a*b, MPI_DOUBLE_PRECISION, rank, MPI_COMM_WORLD, IERR)
   else
      call mpi_bcast(val(:,:,ic,id), a*b, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, IERR)
   endif
   enddo
   enddo

 end subroutine pp_bcast_d4d

!-------------------------------------------------------------------!
!       character
!-------------------------------------------------------------------!

 subroutine pp_bcast_c(val, rank)
   implicit none
   character(LEN=*),           intent(inout) :: val
   integer,          optional, intent(in)    :: rank

   integer :: k

   k = len(val)

   if (present(rank)) then
      call mpi_bcast(val, k, MPI_CHARACTER, rank, MPI_COMM_WORLD, IERR)
   else
      call mpi_bcast(val, k, MPI_CHARACTER, root, MPI_COMM_WORLD, IERR)
   endif

 end subroutine pp_bcast_c

!-------------------------------------------------------------------!

 subroutine pp_bcast_c1d(val, n, rank)
   implicit none
   integer,                        intent(in)    :: n
   character(LEN=*), dimension(n), intent(inout) :: val
   integer,          optional,     intent(in)    :: rank

   integer :: k

   k = len(val)

   if (present(rank)) then
      call mpi_bcast(val, k*n, MPI_CHARACTER, rank, MPI_COMM_WORLD, IERR)
   else
      call mpi_bcast(val, k*n, MPI_CHARACTER, root, MPI_COMM_WORLD, IERR)
   endif

 end subroutine pp_bcast_c1d

!===================================================================!
!       mpi_send / mpi_recv
!===================================================================!
!-------------------------------------------------------------------!
!       send : real(kind=dp)
!-------------------------------------------------------------------!

 subroutine pp_send_d(val, rank)
   implicit none
   real(kind=dp), intent(in) :: val
   integer,       intent(in) :: rank

   integer, parameter :: tag = 1

   call mpi_send(val, 1, MPI_DOUBLE_PRECISION, rank, tag, MPI_COMM_WORLD, IERR)

 end subroutine pp_send_d

!-------------------------------------------------------------------!

 subroutine pp_send_d1d(val, rank, a)
   implicit none
   integer,                     intent(in) :: a
   real(kind=dp), dimension(a), intent(in) :: val
   integer,                     intent(in) :: rank

   integer, parameter :: tag = 1

   call mpi_send(val, a, MPI_DOUBLE_PRECISION, rank, tag, MPI_COMM_WORLD, IERR)

 end subroutine pp_send_d1d

!-------------------------------------------------------------------!

 subroutine pp_send_d4d(val, rank, a, b, c, d)
   implicit none
   integer,                           intent(in) :: a, b, c, d
   real(kind=dp), dimension(a,b,c,d), intent(in) :: val
   integer,                           intent(in) :: rank

   integer, parameter :: tag = 1

   integer :: ii, ij

   do ii = 1, d
   do ij = 1, c
      call mpi_send(val(:,:,ij,ii), a*b, MPI_DOUBLE_PRECISION, rank, tag, MPI_COMM_WORLD, IERR)
   enddo
   enddo


 end subroutine pp_send_d4d

!-------------------------------------------------------------------!
!       recv : real(kind=dp)
!-------------------------------------------------------------------!

 subroutine pp_recv_d(val, rank)
   implicit none
   real(kind=dp), intent(out) :: val
   integer,       intent(in)  :: rank

   integer, parameter :: tag = 1

   call mpi_recv(val, 1, MPI_DOUBLE_PRECISION, rank, tag, MPI_COMM_WORLD, status, IERR)

 end subroutine pp_recv_d

!-------------------------------------------------------------------!

 subroutine pp_recv_d1d(val, rank, a)
   implicit none
   integer,                     intent(in)  :: a
   real(kind=dp), dimension(a), intent(out) :: val
   integer,                     intent(in)  :: rank

   integer, parameter :: tag = 1

   call mpi_recv(val, a, MPI_DOUBLE_PRECISION, rank, tag, MPI_COMM_WORLD, status, IERR)

 end subroutine pp_recv_d1d

!-------------------------------------------------------------------!

 subroutine pp_recv_d4d(val, rank, a, b, c, d)
   implicit none
   integer,                           intent(in)  :: a, b, c, d
   real(kind=dp), dimension(a,b,c,d), intent(out) :: val
   integer,                           intent(in)  :: rank

   integer :: ii, ij
   integer, parameter :: tag = 1

   do ii = 1, d
   do ij = 1, c
      call mpi_recv(val(:,:,ij,ii), a*b, MPI_DOUBLE_PRECISION, rank, tag, MPI_COMM_WORLD, status, IERR)
   enddo
   enddo

 end subroutine pp_recv_d4d

!-------------------------------------------------------------------!

 subroutine pp_barrier()
   implicit none

   call mpi_barrier(MPI_COMM_WORLD, IERR)

 end subroutine pp_barrier
end module parallel
