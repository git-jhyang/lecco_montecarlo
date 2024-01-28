module vasp

use constants,  only :  dp, LINELEN, CHKEOF, CHKOK

use util,       only :  cart_to_dir

use iotool,     only :  io_checkdim,    &
                        io_linesplit,   &
                        io_printerr

implicit none

type, public :: poscar
   character(LEN=LINELEN)                        :: Des
   real(kind=dp),    dimension(3,3)              :: LattVec
   integer                                       :: nSpec = 0
   integer                                       :: nAtom
   integer,          dimension(:),   allocatable :: nAtomSpec
   character(LEN=2), dimension(:),   allocatable :: NameAtoms
   logical                                       :: isCart = .false.
   logical                                       :: isSelc = .false.
   real(kind=dp),    dimension(:,:), allocatable :: coo
   logical,          dimension(:,:), allocatable :: selc
   real(kind=dp),    dimension(:,:), allocatable :: force
end type poscar

type, public :: chgcar
   integer,       dimension(3)                  :: nGrid
   logical                                      :: isLong
   real(kind=dp), dimension(:,:,:), allocatable :: chg
   real(kind=dp)                                :: eTime
end type chgcar

contains

!-----------------------------------------------------------------------------!

 subroutine chg_get(u, DatChg)
   implicit none
   integer,      intent(in)    :: u
   type(chgcar), intent(inout) :: DatChg

   character(LEN=LINELEN) :: line
   integer                :: CHK, nLine, nCol, nPnt, nRem
   integer                :: ii, ix1, ix2, ix3, ix4, iy1, iy2, iz1, iz2
   integer                :: t1, t2, trate
   logical                :: isOpen

   nCol = 0

   call system_clock(t1, trate)
 
   read (u,'(A)') line
    call io_checkdim(line, 1, CHK, nMax=nCol)
    if (nCol .ne. 5) DatChg%isLong = .false.
    if (nCol .ne. 10) DatChg%isLong = .true.

   nPnt  = DatChg%nGrid(1)*DatChg%nGrid(2)*DatChg%nGrid(3)
   nLine = ceiling(float(nPnt)/float(nCol))

   ix1 = 1
   ix3 = 1
   iy1 = 1
   iz1 = 1
   iz2 = 1
   do ii = 2, nLine
      ix2 = ix1 + nCol - 1
      if (ix2 .le. DatChg%nGrid(1)) then
         read (line, *) DatChg%chg(ix1:ix2,iy1,iz1)
      else
         ix4 = ix2 - DatChg%nGrid(1) 
         ix2 = DatChg%nGrid(1)
         iy2 = iy1 + 1
         if (iy1 .eq. DatChg%nGrid(2)) then
            iy2 = 1
            iz2 = iz1 + 1
         endif
         read (line, *) DatChg%chg(ix1:ix2,iy1,iz1), DatChg%chg(ix3:ix4,iy2,iz2)
         ix2 = ix4
         iy1 = iy2
         iz1 = iz2
      endif
      ix1 = ix2 + 1
      if (ix2 .eq. DatChg%nGrid(1)) then
         ix1 = 1
         iy1 = iy1 + 1
         if (iy1 .gt. DatChg%nGrid(2)) then
            iy1 = 1
            iz1 = iz1 + 1
            iz2 = iz1
         endif
      endif
      read (u,'(A)') line
   enddo
   read (line, *) DatChg%chg(ix1:DatChg%nGrid(1),iy1,iz1)

   call system_clock(t2, trate)
   DatChg%eTime = float(t2-t1)/float(trate)

 end subroutine chg_get

!-----------------------------------------------------------------------------!

 subroutine chg_write(u, DatChg)
   implicit none
   integer,           intent(in)    :: u
   type(chgcar),      intent(inout) :: DatChg

   integer :: t1, t2, trate
   integer :: ix, iy, iz
   logical :: isLong
   character(LEN=LINELEN) :: form
 
   call system_clock(t1, trate)

   if (DatChg%islong) then
      form = '(5(1X,E17.11))'
   else
      form = '(10(1X,G11.5))'
   endif

   write (u,'(3I5)') DatChg%nGrid
   write (u, form) (((DatChg%chg(ix,iy,iz), ix=1,DatChg%nGrid(1)), &
                   iy=1,DatChg%nGrid(2)), iz=1,DatChg%nGrid(3))

   call system_clock(t2, trate)
   DatChg%eTime = float(t2-t1)/float(trate)

 end subroutine chg_write

!-----------------------------------------------------------------------------!

 subroutine chg_sum(DatChg1, DatChg2, DatChgOut, Mult)
   implicit none
   type(chgcar),                          intent(in)    :: DatChg1
   type(chgcar),                          intent(in)    :: DatChg2
   type(chgcar),                          intent(inout) :: DatChgOut
   real(kind=dp), dimension(2), optional, intent(in)    :: Mult

   real(kind=dp), dimension(2) :: M
   integer                     :: ix, iy, iz
   integer                     :: t1, t2, trate

   call system_clock(t1, trate)

   M(:) = 1d0

   if (present(Mult)) M = Mult

   do iz = 1, DatChgOut%nGrid(3)
   do iy = 1, DatChgOut%nGrid(2)
   do ix = 1, DatChgOut%nGrid(1)
      DatChgOut%chg(ix,iy,iz) = M(1)*DatChg1%chg(ix,iy,iz) + M(2)*DatChg2%chg(ix,iy,iz) 
   enddo
   enddo
   enddo

   call system_clock(t2, trate)
   DatChgOut%eTime = float(t2-t1)/float(trate)

 end subroutine chg_sum

!-----------------------------------------------------------------------------!

 subroutine poscar_read(u, DatPos)
   implicit none
   integer,      intent(in)  :: u
   type(poscar), intent(out) :: DatPos

   integer                              :: ii, IOCHK
   logical                              :: iserr, iscart
   character(LEN=LINELEN)               :: line
   character(LEN=LINELEN), dimension(6) :: line6
   real(kind=dp)                        :: LattCont

   CHK = CHKOK

   read (u,'(A)', iostat=IOCHK, end=100) DatPos%Des ! file comment
   read (u,'(A)', iostat=IOCHK, end=100) line
    call io_checkdim(line, 1, CHK)
    if (CHK .ne. CHKOK) goto 100
    call io_linesplit(line, 1, LattCont)

   do ii = 1, 3
      read (u,'(A)', iostat=IOCHK, end=100) line
       call io_checkdim(line, 3, CHK)
       if (CHK .ne. CHKOK) goto 100
       call io_linesplit(line, 3, DatPos%LattVec(:,ii))
       DatPos%LattVec(:,ii) = DatPos%LattVec(:,ii)*LattCont
   enddo

   read (u,'(A)', iostat=IOCHK, end=100) line ! name of atoms, vasp 5 style
    call io_checkdim(line, 1, CHK, nMax=DatPos%nSpec)
    if (CHK .ne. CHKOK) goto 100
    allocate (DatPos%NameAtoms(DatPos%nSpec), DatPos%nAtomSpec(DatPos%nSpec))
    call io_linesplit(line, DatPos%nSpec, DatPos%NameAtoms)

   read(u,'(A)', iostat=IOCHK, end=100) line  ! number of atoms of each species
    call io_checkdim(line, 1, CHK, nMax=DatPos%nSpec)
    if (CHK .ne. CHKOK) goto 100
    call io_linesplit(line, DatPos%nSpec, DatPos%nAtomSpec)
    DatPos%nAtom = sum(DatPos%nAtomSpec)
    allocate (DatPos%coo(3,DatPos%nAtom))
    allocate (DatPos%force(3,DatPos%nAtom))

   read(u,'(A)', iostat=IOCHK, end=100) line  ! read coordination type
    line = io_getcapital(line)
    if (line(1:1) .eq. 'S') then
       DatPos%isSelc = .true.
       read (u, '(A)', iostat=IOCHK, end=100) line
       line = io_getcapital(line)
       allocate (DatPos%selc(3,DatPos%nAtom))
    endif
    if (line(1:1) .eq. 'C') DatPos%isCart = .true.

   do ii = 1, DatPos%nAtom      ! read coordination
      read(u,'(A)', iostat=IOCHK, end=100) line
       if (DatPos%isSelc) then
          call io_checkdim(line, 6, CHK)
          if (CHK .ne. CHKOK) goto 100
          call io_linesplit(line, 6, line6)
          line = trim(line6(1))//' '//trim(line6(2))//' '//trim(line6(3))
          call io_linesplit(line, 3, DatPos%coo(:,ii))
          line = trim(line6(4))//' '//trim(line6(5))//' '//trim(line6(6))
          call io_linesplit(line, 3, DatPos%selc(:,ii))
       else
          call io_checkdim(line, 3, CHK)
          if (CHK .ne. CHKOK) goto 100
          call io_linesplit(line, 3, DatPos%coo(:,ii))
       endif
       if (DatPos%isCart) call cart_to_dir(DatPos%LattVec, DatPos%coo(:,ii))
   enddo
   DatPos%force(:,:) = 0d0

   read (u, *, iostat=IOCHK, end=101)
   do ii = 1, DatPos%nAtom
      read (u, '(A)', iostat=IOCHK, end=101) line
       call io_checkdim(line, 3, CHK)
       if (CHK .ne. CHKOK) goto 100
       call io_linesplit(line, 3, DatPos%force(:,ii))
   enddo

   return

100 continue
   if (IOCHK .le. CHKEOF) call io_printerr('poscar_read - I/O error','touched end of file')
   if (IOCHK .ge. CHKERR) call io_printerr('poscar_read - I/O error','variable error',line)
   if (CHK .ge. CHKERR) call io_printerr('poscar_read - io_checkdim','dimension is not coorrect',line)
   return

101 continue
   return
 end subroutine poscar_read

!-----------------------------------------------------------------------------!

 subroutine poscar_write(u, DatPos, islong_in)
   implicit none
   integer,           intent(in) :: u
   type(poscar),      intent(in) :: DatPos
   logical, optional, intent(in) :: islong_in

   integer                     :: ii
   logical                     :: islong
   character(LEN=LINELEN)      :: form
   real(kind=dp), dimension(3) :: coo

   islong = .true.
   if (present(islong_in)) islong = islong_in

   write (u,'(A)') trim(adjustl(DatPos%Des))
   write (u,'(A)') '1.00000000000'

   if (islong) then
      form='(1x,3F22.16)'
   else
      form='(1x,3F12.6)'
   endif

   write (u, form) (DatPos%LattVec(:,ii),ii=1,3)

   do ii = 1, DatPos%nSpec
      write (u,'(A4)',advance='no') DatPos%NameAtoms(ii)
   enddo
   write (u,*)
   do ii = 1, DatPos%nSpec
      write (u,'(I4)',advance='no') DatPos%nAtomSpec(ii)
   enddo
   write (u,*)
   write (u,'(A)') 'Direct'

   if (islong) then
      form='(3F20.16)'
   else
      form='(3F10.6)'
   endif

   do ii = 1, DatPos%nAtom
      coo = DatPos%coo(:,ii)
      if (DatPos%isCart) call cart_to_dir(DatPos%LattVec, coo)
      write (u, form) coo
   enddo

 end subroutine poscar_write

!-----------------------------------------------------------------------------!

end module vasp

