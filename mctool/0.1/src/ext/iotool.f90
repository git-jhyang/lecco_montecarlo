module iotool

use constants,  only:   dp,             &
                        LINELEN,        &
                        CHKOK,          &
                        CHKERR,         &
                        CHKEOF

implicit none

!--------------------------------------------------------------------------!
!       function and subroutines
!--------------------------------------------------------------------------!
!       = functions =
!       io_getunit            : get integer for file index, return 'u'
!       io_getcapital(chr_in) : capitalize alphabet, return 'chr_out'
!       io_center(chr_in,n)   : locate 'chr_in' at the center of 'chr_out'
!                                 with whole length of 'n', return 'chr_out'
!       io_isodd(a)           : if integer 'a' is odd, return true
!       io_assert_file(file)  : check 'file' exists, return true
!
!       = subroutines =
!       io_findline(u,a,b,i)  : find line that contains keyword 'a' in 
!                                 the file with index 'u' and return 'b' 
!                                 that whole line except keyword 'a'.
!                                 'i' is error code
!       io_printhead(a,b,u)   : print headline with keyword 'a' decorated
!                                 with symbol 'b' in file indexed by 'u'.
!                                 'b' and 'u' are optional
!       io_printerr(a,b,c)    : print error message about action 'a' with 
!                                 some detail 'b', and more detail 'c' and 
!                                 terminate program
!--------------------------------------------------------------------------!

public ::       io_getunit,           io_getcapital,        &
                io_center,            io_isodd,             &
                io_assert_file,       io_findline,          &
                io_printhead,         io_printerr,          &
                io_findstr,           io_stripline,         &
                io_checkdim,          io_linesplit

interface io_linesplit
   module procedure     io_linesplit_i,         io_linesplit_in,        &
                        io_linesplit_l,         io_linesplit_ln,        &
                        io_linesplit_r,         io_linesplit_rn,        &
                        io_linesplit_d,         io_linesplit_dn,        &
                        io_linesplit_c,         io_linesplit_cn

end interface io_linesplit

contains

!--------------------------------------------------------------------------!

 function io_getunit(u_try) result(u)
   implicit none
   integer, optional, intent(in) :: u_try

   integer :: u
   integer, parameter :: u_ini = 10
   logical :: isfile, isopen

   if (present(u_try)) then
      u = u_try
   else
      u = u_ini
   endif

   search : do
      inquire(unit = u, exist = isfile)
      if (isfile) then
         inquire (unit = u, opened = isopen)
         if (.not. isopen) exit search
      else
         exit search
      endif
      u = u + 1
      if (u .ge. 100) then
         call io_printerr('programming error','more than 90 files are opend')
         stop
      endif
   enddo search
 end function io_getunit

!--------------------------------------------------------------------------!

 function io_getcapital(chr_in) result(chr_out)
   implicit none
   character(LEN=*), intent(in) :: chr_in  
   
   character(LEN=len(chr_in))   :: chr_out
   integer                      :: ii, ij

   chr_out = chr_in
   do ii = 1, len_trim(chr_in)
      select case(iachar(chr_in(ii:ii)))
         case(97:122)
            chr_out(ii:ii) = achar(iachar(chr_in(ii:ii))-32)
      end select
   enddo
 end function io_getcapital

!--------------------------------------------------------------------------!

 function io_center(chr, n) result(chr_out)
   implicit none

   character(LEN=*), intent(in) :: chr
   integer,          intent(in) :: n

   character(LEN=n) :: chr_out
   integer          :: ii

   ii = len_trim(adjustl(chr)) 
   if (ii .ge. n) then
      chr_out = trim(adjustl(chr))
      return
   endif

   ii = (n-ii)/2

   chr_out = repeat(' ',ii)//trim(adjustl(chr))

 end function io_center

!--------------------------------------------------------------------------!

 function io_isodd(val) result(isodd)
   implicit none
   integer, intent(in) :: val

   logical :: isodd

   isodd = .false.
   if (mod(val,2) .eq. 1) isodd = .true.

 end function

!--------------------------------------------------------------------------!

 function io_assert_file(fname) result(isfile)
   implicit none
   character(LEN=*), intent(in)  :: fname

   logical :: isfile

   inquire(file = trim(adjustl(fname)), exist = isfile)
   if (.not. isfile) call io_printerr('File does not exists',trim(adjustl(fname)))

 end function io_assert_file

!--------------------------------------------------------------------------!

 subroutine io_findline(u_id, chr_in, chr_out, ierr)
   implicit none
   integer,          intent(in)  :: u_id
   character(LEN=*), intent(in)  :: chr_in
   character(LEN=*), intent(out) :: chr_out
   integer,          intent(out) :: ierr

   integer                :: ii, chk, id01, id02, id03
   character(LEN=LINELEN) :: line, line_head
   character(LEN=LINELEN) :: chr_tmp

   rewind(u_id)
   chr_tmp = trim(io_getcapital(adjustl(chr_in)))
   do ii = 1, LINELEN
      if (chr_tmp(ii:ii) .eq. ' ') then
         id02 = ii-1
         exit
      endif
   enddo
   if (id02 .eq. 0) then
      call io_printerr('io_findline','KEYWORD is not given')
      ierr = CHKERR
      return    ! keyword (chr_in) is not given, return 
   endif
   search : do
      read (u_id, '(A)', iostat = chk) line
      if (chk .eq. CHKOK) then
         if ((line(1:1) .eq. '#') .or. (line(1:1) .eq. '!')) cycle search
         line = trim(io_stripline(line))
         line_head = io_getcapital(line(1:id02))
                ! keyword region saved in 'line_head' (capitalized)
         if (line_head(1:id02) .ne. chr_tmp(1:id02)) cycle
                ! got keyword and save line w/o keyword
         write (chr_out,'(A)') trim(line(id02+1:))
         ierr = CHKOK
         return
      elseif (chk .ge. CHKERR) then
                ! error while read line, rarely occur
         call io_printerr('io_findline','error while finding KEYWORD '//trim(chr_tmp),trim(adjustl(line)))
         ierr = CHKERR
         return
      elseif (chk .le. CHKEOF) then
                ! touched end of file and fail to find keyword
         chr_out = ''
         ierr = CHKEOF
         return
      endif
   enddo search

 end subroutine io_findline

!--------------------------------------------------------------------------!

 function io_findstr(line, str) result(idx)
   implicit none
   character(LEN=*), intent(in) :: line
   character(LEN=*), intent(in) :: str

   character(LEN=len(line)) :: line_tmp
   integer :: iline, idx, lline, lstr

   line_tmp = io_stripline(line)//' '
   lline = len(trim(adjustl(line_tmp)))
   lstr  = len(trim(adjustl(str)))

   idx = 0

   if (lline .le. lstr) return
   if (lstr .eq. 0) return

   do iline = 1, lline-lstr+1
      if (line_tmp(iline:iline+lstr-1) .eq. str) then
         idx = iline
         return
      endif
   enddo

 end function io_findstr

!--------------------------------------------------------------------------!

 function io_stripline(line_in) result (line_out)
   implicit none
   character(LEN=*), intent(in) :: line_in
   character(LEN=len(line_in))  :: line_out

   character(LEN=len(line_in)) :: line_tmp
   integer :: i, iEnd

   line_tmp = trim(adjustl(line_in))
   iEnd = len(line_tmp)
   search : do i = 1, len(line_tmp)
      if ((line_tmp(i:i) .eq. '#') .or. (line_tmp(i:i) .eq. '!')) then
         iEnd = i-1
         exit search
      endif
   enddo search
   if (iEnd .eq. 0) return
   write (line_out,'(A)') line_tmp(1:iEnd)

 end function io_stripline

!--------------------------------------------------------------------------!

 subroutine io_checkdim(line, n, CHK, nMax)
   implicit none
   character(LEN=*),  intent(in)    :: line
   integer,           intent(in)    :: n
   integer,           intent(out)   :: CHK
   integer, optional, intent(inout) :: nMax

   character(LEN=len(line)) :: line_tmp
   integer                  :: i, iCnt

   line_tmp = io_stripline(line)
   iCnt     = 0

   do i = 1, len(line_tmp//' ')-1
      if ((line_tmp(i:i) .ne. ' ') .and. (line_tmp(i+1:i+1) .eq. ' ')) iCnt = iCnt + 1
   enddo

   CHK = CHKOK
   if (n .gt. iCnt) CHK = CHKERR
   if (present(nMax)) nMax = iCnt

 end subroutine io_checkdim

!--------------------------------------------------------------------------!

 subroutine io_linesplit_i(line_in, line_out, CHK, PrintErr_in)
   implicit none
   character(LEN=*), intent(in)  :: line_in
   integer,          intent(out) :: line_out
   integer,          intent(out) :: CHK
   logical, optional, intent(in) :: PrintErr_in

   character(LEN=len(line_in)) :: line_tmp
   character(LEN=len(line_in)) :: line_val
   integer                     :: i
   logical                     :: PrintErr

   PrintErr = .true.
   if (present(PrintErr_in)) PrintErr = PrintErr_in 

   line_tmp = trim(adjustl(io_stripline(line_in)))
   CHK      = CHKOK

   do i = 1, len(line_tmp//' ')-1
      if ((line_tmp(i:i) .ne. ' ') .and. (line_tmp(i+1:i+1) .eq. ' ')) then
         line_val = trim(adjustl(line_tmp(1:i)))
         read (line_val, *, err=100) line_out
         return
      endif
   enddo

100 continue
   if (PrintErr) call io_printerr('read error - io_linesplit',&
        trim(adjustl(line_val))//' is not INTEGER')
   CHK = CHKERR

 end subroutine io_linesplit_i

!--------------------------------------------------------------------------!

 subroutine io_linesplit_in(line_in, n, line_out, CHK, PrintErr_in)
   implicit none
   character(LEN=*),      intent(in)  :: line_in
   integer,               intent(in)  :: n
   integer, dimension(n), intent(out) :: line_out
   integer,               intent(out) :: CHK
   logical, optional,     intent(in)  :: PrintErr_in

   character(LEN=len(line_in)) :: line_tmp
   character(LEN=len(line_in)) :: line_val
   integer                     :: i, iCnt, id01, id02
   logical                     :: PrintErr

   PrintErr = .true.
   if (present(PrintErr_in)) PrintErr = PrintErr_in 
 
   line_tmp = trim(adjustl(io_stripline(line_in)))

   id01 = 1
   iCnt = 0
   CHK  = CHKOK

   do i = 1, len(line_tmp//' ')-1
      if ((line_tmp(i:i) .ne. ' ') .and. (line_tmp(i+1:i+1) .eq. ' ')) then
         iCnt = iCnt + 1
         id02 = i
         line_val = trim(adjustl(line_tmp(id01:id02)))
         read (line_val, *, err=100) line_out(iCnt)
         if (iCnt .eq. n) return
         id01 = id02 + 1
      endif
   enddo

100 continue
   if (PrintErr) call io_printerr('read error - io_linesplit',&
        trim(adjustl(line_val))//' is not INTEGER')
   CHK = CHKERR

 end subroutine io_linesplit_in

!--------------------------------------------------------------------------!

 subroutine io_linesplit_l(line_in, line_out, CHK, PrintErr_in)
   implicit none
   character(LEN=*), intent(in)  :: line_in
   logical,          intent(out) :: line_out
   integer,          intent(out) :: CHK
   logical, optional, intent(in) :: PrintErr_in

   character(LEN=len(line_in)) :: line_tmp
   character(LEN=len(line_in)) :: line_val
   integer                     :: i
   logical                     :: PrintErr

   PrintErr = .true.
   if (present(PrintErr_in)) PrintErr = PrintErr_in 
 
   line_tmp = trim(adjustl(io_stripline(line_in)))
   CHK      = CHKOK

   do i = 1, len(line_tmp//' ')-1
      if ((line_tmp(i:i) .ne. ' ') .and. (line_tmp(i+1:i+1) .eq. ' ')) then
         line_val = trim(adjustl(line_tmp(1:i)))
         read (line_val, *, err=100) line_out
         return
      endif
   enddo

100 continue
   if (PrintErr) call io_printerr('read error - io_linesplit',&
        trim(adjustl(line_val))//' is not BOOLEAN')
   CHK = CHKERR

 end subroutine io_linesplit_l

!--------------------------------------------------------------------------!

 subroutine io_linesplit_ln(line_in, n, line_out, CHK, PrintErr_in)
   implicit none
   character(LEN=*),      intent(in)  :: line_in
   integer,               intent(in)  :: n
   logical, dimension(n), intent(out) :: line_out
   integer,               intent(out) :: CHK
   logical, optional,     intent(in)  :: PrintErr_in

   character(LEN=len(line_in)) :: line_tmp
   character(LEN=len(line_in)) :: line_val
   integer                     :: i, iCnt, id01, id02
   logical                     :: PrintErr

   PrintErr = .true.
   if (present(PrintErr_in)) PrintErr = PrintErr_in 
 
   line_tmp = trim(adjustl(io_stripline(line_in)))

   id01 = 1
   iCnt = 0
   CHK  = CHKOK

   do i = 1, len(line_tmp//' ')-1
      if ((line_tmp(i:i) .ne. ' ') .and. (line_tmp(i+1:i+1) .eq. ' ')) then
         iCnt = iCnt + 1
         id02 = i
         line_val = trim(adjustl(line_tmp(id01:id02)))
         read (line_val, *, err=100) line_out(iCnt)
         if (iCnt .eq. n) return
         id01 = id02 + 1
      endif
   enddo

100 continue
   if (PrintErr) call io_printerr('read error - io_linesplit',&
        trim(adjustl(line_val))//' is not BOOLEAN')
   CHK = CHKERR

 end subroutine io_linesplit_ln

!--------------------------------------------------------------------------!

 subroutine io_linesplit_r(line_in, line_out, CHK, PrintErr_in)
   implicit none
   character(LEN=*), intent(in)  :: line_in
   real,             intent(out) :: line_out
   integer,          intent(out) :: CHK
   logical, optional, intent(in) :: PrintErr_in

   character(LEN=len(line_in)) :: line_tmp
   character(LEN=len(line_in)) :: line_val
   integer                     :: i
   logical                     :: PrintErr

   PrintErr = .true.
   if (present(PrintErr_in)) PrintErr = PrintErr_in 
 
   line_tmp = trim(adjustl(io_stripline(line_in)))
   CHK      = CHKOK

   do i = 1, len(line_tmp//' ')-1
      if ((line_tmp(i:i) .ne. ' ') .and. (line_tmp(i+1:i+1) .eq. ' ')) then
         line_val = trim(adjustl(line_tmp(1:i)))
         read (line_val, *, err=100) line_out
         return
      endif
   enddo

100 continue
   if (PrintErr) call io_printerr('read error - io_linesplit',&
        trim(adjustl(line_val))//' is not REAL')
   CHK = CHKERR

 end subroutine io_linesplit_r

!--------------------------------------------------------------------------!

 subroutine io_linesplit_rn(line_in, n, line_out, CHK, PrintErr_in)
   implicit none
   character(LEN=*),   intent(in)  :: line_in
   integer,            intent(in)  :: n
   real, dimension(n), intent(out) :: line_out
   integer,            intent(out) :: CHK
   logical, optional,  intent(in)  :: PrintErr_in

   character(LEN=len(line_in)) :: line_tmp
   character(LEN=len(line_in)) :: line_val
   integer                     :: i, iCnt, id01, id02
   logical                     :: PrintErr

   PrintErr = .true.
   if (present(PrintErr_in)) PrintErr = PrintErr_in 
 
   line_tmp = trim(adjustl(io_stripline(line_in)))

   id01 = 1
   iCnt = 0
   CHK  = CHKOK

   do i = 1, len(line_tmp//' ')-1
      if ((line_tmp(i:i) .ne. ' ') .and. (line_tmp(i+1:i+1) .eq. ' ')) then
         iCnt = iCnt + 1
         id02 = i
         line_val = trim(adjustl(line_tmp(id01:id02)))
         read (line_val, *, err=100) line_out(iCnt)
         if (iCnt .eq. n) return
         id01 = id02 + 1
      endif
   enddo

100 continue
   if (PrintErr) call io_printerr('read error - io_linesplit',&
        trim(adjustl(line_val))//' is not REAL')
   CHK = CHKERR

 end subroutine io_linesplit_rn

!--------------------------------------------------------------------------!

 subroutine io_linesplit_d(line_in, line_out, CHK, PrintErr_in)
   implicit none
   character(LEN=*), intent(in)  :: line_in
   real(kind=dp),    intent(out) :: line_out
   integer,          intent(out) :: CHK
   logical, optional, intent(in) :: PrintErr_in

   character(LEN=len(line_in)) :: line_tmp
   character(LEN=len(line_in)) :: line_val
   integer                     :: i
   logical                     :: PrintErr

   PrintErr = .true.
   if (present(PrintErr_in)) PrintErr = PrintErr_in 
 
   line_tmp = trim(adjustl(io_stripline(line_in)))
   CHK      = CHKOK

   do i = 1, len(line_tmp//' ')-1
      if ((line_tmp(i:i) .ne. ' ') .and. (line_tmp(i+1:i+1) .eq. ' ')) then
         line_val = trim(adjustl(line_tmp(1:i)))
         read (line_val, *, err=100) line_out
         return
      endif
   enddo

100 continue
   if (PrintErr) call io_printerr('read error - io_linesplit',&
        trim(adjustl(line_val))//' is not DOUBLE')
   CHK = CHKERR

 end subroutine io_linesplit_d

!--------------------------------------------------------------------------!

 subroutine io_linesplit_dn(line_in, n, line_out, CHK, PrintErr_in)
   implicit none
   character(LEN=*),            intent(in)  :: line_in
   integer,                     intent(in)  :: n
   real(kind=dp), dimension(n), intent(out) :: line_out
   integer,                     intent(out) :: CHK
   logical, optional,           intent(in)  :: PrintErr_in

   character(LEN=len(line_in)) :: line_tmp
   character(LEN=len(line_in)) :: line_val
   integer                     :: i, iCnt, id01, id02
   logical                     :: PrintErr

   PrintErr = .true.
   if (present(PrintErr_in)) PrintErr = PrintErr_in 
 
   line_tmp = trim(adjustl(io_stripline(line_in)))

   id01 = 1
   iCnt = 0
   CHK  = CHKOK

   do i = 1, len(line_tmp//' ')-1
      if ((line_tmp(i:i) .ne. ' ') .and. (line_tmp(i+1:i+1) .eq. ' ')) then
         iCnt = iCnt + 1
         id02 = i
         line_val = trim(adjustl(line_tmp(id01:id02)))
         read (line_val, *, err=100) line_out(iCnt)
         if (iCnt .eq. n) return
         id01 = id02 + 1
      endif
   enddo

100 continue
   if (PrintErr) call io_printerr('read error - io_linesplit',&
        trim(adjustl(line_val))//' is not DOUBLE')
   CHK = CHKERR

 end subroutine io_linesplit_dn

!--------------------------------------------------------------------------!

 subroutine io_linesplit_c(line_in, line_out, CHK)
   implicit none
   character(LEN=*), intent(in)  :: line_in
   character(LEN=*), intent(out) :: line_out
   integer,          intent(out) :: CHK

   character(LEN=len(line_in)) :: line_tmp
   character(LEN=len(line_in)) :: line_val
   integer                     :: i
 
   line_tmp = trim(adjustl(io_stripline(line_in)))
   CHK      = CHKOK

   do i = 1, len(line_tmp//' ')-1
      if ((line_tmp(i:i) .ne. ' ') .and. (line_tmp(i+1:i+1) .eq. ' ')) then
         line_val = trim(adjustl(line_tmp(1:i)))
         read (line_val, '(A)') line_out
         return
      endif
   enddo

 end subroutine io_linesplit_c

!--------------------------------------------------------------------------!

 subroutine io_linesplit_cn(line_in, n, line_out, CHK)
   implicit none
   character(LEN=*),               intent(in)  :: line_in
   integer,                        intent(in)  :: n
   character(LEN=*), dimension(n), intent(out) :: line_out
   integer,                        intent(out) :: CHK

   character(LEN=len(line_in)) :: line_tmp
   character(LEN=len(line_in)) :: line_val
   integer                     :: i, iCnt, id01, id02
 
   line_tmp = trim(adjustl(io_stripline(line_in)))

   id01 = 1
   iCnt = 0
   CHK  = CHKOK

   do i = 1, len(line_tmp//' ')-1
      if ((line_tmp(i:i) .ne. ' ') .and. (line_tmp(i+1:i+1) .eq. ' ')) then
         iCnt = iCnt + 1
         id02 = i
         line_val = trim(adjustl(line_tmp(id01:id02)))
         read (line_val, '(A)') line_out(iCnt)
         if (iCnt .eq. n) return
         id01 = id02 + 1
      endif
   enddo

 end subroutine io_linesplit_cn

!--------------------------------------------------------------------------!

 subroutine io_printhead(chr_in, symb, u)
   implicit none
   character(LEN=*), intent(in)           :: chr_in
   character,        intent(in), optional :: symb
   integer,          intent(in), optional :: u

   character :: chr

   chr = '='
   if (present(symb)) chr = symb
 
   if (present(u)) then
      write (u,*)
      write (u,*) repeat(chr,75)
      write (u,*) io_center(chr_in, 75)
      write (u,*) repeat(chr,75)
      write (u,*)
   else
      write (*,*)
      write (*,*) repeat(chr,75)
      write (*,*) io_center(chr_in, 75)
      write (*,*) repeat(chr,75)
      write (*,*)
   endif

 end subroutine io_printhead

!--------------------------------------------------------------------------!

 subroutine io_printerr(str1, str2, line)
   implicit none

   character(LEN=*),           intent(in) :: str1
   character(LEN=*),           intent(in) :: str2
   character(LEN=*), optional, intent(in) :: line

   write (*,'(1x,"!!! ERROR : ",A," : ",A)') trim(adjustl(str1)), trim(adjustl(str2))
   if (present(line)) write (*,'(13x,A)') trim(adjustl(line))

 end subroutine io_printerr

!--------------------------------------------------------------------------!

end module iotool

