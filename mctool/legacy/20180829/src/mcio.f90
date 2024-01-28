module mcio

use constants,  only:   LINELEN,        &
                        CHKOK,          &
                        CHKERR,         &
                        CHKEOF

implicit none

!--------------------------------------------------------------------------!
!       function and subroutines
!--------------------------------------------------------------------------!
!       = functions =
!       mcio_getunit            : get integer for file index, return 'u'
!       mcio_getcapital(chr_in) : capitalize alphabet, return 'chr_out'
!       mcio_center(chr_in,n)   : locate 'chr_in' at the center of 'chr_out'
!                                 with whole length of 'n', return 'chr_out'
!       mcio_isodd(a)           : if integer 'a' is odd, return true
!       mcio_assert_file(file)  : check 'file' exists, return true
!
!       = subroutines =
!       mcio_findline(u,a,b,i)  : find line that contains keyword 'a' in 
!                                 the file with index 'u' and return 'b' 
!                                 that whole line except keyword 'a'.
!                                 'i' is error code
!       mcio_printhead(a,b,u)   : print headline with keyword 'a' decorated
!                                 with symbol 'b' in file indexed by 'u'.
!                                 'b' and 'u' are optional
!       mcio_printerr(a,b,c)    : print error message about action 'a' with 
!                                 some detail 'b', and more detail 'c' and 
!                                 terminate program
!--------------------------------------------------------------------------!

public ::       mcio_getunit,           mcio_getcapital,        &
                mcio_center,            mcio_isodd,             &
                mcio_assert_file,       mcio_findline,          &
                mcio_printhead,         mcio_printerr,          &
                mcio_findstr

contains

!--------------------------------------------------------------------------!

 function mcio_getunit(u_try) result(u)
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
         write (*,*) " !!! ERROR : Fail to fine unused 'unit' (u < 100) !!! "
         stop
      endif
   enddo search
 end function mcio_getunit

!--------------------------------------------------------------------------!

 function mcio_getcapital(chr_in) result(chr_out)
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
 end function mcio_getcapital

!--------------------------------------------------------------------------!

 function mcio_center(chr, n) result(chr_out)
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

   chr_out = repeat(' ',ii) // trim(adjustl(chr))

 end function mcio_center

!--------------------------------------------------------------------------!

 function mcio_isodd(val) result(isodd)
   implicit none
   integer, intent(in) :: val

   logical :: isodd

   isodd = .false.
   if (mod(val,2) .eq. 1) isodd = .true.

 end function

!--------------------------------------------------------------------------!

 function mcio_assert_file(fname) result(isfile)
   implicit none
   character(LEN=*), intent(in)  :: fname

   logical :: isfile

   inquire(file = trim(adjustl(fname)), exist = isfile)
   if (.not. isfile) then
      write (*,*) " !!! ERROR : File does not exist : ", trim(adjustl(fname))
   endif

 end function mcio_assert_file

!--------------------------------------------------------------------------!

 subroutine mcio_findline(u_id, chr_in, chr_out, ierr)
   implicit none
   integer,          intent(in)  :: u_id
   character(LEN=*), intent(in)  :: chr_in
   character(LEN=*), intent(out) :: chr_out
   integer,          intent(out) :: ierr

   integer                :: ii, chk, id01, id02
   character(LEN=LINELEN) :: line, line_head
   character(LEN=LINELEN) :: chr_tmp

   rewind(u_id)
   chr_tmp = adjustl(mcio_getcapital(chr_in))
   do ii = 1, LINELEN
      if (chr_tmp(ii:ii) .eq. ' ') then
         id02 = ii-1
         exit
      endif
   enddo
   if (id02 .eq. 0) then
      write (*,*) " !!! ERROR : 'chr_in' is missing (mcio_findline) !!!"
      line = 'ERROR - ERROR - ERROR'
      ierr = CHKERR
      return    ! keyword (chr_in) is not given, return 
   endif
   search : do
      read(u_id, '(A)', iostat = chk) line
      if (line(1:1) .eq. '#' .or. line(1:1) .eq. '!') cycle
      if (chk .eq. CHKOK) then
         line = trim(adjustl(line))
         line_head = mcio_getcapital(line(1:id02))
                ! keyword region saved in 'line_head' (capitalized)
         if (line_head(1:id02) .ne. chr_tmp(1:id02)) cycle
         if (line_head(1:id02) .eq. chr_tmp(1:id02)) then
                ! got keyword and save line w/o keyword
            write (chr_out,'(A)') trim(adjustl(line(id02+1:)))
            ierr = CHKOK
            exit search
         endif
      elseif (chk .ge. 1) then
                ! error while read line, rarely occur
         write (*,*) ' !!! ERROR : error while reading (mcio_findline) : ', &
                        trim(adjustl(chr_tmp))
         write (*,*) '           ',trim(adjustl(line))
         ierr = CHKERR
         exit search
      elseif (chk .le. -1) then
                ! touched end of file and fail to find keyword
         ierr = CHKEOF
         exit search
      endif
   enddo search

 end subroutine mcio_findline

!--------------------------------------------------------------------------!

 function mcio_findstr(line, str) result(idx)
   implicit none
   character(LEN=*), intent(in) :: line
   character(LEN=*), intent(in) :: str

   integer :: iline, idx, lline, lstr

   lline = len(adjustl(trim(line)))
   lstr  = len(adjustl(trim(str)))

   idx = 0

   if (lline .le. lstr) return
   if (lstr .eq. 0) return

   do iline = 1, lline-lstr+1
      if (line(iline:iline+lstr-1) .eq. str) then
         idx = iline
         return
      endif
   enddo

 end function mcio_findstr

!--------------------------------------------------------------------------!

 subroutine mcio_printhead(chr_in, symb, u)
   implicit none
   character(LEN=*), intent(in)           :: chr_in
   character,        intent(in), optional :: symb
   integer,          intent(in), optional :: u

   character :: chr

   chr = '='
   if (present(symb)) chr    = symb
 
   if (present(u)) then
      write (u,*)
      write (u,*) repeat(chr,75)
      write (u,*) mcio_center(chr_in, 75)
      write (u,*) repeat(chr,75)
      write (u,*)
   else
      write (*,*)
      write (*,*) repeat(chr,75)
      write (*,*) mcio_center(chr_in, 75)
      write (*,*) repeat(chr,75)
      write (*,*)
   endif

 end subroutine mcio_printhead

!--------------------------------------------------------------------------!

 subroutine mcio_printerr(act, chr_in, line)
   implicit none

   character(LEN=*),           intent(in) :: act
   character(LEN=*),           intent(in) :: chr_in
   character(LEN=*), optional, intent(in) :: line

   character(LEN=LINELEN) :: act_adj

   act_adj = trim(adjustl(act))

   write (*,200) ' !!! ERROR : ',trim(act_adj), ' error : ', trim(adjustl(chr_in))
   if (present(line)) write (*,*) trim(line)

200 format(1x,4A)

 end subroutine mcio_printerr

!--------------------------------------------------------------------------!

end module mcio

