program main_gen

use constants,  only :  PATHLEN, VARLEN,        &
                        CHKOK, CHKERR,          &
                        dp, kb_ev

use iotool,     only :  io_getunit,           &
                        io_assert_file,       &
                        io_printhead,         &
                        io_printerr

use parallel,   only :  pMast,          &
                        nprc,           &
                        pp_init,        &
                        pp_fin,         &
                        pp_bcast

use input,      only :  Var_Input,      &
                        inp_init,       &
                        inp_fin

use coordinate, only :  Var_Coord,      &
                        coo_init,       &
                        coo_fin,        &
                        coo_writeposcar

use generator,  only :  Var_Gen,        &
                        gen_init,       &
                        gen_load_coo
!                        gen_fin,        &

implicit none

!-------------------------------------------------------------------!
!       Variables
!-------------------------------------------------------------------!

type(Var_Input) :: vinp
type(Var_Coord) :: vcoo
type(Var_Gen)   :: vgen

integer :: u, t_start, t_end, t_rate, t_max
integer :: t_h, t_m, t_s, t_ms

logical :: main_isinit = .false.
logical :: main_iserr  = .true.

real(kind=dp), dimension(3) :: time_tot

character(LEN=PATHLEN)          :: f_name
character(LEN=VARLEN)           :: f_comm
character(LEN=14), dimension(3) :: time_form 

!-------------------------------------------------------------------!

call initialize()

if (.not. main_isinit) then
   if (pMast) call io_printerr('main initialize','error occurred')
   main_iserr = .true.
   call finalize()
endif

do
   call gen_load_coo(vgen)
   if (vgen%isfin) exit
   if (vgen%iserr) then
      call io_printerr('generator failed','error occured')
      main_iserr = .true.
      exit
   endif
   write (f_name,'("outputs/",A,"_GEN_",I6.6)') trim(vinp%fname_coo_out), vgen%iGrid
   write (f_comm,'("GEN_",I6.6)') vgen%iGrid
   call coo_writeposcar(f_name, vcoo%nSpec, vgen%nAtom, f_comm, vgen%LattVec, &
                        vcoo%NameAtom, vgen%nAtomSpec, .true., vgen%coo)
enddo

call main_calctime()

write (*,'(3x,"Total Number of Coordinate     : ",I8)') vgen%Cnt_all
write (*,'(3x,"Number of Generated Coordinate : ",I8)') vgen%iGrid

call finalize()

contains

!-------------------------------------------------------------------!

 subroutine initialize()
   implicit none

   if (main_isinit) return

   call pp_init()
   if (nprc .ne. 1) then
      call io_printerr('parallelization','program designed for serial run')
      return
   endif

   if (pMast) then      ! initialize program timer
      time_tot(:) = 0d0
      call system_clock(t_start, t_rate, t_max)
   endif

   vinp = inp_init('GEN')
   if (vinp%iserr) return

   vcoo = coo_init(vinp)
   if (vcoo%iserr) return

   vgen = gen_init(vinp, vcoo)
   if (vgen%iserr) return

   if (pMast) then
      call main_calctime()
      write (*,'(2x,"Elapsed time for initialization : ",A14)') time_form(3)
   endif

   main_iserr  = .false.
   main_isinit = .true.

 end subroutine initialize

!-------------------------------------------------------------------!

 subroutine finalize()
   implicit none

   if (.not. vinp%iserr) call inp_fin(vinp)
   if (.not. vcoo%iserr) call coo_fin(vcoo)
!   if (.not. vopt%iserr) call opt_fin(vgen)

   if (pMast) then
      if (main_iserr) then
         call io_printhead('ERROR - GENERATOR TERMINATED - ERROR')
      else
         write (*,'(2x,"Elapsed Time : ",A13)') time_form(3)
         call io_printhead('GENERATION FINISHED')
      endif
   endif
   call pp_fin()
   stop

 end subroutine finalize

!-------------------------------------------------------------------!

 subroutine main_calctime()
   implicit none
   integer :: ii

   call system_clock(t_end, t_rate, t_max)
   time_tot(1) = real(t_end - t_start)/real(t_rate)
   if (time_tot(1) .lt. 0d0) time_tot(1) = time_tot(1) + real(t_max)/real(t_rate)
   time_tot(2) = time_tot(2) + time_tot(1)
   time_tot(3) = time_tot(3) + time_tot(2)

   do ii = 1, 3
      t_h  = floor(time_tot(ii)/3600)
      t_m  = floor(time_tot(ii)/60 - t_h*60)
      t_s  = floor(time_tot(ii) - t_m*60 - t_h*360)
      t_ms = floor((time_tot(ii) - floor(time_tot(ii)))*100)
      if (t_h .ne. 0) then
         write (time_form(ii),'(I4,"h",I2.2,"m",I2.2,".",I2.2,"s")') t_h, t_m, t_s, t_ms
      elseif (t_m .ne. 0) then
         write (time_form(ii),'(5x, I2,"m",I2.2,".",I2.2,"s")') t_m, t_s, t_ms
      else
         write (time_form(ii),'(8x, I2,".",I2.2,"s")') t_s, t_ms
      endif
   enddo
   call system_clock(t_start, t_rate, t_max)

 end subroutine main_calctime

!-------------------------------------------------------------------!

end program main_gen
