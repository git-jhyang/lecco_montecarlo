program main_opt

use constants,  only :  PATHLEN, VARLEN,        &
                        CHKOK, CHKERR,          &
                        dp, kb_ev

use mt_serial,  only :  grnd

use mcio,       only :  mcio_getunit,           &
                        mcio_assert_file,       &
                        mcio_printhead,         &
                        mcio_printerr

use parallel,   only :  pMast,          &
                        pp_mtserial,    &
                        pp_init,        &
                        pp_fin,         &
                        pp_bcast

use input,      only :  Var_Input,      &
                        inp_init,       &
                        inp_fin

use lattice,    only :  Var_Latt,       &
                        lat_init,       &
                        lat_fin

use coordinate, only :  Var_Coord,      &
                        coo_init,       &
                        coo_fin,        &
                        coo_writeposcar

use optimizer,  only :  Var_Opt,        &
                        opt_init,       &
                        opt_fin,        &
                        opt_ranwalk,    &
                        opt_coo_gather, &
                        opt_coo_update, &
                        opt_coo_rollback

use energy,     only :  ene_init,    &
                        ene_fin,     &
                        calc_energy

implicit none

!-------------------------------------------------------------------!
!       Variables
!-------------------------------------------------------------------!

type(Var_Input) :: vinp
type(Var_Coord) :: vcoo
type(Var_Latt)  :: vlat
type(Var_Opt)   :: vopt

integer :: iOptSite, iOptMig, iOptTry
integer :: rSite, rUnit, iSite
integer :: u, t_start, t_end, t_rate, t_max

logical :: main_isinit = .false.
logical :: main_iserr  = .true.
logical :: ene_iserr   = .false.

logical :: accept_ranwalk

real(kind=dp), dimension(3) :: time_tot
real(kind=dp), dimension(3) :: ene_coh, ene_tot

character(LEN=PATHLEN)          :: f_name
character(LEN=VARLEN)           :: f_comm
character(LEN=14), dimension(3) :: time_form 

!-------------------------------------------------------------------!

call initialize()

if (.not. main_isinit) then
   if (pMast) call mcio_printerr('main initialize','error occurred')
   main_iserr = .true.
   call finalize()
endif

if (pMast) then
   call mcio_printhead('OPTIMIZATION START')
   call mcio_printhead('OPTIMIZATION START',u=u)
endif

loop_opt : do iOptSite = 1, vinp%nOptSite
   if (pMast) then
      time_tot(2) = 0d0
      if (vinp%selective_site) rSite = floor(vlat%nSite*grnd()) + 1
      if (vinp%selective_unit) then
         do
            rUnit = floor(vinp%nUnit*grnd()) + 1
            if (vcoo%isMove(rUnit,mod(rSite,2))) exit
         enddo
      endif
      call main_printinfo_init()
   endif
   loop_mig : do iOptMig = 1, vinp%nOptMig
      loop_try : do iOptTry = 1, vinp%nOptTry
         if (pMast) then
            if (vinp%selective_site) then
               if (vinp%selective_unit) then
                  call opt_ranwalk(vinp%rate_dis, vinp%rate_ang, rSite, rUnit)
               else
                  call opt_ranwalk(vinp%rate_dis, vinp%rate_ang, rSite)
               endif
               call opt_coo_gather(vcoo, vopt, rSite)
            else
               do iSite = 1, vlat%nSite
                  call opt_ranwalk(vinp%rate_dis, vinp%rate_ang, iSite)
                  call opt_coo_gather(vcoo, vopt, iSite)
               enddo
            endif
         endif
         call pp_bcast(vopt%coo, 3, vopt%nAtom)
         ene_iserr = calc_energy(vcoo%nSpec, vopt%nAtom, vopt%LattVec, vcoo%NameAtom, &
                                 vopt%TypeAtoms, vinp%ispbc, vopt%coo, ene_coh(3), ene_tot(3))
         if (ene_iserr) then
            if (pMast) write (f_comm,'("iOpt iMig iTry - ",3I7)') iOptSite, iOptMig, iOptTry
            if (pMast) call mcio_printerr('main - energy calculation',f_comm)
            main_iserr = .true.
            call finalize()
         endif
         if (pMast) accept_ranwalk = isaccept_montecarlo(ene_tot(2), ene_tot(3))
         call pp_bcast(accept_ranwalk)
         if (accept_ranwalk) exit loop_try
      enddo loop_try
      if (pMast) then
         if (accept_ranwalk) then
            if (vinp%selective_site) then
               call opt_coo_update(rSite)
            else
               do iSite = 1, vlat%nSite
                  call opt_coo_update(iSite)
               enddo
            endif
         else
            if (vinp%selective_site) then
               call opt_coo_rollback(rSite)
               call opt_coo_gather(vcoo, vopt, rSite)
            else
               do iSite = 1, vlat%nSite
                  call opt_coo_rollback(iSite)
                  call opt_coo_gather(vcoo, vopt, iSite)
               enddo
            endif
         endif
         call main_calctime()
         call main_printinfo_mig(accept_ranwalk)
      endif
      if (.not. accept_ranwalk) exit loop_mig
   enddo loop_mig
   if (pMast) then
      call main_printinfo_site(accept_ranwalk)
      if (.not. (accept_ranwalk .or. (iOptMig .ne. 1))) cycle loop_opt
      ene_coh(1) = ene_coh(2)
      ene_tot(1) = ene_tot(2)
      write (f_name,'(A,"_OUT")') trim(vinp%fname_coo_out)
      write (f_comm,'("L_SITE=",I6.6," E_TOTAL=",ES16.9E2," eV")') iOptSite, ene_tot(1)
      call coo_writeposcar(f_name, vcoo%nSpec, vopt%nAtom, f_comm, vopt%LattVec, &
                           vcoo%NameAtom, vopt%nAtomSpec, .true., vopt%coo)
      if (vinp%isTrace) then
         write (f_name,'("trjs/",A,"_TRJ_",I6.6)') trim(vinp%fname_coo_out), iOptSite
         call coo_writeposcar(f_name, vcoo%nSpec, vopt%nAtom, f_comm, vopt%LattVec, &
                              vcoo%NameAtom, vopt%nAtomSpec, .true., vopt%coo)
      endif
   endif
enddo loop_opt

call finalize()

contains

!-------------------------------------------------------------------!

 subroutine initialize()
   implicit none

   if (main_isinit) return

   call pp_init()

   if (pMast) then      ! initialize program timer
      time_tot(:) = 0d0
      call system_clock(t_start, t_rate, t_max)
   endif

   vinp = inp_init('OPT')
   if (.not. vinp%pp) call mcio_printerr('main - input bcast','input is not bcasted')
   if (vinp%iserr) return

   if (      vinp%lseed) call pp_mtserial(vinp%seed)
   if (.not. vinp%lseed) call pp_mtserial()

   vlat = lat_init(vinp)
   if (.not. vlat%pp) call mcio_printerr('main - lattice bcast','lattice is not bcasted')
   if (vlat%iserr) return

   vcoo = coo_init(vinp)
   if (.not. vcoo%pp) call mcio_printerr('main - coord bcast','coord is not bcasted')
   if (vcoo%iserr) return

   vopt = opt_init(vinp, vcoo, vlat)
   if (.not. vopt%pp) call mcio_printerr('main - optimizer bcast','optimizer is not bcasted')
   if (vopt%iserr) return

   if (.not. (vinp%pp .and. vcoo%pp .and. vlat%pp .and. vopt%pp)) return

   ene_iserr = ene_init(vinp, vcoo)
   ene_iserr = calc_energy(vcoo%nSpec, vopt%nAtom, vopt%LattVec, vcoo%NameAtom, &
                           vopt%TypeAtoms, vinp%ispbc, vopt%coo, ene_coh(1), ene_tot(1))
   if (ene_iserr) return

   ene_coh(2) = ene_coh(1)
   ene_tot(2) = ene_tot(1)

   if (pMast) then
      call main_calctime()
      write (*,'(2x,"Elapsed time for initialization : ",A14)') time_form(3)
      u = mcio_getunit()
      open (u, name='opt.log', action='write')
   endif

   main_iserr  = .false.
   main_isinit = .true.

 end subroutine initialize

!-------------------------------------------------------------------!

 subroutine finalize()
   implicit none

   if (.not. vinp%iserr) call inp_fin(vinp)
   if (.not. vcoo%iserr) call coo_fin(vcoo)
   if (.not. vlat%iserr) call lat_fin(vlat)
   if (.not. vopt%iserr) call opt_fin(vopt)
   if (.not. ene_iserr)  call ene_fin(ene_iserr)

   if (pMast) then
      if (main_iserr) then
         call mcio_printhead('ERROR - OPTIMIZATION TERMINATED - ERROR')
      else
         write (*,'("2x,Elapsed Time : ",A13)') time_form(3)
         call mcio_printhead('OPTIMIZATION FINISHED')
         call mcio_printhead('OPTIMIZATION FINISHED',u=u)
      endif
      close (u)
   endif
   call pp_fin()
   stop

 end subroutine finalize

!-------------------------------------------------------------------!

 function isaccept_montecarlo(ene_1, ene_2) result(is_accept)
   implicit none
   real(kind=dp), intent(in) :: ene_1 ! old energy
   real(kind=dp), intent(in) :: ene_2 ! new energy

   logical       :: is_accept
   real(kind=dp) :: prob

   is_accept = .true.
   prob = exp((ene_1 - ene_2)/(kb_ev*vinp%temperature))
   if (prob .ge. 1d0) return 
   if (prob .ge. grnd()) return

   is_accept = .false.

 end function isaccept_montecarlo
 
!-------------------------------------------------------------------!

 subroutine main_printinfo_init()
   implicit none

   if (vinp%selective_site) then
      if (vinp%selective_unit) then
         write (*,'(2x,"iL_SITE= ",I7," (OF ",I7," )  SITE NUMBER : ",I5,"  UNIT NUMBER : ",I5)') &
                iOptSite, vinp%nOptSite, rSite, rUnit
      else
         write (*,'(2x,"iL_SITE= ",I7," (OF ",I7," )  SITE NUMBER : ",I5)') &
                iOptSite, vinp%nOptSite, rSite
      endif
   else
      write (*,'(2x,"iL_SITE= ",I7," (OF ",I7," )  SITE NUMBER : ALL SITES")') &
                iOptSite, vinp%nOptSite
   endif

 end subroutine main_printinfo_init

!-------------------------------------------------------------------!

 subroutine main_printinfo_mig(is_accept)
   implicit none
   logical, intent(in) :: is_accept

   if (is_accept) then
      write (*,'(5x,"iL_MIG= ",I5,"  iL_TRY= ",I5,"  dE= ",ES13.6E2, " eV  dT= ",A14)') &
                iOptMig, iOptTry, ene_tot(3)-ene_tot(2), adjustl(time_form(1))
      ene_tot(2) = ene_tot(3)
      ene_coh(2) = ene_coh(3)
   else
      write (*,'(5x,"iL_MIG= ",I5,"  iL_TRY= ",I5,"  dE= None (mig. fail)  dT= ",A14)') &
                iOptMig, iOptTry-1, adjustl(time_form(1))
   endif

 end subroutine main_printinfo_mig

!-------------------------------------------------------------------!

 subroutine main_printinfo_site(is_accept)
   implicit none
   logical, intent(in) :: is_accept

  if (is_accept .or. (iOptMig .ne. 1)) then
   write (*,'(3x,"E_COHENSIVE  = ",F25.12," eV")') ene_coh(2)
   write (*,'(3x,"E_TOTAL      = ",F25.12," eV")') ene_tot(2)
   write (*,'(3x,"dE_SITE      = ",F25.12," eV")') ene_tot(2) - ene_tot(1)
  endif
   write (*,'(3x,"Elapsed Time =   ",A14)') time_form(2)
  if (is_accept .or. (iOptMig .ne. 1)) then
   write (*,'(2x,"iLS= ",I6,"  ET= ",ES13.6E2,"  EC= ",ES13.6E2,"  dE= ",ES10.3E2,"  T= ",A14)') &
        iOptSite, ene_tot(2), ene_coh(2), ene_tot(2)-ene_tot(1), adjustl(time_form(3))
  else
   write (*,'(2x,"iLS= ",I6,"  ET= ",A13,     "  EC= ",A13,     "  dT= ",A10,     "  T= ",A14)') &
        iOptSite, "None(fail)", "None(fail)", "None(fail)", adjustl(time_form(3))
  endif
   write (*,*)

 end subroutine main_printinfo_site

!-------------------------------------------------------------------!

 subroutine main_calctime()
   implicit none

   integer :: t_h, t_m, t_s, t_ms
   integer :: ii

   call system_clock(t_end, t_rate, t_max)
   time_tot(1) = real(t_end - t_start)/real(t_rate)
   if (time_tot(1) .lt. 0d0) time_tot(1) = time_tot(1) + real(t_max)/real(t_rate)
   time_tot(2) = time_tot(2) + time_tot(1)
   time_tot(3) = time_tot(3) + time_tot(1)

   do ii = 1, 3
      t_h  = floor(time_tot(ii)/3600)
      t_m  = floor(time_tot(ii)/60 - t_h*60)
      t_s  = floor(time_tot(ii) - t_m*60 - t_h*3600)
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

end program main_opt
