module energy

use constants,  only :  dp,      &
                        VARLEN,  &
                        PATHLEN, &
                        CHKERR,  &
                        CHKOK

use util,       only :  cart_to_dir

use mcio,       only :  mcio_printhead, &
                        mcio_printerr

use parallel,   only :  root,           &
                        nrnk,           &
                        nprc,           &
                        pMast,          &
                        pp_fin,         &
                        pp_bcast,       &
                        pp_barrier,     &
                        pp_send,        &
                        pp_recv

use input,      only :  Var_Input,      &
                        inp_read_network

use coordinate, only :  Var_Coord

use aenet,      only :  aenet_init, aenet_final,        &
                        aenet_load_potential,           &
                        aenet_print_info,               &
                        aenet_Rc_min, aenet_Rc_max,     &
                        aenet_atomic_energy,            &
                        aenet_atomic_energy_and_forces, &
                        aenet_free_atom_energy,         &
                        aenet_nnb_max

use lclist,     only :  lcl_init,        &
                        lcl_final,       &
                        lcl_nmax_nbdist, &
                        lcl_nbdist_cart 

implicit none

character(LEN=VARLEN), private :: method
logical,               private :: ene_isinit   = .false.
logical,               private :: aenet_isinit = .false.
logical,               private :: amber_isinit = .false.

contains 

!-------------------------------------------------------------------!

 function ene_init(vinp, vcoo) result (iserr)
   implicit none
   type(Var_Input), intent(in)  :: vinp
   type(Var_Coord), intent(in)  :: vcoo

   integer :: iSpec, stat
   logical :: iserr
   character(LEN=PATHLEN), dimension(vcoo%nSpec) :: aenet_networks

   iserr  = .false.
   method = trim(adjustl(vinp%method))

   if (pMast) then
      call mcio_printhead('Energy calculation')
      write (*,'(1x,A,A)') ' Calculation method              : ', method
   endif

   select case (method)
   case ('AENET')
      if (pMast) then
         call mcio_printhead('AENET Initialization')
         call inp_read_network(vinp%fname_in, vcoo%nSpec, vcoo%NameAtom, aenet_networks, iserr)
      endif
      call pp_bcast(iserr)
      if (iserr) then
         if (pMast) call mcio_printerr('energy initalize','fail to read aenet network files')
         return
      endif
      call pp_bcast(aenet_networks, vcoo%nSpec)

      call aenet_init(vcoo%NameAtom, stat)
      if (stat .ne. CHKOK) then
         if (pMast) call mcio_printerr('aenet initialize','fail to initialize aenet')
         iserr = .true.
         return
      endif
      loop_spec : do iSpec = 1, vcoo%nSpec
         call aenet_load_potential(iSpec, aenet_networks(iSpec), stat)
         if (stat .ne. CHKOK) then
            if (pMast) call mcio_printerr('aenet initialize - network loading',vcoo%NameAtom(iSpec))
            iserr = .true.
            return
         endif
      enddo loop_spec
      if (pMast) then
         call aenet_print_info()
         call mcio_printhead('AENET initialization finished')
      endif
      aenet_isinit = .true.
   end select

   call pp_barrier()
   ene_isinit = .true.

 end function ene_init

!-------------------------------------------------------------------!

 function calc_energy(nSpec, nAtom, LattVec, NameAtom, TypeAtoms, &
                      ispbc, coo_in, E_coh, E_tot, force) result (iserr)
   implicit none
   integer,                                        intent(in)  :: nSpec
   integer,                                        intent(in)  :: nAtom
   real(kind=dp),    dimension(3,3),               intent(in)  :: LattVec
   character(LEN=*), dimension(nSpec),             intent(in)  :: NameAtom
   integer,          dimension(nAtom),             intent(in)  :: TypeAtoms
   logical,                                        intent(in)  :: ispbc
   real(kind=dp),    dimension(3,nAtom),           intent(in)  :: coo_in
   real(kind=dp),                                  intent(out) :: E_coh
   real(kind=dp),                                  intent(out) :: E_tot
   real(kind=dp),    dimension(3,nAtom), optional, intent(out) :: force

   logical :: iserr
   logical :: lforce = .false.

   iserr = .false.

   if (present(force)) lforce = .true.

   if (.not. ene_isinit) then
      if (pMast) call mcio_printerr('energy initialize', method)
      call ene_fin(iserr)
      return
   endif

   E_coh = 0d0
   E_tot = 0d0

   select case(method)
   case ('AENET')
      if (.not. aenet_isinit) then
         if (pMast) call mcio_printerr('energy initialize','AENET is not initialized')
         iserr = .true.
         return
      endif
      if (lforce) then
         call get_energy_aenet(nAtom, LattVec, TypeAtoms, ispbc, coo_in, E_coh, E_tot, force)
      else
         call get_energy_aenet(nAtom, LattVec, TypeAtoms, ispbc, coo_in, E_coh, E_tot)
      endif
   end select

 end function calc_energy

!-------------------------------------------------------------------!

 subroutine ene_fin(iserr) 
   implicit none
   logical,          intent(inout) :: iserr

   integer :: stat

   if (.not. ene_isinit) then
      if (pMast) call mcio_printerr('energy finalize','energy is not initialized')
      iserr = .true.
      return
   endif

   if (aenet_isinit) call aenet_final(stat)


 end subroutine ene_fin

!-------------------------------------------------------------------!
!       get_energy_aenet
!       come from 'predict.F90' of aenet package
!-------------------------------------------------------------------!

 subroutine get_energy_aenet(nAtom, LattVec, TypeAtoms, ispbc, coo_in, &
                             E_coh, E_tot, force)

   implicit none
   integer,                                     intent(in)    :: nAtom
   real(kind=dp), dimension(3,3),               intent(in)    :: LattVec
   integer,       dimension(nAtom),             intent(in)    :: TypeAtoms
   logical,                                     intent(in)    :: ispbc
   real(kind=dp), dimension(3,nAtom),           intent(in)    :: coo_in
   real(kind=dp),                               intent(inout) :: E_coh
   real(kind=dp),                               intent(inout) :: E_tot
   real(kind=dp), dimension(3,nAtom), optional, intent(out)   :: force

   integer                                   :: nnb
   integer,       dimension(aenet_nnb_max)   :: nblist
   real(kind=dp), dimension(3,aenet_nnb_max) :: nbcoo
   real(kind=dp), dimension(aenet_nnb_max)   :: nbdist
   integer,       dimension(aenet_nnb_max)   :: nbtype

   integer                            :: iAtom, iprc, stat
   logical                            :: lforce = .false.
   real(kind=dp)                      :: E_i, E_tmp_tot, E_tmp_coh
   real(kind=dp), dimension(3,nAtom)  :: coo_frac

!write (*,*) 'enter'
   loop_atom_1 : do iAtom = 1, nAtom
      if (mod(iAtom-1,nprc) .ne. nrnk) cycle loop_atom_1 ! parallel run
      coo_frac(1:3,iAtom) = coo_in(1:3,iAtom)
      call cart_to_dir(LattVec, coo_frac(1:3,iAtom))
!write (*,*) nrnk, iAtom, 'conv'
   enddo loop_atom_1

   loop_atom_2 : do iAtom = 1, nAtom
      if (nprc .eq. 1) exit loop_atom_2
      if (mod(iAtom-1,nprc) .ne. nrnk) cycle loop_atom_2 ! parallel run
      if (.not. pMast) then
         call pp_send(coo_frac(1:3,iAtom), root, 3)
      else
         loop_prc_1 : do iprc = 1, nprc-1
            if (iAtom+iprc .gt. nAtom) exit loop_prc_1
            call pp_recv(coo_frac(1:3,iAtom+iprc), iprc, 3)
         enddo loop_prc_1
      endif
   enddo loop_atom_2

   call pp_bcast(coo_frac, 3, nAtom)

   if (present(force)) lforce = .true.
!write (*,*) nrnk, 'ene calc start'
!if (pMast) then
!do iAtom = 1, nAtom
!   write (*,'(2I4,6F12.5)') nrnk, TypeAtoms(iAtom), coo_frac(:,iAtom), coo_in(:,iAtom)
!enddo
!endif
!stop
   call lcl_init(aenet_Rc_min, aenet_Rc_max, LattVec, nAtom, TypeAtoms, &
                 coo_frac, ispbc)
!write (*,*) nrnk,'lcl init'

   loop_atom_3 : do iAtom = 1, nAtom
      if (mod(iAtom-1,nprc) .ne. nrnk) cycle loop_atom_3 ! parallel run
      nnb = aenet_nnb_max
!write (*,*) nrnk, iAtom, 'calc start', nnb
      call lcl_nbdist_cart(iAtom, nnb, nbcoo, nbdist, aenet_Rc_max, &
                           nblist=nblist, nbtype=nbtype)
!write (*,*) nrnk, iAtom, 'lcl load', nnb
      if (lforce) then
         call aenet_atomic_energy_and_forces(coo_in(1:3,iAtom),    &
                        TypeAtoms(iAtom), iAtom, nnb, nbcoo, nbtype, &
                        nblist, nAtom, E_i, force, stat)
      else
         call aenet_atomic_energy(coo_in(1:3,iAtom), TypeAtoms(iAtom), &
                        nnb, nbcoo, nbtype, E_i, stat)
      endif
      E_tot = E_tot + E_i
      E_coh = E_coh + E_i - aenet_free_atom_energy(TypeAtoms(iAtom))
!write (*,*) nrnk, iAtom, nAtom, E_tot, E_coh
   enddo loop_atom_3

   call lcl_final()
   
!write (*,*) nrnk, 'calc finished'
   if (.not. pMast) then
      call pp_send(E_tot, root)
      call pp_send(E_coh, root) 
   elseif (nprc .ne. 1) then
      loop_prc_2 : do iprc = 1, nprc-1
         call pp_recv(E_tmp_tot, iprc)
         call pp_recv(E_tmp_coh, iprc)
         E_tot = E_tot + E_tmp_tot
         E_coh = E_coh + E_tmp_coh
      enddo loop_prc_2
   endif
!write (*,*) nrnk, 'comm fin'
   if (lforce) then
      loop_atom_4 : do iAtom = 1, nAtom
         if (nprc .eq. 1) exit loop_atom_4
         if (mod(iAtom-1,nprc) .ne. nrnk) cycle loop_atom_4 ! parallel run
         if (.not. pMast) then
            call pp_send(force(1:3,iAtom), root, 3)
         else
            loop_prc_3 : do iprc = 1, nprc-1
               if (iAtom+iprc .gt. nAtom) exit loop_prc_3
               call pp_recv(force(1:3,iAtom+iprc), iprc, 3)
            enddo loop_prc_3
         endif
      enddo loop_atom_4
   endif

 end subroutine get_energy_aenet

!-------------------------------------------------------------------!

end module energy
