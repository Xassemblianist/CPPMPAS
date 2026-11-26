! SPDX-License-Identifier: GPL-3.0-or-later
! Column driver of MPAS's own surface layer, YSU PBL, physics interface and physics tendencies (build_ref.sh): reads
! the pool arrays from in.bin (records, tests/pbl_column/columns.py), sets the module variables physics_run_init sets
! (one rank, one thread: its..ite = 1..nCellsSolve, kts..kte = 1..nVertLevels), calls per step what physics_driver
! (mpas_atmphys_driver.F:202-330) and the dynamics call for these schemes: allocate_forall_physics, MPAS_to_physics,
! allocate / driver / deallocate_sfclayer, allocate / driver / deallocate_pbl, then physics_get_tend with rk_step 1
! and 2, and writes the snapshots of columns.py to out.bin. Compiled without -fconvert=big-endian (native-endian files).
module recio
  use mpas_kind_types
  use mpas_pool_routines
  implicit none
  integer, parameter :: maxrec = 200
  type rec_t
    character(len=32) :: name
    integer(4) :: kind, nd, dims(4)
    real(kind=RKIND), pointer :: r(:) => null()
    integer, pointer :: iv(:) => null()
  end type
  type(rec_t), target :: recs(maxrec)
  integer :: nrec = 0
contains
  subroutine read_all(u)
    integer, intent(in) :: u
    character(len=32) :: nm
    integer(4) :: meta(6)
    integer :: ios, n
    real(kind=RKIND), pointer :: p1(:), p2(:,:), p3(:,:,:)
    integer, pointer :: q1(:), q2(:,:)
    do
      read(u, iostat=ios) nm
      if (ios /= 0) exit
      read(u) meta
      nrec = nrec + 1
      recs(nrec)%name = nm
      recs(nrec)%kind = meta(1)
      recs(nrec)%nd = meta(2)
      recs(nrec)%dims = meta(3:6)
      n = product(meta(3:2+meta(2)))
      if (meta(1) == 0) then
        allocate(recs(nrec)%r(n))
        read(u) recs(nrec)%r
        select case (meta(2))
        case (1)
          p1(1:meta(3)) => recs(nrec)%r
          call reg_r1(trim(nm), p1)
        case (2)
          p2(1:meta(3), 1:meta(4)) => recs(nrec)%r
          call reg_r2(trim(nm), p2)
        case (3)
          p3(1:meta(3), 1:meta(4), 1:meta(5)) => recs(nrec)%r
          call reg_r3(trim(nm), p3)
        end select
      else
        allocate(recs(nrec)%iv(n))
        read(u) recs(nrec)%iv
        select case (meta(2))
        case (1)
          q1(1:meta(3)) => recs(nrec)%iv
          call reg_i1(trim(nm), q1)
        case (2)
          q2(1:meta(3), 1:meta(4)) => recs(nrec)%iv
          call reg_i2(trim(nm), q2)
        end select
      end if
    end do
  end subroutine
  integer function rec_index(name)
    character(len=*), intent(in) :: name
    do rec_index = 1, nrec
      if (trim(recs(rec_index)%name) == name) return
    end do
    print *, 'no record ', name
    stop 1
  end function
  subroutine put_r(u, name, nd, dims, a)
    integer, intent(in) :: u, nd, dims(4)
    character(len=*), intent(in) :: name
    real(kind=RKIND), intent(in) :: a(:)
    character(len=32) :: nm
    nm = name
    write(u) nm, 0_4, int(nd, 4), int(dims, 4), a
  end subroutine
  subroutine marker(u, name, step)
    integer, intent(in) :: u, step
    character(len=*), intent(in) :: name
    character(len=32) :: nm
    nm = name
    write(u) nm, 1_4, 1_4, 1_4, 1_4, 1_4, 1_4, int(step, 4)
  end subroutine
  subroutine dump_all(u)  ! every pool array, in the order of in.bin
    integer, intent(in) :: u
    integer :: r
    do r = 1, nrec
      if (trim(recs(r)%name) == 'hdr' .or. trim(recs(r)%name) == 'dtl') cycle
      if (recs(r)%kind == 0) then
        write(u) recs(r)%name, recs(r)%kind, recs(r)%nd, recs(r)%dims, recs(r)%r
      else
        write(u) recs(r)%name, recs(r)%kind, recs(r)%nd, recs(r)%dims, recs(r)%iv
      end if
    end do
  end subroutine
end module recio

subroutine no_halo_exchange(domain, halo_group, withGPUAwareMPI, ierr)
  use mpas_derived_types, only : domain_type
  type (domain_type), intent(inout) :: domain
  character(len=*), intent(in) :: halo_group
  logical, intent(in), optional :: withGPUAwareMPI
  integer, intent(out), optional :: ierr
  ! one rank: the halo cells hold what in.bin gives them
end subroutine

program pbl_column_ref
  use mpas_kind_types
  use mpas_pool_routines
  use mpas_derived_types
  use mpas_atm_dimensions, only: nVertLevels, num_scalars
  use mpas_atmphys_vars
  use mpas_atmphys_interface, only: allocate_forall_physics, deallocate_forall_physics, MPAS_to_physics
  use mpas_atmphys_driver_sfclayer, only: init_sfclayer, allocate_sfclayer, deallocate_sfclayer, driver_sfclayer
  use mpas_atmphys_driver_pbl, only: init_pbl, allocate_pbl, deallocate_pbl, driver_pbl
  use mpas_atmphys_todynamics, only: physics_get_tend
  use recio
  implicit none
  interface
    subroutine get_sfclayrev_tables(t) bind(C, name='get_sfclayrev_tables')
      real(8), intent(out) :: t(1001, 4)
    end subroutine
    subroutine no_halo_exchange(domain, halo_group, withGPUAwareMPI, ierr)
      use mpas_derived_types, only : domain_type
      type (domain_type), intent(inout) :: domain
      character(len=*), intent(in) :: halo_group
      logical, intent(in), optional :: withGPUAwareMPI
      integer, intent(out), optional :: ierr
    end subroutine
  end interface
  type(mpas_pool_type) :: configs, mesh, state, diag, diag_physics, tend_physics, tend, sfc_input
  type(block_type) :: block
  type(domain_type), target :: domain
  character(len=512) :: fin, fout
  integer :: hdr(10), nC, nCt, nEs, nE, nV, nsteps, step, rk
  integer, target :: d_nCells, d_nEdges, d_nCellsSolve, d_nEdgesSolve
  integer, target :: iqv = 1, iqc = 2, iqr = 3, iqi = 4, iqs = 5, iqg = 6, inone = -1
  real(kind=RKIND), target :: len_disp
  logical, target :: frac_seaice, ysu_pblmix_cfg, do_restart = .false., do_dacycling = .false.
  character(len=StrKIND), target :: microp = 'mp_wsm6', pbl = 'bl_ysu', sfclayer = 'sf_monin_obukhov_rev', &
    convection, radt_lw, radt_sw
  real(kind=RKIND), pointer :: tend_ru(:,:), tend_rth(:,:), tend_rho(:,:)
  real(8) :: tables(1001, 4)

  call get_command_argument(1, fin)
  call get_command_argument(2, fout)
  open(10, file=trim(fin), access='stream', form='unformatted', status='old')
  call read_all(10)
  close(10)
  hdr = recs(rec_index('hdr'))%iv(1:10)
  nC = hdr(1); nCt = hdr(2); nEs = hdr(3); nE = hdr(4); nV = hdr(5); nsteps = hdr(6)
  frac_seaice = hdr(7) == 1
  ysu_pblmix_cfg = hdr(8) == 1
  select case (hdr(9))
  case (0); convection = 'off'
  case (1); convection = 'cu_ntiedtke'
  case (2); convection = 'cu_kain_fritsch'
  case (3); convection = 'cu_grell_freitas'
  end select
  radt_lw = 'off'; radt_sw = 'off'
  if (hdr(10) == 1) then
    radt_lw = 'rrtmg_lw'; radt_sw = 'rrtmg_sw'
  end if
  len_disp = recs(rec_index('dtl'))%r(2)

  d_nCells = nCt; d_nEdges = nE; d_nCellsSolve = nC; d_nEdgesSolve = nEs
  call reg_i0('nCells', d_nCells); call reg_i0('nEdges', d_nEdges)
  call reg_i0('nCellsSolve', d_nCellsSolve); call reg_i0('nEdgesSolve', d_nEdgesSolve)
  call reg_i0('index_qv', iqv); call reg_i0('index_qc', iqc); call reg_i0('index_qr', iqr)
  call reg_i0('index_qi', iqi); call reg_i0('index_qs', iqs); call reg_i0('index_qg', iqg)
  call reg_i0('index_nc', inone); call reg_i0('index_ni', inone); call reg_i0('index_nifa', inone)
  call reg_i0('index_nwfa', inone)
  call reg_cc('config_microp_scheme', microp); call reg_cc('config_pbl_scheme', pbl)
  call reg_cc('config_sfclayer_scheme', sfclayer); call reg_cc('config_convection_scheme', convection)
  call reg_cc('config_radt_lw_scheme', radt_lw); call reg_cc('config_radt_sw_scheme', radt_sw)
  call reg_cl('config_frac_seaice', frac_seaice); call reg_cl('config_ysu_pblmix', ysu_pblmix_cfg)
  call reg_cl('config_do_restart', do_restart); call reg_cl('config_do_DAcycling', do_dacycling)
  call reg_cr('config_len_disp', len_disp)
  call mpas_pool_get_array(diag, 'tend_ru_physics', tend_ru)
  call mpas_pool_get_array(diag, 'tend_rtheta_physics', tend_rth)
  call mpas_pool_get_array(diag, 'tend_rho_physics', tend_rho)
  block%domain => domain

  ! physics_run_init (mpas_atmphys_manager.F:540-544, 661-669, 740-744) and init_physics_flags (mpas_atmphys_init.F)
  nVertLevels = nV; num_scalars = 6
  ims = 1; ime = nC; jms = 1; jme = 1; kms = 1; kme = nV + 1
  ids = ims; ide = ime + 1; jds = jms; jde = jme + 1; kds = kms; kde = kme
  its = ims; ite = ime; jts = jms; jte = jme; kts = kms; kte = kme - 1
  dt_pbl = recs(rec_index('dtl'))%r(1)
  xice_threshold = 0.5
  if (frac_seaice) xice_threshold = 0.02
  f_qc = .true.; f_qr = .true.; f_qi = .true.; f_qs = .true.; f_qg = .true.; f_qoz = .false.
  f_nc = .false.; f_ni = .false.; f_nifa = .false.; f_nwfa = .false.; f_nbca = .false.

  call init_sfclayer(configs)
  call init_pbl(configs)

  open(11, file=trim(fout), access='stream', form='unformatted', status='replace')
  call get_sfclayrev_tables(tables)
  call put_r(11, 'psim_stab', 1, [1001, 1, 1, 1], tables(:, 1))
  call put_r(11, 'psim_unstab', 1, [1001, 1, 1, 1], tables(:, 2))
  call put_r(11, 'psih_stab', 1, [1001, 1, 1, 1], tables(:, 3))
  call put_r(11, 'psih_unstab', 1, [1001, 1, 1, 1], tables(:, 4))
  do step = 1, nsteps
    call allocate_forall_physics(configs)
    call MPAS_to_physics(configs, mesh, state, 1, diag, diag_physics, 1, nC)
    call marker(11, '#A', step)
    call dump_all(11)
    call dump_vars(11)
    call allocate_sfclayer(configs)
    call driver_sfclayer(step, configs, mesh, diag_physics, sfc_input, 1, nC)
    call deallocate_sfclayer(configs)
    call marker(11, '#B', step)
    call dump_all(11)
    call allocate_pbl(configs)
    call driver_pbl(step, configs, mesh, sfc_input, diag_physics, tend_physics, 1, nC)
    call deallocate_pbl(configs)
    call marker(11, '#C', step)
    call dump_all(11)
    do rk = 1, 2
      call physics_get_tend(block, mesh, state, diag, tend, tend_physics, configs, rk, 1, tend_ru, tend_rth, &
                            tend_rho, no_halo_exchange)
      call marker(11, merge('#D', '#E', rk == 1), step)
      call dump_all(11)
    end do
    call deallocate_forall_physics(configs)
  end do
  close(11)
contains
  subroutine dump_vars(u)  ! the mpas_atmphys_vars arrays MPAS_to_physics defines (1:nC, levels it sets)
    integer, intent(in) :: u
    call dumpv(u, '@qv_p', qv_p, 1, nV); call dumpv(u, '@qc_p', qc_p, 1, nV); call dumpv(u, '@qr_p', qr_p, 1, nV)
    call dumpv(u, '@qi_p', qi_p, 1, nV); call dumpv(u, '@qs_p', qs_p, 1, nV); call dumpv(u, '@qg_p', qg_p, 1, nV)
    call dumpv(u, '@u_p', u_p, 1, nV); call dumpv(u, '@v_p', v_p, 1, nV); call dumpv(u, '@zz_p', zz_p, 1, nV)
    call dumpv(u, '@rho_p', rho_p, 1, nV); call dumpv(u, '@th_p', th_p, 1, nV); call dumpv(u, '@t_p', t_p, 1, nV)
    call dumpv(u, '@pi_p', pi_p, 1, nV); call dumpv(u, '@pres_p', pres_p, 1, nV); call dumpv(u, '@zmid_p', zmid_p, 1, nV)
    call dumpv(u, '@dz_p', dz_p, 1, nV); call dumpv(u, '@znu_p', znu_p, 1, nV)
    call dumpv(u, '@w_p', w_p, 1, nV + 1); call dumpv(u, '@z_p', z_p, 1, nV + 1)
    call dumpv(u, '@fzm_p', fzm_p, 2, nV); call dumpv(u, '@fzp_p', fzp_p, 2, nV)
    call dumpv(u, '@t2_p', t2_p, 1, nV + 1); call dumpv(u, '@pres2_p', pres2_p, 1, nV + 1)
    call dumpv(u, '@pres2_hyd_p', pres2_hyd_p, 1, nV + 1); call dumpv(u, '@pres2_hydd_p', pres2_hydd_p, 1, nV + 1)
    call dumpv(u, '@pres_hyd_p', pres_hyd_p, 1, nV); call dumpv(u, '@pres_hydd_p', pres_hydd_p, 1, nV)
    call dumpv(u, '@znu_hyd_p', znu_hyd_p, 1, nV)
    call put_r(u, '@psfc_p', 1, [nC, 1, 1, 1], psfc_p(1:nC, 1))
    call put_r(u, '@psfc_hyd_p', 1, [nC, 1, 1, 1], psfc_hyd_p(1:nC, 1))
    call put_r(u, '@psfc_hydd_p', 1, [nC, 1, 1, 1], psfc_hydd_p(1:nC, 1))
  end subroutine
  subroutine dumpv(u, name, a, k0, k1)
    integer, intent(in) :: u, k0, k1
    character(len=*), intent(in) :: name
    real(kind=RKIND), intent(in) :: a(:,:,:)
    call put_r(u, name, 2, [nC, k1 - k0 + 1, 1, 1], reshape(a(1:nC, k0:k1, 1), [nC * (k1 - k0 + 1)]))
  end subroutine
end program pbl_column_ref
