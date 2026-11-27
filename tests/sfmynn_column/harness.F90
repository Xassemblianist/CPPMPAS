! SPDX-License-Identifier: GPL-3.0-or-later
! Column driver of MPAS's own MYNN surface layer and physics interface (build_ref.sh): reads the pool arrays from
! in.bin (records, tests/sfmynn_column/columns.py), sets the module variables physics_run_init sets (one rank, one
! thread: its..ite = 1..nCellsSolve, kts..kte = 1..nVertLevels), calls init_sfclayer (sf_mynn_init) once, then per
! step what physics_driver (mpas_atmphys_driver.F:202-277) calls for the surface layer: allocate_forall_physics,
! MPAS_to_physics, allocate / driver / deallocate_sfclayer, and writes the snapshots of columns.py to out.bin.
! config_pbl_scheme is 'bl_ysu' here: it selects only what allocate_forall_physics / MPAS_to_physics add for bl_mynn
! (the number concentrations nc_p, ni_p, nifa_p, nwfa_p), which the surface layer does not read; physics_control's
! check of the pair (sf_mynn goes with bl_mynn) is not run. Compiled without -fconvert=big-endian (native-endian files).
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

program sfmynn_column_ref
  use mpas_kind_types
  use mpas_pool_routines
  use mpas_atm_dimensions, only: nVertLevels, num_scalars
  use mpas_atmphys_vars
  use mpas_atmphys_interface, only: allocate_forall_physics, deallocate_forall_physics, MPAS_to_physics
  use mpas_atmphys_driver_sfclayer, only: init_sfclayer, allocate_sfclayer, deallocate_sfclayer, driver_sfclayer
  use recio
  implicit none
  interface
    subroutine get_sfmynn_tables(t) bind(C, name='get_sfmynn_tables')
      real(8), intent(out) :: t(1001, 4)
    end subroutine
  end interface
  type(mpas_pool_type) :: configs, mesh, state, diag, diag_physics, sfc_input
  character(len=512) :: fin, fout
  integer :: hdr(8), nC, nCt, nV, nsteps, step
  integer, target :: d_nCells, d_nCellsSolve
  integer, target :: iqv = 1, iqc = 2, iqr = 3, iqi = 4, iqs = 5, iqg = 6, inone = -1
  real(kind=RKIND), target :: len_disp
  logical, target :: frac_seaice, do_restart, do_dacycling = .false.
  character(len=StrKIND), target :: microp = 'mp_wsm6', pbl = 'bl_ysu', sfclayer = 'sf_mynn'
  real(8) :: tables(1001, 4)

  call get_command_argument(1, fin)
  call get_command_argument(2, fout)
  open(10, file=trim(fin), access='stream', form='unformatted', status='old')
  call read_all(10)
  close(10)
  hdr = recs(rec_index('hdr'))%iv(1:8)
  nC = hdr(1); nCt = hdr(2); nV = hdr(3); nsteps = hdr(4)
  frac_seaice = hdr(5) == 1
  do_restart = hdr(6) == 1
  len_disp = recs(rec_index('dtl'))%r(2)

  d_nCells = nCt; d_nCellsSolve = nC
  call reg_i0('nCells', d_nCells); call reg_i0('nCellsSolve', d_nCellsSolve)
  call reg_i0('index_qv', iqv); call reg_i0('index_qc', iqc); call reg_i0('index_qr', iqr)
  call reg_i0('index_qi', iqi); call reg_i0('index_qs', iqs); call reg_i0('index_qg', iqg)
  call reg_i0('index_nc', inone); call reg_i0('index_ni', inone); call reg_i0('index_nifa', inone)
  call reg_i0('index_nwfa', inone)
  call reg_cc('config_microp_scheme', microp); call reg_cc('config_pbl_scheme', pbl)
  call reg_cc('config_sfclayer_scheme', sfclayer)
  call reg_cl('config_frac_seaice', frac_seaice); call reg_cl('config_do_restart', do_restart)
  call reg_cl('config_do_DAcycling', do_dacycling)
  call reg_cr('config_len_disp', len_disp)

  ! physics_run_init (mpas_atmphys_manager.F:540-544, 661-669, 740-744) and init_physics_flags (mpas_atmphys_init.F)
  nVertLevels = nV; num_scalars = 6
  ims = 1; ime = nC; jms = 1; jme = 1; kms = 1; kme = nV + 1
  ids = ims; ide = ime + 1; jds = jms; jde = jme + 1; kds = kms; kde = kme
  its = ims; ite = ime; jts = jms; jte = jme; kts = kms; kte = kme - 1
  xice_threshold = 0.5
  if (frac_seaice) xice_threshold = 0.02
  f_qc = .true.; f_qr = .true.; f_qi = .true.; f_qs = .true.; f_qg = .true.; f_qoz = .false.
  f_nc = .false.; f_ni = .false.; f_nifa = .false.; f_nwfa = .false.; f_nbca = .false.

  call init_sfclayer(configs)

  open(11, file=trim(fout), access='stream', form='unformatted', status='replace')
  call get_sfmynn_tables(tables)
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
    call marker(11, '#B', step)
    call dump_all(11)
    if (frac_seaice) then  ! the *_sea arrays deallocate_sfclayer keeps for driver_seaice
      call put2(11, '@chs_sea', chs_sea); call put2(11, '@chs2_sea', chs2_sea); call put2(11, '@cqs2_sea', cqs2_sea)
      call put2(11, '@cpm_sea', cpm_sea); call put2(11, '@hfx_sea', hfx_sea); call put2(11, '@qfx_sea', qfx_sea)
      call put2(11, '@qgh_sea', qgh_sea); call put2(11, '@qsfc_sea', qsfc_sea); call put2(11, '@lh_sea', lh_sea)
    end if
    call deallocate_sfclayer(configs)
    call deallocate_forall_physics(configs)
  end do
  close(11)
contains
  subroutine put2(u, name, a)
    integer, intent(in) :: u
    character(len=*), intent(in) :: name
    real(kind=RKIND), intent(in) :: a(:,:)
    call put_r(u, name, 1, [nC, 1, 1, 1], a(1:nC, 1))
  end subroutine
  subroutine dump_vars(u)  ! the mpas_atmphys_vars arrays MPAS_to_physics defines (1:nC, levels it sets)
    integer, intent(in) :: u
    call dumpv(u, '@qv_p', qv_p, 1, nV); call dumpv(u, '@qc_p', qc_p, 1, nV)
    call dumpv(u, '@u_p', u_p, 1, nV); call dumpv(u, '@v_p', v_p, 1, nV)
    call dumpv(u, '@rho_p', rho_p, 1, nV); call dumpv(u, '@th_p', th_p, 1, nV); call dumpv(u, '@t_p', t_p, 1, nV)
    call dumpv(u, '@pi_p', pi_p, 1, nV); call dumpv(u, '@dz_p', dz_p, 1, nV)
    call dumpv(u, '@pres_hyd_p', pres_hyd_p, 1, nV)
    call put_r(u, '@psfc_p', 1, [nC, 1, 1, 1], psfc_p(1:nC, 1))
  end subroutine
  subroutine dumpv(u, name, a, k0, k1)
    integer, intent(in) :: u, k0, k1
    character(len=*), intent(in) :: name
    real(kind=RKIND), intent(in) :: a(:,:,:)
    call put_r(u, name, 2, [nC, k1 - k0 + 1, 1, 1], reshape(a(1:nC, k0:k1, 1), [nC * (k1 - k0 + 1)]))
  end subroutine
end program sfmynn_column_ref
