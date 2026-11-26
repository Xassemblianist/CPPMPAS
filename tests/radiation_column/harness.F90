! SPDX-License-Identifier: GPL-3.0-or-later
! Column driver of MPAS's own cloudiness, RRTMG radiation and ozone climatology code (build_ref.sh): reads the pool
! arrays and the mpas_atmphys_vars arrays MPAS_to_physics fills (names starting with "@") from in.bin (records,
! tests/radiation_column/columns.py), sets the module variables physics_run_init sets (one rank, one thread:
! its..ite = 1..nCellsSolve, kts..kte = 1..nVertLevels), runs the initialisation of mpas_atmphys_init.F:397-449
! (init_o3climatology, init_radiation_sw, init_radiation_lw, reading OZONE_*.TBL and RRTMG_*_DATA.DBL from the
! working directory), then per step what physics_timetracker (mpas_atmphys_manager.F:199-212, 302-312) and
! physics_driver (mpas_atmphys_driver.F:202-246) do for radiation: the time variables, o3climatology_from_MPAS,
! allocate / driver / deallocate of the cloudiness, SW and LW drivers, and writes the snapshots of columns.py to
! out.bin. Compiled with -fconvert=big-endian as MPAS (the RRTMG_*_DATA.DBL reads); in.bin / out.bin are opened
! little-endian.
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
  subroutine marker(u, name, step)
    integer, intent(in) :: u, step
    character(len=*), intent(in) :: name
    character(len=32) :: nm
    nm = name
    write(u) nm, 1_4, 1_4, 1_4, 1_4, 1_4, 1_4, int(step, 4)
  end subroutine
  subroutine dump_all(u)  ! every pool array, in the order of in.bin (not hdr, dtl, time and the "@" arrays)
    integer, intent(in) :: u
    integer :: r
    do r = 1, nrec
      if (trim(recs(r)%name) == 'hdr' .or. trim(recs(r)%name) == 'dtl' .or. trim(recs(r)%name) == 'time') cycle
      if (recs(r)%name(1:1) == '@') cycle
      if (recs(r)%kind == 0) then
        write(u) recs(r)%name, recs(r)%kind, recs(r)%nd, recs(r)%dims, recs(r)%r
      else
        write(u) recs(r)%name, recs(r)%kind, recs(r)%nd, recs(r)%dims, recs(r)%iv
      end if
    end do
  end subroutine
end module recio

program radiation_column_ref
  use mpas_kind_types
  use mpas_pool_routines
  use mpas_derived_types
  use mpas_atmphys_vars
  use mpas_atmphys_manager, only: gmt, curr_julday, julday, year
  use mpas_atmphys_o3climatology, only: init_o3climatology, o3climatology_from_MPAS
  use mpas_atmphys_driver_cloudiness, only: allocate_cloudiness, deallocate_cloudiness, driver_cloudiness
  use mpas_atmphys_driver_radiation_sw, only: allocate_radiation_sw, deallocate_radiation_sw, driver_radiation_sw, &
                                              init_radiation_sw
  use mpas_atmphys_driver_radiation_lw, only: allocate_radiation_lw, deallocate_radiation_lw, driver_radiation_lw, &
                                              init_radiation_lw
  use recio
  implicit none
  type(mpas_pool_type) :: configs, mesh, state, diag_physics, atm_input, sfc_input, tend_physics
  type(dm_info) :: dminfo
  character(len=512) :: fin, fout
  integer :: hdr(10), nC, nV, nsteps, step
  integer, target :: d_nCells, d_nCellsSolve, d_nVertLevels, d_nOzn, d_nMonths
  real(kind=RKIND), target :: len_disp
  logical, target :: o3clim_cfg, microp_re = .false.
  character(len=StrKIND), target :: microp = 'mp_wsm6', radt_lw = 'rrtmg_lw', radt_sw = 'rrtmg_sw', &
    radt_cld = 'cld_fraction', overlap, dcorrlen
  real(kind=RKIND), pointer :: tm(:,:)
  real(kind=RKIND) :: xtime_s

  call get_command_argument(1, fin)
  call get_command_argument(2, fout)
  open(10, file=trim(fin), access='stream', form='unformatted', status='old', convert='little_endian')
  call read_all(10)
  close(10)
  hdr = recs(rec_index('hdr'))%iv(1:10)
  nC = hdr(1); nV = hdr(2); nsteps = hdr(3)
  select case (hdr(6))
  case (0); overlap = 'none'
  case (1); overlap = 'random'
  case (2); overlap = 'maximum_random'
  case (3); overlap = 'maximum'
  case (4); overlap = 'exponential'
  case (5); overlap = 'exponential_random'
  end select
  dcorrlen = merge('latitude_varying', 'constant        ', hdr(7) == 1)
  o3clim_cfg = hdr(8) == 1
  len_disp = recs(rec_index('dtl'))%r(3)
  tm(1:4, 1:nsteps) => recs(rec_index('time'))%r

  d_nCells = nC + 1; d_nCellsSolve = nC; d_nVertLevels = nV; d_nOzn = hdr(4); d_nMonths = hdr(5)
  call reg_i0('nCells', d_nCells); call reg_i0('nCellsSolve', d_nCellsSolve)
  call reg_i0('nVertLevels', d_nVertLevels); call reg_i0('nOznLevels', d_nOzn); call reg_i0('nMonths', d_nMonths)
  call reg_cc('config_microp_scheme', microp); call reg_cc('config_radt_lw_scheme', radt_lw)
  call reg_cc('config_radt_sw_scheme', radt_sw); call reg_cc('config_radt_cld_scheme', radt_cld)
  call reg_cc('config_radt_cld_overlap', overlap); call reg_cc('config_radt_cld_dcorrlen', dcorrlen)
  call reg_cl('config_o3climatology', o3clim_cfg); call reg_cl('config_microp_re', microp_re)
  call reg_cr('config_len_disp', len_disp)

  ! physics_run_init (mpas_atmphys_manager.F:466-490, 661-669, 737, 747-760)
  ims = 1; ime = nC; jms = 1; jme = 1; kms = 1; kme = nV + 1
  ids = ims; ide = ime + 1; jds = jms; jde = jme + 1; kds = kms; kde = kme
  its = ims; ite = ime; jts = jms; jte = jme; kts = kms; kte = kme - 1
  dt_radtsw = recs(rec_index('dtl'))%r(1)
  dt_radtlw = recs(rec_index('dtl'))%r(2)
  num_months = d_nMonths
  num_oznlevels = d_nOzn
  has_reqc = 0; has_reqi = 0; has_reqs = 0
  f_qc = .true.; f_qr = .true.; f_qi = .true.; f_qs = .true.; f_qg = .true.

  ! the arrays of allocate_forall_physics (mpas_atmphys_interface.F) the radiation code reads, as MPAS_to_physics
  ! leaves them
  allocate(pres_hyd_p(ims:ime,kms:kme,jms:jme), pres2_hyd_p(ims:ime,kms:kme,jms:jme), pi_p(ims:ime,kms:kme,jms:jme))
  allocate(t_p(ims:ime,kms:kme,jms:jme), t2_p(ims:ime,kms:kme,jms:jme), dz_p(ims:ime,kms:kme,jms:jme))
  allocate(pres_p(ims:ime,kms:kme,jms:jme), qv_p(ims:ime,kms:kme,jms:jme), qc_p(ims:ime,kms:kme,jms:jme))
  allocate(qi_p(ims:ime,kms:kme,jms:jme), qs_p(ims:ime,kms:kme,jms:jme))
  pres_hyd_p = 0; pres2_hyd_p = 0; pi_p = 0; t_p = 0; t2_p = 0; dz_p = 0; pres_p = 0
  qv_p = 0; qc_p = 0; qi_p = 0; qs_p = 0
  call setv('@pres_hyd_p', pres_hyd_p, nV); call setv('@pres2_hyd_p', pres2_hyd_p, nV + 1)
  call setv('@pi_p', pi_p, nV); call setv('@t_p', t_p, nV); call setv('@t2_p', t2_p, nV + 1)
  call setv('@dz_p', dz_p, nV); call setv('@pres_p', pres_p, nV)
  call setv('@qv_p', qv_p, nV); call setv('@qc_p', qc_p, nV); call setv('@qi_p', qi_p, nV); call setv('@qs_p', qs_p, nV)

  ! physics_init (mpas_atmphys_init.F:397, 428-449)
  if (o3clim_cfg) call init_o3climatology(mesh, atm_input)
  call init_radiation_sw(dminfo, configs)
  call init_radiation_lw(dminfo, configs)

  open(11, file=trim(fout), access='stream', form='unformatted', status='replace', convert='little_endian')
  call marker(11, '#I', 0)
  call dump_all(11)
  do step = 1, nsteps
    ! physics_timetracker (mpas_atmphys_manager.F:205-212, 302-312); gmt is the start time's (physics_run_init)
    xtime_s = tm(1, step)
    julday = nint(tm(2, step)); curr_julday = tm(3, step); gmt = tm(4, step); year = 2026
    if (o3clim_cfg) call o3climatology_from_MPAS(curr_julday, mesh, atm_input, diag_physics)
    call allocate_cloudiness
    call driver_cloudiness(configs, mesh, diag_physics, sfc_input, 1, nC)
    call allocate_radiation_sw(configs, xtime_s)
    call driver_radiation_sw(step, configs, mesh, state, 1, diag_physics, atm_input, sfc_input, tend_physics, &
                             xtime_s, 1, nC)
    call allocate_radiation_lw(configs, xtime_s)
    call driver_radiation_lw(xtime_s, configs, mesh, state, 1, diag_physics, atm_input, sfc_input, tend_physics, &
                             1, nC)
    call deallocate_cloudiness
    call deallocate_radiation_sw(configs)
    call deallocate_radiation_lw(configs)
    call marker(11, '#S', step)
    call dump_all(11)
  end do
  close(11)
contains
  subroutine setv(name, a, nk)
    character(len=*), intent(in) :: name
    real(kind=RKIND), intent(inout) :: a(:,:,:)
    integer, intent(in) :: nk
    integer :: r
    r = rec_index(name)
    a(1:nC, 1:nk, 1) = reshape(recs(r)%r, [nC, nk])
  end subroutine
end program radiation_column_ref
