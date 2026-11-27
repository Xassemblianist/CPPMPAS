! SPDX-License-Identifier: GPL-3.0-or-later
! Column driver of MPAS's own Thompson path (build_ref.sh): reads the pool arrays from in.bin, sets the module
! variables physics_run_init sets (one rank: its..ite = 1..nCells, kts..kte = 1..nVertLevels; dt_microp and n_microp
! as mpas_atmphys_manager.F:684-694 sets them for mp_thompson), calls init_microphysics (thompson_init with
! l_mp_tables = .true., reading the MP_THOMPSON_*_DATA.DBL files of the current directory, and
! init_thompson_clouddroplets_forMPAS), writes every module variable thompson_init sets (dump.c) to init.bin, then
! calls driver_microphysics nsteps times and writes the pool arrays it changed after every step to out.bin.
! This main program is compiled with -fconvert=big-endian like MPAS's (gfortran applies -fconvert from the main
! program: the tables are read as the reference reads them); in.bin and out.bin are opened little-endian.
! Layouts: tests/thompson_column/columns.py. A fourth argument 'l_diags' sets l_diags = .true. and appends
! refl10cm, refl10cm_max, refl10cm_1km, refl10cm_1km_max and relhum to every step's output.
program thompson_column_ref
  use iso_c_binding, only: c_char, c_null_char
  use mpas_kind_types
  use mpas_derived_types, only: dm_info
  use mpas_pool_routines
  use mpas_atmphys_vars
  use mpas_atmphys_driver_microphysics, only: driver_microphysics, init_microphysics
  implicit none
  interface
    subroutine dump_thompson_init(path) bind(C, name='dump_thompson_init')
      import :: c_char
      character(kind=c_char), intent(in) :: path(*)
    end subroutine
  end interface
  integer(4) :: hdr(5)
  integer :: nC, nV, nsteps, step, nscal
  integer, target :: mp_top, iqv = 1, iqc = 2, iqr = 3, iqi = 4, iqs = 5, iqg = 6, ini = 7, inr = 8, ncs, nvl
  real(kind=RKIND) :: dt
  real(kind=RKIND), allocatable, target :: zgrid(:,:), zz(:,:), rho_zz(:,:), theta_m(:,:), w(:,:), scalars(:,:,:), &
    exner(:,:), pressure_base(:,:), pressure_p(:,:), exner_base(:,:), rtheta_base(:,:), rtheta_p(:,:), &
    surface_pressure(:), rainprod(:,:), evapprod(:,:), re_cloud(:,:), re_ice(:,:), re_snow(:,:), rainnc(:), snownc(:), &
    graupelnc(:), rainncv(:), snowncv(:), graupelncv(:), sr(:), precipw(:), tend_sfc_pressure(:), rt_diabatic_tend(:,:), &
    dtheta_dt_mp(:,:), refl10cm(:,:), rthmpten(:,:), rqvmpten(:,:), rqcmpten(:,:), rqrmpten(:,:), rqimpten(:,:), &
    rqsmpten(:,:), rqgmpten(:,:), rnimpten(:,:), rnrmpten(:,:), refl10cm_max(:), refl10cm_1km(:), refl10cm_1km_max(:), &
    relhum(:,:), nt_c(:), mu_c(:), rlandmask(:)
  logical :: diags
  character(len=16) :: a4
  integer, allocatable, target :: i_rainnc(:), landmask(:)
  type(mpas_pool_type) :: configs, mesh, state, diag, diag_physics, tend_physics, tend, sfc_input
  type(dm_info) :: dminfo
  character(len=512) :: fin, fout, finit

  call get_command_argument(1, fin)
  call get_command_argument(2, fout)
  call get_command_argument(3, finit)
  call get_command_argument(4, a4)
  diags = trim(a4) == 'l_diags'
  open(10, file=trim(fin), access='stream', form='unformatted', status='old', convert='little_endian')
  read(10) hdr
  nC = hdr(1); nV = hdr(2); nsteps = hdr(3); mp_top = hdr(4); nscal = hdr(5)
  read(10) dt
  allocate(zgrid(nV+1,nC), zz(nV,nC), rho_zz(nV,nC), theta_m(nV,nC), w(nV+1,nC), scalars(nscal,nV,nC), &
    exner(nV,nC), pressure_base(nV,nC), pressure_p(nV,nC), exner_base(nV,nC), rtheta_base(nV,nC), rtheta_p(nV,nC), &
    surface_pressure(nC), rainprod(nV,nC), evapprod(nV,nC), re_cloud(nV,nC), re_ice(nV,nC), re_snow(nV,nC), &
    rainnc(nC), snownc(nC), graupelnc(nC), rainncv(nC), snowncv(nC), graupelncv(nC), sr(nC), precipw(nC), &
    tend_sfc_pressure(nC), rt_diabatic_tend(nV,nC), dtheta_dt_mp(nV,nC), refl10cm(nV,nC), rthmpten(nV,nC), &
    rqvmpten(nV,nC), rqcmpten(nV,nC), rqrmpten(nV,nC), rqimpten(nV,nC), rqsmpten(nV,nC), rqgmpten(nV,nC), &
    rnimpten(nV,nC), rnrmpten(nV,nC), i_rainnc(nC), refl10cm_max(nC), refl10cm_1km(nC), refl10cm_1km_max(nC), &
    relhum(nV,nC), nt_c(nC), mu_c(nC), rlandmask(nC), landmask(nC))
  read(10) zgrid, zz, rho_zz, theta_m, w, scalars, exner, pressure_base, pressure_p, exner_base, rtheta_base, &
    rtheta_p, surface_pressure, rainprod, evapprod, re_cloud, re_ice, re_snow, rainnc, snownc, graupelnc, rlandmask, &
    nt_c, mu_c
  close(10)
  landmask = nint(rlandmask)
  ! values the driver must overwrite before reading them
  rainncv = -999.; snowncv = -999.; graupelncv = -999.; sr = -999.; precipw = -999.; tend_sfc_pressure = -999.
  rt_diabatic_tend = -999.; dtheta_dt_mp = -999.; i_rainnc = 0
  rthmpten = -999.; rqvmpten = -999.; rqcmpten = -999.; rqrmpten = -999.; rqimpten = -999.; rqsmpten = -999.
  rqgmpten = -999.; rnimpten = -999.; rnrmpten = -999.
  ! diagnostics: overwritten by l_diags except the running maximum refl10cm_1km_max (its value since the last reset)
  refl10cm = -999.; refl10cm_max = -999.; refl10cm_1km = -999.; refl10cm_1km_max = 3.; relhum = -999.
  ! every third cell starts from a quiet NaN: MAX(refl10cm_1km_max(i), refl10cm_1km(i)) returns its first argument
  ! unless the second is greater (maxsd), so this shows the operand order
  refl10cm_1km_max(2::3) = transfer(int(z'7FF8000000000000', 8), 1.)

  ncs = nC; nvl = nV
  call reg_r2('zgrid', zgrid); call reg_r2('zz', zz); call reg_r2('rho_zz', rho_zz); call reg_r2('theta_m', theta_m)
  call reg_r2('w', w); call reg_r3('scalars', scalars); call reg_r2('exner', exner)
  call reg_r2('pressure_base', pressure_base); call reg_r2('pressure_p', pressure_p)
  call reg_r2('exner_base', exner_base); call reg_r2('rtheta_base', rtheta_base); call reg_r2('rtheta_p', rtheta_p)
  call reg_r1('surface_pressure', surface_pressure); call reg_r1('tend_sfc_pressure', tend_sfc_pressure)
  call reg_r2('rt_diabatic_tend', rt_diabatic_tend); call reg_r2('dtheta_dt_mp', dtheta_dt_mp)
  call reg_r2('rainprod', rainprod); call reg_r2('evapprod', evapprod); call reg_r2('re_cloud', re_cloud)
  call reg_r2('re_ice', re_ice); call reg_r2('re_snow', re_snow); call reg_r2('refl10cm', refl10cm)
  call reg_r1('rainnc', rainnc); call reg_r1('rainncv', rainncv); call reg_r1('snownc', snownc)
  call reg_r1('snowncv', snowncv); call reg_r1('graupelnc', graupelnc); call reg_r1('graupelncv', graupelncv)
  call reg_r1('sr', sr); call reg_r1('precipw', precipw); call reg_i1('i_rainnc', i_rainnc)
  call reg_r2('rthmpten', rthmpten); call reg_r2('rqvmpten', rqvmpten); call reg_r2('rqcmpten', rqcmpten)
  call reg_r2('rqrmpten', rqrmpten); call reg_r2('rqimpten', rqimpten); call reg_r2('rqsmpten', rqsmpten)
  call reg_r2('rqgmpten', rqgmpten); call reg_r2('rnimpten', rnimpten); call reg_r2('rnrmpten', rnrmpten)
  call reg_i0('mp_top_level', mp_top)
  call reg_r1('refl10cm_max', refl10cm_max); call reg_r1('refl10cm_1km', refl10cm_1km)
  call reg_r1('refl10cm_1km_max', refl10cm_1km_max); call reg_r2('relhum', relhum)
  call reg_r1('nt_c', nt_c); call reg_r1('mu_c', mu_c); call reg_i1('landmask', landmask)
  call reg_i0('nCellsSolve', ncs); call reg_i0('nVertLevels', nvl)
  call reg_i0('index_qv', iqv); call reg_i0('index_qc', iqc); call reg_i0('index_qr', iqr)
  call reg_i0('index_qi', iqi); call reg_i0('index_qs', iqs); call reg_i0('index_qg', iqg)
  call reg_i0('index_ni', ini); call reg_i0('index_nr', inr)

  ! physics_run_init (mpas_atmphys_manager.F:653-694), one rank
  ims = 1; ime = nC; jms = 1; jme = 1; kms = 1; kme = nV + 1
  ids = ims; ide = ime + 1; jds = jms; jde = jme + 1; kds = kms; kde = kme
  its = ims; ite = ime; jts = jms; jte = jme; kts = kms; kte = kme - 1
  dt_dyn = dt
  dt_microp = dt_dyn
  n_microp  = 1
  dt_microp = 90._RKIND
  n_microp  = max(nint(dt_dyn/dt_microp),1)
  dt_microp = dt_dyn / n_microp
  if(dt_dyn <= dt_microp) dt_microp = dt_dyn
  l_diags = diags; l_acrain = .false.
  has_reqc = 0; has_reqi = 0; has_reqs = 0
  ! physics_tables_init (mpas_atmphys_control.F:406)
  l_mp_tables = .true.
  ! init_microphysics (DRV:268-314): thompson_init, init_thompson_clouddroplets_forMPAS
  call init_microphysics(dminfo, configs, mesh, state, 2, sfc_input, diag_physics)
  call dump_thompson_init(trim(finit)//c_null_char)

  open(11, file=trim(fout), access='stream', form='unformatted', status='replace', convert='little_endian')
  write(11) nt_c, mu_c, dt_microp, real(n_microp, RKIND)
  do step = 1, nsteps
    call driver_microphysics(configs, mesh, state, 2, diag, diag_physics, tend_physics, tend, step, 1, nC)
    write(11) scalars, theta_m, rtheta_p, exner, pressure_p, surface_pressure, tend_sfc_pressure, rt_diabatic_tend, &
      dtheta_dt_mp, precipw, rainnc, rainncv, snownc, snowncv, graupelnc, graupelncv, sr, rainprod, evapprod, &
      re_cloud, re_ice, re_snow, rthmpten, rqvmpten, rqcmpten, rqrmpten, rqimpten, rqsmpten, rqgmpten, rnimpten, &
      rnrmpten
    if (diags) write(11) refl10cm, refl10cm_max, refl10cm_1km, refl10cm_1km_max, relhum
  end do
  close(11)
end program thompson_column_ref
