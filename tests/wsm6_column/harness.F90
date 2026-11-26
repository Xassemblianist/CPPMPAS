! SPDX-License-Identifier: GPL-3.0-or-later
! Column driver of MPAS's own WSM6 path (build_ref.sh): reads the pool arrays from in.bin, sets the module
! variables physics_run_init sets (one rank: its..ite = 1..nCells, kts..kte = 1..nVertLevels), calls
! driver_microphysics nsteps times and writes the pool arrays it changed after every step to out.bin.
! Both files are native-endian streams (this file is compiled without -fconvert=big-endian). Layout:
! tests/wsm6_column/columns.py.
program wsm6_column_ref
  use mpas_kind_types
  use mpas_pool_routines
  use mpas_atmphys_vars
  use drv_extract, only: driver_microphysics
  implicit none
  integer(4) :: hdr(5)
  integer :: nC, nV, nsteps, step, nscal
  integer, target :: mp_top, iqv = 1, iqc = 2, iqr = 3, iqi = 4, iqs = 5, iqg = 6
  real(kind=RKIND) :: dt
  real(kind=RKIND), allocatable, target :: zgrid(:,:), zz(:,:), rho_zz(:,:), theta_m(:,:), w(:,:), scalars(:,:,:), &
    exner(:,:), pressure_base(:,:), pressure_p(:,:), exner_base(:,:), rtheta_base(:,:), rtheta_p(:,:), &
    surface_pressure(:), rainprod(:,:), evapprod(:,:), re_cloud(:,:), re_ice(:,:), re_snow(:,:), rainnc(:), snownc(:), &
    graupelnc(:), rainncv(:), snowncv(:), graupelncv(:), sr(:), precipw(:), tend_sfc_pressure(:), rt_diabatic_tend(:,:), &
    dtheta_dt_mp(:,:), refl10cm(:,:), rthmpten(:,:), rqvmpten(:,:), rqcmpten(:,:), rqrmpten(:,:), rqimpten(:,:), &
    rqsmpten(:,:), rqgmpten(:,:)
  integer, allocatable, target :: i_rainnc(:)
  type(mpas_pool_type) :: configs, mesh, state, diag, diag_physics, tend_physics, tend
  character(len=512) :: fin, fout

  call get_command_argument(1, fin)
  call get_command_argument(2, fout)
  open(10, file=trim(fin), access='stream', form='unformatted', status='old')
  read(10) hdr
  nC = hdr(1); nV = hdr(2); nsteps = hdr(3); mp_top = hdr(4); nscal = hdr(5)
  read(10) dt
  allocate(zgrid(nV+1,nC), zz(nV,nC), rho_zz(nV,nC), theta_m(nV,nC), w(nV+1,nC), scalars(nscal,nV,nC), &
    exner(nV,nC), pressure_base(nV,nC), pressure_p(nV,nC), exner_base(nV,nC), rtheta_base(nV,nC), rtheta_p(nV,nC), &
    surface_pressure(nC), rainprod(nV,nC), evapprod(nV,nC), re_cloud(nV,nC), re_ice(nV,nC), re_snow(nV,nC), &
    rainnc(nC), snownc(nC), graupelnc(nC), rainncv(nC), snowncv(nC), graupelncv(nC), sr(nC), precipw(nC), &
    tend_sfc_pressure(nC), rt_diabatic_tend(nV,nC), dtheta_dt_mp(nV,nC), refl10cm(nV,nC), rthmpten(nV,nC), &
    rqvmpten(nV,nC), rqcmpten(nV,nC), rqrmpten(nV,nC), rqimpten(nV,nC), rqsmpten(nV,nC), rqgmpten(nV,nC), i_rainnc(nC))
  read(10) zgrid, zz, rho_zz, theta_m, w, scalars, exner, pressure_base, pressure_p, exner_base, rtheta_base, &
    rtheta_p, surface_pressure, rainprod, evapprod, re_cloud, re_ice, re_snow, rainnc, snownc, graupelnc
  close(10)
  ! values the driver must overwrite before reading them
  rainncv = -999.; snowncv = -999.; graupelncv = -999.; sr = -999.; precipw = -999.; tend_sfc_pressure = -999.
  rt_diabatic_tend = -999.; dtheta_dt_mp = -999.; refl10cm = 0.; i_rainnc = 0
  rthmpten = -999.; rqvmpten = -999.; rqcmpten = -999.; rqrmpten = -999.; rqimpten = -999.; rqsmpten = -999.
  rqgmpten = -999.

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
  call reg_r2('rqgmpten', rqgmpten); call reg_i0('mp_top_level', mp_top)
  call reg_i0('index_qv', iqv); call reg_i0('index_qc', iqc); call reg_i0('index_qr', iqr)
  call reg_i0('index_qi', iqi); call reg_i0('index_qs', iqs); call reg_i0('index_qg', iqg)

  ! physics_run_init (mpas_atmphys_manager.F:653-708), one rank
  ims = 1; ime = nC; jms = 1; jme = 1; kms = 1; kme = nV + 1
  ids = ims; ide = ime + 1; jds = jms; jde = jme + 1; kds = kms; kde = kme
  its = ims; ite = ime; jts = jms; jte = jme; kts = kms; kte = kme - 1
  dt_dyn = dt; dt_microp = dt_dyn; n_microp = 1
  l_diags = .false.; l_acrain = .false.
  has_reqc = 0; has_reqi = 0; has_reqs = 0
  ! init_microphysics (DRV:305-307)
  call init_wsm6()

  open(11, file=trim(fout), access='stream', form='unformatted', status='replace')
  do step = 1, nsteps
    call driver_microphysics(configs, mesh, state, 2, diag, diag_physics, tend_physics, tend, step, 1, nC)
    write(11) scalars, theta_m, rtheta_p, exner, pressure_p, surface_pressure, tend_sfc_pressure, rt_diabatic_tend, &
      dtheta_dt_mp, precipw, rainnc, rainncv, snownc, snowncv, graupelnc, graupelncv, sr, rainprod, evapprod, &
      re_cloud, re_ice, re_snow, rthmpten, rqvmpten, rqcmpten, rqrmpten, rqimpten, rqsmpten, rqgmpten
  end do
  close(11)
contains
  subroutine init_wsm6()
    use mp_wsm6, only: mp_wsm6_init
    use mpas_atmphys_constants
    character(len=512) :: errmsg
    integer :: errflg
    integer, parameter :: hail_opt = 0  ! DRV:91
    call mp_wsm6_init(den0=rho_a,denr=rho_r,dens=rho_s,cl=cliq,cpv=cpv,hail_opt=hail_opt,errmsg=errmsg,errflg=errflg)
  end subroutine
end program wsm6_column_ref
