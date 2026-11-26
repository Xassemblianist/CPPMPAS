! SPDX-License-Identifier: GPL-3.0-or-later
! Column driver of MPAS's own Noah path (build_ref.sh): reads the pool arrays from in.bin, sets the module variables
! physics_run_init sets (one rank: its..ite = ims..ime = 1..nCells, kms = 1, kme = nV+1), calls init_lsm, then per
! step what physics_driver calls for config_lsm_scheme = 'sf_noah' (mpas_atmphys_driver.F:283-315): allocate_lsm,
! driver_lsm, deallocate_lsm, allocate_seaice, driver_seaice, deallocate_seaice, atmphys_sfc_diagnostics; the *_sea
! arrays of the surface layer are given (the sfclayer driver allocates them every step; deallocate_seaice frees
! them). Writes the tables and the pool arrays after init_lsm and after every step to out.bin. Native-endian
! streams (this file is compiled without -fconvert=big-endian). Layout: tests/noah_column/columns.py.
program noah_column_ref
  use mpas_kind_types
  use mpas_pool_routines
  use mpas_atmphys_vars
  use mpas_atmphys_landuse, only: isurban
  use mpas_atmphys_driver_lsm, only: init_lsm, allocate_lsm, deallocate_lsm, driver_lsm
  use mpas_atmphys_driver_seaice, only: allocate_seaice, deallocate_seaice, driver_seaice
  use mpas_atmphys_sfc_diagnostics, only: atmphys_sfc_diagnostics
  use module_sf_noahlsm
  implicit none
  integer, parameter :: NR1 = 53, NS2 = 5, NA = 4, NSEA = 9
  character(len=24), parameter :: r1names(NR1) = [character(len=24) :: &
    'acsnom', 'acsnow', 'canwat', 'chs', 'chs2', 'chklowq', 'cpm', 'cqs2', 'glw', 'grdflx', 'gsw', 'hfx', 'lai', &
    'lh', 'noahres', 'potevp', 'qfx', 'qgh', 'qsfc', 'br', 'sfc_albedo', 'sfc_emibck', 'sfc_emiss', 'sfcrunoff', &
    'smstav', 'smstot', 'snotime', 'snopcx', 'swddif', 'swddir', 'udrunoff', 'z0', 'znt', 'sr', 'rainncv', 'raincv', &
    'q2', 't2m', 'th2m', 'shdmin', 'shdmax', 'snoalb', 'sfc_albbck', 'snow', 'snowc', 'snowh', 'tmn', 'skintemp', &
    'vegfra', 'xice', 'xland', 'surface_pressure', 'spare']
  character(len=24), parameter :: s2names(NS2) = [character(len=24) :: 'dzs', 'sh2o', 'smcrel', 'smois', 'tslb']
  integer(4) :: hdr(10)
  integer :: nC, nsoil, nsteps, nV, step, k
  integer, target :: nCells_t, nSoil_t, nCellsSolve_t, isice_t
  real(kind=RKIND) :: dt
  real(kind=RKIND), allocatable, target :: r1(:,:), s2(:,:,:), a3(:,:,:), sea(:,:)
  integer, allocatable, target :: ivgtyp(:), isltyp(:)
  character(len=StrKIND), target :: mminlu, c_lsm, c_mp, c_cu, c_soil
  logical, target :: l_albedo, l_snowalb, l_restart
  character(len=64) :: mminlu_in
  type(mpas_pool_type) :: configs, mesh, diag, diag_physics, sfc_input, output_noahmp
  type(dm_info) :: dminfo
  character(len=512) :: fin, fout

  call get_command_argument(1, fin)
  call get_command_argument(2, fout)
  open(10, file=trim(fin), access='stream', form='unformatted', status='old')
  read(10) hdr
  nC = hdr(1); nsoil = hdr(2); nsteps = hdr(3); nV = hdr(4)
  read(10) dt, mminlu_in
  allocate(ivgtyp(nC), isltyp(nC), r1(nC,NR1), s2(nsoil,nC,NS2), a3(nC,nV+1,NA), sea(nC,NSEA))
  read(10) ivgtyp, isltyp, r1, s2, a3, sea
  close(10)
  isurban = hdr(5); isice_t = hdr(6)
  l_albedo = hdr(7) /= 0; l_snowalb = hdr(8) /= 0; l_restart = hdr(9) /= 0
  nCells_t = nC; nCellsSolve_t = nC; nSoil_t = nsoil
  mminlu = mminlu_in; c_lsm = 'sf_noah'; c_mp = 'mp_wsm6'; c_cu = 'cu_ntiedtke'; c_soil = 'STAS'

  do k = 1, NR1
    call reg_r1(trim(r1names(k)), r1(:,k))
  end do
  do k = 1, NS2
    call reg_r2(trim(s2names(k)), s2(:,:,k))
  end do
  call reg_i1('ivgtyp', ivgtyp); call reg_i1('isltyp', isltyp)
  call reg_i0('nCells', nCells_t); call reg_i0('nCellsSolve', nCellsSolve_t); call reg_i0('nSoilLevels', nSoil_t)
  call reg_i0('isice', isice_t); call reg_c0('mminlu', mminlu)
  call reg_c0('config_lsm_scheme', c_lsm); call reg_c0('config_microp_scheme', c_mp)
  call reg_c0('config_convection_scheme', c_cu); call reg_c0('input_soil_data', c_soil)
  call reg_l0('config_sfc_albedo', l_albedo); call reg_l0('config_sfc_snowalbedo', l_snowalb)
  call reg_l0('config_do_restart', l_restart)

  ! physics_run_init (mpas_atmphys_manager.F:653-663, 678; 542; 740-744 with config_frac_seaice = hdr(10))
  ims = 1; ime = nC; jms = 1; jme = 1; kms = 1; kme = nV + 1
  ids = ims; ide = ime + 1; jds = jms; jde = jme + 1; kds = kms; kde = kme
  its = ims; ite = ime; jts = jms; jte = jme; kts = kms; kte = kme - 1
  num_soils = nsoil
  dt_pbl = dt
  if (hdr(10) == 0) then
    xice_threshold = 0.5
  else
    xice_threshold = 0.02
  end if
  ! MPAS_to_physics (mpas_atmphys_interface.F): the physics arrays, given
  allocate(dz_p(ims:ime,kms:kme,jms:jme), t_p(ims:ime,kms:kme,jms:jme), qv_p(ims:ime,kms:kme,jms:jme), &
           pres2_hyd_p(ims:ime,kms:kme,jms:jme))
  dz_p(:,:,1) = a3(:,:,1); t_p(:,:,1) = a3(:,:,2); qv_p(:,:,1) = a3(:,:,3); pres2_hyd_p(:,:,1) = a3(:,:,4)

  call init_lsm(dminfo, mesh, configs, diag_physics, sfc_input)
  open(11, file=trim(fout), access='stream', form='unformatted', status='replace')
  write(11) lucats, bare, natural, lcz_1, lcz_2, lcz_3, lcz_4, lcz_5, lcz_6, lcz_7, lcz_8, lcz_9, lcz_10, lcz_11, &
    slcats, slpcats, nrotbl
  write(11) shdtbl, rstbl, rgltbl, hstbl, snuptbl, maxalb, laimintbl, laimaxtbl, emissmintbl, emissmaxtbl, &
    albedomintbl, albedomaxtbl, z0mintbl, z0maxtbl, ztopvtbl, zbotvtbl, topt_data, cmcmax_data, cfactr_data, &
    rsmax_data, bb, drysmc, f11, maxsmc, refsmc, satpsi, satdk, satdw, wltsmc, qtz, slope_data, sbeta_data, &
    fxexp_data, csoil_data, salp_data, refdk_data, refkdt_data, frzk_data, zbot_data, czil_data, smlow_data, &
    smhigh_data, lvcoef_data
  write(11) r1, s2
  do step = 1, nsteps
    call allocate_lsm
    call driver_lsm(step, configs, mesh, diag_physics, sfc_input, its, ite)
    call deallocate_lsm
    ! the surface layer's open-water values (mpas_atmphys_driver_sfclayer.F allocates them every step)
    allocate(cpm_sea(ims:ime,jms:jme), chs_sea(ims:ime,jms:jme), chs2_sea(ims:ime,jms:jme), &
             cqs2_sea(ims:ime,jms:jme), hfx_sea(ims:ime,jms:jme), qfx_sea(ims:ime,jms:jme), &
             qgh_sea(ims:ime,jms:jme), qsfc_sea(ims:ime,jms:jme), lh_sea(ims:ime,jms:jme))
    cpm_sea(:,1) = sea(:,1); chs_sea(:,1) = sea(:,2); chs2_sea(:,1) = sea(:,3); cqs2_sea(:,1) = sea(:,4)
    hfx_sea(:,1) = sea(:,5); qfx_sea(:,1) = sea(:,6); qgh_sea(:,1) = sea(:,7); qsfc_sea(:,1) = sea(:,8)
    lh_sea(:,1) = sea(:,9)
    call allocate_seaice(configs)
    call driver_seaice(configs, diag_physics, sfc_input, its, ite)
    call deallocate_seaice(configs)
    call atmphys_sfc_diagnostics(configs, mesh, diag, diag_physics, sfc_input, output_noahmp, its, ite)
    write(11) r1, s2
  end do
  close(11)
end program noah_column_ref
