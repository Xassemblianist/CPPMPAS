! SPDX-License-Identifier: GPL-3.0-or-later
! Column driver of MPAS's own Grell-Freitas convection, physics interface and physics tendencies (build_ref.sh): reads
! the pool arrays from in.bin (records, tests/gf_column/columns.py), sets the module variables physics_run_init sets
! (one rank, one thread: its..ite = 1..nCellsSolve, kts..kte = 1..nVertLevels), calls per step what physics_driver
! (mpas_atmphys_driver.F:202-370) and the dynamics call for the convection: allocate_forall_physics, MPAS_to_physics,
! update_convection_step1, allocate / driver / deallocate_convection, update_convection_step2, then physics_get_tend
! with rk_step 1 and 2, and writes the snapshots of columns.py to out.bin. The module recio is
! tests/cu_gwdo_column/harness.F90's. Compiled without -fconvert=big-endian (native-endian files).
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

program gf_column_ref
  use mpas_kind_types
  use mpas_pool_routines
  use mpas_derived_types
  use mpas_atm_dimensions, only: nVertLevels, num_scalars
  use mpas_atmphys_vars
  use mpas_atmphys_interface, only: allocate_forall_physics, deallocate_forall_physics, MPAS_to_physics
  use mpas_atmphys_driver_convection, only: allocate_convection, deallocate_convection, driver_convection, &
    update_convection_step1, update_convection_step2
  use mpas_atmphys_todynamics, only: physics_get_tend
  use recio
  implicit none
  interface
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
  integer, target :: iqv = 1, iqc = 2, iqr = 3, iqi = 4, iqs = 5, iqg = 6, inone = -1, gf_deep, gf_shallow
  real(kind=RKIND), target :: len_disp, bucket_rainc
  logical, target :: do_restart, do_dacycling = .false.
  character(len=StrKIND), target :: microp = 'mp_wsm6', pbl = 'bl_ysu', convection = 'cu_grell_freitas', &
    radt_lw = 'rrtmg_lw', radt_sw = 'rrtmg_sw'
  real(kind=RKIND), pointer :: tend_ru(:,:), tend_rth(:,:), tend_rho(:,:)

  call get_command_argument(1, fin)
  call get_command_argument(2, fout)
  open(10, file=trim(fin), access='stream', form='unformatted', status='old')
  call read_all(10)
  close(10)
  hdr = recs(rec_index('hdr'))%iv(1:10)
  nC = hdr(1); nCt = hdr(2); nEs = hdr(3); nE = hdr(4); nV = hdr(5); nsteps = hdr(6)
  do_restart = hdr(8) == 1
  gf_deep = hdr(9); gf_shallow = hdr(10)
  len_disp = recs(rec_index('dtl'))%r(2)
  bucket_rainc = recs(rec_index('dtl'))%r(3)

  d_nCells = nCt; d_nEdges = nE; d_nCellsSolve = nC; d_nEdgesSolve = nEs
  call reg_i0('nCells', d_nCells); call reg_i0('nEdges', d_nEdges)
  call reg_i0('nCellsSolve', d_nCellsSolve); call reg_i0('nEdgesSolve', d_nEdgesSolve)
  call reg_i0('index_qv', iqv); call reg_i0('index_qc', iqc); call reg_i0('index_qr', iqr)
  call reg_i0('index_qi', iqi); call reg_i0('index_qs', iqs); call reg_i0('index_qg', iqg)
  call reg_i0('index_nc', inone); call reg_i0('index_ni', inone); call reg_i0('index_nifa', inone)
  call reg_i0('index_nwfa', inone)
  call reg_i0('config_gfconv_closure_deep', gf_deep); call reg_i0('config_gfconv_closure_shallow', gf_shallow)
  call reg_cc('config_microp_scheme', microp); call reg_cc('config_pbl_scheme', pbl)
  call reg_cc('config_convection_scheme', convection)
  call reg_cc('config_radt_lw_scheme', radt_lw); call reg_cc('config_radt_sw_scheme', radt_sw)
  call reg_cl('config_do_restart', do_restart); call reg_cl('config_do_DAcycling', do_dacycling)
  call reg_cr('config_len_disp', len_disp); call reg_cr('config_bucket_rainc', bucket_rainc)
  call mpas_pool_get_array(diag, 'tend_ru_physics', tend_ru)
  call mpas_pool_get_array(diag, 'tend_rtheta_physics', tend_rth)
  call mpas_pool_get_array(diag, 'tend_rho_physics', tend_rho)
  block%domain => domain

  ! physics_run_init (mpas_atmphys_manager.F:540-544, 661-669, 700-702) and init_physics_flags (mpas_atmphys_init.F)
  nVertLevels = nV; num_scalars = 6
  ims = 1; ime = nC; jms = 1; jme = 1; kms = 1; kme = nV + 1
  ids = ims; ide = ime + 1; jds = jms; jde = jme + 1; kds = kms; kde = kme
  its = ims; ite = ime; jts = jms; jte = jme; kts = kms; kte = kme - 1
  dt_dyn = recs(rec_index('dtl'))%r(1)
  dt_pbl = dt_dyn
  n_cu = 1
  dt_cu = dt_dyn
  l_acrain = hdr(7) == 1
  l_conv = .true.
  f_qc = .true.; f_qr = .true.; f_qi = .true.; f_qs = .true.; f_qg = .true.; f_qoz = .false.
  f_nc = .false.; f_ni = .false.; f_nifa = .false.; f_nwfa = .false.; f_nbca = .false.

  open(11, file=trim(fout), access='stream', form='unformatted', status='replace')
  do step = 1, nsteps
    call allocate_forall_physics(configs)
    call MPAS_to_physics(configs, mesh, state, 1, diag, diag_physics, 1, nC)
    call marker(11, '#A', step)
    call dump_all(11)
    call update_convection_step1(configs, diag_physics, tend_physics, 1, nC)
    if (l_conv) then
      call allocate_convection(configs)
      call driver_convection(step, configs, mesh, sfc_input, diag_physics, tend_physics, 1, nC)
      call deallocate_convection(configs)
    end if
    call update_convection_step2(configs, diag_physics, 1, nC)
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
end program gf_column_ref
