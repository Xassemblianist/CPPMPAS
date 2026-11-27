! SPDX-License-Identifier: GPL-3.0-or-later
! Stand-ins for the MPAS framework modules that the MPAS radiation sources of the radiation column check use (no
! arithmetic): a pool that hands out pointers the harness registered by name (arrays, dimensions and namelist
! options), empty timers and log, dm_info / IO_NODE and one-rank broadcasts (no-ops), Fortran units, an empty
! Noah-MP type, and the four time variables of mpas_atmphys_manager the radiation drivers use (the harness sets them
! as physics_run_init / physics_timetracker do, mpas_atmphys_manager.F:205-212, 457-463).
module mpas_log
contains
  subroutine mpas_log_write(msg, intArgs, realArgs, logicArgs, messageType)
    character(len=*), intent(in) :: msg
    integer, intent(in), optional :: intArgs(:), messageType
    real(kind=8), intent(in), optional :: realArgs(:)
    logical, intent(in), optional :: logicArgs(:)
    print *, 'mpas_log_write: ', trim(msg)
    ! MPAS_LOG_CRIT (= 4, mpas_log_types.inc:8) aborts the run as MPAS's mpas_log_write does; a nonzero exit
    ! status, so the check cannot take the run for a finished one
    if (present(messageType)) then
      if (messageType == 4) error stop 'MPAS_LOG_CRIT'
    end if
  end subroutine
end module mpas_log

module mpas_timer
contains
  subroutine mpas_timer_start(n)
    character(len=*), intent(in) :: n
  end subroutine
  subroutine mpas_timer_stop(n)
    character(len=*), intent(in) :: n
  end subroutine
end module mpas_timer

module NoahmpIOVarType
  type NoahmpIO_type
    integer :: dummy = 0
  end type
end module NoahmpIOVarType

module mpas_derived_types
  integer, parameter :: MPAS_LOG_ERR = 3, MPAS_LOG_CRIT = 4, IO_NODE = 0
  type dm_info
    integer :: nprocs = 1, my_proc_id = 0
  end type
  type domain_type
    integer :: dummy = 0
  end type
  type block_type
    type(domain_type), pointer :: domain => null()
  end type
end module mpas_derived_types

module mpas_dmpar
  use mpas_kind_types
  use mpas_derived_types
contains
  subroutine mpas_dmpar_bcast_int(dminfo, i, proc)
    type(dm_info), intent(in) :: dminfo
    integer, intent(inout) :: i
    integer, intent(in), optional :: proc
  end subroutine
  subroutine mpas_dmpar_bcast_real(dminfo, r, proc)
    type(dm_info), intent(in) :: dminfo
    real(kind=RKIND), intent(inout) :: r
    integer, intent(in), optional :: proc
  end subroutine
  subroutine mpas_dmpar_bcast_reals(dminfo, n, r, proc)
    type(dm_info), intent(in) :: dminfo
    integer, intent(in) :: n
    real(kind=RKIND), intent(inout) :: r(n)
    integer, intent(in), optional :: proc
  end subroutine
  subroutine mpas_dmpar_bcast_doubles(dminfo, n, r, proc)
    type(dm_info), intent(in) :: dminfo
    integer, intent(in) :: n
    real(kind=R8KIND), intent(inout) :: r(n)
    integer, intent(in), optional :: proc
  end subroutine
  subroutine mpas_dmpar_max_int(dminfo, i, imax)
    type(dm_info), intent(in) :: dminfo
    integer, intent(in) :: i
    integer, intent(out) :: imax
    imax = i
  end subroutine
end module mpas_dmpar

module mpas_io_units
contains
  subroutine mpas_new_unit(newUnit, unformatted)
    integer, intent(inout) :: newUnit
    logical, optional, intent(in) :: unformatted
    logical :: opened
    do newUnit = 20, 99
      inquire(unit=newUnit, opened=opened)
      if (.not. opened) return
    end do
    newUnit = -1
  end subroutine
  subroutine mpas_release_unit(unit)
    integer, intent(in) :: unit
  end subroutine
end module mpas_io_units

module mpas_atmphys_manager
  use mpas_kind_types
  integer :: julday, year
  real(kind=RKIND) :: gmt, curr_julday
end module mpas_atmphys_manager

module mpas_atm_dimensions
  integer :: nVertLevels, maxEdges, maxEdges2, num_scalars
end module mpas_atm_dimensions

module mpas_pool_routines
  use mpas_kind_types
  use mpas_log
  use mpas_derived_types
  implicit none
  type mpas_pool_type
    integer :: dummy = 0
  end type
  integer, parameter :: nmax = 256
  type r1p
    real(kind=RKIND), pointer :: p(:) => null()
  end type
  type r2p
    real(kind=RKIND), pointer :: p(:,:) => null()
  end type
  type r3p
    real(kind=RKIND), pointer :: p(:,:,:) => null()
  end type
  type i1p
    integer, pointer :: p(:) => null()
  end type
  type i2p
    integer, pointer :: p(:,:) => null()
  end type
  type i0p
    integer, pointer :: p => null()
  end type
  type rcp_
    real(kind=RKIND), pointer :: p => null()
  end type
  type lcp
    logical, pointer :: p => null()
  end type
  type ccp
    character(len=StrKIND), pointer :: p => null()
  end type
  character(len=64) :: n1(nmax) = '', n2(nmax) = '', n3(nmax) = '', ni1(nmax) = '', ni2(nmax) = '', ni0(nmax) = '', &
    nrc(nmax) = '', nlc(nmax) = '', ncc(nmax) = ''
  type(r1p) :: v1(nmax)
  type(r2p) :: v2(nmax)
  type(r3p) :: v3(nmax)
  type(i1p) :: vi1(nmax)
  type(i2p) :: vi2(nmax)
  type(i0p) :: vi0(nmax)
  type(rcp_) :: vrc(nmax)
  type(lcp) :: vlc(nmax)
  type(ccp) :: vcc(nmax)
  interface mpas_pool_get_array
    module procedure get_r1, get_r2, get_r3, get_i1, get_i2, get_i0
  end interface
  interface mpas_pool_get_config
    module procedure cfg_c, cfg_r, cfg_l, cfg_i
  end interface
contains
  integer function slot(names, name)
    character(len=64), intent(inout) :: names(nmax)
    character(len=*), intent(in) :: name
    do slot = 1, nmax
      if (trim(names(slot)) == trim(name)) return
      if (names(slot) == '') then
        names(slot) = name
        return
      end if
    end do
    error stop 'pool stub full'
  end function
  integer function find(names, name)
    character(len=64), intent(in) :: names(nmax)
    character(len=*), intent(in) :: name
    do find = 1, nmax
      if (trim(names(find)) == trim(name)) return
    end do
    print *, 'pool stub: no entry ', trim(name)
    stop 1
  end function
  subroutine reg_r1(name, a)
    character(len=*), intent(in) :: name
    real(kind=RKIND), target :: a(:)
    v1(slot(n1, name))%p => a
  end subroutine
  subroutine reg_r2(name, a)
    character(len=*), intent(in) :: name
    real(kind=RKIND), target :: a(:,:)
    v2(slot(n2, name))%p => a
  end subroutine
  subroutine reg_r3(name, a)
    character(len=*), intent(in) :: name
    real(kind=RKIND), target :: a(:,:,:)
    v3(slot(n3, name))%p => a
  end subroutine
  subroutine reg_i1(name, a)
    character(len=*), intent(in) :: name
    integer, target :: a(:)
    vi1(slot(ni1, name))%p => a
  end subroutine
  subroutine reg_i2(name, a)
    character(len=*), intent(in) :: name
    integer, target :: a(:,:)
    vi2(slot(ni2, name))%p => a
  end subroutine
  subroutine reg_i0(name, a)
    character(len=*), intent(in) :: name
    integer, target :: a
    vi0(slot(ni0, name))%p => a
  end subroutine
  subroutine reg_cr(name, a)
    character(len=*), intent(in) :: name
    real(kind=RKIND), target :: a
    vrc(slot(nrc, name))%p => a
  end subroutine
  subroutine reg_cl(name, a)
    character(len=*), intent(in) :: name
    logical, target :: a
    vlc(slot(nlc, name))%p => a
  end subroutine
  subroutine reg_cc(name, a)
    character(len=*), intent(in) :: name
    character(len=StrKIND), target :: a
    vcc(slot(ncc, name))%p => a
  end subroutine
  subroutine get_r1(pool, name, p, tl)
    type(mpas_pool_type), intent(in) :: pool
    character(len=*), intent(in) :: name
    real(kind=RKIND), pointer :: p(:)
    integer, intent(in), optional :: tl
    p => v1(find(n1, name))%p
  end subroutine
  subroutine get_r2(pool, name, p, tl)
    type(mpas_pool_type), intent(in) :: pool
    character(len=*), intent(in) :: name
    real(kind=RKIND), pointer :: p(:,:)
    integer, intent(in), optional :: tl
    p => v2(find(n2, name))%p
  end subroutine
  subroutine get_r3(pool, name, p, tl)
    type(mpas_pool_type), intent(in) :: pool
    character(len=*), intent(in) :: name
    real(kind=RKIND), pointer :: p(:,:,:)
    integer, intent(in), optional :: tl
    p => v3(find(n3, name))%p
  end subroutine
  subroutine get_i1(pool, name, p, tl)
    type(mpas_pool_type), intent(in) :: pool
    character(len=*), intent(in) :: name
    integer, pointer :: p(:)
    integer, intent(in), optional :: tl
    p => vi1(find(ni1, name))%p
  end subroutine
  subroutine get_i2(pool, name, p, tl)
    type(mpas_pool_type), intent(in) :: pool
    character(len=*), intent(in) :: name
    integer, pointer :: p(:,:)
    integer, intent(in), optional :: tl
    p => vi2(find(ni2, name))%p
  end subroutine
  subroutine get_i0(pool, name, p, tl)
    type(mpas_pool_type), intent(in) :: pool
    character(len=*), intent(in) :: name
    integer, pointer :: p
    integer, intent(in), optional :: tl
    p => vi0(find(ni0, name))%p
  end subroutine
  subroutine mpas_pool_get_dimension(pool, name, p)
    type(mpas_pool_type), intent(in) :: pool
    character(len=*), intent(in) :: name
    integer, pointer :: p
    p => vi0(find(ni0, name))%p
  end subroutine
  subroutine cfg_c(pool, name, p)
    type(mpas_pool_type), intent(in) :: pool
    character(len=*), intent(in) :: name
    character(len=StrKIND), pointer :: p
    p => vcc(find(ncc, name))%p
  end subroutine
  subroutine cfg_r(pool, name, p)
    type(mpas_pool_type), intent(in) :: pool
    character(len=*), intent(in) :: name
    real(kind=RKIND), pointer :: p
    p => vrc(find(nrc, name))%p
  end subroutine
  subroutine cfg_l(pool, name, p)
    type(mpas_pool_type), intent(in) :: pool
    character(len=*), intent(in) :: name
    logical, pointer :: p
    p => vlc(find(nlc, name))%p
  end subroutine
  subroutine cfg_i(pool, name, p)
    type(mpas_pool_type), intent(in) :: pool
    character(len=*), intent(in) :: name
    integer, pointer :: p
    print *, 'pool stub: integer config ', trim(name), ' not provided'
    stop 1
  end subroutine
end module mpas_pool_routines
