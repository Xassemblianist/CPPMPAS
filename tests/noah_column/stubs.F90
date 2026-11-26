! SPDX-License-Identifier: GPL-3.0-or-later
! Stand-ins for the MPAS framework modules that the extracted MPAS physics sources use (no arithmetic): derived types
! (pool, dm_info, log levels), a log that stops on MPAS_LOG_CRIT as MPAS's does, empty timers, one-rank dmpar
! broadcasts (no-ops), file units, the isurban module variable of mpas_atmphys_landuse (set by the harness, as
! landuse_init_forMPAS does), an empty Noah-MP type, and a pool that hands out pointers the harness registered by name.
module mpas_derived_types
  implicit none
  integer, parameter :: MPAS_LOG_OUT = 1, MPAS_LOG_WARN = 2, MPAS_LOG_ERR = 3, MPAS_LOG_CRIT = 4
  type mpas_pool_type
    integer :: dummy = 0
  end type
  type dm_info
    integer :: nprocs = 1, my_proc_id = 0
  end type
end module mpas_derived_types

module mpas_log
  use mpas_derived_types, only: MPAS_LOG_CRIT
contains
  subroutine mpas_log_write(msg, intArgs, realArgs, logicArgs, messageType)
    character(len=*), intent(in) :: msg
    integer, intent(in), optional :: intArgs(:), messageType
    real(kind=8), intent(in), optional :: realArgs(:)
    logical, intent(in), optional :: logicArgs(:)
    if (present(messageType)) then
      if (messageType /= 1) print *, trim(msg)
      if (messageType == MPAS_LOG_CRIT) stop 1
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

module mpas_dmpar
  use mpas_derived_types
  implicit none
  integer, parameter :: IO_NODE = 0
contains
  subroutine mpas_dmpar_bcast_char(dminfo, c)
    type(dm_info), intent(in) :: dminfo
    character(len=*), intent(inout) :: c
  end subroutine
  subroutine mpas_dmpar_bcast_int(dminfo, i)
    type(dm_info), intent(in) :: dminfo
    integer, intent(inout) :: i
  end subroutine
  subroutine mpas_dmpar_bcast_ints(dminfo, n, a)
    type(dm_info), intent(in) :: dminfo
    integer, intent(in) :: n
    integer, intent(inout) :: a(n)
  end subroutine
  subroutine mpas_dmpar_bcast_real(dminfo, r)
    type(dm_info), intent(in) :: dminfo
    real(kind=8), intent(inout) :: r
  end subroutine
  subroutine mpas_dmpar_bcast_reals(dminfo, n, a)
    type(dm_info), intent(in) :: dminfo
    integer, intent(in) :: n
    real(kind=8), intent(inout) :: a(n)
  end subroutine
end module mpas_dmpar

module mpas_io_units
contains
  subroutine mpas_new_unit(newUnit, unformatted)
    integer, intent(inout) :: newUnit
    logical, optional, intent(in) :: unformatted
    newUnit = 41
  end subroutine
  subroutine mpas_release_unit(releasedUnit)
    integer, intent(in) :: releasedUnit
  end subroutine
end module mpas_io_units

module mpas_atmphys_landuse
  implicit none
  integer :: isurban
end module mpas_atmphys_landuse

module mpas_pool_routines
  use mpas_kind_types
  use mpas_derived_types
  use mpas_log
  implicit none
  integer, parameter :: nmax = 128
  type r1p
    real(kind=RKIND), pointer :: p(:) => null()
  end type
  type r2p
    real(kind=RKIND), pointer :: p(:,:) => null()
  end type
  type i1p
    integer, pointer :: p(:) => null()
  end type
  type i0p
    integer, pointer :: p => null()
  end type
  type c0p
    character(len=StrKIND), pointer :: p => null()
  end type
  type l0p
    logical, pointer :: p => null()
  end type
  character(len=64) :: n1(nmax) = '', n2(nmax) = '', ni1(nmax) = '', ni0(nmax) = '', nc0(nmax) = '', nl0(nmax) = ''
  type(r1p) :: v1(nmax)
  type(r2p) :: v2(nmax)
  type(i1p) :: vi1(nmax)
  type(i0p) :: vi0(nmax)
  type(c0p) :: vc0(nmax)
  type(l0p) :: vl0(nmax)
  interface mpas_pool_get_array
    module procedure get_r1, get_r2, get_i1, get_i0, get_c0
  end interface
  interface mpas_pool_get_config
    module procedure get_c0, get_l0
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
    stop 'pool stub full'
  end function
  integer function find(names, name)
    character(len=64), intent(in) :: names(nmax)
    character(len=*), intent(in) :: name
    do find = 1, nmax
      if (trim(names(find)) == trim(name)) return
    end do
    print *, 'pool stub: no field ', trim(name)
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
  subroutine reg_i1(name, a)
    character(len=*), intent(in) :: name
    integer, target :: a(:)
    vi1(slot(ni1, name))%p => a
  end subroutine
  subroutine reg_i0(name, a)
    character(len=*), intent(in) :: name
    integer, target :: a
    vi0(slot(ni0, name))%p => a
  end subroutine
  subroutine reg_c0(name, a)
    character(len=*), intent(in) :: name
    character(len=StrKIND), target :: a
    vc0(slot(nc0, name))%p => a
  end subroutine
  subroutine reg_l0(name, a)
    character(len=*), intent(in) :: name
    logical, target :: a
    vl0(slot(nl0, name))%p => a
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
  subroutine get_i1(pool, name, p, tl)
    type(mpas_pool_type), intent(in) :: pool
    character(len=*), intent(in) :: name
    integer, pointer :: p(:)
    integer, intent(in), optional :: tl
    p => vi1(find(ni1, name))%p
  end subroutine
  subroutine get_i0(pool, name, p, tl)
    type(mpas_pool_type), intent(in) :: pool
    character(len=*), intent(in) :: name
    integer, pointer :: p
    integer, intent(in), optional :: tl
    p => vi0(find(ni0, name))%p
  end subroutine
  subroutine get_c0(pool, name, p, tl)
    type(mpas_pool_type), intent(in) :: pool
    character(len=*), intent(in) :: name
    character(len=StrKIND), pointer :: p
    integer, intent(in), optional :: tl
    p => vc0(find(nc0, name))%p
  end subroutine
  subroutine get_l0(pool, name, p)
    type(mpas_pool_type), intent(in) :: pool
    character(len=*), intent(in) :: name
    logical, pointer :: p
    p => vl0(find(nl0, name))%p
  end subroutine
  subroutine mpas_pool_get_dimension(pool, name, p)
    type(mpas_pool_type), intent(in) :: pool
    character(len=*), intent(in) :: name
    integer, pointer :: p
    p => vi0(find(ni0, name))%p
  end subroutine
end module mpas_pool_routines
