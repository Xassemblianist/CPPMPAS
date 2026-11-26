! SPDX-License-Identifier: GPL-3.0-or-later
! Stand-ins for the MPAS framework modules that the extracted MPAS physics sources use (no arithmetic): a pool that
! hands out pointers the harness registered by name, empty timers and log, and an empty Noah-MP type.
module mpas_log
contains
  subroutine mpas_log_write(msg, intArgs, realArgs, logicArgs, messageType)
    character(len=*), intent(in) :: msg
    integer, intent(in), optional :: intArgs(:), messageType
    real(kind=8), intent(in), optional :: realArgs(:)
    logical, intent(in), optional :: logicArgs(:)
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

module mpas_pool_routines
  use mpas_kind_types
  use mpas_log
  implicit none
  type mpas_pool_type
    integer :: dummy = 0
  end type
  integer, parameter :: nmax = 128
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
  type i0p
    integer, pointer :: p => null()
  end type
  character(len=64) :: n1(nmax) = '', n2(nmax) = '', n3(nmax) = '', ni1(nmax) = '', ni0(nmax) = ''
  type(r1p) :: v1(nmax)
  type(r2p) :: v2(nmax)
  type(r3p) :: v3(nmax)
  type(i1p) :: vi1(nmax)
  type(i0p) :: vi0(nmax)
  character(len=StrKIND), target :: cfg_microp_scheme = 'mp_wsm6'
  real(kind=RKIND), target :: cfg_bucket_rainnc = 0.
  logical, target :: cfg_do_restart = .false.
  interface mpas_pool_get_array
    module procedure get_r1, get_r2, get_r3, get_i1, get_i0
  end interface
  interface mpas_pool_get_config
    module procedure cfg_c, cfg_r, cfg_l
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
  subroutine reg_i0(name, a)
    character(len=*), intent(in) :: name
    integer, target :: a
    vi0(slot(ni0, name))%p => a
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
    if (name /= 'config_microp_scheme') stop 'pool stub: config'
    p => cfg_microp_scheme
  end subroutine
  subroutine cfg_r(pool, name, p)
    type(mpas_pool_type), intent(in) :: pool
    character(len=*), intent(in) :: name
    real(kind=RKIND), pointer :: p
    if (name /= 'config_bucket_rainnc') stop 'pool stub: config'
    p => cfg_bucket_rainnc
  end subroutine
  subroutine cfg_l(pool, name, p)
    type(mpas_pool_type), intent(in) :: pool
    character(len=*), intent(in) :: name
    logical, pointer :: p
    if (name /= 'config_do_restart') stop 'pool stub: config'
    p => cfg_do_restart
  end subroutine
end module mpas_pool_routines
