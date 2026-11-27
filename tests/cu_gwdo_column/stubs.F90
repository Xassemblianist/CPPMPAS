! SPDX-License-Identifier: GPL-3.0-or-later
! Stand-ins for the MPAS framework modules that the MPAS physics sources of the convection / GWDO column check use (no
! arithmetic): tests/pbl_column/stubs.F90 (a pool that hands out pointers the harness registered by name, empty timers
! and log, the derived types todynamics names, the inner dimensions of mpas_atm_dimensions.F, an empty Noah-MP type)
! with the log levels mpas_atmphys_utilities.F names (MPAS_LOG_CRIT stops), integer namelist options taken from the
! registered dimensions (reg_i0; driver_convection reads config_gfconv_closure_*), and two modules
! mpas_atmphys_driver_gwdo.F uses only outside the bl_ysu_gwdo path: mpas_atmphys_manager (curr_julday, read by the
! bl_ugwp_gwdo call) and module_bl_ugwp_gwdo (gwdo_ugwp, never called: config_gwdo_scheme = 'bl_ugwp_gwdo' stops the
! harness).
module mpas_log
contains
  subroutine mpas_log_write(msg, intArgs, realArgs, logicArgs, messageType)
    character(len=*), intent(in) :: msg
    integer, intent(in), optional :: intArgs(:), messageType
    real(kind=8), intent(in), optional :: realArgs(:)
    logical, intent(in), optional :: logicArgs(:)
    print *, 'mpas_log_write: ', trim(msg)
    if (present(messageType)) then
      if (messageType == 4) stop 1  ! MPAS_LOG_CRIT aborts
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
  integer, parameter :: MPAS_LOG_ERR = 3, MPAS_LOG_CRIT = 4
  type domain_type
    integer :: dummy = 0
  end type
  type block_type
    type(domain_type), pointer :: domain => null()
  end type
end module mpas_derived_types

module mpas_dmpar
end module mpas_dmpar

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
    stop 'pool stub full'
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
    p => vi0(find(ni0, name))%p  ! integer options are registered with the dimensions (reg_i0)
  end subroutine
end module mpas_pool_routines

module mpas_atmphys_manager
  use mpas_kind_types
  real(kind=RKIND) :: curr_julday = 0.
end module mpas_atmphys_manager

module module_bl_ugwp_gwdo
  ! the interface of physics_wrf/module_bl_ugwp_gwdo.F:16-169 (its declarations copied), body replaced by a stop
  use mpas_kind_types, only: kind_phys => RKIND
  implicit none
contains
 subroutine gwdo_ugwp(u3d,v3d,t3d,qv3d,p3d,p3di,pi3d,z,pblh,kpbl2d,br1,xland, &
                 rublten,rvblten,rthblten,                                    &
                 dtaux3d,dtauy3d,dusfcg,dvsfcg,                               &
                 dtaux3d_ls,dtauy3d_ls,dtaux3d_bl,dtauy3d_bl,                 &
                 dtaux3d_ss,dtauy3d_ss,dtaux3d_fd,dtauy3d_fd,                 &
                 dusfc_ls,dvsfc_ls,dusfc_bl,dvsfc_bl,dusfc_ss,dvsfc_ss,       &
                 dusfc_fd,dvsfc_fd,ugwp_diags,ngw_scheme,xlatd,               &
                 jindx1_tau,jindx2_tau,ddy_j1tau,ddy_j2tau,r_DoY,             &
                 raincv,rainncv,ntau_d1y,ntau_d2t,days_limb,tau_limb,         &
                 dudt_ngw,dvdt_ngw,dtdt_ngw,                                  &
                 var2dls,oc12dls,oa2d1ls,oa2d2ls,oa2d3ls,oa2d4ls,ol2d1ls,     &
                 ol2d2ls,ol2d3ls,ol2d4ls,var2dss,oc12dss,oa2d1ss,oa2d2ss,     &
                 oa2d3ss,oa2d4ss,ol2d1ss,ol2d2ss,ol2d3ss,ol2d4ss,             &
                 sina,cosa,zi,dz,znu,znw,p_top,                               &
                 cp,g,rd,rv,ep1,pi,                                           &
                 dt,dx,itimestep,                                             &
                 ids,ide, jds,jde, kds,kde,                                   &
                 ims,ime, jms,jme, kms,kme,                                   &
                 its,ite, jts,jte, kts,kte,                                   &
                 errmsg,errflg                                                &
                )
!=================================================================================================================
!                                                                       
!-- u3d         3d u-velocity interpolated to theta points (m/s)
!-- v3d         3d v-velocity interpolated to theta points (m/s)
!-- t3d         temperature (k)
!-- qv3d        3d water vapor mixing ratio (kg/kg)
!-- p3d         3d pressure (pa)
!-- p3di        3d pressure (pa) at interface level
!-- pi3d        3d exner function (dimensionless)
!-- pblh        PBL height (m)
!-- br1         bulk Richardson number at lowest model level
!-- kpbl2d      index level of PBL top 
!-- xland       land mask (1 for land, 2 for water)
!-- rublten     u tendency due to pbl parameterization (m/s/s) 
!-- rvblten     v tendency due to pbl parameterization (m/s/s)
!-- rthblten    potential temperature tendency due to pbl parameterization (K/s)
!-- sina        sine rotation angle
!-- cosa        cosine rotation angle
!-- znu         eta values (sigma values)
!-- cp          heat capacity at constant pressure for dry air (j/kg/k)
!-- g           acceleration due to gravity (m/s^2)
!-- rd          gas constant for dry air (j/kg/k)
!-- z           height above sea level of layer centers (m)
!-- zi          height above sea level of layer interfaces (m)
!-- dz          layer thickness (m)
!-- rv          gas constant for water vapor (j/kg/k)
!-- dt          time step (s)
!-- dx          model grid interval (m)
!-- ep1         constant for virtual temperature (r_v/r_d - 1) (dimensionless)
!-- ids         start index for i in domain
!-- ide         end index for i in domain
!-- jds         start index for j in domain
!-- jde         end index for j in domain
!-- kds         start index for k in domain
!-- kde         end index for k in domain
!-- ims         start index for i in memory
!-- ime         end index for i in memory
!-- jms         start index for j in memory
!-- jme         end index for j in memory
!-- kms         start index for k in memory
!-- kme         end index for k in memory
!-- its         start index for i in tile
!-- ite         end index for i in tile
!-- jts         start index for j in tile
!-- jte         end index for j in tile
!-- kts         start index for k in tile
!-- kte         end index for k in tile
!
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: ids,ide,jds,jde,kds,kde, &
                      ims,ime,jms,jme,kms,kme, &
                      its,ite,jts,jte,kts,kte
 integer,intent(in):: itimestep

 integer,intent(in),dimension(ims:ime,jms:jme):: kpbl2d

 real(kind=kind_phys),intent(in):: dt,cp,g,rd,rv,ep1,pi
 real(kind=kind_phys),intent(in),optional:: p_top

 real(kind=kind_phys),intent(in),dimension(kms:kme),optional::        &
                                                                 znu, &
                                                                 znw

 real(kind=kind_phys),intent(in),dimension(ims:ime,jms:jme)::         &
                                                                  dx, &
                                                           sina,cosa, &
                                                 pblh,br1,xland

 real(kind=kind_phys),intent(in),dimension(ims:ime,jms:jme),optional::&
                            xlatd,raincv,rainncv,ddy_j1tau,ddy_j2tau
 integer,intent(in),dimension(ims:ime,jms:jme),optional::             &
                                               jindx1_tau,jindx2_tau
 integer,intent(in):: ntau_d1y,ntau_d2t
 real(kind=kind_phys),intent(in),dimension(ntau_d2t),optional::          &
                            days_limb
 real(kind=kind_phys),intent(in),dimension(ntau_d1y,ntau_d2t),optional:: &
                            tau_limb
 real(kind=kind_phys),intent(in)                           :: r_DoY

 real(kind=kind_phys),intent(in),dimension(ims:ime,jms:jme):: &
                                                     var2dls,var2dss, &
                                                     oc12dls,oc12dss, & 
                                     oa2d1ls,oa2d2ls,oa2d3ls,oa2d4ls, &
                                     oa2d1ss,oa2d2ss,oa2d3ss,oa2d4ss, &
                                     ol2d1ls,ol2d2ls,ol2d3ls,ol2d4ls, &
                                     ol2d1ss,ol2d2ss,ol2d3ss,ol2d4ss


 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme):: &
                                                                qv3d, &
                                                                 p3d, &
                                                                pi3d, &
                                                                 t3d, &
                                                                 u3d, &
                                                                 v3d, &
                                                                   z, &
                                                                  zi, &
                                                                  dz

 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme):: &
                                                                 p3di

 logical,intent(in)::                           ugwp_diags,ngw_scheme

!--- output arguments:
 character(len=*),intent(out)::  errmsg

 integer,intent(out):: errflg

 real(kind=kind_phys),intent(out),dimension(ims:ime,jms:jme)::    &
                                                   dusfcg,dvsfcg

 real(kind=kind_phys),intent(out),dimension(:,:),optional::       &
                             dusfc_ls,dusfc_bl,dusfc_ss,dusfc_fd, &
                             dvsfc_ls,dvsfc_bl,dvsfc_ss,dvsfc_fd

 real(kind=kind_phys),intent(out),dimension(ims:ime,kms:kme,jms:jme ):: &
                                                       dtaux3d,dtauy3d

 real(kind=kind_phys),intent(out),dimension(:,:,:),optional::           &
                           dtaux3d_ls,dtaux3d_bl,dtaux3d_ss,dtaux3d_fd, &
                           dtauy3d_ls,dtauy3d_bl,dtauy3d_ss,dtauy3d_fd

 real(kind=kind_phys),intent(out),dimension(:,:,:),optional::           &
                           dudt_ngw,dvdt_ngw,dtdt_ngw

!--- inout arguments:
 real(kind=kind_phys),intent(inout),dimension(ims:ime,kms:kme,jms:jme):: &
                           rublten,rvblten,rthblten

!--- local variables and arrays:
  stop 'gwdo_ugwp: not part of the bl_ysu_gwdo check'
  end subroutine gwdo_ugwp
end module module_bl_ugwp_gwdo
