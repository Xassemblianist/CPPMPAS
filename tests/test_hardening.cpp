// SPDX-License-Identifier: GPL-3.0-or-later
// Direct calls that reach the situations of the GPU audit (MPAS leaves them undefined, or the port's capacity ends
// there), on the CPU build (test_hardening) and, the same source compiled by nvcc, on the GPU (test_gpu_hardening):
// the kernels below run as one-thread launches on the device there. Each check prints ok / FAIL lines.
//   views       the index arithmetic of YSU's, Kessler's and Noah's View3 / View2 at the audit's mesh sizes (offsets
//               above 2^31 - 1): every offset the exact 64-bit value, in a kernel
//   sfclayrev   the revised surface layer's psim / psih lookups: a z/L whose table index MPAS's int() leaves out of
//               the tables (s = 100 z/L <= -1, >= 2^31 - 1, NaN) is refused; the edges of the defined range run, the
//               kernel's results the host's bits
//   ntclamp     New Tiedtke's MIN / MAX as minsd / maxsd: cuadjtqn (kcall 2) on a -0.0 humidity at 8 K, where qsat
//               underflows to 0 and zcond = min(-0.0, 0.) is +0.0: pq stays -0.0 (CUDA's min.f64 gave -0.0 and pq
//               +0.0 before), and the clamp helpers on +-0 and NaN, kernel against host
//   raw         the column routines called directly with more levels than their arrays (and NTiedtke's, YSU's,
//               WSM6's, Kessler's, GWDO's lower bounds): refused before an array is touched
//   thin        the drivers on columns thinner than MPAS's code reads: driver_gwdo with 1 level, WSM6 and Kessler's
//               microphysics drivers with 1 / 2 levels, Noah with one soil layer (init_lsm, lsm, seaice_noah)
//   seaice      config_frac_seaice = false is refused by check_physics_supported, and driver_seaice refuses unset
//               open-water arrays
//   capacity    the column, scalar and RBF capacities refused on this build as on the other (GPU work arrays): 257
//               levels, 9 scalars, mpas_legs N = 19 and 17 RBF points (the latter two in a kernel)
//   radius      RRTMG's ice and snow sizes outside MPAS's ranges (module_ra_rrtmg_sw.F:2260, 2283, 2332,
//               module_ra_rrtmg_lw.F:2951, 2966, 2988) stop, as MPAS does; the range ends run
//   prep        the standalone RRTMG prep stages refuse more than 128 layers
//
// usage: test_hardening <check>...   (or all)    exit 0 ok, 1 failure, 77 (GPU build) no CUDA device
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "core/mesh.hpp"
#include "core_atmosphere/atm_advance_scalars.hpp"
#include "core_atmosphere/atm_microphysics.hpp"
#include "operators/vector_operations.hpp"
#include "physics/bl_gwdo.hpp"
#include "physics/bl_ysu.hpp"
#include "physics/cu_ntiedtke.hpp"
#include "physics/module_mp_wsm6.hpp"
#include "physics/module_sf_noahdrv.hpp"
#include "physics/mp_kessler.hpp"
#include "physics/mp_wsm6.hpp"
#include "physics/mpas_atmphys_constants.hpp"
#include "physics/mpas_atmphys_control.hpp"
#include "physics/mpas_atmphys_driver_lsm.hpp"
#include "physics/mpas_atmphys_gwdo.hpp"
#include "physics/mpas_atmphys_interface.hpp"
#include "physics/mpas_atmphys_microphysics.hpp"
#include "physics/rrtmg/rrtmg_lw_cldprmc.hpp"
#include "physics/rrtmg/rrtmg_lw_stages.hpp"
#include "physics/rrtmg/rrtmg_sw_cldprmc.hpp"
#include "physics/rrtmg/rrtmg_sw_stages.hpp"
#include "physics/sf_noah.hpp"
#include "physics/sf_sfclayrev.hpp"
#if defined(CPPMPAS_USE_CUDA)
#include "core/gpu_probe.hpp"
#endif

using namespace cppmpas;

namespace {

int failures = 0;
void check(bool ok, const std::string& what) {
  std::printf("%s %s\n", ok ? "ok  " : "FAIL", what.c_str());
  if (!ok) ++failures;
}
template <class F> std::string error_of(const F& f) {
  try {
    f();
  } catch (const Error& e) {
    return e.what();
  }
  return "";
}
void expect_error(const std::string& got, const std::string& want, const std::string& what) {
  const bool ok = got.find(want) != std::string::npos;
  check(ok, what + (ok ? "" : ": got '" + got + "', want '" + want + "'"));
}
std::uint64_t bits(double v) {
  std::uint64_t b;
  std::memcpy(&b, &v, 8);
  return b;
}
const char* kBuild =
#if defined(CPPMPAS_USE_CUDA)
    "GPU";
#else
    "CPU";
#endif

// ---- M1: offsets in std::ptrdiff_t ----
struct OffCase {
  int i, j, k, sj, sk;
};
void check_views() {
  // x1.65536002 at 55 levels (YSU's interface level k = 55, Kessler's k = 54; the audit's first overflowing index
  // i = 50'331'584 at k = 32), 835'586 cells at 100 levels, Noah's 8 soil levels of 268'435'457 points, a j plane
  constexpr int nc = 6;
  const OffCase cs[nc] = {{50331584, 0, 32, 0, 65536002}, {65536001, 0, 55, 0, 65536002}, {65536001, 0, 54, 0, 65536002},
                          {835585, 0, 100, 0, 835586},     {268435456, 0, 7, 0, 268435457}, {7, 3, 1000000, 40000000, 3000}};
  FArray<real> out({4 * nc}, 0.0);
  parallel_for(1, 1, [=, o = out.view()] CPPMPAS_LAMBDA(int) {
    for (int n = 0; n < nc; ++n) {
      const OffCase c = cs[n];
      const phys::ysu::View3<real> y{nullptr, c.sj, c.sk};
      const phys::kessler::View3<real> ke{nullptr, c.sj, c.sk};
      const phys::noah::View3<real> no{nullptr, c.sj, c.sk};
      const phys::noah::View2<real> n2{nullptr, c.sk};
      o(4 * n + 1) = real(y.offset(c.i, c.j, c.k));
      o(4 * n + 2) = real(ke.offset(c.i, c.j, c.k));
      o(4 * n + 3) = real(no.offset(c.i, c.j, c.k));
      o(4 * n + 4) = real(n2.offset(c.i, c.k));
    }
  });
  for (int n = 0; n < nc; ++n) {
    const OffCase c = cs[n];
    const long long want = (long long)c.i + (long long)c.j * c.sj + (long long)c.k * c.sk;
    const long long want2 = (long long)c.i + (long long)c.k * c.sk;
    const bool ok = out(4 * n + 1) == double(want) && out(4 * n + 2) == double(want) &&
                    out(4 * n + 3) == double(want) && out(4 * n + 4) == double(want2);
    check(ok, std::string("views: offset of (") + std::to_string(c.i) + ", " + std::to_string(c.j) + ", " +
                  std::to_string(c.k) + ") = " + std::to_string(want) + " in YSU's, Kessler's and Noah's views (" +
                  kBuild + ")");
  }
}

// ---- G3: the revised surface layer's lookups ----
CPPMPAS_HD real sfc_lookup(int f, real z, const phys::sfclayrev::SFClayRevTables<real>& t) {
  namespace s = phys::sfclayrev;
  return f == 0 ? s::psim_stable(z, t) : f == 1 ? s::psih_stable(z, t) : f == 2 ? s::psim_unstable(z, t)
                                                                                 : s::psih_unstable(z, t);
}
void check_sfclayrev() {
  auto tbl = std::make_unique<phys::sfclayrev::SFClayRevTables<real>>();
  phys::sfclayrev::init_sfclayrev_tables_o3(*tbl);
  phys::KernelCopy<phys::sfclayrev::SFClayRevTables<real>> keep;
  const phys::sfclayrev::SFClayRevTables<real>* tp = phys::kernel_ptr(*tbl, keep);
  const char* fn[4] = {"psim_stable", "psih_stable", "psim_unstable", "psih_unstable"};
  const real nan = std::numeric_limits<real>::quiet_NaN();
  // refused: the stable lookups of z/L <= -0.01 or >= 21474836.47 (int() of s = 100 z/L: a negative index, or
  // -2^31 from cvttsd2si), the unstable ones mirrored, and NaN
  struct Bad {
    int f;
    real z;
  };
  const Bad bad[] = {{0, 2.2e7}, {0, -0.02}, {1, -0.01}, {1, 21474836.47}, {2, 0.02}, {2, -2.2e7},
                     {3, 0.01},  {3, nan},   {0, nan},   {2, 1e300},       {1, -1e300}};
  for (const Bad& b : bad) {
    const int f = b.f;
    const real z = b.z;
    const std::string e = error_of([&] {
      parallel_for(1, 1, [=] CPPMPAS_LAMBDA(int) {
        const real r = sfc_lookup(f, z, *tp);
        (void)r;
      });
    });
    expect_error(e, "sf_sfclayrev: z/L outside the psim / psih lookup",
                 std::string("sfclayrev: ") + fn[f] + "(" + std::to_string(z) + ") refused (" + kBuild + ")");
  }
  // accepted: s just above -1 (index 0, extrapolated as MPAS does), 0, s = 999 (the full function), large z/L
  const real good[] = {-0.0099999, -0.0, 0.0, 0.0567, 9.99, 9.98, 2.0e7};
  constexpr int ng = sizeof(good) / sizeof(good[0]);
  FArray<real> out({4 * ng}, 0.0);
  for (int n = 0; n < ng; ++n)
    for (int f = 0; f < 4; ++f) {
      const real z = (f >= 2 ? -1.0 : 1.0) * good[n];
      const int slot = 4 * n + f + 1;
      const std::string e = error_of([&] {
        parallel_for(1, 1, [=, o = out.view()] CPPMPAS_LAMBDA(int) { o(slot) = sfc_lookup(f, z, *tp); });
      });
      const real host = sfc_lookup(f, z, *tbl);
      check(e.empty() && bits(out(slot)) == bits(host) && std::isfinite(host),
            std::string("sfclayrev: ") + fn[f] + "(" + std::to_string(z) + ") runs, the kernel's value the host's (" +
                kBuild + ")" + (e.empty() ? "" : ": " + e));
    }
}

// ---- G2: New Tiedtke's clamps ----
void check_ntclamp() {
  namespace nt = phys::cu_ntiedtke;
  namespace pc = physics_constants;
  const nt::NtConst<real> c = nt::make_ntconst<real>(pc::cp, pc::R_d, pc::rv, pc::xlv, pc::xls, pc::xlf, pc::gravity);
  FArray<real> out({12}, 0.0);
  auto run = [&](real* o) {
    real pt[nt::KM] = {}, pq[nt::KM] = {};
    pt[1] = 8.0;
    pq[1] = -0.0;
    nt::cuadjtqn(c, 1, real(5.e4), pt, pq, true, 2);
    o[0] = pt[1];
    o[1] = pq[1];
  };
  parallel_for(1, 1, [=, o = out.view()] CPPMPAS_LAMBDA(int) {
    real pt[nt::KM] = {}, pq[nt::KM] = {};
    pt[1] = 8.0;
    pq[1] = -0.0;
    nt::cuadjtqn(c, 1, real(5.e4), pt, pq, true, 2);
    o(1) = pt[1];
    o(2) = pq[1];
    const real nz = -0.0, pz = 0.0, nan = __builtin_nan("");
    o(3) = nt::cmin(nz, real(0.));  // MIN(-0.0, 0.): +0.0
    o(4) = nt::cmax(nz, real(0.));  // MAX(-0.0, 0.): +0.0
    o(5) = nt::vmin(nz, pz);        // MIN(a, b) = b < a ? b : a: a tie gives a
    o(6) = nt::vmax(pz, nz);
    o(7) = nt::cmin(nan, real(0.5));  // minsd: a NaN operand gives the source operand (here the constant)
    o(8) = nt::vmax(nan, real(1.0));  // MAX(a, b) = maxsd b, a: the source operand a
    o(9) = nt::vmax(real(1.0), nan);
  });
  real h[2];
  run(h);
  check(out(1) == 8.0 && bits(h[0]) == bits(8.0), std::string("ntclamp: cuadjtqn at 8 K leaves pt (") + kBuild + ")");
  check(bits(out(2)) == 0x8000000000000000ull && bits(h[1]) == 0x8000000000000000ull,
        std::string("ntclamp: cuadjtqn of pq = -0.0 at 8 K: zcond = min(-0.0, 0.) = +0.0, pq stays -0.0 (") + kBuild +
            ", kernel and host)");
  check(bits(out(3)) == 0 && bits(out(4)) == 0, "ntclamp: MIN / MAX(-0.0, 0.) are +0.0");
  check(bits(out(5)) == 0x8000000000000000ull && bits(out(6)) == 0,
        "ntclamp: MIN / MAX of two variables tie to the first argument");
  check(out(7) == 0.5 && std::isnan(out(8)) && out(9) == 1.0,
        "ntclamp: a NaN operand gives the source operand of minsd / maxsd");
}

// ---- M5: direct calls beyond the column arrays ----
void check_raw() {
  namespace nt = phys::cu_ntiedtke;
  for (int km : {101, 0}) {
    const std::string e = error_of([&] {
      parallel_for(1, 1, [=] CPPMPAS_LAMBDA(int) {
        const nt::NtConst<real> c{};
        real zprecc = 0;
        nt::cu_ntiedtke_run<real>(c, km, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                                  nullptr, nullptr, nullptr, nullptr, nullptr, 0., 0., zprecc, 1, 600., 1.e4,
                                  nt::ColumnSlots{});
      });
    });
    expect_error(e, "cu_ntiedtke_run: km outside 1..100",
                 "raw: cu_ntiedtke_run with km = " + std::to_string(km) + " refused (" + kBuild + ")");
  }
  {
    phys::ysu::YSUKernel k{};
    k.nz = 129;
    const std::string e = error_of([&] {
      parallel_for(1, 1, [=] CPPMPAS_LAMBDA(int) { k.compute<128>(0, 0); });
    });
    expect_error(e, "ysu: nz outside 1..MaxLevels", std::string("raw: YSU compute<128> with nz = 129 refused (") + kBuild + ")");
  }
  {
    const auto p = phys::module_mp_wsm6::mpas_wsm6_params();
    const std::string e = error_of([&] {
      parallel_for(1, 1, [=] CPPMPAS_LAMBDA(int) {
        real r = 0, rn = 0, s = 0, sn = 0, g = 0, gn = 0, sr = 0;
        p.operator()<128>(129, 60., nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                          nullptr, nullptr, r, rn, s, sn, g, gn, sr);
      });
    });
    expect_error(e, "mp_wsm6: nz outside 1..MaxLevels", std::string("raw: WSM6 operator()<128> with nz = 129 refused (") + kBuild + ")");
  }
  for (int nz : {129, 1}) {
    phys::kessler::KesslerKernel k{};
    k.nz = nz;
    const std::string e = error_of([&] { parallel_for(1, 1, [=] CPPMPAS_LAMBDA(int) { k(0, 0); }); });
    expect_error(e, "kessler: nz outside 2..128", "raw: Kessler with nz = " + std::to_string(nz) + " refused (" + kBuild + ")");
  }
  for (int kte : {129, 1}) {
    const std::string e = error_of([&] {
      parallel_for(1, 1, [=] CPPMPAS_LAMBDA(int) {
        phys::bl_gwdo::GwdoColumn<real> col{};
        phys::bl_gwdo::bl_gwdo_run(col, kte, 9.80616, 1004.5, 287., 0.608, 3.14159265358979, 600.);
      });
    });
    expect_error(e, "bl_gwdo_run: kte outside 2..128", "raw: bl_gwdo_run with kte = " + std::to_string(kte) + " refused (" + kBuild + ")");
  }
}

// ---- M4: thin columns through the drivers ----
void check_thin() {
  {
    phys::GwdoPools g;
    g.gwdo_scheme = "bl_ysu_gwdo";
    g.nVertLevels = 1;
    const phys::PhysicsVars v;
    expect_error(error_of([&] { phys::driver_gwdo(1, g, v, 1, 1); }), "driver_gwdo: fewer than 2 levels",
                 "thin: driver_gwdo with 1 level refused");
  }
  {
    phys::MicrophysicsPools m;
    m.index_qv = 1, m.index_qc = 2, m.index_qr = 3, m.index_qi = 4, m.index_qs = 5, m.index_qg = 6;
    m.nVertLevels = 1;
    m.mp_top_level = 1;
    expect_error(error_of([&] { phys::driver_microphysics_wsm6(phys::module_mp_wsm6::mpas_wsm6_params(), m, 1, 1); }),
                 "driver_microphysics: nVertLevels < 2", "thin: WSM6's driver_microphysics with 1 level refused");
  }
  {
    Mesh m;
    m.nVertLevels = 2;
    State s;
    Diag d;
    FArray<real> a, b, c, e;
    MicrophysicsFields mp;
    mp.active = true;
    expect_error(error_of([&] { driver_microphysics(m, s, a, b, c, d, e, mp); }),
                 "driver_microphysics: nVertLevels < 3", "thin: Kessler's driver_microphysics with 2 levels refused");
  }
  {
    phys::LsmConfig cfg;
    phys::LsmPools p;
    p.nSoilLevels = 1;
    phys::noah::NoahTables t;
    expect_error(error_of([&] { phys::init_lsm(cfg, p, "/nonexistent", t); }), "nSoilLevels = 1: the Noah port takes 2..",
                 "thin: init_lsm with one soil layer refused");
    phys::noahdrv::LsmArgs a;
    a.num_soil_layers = 1;
    expect_error(error_of([&] { phys::noahdrv::lsm(t, a); }), "num_soil_layers = 1 is outside 2..",
                 "thin: module_sf_noahdrv lsm with one soil layer refused");
    phys::noahdrv::SeaiceArgs sa;
    sa.num_soil_layers = 1;
    expect_error(error_of([&] { phys::noahdrv::seaice_noah(sa); }), "num_soil_layers = 1 is outside 2..",
                 "thin: seaice_noah with one soil layer refused");
  }
}

// ---- M2: fractional sea ice ----
void check_seaice() {
  phys::PhysicsConfig c;
  c.config_physics_suite = "mesoscale_reference";
  c.config_microp_scheme = "mp_wsm6";
  c.config_convection_scheme = "cu_ntiedtke";
  c.config_pbl_scheme = "bl_ysu";
  c.config_gwdo_scheme = "bl_ysu_gwdo";
  c.config_radt_lw_scheme = "rrtmg_lw";
  c.config_radt_sw_scheme = "rrtmg_sw";
  c.config_radt_cld_scheme = "cld_fraction";
  c.config_sfclayer_scheme = "sf_monin_obukhov_rev";
  c.config_lsm_scheme = "sf_noah";
  c.config_les_model = "none";
  c.config_bucket_update = c.config_conv_interval = c.config_pbl_interval = c.config_sounding_interval = "none";
  check(error_of([&] { phys::check_physics_supported(c); }).empty(),
        "seaice: check_physics_supported accepts mesoscale_reference with config_frac_seaice = true");
  c.config_frac_seaice = false;
  expect_error(error_of([&] { phys::check_physics_supported(c); }), "config_frac_seaice = false has no verified port",
               "seaice: config_frac_seaice = false refused");
  phys::LsmConfig cfg;
  phys::LsmPools p;
  p.nSoilLevels = 4;
  phys::LsmVars v;
  v.num_soils = 4;
  expect_error(error_of([&] { phys::driver_seaice(cfg, p, v, 1, 1); }),
               "driver_seaice: the surface layer's open-water arrays (*_sea) are not set",
               "seaice: driver_seaice without the open-water arrays refused");
}

// ---- G5: one range of inputs for both builds ----
void check_capacity() {
  expect_error(error_of([] { fits_column(257, "test_hardening"); }), "is more than the 256 levels",
               std::string("capacity: 257 levels refused (") + kBuild + ")");
  check(error_of([] { fits_column(256, "test_hardening"); }).empty(), std::string("capacity: 256 levels accepted (") + kBuild + ")");
  expect_error(error_of([] { check_scalar_capacity(9); }), "more than 8 scalars",
               std::string("capacity: 9 scalars refused (") + kBuild + ")");
  check(error_of([] { check_scalar_capacity(8); }).empty(), std::string("capacity: 8 scalars accepted (") + kBuild + ")");
  FArray<real> a({19 * 19}, 0.0), b({19}, 0.0), x({19}, 0.0);
  FArray<int> ix({19}, 0);
  const std::string legs = error_of([&] {
    parallel_for(1, 1, [=, A = a.view(), B = b.view(), X = x.view(), I = ix.view()] CPPMPAS_LAMBDA(int) {
      for (int i = 0; i < 19; ++i) {
        A[i + 19 * i] = 2.0;
        B[i] = 2.0;
      }
      mpas_legs(&A[0], 19, &B[0], &X[0], &I[0]);
    });
  });
  expect_error(legs, "mpas_legs: N > kMaxRbfPoints + 2", std::string("capacity: mpas_legs N = 19 refused in a kernel (") + kBuild + ")");
  FArray<real> src({17 * 3}, 0.0), uv({17 * 3}, 0.0), co({17 * 3}, 0.0);
  const std::string rbf = error_of([&] {
    parallel_for(1, 1, [=, S = src.view(), U = uv.view(), C = co.view()] CPPMPAS_LAMBDA(int) {
      const real dest[3] = {0.0, 0.0, 1.0};
      const real basis[6] = {1.0, 0.0, 0.0, 1.0, 0.0, 0.0};
      mpas_rbf_interp_func_3D_plane_vec_const_dir_comp_coeffs(17, &S[0], &U[0], dest, 1.0, basis, &C[0]);
    });
  });
  expect_error(rbf, "pointCount > kMaxRbfPoints", std::string("capacity: 17 RBF points refused in a kernel (") + kBuild + ")");
}

// ---- L1.1: RRTMG's cloud particle sizes ----
void check_radius() {
  namespace rr = phys::rrtmg;
  auto sw = std::make_unique<rr::RrtmgSwTables>();  // zero tables: the values do not matter, only the stops
  auto lw = std::make_unique<rr::RrtmgLwTables>();
  for (int ig = 0; ig < 112; ++ig) sw->ngb[ig] = 16 + ig % 14;  // the band of each g-point: 16..29 (SW), 1..16 (LW)
  for (int ig = 0; ig < 140; ++ig) lw->ngb[ig] = 1 + ig % 16;
  phys::KernelCopy<rr::RrtmgSwTables> ksw;
  phys::KernelCopy<rr::RrtmgLwTables> klw;
  const rr::RrtmgSwTables* swp = phys::kernel_ptr(*sw, ksw);
  const rr::RrtmgLwTables* lwp = phys::kernel_ptr(*lw, klw);
  struct Case {
    int iceflag;
    real radice, radsno;
    const char* want;  // "" runs
  };
  const Case cs[] = {{2, 134.0, 10.0, "ICE RADIUS OUT OF BOUNDS"},
                     {2, 4.9, 10.0, "ICE RADIUS OUT OF BOUNDS"},
                     {2, 131.0, 10.0, ""},
                     {2, 5.0, 10.0, ""},
                     {3, 140.5, 10.0, "ERROR: ICE GENERALIZED EFFECTIVE SIZE OUT OF BOUNDS"},
                     {3, 140.0, 10.0, ""},
                     {5, 20.0, 141.0, "ERROR: SNOW GENERALIZED EFFECTIVE SIZE OUT OF BOUNDS"},
                     {5, 20.0, 4.0, "ERROR: SNOW GENERALIZED EFFECTIVE SIZE OUT OF BOUNDS"},
                     {5, 20.0, 140.0, ""}};
  FArray<real> buf({8 * 140}, 0.0);  // LW: cldfmc, ciwp, clwp, cswp, reicmc, relqmc, resnmc (1 layer), taucmc
  for (const Case& c : cs) {
    const int iceflag = c.iceflag;
    const real radice = c.radice, radsno = c.radsno;
    const std::string what = " iceflag " + std::to_string(iceflag) + ", radice " + std::to_string(radice) +
                             ", radsno " + std::to_string(radsno) + " (" + kBuild + ")";
    const std::string esw = error_of([&] {
      parallel_for(1, 1, [=] CPPMPAS_LAMBDA(int) {
        rr::rb taorm = 0, taucm = 0, ssacm = 1, asmcm = 0;
        for (int ig = 0; ig < 112; ++ig)
          rr::cldprmc_sw_elem(2, iceflag, 1, ig, 1.0, 1.0, 0.0, iceflag == 5 ? 1.0 : 0.0, radice, 10.0, radsno, 0.0,
                              taorm, taucm, ssacm, asmcm, *swp);
      });
    });
    if (*c.want) expect_error(esw, c.want, std::string("radius: SW") + what + " stops");
    else check(esw.empty(), std::string("radius: SW") + what + " runs" + (esw.empty() ? "" : ": " + esw));
    // the inputs written on the host (the CSD views read through __ldg: nothing the kernel writes)
    for (int ig = 0; ig < 140; ++ig) {
      buf(ig + 1) = 1.0;            // cldfmc[ig][0]
      buf(7 * 140 + ig + 1) = 0.0;  // taucmc[ig][0]
    }
    buf(140 + 1) = 1.0;                        // ciwp
    buf(280 + 1) = 0.0;                        // clwp
    buf(420 + 1) = iceflag == 5 ? 1.0 : 0.0;   // cswp
    buf(560 + 1) = radice;
    buf(700 + 1) = 10.0;
    buf(840 + 1) = radsno;
    const std::string elw = error_of([&] {
      parallel_for(1, 1, [=, B = buf.view()] CPPMPAS_LAMBDA(int) {
        real* d = &B[0];
        rr::rrtmg_lw_cldprmc(1, 2, iceflag, 1, rr::CSD2{d, 1, 1}, rr::CSD{d + 140, 1}, rr::CSD{d + 280, 1},
                             rr::CSD{d + 420, 1}, rr::CSD{d + 560, 1}, rr::CSD{d + 700, 1}, rr::CSD{d + 840, 1},
                             rr::SD2{d + 7 * 140, 1, 1}, *lwp, 0, 1);
      });
    });
    if (*c.want) expect_error(elw, c.want, std::string("radius: LW") + what + " stops");
    else check(elw.empty(), std::string("radius: LW") + what + " runs" + (elw.empty() ? "" : ": " + elw));
  }
}

// ---- L1.2: the standalone prep stages ----
void check_prep() {
  namespace rr = phys::rrtmg;
  auto sw = std::make_unique<rr::RrtmgSwTables>();
  auto lw = std::make_unique<rr::RrtmgLwTables>();
  const rr::RrtmgLwConsts kc{};
  rr::RrtmgLwCallConfig lc{};
  lc.nlayers = 129;
  lc.num_levels = 128;
  expect_error(error_of([&] { rr::rrtmg_lw_prep(lc, rr::LwColumnIn{}, rr::LwViews{}, *lw, kc); }),
               "rrtmg_lw_prep: nlayers outside 1..128", "prep: rrtmg_lw_prep with 129 layers refused");
  rr::RrtmgSwCallConfig sc{};
  sc.nlayers = 129;
  sc.num_levels = 128;
  rr::SwColumnIn si{};
  si.coszen = 0.5;
  int iscal[4] = {};
  rr::SwViews sv{};
  sv.iscal = rr::SI{iscal, 1};
  expect_error(error_of([&] { rr::rrtmg_sw_prep(sc, si, sv, *sw, kc); }), "rrtmg_sw_prep: nlayers outside 1..128",
               "prep: rrtmg_sw_prep with 129 layers refused");
}

}  // namespace

int main(int argc, char** argv) {
#if defined(CPPMPAS_USE_CUDA)
  std::string why;
  if (!gpu::init(why)) {
    std::printf("SKIP: no usable CUDA device (%s)\n", why.c_str());
    return gpu::kNoDeviceExit;
  }
#endif
  if (argc < 2) {
    std::fprintf(stderr, "usage: test_hardening <check>... | all\n");
    return 2;
  }
  const char* all[] = {"views", "sfclayrev", "ntclamp", "raw", "thin", "seaice", "capacity", "radius", "prep"};
  std::vector<std::string> run(argv + 1, argv + argc);
  if (run.size() == 1 && run[0] == "all") run.assign(std::begin(all), std::end(all));
  for (const std::string& r : run) {
    if (r == "views") check_views();
    else if (r == "sfclayrev") check_sfclayrev();
    else if (r == "ntclamp") check_ntclamp();
    else if (r == "raw") check_raw();
    else if (r == "thin") check_thin();
    else if (r == "seaice") check_seaice();
    else if (r == "capacity") check_capacity();
    else if (r == "radius") check_radius();
    else if (r == "prep") check_prep();
    else {
      std::fprintf(stderr, "unknown check %s\n", r.c_str());
      return 2;
    }
  }
  std::printf("%s\n", failures ? "FAILED" : "all checks passed");
  return failures ? 1 : 0;
}
