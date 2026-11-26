// SPDX-License-Identifier: GPL-3.0-or-later
// Unit test of atm_compute_solve_diagnostics (mpas_atm_time_integration.F:6967-7437) and mpas_reconstruct_2d
// (mpas_vector_reconstruction.F:205-339) as the dycore will call them: on arrays that already hold values.
//  * Reference values: with every output poisoned (NaN, garbage columns included) beforehand, the two routines on
//    the init state give again the arrays atm_core_init gave (those are the reference's initdump, tolerance 0, in
//    the initdump_*_* tests), and leave the garbage columns alone (MPAS never writes them).
//  * Synthetic, non-degenerate data (u, rho_zz and fVertex perturbed): every output, after NaN poison, against an
//    independent straightforward implementation of the Fortran (ref_solve_diag / ref_reconstruct below), in a
//    sequence of calls as atm_srk3 makes them: rk_step absent, 1, 2, 3. At rk_step 1 and 2 v is not recomputed
//    (TI:7307-7310) and pv_edge uses the v it holds (TI:7428).
//  * mpas_reconstruct_2d's sincos (VR:310-313, two sincos@plt in the reference): latCell / lonCell set to values
//    where glibc's sincos and sin / cos differ by one ulp.
//  * Fortran's -a*b + c is -(a*b) + c, evaluated as c - a*b by the reference (no negation): a NaN of u keeps its sign
//    in uReconstructZonal / uReconstructMeridional.
//  * The refusals: rk_step outside 1..3, ranges other than one thread's, arrays of other extents, kiteForCell
//    outside 1..3, vertexDegree /= 3, config_apvm_upwinding <= 0 or NaN, nVertLevels < 1.
//
// usage: test_solve_diag {<init file> <namelist.atmosphere> <scalars>}...
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include "core/error.hpp"
#include "core/namelist.hpp"
#include "core_atmosphere/atm_init.hpp"
#include "core_atmosphere/atm_solve_diagnostics.hpp"
#include "io/mpas_reader.hpp"
#include "operators/vector_reconstruct.hpp"

extern "C" void sincos(double, double*, double*);

using namespace cppmpas;

// glibc's own entries, through pointers so that the compiler neither folds nor fuses them (sin and cos of one
// argument would become one sincos call, as gfortran makes in VR:310-313)
static double (*volatile c_sin)(double) = ::sin;
static double (*volatile c_cos)(double) = ::cos;
static void (*volatile c_sincos)(double, double*, double*) = ::sincos;

static int fails = 0;
static void expect(bool ok, const std::string& what) {
  if (!ok) {
    std::printf("FAIL %s\n", what.c_str());
    ++fails;
  }
}
template <class F> static bool throws_with(F f, const std::string& text) {
  try {
    f();
  } catch (const Error& e) {
    if (std::string(e.what()).find(text) != std::string::npos) return true;
    std::printf("  unexpected error: %s\n", e.what());
    return false;
  }
  return false;
}
static std::uint64_t bits(double x) {
  std::uint64_t b;
  std::memcpy(&b, &x, 8);
  return b;
}
static double from_bits(std::uint64_t b) {
  double x;
  std::memcpy(&x, &b, 8);
  return x;
}

// deterministic values in [-1, 1)
struct Lcg {
  std::uint64_t s;
  double next() {
    s = s * 6364136223846793005ull + 1442695040888963407ull;
    return double(s >> 11) * 0x1p-52 - 1.0;
  }
};

// the 15 arrays the two routines write
struct Out {
  FArray<real> v, rho_edge, pv_edge, gradPVt, gradPVn, vorticity, pv_vertex, divergence, ke, pv_cell, X, Y, Z, zonal,
      merid;
};
static std::vector<std::pair<const char*, FArray<real>*>> fields(Diag& d) {
  return {{"v", &d.v},
          {"rho_edge", &d.rho_edge},
          {"pv_edge", &d.pv_edge},
          {"gradPVt", &d.gradPVt},
          {"gradPVn", &d.gradPVn},
          {"vorticity", &d.vorticity},
          {"pv_vertex", &d.pv_vertex},
          {"divergence", &d.divergence},
          {"ke", &d.ke},
          {"pv_cell", &d.pv_cell},
          {"uReconstructX", &d.uReconstructX},
          {"uReconstructY", &d.uReconstructY},
          {"uReconstructZ", &d.uReconstructZ},
          {"uReconstructZonal", &d.uReconstructZonal},
          {"uReconstructMeridional", &d.uReconstructMeridional}};
}
static std::vector<FArray<real>*> fields(Out& o) {
  return {&o.v,  &o.rho_edge,  &o.pv_edge,    &o.gradPVt, &o.gradPVn, &o.vorticity, &o.pv_vertex, &o.divergence,
          &o.ke, &o.pv_cell,   &o.X,          &o.Y,       &o.Z,       &o.zonal,     &o.merid};
}

// ---------------------------------------------------------------------------------------------------------------
// Independent implementation of the Fortran, written from mpas_atm_time_integration.F:7136-7430 directly: each
// accumulation in a local variable that starts at 0, one statement per Fortran statement.
static void ref_solve_diag(const Mesh& m, real dt, real apvm, const FArray<real>& h, const FArray<real>& u,
                           const FArray<real>& v_prior, bool reconstruct_v, Out& o) {
  const int nL = m.nVertLevels, nC = m.nCells, nE = m.nEdges, nV = m.nVertices;
  for (auto* a : {&o.v, &o.rho_edge, &o.pv_edge, &o.gradPVt, &o.gradPVn}) a->allocate({nL, nE + 1}, 0.0);
  for (auto* a : {&o.vorticity, &o.pv_vertex}) a->allocate({nL, nV + 1}, 0.0);
  for (auto* a : {&o.divergence, &o.ke, &o.pv_cell}) a->allocate({nL, nC + 1}, 0.0);
  FArray<real> ke_edge({nL, nE + 1}, 0.0), ke_vertex({nL, nV + 1}, 0.0);
  for (int e = 1; e <= nE; ++e)
    for (int k = 1; k <= nL; ++k) {
      o.rho_edge(k, e) = 0.5 * (h(k, m.cellsOnEdge(1, e)) + h(k, m.cellsOnEdge(2, e)));            // :7141
      ke_edge(k, e) = (m.dcEdge(e) * m.dvEdge(e)) * (u(k, e) * u(k, e));                          // :7147, :7151
    }
  for (int iv = 1; iv <= nV; ++iv)
    for (int k = 1; k <= nL; ++k) {
      real acc = 0.0;  // :7163-7176
      for (int i = 1; i <= m.vertexDegree; ++i) {
        const int e = m.edgesOnVertex(i, iv);
        acc = acc + (m.edgesOnVertex_sign(i, iv) * m.dcEdge(e)) * u(k, e);
      }
      o.vorticity(k, iv) = acc * m.invAreaTriangle(iv);  // :7179
    }
  for (int c = 1; c <= nC; ++c)
    for (int k = 1; k <= nL; ++k) {
      real dv = 0.0, kc = 0.0;  // :7190-7203, :7223-7236
      for (int i = 1; i <= m.nEdgesOnCell(c); ++i) {
        const int e = m.edgesOnCell(i, c);
        dv = dv + (m.edgesOnCell_sign(i, c) * m.dvEdge(e)) * u(k, e);
        kc = kc + 0.25 * ke_edge(k, e);
      }
      o.divergence(k, c) = dv * m.invAreaCell(c);  // :7206
      o.ke(k, c) = kc * m.invAreaCell(c);          // :7238
    }
  for (int iv = 1; iv <= nV; ++iv)  // :7263
    for (int k = 1; k <= nL; ++k)
      ke_vertex(k, iv) =
          (ke_edge(k, m.edgesOnVertex(1, iv)) + ke_edge(k, m.edgesOnVertex(2, iv)) + ke_edge(k, m.edgesOnVertex(3, iv))) *
          (0.25 * m.invAreaTriangle(iv));
  for (int c = 1; c <= nC; ++c)
    for (int k = 1; k <= nL; ++k) {
      real kc = 0.625 * o.ke(k, c);  // :7280, ke_fact = 1.0 - .375
      for (int i = 1; i <= m.nEdgesOnCell(c); ++i)  // :7287-7296
        kc = kc + 0.375 * m.kiteAreasOnVertex(m.kiteForCell(i, c), m.verticesOnCell(i, c)) *
                      ke_vertex(k, m.verticesOnCell(i, c)) * m.invAreaCell(c);
      o.ke(k, c) = kc;
    }
  for (int e = 1; e <= nE; ++e)
    for (int k = 1; k <= nL; ++k) {
      if (!reconstruct_v) {
        o.v(k, e) = v_prior(k, e);
        continue;
      }
      real acc = 0.0;  // :7318-7327
      for (int i = 1; i <= m.nEdgesOnEdge(e); ++i) acc = acc + m.weightsOnEdge(i, e) * u(k, m.edgesOnEdge(i, e));
      o.v(k, e) = acc;
    }
  for (int iv = 1; iv <= nV; ++iv)
    for (int k = 1; k <= nL; ++k) o.pv_vertex(k, iv) = m.fVertex(iv) + o.vorticity(k, iv);  // :7351
  for (int c = 1; c <= nC; ++c)
    for (int k = 1; k <= nL; ++k) {
      real acc = 0.0;  // :7384-7393
      for (int i = 1; i <= m.nEdgesOnCell(c); ++i)
        acc = acc + m.kiteAreasOnVertex(m.kiteForCell(i, c), m.verticesOnCell(i, c)) *
                        o.pv_vertex(k, m.verticesOnCell(i, c)) * m.invAreaCell(c);
      o.pv_cell(k, c) = acc;
    }
  for (int e = 1; e <= nE; ++e)
    for (int k = 1; k <= nL; ++k) {
      const real p1 = o.pv_vertex(k, m.verticesOnEdge(1, e)), p2 = o.pv_vertex(k, m.verticesOnEdge(2, e));
      const real pe = 0.5 * (p1 + p2);                                                                     // :7367
      o.gradPVt(k, e) = (p2 - p1) * (1.0 * m.invDvEdge(e));                                               // :7426
      o.gradPVn(k, e) = (o.pv_cell(k, m.cellsOnEdge(2, e)) - o.pv_cell(k, m.cellsOnEdge(1, e))) * (1.0 * m.invDcEdge(e));
      o.pv_edge(k, e) = pe - (apvm * dt) * (o.v(k, e) * o.gradPVt(k, e) + u(k, e) * o.gradPVn(k, e));    // :7428
    }
}

// mpas_vector_reconstruction.F:274-334 with glibc's own sincos, the evaluation of the reference executable
// (c - a*b, as at 0x589bf0 / 0x589c18; equal to -(a*b) + c for non-NaN values)
static void ref_reconstruct(const Mesh& m, const FArray<real>& u, Out& o) {
  const int nL = m.nVertLevels, nC = m.nCells;
  for (auto* a : {&o.X, &o.Y, &o.Z, &o.zonal, &o.merid}) a->allocate({nL, nC + 1}, 0.0);
  for (int c = 1; c <= nC; ++c) {
    double slat, clat, slon, clon;
    c_sincos(m.latCell(c), &slat, &clat);
    c_sincos(m.lonCell(c), &slon, &clon);
    for (int k = 1; k <= nL; ++k) {
      real x = 0.0, y = 0.0, z = 0.0;
      for (int i = 1; i <= m.nEdgesOnCell(c); ++i) {
        const real ue = u(k, m.edgesOnCell(i, c));
        x = x + m.coeffs_reconstruct(1, i, c) * ue;
        y = y + m.coeffs_reconstruct(2, i, c) * ue;
        z = z + m.coeffs_reconstruct(3, i, c) * ue;
      }
      o.X(k, c) = x, o.Y(k, c) = y, o.Z(k, c) = z;
      o.zonal(k, c) = m.on_a_sphere ? y * clon - x * slon : x;
      o.merid(k, c) = m.on_a_sphere ? z * clat - (x * clon + y * slon) * slat : y;
    }
  }
}

// every element of columns 1..n of each output equal, bit for bit; the garbage columns (n+1) untouched
static void compare(Diag& d, Out& o, const Diag* before, const std::string& what) {
  auto df = fields(d);
  auto of = fields(o);
  for (std::size_t f = 0; f < df.size(); ++f) {
    const FArray<real>& a = *df[f].second;
    const FArray<real>& b = *of[f];
    const int nL = a.extent(0), n = a.extent(1) - 1;
    long bad = 0;
    for (int j = 1; j <= n; ++j)
      for (int k = 1; k <= nL; ++k) bad += bits(a(k, j)) != bits(b(k, j));
    expect(bad == 0, what + ": " + df[f].first + " differs in " + std::to_string(bad) + " elements");
    if (before) {
      const FArray<real>& p = *fields(const_cast<Diag&>(*before))[f].second;
      long g = 0;
      for (int k = 1; k <= nL; ++k) g += bits(a(k, n + 1)) != bits(p(k, n + 1));
      expect(g == 0, what + ": " + df[f].first + " garbage column written");
    }
  }
}

static void poison(Diag& d, std::uint64_t pattern) {
  for (auto& [n, a] : fields(d)) {
    (void)n;
    const double x = from_bits(pattern);
    for (std::size_t i = 0; i < a->size(); ++i) a->data()[i] = x;
  }
}
static void poison_work(SolveDiagWork& w, std::uint64_t pattern) {  // 1..n; MPAS zeroes column n+1 (CORE:492-499)
  for (auto* a : {&w.ke_vertex, &w.ke_edge})
    for (int j = 1; j < a->extent(1); ++j)
      for (int k = 1; k <= a->extent(0); ++k) (*a)(k, j) = from_bits(pattern);
}

static void run_both(const Mesh& m, const InitConfig& cfg, const FArray<real>& h, const FArray<real>& u, Diag& d,
                     SolveDiagWork& w, std::optional<int> rk) {
  atm_compute_solve_diagnostics(cfg.config_dt, cfg.config_apvm_upwinding, m, h, u, d, w, 1, m.nCells, 1, m.nVertices,
                                1, m.nEdges, rk);
  mpas_reconstruct_2d(m, u, d.uReconstructX, d.uReconstructY, d.uReconstructZ, d.uReconstructZonal,
                      d.uReconstructMeridional);
}
static void ref_both(const Mesh& m, const InitConfig& cfg, const FArray<real>& h, const FArray<real>& u,
                     const FArray<real>& v_prior, bool reconstruct_v, Out& o) {
  ref_solve_diag(m, cfg.config_dt, cfg.config_apvm_upwinding, h, u, v_prior, reconstruct_v, o);
  ref_reconstruct(m, u, o);
}

static void perturb(const FArray<real>& a, FArray<real>& out, Lcg& g, real abs_scale, real rel_scale) {
  out = a;
  const int nL = a.extent(0), n = a.extent(1) - 1;
  for (int j = 1; j <= n; ++j)
    for (int k = 1; k <= nL; ++k) out(k, j) = a(k, j) * (1.0 + rel_scale * g.next()) + abs_scale * g.next();
}

static void one_case(const std::string& init, const std::string& nml, const std::string& scalars_csv) {
  std::vector<std::string> scalars;
  for (std::size_t p = 0; p <= scalars_csv.size();) {
    const std::size_t q = std::min(scalars_csv.find(',', p), scalars_csv.size());
    scalars.push_back(scalars_csv.substr(p, q - p));
    p = q + 1;
  }
  MpasInput in = read_mpas_init(init, scalars);
  InitConfig cfg = read_init_config(nml);
  std::vector<std::string> log;
  atm_core_init(in.mesh, in.state, in.diag, cfg, log);
  Mesh& m = in.mesh;
  const std::string tag = init.substr(init.find_last_of('/') + 1);
  const int nL = m.nVertLevels;

  // --- allocation: MPAS's extents, the garbage column of the work arrays zero (CORE:492-499)
  SolveDiagWork w;
  allocate_solve_diag_work(m, w);
  expect(w.ke_vertex.extents() == std::vector<int>{nL, m.nVertices + 1} &&
             w.ke_edge.extents() == std::vector<int>{nL, m.nEdges + 1},
         tag + ": ke_vertex / ke_edge are (nVertLevels, n+1)");
  {
    Diag e;
    allocate_solve_diag_fields(m, e);
    bool ok = true;
    for (auto& [n, a] : fields(e)) {
      const std::string s = n;
      const int n1 = (s == "v" || s == "rho_edge" || s == "pv_edge" || s == "gradPVt" || s == "gradPVn") ? m.nEdges + 1
                     : (s == "vorticity" || s == "pv_vertex")                                            ? m.nVertices + 1
                                                                                                          : m.nCells + 1;
      ok = ok && a->extents() == std::vector<int>{nL, n1};
    }
    expect(ok, tag + ": allocate_solve_diag_fields gives MPAS's extents");
  }

  // --- reference values: the init state again, on poisoned arrays
  const Diag d0 = in.diag;
  {
    Out o;
    ref_both(m, cfg, in.state.rho_zz, in.state.u, d0.v, true, o);
    Diag dd = d0;
    compare(dd, o, nullptr, tag + " independent implementation on the init state = atm_core_init (initdump)");
  }
  for (std::uint64_t pat : {0x7ff8000000000000ull, 0xfff80000000abcdeull, 0x8000000000000000ull}) {
    Diag d = in.diag;
    poison(d, pat);
    poison_work(w, pat);
    const Diag before = d;
    run_both(m, cfg, in.state.rho_zz, in.state.u, d, w, std::nullopt);
    Out o;
    ref_both(m, cfg, in.state.rho_zz, in.state.u, d0.v, true, o);
    char p[32];
    std::snprintf(p, sizeof p, "%016llx", static_cast<unsigned long long>(pat));
    compare(d, o, &before, tag + " init state on arrays poisoned with " + p);
  }

  // --- synthetic data, the atm_srk3 sequence (rk_step absent, then 1, 2, 3), each call on the previous outputs
  Lcg g{0x9e3779b97f4a7c15ull ^ std::uint64_t(m.nCells)};
  for (int iv = 1; iv <= m.nVertices; ++iv) m.fVertex(iv) = m.fVertex(iv) + 1.0e-4 * g.next();
  Diag d = in.diag;
  poison(d, 0x7ff8000000000001ull);
  poison_work(w, 0x7ff8000000000001ull);
  FArray<real> h, u;
  for (std::optional<int> rk : {std::optional<int>{}, std::optional<int>{1}, std::optional<int>{2},
                                std::optional<int>{3}, std::optional<int>{1}}) {
    perturb(in.state.rho_zz, h, g, 0.0, 0.05);
    perturb(in.state.u, u, g, 5.0, 0.2);
    const Diag before = d;
    const bool recon = !rk || *rk == 3;
    run_both(m, cfg, h, u, d, w, rk);
    Out o;
    ref_both(m, cfg, h, u, before.v, recon, o);
    const std::string what = tag + " synthetic, rk_step " + (rk ? std::to_string(*rk) : std::string("absent"));
    compare(d, o, &before, what);
    if (!recon) {
      long changed = 0;
      for (std::size_t i = 0; i < d.v.size(); ++i) changed += bits(d.v.data()[i]) != bits(before.v.data()[i]);
      expect(changed == 0, what + ": v unchanged");
    }
  }

  // --- the sincos of VR:310-313. glibc 2.41's sincos and sin differ only in the sine, on 0.855469 <= |x| < 2.426265
  // (s_sincos.c:66-72 against s_sin.c: do_cos(a, da) with a + da = (hp0 - |x|) + hp1 split, against
  // do_cos(hp0 - |x|, hp1)); the cosines are the same function. Of the candidates y = hp0 - |x| = +-(2m+1)/256
  // (+-4 ulps) only x = hp0 + 27/256 = 1.6762650767948966 differs, a longitude; no latitude (|x| < pi/2) differs
  // there, nor in 2e7 random latitudes of the strip, so the latitude's sincos is shown by the object code alone.
  // lonCell of cell 2 is set to that value (and latCell of cell 1, should a latitude be found).
  if (m.on_a_sphere) {
    const double hp0 = 0x1.921fb54442d18p0, inf = INFINITY;
    double xlat = 0.0, xlon = 0.0;
    for (int mm = 0; mm <= 109; ++mm)
      for (double x0 : {hp0 - (2 * mm + 1) / 256.0, hp0 + (2 * mm + 1) / 256.0}) {
        if (x0 < 0.855469 || x0 >= 2.426265) continue;
        double x = x0;
        for (int k = 0; k < 4; ++k) x = std::nextafter(x, 0.0);
        for (int k = -4; k <= 4; ++k, x = std::nextafter(x, inf)) {
          double s, c;
          c_sincos(x, &s, &c);
          if (bits(s) == bits(c_sin(x))) continue;
          if (x < hp0 && xlat == 0.0) xlat = x;
          if (x > hp0 && xlon == 0.0) xlon = x;
        }
      }
    expect(xlon != 0.0, tag + ": a longitude where sincos's sine and sin differ");
    if (xlat != 0.0) m.latCell(1) = xlat;
    m.lonCell(2) = xlon;
    std::printf("%s: sincos strip: latCell(1) = %.17g, lonCell(2) = %.17g\n", tag.c_str(), m.latCell(1), xlon);
    Diag e = d;
    poison(e, 0x7ff8000000000002ull);
    mpas_reconstruct_2d(m, u, e.uReconstructX, e.uReconstructY, e.uReconstructZ, e.uReconstructZonal,
                        e.uReconstructMeridional);
    Out o;
    ref_reconstruct(m, u, o);
    long bad = 0, sens = 0;
    for (int c : {1, 2})
      for (int k = 1; k <= nL; ++k) {
        bad += bits(e.uReconstructZonal(k, c)) != bits(o.zonal(k, c)) ||
               bits(e.uReconstructMeridional(k, c)) != bits(o.merid(k, c));
        // the same with sin / cos: the strip values must change the result, or this check could not fail
        const double slat = c_sin(m.latCell(c)), clat = c_cos(m.latCell(c));
        const double slon = c_sin(m.lonCell(c)), clon = c_cos(m.lonCell(c));
        const double x = o.X(k, c), y = o.Y(k, c), z = o.Z(k, c);
        sens += bits(y * clon - x * slon) != bits(o.zonal(k, c)) ||
                bits(z * clat - (x * clon + y * slon) * slat) != bits(o.merid(k, c));
      }
    expect(bad == 0, tag + ": reconstruct on the sincos strip");
    expect(sens > 0, tag + ": the sincos strip changes the reconstruction (sin / cos would give other bits)");
  }

  // --- NaN signs: -(a*b) + c is c - a*b in the reference, so a NaN of u reaches Zonal / Meridional unnegated
  {
    FArray<real> un = u;
    const int e1 = m.edgesOnCell(1, 3), e2 = m.edgesOnCell(1, 7);
    const std::uint64_t q1 = 0x7ff8000000000123ull, q2 = 0xfff8000000000456ull;
    un(1, e1) = from_bits(q1);
    un(nL, e2) = from_bits(q2);
    Diag e = d;
    mpas_reconstruct_2d(m, un, e.uReconstructX, e.uReconstructY, e.uReconstructZ, e.uReconstructZonal,
                        e.uReconstructMeridional);
    bool ok = true;
    for (auto [c, k, q] : {std::tuple{3, 1, q1}, std::tuple{7, nL, q2}})
      for (FArray<real>* a : {&e.uReconstructX, &e.uReconstructY, &e.uReconstructZ, &e.uReconstructZonal,
                              &e.uReconstructMeridional})
        ok = ok && bits((*a)(k, c)) == q;
    expect(ok, tag + ": a NaN of u keeps its sign and payload in uReconstruct* (no negation)");
  }

  // --- refusals, each before anything is written
  auto call = [&](Diag& dd, SolveDiagWork& ww, int cs, int ce, int vs, int ve, int es, int ee, std::optional<int> rk,
                  real apvm) {
    atm_compute_solve_diagnostics(cfg.config_dt, apvm, m, in.state.rho_zz, in.state.u, dd, ww, cs, ce, vs, ve, es, ee,
                                  rk);
  };
  const int nC = m.nCells, nV = m.nVertices, nE = m.nEdges;
  const real a = cfg.config_apvm_upwinding;
  {
    Diag dd = d;
    for (int rk : {0, 4, -1})
      expect(throws_with([&] { call(dd, w, 1, nC, 1, nV, 1, nE, rk, a); }, "rk_step = " + std::to_string(rk)),
             tag + ": rk_step " + std::to_string(rk) + " refused");
    const char* ranges = "ranges must be 1..nCells";
    expect(throws_with([&] { call(dd, w, 1, nC - 1, 1, nV, 1, nE, {}, a); }, ranges), tag + ": cell subrange refused");
    expect(throws_with([&] { call(dd, w, 2, nC, 1, nV, 1, nE, {}, a); }, ranges), tag + ": cellStart 2 refused");
    expect(throws_with([&] { call(dd, w, 0, nC, 1, nV, 1, nE, {}, a); }, ranges), tag + ": cellStart 0 refused");
    expect(throws_with([&] { call(dd, w, 1, nC, 1, nV + 1, 1, nE, {}, a); }, ranges), tag + ": vertexEnd n+1 refused");
    expect(throws_with([&] { call(dd, w, 1, nC, 0, nV, 1, nE, {}, a); }, ranges), tag + ": vertexStart 0 refused");
    expect(throws_with([&] { call(dd, w, 1, nC, 1, nV, 1, nE + 1, {}, a); }, ranges), tag + ": edgeEnd n+1 refused");
    expect(throws_with([&] { call(dd, w, 1, nC, 1, nV, 2, nE, {}, a); }, ranges), tag + ": edgeStart 2 refused");
    for (real bad : {0.0, -0.5, std::nan("")})
      expect(throws_with([&] { call(dd, w, 1, nC, 1, nV, 1, nE, {}, bad); }, "config_apvm_upwinding <= 0"),
             tag + ": config_apvm_upwinding " + std::to_string(bad) + " refused");
    SolveDiagWork none;
    expect(throws_with([&] { call(dd, none, 1, nC, 1, nV, 1, nE, {}, a); }, "ke_edge is not allocated"),
           tag + ": unallocated work arrays refused");
    Diag dv = d;
    dv.pv_cell.allocate({nL, nC}, 0.0);
    expect(throws_with([&] { call(dv, w, 1, nC, 1, nV, 1, nE, {}, a); }, "pv_cell is not allocated"),
           tag + ": pv_cell without its garbage column refused");
    const int kf = m.kiteForCell(2, 5);
    m.kiteForCell(2, 5) = 0;
    expect(throws_with([&] { call(dd, w, 1, nC, 1, nV, 1, nE, {}, a); }, "kiteForCell(2, 5) = 0"),
           tag + ": kiteForCell 0 refused");
    m.kiteForCell(2, 5) = 4;
    expect(throws_with([&] { call(dd, w, 1, nC, 1, nV, 1, nE, {}, a); }, "kiteForCell(2, 5) = 4"),
           tag + ": kiteForCell 4 refused");
    m.kiteForCell(2, 5) = kf;
    m.vertexDegree = 4;
    expect(throws_with([&] { call(dd, w, 1, nC, 1, nV, 1, nE, {}, a); }, "vertexDegree /= 3"),
           tag + ": vertexDegree 4 refused");
    m.vertexDegree = 3;
    long changed = 0;
    for (std::size_t f = 0; f < fields(dd).size(); ++f) {
      const FArray<real>& x = *fields(dd)[f].second;
      const FArray<real>& y = *fields(d)[f].second;
      for (std::size_t i = 0; i < x.size(); ++i) changed += bits(x.data()[i]) != bits(y.data()[i]);
    }
    expect(changed == 0, tag + ": a refused call writes nothing");
    Diag dr = d;
    dr.uReconstructMeridional.allocate({nL, nC}, 0.0);
    expect(throws_with([&] { mpas_reconstruct_2d(m, u, dr.uReconstructX, dr.uReconstructY, dr.uReconstructZ,
                                                 dr.uReconstructZonal, dr.uReconstructMeridional); },
                       "uReconstructMeridional is not allocated"),
           tag + ": reconstruct refuses an output without its garbage column");
    long rchanged = 0;
    for (std::size_t i = 0; i < dr.uReconstructX.size(); ++i)
      rchanged += bits(dr.uReconstructX.data()[i]) != bits(d.uReconstructX.data()[i]);
    expect(rchanged == 0, tag + ": a refused reconstruct writes nothing");
    Mesh nocoef = m;
    nocoef.coeffs_reconstruct = FArray<real>();
    expect(throws_with([&] { mpas_reconstruct_2d(nocoef, u, dr.uReconstructX, dr.uReconstructY, dr.uReconstructZ,
                                                 dr.uReconstructZonal, dr.uReconstructZ); },
                       "coeffs_reconstruct is not"),
           tag + ": reconstruct refuses a mesh without coeffs_reconstruct");
    Mesh flat = m;
    flat.nVertLevels = 0;
    SolveDiagWork w0;
    Diag d0f;
    expect(throws_with([&] { allocate_solve_diag_work(flat, w0); }, "nVertLevels < 1"),
           tag + ": allocate_solve_diag_work refuses nVertLevels 0");
    expect(throws_with([&] { allocate_solve_diag_fields(flat, d0f); }, "nVertLevels < 1"),
           tag + ": allocate_solve_diag_fields refuses nVertLevels 0");
  }
}

int main(int argc, char** argv) {
  if (argc < 4 || (argc - 1) % 3 != 0) {
    std::printf("usage: test_solve_diag {<init file> <namelist.atmosphere> <scalars>}...\n");
    return 2;
  }
  try {
    for (int i = 1; i + 2 < argc; i += 3) one_case(argv[i], argv[i + 1], argv[i + 2]);
  } catch (const Error& e) {
    std::printf("FAIL unexpected error: %s\n", e.what());
    return 1;
  }
  if (fails == 0) std::printf("OK solve diagnostics / reconstruct unit test\n");
  return fails ? 1 : 0;
}
