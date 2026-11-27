// SPDX-License-Identifier: GPL-3.0-or-later
// See mp_thompson_init.hpp. MT = src/core_atmosphere/physics/physics_wrf/module_mp_thompson.F (MPAS v8.4.2).
#include "physics/mp_thompson_init.hpp"

#include <cmath>
#include <limits>
#include <cstdio>
#include <cstring>
#include <fstream>

#include "core/error.hpp"
#include "core/libm.hpp"
#include "physics/mpas_atmphys_functions.hpp"

namespace cppmpas::phys::thompson {

namespace {

using detail::c_max;
using detail::c_min;

// r_c (MT:274-280): the cloud water contents of the tables (kg/m**3)
constexpr double r_c[ntb_c] = {1.e-6, 2.e-6, 3.e-6, 4.e-6, 5.e-6, 6.e-6, 7.e-6, 8.e-6, 9.e-6, 1.e-5, 2.e-5, 3.e-5, 4.e-5,
                               5.e-5, 6.e-5, 7.e-5, 8.e-5, 9.e-5, 1.e-4, 2.e-4, 3.e-4, 4.e-4, 5.e-4, 6.e-4, 7.e-4, 8.e-4,
                               9.e-4, 1.e-3, 2.e-3, 3.e-3, 4.e-3, 5.e-3, 6.e-3, 7.e-3, 8.e-3, 9.e-3, 1.e-2};

// A power of constants gfortran folds at compile time (MPFR, correctly rounded): glibc's pow gives the same
// double for every one of them (checked against the reference's thompson_init, tests/thompson_column)
double folded_pow(double a, double b) { return ::cppmpas::libm::pow(a, b); }
double folded_log(double a) { return ::cppmpas::libm::log(a); }

// MT:680-736 bin edges: xDx(n) = DEXP(REAL(n-1)/REAL(nb) * DLOG(xDx(nb+1)/xDx(1)) + DLOG(xDx(1))), n = 2..nb.
// gfortran folds the logs and the whole terms n = 2, nb - 1 and nb at compile time (MPFR, correctly rounded: the
// peeled first iteration and the two of the epilogue) and runs n = 3 .. nb - 2 through libmvec's 2-lane exp (2
// _ZGVbN2v_exp per loop). The folded exps: expl rounded to double (CPPWRF 7fac5bf src/phys/mp_thompson/
// thompson_tables.cpp:47-53, the same rule in WRF's object), checked against the reference's xDx.
// gfortran folds the edges n = 2, nb-1 and nb with MPFR (correctly rounded); x87 expl rounded to double gives the same
// 15 values (checked by thompson_column's init records). This needs the x86-64 80-bit long double.
static_assert(std::numeric_limits<long double>::digits == 64, "bin_edges: the folded edges need x87 long double");
void bin_edges(double* xDx, int nb) {
  const double l1 = folded_log(xDx[nb] / xDx[0]);
  const double l0 = folded_log(xDx[0]);
  for (int n = 2; n <= nb; ++n) {
    const double a = double(n - 1) / double(nb) * l1 + l0;
    xDx[n - 1] = (n == 2 || n >= nb - 1) ? double(expl(static_cast<long double>(a))) : ::cppmpas::libm::mvec_exp(a);
  }
}

// table_Efrw (MT:4068-4124)
void table_Efrw(const ThompsonView& v, double* t_Efrw) {
  for (int j = 1; j <= nbc; ++j) {
    for (int i = 1; i <= nbr; ++i) {
      double Ef_rw = 0.0;
      const double Dcj = v.Dc[j - 1], Dri = v.Dr[i - 1];
      const double p = Dcj / Dri;
      double& t = t_Efrw[(i - 1) + std::size_t(nbr) * (j - 1)];
      if (Dri < 50.E-6 || Dcj < 3.E-6) {
        t = 0.0;
      } else if (p > 0.25) {
        const double X = Dcj * 1.E6;
        if (Dri < 75.e-6) {
          Ef_rw = 0.026794 * X - 0.20604;
        } else if (Dri < 125.e-6) {
          Ef_rw = -0.00066842 * X * X + 0.061542 * X - 0.37089;
        } else if (Dri < 175.e-6) {
          Ef_rw = 4.091e-06 * X * X * X * X - 0.00030908 * X * X * X + 0.0066237 * X * X - 0.0013687 * X - 0.073022;
        } else if (Dri < 250.e-6) {
          Ef_rw = 9.6719e-5 * X * X * X - 0.0068901 * X * X + 0.17305 * X - 0.65988;
        } else if (Dri < 350.e-6) {
          Ef_rw = 9.0488e-5 * X * X * X - 0.006585 * X * X + 0.16606 * X - 0.56125;
        } else {
          Ef_rw = 0.00010721 * X * X * X - 0.0072962 * X * X + 0.1704 * X - 0.46929;
        }
      } else {
        const double vtr = -0.1021 + 4.932E3 * Dri - 0.9551E6 * Dri * Dri + 0.07934E9 * Dri * Dri * Dri -
                           0.002362E12 * Dri * Dri * Dri * Dri;
        const double stokes = Dcj * Dcj * vtr * rho_w / (9. * 1.718E-5 * Dri);
        const double reynolds = 9. * stokes / (p * p * rho_w);
        const double F = ::cppmpas::libm::log(reynolds);
        const double G = -0.1007 - 0.358 * F + 0.0261 * F * F;
        const double K0 = ::cppmpas::libm::exp(G);
        const double z = ::cppmpas::libm::log(stokes / (K0 + 1.E-15));
        const double H = 0.1465 + 1.302 * z - 0.607 * z * z + 0.293 * z * z * z;
        const double yc0 = 2.0 / PI * ::cppmpas::libm::atan(H);
        Ef_rw = (yc0 + p) * (yc0 + p) / ((1. + p) * (1. + p));
      }
      t = c_max(0.0, c_min(Ef_rw, 0.95));
    }
  }
}

// table_Efsw (MT:4132-4168)
void table_Efsw(const ThompsonView& v, double* t_Efsw) {
  for (int j = 1; j <= nbc; ++j) {
    const double Dcj = v.Dc[j - 1];
    const double vtc = 1.19E4 * (1.0E4 * Dcj * Dcj * 0.25);
    for (int i = 1; i <= nbs; ++i) {
      const double Dsi = v.Ds[i - 1];
      double& t = t_Efsw[(i - 1) + std::size_t(nbs) * (j - 1)];
      const double vts = av_s * ::cppmpas::libm::pow(Dsi, bv_s) * ::cppmpas::libm::exp(-fv_s * Dsi) - vtc;
      const double Ds_m = ::cppmpas::libm::pow(am_s * (Dsi * Dsi) / am_r, v.obmr);  // Ds(i)**bm_s, bm_s = 2.0
      const double p = Dcj / Ds_m;
      if (p > 0.25 || Dsi < D0s || Dcj < 6.E-6 || vts < 1.E-3) {
        t = 0.0;
      } else {
        const double stokes = Dcj * Dcj * vts * rho_w / (9. * 1.718E-5 * Ds_m);
        const double reynolds = 9. * stokes / (p * p * rho_w);
        const double F = ::cppmpas::libm::log(reynolds);
        const double G = -0.1007 - 0.358 * F + 0.0261 * F * F;
        const double K0 = ::cppmpas::libm::exp(G);
        const double z = ::cppmpas::libm::log(stokes / (K0 + 1.E-15));
        const double H = 0.1465 + 1.302 * z - 0.607 * z * z + 0.293 * z * z * z;
        const double yc0 = 2.0 / PI * ::cppmpas::libm::atan(H);
        const double Ef_sw = (yc0 + p) * (yc0 + p) / ((1. + p) * (1. + p));
        t = c_max(0.0, c_min(Ef_sw, 0.95));
      }
    }
  }
}

// table_dropEvap (MT:4222-4310)
void table_dropEvap(const ThompsonView& v, double* tpc_wev, double* tnc_wev) {
  double N_c[nbc], massc[nbc];
  for (int n = 1; n <= nbc; ++n)  // MT:4234-4236, _ZGVbN2vv_pow in pairs (100 values)
    massc[n - 1] = am_r * detail::vpow(v.Dc[n - 1], bm_r, n - 1, nbc);

  for (int k = 1; k <= nbc; ++k) {
    const int nu_c = detail::i_min(15, detail::nint(1000.E6 / v.t_Nc[k - 1]) + 2);
    for (int j = 1; j <= ntb_c; ++j) {
      const double lamc =
          ::cppmpas::libm::pow(v.t_Nc[k - 1] * am_r * v.CCG(2, nu_c) * v.ocg1[nu_c - 1] / r_c[j - 1], v.obmr);
      const double N0_c = v.t_Nc[k - 1] * v.ocg1[nu_c - 1] * ::cppmpas::libm::pow(lamc, v.CCE(1, nu_c));
      for (int i = 1; i <= nbc; ++i) {
        N_c[i - 1] = N0_c * detail::powi(v.Dc[i - 1], nu_c) * ::cppmpas::libm::exp(-lamc * v.Dc[i - 1]) * v.dtc[i - 1];
        double summ = 0., summ2 = 0.;
        for (int n = 1; n <= i; ++n) {
          summ = summ + massc[n - 1] * N_c[n - 1];
          summ2 = summ2 + N_c[n - 1];
        }
        const std::size_t ix = (i - 1) + std::size_t(nbc) * ((j - 1) + std::size_t(ntb_c) * (k - 1));
        tpc_wev[ix] = summ;
        tnc_wev[ix] = summ2;
      }
    }
  }
}

// One table file of MT:850-965: big-endian Fortran unformatted sequential records, one per array
struct DblFile {
  std::ifstream f;
  std::string name;
  DblFile(const std::string& dir, const std::string& file, const char* fatal) : f(dir + "/" + file, std::ios::binary), name(file) {
    if (!f) throw Error(std::string("subroutine thompson_init: ") + fatal);  // physics_error_fatal
  }
  static std::uint32_t be32(const unsigned char* b) {
    return (std::uint32_t(b[0]) << 24) | (std::uint32_t(b[1]) << 16) | (std::uint32_t(b[2]) << 8) | b[3];
  }
  void read(FStorage<double>& a, std::size_t n) {
    a.assign(n, 0.0);
    unsigned char m[4];
    if (!f.read(reinterpret_cast<char*>(m), 4) || be32(m) != n * 8)
      throw Error("thompson_init: " + name + ": record length is not " + std::to_string(n * 8) + " bytes");
    std::vector<unsigned char> b(n * 8);
    if (!f.read(reinterpret_cast<char*>(b.data()), std::streamsize(b.size())))
      throw Error("thompson_init: " + name + ": short record");
    for (std::size_t i = 0; i < n; ++i) {
      std::uint64_t u = 0;
      for (int c = 0; c < 8; ++c) u = (u << 8) | b[i * 8 + c];
      std::memcpy(&a[i], &u, 8);
    }
    if (!f.read(reinterpret_cast<char*>(m), 4) || be32(m) != n * 8)
      throw Error("thompson_init: " + name + ": record end marker");
  }
};

}  // namespace

void ThompsonTables::point() {
  v.tcg_racg = tcg_racg.data(); v.tmr_racg = tmr_racg.data(); v.tcr_gacr = tcr_gacr.data();
  v.tmg_gacr = tmg_gacr.data(); v.tnr_racg = tnr_racg.data(); v.tnr_gacr = tnr_gacr.data();
  v.tcs_racs1 = tcs_racs1.data(); v.tmr_racs1 = tmr_racs1.data(); v.tcs_racs2 = tcs_racs2.data();
  v.tmr_racs2 = tmr_racs2.data(); v.tcr_sacr1 = tcr_sacr1.data(); v.tms_sacr1 = tms_sacr1.data();
  v.tcr_sacr2 = tcr_sacr2.data(); v.tms_sacr2 = tms_sacr2.data(); v.tnr_racs1 = tnr_racs1.data();
  v.tnr_racs2 = tnr_racs2.data(); v.tnr_sacr1 = tnr_sacr1.data(); v.tnr_sacr2 = tnr_sacr2.data();
  v.tpi_qrfz = tpi_qrfz.data(); v.tni_qrfz = tni_qrfz.data(); v.tpg_qrfz = tpg_qrfz.data();
  v.tnr_qrfz = tnr_qrfz.data(); v.tpi_qcfz = tpi_qcfz.data(); v.tni_qcfz = tni_qcfz.data();
  v.tpi_ide = tpi_ide.data(); v.tps_iaus = tps_iaus.data(); v.tni_iaus = tni_iaus.data();
  v.t_Efrw = t_Efrw.data(); v.t_Efsw = t_Efsw.data();
}

void thompson_init(const std::string& dir, ThompsonTables& t) {
  ThompsonView& v = t.v;
  v = ThompsonView{};

  v.Sc3 = folded_pow(Sc, 1. / 3.);  // MT:501-511
  v.D0i = folded_pow(xm0i / am_i, 1. / bm_i);
  v.xm0s = am_s * (D0s * D0s);  // D0s**bm_s, bm_s = 2.0
  v.xm0g = am_g * folded_pow(D0g, bm_g);

  for (int n = 1; n <= 15; ++n) {  // MT:515-527
    double* c = &v.cce[(n - 1) * 5];
    double* g = &v.ccg[(n - 1) * 5];
    c[0] = n + 1.;
    c[1] = bm_r + n + 1.;
    c[2] = bm_r + n + 4.;
    c[3] = n + bv_c + 1.;
    c[4] = bm_r + n + bv_c + 1.;
    for (int i = 0; i < 5; ++i) g[i] = wgamma(c[i]);
    v.ocg1[n - 1] = 1. / g[0];
    v.ocg2[n - 1] = 1. / g[1];
  }

  v.cie[0] = mu_i + 1.;  // MT:529-546
  v.cie[1] = bm_i + mu_i + 1.;
  v.cie[2] = bm_i + mu_i + bv_i + 1.;
  v.cie[3] = mu_i + bv_i + 1.;
  v.cie[4] = mu_i + 2.;
  v.cie[5] = bm_i * 0.5 + mu_i + bv_i + 1.;
  v.cie[6] = bm_i * 0.5 + mu_i + 1.;
  for (int n = 0; n < 7; ++n) v.cig[n] = wgamma(v.cie[n]);
  v.oig1 = 1. / v.cig[0];
  v.oig2 = 1. / v.cig[1];
  v.obmi = 1. / bm_i;

  v.cre[0] = bm_r + 1.;  // MT:548-568
  v.cre[1] = mu_r + 1.;
  v.cre[2] = bm_r + mu_r + 1.;
  v.cre[3] = bm_r * 2. + mu_r + 1.;
  v.cre[4] = mu_r + bv_r + 1.;
  v.cre[5] = bm_r + mu_r + bv_r + 1.;
  v.cre[6] = bm_r * 0.5 + mu_r + bv_r + 1.;
  v.cre[7] = bm_r + mu_r + bv_r + 3.;
  v.cre[8] = mu_r + bv_r + 3.;
  v.cre[9] = mu_r + 2.;
  v.cre[10] = 0.5 * (bv_r + 5. + 2. * mu_r);
  v.cre[11] = bm_r * 0.5 + mu_r + 1.;
  v.cre[12] = bm_r * 2. + mu_r + bv_r + 1.;
  for (int n = 0; n < 13; ++n) v.crg[n] = wgamma(v.cre[n]);
  v.obmr = 1. / bm_r;
  v.ore1 = 1. / v.cre[0];
  v.org1 = 1. / v.crg[0];
  v.org2 = 1. / v.crg[1];
  v.org3 = 1. / v.crg[2];

  v.cse[0] = bm_s + 1.;  // MT:570-593
  v.cse[1] = bm_s + 2.;
  v.cse[2] = bm_s * 2.;
  v.cse[3] = bm_s + bv_s + 1.;
  v.cse[4] = bm_s * 2. + bv_s + 1.;
  v.cse[5] = bm_s * 2. + 1.;
  v.cse[6] = bm_s + mu_s + 1.;
  v.cse[7] = bm_s + mu_s + 2.;
  v.cse[8] = bm_s + mu_s + 3.;
  v.cse[9] = bm_s + mu_s + bv_s + 1.;
  v.cse[10] = bm_s * 2. + mu_s + bv_s + 1.;
  v.cse[11] = bm_s * 2. + mu_s + 1.;
  v.cse[12] = bv_s + 2.;
  v.cse[13] = bm_s + bv_s;
  v.cse[14] = mu_s + 1.;
  v.cse[15] = 1.0 + (1.0 + bv_s) / 2.;
  v.cse[16] = v.cse[15] + mu_s + 1.;
  v.cse[17] = bv_s + mu_s + 3.;
  for (int n = 0; n < 18; ++n) v.csg[n] = wgamma(v.cse[n]);
  v.oams = 1. / am_s;
  v.obms = 1. / bm_s;
  v.ocms = folded_pow(v.oams, v.obms);

  v.cge[0] = bm_g + 1.;  // MT:595-617
  v.cge[1] = mu_g + 1.;
  v.cge[2] = bm_g + mu_g + 1.;
  v.cge[3] = bm_g * 2. + mu_g + 1.;
  v.cge[4] = bm_g * 2. + mu_g + bv_g + 1.;
  v.cge[5] = bm_g + mu_g + bv_g + 1.;
  v.cge[6] = bm_g + mu_g + bv_g + 2.;
  v.cge[7] = bm_g + mu_g + bv_g + 3.;
  v.cge[8] = mu_g + bv_g + 3.;
  v.cge[9] = mu_g + 2.;
  v.cge[10] = 0.5 * (bv_g + 5. + 2. * mu_g);
  v.cge[11] = 0.5 * (bv_g + 5.) + mu_g;
  for (int n = 0; n < 12; ++n) v.cgg[n] = wgamma(v.cge[n]);
  v.oamg = 1. / am_g;
  v.obmg = 1. / bm_g;
  v.ocmg = folded_pow(v.oamg, v.obmg);
  v.oge1 = 1. / v.cge[0];
  v.ogg1 = 1. / v.cgg[0];
  v.ogg2 = 1. / v.cgg[1];
  v.ogg3 = 1. / v.cgg[2];

  // MT:623-657: rate-equation constants
  v.t1_qr_qc = PI * .25 * av_r * v.crg[8];
  v.t1_qr_qi = PI * .25 * av_r * v.crg[8];
  v.t2_qr_qi = PI * .25 * am_r * av_r * v.crg[7];
  v.t1_qg_qc = PI * .25 * av_g * v.cgg[8];
  v.t1_qs_qc = PI * .25 * av_s;
  v.t1_qs_qi = PI * .25 * av_s;
  v.t1_qr_ev = 0.78 * v.crg[9];
  v.t2_qr_ev = 0.308 * v.Sc3 * std::sqrt(av_r) * v.crg[10];
  v.t1_qs_sd = 0.86;
  v.t2_qs_sd = 0.28 * v.Sc3 * std::sqrt(av_s);
  v.t1_qs_me = PI * 4. * C_sqrd * olfus * 0.86;
  v.t2_qs_me = PI * 4. * C_sqrd * olfus * 0.28 * v.Sc3 * std::sqrt(av_s);
  v.t1_qg_sd = 0.86 * v.cgg[9];
  v.t2_qg_sd = 0.28 * v.Sc3 * std::sqrt(av_g) * v.cgg[10];
  v.t1_qg_me = PI * 4. * C_cube * olfus * 0.86 * v.cgg[9];
  v.t2_qg_me = PI * 4. * C_cube * olfus * 0.28 * v.Sc3 * std::sqrt(av_g) * v.cgg[10];

  // MT:660-668: table index offsets, folded NINT(ALOG10(..)) of the first table values
  auto n10 = [](double x) { return detail::nint(::cppmpas::libm::log10(x)); };
  v.nic2 = n10(r_c1);
  v.nii2 = n10(r_i1);
  v.nii3 = n10(Nt_i1);
  v.nir2 = n10(r_r1);
  v.nir3 = n10(1.e6);  // N0r_exp(1)
  v.nis2 = n10(r_s1);
  v.nig2 = n10(r_g1);
  v.nig3 = n10(1.e4);  // N0g_exp(1)
  v.niIN2 = n10(Nt_IN1);

  // Bins (MT:671-745)
  v.Dc[0] = D0c * 1.0;
  v.dtc[0] = D0c * 1.0;
  for (int n = 2; n <= nbc; ++n) {
    v.Dc[n - 1] = v.Dc[n - 2] + 1.0E-6;
    v.dtc[n - 1] = (v.Dc[n - 1] - v.Dc[n - 2]);
  }
  double* xDx = t.xDx;
  auto edges = [&](double x1, double xtop, int nb, double* D, double* dt) {
    xDx[0] = x1;
    xDx[nb] = xtop;
    bin_edges(xDx, nb);
    for (int n = 1; n <= nb; ++n) {
      D[n - 1] = std::sqrt(xDx[n - 1] * xDx[n]);
      dt[n - 1] = xDx[n] - xDx[n - 1];
    }
  };
  edges(v.D0i * 1.0, 5.0 * D0s, nbi, v.Di, v.dti);
  edges(D0r * 1.0, 0.005, nbr, v.Dr, v.dtr);
  edges(D0s * 1.0, 0.02, nbs, v.Ds, v.dts);
  edges(D0g * 1.0, 0.05, nbg, v.Dg, v.dtg);
  xDx[0] = 1.0;  // cloud droplet number (MT:737-745)
  xDx[nbc] = 3000.0;
  bin_edges(xDx, nbc);
  for (int n = 1; n <= nbc; ++n) v.t_Nc[n - 1] = std::sqrt(xDx[n - 1] * xDx[n]) * 1.E6;
  v.nic1 = static_cast<int>(::cppmpas::libm::log(v.t_Nc[nbc - 1] / v.t_Nc[0]));

  // Tables (MT:751-843 set them to 0; l_mp_tables = .true.)
  t.t_Efrw.assign(std::size_t(nbr) * nbc, 0.0);
  t.t_Efsw.assign(std::size_t(nbs) * nbc, 0.0);
  t.tpc_wev.assign(std::size_t(nbc) * ntb_c * nbc, 0.0);
  t.tnc_wev.assign(std::size_t(nbc) * ntb_c * nbc, 0.0);
  table_Efrw(v, t.t_Efrw.data());  // MT:849-850
  table_Efsw(v, t.t_Efsw.data());
  table_dropEvap(v, t.tpc_wev.data(), t.tnc_wev.data());  // MT:854

  const std::size_t n_gr = std::size_t(ntb_g1) * ntb_g * ntb_r1 * ntb_r;
  const std::size_t n_sr = std::size_t(ntb_s) * ntb_t * ntb_r1 * ntb_r;
  const std::size_t n_rfz = std::size_t(ntb_r) * ntb_r1 * 45 * ntb_IN;
  const std::size_t n_cfz = std::size_t(ntb_c) * nbc * 45 * ntb_IN;
  const std::size_t n_ia = std::size_t(ntb_i) * ntb_i1;
  {  // MT:862-874
    DblFile f(dir, "MP_THOMPSON_QRacrQG_DATA.DBL", "failure opening MP_THOMPSON_QRacrQG.DBL");
    f.read(t.tcg_racg, n_gr); f.read(t.tmr_racg, n_gr); f.read(t.tcr_gacr, n_gr);
    f.read(t.tmg_gacr, n_gr); f.read(t.tnr_racg, n_gr); f.read(t.tnr_gacr, n_gr);
  }
  {  // MT:890-907
    DblFile f(dir, "MP_THOMPSON_QRacrQS_DATA.DBL", "failure opening MP_THOMPSON_QRacrQS.DBL");
    f.read(t.tcs_racs1, n_sr); f.read(t.tmr_racs1, n_sr); f.read(t.tcs_racs2, n_sr); f.read(t.tmr_racs2, n_sr);
    f.read(t.tcr_sacr1, n_sr); f.read(t.tms_sacr1, n_sr); f.read(t.tcr_sacr2, n_sr); f.read(t.tms_sacr2, n_sr);
    f.read(t.tnr_racs1, n_sr); f.read(t.tnr_racs2, n_sr); f.read(t.tnr_sacr1, n_sr); f.read(t.tnr_sacr2, n_sr);
  }
  {  // MT:933-944
    DblFile f(dir, "MP_THOMPSON_freezeH2O_DATA.DBL", "failure opening MP_THOMPSON_freezeH2O.DBL");
    f.read(t.tpi_qrfz, n_rfz); f.read(t.tni_qrfz, n_rfz); f.read(t.tpg_qrfz, n_rfz); f.read(t.tnr_qrfz, n_rfz);
    f.read(t.tpi_qcfz, n_cfz); f.read(t.tni_qcfz, n_cfz);
  }
  {  // MT:960-968
    DblFile f(dir, "MP_THOMPSON_QIautQS_DATA.DBL", "failure opening MP_THOMPSON_QIautQS.DBL");
    f.read(t.tpi_ide, n_ia); f.read(t.tps_iaus, n_ia); f.read(t.tni_iaus, n_ia);
  }

  // MT:978-988: radar_init with xam_r = am_r, xbm_r = bm_r, xmu_r = mu_r, ... of Thompson
  v.radar = radar::radar_init(am_r, bm_r, mu_r, am_s, bm_s, mu_s, am_g, bm_g, mu_g);
  t.point();
}

// The dump of the column check (tests/thompson_column/dump.c writes the same records from MPAS's module variables):
// per entry a 32-byte name, an int64 count and count doubles; an int64 count of -1 and two doubles (count, FNV-1a
// hash of the bytes) for the tables over 100000 values.
void dump_thompson_init(const ThompsonTables& t, const std::string& path) {
  std::FILE* f = std::fopen(path.c_str(), "wb");
  if (!f) throw Error("dump_thompson_init: cannot open " + path);
  auto put = [&](const char* name, const double* a, std::size_t n) {
    char nm[32] = {};
    std::strncpy(nm, name, 31);
    std::fwrite(nm, 1, 32, f);
    if (n > 100000) {
      std::uint64_t h = 1469598103934665603ULL;
      const unsigned char* b = reinterpret_cast<const unsigned char*>(a);
      for (std::size_t i = 0; i < n * 8; ++i) h = (h ^ b[i]) * 1099511628211ULL;
      const std::int64_t m = -1;
      std::fwrite(&m, 8, 1, f);
      const double cnt = double(n);
      double hd;
      std::memcpy(&hd, &h, 8);
      std::fwrite(&cnt, 8, 1, f);
      std::fwrite(&hd, 8, 1, f);
    } else {
      const std::int64_t m = std::int64_t(n);
      std::fwrite(&m, 8, 1, f);
      std::fwrite(a, 8, n, f);
    }
  };
  auto puti = [&](const char* name, int i) {
    const double d = i;
    put(name, &d, 1);
  };
  const ThompsonView& v = t.v;
  put("cce", v.cce, 75); put("ccg", v.ccg, 75); put("ocg1", v.ocg1, 15); put("ocg2", v.ocg2, 15);
  put("cie", v.cie, 7); put("cig", v.cig, 7); put("oig1", &v.oig1, 1); put("oig2", &v.oig2, 1);
  put("obmi", &v.obmi, 1);
  put("cre", v.cre, 13); put("crg", v.crg, 13); put("ore1", &v.ore1, 1); put("org1", &v.org1, 1);
  put("org2", &v.org2, 1); put("org3", &v.org3, 1); put("obmr", &v.obmr, 1);
  put("cse", v.cse, 18); put("csg", v.csg, 18); put("oams", &v.oams, 1); put("obms", &v.obms, 1);
  put("ocms", &v.ocms, 1);
  put("cge", v.cge, 12); put("cgg", v.cgg, 12); put("oge1", &v.oge1, 1); put("ogg1", &v.ogg1, 1);
  put("ogg2", &v.ogg2, 1); put("ogg3", &v.ogg3, 1); put("oamg", &v.oamg, 1); put("obmg", &v.obmg, 1);
  put("ocmg", &v.ocmg, 1);
  put("t1_qr_qc", &v.t1_qr_qc, 1); put("t1_qr_qi", &v.t1_qr_qi, 1); put("t2_qr_qi", &v.t2_qr_qi, 1);
  put("t1_qg_qc", &v.t1_qg_qc, 1); put("t1_qs_qc", &v.t1_qs_qc, 1); put("t1_qs_qi", &v.t1_qs_qi, 1);
  put("t1_qr_ev", &v.t1_qr_ev, 1); put("t2_qr_ev", &v.t2_qr_ev, 1); put("t1_qs_sd", &v.t1_qs_sd, 1);
  put("t2_qs_sd", &v.t2_qs_sd, 1); put("t1_qg_sd", &v.t1_qg_sd, 1); put("t2_qg_sd", &v.t2_qg_sd, 1);
  put("t1_qs_me", &v.t1_qs_me, 1); put("t2_qs_me", &v.t2_qs_me, 1); put("t1_qg_me", &v.t1_qg_me, 1);
  put("t2_qg_me", &v.t2_qg_me, 1);
  put("sc3", &v.Sc3, 1); put("d0i", &v.D0i, 1); put("xm0s", &v.xm0s, 1); put("xm0g", &v.xm0g, 1);
  puti("nic1", v.nic1); puti("nic2", v.nic2); puti("nii2", v.nii2); puti("nii3", v.nii3); puti("nir2", v.nir2);
  puti("nir3", v.nir3); puti("nis2", v.nis2); puti("nig2", v.nig2); puti("nig3", v.nig3); puti("niin2", v.niIN2);
  put("dc", v.Dc, nbc); put("dtc", v.dtc, nbc); put("di", v.Di, nbi); put("dti", v.dti, nbi);
  put("dr", v.Dr, nbr); put("dtr", v.dtr, nbr); put("ds", v.Ds, nbs); put("dts", v.dts, nbs);
  put("dg", v.Dg, nbg); put("dtg", v.dtg, nbg); put("t_nc", v.t_Nc, nbc); put("xdx", t.xDx, nbins + 1);
  put("t_efrw", t.t_Efrw.data(), t.t_Efrw.size()); put("t_efsw", t.t_Efsw.data(), t.t_Efsw.size());
  put("tpc_wev", t.tpc_wev.data(), t.tpc_wev.size()); put("tnc_wev", t.tnc_wev.data(), t.tnc_wev.size());
  const FStorage<double>* tabs[] = {&t.tcg_racg, &t.tmr_racg, &t.tcr_gacr, &t.tmg_gacr, &t.tnr_racg,
                                       &t.tnr_gacr, &t.tcs_racs1, &t.tmr_racs1, &t.tcs_racs2, &t.tmr_racs2,
                                       &t.tcr_sacr1, &t.tms_sacr1, &t.tcr_sacr2, &t.tms_sacr2, &t.tnr_racs1,
                                       &t.tnr_racs2, &t.tnr_sacr1, &t.tnr_sacr2, &t.tpi_qrfz, &t.tni_qrfz,
                                       &t.tpg_qrfz, &t.tnr_qrfz, &t.tpi_qcfz, &t.tni_qcfz, &t.tpi_ide,
                                       &t.tps_iaus, &t.tni_iaus};
  const char* tnames[] = {"tcg_racg", "tmr_racg", "tcr_gacr", "tmg_gacr", "tnr_racg", "tnr_gacr", "tcs_racs1",
                          "tmr_racs1", "tcs_racs2", "tmr_racs2", "tcr_sacr1", "tms_sacr1", "tcr_sacr2",
                          "tms_sacr2", "tnr_racs1", "tnr_racs2", "tnr_sacr1", "tnr_sacr2", "tpi_qrfz", "tni_qrfz",
                          "tpg_qrfz", "tnr_qrfz", "tpi_qcfz", "tni_qcfz", "tpi_ide", "tps_iaus", "tni_iaus"};
  for (int i = 0; i < 27; ++i) put(tnames[i], tabs[i]->data(), tabs[i]->size());
  const radar::RadarTables& r = v.radar;
  put("xxds", r.xxDs, radar::nrbins); put("xdts", r.xdts, radar::nrbins); put("xxdg", r.xxDg, radar::nrbins);
  put("xdtg", r.xdtg, radar::nrbins); put("simpson", r.simpson, radar::nrbins + 1);
  put("k_w", &r.k_w, 1); put("pi5", &r.pi5, 1); put("lamda4", &r.lamda4, 1);
  const double mw[2] = {r.m_w_0.re, r.m_w_0.im}, mi[2] = {r.m_i_0.re, r.m_i_0.im};
  put("m_w_0", mw, 2); put("m_i_0", mi, 2);
  put("xcre", r.xcre, 4); put("xcse", r.xcse, 4); put("xcge", r.xcge, 4); put("xcrg", r.xcrg, 4);
  put("xcsg", r.xcsg, 4); put("xcgg", r.xcgg, 4);
  put("xam_r", &r.xam_r, 1); put("xbm_r", &r.xbm_r, 1); put("xmu_r", &r.xmu_r, 1); put("xobmr", &r.xobmr, 1);
  put("xam_s", &r.xam_s, 1); put("xbm_s", &r.xbm_s, 1); put("xmu_s", &r.xmu_s, 1); put("xoams", &r.xoams, 1);
  put("xobms", &r.xobms, 1); put("xocms", &r.xocms, 1);
  put("xam_g", &r.xam_g, 1); put("xbm_g", &r.xbm_g, 1); put("xmu_g", &r.xmu_g, 1); put("xoamg", &r.xoamg, 1);
  put("xobmg", &r.xobmg, 1); put("xocmg", &r.xocmg, 1);
  put("xorg2", &r.xorg2, 1); put("xosg2", &r.xosg2, 1); put("xogg2", &r.xogg2, 1);
  std::fclose(f);
}

}  // namespace cppmpas::phys::thompson
