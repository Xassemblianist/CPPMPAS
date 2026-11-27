/* SPDX-License-Identifier: GPL-3.0-or-later */
/* dump_thompson_init(path): the module variables MPAS's thompson_init leaves (module_mp_thompson and mp_radar), read
 * through their symbols (gfortran: __<module>_MOD_<name>; an ALLOCATABLE array is a descriptor whose first member is
 * the data address), written in the record format and order of CPPMPAS's dump_thompson_init
 * (src/physics/mp_thompson_init.cpp): a 32-byte name, an int64 count and the values; an int64 -1, the count and the
 * FNV-1a hash of the bytes for arrays over 100000 values. Integers are written as doubles. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define M(x) __module_mp_thompson_MOD_##x
#define RD(x) __mp_radar_MOD_##x
struct desc { double *base; };
extern double M(cce)[75], M(ccg)[75], M(ocg1)[15], M(ocg2)[15], M(cie)[7], M(cig)[7], M(oig1), M(oig2), M(obmi);
extern double M(cre)[13], M(crg)[13], M(ore1), M(org1), M(org2), M(org3), M(obmr);
extern double M(cse)[18], M(csg)[18], M(oams), M(obms), M(ocms);
extern double M(cge)[12], M(cgg)[12], M(oge1), M(ogg1), M(ogg2), M(ogg3), M(oamg), M(obmg), M(ocmg);
extern double M(t1_qr_qc), M(t1_qr_qi), M(t2_qr_qi), M(t1_qg_qc), M(t1_qs_qc), M(t1_qs_qi), M(t1_qr_ev), M(t2_qr_ev);
extern double M(t1_qs_sd), M(t2_qs_sd), M(t1_qg_sd), M(t2_qg_sd), M(t1_qs_me), M(t2_qs_me), M(t1_qg_me), M(t2_qg_me);
extern double M(sc3), M(d0i), M(xm0s), M(xm0g);
extern int M(nic1), M(nic2), M(nii2), M(nii3), M(nir2), M(nir3), M(nis2), M(nig2), M(nig3), M(niin2);
extern double M(dc)[100], M(dtc)[100], M(di)[100], M(dti)[100], M(dr)[100], M(dtr)[100], M(ds)[100], M(dts)[100];
extern double M(dg)[100], M(dtg)[100], M(t_nc)[100], M(xdx)[101];
extern struct desc M(t_efrw), M(t_efsw), M(tpc_wev), M(tnc_wev);
extern struct desc M(tcg_racg), M(tmr_racg), M(tcr_gacr), M(tmg_gacr), M(tnr_racg), M(tnr_gacr);
extern struct desc M(tcs_racs1), M(tmr_racs1), M(tcs_racs2), M(tmr_racs2), M(tcr_sacr1), M(tms_sacr1);
extern struct desc M(tcr_sacr2), M(tms_sacr2), M(tnr_racs1), M(tnr_racs2), M(tnr_sacr1), M(tnr_sacr2);
extern struct desc M(tpi_qrfz), M(tni_qrfz), M(tpg_qrfz), M(tnr_qrfz), M(tpi_qcfz), M(tni_qcfz);
extern struct desc M(tpi_ide), M(tps_iaus), M(tni_iaus);
extern double RD(xxds)[50], RD(xdts)[50], RD(xxdg)[50], RD(xdtg)[50], RD(simpson)[51], RD(k_w), RD(pi5), RD(lamda4);
extern double RD(m_w_0)[2], RD(m_i_0)[2], RD(xcre)[4], RD(xcse)[4], RD(xcge)[4], RD(xcrg)[4], RD(xcsg)[4], RD(xcgg)[4];
extern double RD(xam_r), RD(xbm_r), RD(xmu_r), RD(xobmr), RD(xam_s), RD(xbm_s), RD(xmu_s), RD(xoams), RD(xobms);
extern double RD(xocms), RD(xam_g), RD(xbm_g), RD(xmu_g), RD(xoamg), RD(xobmg), RD(xocmg), RD(xorg2), RD(xosg2);
extern double RD(xogg2);

static FILE *f;
static void put(const char *name, const double *a, size_t n) {
  char nm[32] = {0};
  strncpy(nm, name, 31);
  fwrite(nm, 1, 32, f);
  if (n > 100000) {
    uint64_t h = 1469598103934665603ULL;
    const unsigned char *b = (const unsigned char *)a;
    for (size_t i = 0; i < n * 8; ++i) h = (h ^ b[i]) * 1099511628211ULL;
    int64_t m = -1;
    double cnt = (double)n, hd;
    memcpy(&hd, &h, 8);
    fwrite(&m, 8, 1, f);
    fwrite(&cnt, 8, 1, f);
    fwrite(&hd, 8, 1, f);
  } else {
    int64_t m = (int64_t)n;
    fwrite(&m, 8, 1, f);
    fwrite(a, 8, n, f);
  }
}
static void puti(const char *name, int i) {
  double d = i;
  put(name, &d, 1);
}
#define P1(x) put(#x, &M(x), 1)
#define PA(x, n) put(#x, M(x), n)
#define PD(x, n) put(#x, M(x).base, n)
#define R1(x) put(#x, &RD(x), 1)
#define RA(x, n) put(#x, RD(x), n)

void dump_thompson_init(const char *path) {
  f = fopen(path, "wb");
  if (!f) {
    fprintf(stderr, "dump_thompson_init: cannot open %s\n", path);
    exit(1);
  }
  PA(cce, 75); PA(ccg, 75); PA(ocg1, 15); PA(ocg2, 15); PA(cie, 7); PA(cig, 7); P1(oig1); P1(oig2); P1(obmi);
  PA(cre, 13); PA(crg, 13); P1(ore1); P1(org1); P1(org2); P1(org3); P1(obmr);
  PA(cse, 18); PA(csg, 18); P1(oams); P1(obms); P1(ocms);
  PA(cge, 12); PA(cgg, 12); P1(oge1); P1(ogg1); P1(ogg2); P1(ogg3); P1(oamg); P1(obmg); P1(ocmg);
  P1(t1_qr_qc); P1(t1_qr_qi); P1(t2_qr_qi); P1(t1_qg_qc); P1(t1_qs_qc); P1(t1_qs_qi); P1(t1_qr_ev); P1(t2_qr_ev);
  P1(t1_qs_sd); P1(t2_qs_sd); P1(t1_qg_sd); P1(t2_qg_sd); P1(t1_qs_me); P1(t2_qs_me); P1(t1_qg_me); P1(t2_qg_me);
  P1(sc3); P1(d0i); P1(xm0s); P1(xm0g);
  puti("nic1", M(nic1)); puti("nic2", M(nic2)); puti("nii2", M(nii2)); puti("nii3", M(nii3)); puti("nir2", M(nir2));
  puti("nir3", M(nir3)); puti("nis2", M(nis2)); puti("nig2", M(nig2)); puti("nig3", M(nig3)); puti("niin2", M(niin2));
  PA(dc, 100); PA(dtc, 100); PA(di, 100); PA(dti, 100); PA(dr, 100); PA(dtr, 100); PA(ds, 100); PA(dts, 100);
  PA(dg, 100); PA(dtg, 100); PA(t_nc, 100); PA(xdx, 101);
  PD(t_efrw, 10000); PD(t_efsw, 10000); PD(tpc_wev, 370000); PD(tnc_wev, 370000);
  const size_t gr = 28 * 28 * 37 * 37, sr = 28 * 9 * 37 * 37, rfz = 37 * 37 * 45 * 55, cfz = 37 * 100 * 45 * 55,
               ia = 64 * 55;
  PD(tcg_racg, gr); PD(tmr_racg, gr); PD(tcr_gacr, gr); PD(tmg_gacr, gr); PD(tnr_racg, gr); PD(tnr_gacr, gr);
  PD(tcs_racs1, sr); PD(tmr_racs1, sr); PD(tcs_racs2, sr); PD(tmr_racs2, sr); PD(tcr_sacr1, sr); PD(tms_sacr1, sr);
  PD(tcr_sacr2, sr); PD(tms_sacr2, sr); PD(tnr_racs1, sr); PD(tnr_racs2, sr); PD(tnr_sacr1, sr); PD(tnr_sacr2, sr);
  PD(tpi_qrfz, rfz); PD(tni_qrfz, rfz); PD(tpg_qrfz, rfz); PD(tnr_qrfz, rfz); PD(tpi_qcfz, cfz); PD(tni_qcfz, cfz);
  PD(tpi_ide, ia); PD(tps_iaus, ia); PD(tni_iaus, ia);
  RA(xxds, 50); RA(xdts, 50); RA(xxdg, 50); RA(xdtg, 50); RA(simpson, 51); R1(k_w); R1(pi5); R1(lamda4);
  RA(m_w_0, 2); RA(m_i_0, 2);
  RA(xcre, 4); RA(xcse, 4); RA(xcge, 4); RA(xcrg, 4); RA(xcsg, 4); RA(xcgg, 4);
  R1(xam_r); R1(xbm_r); R1(xmu_r); R1(xobmr); R1(xam_s); R1(xbm_s); R1(xmu_s); R1(xoams); R1(xobms); R1(xocms);
  R1(xam_g); R1(xbm_g); R1(xmu_g); R1(xoamg); R1(xobmg); R1(xocmg); R1(xorg2); R1(xosg2); R1(xogg2);
  if (ferror(f) || fclose(f) != 0) {
    fprintf(stderr, "dump_thompson_init: write error on %s\n", path);
    exit(1);
  }
}
