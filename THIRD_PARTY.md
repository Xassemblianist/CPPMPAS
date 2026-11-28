# Third-party components and notices

CPPMPAS's project-owned code is available under GPL-3.0-or-later. The terms below govern material from other projects. Keep the applicable notices with source distributions and in accompanying materials for binary distributions where the upstream terms require them.

The legacy RRTMG and Numerical Recipes entries identify unresolved redistribution rights in the current source tree. Confirm the relevant grants or replace the affected code and data before distributing a combined GPL or paid build.

## Source components

| Component | CPPMPAS paths and origin | Terms |
| --- | --- | --- |
| MPAS-Atmosphere v8.4.2 and MPAS framework | Ports under `src/core_atmosphere/`, `src/operators/`, and parts of `src/core/`, `src/io/` and `src/physics/`. The two generated files `src/core/namelist_table.inc` and `src/io/registry_atm_table.inc` include MPAS Registry material. | LANS/UCAR BSD-style terms, including U.S. Government rights, modification marking, notice retention and no endorsement without permission. See the full notice below. [MPAS v8.4.2 license](https://github.com/MPAS-Dev/MPAS-Model/blob/v8.4.2/LICENSE). |
| NCAR MMM-physics | The reference pin is `a4baf7f` (`ci/data.tsv:11`). Ports include GWDO, MYNN PBL and surface layer, YSU, New Tiedtke, WSM6, radar and revised MM5 surface-layer code under `src/physics/`. | BSD-3-Clause, Copyright (c) 2022, National Center for Atmospheric Research. See the full notice below. The gamma calculation in `mp_radar.hpp` also has Numerical Recipes ancestry. [MMM-physics license](https://github.com/NCAR/MMM-physics/blob/main/LICENSE). |
| Weather Research and Forecasting model (WRF) physics | MPAS v8.4.2 `physics_wrf` modules and WRF-lineage ports carried over from CPPWRF. Portions include Kessler, Noah land/ice and drivers, WRF-lineage physics wrappers, MPAS's Grell–Freitas and MPAS's WRF-v4.1.4 Thompson. See `SHARED_WITH_CPPWRF.md:18–19,22–23,26–40`. | NCAR/UCAR's WRF Public Domain Notice and disclaimer below. Separately licensed code has its own terms, including AER RRTMG and Numerical Recipes. Separate CAM terms for the exact `getfactors` port are not checked. [UCAR WRF Public Domain Notice](https://www2.mmm.ucar.edu/wrf/site/access_code/public_domain_notice.html). |
| Shared CPPWRF code | `SHARED_WITH_CPPWRF.md` records the source commits, shared owner-written GPU/build helpers and physics kernels. | The owner-written portions follow CPPMPAS's project license here. The license of each upstream portion follows its entry in this document. |
| GNU C Library `libm` and `libmvec` | GNU C Library 2.41 and Debian 13 x86-64 builds. `src/core/{glibc_math.hpp,glibc_math_tables.hpp,libm.hpp,libmvec_x86.hpp,libmvec_tables.hpp}` contain modified source or machine-code translations, extracted tables and C++ wrappers. `libmvec_tables.hpp` also includes reciprocal outputs measured on an Intel Core i7-10750H. | LGPL-2.1-or-later for the covered GNU C Library portions. Retain the original notices and mark the changes and dates. Include a copy of the LGPL, glibc's applicable notices, corresponding source and any materials needed for the applicable relinking or replacement route when distributing covered binaries. The embedded CPU/GPU math is subject to these terms. [glibc 2.41 LGPL text](https://github.com/bminor/glibc/blob/glibc-2.41/COPYING.LIB), [glibc 2.41 third-party notices](https://github.com/bminor/glibc/blob/glibc-2.41/LICENSES), [LGPL 2.1 §§1–6](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html). |
| Sun fdlibm inside glibc | Portions of `src/core/glibc_math.hpp` port sources with Sun notices, including `e_log10.c`, `s_tanh.c`, `s_expm1.c`, `e_sinh.c`, `s_erf.c`, `e_lgamma_r.c` and `s_log1p.c`. | Preserve the SunPro permission notice below. [Example upstream file](https://github.com/bminor/glibc/blob/glibc-2.41/sysdeps/ieee754/dbl-64/e_log10.c#L4-L9). |
| Atmospheric & Environmental Research (AER) RRTMG LW/SW | RRTMG stage and kernel headers under `src/physics/rrtmg/`, the RRTMG-lineage portions of `src/physics/module_ra_rrtmg.*`, and the external RRTMG coefficients and CPPWRF table dumps. See `SHARED_WITH_CPPWRF.md:25–29`. `kissvec.hpp` includes the KISS PRNG marked public domain in the upstream Fortran, alongside McICA code. | Legacy WRF/MPAS sources give AER notices with a no-sale condition. Retain each applicable notice on each copy. Scope of the legacy grant to the separate driver/loader files and dumps is subject to upstream verification. AER also publishes BSD-3-Clause RRTMG LW/SW releases. Applicability of those grants to the exact WRF/MPAS portions and input tables is not checked. [AER RRTMG LW](https://github.com/AER-RC/RRTMG_LW/blob/master/LICENSE.md), [AER RRTMG SW](https://github.com/AER-RC/RRTMG_SW/blob/master/LICENSE.md), [MPAS RRTMG LW source](https://github.com/MPAS-Dev/MPAS-Model/blob/v8.4.2/src/core_atmosphere/physics/physics_wrf/module_ra_rrtmg_lw.F#L2056-L2066), [MPAS RRTMG SW source](https://github.com/MPAS-Dev/MPAS-Model/blob/v8.4.2/src/core_atmosphere/physics/physics_wrf/module_ra_rrtmg_sw.F#L1360-L1370). |
| Numerical Recipes Software 2.02 | `src/physics/mpas_atmphys_functions.hpp:13–31` and `src/physics/mp_radar.hpp:86–102`. Upstream MPAS and MMM-physics mark the translated gamma routine. | (C) Copr. 1986-92 Numerical Recipes Software 2.02. Permission to redistribute these translations is unverified. [MPAS source](https://github.com/MPAS-Dev/MPAS-Model/blob/v8.4.2/src/core_atmosphere/physics/mpas_atmphys_functions.F#L96-L120), [Publisher's licensing terms](https://www.numerical.recipes/licenseinfo.html). |

George H. Bryan, NCAR Mesoscale and Microscale Meteorology Division, authored the ported `getcape` routine in the MPAS convective diagnostics. The upstream routine is provided without warranty. [MPAS source attribution](https://github.com/MPAS-Dev/MPAS-Model/blob/v8.4.2/src/core_atmosphere/diagnostics/mpas_convective_diagnostics.F#L604-L615).

## Cases, tables and external build inputs

- NCAR's Schaer mountain-wave, supercell and Jablonowski–Williamson baroclinic-wave inputs are fetched for tests. The tracked MW restart namelist adapts the NCAR forecast namelist. The tracked RGR restart namelist comes from `mpas842-dbl-rgr-001`. The redistribution terms for the case archives and the underlying RGR configuration are not checked. [NCAR idealized test cases](https://www2.mmm.ucar.edu/projects/mpas/site/downloads/idealized.html).
- `docs/REFERENCES.md`, `ci/data.tsv` and `tools/colab/drive_manifest.tsv` identify reference releases, runtime tables and their checksums. The separate input and build releases include MPAS source, MPAS-Data tables, Thompson and AER RRTMG coefficients, and mesh and forcing inputs. The hashes, filenames and links do not establish redistribution rights. Check each provider's terms before redistributing a table, archive or reference output. [MPAS-Data](https://github.com/MPAS-Dev/MPAS-Data).
- The reference build also has a NOAA-GSL UGWP dependency. The CPPMPAS test stub at `tests/cu_gwdo_column/stubs.F90:255–285` copies declarations and comments from the MPAS wrapper. The license for any separate NOAA material in those declarations or in distributed upstream archives is not checked. [NOAA-GSL UGWP](https://github.com/NOAA-GSL/UGWP).
- CPPMPAS uses an externally installed NetCDF-C library. Builds can also use an OpenMP runtime and the NVIDIA CUDA Toolkit. These retain their own terms. A distributed CI image or binary bundle can include Debian package contents, additional runtime libraries or a separately built glibc-derived test preload. Include the applicable notices and fulfill any source or redistribution terms for what the bundle actually contains. [NetCDF-C copyright](https://docs.unidata.ucar.edu/netcdf-c/current/copyright.html), [GCC Runtime Library Exception](https://gcc.gnu.org/runtime-library-exception.html), [NVIDIA CUDA Toolkit terms](https://docs.nvidia.com/cuda/eula/contents.html).

## License and attribution texts

### MPAS v8.4.2

```text
Copyright (c) 2013-2019,  Los Alamos National Security, LLC (LANS) (Ocean: LA-CC-13-047;
Land Ice: LA-CC-13-117) and the University Corporation for Atmospheric Research (UCAR).

All rights reserved. 

LANS is the operator of the Los Alamos National Laboratory under Contract No.
DE-AC52-06NA25396 with the U.S. Department of Energy.  UCAR manages the National
Center for Atmospheric Research under Cooperative Agreement ATM-0753581 with the
National Science Foundation.  The U.S. Government has rights to use, reproduce,
and distribute this software.  NO WARRANTY, EXPRESS OR IMPLIED IS OFFERED BY
LANS, UCAR OR THE GOVERNMENT AND NONE OF THEM ASSUME ANY LIABILITY FOR THE USE
OF THIS SOFTWARE.  If software is modified to produce derivative works, such
modified software should be clearly marked, so as not to confuse it with the
version available from LANS and UCAR.

Additionally, redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1) Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3) None of the names of LANS, UCAR or the names of its contributors, if any, may
be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

### NCAR MMM-physics

```text
BSD 3-Clause License

Copyright (c) 2022, National Center for Atmospheric Research
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

### WRF Public Domain Notice as supplied with MPAS v8.4.2

```text
------------------------
WRF PUBLIC DOMAIN NOTICE
------------------------

WRF was developed at the National Center for Atmospheric Research
(NCAR) which is operated by the University Corporation for
Atmospheric Research (UCAR). NCAR and UCAR make no proprietary
claims, either statutory or otherwise, to this version and
release of WRF and consider WRF to be in the public domain for
use by any person or entity for any purpose without any fee or
charge. UCAR requests that any WRF user include this notice on
any partial or full copies of WRF. WRF is provided on an "AS
IS" basis and any warranties, either express or implied,
including but not limited to implied warranties of
non-infringement, originality, merchantability and fitness for a
particular purpose, are disclaimed. In no event shall
UCAR be liable for any damages, whatsoever, whether direct,
indirect, consequential or special, that arise out of or in
connection with the access, use or performance of WRF, including
infringement actions.
```

### Legacy AER RRTMG

2002–2008 modules

```text
Copyright 2002-2008, Atmospheric & Environmental Research, Inc. (AER).
This software may be used, copied, or redistributed as long as it is
not sold and this copyright notice is reproduced on each copy made.
This model is provided as is without any express or implied warranties.
(http://www.rtweb.aer.com/)
```

Longwave cloud optical properties and optical depth

```text
Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).
This software may be used, copied, or redistributed as long as it is
not sold and this copyright notice is reproduced on each copy made.
This model is provided as is without any express or implied warranties.
(http://www.rtweb.aer.com/)
```

McICA cloud subcolumns

```text
Copyright 2006-2008, Atmospheric & Environmental Research, Inc. (AER).
This software may be used, copied, or redistributed as long as it is
not sold and this copyright notice is reproduced on each copy made.
This model is provided as is without any express or implied warranties.
(http://www.rtweb.aer.com/)
```

### Sun fdlibm

```text
Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.

Developed at SunPro, a Sun Microsystems, Inc. business.
Permission to use, copy, modify, and distribute this
software is freely granted, provided that this notice
is preserved.
```
