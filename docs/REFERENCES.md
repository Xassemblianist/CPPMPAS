# MPAS reference releases

The tests use double-precision MPAS-Atmosphere 8.4.2 releases from [Xassemblianist/cppmpas-data](https://github.com/Xassemblianist/cppmpas-data).
MPAS was built in Debian 13 with gfortran 14.2 and the flags `-O3`, `-fdefault-real-8` and
`-fdefault-double-8`. Each run used one MPI rank. RG and RG24 use the 240 km x1.10242 mesh.
The Intel `libmvec` preload provides the reference vector `log` and `log10` on either CPU vendor.

The size and SHA-256 are for the complete compressed archive. For a split archive, verify the joined archive
against the checksum in the release's `.sha256` asset. Each archive records its runner and build in `provenance.json`.

| Tag | Case | Test coverage or status | Runner CPU vendor | Preload | Archive bytes | Archive SHA-256 |
| --- | --- | --- | --- | --- | ---: | --- |
| [`mpas842-dbl-mw-001`](https://github.com/Xassemblianist/cppmpas-data/releases/tag/mpas842-dbl-mw-001) | MW mountain wave | Initial state, step files, forecast output and restart checkpoint | AMD | None | 26,920,295 | `264e7da4a84df30762b8802f8db5c0218183e025a722d17ca0f9948855507235` |
| [`mpas842-dbl-mw-003`](https://github.com/Xassemblianist/cppmpas-data/releases/tag/mpas842-dbl-mw-003) | MW mountain wave | Initial-state fields in `initdump.nc` | AMD | None | 27,488,159 | `40decbd87b26ac699236f369c0d4952457f55530ceb65d6d96a8ef799f0ca43d` |
| [`mpas842-dbl-jw-001`](https://github.com/Xassemblianist/cppmpas-data/releases/tag/mpas842-dbl-jw-001) | JW Jablonowski–Williamson baroclinic wave | Initial state, step files, forecast output and restart checkpoint | see `provenance.json` | None | 1,330,376,772 | `f7c81b6c57d3bfaa767015b690c46c4bd81b664889934415810e8a324d430a75` |
| [`mpas842-dbl-jw-002`](https://github.com/Xassemblianist/cppmpas-data/releases/tag/mpas842-dbl-jw-002) | JW Jablonowski–Williamson baroclinic wave | Initial-state fields in `initdump.nc` | AMD | None | 1,651,734,882 | `11d679fa0e9c4ac62d9145e6d04f2a42f018828b9e5d87d32c6560808e61e74f` |
| [`mpas842-dbl-sc-001`](https://github.com/Xassemblianist/cppmpas-data/releases/tag/mpas842-dbl-sc-001) | SC supercell with Kessler microphysics | Initial state, step files and restart checkpoint | see `provenance.json` | None | 403,858,595 | `ea47e824a970973f05e907e0ef1cc912ffe816b60df2e59df3360bceca674f39` |
| [`mpas842-dbl-sc-002`](https://github.com/Xassemblianist/cppmpas-data/releases/tag/mpas842-dbl-sc-002) | SC supercell with Kessler microphysics | Initial-state fields in `initdump.nc` and 30 min forecast output | AMD | None | 406,492,518 | `a44918ffb917237ac8db8c4c1dc0275efedd050da7546824c1fd3e0831493c50` |
| [`mpas842-dbl-rg-004`](https://github.com/Xassemblianist/cppmpas-data/releases/tag/mpas842-dbl-rg-004) | RG real-data forecast | Step files and 2 h `mesoscale_reference` forecast with history, diagnostics and restart | Intel | Intel `libmvec` | 1,622,029,020 | `59aeab7f5bbfda2cb4aa0b48a6f3d9a2c3686d90f54e76d3945272b1e60b2fba` |
| [`mpas842-dbl-rg24-001`](https://github.com/Xassemblianist/cppmpas-data/releases/tag/mpas842-dbl-rg24-001) | RG24 real-data forecast | 24 h `mesoscale_reference` forecast with history, diagnostics and restart | Intel | Intel `libmvec` | 2,210,029,328 | `59e1eaf8d2d42a4b9e17d62567a0283de349260c2fc61abcca3119264331a372` |
| [`mpas842-dbl-rgr-001`](https://github.com/Xassemblianist/cppmpas-data/releases/tag/mpas842-dbl-rgr-001) | RGR restart run | RG's 2 h forecast, then a `mesoscale_reference` restart run from its 01:10 checkpoint to 02:50 with history, diagnostics and restart | Intel | Intel `libmvec` | 1,352,722,570 | `c4720917afbbfed2e2080187fb2378eb07196c1d42562ac2c426f342debc8264` |
| [`mpas842-dbl-rg-001`](https://github.com/Xassemblianist/cppmpas-data/releases/tag/mpas842-dbl-rg-001) | RG real-data forecast | Superseded by `mpas842-dbl-rg-004` | AMD | None | 1,621,879,539 | `5271eaf1f0c5c16a63de5966c14fac33cd87547f69b99a6378cccca7b6c690a9` |
| [`mpas842-dbl-rg-002`](https://github.com/Xassemblianist/cppmpas-data/releases/tag/mpas842-dbl-rg-002) | RG real-data forecast | Superseded by `mpas842-dbl-rg-004` | AMD | None | 1,621,879,453 | `afc69bf9276fd753d8b2d48a090bb8357121337138f6e76a676020902afd029c` |
| [`mpas842-dbl-rg-003`](https://github.com/Xassemblianist/cppmpas-data/releases/tag/mpas842-dbl-rg-003) | RG real-data forecast | Superseded by `mpas842-dbl-rg-004` | Intel | None | 1,621,890,528 | `e683d9641141de4c309b834a4aca8a59f04a3b98ed2b15a924fdc072424e7024` |

The superseded RG releases used the runner's system `libmvec`. Their results depend on the runner CPU vendor.

Extract the releases used by the tests into separate tag directories under `CPPMPAS_REFS_DIR`, as named in
[tests/CMakeLists.txt](../tests/CMakeLists.txt). The data repository documents the download and verification steps.
[Fetch a reference](https://github.com/Xassemblianist/cppmpas-data#fetch-a-reference)
