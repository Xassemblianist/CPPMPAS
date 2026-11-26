! SPDX-License-Identifier: GPL-3.0-or-later
! The RRTMG tables as MPAS's own initialisation leaves them: rrtmg_initlw_forMPAS and rrtmg_initsw_forMPAS
! (mpas_atmphys_rrtmg_lwinit.F:43-58, mpas_atmphys_rrtmg_swinit.F), which read RRTMG_LW_DATA.DBL / RRTMG_SW_DATA.DBL
! from the working directory and call rrtmg_lw_ini(cp) / rrtmg_sw_ini(cp), written as memory images in the layout of
! CPPWRF's dump programs (CPPWRF 7fac5bf tests/phys/ra_rrtmg_lw/dump_tables.f90, ra_rrtmg_sw/dump_sw_tables.f90, the
! WRITE lists copied unchanged): radiation_dump_tables <lw.bin> <sw.bin>. No arithmetic.
PROGRAM radiation_dump_tables
  USE mpas_derived_types, only: dm_info
  USE mpas_atmphys_rrtmg_lwinit, only: rrtmg_initlw_forMPAS
  USE mpas_atmphys_rrtmg_swinit, only: rrtmg_initsw_forMPAS
  IMPLICIT NONE
  TYPE(dm_info) :: dminfo
  CHARACTER(len=512) :: flw, fsw
  CALL get_command_argument(1, flw)
  CALL get_command_argument(2, fsw)
  CALL rrtmg_initlw_forMPAS(dminfo)
  CALL rrtmg_initsw_forMPAS(dminfo)
  CALL write_lw(trim(flw))
  CALL write_sw(trim(fsw))
CONTAINS
  SUBROUTINE write_lw(path)
    USE rrlw_con
    USE rrlw_wvn
    USE rrlw_tbl
    USE rrlw_cld
    USE rrlw_ref
    CHARACTER(len=*), INTENT(IN) :: path
    INTEGER :: uout, ios
    OPEN(NEWUNIT=uout, FILE=path, ACCESS='STREAM', FORM='UNFORMATTED', CONVERT='LITTLE_ENDIAN', STATUS='REPLACE', IOSTAT=ios)
    IF (ios /= 0) STOP 'cannot open the LW tables file'
  ! Lookup tables from rrlw_tbl: tau_tbl(0:10000), exp_tbl(0:10000), tfn_tbl(0:10000)
  WRITE(uout) tau_tbl
  WRITE(uout) exp_tbl
  WRITE(uout) tfn_tbl

  ! Spectral info from rrlw_wvn: totplnk(181, 16), totplk16(181), delwave(16), wavenum1(16), wavenum2(16)
  WRITE(uout) totplnk
  WRITE(uout) totplk16
  WRITE(uout) delwave
  WRITE(uout) wavenum1
  WRITE(uout) wavenum2
  WRITE(uout) ngc
  WRITE(uout) ngs
  WRITE(uout) ngn
  WRITE(uout) ngb
  WRITE(uout) ngm
  WRITE(uout) wt
  WRITE(uout) rwgt

  ! Cloud tables from rrlw_cld: abscld1, absice0, absice1, absice2, absice3, absliq0, absliq1
  WRITE(uout) abscld1
  WRITE(uout) absice0
  WRITE(uout) absice1
  WRITE(uout) absice2
  WRITE(uout) absice3
  WRITE(uout) absliq0
  WRITE(uout) absliq1

  ! Reference profiles from rrlw_ref: preflog(59), tref(59), chi_mls(7, 59)
  WRITE(uout) preflog
  WRITE(uout) tref
  WRITE(uout) chi_mls

  CALL dump_b01(uout)
  CALL dump_b02(uout)
  CALL dump_b03(uout)
  CALL dump_b04(uout)
  CALL dump_b05(uout)
  CALL dump_b06(uout)
  CALL dump_b07(uout)
  CALL dump_b08(uout)
  CALL dump_b09(uout)
  CALL dump_b10(uout)
  CALL dump_b11(uout)
  CALL dump_b12(uout)
  CALL dump_b13(uout)
  CALL dump_b14(uout)
  CALL dump_b15(uout)
  CALL dump_b16(uout)

    CLOSE(uout)
  END SUBROUTINE write_lw
  SUBROUTINE write_sw(path)
    USE rrsw_con
    USE rrsw_wvn
    USE rrsw_tbl
    USE rrsw_cld
    USE rrsw_ref
    USE rrsw_aer
    CHARACTER(len=*), INTENT(IN) :: path
    INTEGER :: uout, ios
    OPEN(NEWUNIT=uout, FILE=path, ACCESS='STREAM', FORM='UNFORMATTED', CONVERT='LITTLE_ENDIAN', STATUS='REPLACE', IOSTAT=ios)
    IF (ios /= 0) STOP 'cannot open the SW tables file'
  ! 1. Lookup tables from rrsw_tbl: exp_tbl(0:10000)
  WRITE(uout) exp_tbl

  ! 2. Spectral info from rrsw_wvn
  WRITE(uout) wavenum1
  WRITE(uout) wavenum2
  WRITE(uout) delwave
  WRITE(uout) ngc
  WRITE(uout) ngs
  WRITE(uout) ngn
  WRITE(uout) ngb
  WRITE(uout) ngm
  WRITE(uout) wt
  WRITE(uout) rwgt

  ! 3. Cloud tables from rrsw_cld
  WRITE(uout) extliq1
  WRITE(uout) ssaliq1
  WRITE(uout) asyliq1
  WRITE(uout) extice2
  WRITE(uout) ssaice2
  WRITE(uout) asyice2
  WRITE(uout) extice3
  WRITE(uout) ssaice3
  WRITE(uout) asyice3
  WRITE(uout) fdlice3
  WRITE(uout) abari
  WRITE(uout) bbari
  WRITE(uout) cbari
  WRITE(uout) dbari
  WRITE(uout) ebari
  WRITE(uout) fbari

  ! 4. Reference profiles from rrsw_ref: preflog(59), tref(59)
  WRITE(uout) preflog
  WRITE(uout) tref

  ! 5. Aerosol properties from rrsw_aer
  WRITE(uout) rsrtaua
  WRITE(uout) rsrpiza
  WRITE(uout) rsrasya

  ! 6. Band absorption tables
  CALL dump_s16(uout)
  CALL dump_s17(uout)
  CALL dump_s18(uout)
  CALL dump_s19(uout)
  CALL dump_s20(uout)
  CALL dump_s21(uout)
  CALL dump_s22(uout)
  CALL dump_s23(uout)
  CALL dump_s24(uout)
  CALL dump_s25(uout)
  CALL dump_s26(uout)
  CALL dump_s27(uout)
  CALL dump_s28(uout)
  CALL dump_s29(uout)

    CLOSE(uout)
  END SUBROUTINE write_sw
  SUBROUTINE dump_b01(u)
    USE rrlw_kg01
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, ka_mn2, kb_mn2, selfref, forref
  END SUBROUTINE dump_b01

  SUBROUTINE dump_b02(u)
    USE rrlw_kg02
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, selfref, forref
  END SUBROUTINE dump_b02

  SUBROUTINE dump_b03(u)
    USE rrlw_kg03
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, ka_mn2o, kb_mn2o, selfref, forref
  END SUBROUTINE dump_b03

  SUBROUTINE dump_b04(u)
    USE rrlw_kg04
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, selfref, forref
  END SUBROUTINE dump_b04

  SUBROUTINE dump_b05(u)
    USE rrlw_kg05
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, ka_mo3, ccl4, selfref, forref
  END SUBROUTINE dump_b05

  SUBROUTINE dump_b06(u)
    USE rrlw_kg06
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, ka, ka_mco2, cfc11adj, cfc12, selfref, forref
  END SUBROUTINE dump_b06

  SUBROUTINE dump_b07(u)
    USE rrlw_kg07
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, ka_mco2, kb_mco2, selfref, forref
  END SUBROUTINE dump_b07

  SUBROUTINE dump_b08(u)
    USE rrlw_kg08
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, ka_mco2, kb_mco2, ka_mn2o, kb_mn2o, ka_mo3, cfc12, cfc22adj, selfref, forref
  END SUBROUTINE dump_b08

  SUBROUTINE dump_b09(u)
    USE rrlw_kg09
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, ka_mn2o, kb_mn2o, selfref, forref
  END SUBROUTINE dump_b09

  SUBROUTINE dump_b10(u)
    USE rrlw_kg10
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, selfref, forref
  END SUBROUTINE dump_b10

  SUBROUTINE dump_b11(u)
    USE rrlw_kg11
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, ka_mo2, kb_mo2, selfref, forref
  END SUBROUTINE dump_b11

  SUBROUTINE dump_b12(u)
    USE rrlw_kg12
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, ka, selfref, forref
  END SUBROUTINE dump_b12

  SUBROUTINE dump_b13(u)
    USE rrlw_kg13
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, ka_mco2, ka_mco, kb_mo3, selfref, forref
  END SUBROUTINE dump_b13

  SUBROUTINE dump_b14(u)
    USE rrlw_kg14
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, selfref, forref
  END SUBROUTINE dump_b14

  SUBROUTINE dump_b15(u)
    USE rrlw_kg15
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, ka, ka_mn2, selfref, forref
  END SUBROUTINE dump_b15

  SUBROUTINE dump_b16(u)
    USE rrlw_kg16
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, selfref, forref
  END SUBROUTINE dump_b16

  SUBROUTINE dump_s16(u)
    USE rrsw_kg16
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, kb, selfref, forref, sfluxref, rayl, strrat1, layreffr
  END SUBROUTINE dump_s16

  SUBROUTINE dump_s17(u)
    USE rrsw_kg17
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, kb, selfref, forref, sfluxref, rayl, strrat, layreffr
  END SUBROUTINE dump_s17

  SUBROUTINE dump_s18(u)
    USE rrsw_kg18
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, kb, selfref, forref, sfluxref, rayl, strrat, layreffr
  END SUBROUTINE dump_s18

  SUBROUTINE dump_s19(u)
    USE rrsw_kg19
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, kb, selfref, forref, sfluxref, rayl, strrat, layreffr
  END SUBROUTINE dump_s19

  SUBROUTINE dump_s20(u)
    USE rrsw_kg20
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, kb, selfref, forref, sfluxref, absch4, rayl, layreffr
  END SUBROUTINE dump_s20

  SUBROUTINE dump_s21(u)
    USE rrsw_kg21
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, kb, selfref, forref, sfluxref, rayl, strrat, layreffr
  END SUBROUTINE dump_s21

  SUBROUTINE dump_s22(u)
    USE rrsw_kg22
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, kb, selfref, forref, sfluxref, rayl, strrat, layreffr
  END SUBROUTINE dump_s22

  SUBROUTINE dump_s23(u)
    USE rrsw_kg23
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, selfref, forref, sfluxref, rayl, givfac, layreffr
  END SUBROUTINE dump_s23

  SUBROUTINE dump_s24(u)
    USE rrsw_kg24
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, kb, selfref, forref, sfluxref, abso3a, abso3b, rayla, raylb, strrat, layreffr
  END SUBROUTINE dump_s24

  SUBROUTINE dump_s25(u)
    USE rrsw_kg25
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, sfluxref, abso3a, abso3b, rayl, layreffr
  END SUBROUTINE dump_s25

  SUBROUTINE dump_s26(u)
    USE rrsw_kg26
    INTEGER, INTENT(IN) :: u
    WRITE(u) sfluxref, rayl
  END SUBROUTINE dump_s26

  SUBROUTINE dump_s27(u)
    USE rrsw_kg27
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, kb, sfluxref, rayl, scalekur, layreffr
  END SUBROUTINE dump_s27

  SUBROUTINE dump_s28(u)
    USE rrsw_kg28
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, kb, sfluxref, rayl, strrat, layreffr
  END SUBROUTINE dump_s28

  SUBROUTINE dump_s29(u)
    USE rrsw_kg29
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, kb, selfref, forref, sfluxref, absh2o, absco2, rayl, layreffr
  END SUBROUTINE dump_s29

END PROGRAM radiation_dump_tables
